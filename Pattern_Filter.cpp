#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <thread>
#include <future>
#include <numeric> // For std::accumulate
#include <algorithm>
#include <stdexcept>
#include <filesystem> // Requires C++17
#include <random>
#include <unistd.h> // For getpid()
#include <cstdio>   // For std::remove, std::rename
#include <cstdlib>  // For std::mkdtemp, std::getenv
#include <cstring>  // For strerror, strcpy, std::strlen
#include <cerrno>   // For errno
#include <queue>
#include <atomic>
#include <chrono>
#include <iomanip>
#include <memory>   // For unique_ptr
#include <execution> // For potential std::execution::par (C++17)
#include <map>
#include <unordered_set> // For R2 extraction ID lookup
#include <mutex> // For thread safety during result collection
#include <cmath> // For std::ceil

// Assuming gzstream library is available and header is in include path
#include "gzstream.h"

// For Zlib operations during sorting phase and R2 extraction/concatenation
#include <zlib.h>

// For argument parsing (assuming cxxopts.hpp is available)
#include "cxxopts.hpp"

namespace fs = std::filesystem;

// --- Configuration for Sorting & Extraction ---
const size_t GZ_BUFFER_SIZE = 256 * 1024;
const size_t CONCAT_BUFFER_SIZE = 256 * 1024; // Buffer for concatenating gz files
// const unsigned int MAX_MERGE_FILES = 100; // Note: This isn't actively used in the current merge logic but kept for context
// --- End Configuration ---

// =============================================
// == Params Struct Definition (Moved Earlier) ==
// =============================================
struct Params {
    fs::path r1_input;
    fs::path r2_input;
    fs::path r1_output_dir;
    fs::path r2_output_dir;
    fs::path hamming_filterout_dir;
    fs::path base_composition_filterout_dir;
    std::string target_seq;
    int seq_length = 0;     // For Hamming region
    int start_pos = 0;      // For Hamming region (0-based)
    int threshold = 0;      // Hamming distance threshold
    double base_composition_threshold = 0.0; // For base comp region
    int base_comp_start_pos = 0; // MODIFIED: Now a parameter
    int base_comp_length = 0;   // MODIFIED: Now a parameter
    int num_threads = 1;
    long long chunk_read_count = 0;
    fs::path main_temp_dir_base;
    size_t sort_memory_limit_mb = 0;
    int sort_threads = 1;
};


// === Forward Declarations for Helper Functions ===
void handle_incomplete_record(gzFile infile, const fs::path& chunk_path, const std::string& id_line, const std::string& missing_part);
void handle_read_error(igzstream& stream, const fs::path& filepath, long long record_num, const std::string& expected_part);
fs::path extract_r2_reads_parallel(
    const fs::path& sorted_r1_passed_path,
    const fs::path& original_r2_input_path,
    const fs::path& final_r2_output_path,
    const fs::path& main_temp_dir_base,
    const Params& params);
// MODIFIED: Main sorting function now coordinates parallel workers
bool sort_and_merge_chunk_category(
    const std::vector<fs::path>& input_gz_chunks,
    const fs::path& final_output_gz_path,
    const fs::path& sort_temp_dir_base,
    const std::string& category_name,
    size_t sort_memory_limit_bytes,
    int num_sort_threads,
    bool enable_parallel_sort_alg);
// NEW: Worker function for parallel sorting phase 1
struct SortWorkerResult {
    bool success = false;
    std::vector<fs::path> created_temp_files;
    long long records_processed = 0;
};
SortWorkerResult process_sort_sub_chunk(
    const std::vector<fs::path>& sub_input_gz_chunks,
    const fs::path& category_sort_temp_path,
    const std::string& category_name,
    size_t sort_memory_limit_bytes,
    bool enable_parallel_sort_alg,
    int worker_id);
// Merging function remains largely the same (sequential)
bool merge_plain_text_chunks_to_gz(const std::vector<fs::path>& plain_chunk_files, const fs::path& final_output_gz_path);
bool concatenate_gz_files(const std::vector<fs::path>& input_files, const fs::path& output_file);
fs::path process_r2_extraction_chunk(
    const std::vector<std::string>& r2_chunk_records,
    const std::unordered_set<std::string>& r1_passed_ids_cref,
    const fs::path& r2_extract_temp_dir,
    int chunk_index);
// =============================================


// ========================================================================
// == RAII Temporary Directory (Used for Filtering and Sorting Temp Dirs) ==
// ========================================================================
class TempDir {
    fs::path path_;
    bool owns_ = false;

public:
    TempDir(const fs::path& base_dir = fs::current_path(), const std::string& prefix = "fastq_proc_temp_") {
        fs::path effective_base_dir = base_dir;
        if (effective_base_dir.empty()) {
            effective_base_dir = fs::current_path();
             std::cerr << "Warning: TempDir received empty base_dir, using current_path(): " << effective_base_dir << std::endl;
        }

        std::error_code ec;
        if (!fs::exists(effective_base_dir, ec)) {
            if (!fs::create_directories(effective_base_dir, ec)) {
                 throw std::runtime_error("Failed to create specified base directory '" + effective_base_dir.string() + "' for temporary files: " + ec.message());
            }
             std::cout << "Info: Created base directory for temporary files: " << effective_base_dir << std::endl;
        } else if (!fs::is_directory(effective_base_dir, ec)) {
            throw std::runtime_error("Specified base path for temporary files is not a directory: " + effective_base_dir.string());
        }

        std::error_code ec_base_exists;
        if (!fs::exists(effective_base_dir, ec_base_exists) || !fs::is_directory(effective_base_dir, ec_base_exists)) {
             throw std::runtime_error("Base temp directory does not exist or is not a directory after check: " + effective_base_dir.string());
        }

        std::string template_str = (effective_base_dir / (prefix + "XXXXXX")).string();
        std::vector<char> template_cstr(template_str.begin(), template_str.end());
        template_cstr.push_back('\0');

        char* created_path_cstr = mkdtemp(template_cstr.data());
        if (!created_path_cstr) {
             int err_no = errno;
             throw std::runtime_error("Failed to create temporary directory using template '" + template_str + "'. Check permissions and path validity. Error [" + std::to_string(err_no) + "]: " + std::string(strerror(err_no)));
        }
        path_ = created_path_cstr;
        owns_ = true;
        // std::cout << "Debug: TempDir created at: " << path_ << std::endl;
    }

    ~TempDir() {
        if (owns_ && !path_.empty()) {
             // std::cout << "Debug: Cleaning up TempDir: " << path_ << std::endl;
             std::error_code ec_exists;
             if (fs::exists(path_, ec_exists)) {
                try {
                    std::error_code ec_remove;
                    // More robust check for valid temp dir names before removing
                    std::string filename = path_.filename().string();
                    bool safe_to_remove = filename.find("fastq_proc_temp_") != std::string::npos ||
                                          filename.find("filter_sort_main_") != std::string::npos ||
                                          filename.find("sort_temp_") != std::string::npos ||      // Sort category temp dir
                                          filename.find("r2_extract_chunks_") != std::string::npos;

                    if (path_.has_filename() && safe_to_remove)
                    {
                         fs::remove_all(path_, ec_remove);
                         if (ec_remove) {
                             std::cerr << "Warning: Failed to completely remove temporary directory '" << path_.string() << "': " << ec_remove.message() << std::endl;
                         }
                    } else {
                         std::cerr << "Warning: Skipping removal of potential unsafe temporary directory path: " << path_.string() << std::endl;
                    }
                } catch (const fs::filesystem_error& e) {
                    std::cerr << "Warning: Filesystem error during temporary directory cleanup '" << path_.string() << "': " << e.what()
                              << " Path1: " << e.path1() << " Path2: " << e.path2() << std::endl;
                } catch (const std::exception& e) {
                    std::cerr << "Warning: An error occurred during temporary directory cleanup '" << path_.string() << "': " << e.what() << std::endl;
                } catch (...) {
                    std::cerr << "Warning: An unknown error occurred during temporary directory cleanup '" << path_.string() << "'." << std::endl;
                }
            }
        }
    }

    TempDir(const TempDir&) = delete;
    TempDir& operator=(const TempDir&) = delete;
    TempDir(TempDir&& other) noexcept : path_(std::move(other.path_)), owns_(other.owns_) {
        other.owns_ = false;
    }
    TempDir& operator=(TempDir&& other) noexcept {
        if (this != &other) {
            if (owns_ && !path_.empty()) {
                 std::error_code ec;
                 if(fs::exists(path_, ec)) { try { fs::remove_all(path_); } catch (...) { /* Log */ } }
            }
            path_ = std::move(other.path_);
            owns_ = other.owns_;
            other.owns_ = false;
        }
        return *this;
    }
    const fs::path& getPath() const { return path_; }
    void release() { owns_ = false; }
};


// ========================================================================
// == Filtering Stage Structures and Functions (Unchanged)              ==
// ========================================================================

inline int hamming_distance(const std::string& s1, const std::string& s2) {
    int dist = 0;
    size_t len = std::min(s1.length(), s2.length()); // Handle different lengths safely
    for (size_t i = 0; i < len; ++i) {
        if (s1[i] != s2[i]) {
            dist++;
        }
    }
    // If lengths differ, count remaining characters in longer string as mismatches?
    // Standard Hamming distance is only for strings of equal length.
    // Assuming target_seq length == seq_length from params ensures s1 and s2 from substr are same length.
    // If seq1.length() < params.start_pos + params.seq_length, this check won't be reached anyway.
    return dist;
}

inline double calculate_base_composition(const std::string& seq) {
    if (seq.empty()) return 0.0;
    int acgt_count = 0;
    for (char base : seq) {
        char upper_base = std::toupper(base);
        if (upper_base == 'A' || upper_base == 'C' || upper_base == 'G' || upper_base == 'T') {
            acgt_count++;
        }
    }
    return static_cast<double>(acgt_count) / seq.length();
}

struct UnsortedChunkOutputPaths {
    fs::path passed;
    fs::path hamming_rejected;
    fs::path basecomp_rejected;
    bool success = false;
};

UnsortedChunkOutputPaths process_filter_chunk(
    const std::vector<std::string>& reads_r1,
    const fs::path& filter_temp_dir_path, // This is the *main* temp dir path
    int chunk_index,
    const Params& params)
{
    UnsortedChunkOutputPaths output_paths;
    output_paths.success = false;

    try {
        std::error_code ec_exists;
        if (!fs::exists(filter_temp_dir_path, ec_exists) || !fs::is_directory(filter_temp_dir_path, ec_exists)) {
             throw std::runtime_error("Main filter temp dir does not exist or is not a directory: " + filter_temp_dir_path.string());
        }

        output_paths.passed = filter_temp_dir_path / ("passed_part_" + std::to_string(chunk_index) + ".fastq.gz");
        output_paths.hamming_rejected = filter_temp_dir_path / ("hamming_reject_part_" + std::to_string(chunk_index) + ".fastq.gz");
        output_paths.basecomp_rejected = filter_temp_dir_path / ("basecomp_reject_part_" + std::to_string(chunk_index) + ".fastq.gz");

        ogzstream passed_out(output_paths.passed.c_str());
        ogzstream hamming_out(output_paths.hamming_rejected.c_str());
        ogzstream basecomp_out(output_paths.basecomp_rejected.c_str());

        if (!passed_out.good() || !hamming_out.good() || !basecomp_out.good()) {
            throw std::runtime_error("Failed to open one or more temporary chunk output files for chunk " + std::to_string(chunk_index) + " in " + filter_temp_dir_path.string());
        }

        long long records_processed_in_chunk = 0;
        for (size_t i = 0; i < reads_r1.size(); i += 4) {
            if (i + 3 >= reads_r1.size()) {
                 std::cerr << "Warning: Incomplete record at the end of input vector for chunk " << chunk_index << ". Skipping." << std::endl;
                 break;
            }
            const std::string& header1 = reads_r1[i];
            const std::string& seq1 = reads_r1[i+1];
            const std::string& plus1 = reads_r1[i+2];
            const std::string& qual1 = reads_r1[i+3];

             if (header1.empty() || header1[0] != '@') {
                 std::cerr << "Warning: Skipping malformed record (bad header) in chunk " << chunk_index << ": " << header1.substr(0, std::min((size_t)50, header1.length())) << "..." << std::endl;
                 continue;
             }
              if (plus1.empty() || plus1[0] != '+') {
                  std::cerr << "Warning: Skipping malformed record (bad plus line) in chunk " << chunk_index << " for ID " << header1.substr(0, std::min((size_t)50, header1.length())) << "..." << std::endl;
                 continue;
              }
              if (seq1.length() != qual1.length()) {
                  std::cerr << "Warning: Skipping malformed record (seq/qual length mismatch) in chunk " << chunk_index << " for ID " << header1.substr(0, std::min((size_t)50, header1.length())) << "..." << std::endl;
                 continue;
              }

            ogzstream* target_stream = nullptr;
            bool passed_hamming = false;

            // === Hamming Filter Check ===
            size_t required_hamming_len = static_cast<size_t>(params.start_pos + params.seq_length);
            if (seq1.length() >= required_hamming_len) {
                std::string hamming_sub_seq = seq1.substr(params.start_pos, params.seq_length);
                int dist = hamming_distance(hamming_sub_seq, params.target_seq);

                if (dist > params.threshold) {
                    target_stream = &hamming_out;
                } else {
                    passed_hamming = true;
                }
            } else {
                // Read is too short for Hamming region, filter out
                target_stream = &hamming_out;
            }

            // === Base Composition Filter Check (only if passed Hamming) ===
            if (passed_hamming) {
                size_t required_basecomp_len = static_cast<size_t>(params.base_comp_start_pos + params.base_comp_length);
                if (seq1.length() >= required_basecomp_len) {
                    std::string basecomp_sub_seq = seq1.substr(params.base_comp_start_pos, params.base_comp_length);
                    double base_comp = calculate_base_composition(basecomp_sub_seq);

                    if (base_comp < params.base_composition_threshold) {
                        target_stream = &basecomp_out;
                    } else {
                        // Passed both filters
                        target_stream = &passed_out;
                    }
                } else {
                    // Read is too short for Base Composition region, filter out
                    target_stream = &basecomp_out;
                    // Add a warning here
                    std::cerr << "Warning: Read too short (" << seq1.length() << " bp) for base composition region (start "
                              << params.base_comp_start_pos + 1 << ", length " << params.base_comp_length << ") "
                              << "for ID " << header1.substr(0, std::min((size_t)50, header1.length())) << "... Filtering to basecomp_rejected." << std::endl;
                }
            }

            // Write to the determined stream
            if (target_stream) {
                *target_stream << header1 << "\n" << seq1 << "\n" << plus1 << "\n" << qual1 << "\n";
                if (!target_stream->good()) {
                     throw std::runtime_error("Error writing record (starting with " + header1.substr(0, std::min((size_t)20, header1.length())) + "...) "
                                              "to temporary chunk file for chunk " + std::to_string(chunk_index));
                }
            } else {
                 std::cerr << "Critical Error: No output stream determined for record: " << header1 << std::endl;
                 throw std::runtime_error("Internal logic error: Failed to assign an output stream.");
            }
             records_processed_in_chunk++;
        }

        passed_out.close();
        hamming_out.close();
        basecomp_out.close();

        if (!passed_out.good() || !hamming_out.good() || !basecomp_out.good()) {
             throw std::runtime_error("Error occurred during closing of temporary chunk files for chunk " + std::to_string(chunk_index));
        }
        output_paths.success = true;

    } catch (const std::exception& e) {
         std::cerr << "Error in process_filter_chunk (Index " << chunk_index << "): " << e.what() << std::endl;
         output_paths.success = false;
         // Attempt cleanup of potentially incomplete files
         std::error_code ec;
         fs::remove(output_paths.passed, ec);
         fs::remove(output_paths.hamming_rejected, ec);
         fs::remove(output_paths.basecomp_rejected, ec);
    } catch (...) {
         std::cerr << "Unknown error in process_filter_chunk (Index " << chunk_index << ")" << std::endl;
         output_paths.success = false;
         std::error_code ec;
         fs::remove(output_paths.passed, ec);
         fs::remove(output_paths.hamming_rejected, ec);
         fs::remove(output_paths.basecomp_rejected, ec);
    }
    return output_paths;
}


// ========================================================================
// == Sorting Stage Structures and Functions (Partially Modified)       ==
// ========================================================================

struct FastqRecord {
    std::string id;
    std::string seq;
    std::string plus;
    std::string qual;

    std::string get_base_id() const {
        if (id.empty() || id[0] != '@') return id;
        size_t first_delim = id.find_first_of(" \t/");
        if (first_delim == std::string::npos) {
            return id.substr(1);
        }
        return id.substr(1, first_delim - 1);
    }

    // Ensure consistent comparison for sorting
    bool operator<(const FastqRecord& other) const {
         return get_base_id() < other.get_base_id();
    }
     bool operator>(const FastqRecord& other) const {
        return get_base_id() > other.get_base_id();
    }
};

std::string get_base_id_from_string(const std::string& full_id) {
    if (full_id.empty()) return full_id;
    size_t start_pos = (full_id[0] == '@') ? 1 : 0;
    if (start_pos >= full_id.length()) return "";

    size_t first_delim = full_id.find_first_of(" \t/", start_pos);
    if (first_delim == std::string::npos) {
        return full_id.substr(start_pos);
    }
    return full_id.substr(start_pos, first_delim - start_pos);
}


size_t estimate_memory(const std::vector<FastqRecord>& records) {
    size_t total_mem = 0;
    for (const auto& rec : records) {
        total_mem += rec.id.capacity() + rec.seq.capacity() + rec.plus.capacity() + rec.qual.capacity() + 4; // Use capacity for better estimate
    }
    total_mem += records.capacity() * sizeof(FastqRecord); // Vector overhead
    // Add some buffer
    return total_mem * 1.1;
}

bool write_plain_text_chunk(const std::vector<FastqRecord>& records, const fs::path& temp_filename) {
    std::ofstream outfile(temp_filename, std::ios::binary);
    if (!outfile) {
        std::cerr << "Error: Cannot open temporary plain text file for writing: " << temp_filename << " (" << strerror(errno) << ")" << std::endl;
        return false;
    }
    // Buffering handled by OS/fstream, manually setting might not improve much here for sequential write
    // char buffer[GZ_BUFFER_SIZE];
    // outfile.rdbuf()->pubsetbuf(buffer, sizeof(buffer));

    for (const auto& rec : records) {
        outfile.write(rec.id.c_str(), rec.id.size()); outfile.put('\n');
        outfile.write(rec.seq.c_str(), rec.seq.size()); outfile.put('\n');
        outfile.write(rec.plus.c_str(), rec.plus.size()); outfile.put('\n');
        outfile.write(rec.qual.c_str(), rec.qual.size()); outfile.put('\n');
        if (outfile.fail()) {
             std::cerr << "Error: Failed writing record starting with " << rec.id.substr(0, std::min((size_t)30, rec.id.length())) << " to " << temp_filename << std::endl;
             outfile.close();
             return false;
        }
    }
    outfile.close();
    return !outfile.fail();
}

// PlainChunkReader remains the same, used by the sequential merge phase
struct PlainChunkReader {
    std::ifstream file;
    FastqRecord current_record;
    bool eof = false;
    fs::path filename;

    PlainChunkReader(const fs::path& fname) : filename(fname) {
        file.open(fname, std::ios::binary);
        if (!file) {
            std::cerr << "Error: Cannot open temporary chunk file for reading: " << filename << std::endl;
            eof = true;
        } else {
            read_next();
        }
    }

    ~PlainChunkReader() {
        if (file.is_open()) {
            file.close();
        }
         // Optional: Delete the source file when reader is destroyed if it was successfully opened
         // std::error_code ec;
         // fs::remove(filename, ec); // This is handled by the merge function on success
    }

    PlainChunkReader(const PlainChunkReader&) = delete;
    PlainChunkReader& operator=(const PlainChunkReader&) = delete;
    PlainChunkReader(PlainChunkReader&&) = delete;
    PlainChunkReader& operator=(PlainChunkReader&&) = delete;

    bool read_next() {
        if (eof || !file) return false;
        current_record = {}; // Clear previous record

        if (std::getline(file, current_record.id) &&
            std::getline(file, current_record.seq) &&
            std::getline(file, current_record.plus) &&
            std::getline(file, current_record.qual))
        {
             if (current_record.id.empty() || current_record.id[0] != '@') {
                std::cerr << "Warning: Malformed record ID read from plain chunk " << filename << ": " << current_record.id.substr(0, std::min((size_t)50, current_record.id.length())) << "..." << std::endl;
                // Don't necessarily fail the whole merge, but log it.
             }
             if (current_record.plus.empty() || current_record.plus[0] != '+') {
                 if (!current_record.seq.empty()) {
                     std::cerr << "Warning: Malformed plus line read from plain chunk " << filename << " for ID " << current_record.id.substr(0,std::min((size_t)30, current_record.id.length())) << "..." << std::endl;
                 }
             }
             if (current_record.seq.length() != current_record.qual.length()) {
                  std::cerr << "Warning: Seq/Qual length mismatch read from plain chunk " << filename << " for ID " << current_record.id.substr(0,std::min((size_t)30, current_record.id.length())) << "..." << std::endl;
             }
            return true;
        } else {
            if (file.bad()) {
                 std::cerr << "Error: File stream read error on: " << filename << std::endl;
            } else if (!file.eof() && file.fail()) {
                 if (!current_record.id.empty() || !current_record.seq.empty() || !current_record.plus.empty()) {
                     std::cerr << "Warning: Stream failed reading full record (after ID '" << current_record.id.substr(0,std::min((size_t)30, current_record.id.length())) << "...') before EOF in chunk: " << filename << std::endl;
                 }
            }
            eof = true;
            current_record = {};
            return false;
        }
    }

    bool operator>(const PlainChunkReader& other) const {
        if (eof && !other.eof) return true;  // Non-eof is smaller (comes first)
        if (!eof && other.eof) return false; // This reader is smaller
        if (eof && other.eof) return false;  // Both eof, consider equal for heap stability maybe?
        // Primary comparison based on base ID
        return current_record.get_base_id() > other.current_record.get_base_id();
    }
};

struct PlainChunkReaderPtrCompare {
    bool operator()(const PlainChunkReader* a, const PlainChunkReader* b) const {
        return *a > *b; // Use the > operator defined in PlainChunkReader
    }
};

// Merge function (Phase 2) remains sequential
bool merge_plain_text_chunks_to_gz(const std::vector<fs::path>& plain_chunk_files, const fs::path& final_output_gz_path) {
     if (plain_chunk_files.empty()) {
        std::cout << "No sorted plain text chunks to merge into " << final_output_gz_path << ". Creating empty file." << std::endl;
        fs::path parent_dir = final_output_gz_path.parent_path();
        if (!parent_dir.empty()) {
             std::error_code ec; fs::create_directories(parent_dir, ec);
             if(ec) { std::cerr << "Warning: Failed to create parent directory for empty file " << final_output_gz_path << ": " << ec.message() << std::endl; }
        }
        gzFile empty_out = gzopen(final_output_gz_path.c_str(), "wb");
        if (!empty_out) {
            std::cerr << "Error: Cannot create empty final output file: " << final_output_gz_path << " (" << strerror(errno) << ")" << std::endl;
            return false;
        }
        gzclose(empty_out);
        return true;
    }

    std::cout << "Merging " << plain_chunk_files.size() << " sorted plain text chunk(s) into " << final_output_gz_path << "..." << std::endl;

    fs::path parent_dir = final_output_gz_path.parent_path();
    if (!parent_dir.empty()) {
         std::error_code ec; fs::create_directories(parent_dir, ec);
         if(ec) { throw std::runtime_error("Failed to create output directory: " + parent_dir.string() + " " + ec.message()); }
    }

    std::priority_queue<PlainChunkReader*, std::vector<PlainChunkReader*>, PlainChunkReaderPtrCompare> min_heap;
    std::vector<std::unique_ptr<PlainChunkReader>> readers;
    readers.reserve(plain_chunk_files.size());

    bool any_reader_valid = false;
    for (const auto& fname : plain_chunk_files) {
        std::error_code ec_exists;
        if (!fs::exists(fname, ec_exists)) {
             std::cerr << "Warning: [Merge] Chunk file listed but not found: " << fname << ". Skipping." << std::endl;
             continue;
        }
         std::error_code ec_sz; uintmax_t fsize = fs::file_size(fname, ec_sz);
         if(ec_sz || fsize == 0) {
             if (ec_sz) std::cerr << "Warning: [Merge] Could not get size of " << fname << ", skipping: " << ec_sz.message() << std::endl;
             else std::cerr << "Warning: [Merge] Chunk file " << fname << " is empty. Skipping." << std::endl;
             std::error_code ec_rm; fs::remove(fname, ec_rm); // Clean up empty/inaccessible file
             continue;
        }

        auto reader = std::make_unique<PlainChunkReader>(fname);
        if (!reader->eof) {
            min_heap.push(reader.get());
            any_reader_valid = true;
            readers.push_back(std::move(reader));
        } else {
             // Reader failed to open or read first record
             std::cerr << "Warning: [Merge] Failed to initialize reader or read first record for " << fname << ". Skipping." << std::endl;
             reader.reset(); // Destroy unique_ptr
             std::error_code ec_rm; fs::remove(fname, ec_rm); // Clean up unusable file
        }
    }

    if (!any_reader_valid) {
        std::cout << "Warning: No valid data found in any input chunk files for merge to " << final_output_gz_path << ". Creating empty file." << std::endl;
        readers.clear(); // Ensure all unique_ptrs are destroyed
        // Double check deletion of original inputs
        for (const auto& fname : plain_chunk_files) {
             std::error_code ec_exists;
             if(fs::exists(fname, ec_exists)) {
                std::error_code ec_rm; fs::remove(fname, ec_rm);
             }
        }
        gzFile empty_out = gzopen(final_output_gz_path.c_str(), "wb");
         if (!empty_out) {
            std::cerr << "Error: Cannot create empty final output file: " << final_output_gz_path << " (" << strerror(errno) << ")" << std::endl;
            return false;
        }
        gzclose(empty_out);
        return true;
    }

    gzFile outfile_gz = gzopen(final_output_gz_path.c_str(), "wb");
    if (!outfile_gz) {
        std::cerr << "Error: Cannot open final output file for writing: " << final_output_gz_path << " (" << strerror(errno) << ")" << std::endl;
        // Destructors of PlainChunkReader unique_ptrs will close files
        return false;
    }
    gzbuffer(outfile_gz, GZ_BUFFER_SIZE);

    bool write_error = false;
    long long records_merged = 0;
    const long long MERGE_PROGRESS_INTERVAL = 5000000;
    std::string record_buffer;
    record_buffer.reserve(1024); // Pre-allocate buffer for writing

    while (!min_heap.empty()) {
        PlainChunkReader* top = min_heap.top();
        min_heap.pop();

        // Construct the 4 lines efficiently
        record_buffer.clear();
        record_buffer.append(top->current_record.id); record_buffer.push_back('\n');
        record_buffer.append(top->current_record.seq); record_buffer.push_back('\n');
        record_buffer.append(top->current_record.plus); record_buffer.push_back('\n');
        record_buffer.append(top->current_record.qual); record_buffer.push_back('\n');

        // Write the record
        if (gzwrite(outfile_gz, record_buffer.c_str(), static_cast<unsigned int>(record_buffer.length())) != (int)record_buffer.length())
        {
            int errnum = 0;
            const char *error_str = gzerror(outfile_gz, &errnum);
            std::cerr << "\nError: Failed writing to gzipped output file " << final_output_gz_path
                      << ". Error (" << std::to_string(errnum) << "): " << (error_str ? error_str : "Unknown") << std::endl;
            write_error = true;
            break; // Stop merging on write error
        }

        records_merged++;
        if (records_merged % MERGE_PROGRESS_INTERVAL == 0) {
            std::cout << "." << std::flush;
        }

        // Read the next record from the same chunk file
        if (top->read_next()) {
            min_heap.push(top); // Push back into the heap if more records exist
        } else {
            // Check if reading failed due to an error or just EOF
            if (top->file.bad()) {
                 std::cerr << "\nError reading from chunk file during merge: " << top->filename << std::endl;
                 write_error = true; // Treat read error as critical failure
                 break;
            }
             // If just EOF, the reader is exhausted and won't be pushed back.
             // The unique_ptr will clean it up later.
        }
    }
    std::cout << std::endl; // Newline after progress dots

    // Ensure heap is empty (important if loop broke early due to error)
    while (!min_heap.empty()) min_heap.pop();

    // Close the output file
    int close_status = gzclose(outfile_gz);

    if (close_status != Z_OK) {
        // Z_BUF_ERROR might be acceptable if some data was written, but still log warning
        if (close_status != Z_BUF_ERROR || records_merged == 0) {
             std::cerr << "Error: Failed to properly close gzipped output file " << final_output_gz_path
                       << ". Close status: " << close_status << ". Output may be corrupted." << std::endl;
             write_error = true; // Mark as failure if close fails seriously
        } else {
             std::cerr << "Warning: gzclose returned status " << close_status << " on " << final_output_gz_path
                       << ". File might be usable, but check integrity." << std::endl;
        }
    }

    // Cleanup and final status reporting
    bool final_success = !write_error;
    if (final_success) {
        std::cout << "Merge successful for " << final_output_gz_path << ". Merged " << records_merged << " records." << std::endl;
        // Delete temporary plain text chunks ONLY if merge was successful
        for (const auto& ptr : readers) { // Iterate through the unique_ptrs to get filenames
             const fs::path& fname = ptr->filename;
             std::error_code ec_exists;
             if (fs::exists(fname, ec_exists)) {
               std::error_code ec;
               fs::remove(fname, ec);
               if (ec) {
                    std::cerr << "Warning: Failed to delete temporary plain text chunk: " << fname << " (" << ec.message() << ")" << std::endl;
               }
            }
        }
    } else {
         std::cerr << "Error: Merge failed for " << final_output_gz_path << ". Input plain text chunk files NOT deleted in corresponding temp dir: " << final_output_gz_path.parent_path() << std::endl;
         // Attempt to delete the potentially corrupt/incomplete output file
         std::error_code ec;
         fs::remove(final_output_gz_path, ec);
    }

    readers.clear(); // Explicitly clear the vector of unique_ptrs, triggering destructors
    return final_success;
}


// ========================================================================
// == NEW: Parallel Sorting Phase 1 Worker Function                    ==
// ========================================================================
SortWorkerResult process_sort_sub_chunk(
    const std::vector<fs::path>& sub_input_gz_chunks,
    const fs::path& category_sort_temp_path, // Path to the *category specific* temp dir
    const std::string& category_name,
    size_t sort_memory_limit_bytes,
    bool enable_parallel_sort_alg,
    int worker_id)
{
    SortWorkerResult result;
    result.success = false; // Assume failure initially
    long long worker_total_records = 0;
    int plain_chunk_count = 0;
    std::vector<FastqRecord> current_sort_chunk;
    size_t avg_record_size_bytes = 500; // Initial estimate
    size_t estimated_records_per_chunk = std::max((size_t)10000, sort_memory_limit_bytes / avg_record_size_bytes);
    current_sort_chunk.reserve(estimated_records_per_chunk);
    size_t current_memory_estimate = 0;
    char read_buffer[GZ_BUFFER_SIZE]; // Local read buffer per thread

    // const long long WORKER_PROGRESS_INTERVAL = 2000000; // Feedback from worker

    // std::cout << " [Worker " << worker_id << "] Processing " << sub_input_gz_chunks.size() << " chunks. Mem Lim: " << (sort_memory_limit_bytes / (1024*1024)) << "MB. Est Buf: " << estimated_records_per_chunk << std::endl;

    try {
        for(const auto& gz_chunk_path : sub_input_gz_chunks) {
            // Basic checks before opening
            std::error_code ec_exists;
            if (!fs::exists(gz_chunk_path, ec_exists)) {
                std::cerr << "Warning: [Worker " << worker_id << "] Input chunk file not found, skipping: " << gz_chunk_path << std::endl;
                continue;
            }
             std::error_code ec_sz; uintmax_t fsize = fs::file_size(gz_chunk_path, ec_sz);
             if (ec_sz || fsize == 0) {
                 if (ec_sz) std::cerr << "Warning: [Worker " << worker_id << "] Could not get size of " << gz_chunk_path << ", skipping: " << ec_sz.message() << std::endl;
                 else std::cerr << "Warning: [Worker " << worker_id << "] Input chunk file " << gz_chunk_path << " is empty, skipping." << std::endl;
                 if (!ec_sz && fsize == 0) { std::error_code ec_rm; fs::remove(gz_chunk_path, ec_rm); } // Remove empty source chunk
                 continue;
            }

            gzFile infile = gzopen(gz_chunk_path.c_str(), "rb");
            if (!infile) {
                std::cerr << "Warning: [Worker " << worker_id << "] Cannot open input chunk file: " << gz_chunk_path << " (" << strerror(errno) << "). Skipping." << std::endl;
                continue;
            }
            gzbuffer(infile, GZ_BUFFER_SIZE);

            std::string line_id, line_seq, line_plus, line_qual;
            while (true) {
                char* id_ptr = gzgets(infile, read_buffer, sizeof(read_buffer));
                if (id_ptr == nullptr) {
                    if (gzeof(infile)) break; // Normal EOF for this chunk
                    int err; const char* msg = gzerror(infile, &err);
                    if (err != Z_OK && err != Z_STREAM_END) {
                         std::cerr << "\nError: [Worker " << worker_id << "] Reading ID line from chunk " << gz_chunk_path << ": " << (msg ? msg : "Unknown error") << " (err code: " << err << ")" << std::endl;
                         gzclose(infile);
                         throw std::runtime_error("Worker " + std::to_string(worker_id) + ": Error reading gz chunk file (ID line).");
                    }
                     // If gzgets is null but no major error and not EOF, could be truncated file or oddity
                     if (!gzeof(infile)) {
                         std::cerr << "\nWarning: [Worker " << worker_id << "] gzgets returned NULL but not EOF and no major error on " << gz_chunk_path << ". Attempting to continue." << std::endl;
                         continue; // Try reading next line (might recover or fail definitively)
                     }
                    break; // EOF reached
                }
                // Process line_id (remove newline chars)
                line_id = read_buffer;
                while (!line_id.empty() && (line_id.back() == '\n' || line_id.back() == '\r')) {
                    line_id.pop_back();
                }
                 if (line_id.empty() || line_id[0] != '@') {
                     std::cerr << "\nWarning: [Worker " << worker_id << "] Malformed ID line encountered in " << gz_chunk_path << ". Skipping record. Line: '" << line_id.substr(0, std::min((size_t)100, line_id.length())) << "'" << std::endl;
                     // Skip next 3 lines
                     for(int skip=0; skip<3; ++skip) {
                         if (gzgets(infile, read_buffer, sizeof(read_buffer)) == nullptr) {
                             if (!gzeof(infile)) {
                                  int err; gzerror(infile, &err);
                                  std::cerr << "   (Failed to skip subsequent lines, error code: " << err << ")" << std::endl;
                                  gzclose(infile);
                                  throw std::runtime_error("Worker " + std::to_string(worker_id) + ": Error skipping lines after malformed ID.");
                             }
                              break; // EOF while skipping is ok
                         }
                     }
                     continue; // Go to next potential record start
                 }

                // Helper lambda to read next line and handle errors/EOF
                auto read_next_line = [&](std::string& target_line, const char* line_name) -> bool {
                     if (gzgets(infile, read_buffer, sizeof(read_buffer)) == nullptr) {
                         // Use the original handle_incomplete_record which throws an exception
                         handle_incomplete_record(infile, gz_chunk_path, line_id, line_name);
                         return false; // Should not be reached if throws
                     }
                     target_line = read_buffer;
                     while (!target_line.empty() && (target_line.back() == '\n' || target_line.back() == '\r')) {
                          target_line.pop_back();
                     }
                     return true;
                };

                // Read sequence, plus, quality lines
                if (!read_next_line(line_seq, "sequence")) { break; } // Error handled by lambda (throws)
                if (!read_next_line(line_plus, "plus")) { break; }
                if (!read_next_line(line_qual, "quality")) { break; }

                // Basic validation
                if (line_plus.empty() || line_plus[0] != '+') {
                     std::cerr << "\nWarning: [Worker " << worker_id << "] Malformed plus line encountered in " << gz_chunk_path << " for ID " << line_id.substr(0, std::min((size_t)50, line_id.length())) << ". Skipping record." << std::endl;
                     continue;
                 }
                 if (line_seq.length() != line_qual.length()) {
                      std::cerr << "\nWarning: [Worker " << worker_id << "] Sequence and quality length mismatch in " << gz_chunk_path << " for ID " << line_id.substr(0, std::min((size_t)50, line_id.length())) << ". SeqLen=" << line_seq.length() << ", QualLen=" << line_qual.length() << ". Skipping record." << std::endl;
                     continue;
                 }

                // Add record to buffer
                current_sort_chunk.emplace_back(FastqRecord{std::move(line_id), std::move(line_seq), std::move(line_plus), std::move(line_qual)});
                worker_total_records++;

                // Check if buffer needs to be flushed (based on memory or count)
                 bool memory_check_needed = (current_sort_chunk.size() % (estimated_records_per_chunk / 4 + 1)) == 0;
                 if (current_sort_chunk.size() >= estimated_records_per_chunk || (memory_check_needed && current_sort_chunk.size() > 1000) )
                 {
                    current_memory_estimate = estimate_memory(current_sort_chunk);

                    // Flush condition: Memory limit exceeded OR buffer significantly larger than estimate
                    if (current_memory_estimate >= sort_memory_limit_bytes || current_sort_chunk.size() >= estimated_records_per_chunk * 1.5)
                    {
                        // std::cout << (current_memory_estimate >= sort_memory_limit_bytes ? "M" : "S") << std::flush; // Indicate flush reason

                        // --- Parallel Sort Option ---
                        if (enable_parallel_sort_alg) {
                            #ifdef USE_PARALLEL_SORT
                                try {
                                    // std::cout << "p" << std::flush; // Indicate parallel sort attempt
                                    std::sort(std::execution::par, current_sort_chunk.begin(), current_sort_chunk.end());
                                } catch (const std::exception& e) {
                                    std::cerr << "\nWarning: [Worker " << worker_id << "] Parallel sort failed: " << e.what() << ". Falling back to sequential sort for this chunk." << std::endl;
                                    std::sort(current_sort_chunk.begin(), current_sort_chunk.end());
                                } catch (...) {
                                    std::cerr << "\nWarning: [Worker " << worker_id << "] Parallel sort failed with unknown exception. Falling back to sequential sort for this chunk." << std::endl;
                                     std::sort(current_sort_chunk.begin(), current_sort_chunk.end());
                                }
                            #else
                                // std::cout << "s" << std::flush; // Indicate sequential sort (compile flag off)
                                std::sort(current_sort_chunk.begin(), current_sort_chunk.end());
                            #endif
                        } else {
                            // std::cout << "s" << std::flush; // Indicate sequential sort (disabled by flag)
                            std::sort(current_sort_chunk.begin(), current_sort_chunk.end());
                        }
                        // --- End Parallel Sort Option ---

                        // Write the sorted chunk to a plain text temporary file
                        fs::path temp_filename = category_sort_temp_path /
                            ("sorted_" + category_name + "_worker_" + std::to_string(worker_id) + "_chunk_" + std::to_string(plain_chunk_count++) + ".tmp");

                        if (!write_plain_text_chunk(current_sort_chunk, temp_filename)) {
                            gzclose(infile); // Close input file before throwing
                            throw std::runtime_error("Worker " + std::to_string(worker_id) + ": Error writing sorted plain text chunk file: " + temp_filename.string());
                        }
                        result.created_temp_files.push_back(temp_filename);

                        // Clear the buffer and reset memory estimate
                        current_sort_chunk.clear();
                        // Re-reserve space for efficiency
                        current_sort_chunk.reserve(estimated_records_per_chunk);
                        current_memory_estimate = 0;

                        // Recalculate average record size periodically for better estimates
                        if (plain_chunk_count % 5 == 0 && worker_total_records > 0) {
                             size_t approx_mem_last_chunk = current_memory_estimate > 0 ? current_memory_estimate : sort_memory_limit_bytes; // Use estimate if available
                             avg_record_size_bytes = std::max((size_t)50, approx_mem_last_chunk / std::max((size_t)1, estimated_records_per_chunk)); // Avoid division by zero
                             estimated_records_per_chunk = std::max((size_t)10000, sort_memory_limit_bytes / std::max((size_t)1, avg_record_size_bytes));
                             current_sort_chunk.reserve(estimated_records_per_chunk); // Adjust reservation
                             // std::cout << " [W" << worker_id << " Upd. Est. Recs: " << estimated_records_per_chunk << "] " << std::flush;
                        }
                    }
                }
                // Optional: Worker-level progress indicator
                // if (worker_total_records % WORKER_PROGRESS_INTERVAL == 0) {
                //      std::cout << "w" << worker_id << "." << std::flush;
                // }
            } // End while reading records from current gz chunk

            gzclose(infile); // Close the current input chunk file

            // --- Clean up the processed input gzipped chunk ---
            std::error_code ec_exists_del;
             if (fs::exists(gz_chunk_path, ec_exists_del)) {
                std::error_code ec_rm;
                fs::remove(gz_chunk_path, ec_rm);
                if (ec_rm) {
                     // R2 extracted chunks are special intermediate files, less critical to delete immediately
                     if (category_name == "r2_extracted") {
                          // std::cout << "\nNote: [Worker " << worker_id << "] Failed to remove intermediate unsorted R2 chunk " << gz_chunk_path << ": " << ec_rm.message() << std::endl;
                     } else {
                         std::cerr << "\nWarning: [Worker " << worker_id << "] Failed to remove processed input chunk " << gz_chunk_path << ": " << ec_rm.message() << std::endl;
                     }
                }
             }

        } // End loop over all input gz chunks assigned to this worker

        // Process the final remaining buffer chunk for this worker
        if (!current_sort_chunk.empty()) {
            // std::cout << " [Worker " << worker_id << " Final Buf:" << current_sort_chunk.size() << "] " << std::flush;
            // --- Parallel Sort Option ---
            if (enable_parallel_sort_alg) {
                #ifdef USE_PARALLEL_SORT
                    try {
                        // std::cout << "p" << std::flush;
                        std::sort(std::execution::par, current_sort_chunk.begin(), current_sort_chunk.end());
                    } catch (const std::exception& e) {
                         std::cerr << "\nWarning: [Worker " << worker_id << "] Parallel sort failed (final chunk): " << e.what() << ". Falling back to sequential sort." << std::endl;
                        std::sort(current_sort_chunk.begin(), current_sort_chunk.end());
                    } catch (...) {
                         std::cerr << "\nWarning: [Worker " << worker_id << "] Parallel sort failed (final chunk) with unknown exception. Falling back to sequential sort." << std::endl;
                         std::sort(current_sort_chunk.begin(), current_sort_chunk.end());
                    }
                #else
                    // std::cout << "s" << std::flush;
                    std::sort(current_sort_chunk.begin(), current_sort_chunk.end());
                #endif
            } else {
                 // std::cout << "s" << std::flush;
                 std::sort(current_sort_chunk.begin(), current_sort_chunk.end());
            }
            // --- End Parallel Sort Option ---

            fs::path temp_filename = category_sort_temp_path /
                ("sorted_" + category_name + "_worker_" + std::to_string(worker_id) + "_chunk_" + std::to_string(plain_chunk_count++) + ".tmp");
             if (!write_plain_text_chunk(current_sort_chunk, temp_filename)) {
                 throw std::runtime_error("Worker " + std::to_string(worker_id) + ": Error writing final sorted plain text chunk file: " + temp_filename.string());
             }
            result.created_temp_files.push_back(temp_filename);
            current_sort_chunk.clear(); // Clear buffer
             // std::cout << " Done." << std::endl;
        }

        result.success = true; // Mark worker as successful if no exceptions occurred
        result.records_processed = worker_total_records;

    } catch (const std::exception& e) {
        std::cerr << "\nError in sort worker " << worker_id << " for category '" << category_name << "': " << e.what() << std::endl;
        result.success = false;
        // Cleanup any temp files created by this worker before failure
        for(const auto& temp_f : result.created_temp_files) {
            std::error_code ec; fs::remove(temp_f, ec);
        }
        result.created_temp_files.clear();
        result.records_processed = worker_total_records; // Report records processed *before* error
    } catch (...) {
        std::cerr << "\nUnknown error in sort worker " << worker_id << " for category '" << category_name << "'." << std::endl;
        result.success = false;
        for(const auto& temp_f : result.created_temp_files) {
            std::error_code ec; fs::remove(temp_f, ec);
        }
        result.created_temp_files.clear();
        result.records_processed = worker_total_records;
    }

    return result;
}

// ========================================================================
// == MODIFIED: Main Sorting Function (Orchestrator)                   ==
// ========================================================================
bool sort_and_merge_chunk_category(
    const std::vector<fs::path>& input_gz_chunks,
    const fs::path& final_output_gz_path,
    const fs::path& main_temp_dir_base, // Main processing temp dir (e.g., filter_sort_main_XXXXXX)
    const std::string& category_name,
    size_t sort_memory_limit_bytes_per_worker, // Renamed for clarity
    int num_sort_threads,
    bool enable_parallel_sort_alg) // Flag for std::execution::par
{
    if (input_gz_chunks.empty()) {
        std::cout << "No input chunks for category '" << category_name << "', creating empty final file: " << final_output_gz_path << std::endl;
        fs::path parent_dir = final_output_gz_path.parent_path();
        if (!parent_dir.empty()) {
             std::error_code ec; fs::create_directories(parent_dir, ec);
             if(ec) { std::cerr << "Warning: Failed to create parent directory for empty file " << final_output_gz_path << ": " << ec.message() << std::endl; }
        }
        gzFile empty_out = gzopen(final_output_gz_path.c_str(), "wb");
        if (!empty_out) {
            std::cerr << "Error: Cannot create empty final output file: " << final_output_gz_path << " (" << strerror(errno) << ")" << std::endl;
            return false;
        }
        gzclose(empty_out);
        return true;
    }

    std::cout << "\n--- Starting Sort Phase for category '" << category_name << "' ---" << std::endl;
    std::cout << "Input chunks: " << input_gz_chunks.size()
              << ", Output: " << final_output_gz_path
              << ", Sort Workers: " << num_sort_threads
              << ", Mem Limit/Worker: " << (sort_memory_limit_bytes_per_worker / (1024*1024)) << " MB"
              << ", Parallel Sort Alg: " << (enable_parallel_sort_alg ? "Enabled" : "Disabled")
              #ifdef USE_PARALLEL_SORT
              << " (Compiled with USE_PARALLEL_SORT)"
              #else
              << " (Compiled without USE_PARALLEL_SORT)"
              #endif
              << std::endl;

    // Create a *sub-directory* within the main temp dir for this category's sort files
    fs::path category_sort_temp_path = main_temp_dir_base / ("sort_temp_" + category_name);
     // RAII for the category-specific temp directory
    std::unique_ptr<TempDir> category_temp_dir_ptr;
    try {
       category_temp_dir_ptr = std::make_unique<TempDir>(main_temp_dir_base, "sort_temp_" + category_name + "_");
       category_sort_temp_path = category_temp_dir_ptr->getPath(); // Use the path created by TempDir
       std::cout << "Using temporary directory for sorting intermediate files: " << category_sort_temp_path << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error creating temporary directory for sorting category '" << category_name << "': " << e.what() << std::endl;
        return false;
    }


    long long total_records_in_category = 0;
    std::vector<fs::path> all_sorted_plain_chunk_files;
    std::mutex result_mutex; // Mutex to protect access to all_sorted_plain_chunk_files and total_records_in_category
    std::atomic<int> failed_worker_count = 0;

    // --- Phase 1 (Parallel Sort/Write) ---
    std::cout << "Phase 1 (Parallel Sort/Write): Distributing " << input_gz_chunks.size() << " input chunks to " << num_sort_threads << " workers..." << std::endl;
    std::vector<std::future<SortWorkerResult>> sort_futures;
    std::vector<std::vector<fs::path>> worker_inputs(num_sort_threads);

    // Distribute input files relatively evenly among workers
    for (size_t i = 0; i < input_gz_chunks.size(); ++i) {
        worker_inputs[i % num_sort_threads].push_back(input_gz_chunks[i]);
    }

    // Launch worker threads
    for (int i = 0; i < num_sort_threads; ++i) {
        if (!worker_inputs[i].empty()) {
            sort_futures.push_back(std::async(std::launch::async, process_sort_sub_chunk,
                                               std::move(worker_inputs[i]), // Move vector to worker
                                               category_sort_temp_path,    // Pass path by value/const ref
                                               category_name,
                                               sort_memory_limit_bytes_per_worker,
                                               enable_parallel_sort_alg,
                                               i)); // Pass worker ID
        }
    }

    // --- Collect Results from Workers ---
    std::cout << "Waiting for " << sort_futures.size() << " sort workers to complete..." << std::flush;
    size_t completed_workers = 0;
    for (auto& fut : sort_futures) {
        if (!fut.valid()) {
            std::cerr << "\nWarning: Sort worker future invalid before getting result." << std::endl;
            failed_worker_count++;
            completed_workers++;
            continue;
        }
        try {
            SortWorkerResult worker_result = fut.get(); // Wait for worker to finish
            completed_workers++;
            std::cout << "." << std::flush;

            // Lock before modifying shared results
            std::lock_guard<std::mutex> lock(result_mutex);
            if (worker_result.success) {
                all_sorted_plain_chunk_files.insert(all_sorted_plain_chunk_files.end(),
                                                    std::make_move_iterator(worker_result.created_temp_files.begin()),
                                                    std::make_move_iterator(worker_result.created_temp_files.end()));
                total_records_in_category += worker_result.records_processed;
            } else {
                failed_worker_count++;
                // Temporary files from failed worker should have been cleaned by the worker itself
            }
        } catch (const std::exception& e) {
            std::cerr << "\nError retrieving result from sort worker: " << e.what() << std::endl;
            failed_worker_count++;
            completed_workers++; // Count as completed even if failed
        } catch (...) {
            std::cerr << "\nUnknown error retrieving result from sort worker." << std::endl;
            failed_worker_count++;
            completed_workers++;
        }
    }
    std::cout << " Done." << std::endl;

    if (failed_worker_count > 0) {
         std::cerr << "Error: " << failed_worker_count << " sort worker(s) failed for category '" << category_name << "'." << std::endl;
         // Cleanup: remove successfully created temp files as the overall process failed
         for(const auto& temp_f : all_sorted_plain_chunk_files) {
             std::error_code ec; fs::remove(temp_f, ec);
         }
         // category_temp_dir_ptr RAII will handle directory removal
         return false;
    }

    std::cout << "Phase 1 (Parallel Sort/Write): Finished. Total records processed: " << total_records_in_category
                << ". Created " << all_sorted_plain_chunk_files.size() << " sorted plain text chunk(s)." << std::endl;


    // --- Phase 2 (Sequential Merge) ---
    std::cout << "Phase 2 (Sequential Merge): Merging sorted plain text chunks..." << std::endl;
    if (!merge_plain_text_chunks_to_gz(all_sorted_plain_chunk_files, final_output_gz_path)) {
         std::cerr << "Error: Merging of sorted plain text chunks failed for category '" << category_name << "'." << std::endl;
         // merge_plain_text_chunks_to_gz should NOT delete temp files on failure
         // category_temp_dir_ptr RAII will handle directory removal eventually
         return false;
    }

    // If merge was successful, the merge function should have deleted the .tmp files.
    // The category-specific temp directory will be removed by category_temp_dir_ptr destructor.

    std::cout << "--- Sort Phase for category '" << category_name << "' completed successfully. ---" << std::endl;
    return true;
}


// ========================================================================
// == Helper Functions (Error Handling etc. - Unchanged)                ==
// ========================================================================

void handle_incomplete_record(gzFile infile, const fs::path& chunk_path, const std::string& id_line, const std::string& missing_part) {
    int err = 0;
    const char* gz_err_msg = "Unknown (stream closed)";
    bool is_eof = false;
    if (infile) {
        gz_err_msg = gzerror(infile, &err);
        is_eof = gzeof(infile);
    }

    std::string error_base = "Error: Incomplete FASTQ record in chunk " + chunk_path.string() +
                             " (missing " + missing_part + " line for record starting with '" +
                             id_line.substr(0, std::min((size_t)50, id_line.length())) + "...').";

     if (err != Z_OK && err != Z_STREAM_END) {
          error_base += " GZLIB error: " + std::string(gz_err_msg ? gz_err_msg : "N/A") + " (code " + std::to_string(err) + ")";
     } else if (is_eof) {
          error_base += " Reached EOF unexpectedly.";
     } else {
          error_base += " Read failed but not EOF or GZLIB error."; // e.g., premature stream end?
     }
     std::cerr << "\n" << error_base << std::endl;

    // Throw exception to signal critical error during reading
    throw std::runtime_error("Incomplete FASTQ record detected in gz chunk: " + chunk_path.string() + " ID: " + id_line.substr(0,std::min((size_t)50, id_line.length())));
}

void handle_read_error(igzstream& stream, const fs::path& filepath, long long record_num, const std::string& expected_part) {
    std::string error_msg;
    bool was_eof = stream.eof();
    bool was_bad = stream.bad();
    bool was_fail = stream.fail();
    // Don't close the stream here, let the caller manage it or RAII handle it.

    if (was_eof && was_fail && !was_bad) {
         // Most likely case: Clean EOF happened *during* getline, indicating truncated file
         error_msg = "Incomplete FASTQ record found: File " + filepath.string() +
                     " ended prematurely after reading " + std::to_string(record_num) + " records. Missing " + expected_part + " line.";
    } else if (was_eof) {
         error_msg = "EOF reached unexpectedly while reading " + expected_part + " line from input file " + filepath.string() +
                     " after record " + std::to_string(record_num) + ". Stream state (eof,fail,bad): " +
                      std::to_string(was_eof) + "," + std::to_string(was_fail) + "," + std::to_string(was_bad);
    } else if (was_bad) {
         error_msg = "Stream I/O error (badbit set) reading " + expected_part + " line from input file " + filepath.string() +
                     " near record " + std::to_string(record_num + 1) + ".";
    }
    else { // failbit set, but not eof or badbit - might indicate format error?
         error_msg = "Stream logical error (failbit set) reading " + expected_part + " line from input file " + filepath.string() +
                     " near record " + std::to_string(record_num + 1) + ". Stream state (eof,fail,bad): " +
                     std::to_string(was_eof) + "," + std::to_string(was_fail) + "," + std::to_string(was_bad);
    }
     throw std::runtime_error(error_msg);
}


// ========================================================================
// == R2 Extraction Function (PARALLELIZED - Unchanged from previous)   ==
// ========================================================================
// process_r2_extraction_chunk function remains the same...
fs::path process_r2_extraction_chunk(
    const std::vector<std::string>& r2_chunk_records, // Read-only chunk data
    const std::unordered_set<std::string>& r1_passed_ids_cref, // Read-only access to IDs
    const fs::path& r2_extract_temp_dir, // Directory for this chunk's output
    int chunk_index)                     // Unique index for file naming
{
    fs::path temp_chunk_path = r2_extract_temp_dir / ("r2_extracted_part_" + std::to_string(chunk_index) + ".fastq.gz");
    gzFile out_gz = nullptr;
    long long records_written = 0;

    try {
        out_gz = gzopen(temp_chunk_path.c_str(), "wb");
        if (!out_gz) {
            throw std::runtime_error("Failed to open temporary R2 chunk output file: " + temp_chunk_path.string() + " (" + strerror(errno) + ")");
        }
        gzbuffer(out_gz, GZ_BUFFER_SIZE);

        for (size_t i = 0; i < r2_chunk_records.size(); i += 4) {
             if (i + 3 >= r2_chunk_records.size()) {
                 std::cerr << "Warning: Incomplete record at end of R2 extraction chunk " << chunk_index << ". Skipping." << std::endl;
                 break;
            }
            // Use references directly, avoid copying full lines unless needed
            const std::string& header_line = r2_chunk_records[i];
            const std::string& seq_line = r2_chunk_records[i + 1];
            const std::string& plus_line = r2_chunk_records[i + 2];
            const std::string& qual_line = r2_chunk_records[i + 3];

            // Need to strip newline for ID lookup, but write original line
            std::string header_no_nl = header_line;
             while (!header_no_nl.empty() && (header_no_nl.back() == '\n' || header_no_nl.back() == '\r')) {
                header_no_nl.pop_back();
             }

            if (header_no_nl.empty() || header_no_nl[0] != '@') {
                 std::cerr << "Warning: Skipping malformed R2 record (bad header) in chunk " << chunk_index << ": " << header_no_nl.substr(0, std::min((size_t)50, header_no_nl.length())) << "..." << std::endl;
                 continue;
            }

            std::string base_id = get_base_id_from_string(header_no_nl);
            if (r1_passed_ids_cref.count(base_id)) {
                // Write original lines including newline characters
                if (gzputs(out_gz, header_line.c_str()) == -1 ||
                    gzputs(out_gz, seq_line.c_str()) == -1 ||
                    gzputs(out_gz, plus_line.c_str()) == -1 ||
                    gzputs(out_gz, qual_line.c_str()) == -1)
                {
                    int errnum = 0; const char *error_str = gzerror(out_gz, &errnum);
                    throw std::runtime_error("Error writing extracted R2 record for ID " + header_no_nl.substr(0,std::min((size_t)50, header_no_nl.length())) + " to temp chunk " + temp_chunk_path.string() + ". Error (" + std::to_string(errnum) + "): " + (error_str ? error_str : "Unknown"));
                }
                records_written++;
            }
        } // End loop through records in chunk

        int close_status = gzclose(out_gz);
        out_gz = nullptr; // Mark as closed
        if (close_status != Z_OK) {
             // Don't throw, but log warning. The file might still be partially usable or concat might fail later.
             std::cerr << "Warning: Failed to properly close temporary R2 chunk file " << temp_chunk_path
                       << ". Close status: " << close_status << ". File might be corrupted." << std::endl;
             // Decide if we should delete it here or let concatenation try and fail
             // Let's remove it if close failed badly, otherwise keep it for concatenation attempt
             if(close_status != Z_BUF_ERROR) {
                 std::error_code ec_rm; fs::remove(temp_chunk_path, ec_rm);
                 return fs::path{}; // Return empty path indicating failure
             }
        }

        if (records_written == 0) {
             // No matching records, delete the empty file
             std::error_code ec_rm; fs::remove(temp_chunk_path, ec_rm);
             // std::cout << " [Debug R2 Chunk " << chunk_index << ": 0 records written] " << std::flush;
             return fs::path{}; // Return empty path, not an error, just no output
        }
        // std::cout << " [Debug R2 Chunk " << chunk_index << ": " << records_written << " records written] " << std::flush;
        return temp_chunk_path; // Return path to the successfully written chunk

    } catch (const std::exception& e) {
        std::cerr << "Error in process_r2_extraction_chunk (Index " << chunk_index << "): " << e.what() << std::endl;
        if (out_gz) gzclose(out_gz); // Ensure closed if error occurred after open
        std::error_code ec_rm; fs::remove(temp_chunk_path, ec_rm); // Clean up failed file
        return fs::path{}; // Return empty path indicating failure
    } catch (...) {
        std::cerr << "Unknown error in process_r2_extraction_chunk (Index " << chunk_index << ")" << std::endl;
        if (out_gz) gzclose(out_gz);
        std::error_code ec_rm; fs::remove(temp_chunk_path, ec_rm);
        return fs::path{};
    }
}

// concatenate_gz_files function remains the same...
bool concatenate_gz_files(const std::vector<fs::path>& input_files, const fs::path& output_file) {
    if (input_files.empty()) {
        std::cout << "No input files to concatenate for " << output_file << ". Creating empty file." << std::endl;
         fs::path parent_dir = output_file.parent_path();
         if (!parent_dir.empty()) {
             std::error_code ec; fs::create_directories(parent_dir, ec);
             if(ec) { std::cerr << "Warning: Failed to create parent directory for empty concatenated file " << output_file << ": " << ec.message() << std::endl; }
         }
        gzFile empty_out = gzopen(output_file.c_str(), "wb");
        if (!empty_out) {
            std::cerr << "Error: Cannot create empty concatenated output file: " << output_file << " (" << strerror(errno) << ")" << std::endl;
            return false;
        }
        gzclose(empty_out);
        return true;
    }

    std::cout << "Concatenating " << input_files.size() << " gzipped chunk(s) into " << output_file << "..." << std::flush;

    gzFile out_gz = gzopen(output_file.c_str(), "wb");
    if (!out_gz) {
        std::cerr << "\nError: Cannot open final concatenated output file for writing: " << output_file << " (" << strerror(errno) << ")" << std::endl;
        return false;
    }
    gzbuffer(out_gz, CONCAT_BUFFER_SIZE);

    std::vector<char> buffer(CONCAT_BUFFER_SIZE);
    bool success = true;
    long long total_bytes_written = 0;

    for (size_t i = 0; i < input_files.size(); ++i) {
        const auto& in_path = input_files[i];
        gzFile in_gz = gzopen(in_path.c_str(), "rb");
        if (!in_gz) {
            std::cerr << "\nError: Cannot open input chunk file for concatenation: " << in_path << " (" << strerror(errno) << "). Aborting concatenation." << std::endl;
            success = false;
            break;
        }
        // No need to buffer input stream explicitly with gzread loop
        // gzbuffer(in_gz, CONCAT_BUFFER_SIZE);

        int bytes_read;
        while ((bytes_read = gzread(in_gz, buffer.data(), static_cast<unsigned int>(buffer.size()))) > 0) {
            if (gzwrite(out_gz, buffer.data(), static_cast<unsigned int>(bytes_read)) != bytes_read) {
                int errnum = 0;
                const char* error_str = gzerror(out_gz, &errnum);
                std::cerr << "\nError: Failed writing data during concatenation to " << output_file
                          << ". Error (" << std::to_string(errnum) << "): " << (error_str ? error_str : "Unknown") << ". Aborting." << std::endl;
                success = false;
                break;
            }
            total_bytes_written += bytes_read;
        }

        int read_err = Z_OK;
        if (bytes_read < 0) { // Check for read error
            const char* error_str = gzerror(in_gz, &read_err);
            std::cerr << "\nError: Failed reading data from input chunk " << in_path
                      << " during concatenation. Error (" << read_err << "): " << (error_str ? error_str : "Unknown") << ". Aborting." << std::endl;
            success = false;
        }

        gzclose(in_gz); // Close input chunk file

        if (!success) {
            break; // Exit outer loop if write or read error occurred
        }
        if ((i + 1) % 10 == 0 || (i + 1) == input_files.size()) { // Progress indicator
             std::cout << "." << std::flush;
        }
    }

    int close_status = gzclose(out_gz); // Close output file
    std::cout << " Done." << std::endl;

    if (close_status != Z_OK) {
         std::cerr << "Warning: Failed to properly close concatenated output file " << output_file
                   << ". Close status: " << close_status << ". File might be corrupted." << std::endl;
         if (success && total_bytes_written > 0 && close_status != Z_BUF_ERROR) {
             // If write seemed okay but close failed badly, mark as overall failure
             success = false;
         }
    }

    if (success) {
         std::cout << "Concatenation successful. Total bytes written: " << total_bytes_written << std::endl;
         // Clean up the input chunk files *after* successful concatenation and close
         std::cout << "Cleaning up concatenated chunk files..." << std::flush;
         int delete_failures = 0;
         for(const auto& p : input_files) {
             std::error_code ec;
             if (fs::exists(p)) {
                 if (!fs::remove(p, ec)) {
                     std::cerr << "\nWarning: Failed to remove concatenated chunk: " << p << " (" << ec.message() << ")" << std::endl;
                     delete_failures++;
                 }
             }
         }
         std::cout << (delete_failures == 0 ? " Done." : " Done, with " + std::to_string(delete_failures) + " errors.") << std::endl;
    } else {
         std::cerr << "Error: Concatenation failed for " << output_file << ". Input chunk files NOT deleted." << std::endl;
         // Attempt to remove the failed output file
         std::error_code ec_rm; fs::remove(output_file, ec_rm);
         return false;
    }

    return true;
}

// extract_r2_reads_parallel function remains the same...
fs::path extract_r2_reads_parallel(
    const fs::path& sorted_r1_passed_path,
    const fs::path& original_r2_input_path,
    const fs::path& final_r2_output_path, // Used for naming intermediate file
    const fs::path& main_temp_dir_base,   // Base path for R2 extraction temp chunks
    const Params& params)                // To get num_threads, chunk_read_count
{
    std::cout << "\n--- Starting R2 Extraction Phase (Parallel) ---" << std::endl;
    std::cout << "Using R1 passed IDs from: " << sorted_r1_passed_path << std::endl;
    std::cout << "Scanning original R2 file: " << original_r2_input_path << std::endl;
    std::cout << "Using " << params.num_threads << " threads for R2 extraction." << std::endl;

    fs::path intermediate_r2_path;
    gzFile r1_in_gz = nullptr;
    gzFile r2_in_gz = nullptr;
    std::unique_ptr<TempDir> r2_extract_temp_dir_ptr;

    try {
        if (!fs::exists(sorted_r1_passed_path)) {
            throw std::runtime_error("R1 passed file not found: " + sorted_r1_passed_path.string());
        }
        if (!fs::exists(original_r2_input_path)) {
            throw std::runtime_error("Original R2 input file not found: " + original_r2_input_path.string());
        }
         // Check if R1 passed file is empty
         std::error_code ec_r1_sz; uintmax_t r1_size = fs::file_size(sorted_r1_passed_path, ec_r1_sz);
         if (!ec_r1_sz && r1_size < 50) { // Use a small threshold to account for potential headers in empty gz
             std::cout << "Warning: R1 passed file appears empty or very small (" << r1_size << " bytes). Skipping R2 extraction." << std::endl;
             return fs::path{}; // Return empty path
         } else if (ec_r1_sz) {
              std::cerr << "Warning: Could not determine size of R1 passed file (" << ec_r1_sz.message() << "). Proceeding with R2 extraction, but file might be empty." << std::endl;
         }


        fs::path out_dir = final_r2_output_path.parent_path();
        std::string base_name = final_r2_output_path.filename().string();
         // More robust stem extraction for .fastq.gz
         size_t gz_pos = base_name.rfind(".gz");
         if (gz_pos != std::string::npos) {
             base_name = base_name.substr(0, gz_pos);
         }
         size_t fq_pos = base_name.rfind(".fastq");
         if (fq_pos != std::string::npos) {
             base_name = base_name.substr(0, fq_pos);
         } else {
              fq_pos = base_name.rfind(".fq");
              if (fq_pos != std::string::npos) {
                  base_name = base_name.substr(0, fq_pos);
              }
         }
          if (base_name.empty() && final_r2_output_path.has_stem()) {
             base_name = final_r2_output_path.stem().string();
          } else if (base_name.empty()) {
             base_name = "extracted"; // Fallback base name
          }


        intermediate_r2_path = out_dir / (base_name + "_R2_unsorted_temp.fastq.gz"); // Changed intermediate name slightly
        std::cout << "Final intermediate unsorted R2 file will be: " << intermediate_r2_path << std::endl;

        // --- Create a dedicated temp directory for R2 extraction chunks ---
        r2_extract_temp_dir_ptr = std::make_unique<TempDir>(main_temp_dir_base, "r2_extract_chunks_");
        const fs::path& r2_extract_temp_path = r2_extract_temp_dir_ptr->getPath();
        std::cout << "Using temporary directory for R2 extraction chunks: " << r2_extract_temp_path << std::endl;

        // --- Step 1: Read R1 Passed IDs (Single Threaded) ---
        std::unordered_set<std::string> r1_passed_ids;
        long long r1_id_count = 0;
        const long long R1_ID_READ_PROGRESS = 2000000; // Increased progress interval
        std::cout << "Step 1: Reading R1 passed IDs..." << std::flush;

        r1_in_gz = gzopen(sorted_r1_passed_path.c_str(), "rb");
        if (!r1_in_gz) {
            throw std::runtime_error("Failed to open R1 passed file: " + sorted_r1_passed_path.string() + " (" + strerror(errno) + ")");
        }
        gzbuffer(r1_in_gz, GZ_BUFFER_SIZE);
        char buffer[GZ_BUFFER_SIZE]; // Use a reasonable buffer size
        std::string line;

        while (gzgets(r1_in_gz, buffer, sizeof(buffer)) != nullptr) {
            line = buffer;
            // Efficiently strip trailing newline/carriage return
            size_t len = line.length();
            while (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r')) {
                len--;
            }
            line.resize(len); // Adjust string size

            if (!line.empty() && line[0] == '@') {
                 r1_passed_ids.insert(get_base_id_from_string(line));
                 r1_id_count++;
                 if (r1_id_count % R1_ID_READ_PROGRESS == 0) std::cout << "." << std::flush;

                 // Skip next 3 lines efficiently
                 for (int i = 0; i < 3; ++i) {
                     if (gzgets(r1_in_gz, buffer, sizeof(buffer)) == nullptr) {
                         // Check if it's really EOF or an error
                         if (!gzeof(r1_in_gz)) {
                             int err = 0; gzerror(r1_in_gz, &err);
                             gzclose(r1_in_gz);
                             throw std::runtime_error("Incomplete record found while reading R1 passed IDs from " + sorted_r1_passed_path.string() + " after ID " + line.substr(0, std::min((size_t)50, line.length())) + ". GZ Error: " + std::to_string(err));
                         }
                         // If EOF after reading ID and trying to skip, it's an incomplete record at file end
                         if (i < 2) { // If EOF before skipping all 3 lines
                             gzclose(r1_in_gz);
                             throw std::runtime_error("Incomplete record at end of R1 passed file: " + sorted_r1_passed_path.string() + " after ID " + line.substr(0, std::min((size_t)50, line.length())));
                         }
                         goto end_r1_read; // Break outer loop cleanly if EOF occurs after reading qual
                     }
                 }
                 // Continue to next ID line

            } else if (!line.empty()) {
                 // Log non-@ lines but continue, they might be comments or artifacts
                 std::cerr << "\nWarning: Skipping unexpected line in R1 passed file: " << line.substr(0, std::min((size_t)100, line.length())) << std::endl;
            }
        } // End while gzgets

        end_r1_read:; // Label for goto break
        gzclose(r1_in_gz); r1_in_gz = nullptr;
        std::cout << " Done. Read " << r1_id_count << " IDs. Unique IDs in set: " << r1_passed_ids.size() << "." << std::endl;

        if (r1_id_count == 0 || r1_passed_ids.empty()) {
             std::cout << "No passed R1 IDs found or read. Skipping R2 extraction. No R2 output file will be created." << std::endl;
             // Clean up the potentially created (but empty) R2 temp dir
             r2_extract_temp_dir_ptr.reset();
             return fs::path{}; // Return empty path
        }
        if (r1_id_count != (long long)r1_passed_ids.size()) {
             std::cout << "Note: Number of IDs read (" << r1_id_count << ") differs from final set size (" << r1_passed_ids.size() << "). Duplicates existed in the sorted R1 file." << std::endl;
             // r1_id_count = r1_passed_ids.size(); // Use unique count for consistency - not strictly needed here
        }


        // --- Step 2: Scan R2 in Chunks and Dispatch Tasks ---
        std::cout << "Step 2: Reading R2 input and dispatching extraction tasks..." << std::flush;
        r2_in_gz = gzopen(original_r2_input_path.c_str(), "rb");
         if (!r2_in_gz) {
            throw std::runtime_error("Failed to open original R2 input file: " + original_r2_input_path.string() + " (" + strerror(errno) + ")");
        }
        gzbuffer(r2_in_gz, GZ_BUFFER_SIZE);

        std::vector<std::future<fs::path>> r2_extract_futures;
        std::vector<std::string> current_r2_chunk;
        current_r2_chunk.reserve(params.chunk_read_count * 4);
        long long record_count_in_chunk = 0;
        long long total_r2_records_read = 0;
        int chunk_index = 0;
        bool read_error = false;
        char line_buffer[GZ_BUFFER_SIZE]; // Re-use buffer

        while(true) {
             std::string line1, line2, line3, line4;
             std::string current_id_for_error;

             // Read Line 1 (ID)
             if (gzgets(r2_in_gz, line_buffer, sizeof(line_buffer)) == nullptr) {
                 if(gzeof(r2_in_gz)) break; // Normal EOF
                 int err=0; const char* msg = gzerror(r2_in_gz, &err);
                 std::cerr << "\nError reading ID line from R2 input " << original_r2_input_path.string() << ". GZ Error: " << std::to_string(err) << " " << (msg ? msg : "N/A") << std::endl;
                 read_error = true;
                 break;
            }
            line1 = line_buffer; // Copy buffer content
            current_id_for_error = line1; // Keep original line for error msg
            while (!current_id_for_error.empty() && (current_id_for_error.back() == '\n' || current_id_for_error.back() == '\r')) {
                current_id_for_error.pop_back();
            }


            // Read Line 2 (Seq)
            if (gzgets(r2_in_gz, line_buffer, sizeof(line_buffer)) == nullptr) {
                handle_incomplete_record(r2_in_gz, original_r2_input_path, current_id_for_error, "sequence");
                read_error = true; break;
            }
            line2 = line_buffer;

            // Read Line 3 (Plus)
            if (gzgets(r2_in_gz, line_buffer, sizeof(line_buffer)) == nullptr) {
                handle_incomplete_record(r2_in_gz, original_r2_input_path, current_id_for_error, "plus");
                read_error = true; break;
            }
            line3 = line_buffer;

            // Read Line 4 (Qual)
            if (gzgets(r2_in_gz, line_buffer, sizeof(line_buffer)) == nullptr) {
                 handle_incomplete_record(r2_in_gz, original_r2_input_path, current_id_for_error, "quality");
                 read_error = true; break;
             }
             line4 = line_buffer;

             // Now push all 4 lines (including newlines) into the chunk
             current_r2_chunk.push_back(std::move(line1));
             current_r2_chunk.push_back(std::move(line2));
             current_r2_chunk.push_back(std::move(line3));
             current_r2_chunk.push_back(std::move(line4));

             record_count_in_chunk++;
             total_r2_records_read++;
             if (total_r2_records_read % 5000000 == 0) std::cout << "." << std::flush; // Progress

             // Check if chunk is full
             if (record_count_in_chunk >= params.chunk_read_count) {
                 // Launch async task
                 r2_extract_futures.push_back(std::async(std::launch::async, process_r2_extraction_chunk,
                                             std::move(current_r2_chunk), // Move chunk data
                                             std::cref(r1_passed_ids),     // Pass ID set by const ref
                                             r2_extract_temp_path,       // Pass temp path
                                             chunk_index++));              // Pass unique chunk index

                 // Reset chunk buffer
                 current_r2_chunk = {}; // Clear and deallocate
                 current_r2_chunk.reserve(params.chunk_read_count * 4); // Re-reserve
                 record_count_in_chunk = 0;

                 // --- Throttle Task Submission ---
                 // Avoid overwhelming system with too many pending tasks
                 while (r2_extract_futures.size() >= static_cast<size_t>(params.num_threads * 2 + 4)) { // Heuristic limit
                      bool found_ready = false;
                      for(auto it = r2_extract_futures.begin(); it != r2_extract_futures.end(); ) {
                          // Check if future is ready without blocking
                          if (it->valid() && it->wait_for(std::chrono::milliseconds(0)) == std::future_status::ready) {
                              try { it->get(); } catch (...) { /* Error logged later */ } // Get result to potentially release resources/log error early
                              it = r2_extract_futures.erase(it); // Remove completed/failed future
                              found_ready = true;
                          } else if (!it->valid()) {
                               it = r2_extract_futures.erase(it); // Remove invalid future
                          }
                          else {
                              ++it;
                          }
                      }
                      // If no futures were ready, wait briefly before checking again
                      if (!found_ready && !r2_extract_futures.empty()) {
                            std::this_thread::sleep_for(std::chrono::milliseconds(20)); // Short sleep
                      }
                 } // End throttling loop
             } // End if chunk full
        } // end while reading R2

        gzclose(r2_in_gz); r2_in_gz = nullptr; // Close R2 input file

        // Launch task for the last partial chunk if any records were read and no error occurred
        if (!current_r2_chunk.empty() && !read_error) {
             r2_extract_futures.push_back(std::async(std::launch::async, process_r2_extraction_chunk,
                                         std::move(current_r2_chunk),
                                         std::cref(r1_passed_ids),
                                         r2_extract_temp_path,
                                         chunk_index++));
        }
        std::cout << " Done reading R2." << std::endl;
        std::cout << "Finished reading input. Total R2 records read: " << total_r2_records_read << std::endl;

        if(read_error) {
             // Cleanup futures if read failed mid-way
             for(auto& fut : r2_extract_futures) { if(fut.valid()) { fut.wait(); } } // Wait for running tasks to finish/fail
             throw std::runtime_error("Error occurred during R2 file reading. Aborting extraction.");
        }

        // --- Step 3: Collect R2 Extraction Results ---
        std::vector<fs::path> extracted_r2_chunk_files;
        extracted_r2_chunk_files.reserve(r2_extract_futures.size());
        int failed_chunks = 0;
        // long long total_extracted_r2_records = 0; // Approx count based on file sizes later maybe?

        std::cout << "Step 3: Waiting for " << r2_extract_futures.size() << " R2 extraction tasks to complete..." << std::flush;
        size_t total_futures = r2_extract_futures.size();
        size_t completed_count = 0;
        for (auto& fut : r2_extract_futures) {
             if (!fut.valid()) {
                 std::cerr << "\nWarning: An R2 extraction future was found invalid before getting result." << std::endl;
                 failed_chunks++; completed_count++; continue;
             }
             try {
                 fs::path chunk_path = fut.get(); // Wait for future to complete
                 if (!chunk_path.empty()) {
                     // Verify file exists and is not empty before adding
                     std::error_code ec_exists;
                     if (fs::exists(chunk_path, ec_exists)) {
                        std::error_code ec_sz; uintmax_t fsize = fs::file_size(chunk_path, ec_sz);
                        if (!ec_sz && fsize > 0) {
                             extracted_r2_chunk_files.push_back(chunk_path);
                             // std::cout << "+" << std::flush; // Indicate successful chunk file
                        } else {
                            if (!ec_sz && fsize == 0) {
                                 // std::cout << "0" << std::flush; // Indicate empty chunk file
                                 std::error_code ec_rm; fs::remove(chunk_path, ec_rm); // Clean up empty file
                             } else {
                                 std::cerr << "\nWarning: R2 extraction chunk exists but size check failed or is zero: " << chunk_path << " (Error: " << ec_sz.message() << "). Skipping." << std::endl;
                                 failed_chunks++; // Treat unusable file as failure
                             }
                        }
                     } else {
                        std::cerr << "\nWarning: R2 extraction chunk file reported success but not found: " << chunk_path << ". Worker might have failed cleanup." << std::endl;
                        // Don't increment failed_chunks here, as the worker reported success. Log is sufficient.
                     }
                 } else {
                     // Worker explicitly returned empty path, indicating failure or no records written
                     // std::cout << "-" << std::flush; // Indicate failed/empty chunk task
                     // failed_chunks++; // Don't count "no records written" as a failure unless worker reported error
                 }
             } catch (const std::exception& e) {
                  std::cerr << "\nError processing result from an R2 extraction chunk: " << e.what() << std::endl;
                  failed_chunks++;
             } catch (...) {
                  std::cerr << "\nUnknown error processing result from an R2 extraction chunk." << std::endl;
                  failed_chunks++;
             }
              completed_count++;
              if (completed_count % 10 == 0 || completed_count == total_futures) {
                  std::cout << "." << std::flush; // Progress indicator
              }
        }
        r2_extract_futures.clear(); // Clear the vector of futures
        std::cout << " Done." << std::endl;

        if (failed_chunks > 0) {
             std::cerr << "Warning: " << failed_chunks << " R2 extraction chunk task(s) reported errors or produced unusable files." << std::endl;
             // Decide if this is fatal. Let's try to proceed if *some* chunks succeeded.
             // If you want to abort on *any* failure, throw here.
             // throw std::runtime_error("Aborting R2 extraction due to failures in processing chunks.");
        }
        if (extracted_r2_chunk_files.empty()) {
            std::cout << "No R2 records were extracted (no matching IDs found or all tasks failed/produced empty files)." << std::endl;
            r2_extract_temp_dir_ptr.reset(); // Clean up temp dir
            return fs::path{}; // Return empty path
        }

        std::cout << "Collected " << extracted_r2_chunk_files.size() << " non-empty R2 extracted chunk file(s)." << std::endl;

        // --- Step 4: Concatenate Temporary R2 Chunks ---
        std::cout << "Step 4: Concatenating temporary R2 chunk files..." << std::endl;
        if (!concatenate_gz_files(extracted_r2_chunk_files, intermediate_r2_path)) {
             // Concatenation function handles cleanup of its input files on failure
             throw std::runtime_error("Failed to concatenate temporary R2 chunk files into: " + intermediate_r2_path.string());
        }
        // Concatenation function cleans up its input files on success.

        // --- Step 5: Cleanup R2 Extraction Temp Directory ---
        // The TempDir RAII object will handle this when it goes out of scope.
        // For clarity, we can reset the pointer.
        std::cout << "Step 5: Cleaning up R2 extraction temporary directory..." << std::flush;
        r2_extract_temp_dir_ptr.reset();
        std::cout << " Done." << std::endl;

        std::cout << "Successfully created intermediate unsorted R2 file: " << intermediate_r2_path << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "\nError during parallel R2 extraction: " << e.what() << std::endl;
        // Ensure files are closed if open
        if (r1_in_gz) gzclose(r1_in_gz);
        if (r2_in_gz) gzclose(r2_in_gz);
        // Attempt cleanup of intermediate file if it exists
        std::error_code ec; fs::remove(intermediate_r2_path, ec);
        // Let RAII handle the temp dir pointer
        return fs::path{}; // Return empty path on error
    } catch (...) {
         std::cerr << "\nUnknown error during parallel R2 extraction." << std::endl;
         if (r1_in_gz) gzclose(r1_in_gz);
         if (r2_in_gz) gzclose(r2_in_gz);
         std::error_code ec; fs::remove(intermediate_r2_path, ec);
         return fs::path{};
    }

    std::cout << "--- R2 Extraction Phase (Parallel) Completed Successfully ---" << std::endl;
    return intermediate_r2_path; // Return path to the successfully created intermediate file
}


// ========================================================================
// == Main Function (Adjusted for Sort Threads Option and new params)  ==
// ========================================================================
int main(int argc, char* argv[]) {
    cxxopts::Options options("fastq_filter_sort_extract_sortedR2",
        "Filters R1 FASTQ reads sequentially (Hamming -> Base Comp -> Passed),\n"
        "outputs sorted, gzipped files for each R1 category (using parallel sort/write),\n"
        "extracts corresponding R2 reads PARALLELLY, and sorts the extracted R2 reads (using parallel sort/write).");

    const std::string default_sort_mem_mb = "5120"; // 5 GB default, adjusted from 50GB potentially too large for small machines

    // Determine default thread counts
    unsigned int hardware_threads = std::thread::hardware_concurrency();
    int default_filter_extract_threads = (hardware_threads > 1) ? std::max(1u, hardware_threads / 2) : 1;
    int default_sort_threads = (hardware_threads > 1) ? std::max(1u, hardware_threads / 2) : 1; // Default sort threads same as filter/extract
     // Ensure minimum 1 thread
    if (default_filter_extract_threads == 0) default_filter_extract_threads = 1;
    if (default_sort_threads == 0) default_sort_threads = 1;


    options.add_options()
        ("r1_input", "Input R1 FASTQ.gz file", cxxopts::value<std::string>())
        ("r2_input", "Input R2 FASTQ.gz file (for extraction)", cxxopts::value<std::string>())
        ("r1_output_dir", "Output directory for PASSED & SORTED R1 reads", cxxopts::value<std::string>())
        ("r2_output_dir", "Output directory for EXTRACTED & SORTED R2 reads", cxxopts::value<std::string>())
        ("hamming_filterout_dir", "Output directory for Hamming REJECTED & SORTED R1 reads", cxxopts::value<std::string>())
        ("base_composition_filterout_dir", "Output directory for BaseComp REJECTED & SORTED R1 reads", cxxopts::value<std::string>())
        ("target_seq", "Target sequence for Hamming distance comparison", cxxopts::value<std::string>())
        ("seq_length", "Length of sequence region for Hamming comparison", cxxopts::value<int>())
        ("start_pos", "0-based start position in read for Hamming comparison", cxxopts::value<int>())
        ("threshold", "Maximum allowed Hamming distance (reads > threshold are rejected)", cxxopts::value<int>()->default_value("1"))
        ("base_composition_threshold", "Minimum fraction of ACGT in base comp region (reads < threshold are rejected)", cxxopts::value<double>()->default_value("0.9"))
        ("start_pos_base_com", "0-based start position in read for base composition analysis", cxxopts::value<int>()) // NEW PARAM
        ("seq_length_base_com", "Length of sequence region for base composition analysis", cxxopts::value<int>()) // NEW PARAM
        ("t,threads", "Number of parallel threads for R1 filtering AND R2 extraction phases", cxxopts::value<int>()->default_value(std::to_string(default_filter_extract_threads)))
        ("sort-threads", "Number of parallel threads for the SORTING phase (read/sort/write temp files)", cxxopts::value<int>()->default_value(std::to_string(default_sort_threads)))
        ("chunk_records", "Number of FASTQ records per filtering/extraction chunk", cxxopts::value<long long>()->default_value("1000000"))
        ("T,temp_dir", "Base directory for temporary files (default: value of TMPDIR env var, or current directory '.' if unset)", cxxopts::value<std::string>())
        ("sort-mem-mb", "Approx. memory limit (MB) for *each* sorting worker's buffer (R1 and R2)", cxxopts::value<size_t>()->default_value(default_sort_mem_mb))
        ("parallel-sort-alg", "Enable C++17 parallel sort algorithm (std::execution::par) during in-memory sort step. Requires compiler support and USE_PARALLEL_SORT flag during compilation.", cxxopts::value<bool>()->default_value("true")) // Enable by default
        ("h,help", "Print usage");

    std::unique_ptr<TempDir> main_temp_dir_ptr;
    std::string temp_dir_base_str;
    fs::path intermediate_r2_unsorted_path;

    try {
        auto result = options.parse(argc, argv);

        if (result.count("help")) {
            std::cout << options.help() << std::endl;
             std::cout << "\nNote on Base Composition Filter:\n"
                       << "  Uses the specified region (--start_pos_base_com, --seq_length_base_com) for analysis.\n"
                       << "  Reads shorter than this region will be filtered out to the base composition reject directory.\n";
             std::cout << "Note on Temp Dir: Defaults to $TMPDIR environment variable, or the current directory (.) if $TMPDIR is not set.\n";
             std::cout << "Note on Threads:\n"
                       << "  --threads (-t) controls parallelism for R1 filtering and R2 extraction.\n"
                       << "  --sort-threads controls parallelism for the initial read/sort/write phase of sorting (R1 and R2 categories).\n"
                       << "  The final merge step within sorting is sequential.\n";
             std::cout << "Note on Sort Memory (--sort-mem-mb): This limit applies PER sorting worker thread.\n";
             std::cout << "Note on Parallel Sort Alg (--parallel-sort-alg):\n"
                       << "  Requires C++17 compiler and linking with appropriate threading libraries.\n"
                       << "  Must compile with -D USE_PARALLEL_SORT for the option to have an effect.\n" << std::endl;

            return 0;
        }

        // === Argument Validation ===
        std::vector<std::string> required_args = {
            "r1_input", "r2_input", "r1_output_dir", "r2_output_dir",
            "hamming_filterout_dir", "base_composition_filterout_dir",
            "target_seq", "seq_length", "start_pos", // Hamming args
            "start_pos_base_com", "seq_length_base_com" // BaseComp args (NEW)
        };
        bool missing_arg = false;
        for (const auto& arg : required_args) {
            if (result.count(arg) == 0) {
                std::cerr << "Error: Missing required argument: --" << arg << std::endl;
                missing_arg = true;
            }
        }
        if (missing_arg) {
             std::cerr << "\nRun with --help for usage information." << std::endl;
             return 1;
        }

        Params params;
        params.r1_input = result["r1_input"].as<std::string>();
        params.r2_input = result["r2_input"].as<std::string>();
        params.r1_output_dir = result["r1_output_dir"].as<std::string>();
        params.r2_output_dir = result["r2_output_dir"].as<std::string>();
        params.hamming_filterout_dir = result["hamming_filterout_dir"].as<std::string>();
        params.base_composition_filterout_dir = result["base_composition_filterout_dir"].as<std::string>();
        params.target_seq = result["target_seq"].as<std::string>();
        params.seq_length = result["seq_length"].as<int>();
        params.start_pos = result["start_pos"].as<int>();
        params.threshold = result["threshold"].as<int>();
        params.base_composition_threshold = result["base_composition_threshold"].as<double>();
        // Get new base comp params
        params.base_comp_start_pos = result["start_pos_base_com"].as<int>();
        params.base_comp_length = result["seq_length_base_com"].as<int>();

        params.num_threads = result["threads"].as<int>();
        params.sort_threads = result["sort-threads"].as<int>();
        params.chunk_read_count = result["chunk_records"].as<long long>();
        params.sort_memory_limit_mb = result["sort-mem-mb"].as<size_t>();
        bool enable_parallel_sort_alg = result["parallel-sort-alg"].as<bool>();


        // === Determine temporary directory base path ===
        if (result.count("temp_dir")) {
             temp_dir_base_str = result["temp_dir"].as<std::string>();
             std::cout << "User specified temp directory: " << temp_dir_base_str << std::endl;
        } else {
             const char* tmpdir_env = std::getenv("TMPDIR");
             if (tmpdir_env != nullptr && std::strlen(tmpdir_env) > 0) {
                 temp_dir_base_str = tmpdir_env;
                 std::cout << "Using temp directory from $TMPDIR: " << temp_dir_base_str << std::endl;
             } else {
                 temp_dir_base_str = ".";
                 std::cout << "Using default temp directory: current working directory (.)" << std::endl;
             }
        }
        try {
             params.main_temp_dir_base = fs::absolute(temp_dir_base_str);
        } catch (const fs::filesystem_error& e) {
            std::cerr << "Error resolving temporary directory path '" << temp_dir_base_str << "': " << e.what() << std::endl;
            return 1;
        }
        std::cout << "Effective absolute base temp path: " << params.main_temp_dir_base << std::endl;

        // === Parameter Value Validation ===
        if (!fs::exists(params.r1_input)) { throw std::runtime_error("Input R1 file not found: " + params.r1_input.string()); }
        if (!fs::is_regular_file(params.r1_input)) { throw std::runtime_error("Input R1 path is not a regular file: " + params.r1_input.string()); }
        if (!fs::exists(params.r2_input)) { throw std::runtime_error("Input R2 file not found: " + params.r2_input.string()); }
        if (!fs::is_regular_file(params.r2_input)) { throw std::runtime_error("Input R2 path is not a regular file: " + params.r2_input.string()); }

        // Hamming validation
        if (params.target_seq.length() != static_cast<size_t>(params.seq_length)) {
            throw std::runtime_error("target_seq length (" + std::to_string(params.target_seq.length())
                                     + ") must match seq_length (" + std::to_string(params.seq_length) + ").");
        }
         if (params.start_pos < 0) { throw std::runtime_error("start_pos cannot be negative for Hamming filter."); }
         if (params.seq_length <= 0) { throw std::runtime_error("seq_length must be positive for Hamming filter."); }
         if (params.threshold < 0) { throw std::runtime_error("threshold cannot be negative for Hamming filter."); }

        // Base Comp validation
         if (params.base_comp_start_pos < 0) { throw std::runtime_error("start_pos_base_com cannot be negative."); }
         if (params.base_comp_length <= 0) { throw std::runtime_error("seq_length_base_com must be positive."); }
        if (params.base_composition_threshold < 0.0 || params.base_composition_threshold > 1.0) { throw std::runtime_error("base_composition_threshold must be between 0.0 and 1.0."); }


        if (params.num_threads <= 0) { std::cerr << "Warning: num_threads must be positive. Setting to 1." << std::endl; params.num_threads = 1; }
        if (params.sort_threads <= 0) { std::cerr << "Warning: sort_threads must be positive. Setting to 1." << std::endl; params.sort_threads = 1; }
        if (params.chunk_read_count <= 0) { std::cerr << "Warning: chunk_records must be positive. Setting to 1,000,000." << std::endl; params.chunk_read_count = 1000000; }
        if (params.sort_memory_limit_mb == 0) { std::cerr << "Warning: --sort-mem-mb cannot be 0. Setting to default 1024 MB." << std::endl; params.sort_memory_limit_mb = 1024; }


        #ifndef USE_PARALLEL_SORT
        if (enable_parallel_sort_alg) {
             std::cout << "Warning: --parallel-sort-alg=true specified, but code was not compiled with -D USE_PARALLEL_SORT. Parallel algorithm will NOT be used." << std::endl;
             enable_parallel_sort_alg = false; // Override if not compiled correctly
        }
        #endif


        // === Setup ==
        std::cout << "--- Configuration ---" << std::endl;
        std::cout << "Input R1: " << params.r1_input << std::endl;
        std::cout << "Input R2: " << params.r2_input << std::endl;
        std::cout << "Output Passed R1 Dir: " << params.r1_output_dir << std::endl;
        std::cout << "Output Sorted R2 Dir: " << params.r2_output_dir << std::endl;
        std::cout << "Output Hamming Reject Dir: " << params.hamming_filterout_dir << std::endl;
        std::cout << "Output BaseComp Reject Dir: " << params.base_composition_filterout_dir << std::endl;
        std::cout << "Hamming Filter Region (0-based): Pos " << params.start_pos << "-" << (params.start_pos + params.seq_length -1 )
                  << " (Length: " << params.seq_length << ")" << std::endl;
        std::cout << "Hamming Target Seq: " << params.target_seq << ", Max Distance Threshold: " << params.threshold << std::endl;
        std::cout << "Base Comp Filter Region (0-based): Pos " << params.base_comp_start_pos << "-" << (params.base_comp_start_pos + params.base_comp_length - 1)
                  << " (Length: " << params.base_comp_length << ")" << std::endl;
        std::cout << "Base Comp Min ACGT Threshold: " << params.base_composition_threshold << std::endl;
        std::cout << "Threads (Filtering/Extraction): " << params.num_threads << std::endl;
        std::cout << "Threads (Sorting Phase 1): " << params.sort_threads << std::endl;
        std::cout << "Chunk Records (Filtering & Extraction): " << params.chunk_read_count << std::endl;
        std::cout << "Sorting Memory Limit (per Worker): " << params.sort_memory_limit_mb << " MB" << std::endl;
         std::cout << "Use Parallel Sort Algorithm: " << (enable_parallel_sort_alg ? "Yes" : "No") << std::endl;
        std::cout << "Base Temp Directory (Absolute): " << params.main_temp_dir_base << std::endl; // Show absolute path used
        std::cout << "--------------------" << std::endl;

        // Create main temp dir *after* validating the base path and converting to absolute
        main_temp_dir_ptr = std::make_unique<TempDir>(params.main_temp_dir_base, "filter_sort_main_");
        const fs::path& main_processing_temp_path = main_temp_dir_ptr->getPath(); // Base path for filter chunks and sort/extract subdirs
        std::cout << "Using main temporary directory: " << main_processing_temp_path << std::endl;

        // --- Filtering Phase ---
        std::cout << "\n--- Starting R1 Filtering Phase ---" << std::endl;
        std::vector<std::future<UnsortedChunkOutputPaths>> filter_futures;

        igzstream in_r1(params.r1_input.c_str());
        if (!in_r1.good()) {
            throw std::runtime_error("Failed to open input R1 file: " + params.r1_input.string());
        }

        std::vector<std::string> current_chunk_r1;
        const size_t lines_per_record = 4;
        current_chunk_r1.reserve(params.chunk_read_count * lines_per_record);
        std::string line1, line2, line3, line4; // Reuse strings
        long long record_count_in_chunk = 0;
        long long total_records_read = 0;
        int chunk_index = 0;
        auto start_time_filter = std::chrono::high_resolution_clock::now();

        std::cout << "Reading R1 input and dispatching filtering tasks..." << std::flush;

        while (std::getline(in_r1, line1)) {
             if (!std::getline(in_r1, line2)) { handle_read_error(in_r1, params.r1_input, total_records_read, "sequence"); break; }
             if (!std::getline(in_r1, line3)) { handle_read_error(in_r1, params.r1_input, total_records_read, "plus"); break; }
             if (!std::getline(in_r1, line4)) { handle_read_error(in_r1, params.r1_input, total_records_read, "quality"); break; }

             current_chunk_r1.push_back(std::move(line1)); // Move lines into vector
             current_chunk_r1.push_back(std::move(line2));
             current_chunk_r1.push_back(std::move(line3));
             current_chunk_r1.push_back(std::move(line4));

             record_count_in_chunk++;
             total_records_read++;
             if (total_records_read % 2000000 == 0) std::cout << "." << std::flush; // Adjust progress frequency

             if (record_count_in_chunk >= params.chunk_read_count) {
                 filter_futures.push_back(std::async(std::launch::async, process_filter_chunk,
                                           std::move(current_chunk_r1), // Move vector to task
                                           main_processing_temp_path,   // Pass the main temp path
                                           chunk_index++,
                                           std::cref(params)));          // Pass params by const ref

                 current_chunk_r1 = {}; // Clear and deallocate moved-from vector
                 current_chunk_r1.reserve(params.chunk_read_count * lines_per_record); // Re-reserve
                 record_count_in_chunk = 0;

                 // Throttle task submission slightly differently for filtering
                 while (filter_futures.size() >= static_cast<size_t>(params.num_threads * 2 + 4)) { // Similar heuristic
                      bool found_ready = false;
                      for(auto it = filter_futures.begin(); it != filter_futures.end(); ) {
                          if (it->valid() && it->wait_for(std::chrono::milliseconds(0)) == std::future_status::ready) {
                              try {
                                   it->get(); // Get result, check for exceptions implicitly
                              } catch (const std::exception& e) {
                                   std::cerr << "\nAsync Error (throttling check, filter): " << e.what() << std::endl;
                                   // Decide if this error is critical - let's assume it's handled later when collecting results
                              }
                              it = filter_futures.erase(it);
                              found_ready = true;
                          } else if (!it->valid()) {
                               it = filter_futures.erase(it); // Remove invalid future
                          } else {
                              ++it;
                          }
                      }
                      if (!found_ready && !filter_futures.empty()) {
                           std::this_thread::sleep_for(std::chrono::milliseconds(20)); // Short sleep
                      }
                 } // End throttling loop
             }
        }
        // Check stream state after loop
        if (!in_r1.eof() && in_r1.bad()) {
             in_r1.close(); // Close stream before throwing
             throw std::runtime_error("I/O error reading R1 input file: " + params.r1_input.string());
        }
        in_r1.close(); // Close R1 input file

        // Process the last partial chunk
        if (!current_chunk_r1.empty()) {
             filter_futures.push_back(std::async(std::launch::async, process_filter_chunk,
                                       std::move(current_chunk_r1),
                                       main_processing_temp_path,
                                       chunk_index++,
                                       std::cref(params)));
        }
        std::cout << " Done reading R1." << std::endl;
        std::cout << "Finished reading input. Total R1 records read: " << total_records_read << std::endl;

        // --- Collect Filtering Results ---
        std::map<std::string, std::vector<fs::path>> unsorted_chunks_map;
        unsorted_chunks_map["passed"] = {};
        unsorted_chunks_map["hamming_rejected"] = {};
        unsorted_chunks_map["basecomp_rejected"] = {};
        int failed_filter_chunks = 0;

        std::cout << "Waiting for " << filter_futures.size() << " filtering tasks to complete..." << std::flush;
        size_t total_filter_futures = filter_futures.size();
        size_t completed_filter_count = 0;
        for (auto& fut : filter_futures) {
             if (!fut.valid()) {
                 std::cerr << "\nWarning: A filtering future was found invalid before getting result." << std::endl;
                 failed_filter_chunks++; completed_filter_count++; continue;
             }
             try {
                 UnsortedChunkOutputPaths chunk_paths = fut.get(); // Wait for task
                 completed_filter_count++;
                 std::cout << "." << std::flush;

                 if (chunk_paths.success) {
                     // Helper to add valid, non-empty chunk files
                     auto add_if_valid = [&](const fs::path& p, const std::string& cat) {
                         std::error_code ec_exists;
                         if (!fs::exists(p, ec_exists)) return; // Skip if doesn't exist
                         std::error_code ec_sz; uintmax_t fsize = fs::file_size(p, ec_sz);
                         if (!ec_sz && fsize > 0) {
                             unsorted_chunks_map[cat].push_back(p);
                         } else if (!ec_sz && fsize == 0) {
                             std::error_code ec_rm; fs::remove(p, ec_rm); // Remove empty file silently
                         } else if (ec_sz) {
                             std::cerr << "\nWarning: Could not get size of filter chunk " << p << ": " << ec_sz.message() << ". Skipping." << std::endl;
                         }
                     };
                     add_if_valid(chunk_paths.passed, "passed");
                     add_if_valid(chunk_paths.hamming_rejected, "hamming_rejected");
                     add_if_valid(chunk_paths.basecomp_rejected, "basecomp_rejected");
                 } else {
                     failed_filter_chunks++;
                 }
             } catch (const std::exception& e) {
                  std::cerr << "\nError processing result from a filtering chunk: " << e.what() << std::endl;
                  failed_filter_chunks++; completed_filter_count++; // Ensure count advances
             } catch (...) {
                  std::cerr << "\nUnknown error processing result from a filtering chunk." << std::endl;
                  failed_filter_chunks++; completed_filter_count++;
             }
        }
        filter_futures.clear(); // Clear futures

        auto end_time_filter = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> filter_duration = end_time_filter - start_time_filter;
        std::cout << " Done." << std::endl;
        std::cout << "All filtering tasks finished in " << std::fixed << std::setprecision(2) << filter_duration.count() << " seconds." << std::endl;
        if (failed_filter_chunks > 0) {
              std::cerr << "Warning: " << failed_filter_chunks << " filtering chunk(s) failed processing. Output for affected categories might be incomplete." << std::endl;
              // Decide if this is fatal. Let's proceed for now.
        }
        std::cout << "--- R1 Filtering Phase Completed ---" << std::endl;

        // --- Prepare for Sorting ---
        std::cout << "\n--- Starting R1 Sorting and Merging Phase ---" << std::endl;
        auto start_time_sort = std::chrono::high_resolution_clock::now();

        // Extract base filename more robustly
        std::string base_filename = params.r1_input.filename().string();
        size_t gz_pos = base_filename.rfind(".gz");
        if (gz_pos != std::string::npos) {
            base_filename = base_filename.substr(0, gz_pos);
        }
        size_t fq_pos = base_filename.rfind(".fastq");
        if (fq_pos != std::string::npos) {
            base_filename = base_filename.substr(0, fq_pos);
        } else {
             fq_pos = base_filename.rfind(".fq");
             if (fq_pos != std::string::npos) {
                 base_filename = base_filename.substr(0, fq_pos);
             }
        }
        // Handle cases like just "input.gz" -> "input"
         if (base_filename.empty() && params.r1_input.has_stem()) {
            base_filename = params.r1_input.stem().string();
         } else if (base_filename.empty()) {
            base_filename = "output"; // Fallback base name
         }


        // Create output directories
        std::error_code ec_dir;
        auto create_dir = [&](const fs::path& dir) {
            fs::create_directories(dir, ec_dir);
            if(ec_dir) throw fs::filesystem_error("Cannot create directory", dir, ec_dir);
        };
        create_dir(params.r1_output_dir);
        create_dir(params.hamming_filterout_dir);
        create_dir(params.base_composition_filterout_dir);
        create_dir(params.r2_output_dir);

        // Define final output paths
        std::map<std::string, fs::path> final_output_map = {
            {"passed", params.r1_output_dir / (base_filename + "_R1_passed_sorted.fastq.gz")},
            {"hamming_rejected", params.hamming_filterout_dir / (base_filename + "_R1_hamming_rejected_sorted.fastq.gz")},
            {"basecomp_rejected", params.base_composition_filterout_dir / (base_filename + "_R1_basecomp_rejected_sorted.fastq.gz")}
        };
        std::vector<std::string> r1_categories_to_process = {"passed", "hamming_rejected", "basecomp_rejected"};

        bool overall_r1_success = true;
        size_t sort_mem_bytes_per_worker = params.sort_memory_limit_mb * 1024 * 1024;
        bool passed_r1_category_succeeded = false;

        // Process each R1 category (now using parallel sort internally)
        for (const auto& category : r1_categories_to_process) {
            const auto& input_chunks = unsorted_chunks_map[category];
            const auto& final_output_file = final_output_map.at(category);

            // Call the orchestrator function with sort thread count and parallel alg flag
            if (!sort_and_merge_chunk_category(input_chunks,
                                               final_output_file,
                                               main_processing_temp_path, // Pass main temp path
                                               category,
                                               sort_mem_bytes_per_worker,
                                               params.sort_threads, // Use dedicated sort threads
                                               enable_parallel_sort_alg)) // Pass the flag
             {
                  std::cerr << "Error: Sorting and merging failed for R1 category '" << category << "'." << std::endl;
                  overall_r1_success = false;
             }
             else if (category == "passed") {
                 // Check if the output file actually exists and has size > 0 after reported success
                 std::error_code ec_out_ex, ec_out_sz;
                 if (fs::exists(final_output_file, ec_out_ex) && !ec_out_ex) {
                     uintmax_t out_size = fs::file_size(final_output_file, ec_out_sz);
                     // Consider success only if file exists and size calculation didn't error OR if size is > ~0 (for empty valid gz)
                     if (!ec_out_sz || out_size > 28) { // 28 is approx size of empty gz file header/footer
                        passed_r1_category_succeeded = true;
                     } else {
                          std::cerr << "Warning: R1 passed category sorting reported success, but output file '" << final_output_file << "' is missing, empty, or size check failed (Size: " << (ec_out_sz ? -1 : (long long)out_size) << "). Treating as failure for R2 extraction." << std::endl;
                          overall_r1_success = false; // Mark overall as failed if crucial output is bad
                          passed_r1_category_succeeded = false;
                     }
                 } else {
                      std::cerr << "Warning: R1 passed category sorting reported success, but output file '" << final_output_file << "' does not exist. Treating as failure for R2 extraction." << std::endl;
                      overall_r1_success = false;
                      passed_r1_category_succeeded = false;
                 }
             }
        }

        auto end_time_sort = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> sort_duration = end_time_sort - start_time_sort;
         std::cout << "--- R1 Sorting and Merging Phase Completed in " << std::fixed << std::setprecision(2) << sort_duration.count() << " seconds. ---" << std::endl;

        // --- R2 Extraction and Sorting Steps ---
        bool r2_extraction_succeeded = false;
        bool r2_sorting_succeeded = false;
        fs::path final_r2_output_path; // Define outside the conditional block

        if (passed_r1_category_succeeded) {
             std::cout << "\n--- Proceeding to R2 Extraction & Sorting ---" << std::endl;
             // Construct R2 output path using the same base filename
             final_r2_output_path = params.r2_output_dir / (base_filename + "_R2_extracted_sorted.fastq.gz");

             auto start_time_extract_r2 = std::chrono::high_resolution_clock::now();
             intermediate_r2_unsorted_path = extract_r2_reads_parallel(
                 final_output_map["passed"],   // Input: Sorted R1 passed file
                 params.r2_input,              // Input: Original R2 file
                 final_r2_output_path,         // Used for naming intermediate file
                 main_processing_temp_path,    // Base temp dir for R2 extraction sub-temp-dir
                 params                        // Pass all params (contains num_threads for extraction)
             );
             auto end_time_extract_r2 = std::chrono::high_resolution_clock::now();
             std::chrono::duration<double> extract_r2_duration = end_time_extract_r2 - start_time_extract_r2;


             if (!intermediate_r2_unsorted_path.empty()) {
                 r2_extraction_succeeded = true;
                 std::cout << "R2 Extraction successful in " << std::fixed << std::setprecision(2) << extract_r2_duration.count() << " seconds." << std::endl;
                 std::cout << "Intermediate file: " << intermediate_r2_unsorted_path << std::endl;

                 std::cout << "\n--- Starting R2 Sorting Phase ---" << std::endl;
                 auto start_time_sort_r2 = std::chrono::high_resolution_clock::now();

                 // Input for R2 sorting is the single concatenated intermediate file
                 std::vector<fs::path> r2_input_vector = {intermediate_r2_unsorted_path};

                 // Call the same parallelized sort function for the R2 intermediate file
                 r2_sorting_succeeded = sort_and_merge_chunk_category(
                     r2_input_vector,             // Vector containing the single intermediate R2 path
                     final_r2_output_path,        // Final sorted R2 output path
                     main_processing_temp_path,   // Base temp dir for sort sub-temp-dir
                     "r2_extracted",              // Category name for logging/temp files
                     sort_mem_bytes_per_worker,   // Memory limit per worker
                     params.sort_threads,         // Use sort_threads for R2 sorting
                     enable_parallel_sort_alg     // Use parallel sort algorithm if enabled
                 );

                 auto end_time_sort_r2 = std::chrono::high_resolution_clock::now();
                 std::chrono::duration<double> sort_r2_duration = end_time_sort_r2 - start_time_sort_r2;

                 if (r2_sorting_succeeded) {
                     std::cout << "--- R2 Sorting Phase Completed Successfully in " << std::fixed << std::setprecision(2) << sort_r2_duration.count() << " seconds. ---" << std::endl;
                     // sort_and_merge_chunk_category should handle deleting its input (intermediate_r2_unsorted_path) on success.
                     // Double check if it still exists (it shouldn't)
                     std::error_code ec_exists_inter;
                     if (fs::exists(intermediate_r2_unsorted_path, ec_exists_inter)) {
                        std::cerr << "Warning: Intermediate unsorted R2 file still exists after successful sort: " << intermediate_r2_unsorted_path << ". Attempting removal." << std::endl;
                        std::error_code ec_rm_inter; fs::remove(intermediate_r2_unsorted_path, ec_rm_inter);
                     }
                 } else {
                     std::cerr << "Error: R2 Sorting failed after " << std::fixed << std::setprecision(2) << sort_r2_duration.count() << " seconds." << std::endl;
                     // If R2 sort fails, the intermediate file might still exist and be useful
                     if (fs::exists(intermediate_r2_unsorted_path)) {
                        std::cerr << "Intermediate unsorted R2 file (input to failed sort) kept at: " << intermediate_r2_unsorted_path << std::endl;
                     }
                     // Delete the potentially empty/corrupt final R2 output file
                     std::error_code ec_rm_final; fs::remove(final_r2_output_path, ec_rm_final);
                 }

             } else {
                  std::cerr << "Error: R2 extraction failed or produced no reads after "
                            << std::fixed << std::setprecision(2) << extract_r2_duration.count() << " seconds. Skipping R2 sorting." << std::endl;
                  r2_extraction_succeeded = false; // Mark as failed
                  // Delete the potentially empty final R2 output file path target
                  std::error_code ec_rm_final; fs::remove(final_r2_output_path, ec_rm_final);
             }
        } else {
             std::cout << "\n--- Skipping R2 Extraction & Sorting ---" << std::endl;
             std::cout << "Reason: R1 'passed' category processing did not complete successfully or output file was invalid." << std::endl;
              if (failed_filter_chunks > 0 && unsorted_chunks_map["passed"].empty()) {
                  std::cout << "  (No passed R1 chunks were generated due to filtering chunk failures)." << std::endl;
              } else if (!passed_r1_category_succeeded) {
                   std::cout << "  (R1 passed category sorting/merging reported an error or produced an invalid output file)." << std::endl;
              }
             // Define the path just for the final report message
             final_r2_output_path = params.r2_output_dir / (base_filename + "_R2_extracted_sorted.fastq.gz");
        }


        // --- Final Report ---
        std::cout << "\n--- Final Summary ---" << std::endl;
        // Define overall success. Must have R1 success. If R1 passed category was attempted, R2 must also succeed.
        bool final_status_ok = overall_r1_success;
        if (overall_r1_success && passed_r1_category_succeeded) {
             // If R1 passed was expected to succeed, R2 steps must also succeed
             final_status_ok = r2_extraction_succeeded && r2_sorting_succeeded;
        } else if (overall_r1_success && !passed_r1_category_succeeded) {
            // If R1 passed failed (e.g. no passed reads), but other R1 categories were ok, consider it a success *without* R2.
            final_status_ok = true;
        }


        // Report specific failures
        if (failed_filter_chunks > 0) {
             std::cout << "Processing finished with " << failed_filter_chunks << " failed R1 filtering chunk(s)." << std::endl;
        }
         if (!overall_r1_success) { // Catch-all for R1 sort/merge failures
             std::cout << "Processing finished with errors during R1 sorting/merging phase(s)." << std::endl;
        }
        if (passed_r1_category_succeeded) { // Only report R2 issues if R1 passed seemed ok initially
             if (!r2_extraction_succeeded) {
                 std::cout << "Processing finished with errors during R2 extraction." << std::endl;
             } else if (!r2_sorting_succeeded) { // Only report R2 sort failure if extraction succeeded
                 std::cout << "Processing finished with errors during R2 sorting." << std::endl;
             }
        }


        if (final_status_ok) {
             std::cout << "Processing finished successfully!" << std::endl;
        } else {
             std::cout << "Processing finished with errors or warnings. Please check logs and output files." << std::endl;
        }

        // Report final file status
        std::cout << "Final Output Files:" << std::endl;
        auto report_file_status = [](const std::string& label, const fs::path& path, bool expected_to_exist) {
             std::error_code ec_exists;
             bool exists = fs::exists(path, ec_exists);
             if (exists && !ec_exists) {
                  std::error_code ec_sz; uintmax_t fsize = fs::file_size(path, ec_sz);
                  if (!ec_sz) {
                     std::cout << "  " << label << ": " << path << " (Size: " << std::fixed << std::setprecision(2) << fsize / (1024*1024.0) << " MB)" << std::endl;
                  } else {
                      std::cout << "  " << label << ": " << path << " (Exists, but size check failed: " << ec_sz.message() << ")" << std::endl;
                  }
             } else {
                 if (expected_to_exist) {
                     std::cout << "  " << label << ": " << path << " (ERROR: File not created or empty)" << std::endl;
                 } else {
                     std::cout << "  " << label << ": " << path << " (Not created or empty, as expected or due to upstream errors)" << std::endl;
                 }
             }
        };

        report_file_status("R1 Passed", final_output_map["passed"], overall_r1_success && passed_r1_category_succeeded);
        report_file_status("R1 Hamming Reject", final_output_map["hamming_rejected"], overall_r1_success);
        report_file_status("R1 BaseComp Reject", final_output_map["basecomp_rejected"], overall_r1_success);
        report_file_status("R2 Extracted/Sorted", final_r2_output_path, final_status_ok && passed_r1_category_succeeded);


        // Return appropriate exit code
        if (!final_status_ok) return 1;


    } catch (const cxxopts::exceptions::exception& e) {
        std::cerr << "Error parsing options: " << e.what() << std::endl;
        std::cerr << options.help() << std::endl;
        return 1;
    } catch (const fs::filesystem_error& e) {
         std::cerr << "Filesystem Error: " << e.what() << std::endl;
         std::cerr << "  Code: " << e.code().message() << std::endl;
         std::cerr << "  Path1: " << e.path1() << std::endl;
         if (!e.path2().empty()) std::cerr << "  Path2: " << e.path2() << std::endl;
         // Attempt cleanup of known intermediate file if it exists
         if (!intermediate_r2_unsorted_path.empty() && fs::exists(intermediate_r2_unsorted_path)) {
             std::cerr << "Attempting cleanup of intermediate file: " << intermediate_r2_unsorted_path << std::endl;
             std::error_code ec_rm; fs::remove(intermediate_r2_unsorted_path, ec_rm);
         }
         return 1;
    } catch (const std::runtime_error& e) {
        std::cerr << "Runtime Error: " << e.what() << std::endl;
         if (!intermediate_r2_unsorted_path.empty() && fs::exists(intermediate_r2_unsorted_path)) {
             std::cerr << "Attempting cleanup of intermediate file: " << intermediate_r2_unsorted_path << std::endl;
             std::error_code ec_rm; fs::remove(intermediate_r2_unsorted_path, ec_rm);
         }
        return 1;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
         if (!intermediate_r2_unsorted_path.empty() && fs::exists(intermediate_r2_unsorted_path)) {
             std::cerr << "Attempting cleanup of intermediate file: " << intermediate_r2_unsorted_path << std::endl;
             std::error_code ec_rm; fs::remove(intermediate_r2_unsorted_path, ec_rm);
         }
        return 1;
    } catch (...) {
        std::cerr << "An unknown error occurred." << std::endl;
         if (!intermediate_r2_unsorted_path.empty() && fs::exists(intermediate_r2_unsorted_path)) {
             std::cerr << "Attempting cleanup of intermediate file: " << intermediate_r2_unsorted_path << std::endl;
             std::error_code ec_rm; fs::remove(intermediate_r2_unsorted_path, ec_rm);
         }
        return 1;
    }

    std::cout << "\nTemporary directory cleanup initiated..." << std::flush;
    main_temp_dir_ptr.reset(); // Trigger RAII cleanup explicitly before returning
    std::cout << " Done." << std::endl;

    return 0; // Success
}
