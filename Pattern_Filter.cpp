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
#include <ctime>     // for std::time, std::strftime, std::localtime
#include <fstream>   // for std::ofstream (likely already present)
#include <csignal>
#include <cstdlib>
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
// MODIFIED: Simplified parameters for read chopping to match user's command
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
    int base_comp_start_pos = 0; // For base comp region
    int base_comp_length = 0;   // For base comp region
    int num_threads = 1;
    long long chunk_read_count = 0;
    fs::path main_temp_dir_base;
    size_t sort_memory_limit_mb = 0;
    int sort_threads = 1;

    // MODIFIED: Parameters for R1 read chopping simplified to a single region.
    bool chop_reads = false;
    int chop_start = -1;
    int chop_length = -1;
};


// === Forward Declarations for Helper Functions ===
void handle_incomplete_record(gzFile infile, const fs::path& chunk_path, const std::string& id_line, const std::string& missing_part);
void handle_read_error(igzstream& stream, const fs::path& filepath, long long record_num, const std::string& expected_part);
fs::path extract_r2_reads_parallel(
    const fs::path& sorted_r1_passed_path,
    const fs::path& original_r2_input_path,
    const fs::path& final_r2_output_path,
    const fs::path& main_temp_dir_base,
    const Params& params,
    std::unordered_set<std::string>& r2_found_ids_out); // NEW: out param for found R2 IDs
bool sort_and_merge_chunk_category(
    const std::vector<fs::path>& input_gz_chunks,
    const fs::path& final_output_gz_path,
    const fs::path& sort_temp_dir_base,
    const std::string& category_name,
    size_t sort_memory_limit_bytes,
    int num_sort_threads,
    bool enable_parallel_sort_alg);
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
bool merge_plain_text_chunks_to_gz(const std::vector<fs::path>& plain_chunk_files, const fs::path& final_output_gz_path);
bool concatenate_gz_files(const std::vector<fs::path>& input_files, const fs::path& output_file);

// NEW: Struct returned by R2 extraction chunk worker - now contains both the temp path AND the set of IDs found in this chunk.
struct R2ExtractChunkResult {
    fs::path output_path;
    std::unordered_set<std::string> found_ids;
};

R2ExtractChunkResult process_r2_extraction_chunk(
    const std::vector<std::string>& r2_chunk_records,
    const std::unordered_set<std::string>& r1_passed_ids_cref,
    const fs::path& r2_extract_temp_dir,
    int chunk_index);

// NEW: Forward declaration for the R1 filter-by-IDs function (two-way handshake)
bool filter_r1_by_r2_ids(
    const fs::path& r1_input_gz,
    const fs::path& r1_output_gz,
    const std::unordered_set<std::string>& r2_found_ids);
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
    }

    ~TempDir() {
        if (owns_ && !path_.empty()) {
             std::error_code ec_exists;
             if (fs::exists(path_, ec_exists)) {
                try {
                    std::error_code ec_remove;
                    // More robust check for valid temp dir names before removing
                    std::string filename = path_.filename().string();
                    bool safe_to_remove = filename.find("fastq_proc_temp_") != std::string::npos ||
                                          filename.find("filter_sort_main_") != std::string::npos ||
                                          filename.find("sort_temp_") != std::string::npos ||
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
// == Filtering Stage Structures and Functions (MODIFIED)               ==
// ========================================================================

inline int hamming_distance(const std::string& s1, const std::string& s2) {
    int dist = 0;
    size_t len = std::min(s1.length(), s2.length());
    for (size_t i = 0; i < len; ++i) {
        if (s1[i] != s2[i]) {
            dist++;
        }
    }
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

// MODIFIED: This function now contains the simplified R1 chopping logic.
UnsortedChunkOutputPaths process_filter_chunk(
    const std::vector<std::string>& reads_r1,
    const fs::path& filter_temp_dir_path,
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
            // MODIFIED: Make seq1 and qual1 mutable for potential chopping.
            const std::string& header1 = reads_r1[i];
            std::string seq1 = reads_r1[i+1];
            const std::string& plus1 = reads_r1[i+2];
            std::string qual1 = reads_r1[i+3];

             if (header1.empty() || header1[0] != '@' || plus1.empty() || plus1[0] != '+' || seq1.length() != qual1.length()) {
                  std::cerr << "Warning: Skipping malformed record in chunk " << chunk_index << " for ID " << header1.substr(0, std::min((size_t)50, header1.length())) << "..." << std::endl;
                 continue;
              }

            ogzstream* target_stream = nullptr;
            bool is_passed = false;

            // === Step 1: Determine Filter Category (Refactored logic) ===
            // Hamming Filter Check
            size_t required_hamming_len = static_cast<size_t>(params.start_pos + params.seq_length);
            if (seq1.length() < required_hamming_len) {
                target_stream = &hamming_out;
            } else {
                std::string hamming_sub_seq = seq1.substr(params.start_pos, params.seq_length);
                if (hamming_distance(hamming_sub_seq, params.target_seq) > params.threshold) {
                    target_stream = &hamming_out;
                } else {
                    // Base Composition Filter Check (only if passed Hamming)
                    size_t required_basecomp_len = static_cast<size_t>(params.base_comp_start_pos + params.base_comp_length);
                    if (seq1.length() < required_basecomp_len) {
                        target_stream = &basecomp_out;
                        std::cerr << "Warning: Read too short (" << seq1.length() << " bp) for base composition region (start "
                                  << params.base_comp_start_pos + 1 << ", length " << params.base_comp_length << ") "
                                  << "for ID " << header1.substr(0, std::min((size_t)50, header1.length())) << "... Filtering to basecomp_rejected." << std::endl;
                    } else {
                        std::string basecomp_sub_seq = seq1.substr(params.base_comp_start_pos, params.base_comp_length);
                        if (calculate_base_composition(basecomp_sub_seq) < params.base_composition_threshold) {
                            target_stream = &basecomp_out;
                        } else {
                            // Passed both filters
                            target_stream = &passed_out;
                            is_passed = true;
                        }
                    }
                }
            }

            // === Step 2: NEW - Perform R1 Chopping if applicable (with simplified logic) ===
            if (is_passed && params.chop_reads) {
                size_t required_chop_len = static_cast<size_t>(params.chop_start + params.chop_length);

                if (seq1.length() >= required_chop_len) {
                    // Extract the single specified region from sequence and quality strings
                    std::string new_seq = seq1.substr(params.chop_start, params.chop_length);
                    std::string new_qual = qual1.substr(params.chop_start, params.chop_length);

                    // Replace original seq1 and qual1 with the new chopped versions
                    seq1 = std::move(new_seq);
                    qual1 = std::move(new_qual);
                } else {
                    // Edge case: Read passed filters but is too short for the defined chopping region.
                    // This is a structural failure, so we reject it.
                    std::cerr << "Warning: Read passed filters but is too short for chopping (length "
                              << seq1.length() << "). Required length for chopping: " << required_chop_len
                              << ". ID: " << header1.substr(0, std::min((size_t)50, header1.length())) << ". Re-routing to basecomp_rejected." << std::endl;
                    target_stream = &basecomp_out; // Re-route to a rejection bin
                }
            }


            // === Step 3: Write the (potentially modified) record ===
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
// == Sorting Stage Structures and Functions (Unchanged)                ==
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
        total_mem += rec.id.capacity() + rec.seq.capacity() + rec.plus.capacity() + rec.qual.capacity() + 4;
    }
    total_mem += records.capacity() * sizeof(FastqRecord);
    return total_mem * 1.1;
}

bool write_plain_text_chunk(const std::vector<FastqRecord>& records, const fs::path& temp_filename) {
    std::ofstream outfile(temp_filename, std::ios::binary);
    if (!outfile) {
        std::cerr << "Error: Cannot open temporary plain text file for writing: " << temp_filename << " (" << strerror(errno) << ")" << std::endl;
        return false;
    }

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
    }

    PlainChunkReader(const PlainChunkReader&) = delete;
    PlainChunkReader& operator=(const PlainChunkReader&) = delete;
    PlainChunkReader(PlainChunkReader&&) = delete;
    PlainChunkReader& operator=(PlainChunkReader&&) = delete;

    bool read_next() {
        if (eof || !file) return false;
        current_record = {};

        if (std::getline(file, current_record.id) &&
            std::getline(file, current_record.seq) &&
            std::getline(file, current_record.plus) &&
            std::getline(file, current_record.qual))
        {
             if (current_record.id.empty() || current_record.id[0] != '@' || current_record.plus.empty() || current_record.plus[0] != '+' || current_record.seq.length() != current_record.qual.length()) {
                std::cerr << "Warning: Malformed record read from plain chunk " << filename << " for ID " << current_record.id.substr(0, std::min((size_t)50, current_record.id.length())) << "..." << std::endl;
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
        if (eof && !other.eof) return true;
        if (!eof && other.eof) return false;
        if (eof && other.eof) return false;
        return current_record.get_base_id() > other.current_record.get_base_id();
    }
};

struct PlainChunkReaderPtrCompare {
    bool operator()(const PlainChunkReader* a, const PlainChunkReader* b) const {
        return *a > *b;
    }
};

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
             std::error_code ec_rm; fs::remove(fname, ec_rm);
             continue;
        }

        auto reader = std::make_unique<PlainChunkReader>(fname);
        if (!reader->eof) {
            min_heap.push(reader.get());
            any_reader_valid = true;
            readers.push_back(std::move(reader));
        } else {
             std::cerr << "Warning: [Merge] Failed to initialize reader or read first record for " << fname << ". Skipping." << std::endl;
             reader.reset();
             std::error_code ec_rm; fs::remove(fname, ec_rm);
        }
    }

    if (!any_reader_valid) {
        std::cout << "Warning: No valid data found in any input chunk files for merge to " << final_output_gz_path << ". Creating empty file." << std::endl;
        readers.clear();
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
        return false;
    }
    gzbuffer(outfile_gz, GZ_BUFFER_SIZE);

    bool write_error = false;
    long long records_merged = 0;
    const long long MERGE_PROGRESS_INTERVAL = 5000000;
    std::string record_buffer;
    record_buffer.reserve(1024);

    while (!min_heap.empty()) {
        PlainChunkReader* top = min_heap.top();
        min_heap.pop();

        record_buffer.clear();
        record_buffer.append(top->current_record.id); record_buffer.push_back('\n');
        record_buffer.append(top->current_record.seq); record_buffer.push_back('\n');
        record_buffer.append(top->current_record.plus); record_buffer.push_back('\n');
        record_buffer.append(top->current_record.qual); record_buffer.push_back('\n');

        if (gzwrite(outfile_gz, record_buffer.c_str(), static_cast<unsigned int>(record_buffer.length())) != (int)record_buffer.length())
        {
            int errnum = 0;
            const char *error_str = gzerror(outfile_gz, &errnum);
            std::cerr << "\nError: Failed writing to gzipped output file " << final_output_gz_path
                      << ". Error (" << std::to_string(errnum) << "): " << (error_str ? error_str : "Unknown") << std::endl;
            write_error = true;
            break;
        }

        records_merged++;
        if (records_merged % MERGE_PROGRESS_INTERVAL == 0) {
            std::cout << "." << std::flush;
        }

        if (top->read_next()) {
            min_heap.push(top);
        } else {
            if (top->file.bad()) {
                 std::cerr << "\nError reading from chunk file during merge: " << top->filename << std::endl;
                 write_error = true;
                 break;
            }
        }
    }
    std::cout << std::endl;

    while (!min_heap.empty()) min_heap.pop();

    int close_status = gzclose(outfile_gz);

    if (close_status != Z_OK) {
        if (close_status != Z_BUF_ERROR || records_merged == 0) {
             std::cerr << "Error: Failed to properly close gzipped output file " << final_output_gz_path
                       << ". Close status: " << close_status << ". Output may be corrupted." << std::endl;
             write_error = true;
        } else {
             std::cerr << "Warning: gzclose returned status " << close_status << " on " << final_output_gz_path
                       << ". File might be usable, but check integrity." << std::endl;
        }
    }

    bool final_success = !write_error;
    if (final_success) {
        std::cout << "Merge successful for " << final_output_gz_path << ". Merged " << records_merged << " records." << std::endl;
        for (const auto& ptr : readers) {
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
         std::error_code ec;
         fs::remove(final_output_gz_path, ec);
    }

    readers.clear();
    return final_success;
}


// ========================================================================
// == NEW: Parallel Sorting Phase 1 Worker Function (Unchanged)         ==
// ========================================================================
SortWorkerResult process_sort_sub_chunk(
    const std::vector<fs::path>& sub_input_gz_chunks,
    const fs::path& category_sort_temp_path,
    const std::string& category_name,
    size_t sort_memory_limit_bytes,
    bool enable_parallel_sort_alg,
    int worker_id)
{
    SortWorkerResult result;
    result.success = false;
    long long worker_total_records = 0;
    int plain_chunk_count = 0;
    std::vector<FastqRecord> current_sort_chunk;
    size_t avg_record_size_bytes = 500;
    size_t estimated_records_per_chunk = std::max((size_t)10000, sort_memory_limit_bytes / avg_record_size_bytes);
    current_sort_chunk.reserve(estimated_records_per_chunk);
    size_t current_memory_estimate = 0;
    char read_buffer[GZ_BUFFER_SIZE];

    try {
        for(const auto& gz_chunk_path : sub_input_gz_chunks) {
            std::error_code ec_exists;
            if (!fs::exists(gz_chunk_path, ec_exists)) {
                std::cerr << "Warning: [Worker " << worker_id << "] Input chunk file not found, skipping: " << gz_chunk_path << std::endl;
                continue;
            }
             std::error_code ec_sz; uintmax_t fsize = fs::file_size(gz_chunk_path, ec_sz);
             if (ec_sz || fsize == 0) {
                 if (ec_sz) std::cerr << "Warning: [Worker " << worker_id << "] Could not get size of " << gz_chunk_path << ", skipping: " << ec_sz.message() << std::endl;
                 else std::cerr << "Warning: [Worker " << worker_id << "] Input chunk file " << gz_chunk_path << " is empty, skipping." << std::endl;
                 if (!ec_sz && fsize == 0) { std::error_code ec_rm; fs::remove(gz_chunk_path, ec_rm); }
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
                    if (gzeof(infile)) break;
                    int err; const char* msg = gzerror(infile, &err);
                    if (err != Z_OK && err != Z_STREAM_END) {
                         std::cerr << "\nError: [Worker " << worker_id << "] Reading ID line from chunk " << gz_chunk_path << ": " << (msg ? msg : "Unknown error") << " (err code: " << err << ")" << std::endl;
                         gzclose(infile);
                         throw std::runtime_error("Worker " + std::to_string(worker_id) + ": Error reading gz chunk file (ID line).");
                    }
                     if (!gzeof(infile)) {
                         std::cerr << "\nWarning: [Worker " << worker_id << "] gzgets returned NULL but not EOF and no major error on " << gz_chunk_path << ". Attempting to continue." << std::endl;
                         continue;
                     }
                    break;
                }
                line_id = read_buffer;
                while (!line_id.empty() && (line_id.back() == '\n' || line_id.back() == '\r')) {
                    line_id.pop_back();
                }
                 if (line_id.empty() || line_id[0] != '@') {
                     std::cerr << "\nWarning: [Worker " << worker_id << "] Malformed ID line encountered in " << gz_chunk_path << ". Skipping record. Line: '" << line_id.substr(0, std::min((size_t)100, line_id.length())) << "'" << std::endl;
                     for(int skip=0; skip<3; ++skip) {
                         if (gzgets(infile, read_buffer, sizeof(read_buffer)) == nullptr) {
                             if (!gzeof(infile)) {
                                  int err; gzerror(infile, &err);
                                  std::cerr << "   (Failed to skip subsequent lines, error code: " << err << ")" << std::endl;
                                  gzclose(infile);
                                  throw std::runtime_error("Worker " + std::to_string(worker_id) + ": Error skipping lines after malformed ID.");
                             }
                              break;
                         }
                     }
                     continue;
                 }

                auto read_next_line = [&](std::string& target_line, const char* line_name) -> bool {
                     if (gzgets(infile, read_buffer, sizeof(read_buffer)) == nullptr) {
                         handle_incomplete_record(infile, gz_chunk_path, line_id, line_name);
                         return false;
                     }
                     target_line = read_buffer;
                     while (!target_line.empty() && (target_line.back() == '\n' || target_line.back() == '\r')) {
                          target_line.pop_back();
                     }
                     return true;
                };

                if (!read_next_line(line_seq, "sequence")) { break; }
                if (!read_next_line(line_plus, "plus")) { break; }
                if (!read_next_line(line_qual, "quality")) { break; }

                if (line_plus.empty() || line_plus[0] != '+' || line_seq.length() != line_qual.length()) {
                     std::cerr << "\nWarning: [Worker " << worker_id << "] Malformed record encountered in " << gz_chunk_path << " for ID " << line_id.substr(0, std::min((size_t)50, line_id.length())) << ". Skipping." << std::endl;
                     continue;
                 }

                current_sort_chunk.emplace_back(FastqRecord{std::move(line_id), std::move(line_seq), std::move(line_plus), std::move(line_qual)});
                worker_total_records++;

                 bool memory_check_needed = (current_sort_chunk.size() % (estimated_records_per_chunk / 4 + 1)) == 0;
                 if (current_sort_chunk.size() >= estimated_records_per_chunk || (memory_check_needed && current_sort_chunk.size() > 1000) )
                 {
                    current_memory_estimate = estimate_memory(current_sort_chunk);

                    if (current_memory_estimate >= sort_memory_limit_bytes || current_sort_chunk.size() >= estimated_records_per_chunk * 1.5)
                    {
                        if (enable_parallel_sort_alg) {
                            #ifdef USE_PARALLEL_SORT
                                try {
                                    std::sort(std::execution::par, current_sort_chunk.begin(), current_sort_chunk.end());
                                } catch (const std::exception& e) {
                                    std::cerr << "\nWarning: [Worker " << worker_id << "] Parallel sort failed: " << e.what() << ". Falling back to sequential sort for this chunk." << std::endl;
                                    std::sort(current_sort_chunk.begin(), current_sort_chunk.end());
                                } catch (...) {
                                    std::cerr << "\nWarning: [Worker " << worker_id << "] Parallel sort failed with unknown exception. Falling back to sequential sort for this chunk." << std::endl;
                                     std::sort(current_sort_chunk.begin(), current_sort_chunk.end());
                                }
                            #else
                                std::sort(current_sort_chunk.begin(), current_sort_chunk.end());
                            #endif
                        } else {
                            std::sort(current_sort_chunk.begin(), current_sort_chunk.end());
                        }

                        fs::path temp_filename = category_sort_temp_path /
                            ("sorted_" + category_name + "_worker_" + std::to_string(worker_id) + "_chunk_" + std::to_string(plain_chunk_count++) + ".tmp");

                        if (!write_plain_text_chunk(current_sort_chunk, temp_filename)) {
                            gzclose(infile);
                            throw std::runtime_error("Worker " + std::to_string(worker_id) + ": Error writing sorted plain text chunk file: " + temp_filename.string());
                        }
                        result.created_temp_files.push_back(temp_filename);

                        current_sort_chunk.clear();
                        current_sort_chunk.reserve(estimated_records_per_chunk);
                        current_memory_estimate = 0;

                        if (plain_chunk_count % 5 == 0 && worker_total_records > 0) {
                             size_t approx_mem_last_chunk = current_memory_estimate > 0 ? current_memory_estimate : sort_memory_limit_bytes;
                             avg_record_size_bytes = std::max((size_t)50, approx_mem_last_chunk / std::max((size_t)1, estimated_records_per_chunk));
                             estimated_records_per_chunk = std::max((size_t)10000, sort_memory_limit_bytes / std::max((size_t)1, avg_record_size_bytes));
                             current_sort_chunk.reserve(estimated_records_per_chunk);
                        }
                    }
                }
            }

            gzclose(infile);

             std::error_code ec_exists_del;
             if (fs::exists(gz_chunk_path, ec_exists_del)) {
                std::error_code ec_rm;
                fs::remove(gz_chunk_path, ec_rm);
                if (ec_rm) {
                     if (category_name != "r2_extracted") {
                         std::cerr << "\nWarning: [Worker " << worker_id << "] Failed to remove processed input chunk " << gz_chunk_path << ": " << ec_rm.message() << std::endl;
                     }
                }
             }

        }

        if (!current_sort_chunk.empty()) {
            if (enable_parallel_sort_alg) {
                #ifdef USE_PARALLEL_SORT
                    try {
                        std::sort(std::execution::par, current_sort_chunk.begin(), current_sort_chunk.end());
                    } catch (const std::exception& e) {
                         std::cerr << "\nWarning: [Worker " << worker_id << "] Parallel sort failed (final chunk): " << e.what() << ". Falling back to sequential sort." << std::endl;
                        std::sort(current_sort_chunk.begin(), current_sort_chunk.end());
                    } catch (...) {
                         std::cerr << "\nWarning: [Worker " << worker_id << "] Parallel sort failed (final chunk) with unknown exception. Falling back to sequential sort." << std::endl;
                         std::sort(current_sort_chunk.begin(), current_sort_chunk.end());
                    }
                #else
                    std::sort(current_sort_chunk.begin(), current_sort_chunk.end());
                #endif
            } else {
                 std::sort(current_sort_chunk.begin(), current_sort_chunk.end());
            }

            fs::path temp_filename = category_sort_temp_path /
                ("sorted_" + category_name + "_worker_" + std::to_string(worker_id) + "_chunk_" + std::to_string(plain_chunk_count++) + ".tmp");
             if (!write_plain_text_chunk(current_sort_chunk, temp_filename)) {
                 throw std::runtime_error("Worker " + std::to_string(worker_id) + ": Error writing final sorted plain text chunk file: " + temp_filename.string());
             }
            result.created_temp_files.push_back(temp_filename);
            current_sort_chunk.clear();
        }

        result.success = true;
        result.records_processed = worker_total_records;

    } catch (const std::exception& e) {
        std::cerr << "\nError in sort worker " << worker_id << " for category '" << category_name << "': " << e.what() << std::endl;
        result.success = false;
        for(const auto& temp_f : result.created_temp_files) {
            std::error_code ec; fs::remove(temp_f, ec);
        }
        result.created_temp_files.clear();
        result.records_processed = worker_total_records;
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
// == MODIFIED: Main Sorting Function (Orchestrator - Unchanged)        ==
// ========================================================================
bool sort_and_merge_chunk_category(
    const std::vector<fs::path>& input_gz_chunks,
    const fs::path& final_output_gz_path,
    const fs::path& main_temp_dir_base,
    const std::string& category_name,
    size_t sort_memory_limit_bytes_per_worker,
    int num_sort_threads,
    bool enable_parallel_sort_alg)
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

    fs::path category_sort_temp_path = main_temp_dir_base / ("sort_temp_" + category_name);
    std::unique_ptr<TempDir> category_temp_dir_ptr;
    try {
       category_temp_dir_ptr = std::make_unique<TempDir>(main_temp_dir_base, "sort_temp_" + category_name + "_");
       category_sort_temp_path = category_temp_dir_ptr->getPath();
       std::cout << "Using temporary directory for sorting intermediate files: " << category_sort_temp_path << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error creating temporary directory for sorting category '" << category_name << "': " << e.what() << std::endl;
        return false;
    }


    long long total_records_in_category = 0;
    std::vector<fs::path> all_sorted_plain_chunk_files;
    std::mutex result_mutex;
    std::atomic<int> failed_worker_count = 0;

    std::cout << "Phase 1 (Parallel Sort/Write): Distributing " << input_gz_chunks.size() << " input chunks to " << num_sort_threads << " workers..." << std::endl;
    std::vector<std::future<SortWorkerResult>> sort_futures;
    std::vector<std::vector<fs::path>> worker_inputs(num_sort_threads);

    for (size_t i = 0; i < input_gz_chunks.size(); ++i) {
        worker_inputs[i % num_sort_threads].push_back(input_gz_chunks[i]);
    }

    for (int i = 0; i < num_sort_threads; ++i) {
        if (!worker_inputs[i].empty()) {
            sort_futures.push_back(std::async(std::launch::async, process_sort_sub_chunk,
                                               std::move(worker_inputs[i]),
                                               category_sort_temp_path,
                                               category_name,
                                               sort_memory_limit_bytes_per_worker,
                                               enable_parallel_sort_alg,
                                               i));
        }
    }

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
            SortWorkerResult worker_result = fut.get();
            completed_workers++;
            std::cout << "." << std::flush;

            std::lock_guard<std::mutex> lock(result_mutex);
            if (worker_result.success) {
                all_sorted_plain_chunk_files.insert(all_sorted_plain_chunk_files.end(),
                                                    std::make_move_iterator(worker_result.created_temp_files.begin()),
                                                    std::make_move_iterator(worker_result.created_temp_files.end()));
                total_records_in_category += worker_result.records_processed;
            } else {
                failed_worker_count++;
            }
        } catch (const std::exception& e) {
            std::cerr << "\nError retrieving result from sort worker: " << e.what() << std::endl;
            failed_worker_count++;
            completed_workers++;
        } catch (...) {
            std::cerr << "\nUnknown error retrieving result from sort worker." << std::endl;
            failed_worker_count++;
            completed_workers++;
        }
    }
    std::cout << " Done." << std::endl;

    if (failed_worker_count > 0) {
         std::cerr << "Error: " << failed_worker_count << " sort worker(s) failed for category '" << category_name << "'." << std::endl;
         for(const auto& temp_f : all_sorted_plain_chunk_files) {
             std::error_code ec; fs::remove(temp_f, ec);
         }
         return false;
    }

    std::cout << "Phase 1 (Parallel Sort/Write): Finished. Total records processed: " << total_records_in_category
                << ". Created " << all_sorted_plain_chunk_files.size() << " sorted plain text chunk(s)." << std::endl;


    std::cout << "Phase 2 (Sequential Merge): Merging sorted plain text chunks..." << std::endl;
    if (!merge_plain_text_chunks_to_gz(all_sorted_plain_chunk_files, final_output_gz_path)) {
         std::cerr << "Error: Merging of sorted plain text chunks failed for category '" << category_name << "'." << std::endl;
         return false;
    }

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
          error_base += " Read failed but not EOF or GZLIB error.";
     }
     std::cerr << "\n" << error_base << std::endl;

    throw std::runtime_error("Incomplete FASTQ record detected in gz chunk: " + chunk_path.string() + " ID: " + id_line.substr(0,std::min((size_t)50, id_line.length())));
}

void handle_read_error(igzstream& stream, const fs::path& filepath, long long record_num, const std::string& expected_part) {
    std::string error_msg;
    bool was_eof = stream.eof();
    bool was_bad = stream.bad();
    bool was_fail = stream.fail();

    if (was_eof && was_fail && !was_bad) {
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
    else {
         error_msg = "Stream logical error (failbit set) reading " + expected_part + " line from input file " + filepath.string() +
                     " near record " + std::to_string(record_num + 1) + ". Stream state (eof,fail,bad): " +
                     std::to_string(was_eof) + "," + std::to_string(was_fail) + "," + std::to_string(was_bad);
    }
     throw std::runtime_error(error_msg);
}

// =====================================================================
// RAII guard: ensures temp directory is removed on ANY scope exit
// (normal return, exception, std::exit, etc.)
// =====================================================================
class TempDirGuard {
public:
    fs::path dir;
    bool armed = true;
    explicit TempDirGuard(const fs::path& p) : dir(p) {}
    void disarm() { armed = false; }
    ~TempDirGuard() {
        if (armed && !dir.empty() && fs::exists(dir)) {
            std::error_code ec;
            fs::remove_all(dir, ec);
            if (ec) {
                std::cerr << "[Cleanup] Warning: could not remove temp dir "
                          << dir << " (" << ec.message() << ")\n";
            } else {
                std::cerr << "[Cleanup] Removed temp dir: " << dir << "\n";
            }
        }
    }
};

// Global pointer so signal handler can reach it
static fs::path g_temp_dir_for_signal;

static void signal_cleanup_handler(int sig) {
    if (!g_temp_dir_for_signal.empty() && fs::exists(g_temp_dir_for_signal)) {
        std::error_code ec;
        fs::remove_all(g_temp_dir_for_signal, ec);
    }
    // Restore default handler and re-raise so exit status is correct
    std::signal(sig, SIG_DFL);
    std::raise(sig);
}
// ========================================================================
// == R2 Extraction Function (MODIFIED: tracks found IDs for handshake) ==
// ========================================================================
// MODIFIED: now returns R2ExtractChunkResult containing both the output path and the set of IDs that were actually found and written.
R2ExtractChunkResult process_r2_extraction_chunk(
    const std::vector<std::string>& r2_chunk_records,
    const std::unordered_set<std::string>& r1_passed_ids_cref,
    const fs::path& r2_extract_temp_dir,
    int chunk_index)
{
    R2ExtractChunkResult result;
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
            const std::string& header_line = r2_chunk_records[i];
            const std::string& seq_line = r2_chunk_records[i + 1];
            const std::string& plus_line = r2_chunk_records[i + 2];
            const std::string& qual_line = r2_chunk_records[i + 3];

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
                if (gzputs(out_gz, header_line.c_str()) == -1 ||
                    gzputs(out_gz, seq_line.c_str()) == -1 ||
                    gzputs(out_gz, plus_line.c_str()) == -1 ||
                    gzputs(out_gz, qual_line.c_str()) == -1)
                {
                    int errnum = 0; const char *error_str = gzerror(out_gz, &errnum);
                    throw std::runtime_error("Error writing extracted R2 record for ID " + header_no_nl.substr(0,std::min((size_t)50, header_no_nl.length())) + " to temp chunk " + temp_chunk_path.string() + ". Error (" + std::to_string(errnum) + "): " + (error_str ? error_str : "Unknown"));
                }
                records_written++;
                // NEW: Record which IDs were actually found and written for the two-way handshake
                result.found_ids.insert(base_id);
            }
        }

        int close_status = gzclose(out_gz);
        out_gz = nullptr;
        if (close_status != Z_OK) {
             std::cerr << "Warning: Failed to properly close temporary R2 chunk file " << temp_chunk_path
                       << ". Close status: " << close_status << ". File might be corrupted." << std::endl;
             if(close_status != Z_BUF_ERROR) {
                 std::error_code ec_rm; fs::remove(temp_chunk_path, ec_rm);
                 result.output_path = fs::path{};
                 result.found_ids.clear();
                 return result;
             }
        }

        if (records_written == 0) {
             std::error_code ec_rm; fs::remove(temp_chunk_path, ec_rm);
             result.output_path = fs::path{};
             return result;
        }
        result.output_path = temp_chunk_path;
        return result;

    } catch (const std::exception& e) {
        std::cerr << "Error in process_r2_extraction_chunk (Index " << chunk_index << "): " << e.what() << std::endl;
        if (out_gz) gzclose(out_gz);
        std::error_code ec_rm; fs::remove(temp_chunk_path, ec_rm);
        result.output_path = fs::path{};
        result.found_ids.clear();
        return result;
    } catch (...) {
        std::cerr << "Unknown error in process_r2_extraction_chunk (Index " << chunk_index << ")" << std::endl;
        if (out_gz) gzclose(out_gz);
        std::error_code ec_rm; fs::remove(temp_chunk_path, ec_rm);
        result.output_path = fs::path{};
        result.found_ids.clear();
        return result;
    }
}

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
        if (bytes_read < 0) {
            const char* error_str = gzerror(in_gz, &read_err);
            std::cerr << "\nError: Failed reading data from input chunk " << in_path
                      << " during concatenation. Error (" << read_err << "): " << (error_str ? error_str : "Unknown") << ". Aborting." << std::endl;
            success = false;
        }

        gzclose(in_gz);

        if (!success) {
            break;
        }
        if ((i + 1) % 10 == 0 || (i + 1) == input_files.size()) {
             std::cout << "." << std::flush;
        }
    }

    int close_status = gzclose(out_gz);
    std::cout << " Done." << std::endl;

    if (close_status != Z_OK) {
         std::cerr << "Warning: Failed to properly close concatenated output file " << output_file
                   << ". Close status: " << close_status << ". File might be corrupted." << std::endl;
         if (success && total_bytes_written > 0 && close_status != Z_BUF_ERROR) {
             success = false;
         }
    }

    if (success) {
         std::cout << "Concatenation successful. Total bytes written: " << total_bytes_written << std::endl;
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
         std::error_code ec_rm; fs::remove(output_file, ec_rm);
         return false;
    }

    return true;
}

// MODIFIED: now collects the set of R2 IDs that were actually found and returns them via r2_found_ids_out.
fs::path extract_r2_reads_parallel(
    const fs::path& sorted_r1_passed_path,
    const fs::path& original_r2_input_path,
    const fs::path& final_r2_output_path,
    const fs::path& main_temp_dir_base,
    const Params& params,
    std::unordered_set<std::string>& r2_found_ids_out)
{
    std::cout << "\n--- Starting R2 Extraction Phase (Parallel) ---" << std::endl;
    std::cout << "Using R1 passed IDs from: " << sorted_r1_passed_path << std::endl;
    std::cout << "Scanning original R2 file: " << original_r2_input_path << std::endl;
    std::cout << "Using " << params.num_threads << " threads for R2 extraction." << std::endl;

    r2_found_ids_out.clear();

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
         std::error_code ec_r1_sz; uintmax_t r1_size = fs::file_size(sorted_r1_passed_path, ec_r1_sz);
         if (!ec_r1_sz && r1_size < 50) {
             std::cout << "Warning: R1 passed file appears empty or very small (" << r1_size << " bytes). Skipping R2 extraction." << std::endl;
             return fs::path{};
         } else if (ec_r1_sz) {
              std::cerr << "Warning: Could not determine size of R1 passed file (" << ec_r1_sz.message() << "). Proceeding with R2 extraction, but file might be empty." << std::endl;
         }


        fs::path out_dir = final_r2_output_path.parent_path();
        std::string base_name = final_r2_output_path.filename().string();
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
             base_name = "extracted";
          }


        intermediate_r2_path = out_dir / (base_name + "_R2_unsorted_temp.fastq.gz");
        std::cout << "Final intermediate unsorted R2 file will be: " << intermediate_r2_path << std::endl;

        r2_extract_temp_dir_ptr = std::make_unique<TempDir>(main_temp_dir_base, "r2_extract_chunks_");
        const fs::path& r2_extract_temp_path = r2_extract_temp_dir_ptr->getPath();
        std::cout << "Using temporary directory for R2 extraction chunks: " << r2_extract_temp_path << std::endl;

        std::unordered_set<std::string> r1_passed_ids;
        long long r1_id_count = 0;
        const long long R1_ID_READ_PROGRESS = 2000000;
        std::cout << "Step 1: Reading R1 passed IDs..." << std::flush;

        r1_in_gz = gzopen(sorted_r1_passed_path.c_str(), "rb");
        if (!r1_in_gz) {
            throw std::runtime_error("Failed to open R1 passed file: " + sorted_r1_passed_path.string() + " (" + strerror(errno) + ")");
        }
        gzbuffer(r1_in_gz, GZ_BUFFER_SIZE);
        char buffer[GZ_BUFFER_SIZE];
        std::string line;

        while (gzgets(r1_in_gz, buffer, sizeof(buffer)) != nullptr) {
            line = buffer;
            size_t len = line.length();
            while (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r')) {
                len--;
            }
            line.resize(len);

            if (!line.empty() && line[0] == '@') {
                 r1_passed_ids.insert(get_base_id_from_string(line));
                 r1_id_count++;
                 if (r1_id_count % R1_ID_READ_PROGRESS == 0) std::cout << "." << std::flush;

                 for (int i = 0; i < 3; ++i) {
                     if (gzgets(r1_in_gz, buffer, sizeof(buffer)) == nullptr) {
                         if (!gzeof(r1_in_gz)) {
                             int err = 0; gzerror(r1_in_gz, &err);
                             gzclose(r1_in_gz);
                             throw std::runtime_error("Incomplete record found while reading R1 passed IDs from " + sorted_r1_passed_path.string() + " after ID " + line.substr(0, std::min((size_t)50, line.length())) + ". GZ Error: " + std::to_string(err));
                         }
                         if (i < 2) {
                             gzclose(r1_in_gz);
                             throw std::runtime_error("Incomplete record at end of R1 passed file: " + sorted_r1_passed_path.string() + " after ID " + line.substr(0, std::min((size_t)50, line.length())));
                         }
                         goto end_r1_read;
                     }
                 }

            } else if (!line.empty()) {
                 std::cerr << "\nWarning: Skipping unexpected line in R1 passed file: " << line.substr(0, std::min((size_t)100, line.length())) << std::endl;
            }
        }

        end_r1_read:;
        gzclose(r1_in_gz); r1_in_gz = nullptr;
        std::cout << " Done. Read " << r1_id_count << " IDs. Unique IDs in set: " << r1_passed_ids.size() << "." << std::endl;

        if (r1_id_count == 0 || r1_passed_ids.empty()) {
             std::cout << "No passed R1 IDs found or read. Skipping R2 extraction. No R2 output file will be created." << std::endl;
             r2_extract_temp_dir_ptr.reset();
             return fs::path{};
        }
        if (r1_id_count != (long long)r1_passed_ids.size()) {
             std::cout << "Note: Number of IDs read (" << r1_id_count << ") differs from final set size (" << r1_passed_ids.size() << "). Duplicates existed in the sorted R1 file." << std::endl;
        }


        std::cout << "Step 2: Reading R2 input and dispatching extraction tasks..." << std::flush;
        r2_in_gz = gzopen(original_r2_input_path.c_str(), "rb");
         if (!r2_in_gz) {
            throw std::runtime_error("Failed to open original R2 input file: " + original_r2_input_path.string() + " (" + strerror(errno) + ")");
        }
        gzbuffer(r2_in_gz, GZ_BUFFER_SIZE);

        std::vector<std::future<R2ExtractChunkResult>> r2_extract_futures;
        std::vector<std::string> current_r2_chunk;
        current_r2_chunk.reserve(params.chunk_read_count * 4);
        long long record_count_in_chunk = 0;
        long long total_r2_records_read = 0;
        int chunk_index = 0;
        bool read_error = false;
        char line_buffer[GZ_BUFFER_SIZE];

        while(true) {
             std::string line1, line2, line3, line4;
             std::string current_id_for_error;

             if (gzgets(r2_in_gz, line_buffer, sizeof(line_buffer)) == nullptr) {
                 if(gzeof(r2_in_gz)) break;
                 int err=0; const char* msg = gzerror(r2_in_gz, &err);
                 std::cerr << "\nError reading ID line from R2 input " << original_r2_input_path.string() << ". GZ Error: " << std::to_string(err) << " " << (msg ? msg : "N/A") << std::endl;
                 read_error = true;
                 break;
            }
            line1 = line_buffer;
            current_id_for_error = line1;
            while (!current_id_for_error.empty() && (current_id_for_error.back() == '\n' || current_id_for_error.back() == '\r')) {
                current_id_for_error.pop_back();
            }

            if (gzgets(r2_in_gz, line_buffer, sizeof(line_buffer)) == nullptr) {
                handle_incomplete_record(r2_in_gz, original_r2_input_path, current_id_for_error, "sequence");
                read_error = true; break;
            }
            line2 = line_buffer;

            if (gzgets(r2_in_gz, line_buffer, sizeof(line_buffer)) == nullptr) {
                handle_incomplete_record(r2_in_gz, original_r2_input_path, current_id_for_error, "plus");
                read_error = true; break;
            }
            line3 = line_buffer;

            if (gzgets(r2_in_gz, line_buffer, sizeof(line_buffer)) == nullptr) {
                 handle_incomplete_record(r2_in_gz, original_r2_input_path, current_id_for_error, "quality");
                 read_error = true; break;
             }
             line4 = line_buffer;

             current_r2_chunk.push_back(std::move(line1));
             current_r2_chunk.push_back(std::move(line2));
             current_r2_chunk.push_back(std::move(line3));
             current_r2_chunk.push_back(std::move(line4));

             record_count_in_chunk++;
             total_r2_records_read++;
             if (total_r2_records_read % 5000000 == 0) std::cout << "." << std::flush;

             if (record_count_in_chunk >= params.chunk_read_count) {
                 r2_extract_futures.push_back(std::async(std::launch::async, process_r2_extraction_chunk,
                                             std::move(current_r2_chunk),
                                             std::cref(r1_passed_ids),
                                             r2_extract_temp_path,
                                             chunk_index++));

                 current_r2_chunk = {};
                 current_r2_chunk.reserve(params.chunk_read_count * 4);
                 record_count_in_chunk = 0;

                 while (r2_extract_futures.size() >= static_cast<size_t>(params.num_threads * 2 + 4)) {
                      bool found_ready = false;
                      for(auto it = r2_extract_futures.begin(); it != r2_extract_futures.end(); ) {
                          if (it->valid() && it->wait_for(std::chrono::milliseconds(0)) == std::future_status::ready) {
                              try {
                                  R2ExtractChunkResult res = it->get();
                                  // collect IDs as soon as available
                                  for (auto& id : res.found_ids) r2_found_ids_out.insert(std::move(id));
                              } catch (...) { }
                              it = r2_extract_futures.erase(it);
                              found_ready = true;
                          } else if (!it->valid()) {
                               it = r2_extract_futures.erase(it);
                          }
                          else {
                              ++it;
                          }
                      }
                      if (!found_ready && !r2_extract_futures.empty()) {
                            std::this_thread::sleep_for(std::chrono::milliseconds(20));
                      }
                 }
             }
        }

        gzclose(r2_in_gz); r2_in_gz = nullptr;

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
             for(auto& fut : r2_extract_futures) { if(fut.valid()) { fut.wait(); } }
             throw std::runtime_error("Error occurred during R2 file reading. Aborting extraction.");
        }

        std::vector<fs::path> extracted_r2_chunk_files;
        extracted_r2_chunk_files.reserve(r2_extract_futures.size());
        int failed_chunks = 0;

        std::cout << "Step 3: Waiting for " << r2_extract_futures.size() << " R2 extraction tasks to complete..." << std::flush;
        size_t total_futures = r2_extract_futures.size();
        size_t completed_count = 0;
        for (auto& fut : r2_extract_futures) {
             if (!fut.valid()) {
                 std::cerr << "\nWarning: An R2 extraction future was found invalid before getting result." << std::endl;
                 failed_chunks++; completed_count++; continue;
             }
             try {
                 R2ExtractChunkResult res = fut.get();
                 // NEW: merge found IDs from this chunk into the global set
                 for (auto& id : res.found_ids) r2_found_ids_out.insert(std::move(id));

                 fs::path chunk_path = res.output_path;
                 if (!chunk_path.empty()) {
                     std::error_code ec_exists;
                     if (fs::exists(chunk_path, ec_exists)) {
                        std::error_code ec_sz; uintmax_t fsize = fs::file_size(chunk_path, ec_sz);
                        if (!ec_sz && fsize > 0) {
                             extracted_r2_chunk_files.push_back(chunk_path);
                        } else {
                            if (!ec_sz && fsize == 0) {
                                 std::error_code ec_rm; fs::remove(chunk_path, ec_rm);
                             } else {
                                 std::cerr << "\nWarning: R2 extraction chunk exists but size check failed or is zero: " << chunk_path << " (Error: " << ec_sz.message() << "). Skipping." << std::endl;
                                 failed_chunks++;
                             }
                        }
                     } else {
                        std::cerr << "\nWarning: R2 extraction chunk file reported success but not found: " << chunk_path << ". Worker might have failed cleanup." << std::endl;
                     }
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
                  std::cout << "." << std::flush;
              }
        }
        r2_extract_futures.clear();
        std::cout << " Done." << std::endl;

        if (failed_chunks > 0) {
             std::cerr << "Warning: " << failed_chunks << " R2 extraction chunk task(s) reported errors or produced unusable files." << std::endl;
        }
        if (extracted_r2_chunk_files.empty()) {
            std::cout << "No R2 records were extracted (no matching IDs found or all tasks failed/produced empty files)." << std::endl;
            r2_extract_temp_dir_ptr.reset();
            return fs::path{};
        }

        std::cout << "Collected " << extracted_r2_chunk_files.size() << " non-empty R2 extracted chunk file(s)." << std::endl;
        std::cout << "Two-way handshake: collected " << r2_found_ids_out.size() << " unique R2 IDs that matched R1 passed IDs." << std::endl;

        std::cout << "Step 4: Concatenating temporary R2 chunk files..." << std::endl;
        if (!concatenate_gz_files(extracted_r2_chunk_files, intermediate_r2_path)) {
             throw std::runtime_error("Failed to concatenate temporary R2 chunk files into: " + intermediate_r2_path.string());
        }

        std::cout << "Step 5: Cleaning up R2 extraction temporary directory..." << std::flush;
        r2_extract_temp_dir_ptr.reset();
        std::cout << " Done." << std::endl;

        std::cout << "Successfully created intermediate unsorted R2 file: " << intermediate_r2_path << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "\nError during parallel R2 extraction: " << e.what() << std::endl;
        if (r1_in_gz) gzclose(r1_in_gz);
        if (r2_in_gz) gzclose(r2_in_gz);
        std::error_code ec; fs::remove(intermediate_r2_path, ec);
        return fs::path{};
    } catch (...) {
         std::cerr << "\nUnknown error during parallel R2 extraction." << std::endl;
         if (r1_in_gz) gzclose(r1_in_gz);
         if (r2_in_gz) gzclose(r2_in_gz);
         std::error_code ec; fs::remove(intermediate_r2_path, ec);
         return fs::path{};
    }

    std::cout << "--- R2 Extraction Phase (Parallel) Completed Successfully ---" << std::endl;
    return intermediate_r2_path;
}

// ========================================================================
// == NEW: filter_r1_by_r2_ids - Two-Way Handshake R1 Filter            ==
// ========================================================================
// Reads the sorted R1 passed file and outputs a new file containing ONLY those reads whose
// base IDs are present in r2_found_ids. The output file replaces the input atomically.
bool filter_r1_by_r2_ids(
    const fs::path& r1_input_gz,
    const fs::path& r1_output_gz,
    const std::unordered_set<std::string>& r2_found_ids)
{
    std::cout << "\n--- Two-Way Handshake: Filtering R1 by R2 found IDs ---" << std::endl;
    std::cout << "R1 input (to be re-filtered): " << r1_input_gz << std::endl;
    std::cout << "R2 found ID set size: " << r2_found_ids.size() << std::endl;

    gzFile in_gz = gzopen(r1_input_gz.c_str(), "rb");
    if (!in_gz) {
        std::cerr << "Error: Cannot open R1 input for handshake filtering: " << r1_input_gz << " (" << strerror(errno) << ")" << std::endl;
        return false;
    }
    gzbuffer(in_gz, GZ_BUFFER_SIZE);

    gzFile out_gz = gzopen(r1_output_gz.c_str(), "wb");
    if (!out_gz) {
        std::cerr << "Error: Cannot open R1 output for handshake filtering: " << r1_output_gz << " (" << strerror(errno) << ")" << std::endl;
        gzclose(in_gz);
        return false;
    }
    gzbuffer(out_gz, GZ_BUFFER_SIZE);

    char buffer[GZ_BUFFER_SIZE];
    long long records_read = 0;
    long long records_kept = 0;
    long long records_dropped = 0;
    bool error_occurred = false;

    while (true) {
        if (gzgets(in_gz, buffer, sizeof(buffer)) == nullptr) {
            if (gzeof(in_gz)) break;
            int err; const char* msg = gzerror(in_gz, &err);
            if (err != Z_OK && err != Z_STREAM_END) {
                std::cerr << "Error reading R1 file during handshake filtering: " << (msg ? msg : "Unknown") << " (code " << err << ")" << std::endl;
                error_occurred = true;
            }
            break;
        }
        std::string id_line = buffer;

        std::string seq_line, plus_line, qual_line;
        if (gzgets(in_gz, buffer, sizeof(buffer)) == nullptr) { error_occurred = true; std::cerr << "Error: incomplete R1 record during handshake (seq missing)." << std::endl; break; }
        seq_line = buffer;
        if (gzgets(in_gz, buffer, sizeof(buffer)) == nullptr) { error_occurred = true; std::cerr << "Error: incomplete R1 record during handshake (plus missing)." << std::endl; break; }
        plus_line = buffer;
        if (gzgets(in_gz, buffer, sizeof(buffer)) == nullptr) { error_occurred = true; std::cerr << "Error: incomplete R1 record during handshake (qual missing)." << std::endl; break; }
        qual_line = buffer;

        records_read++;

        std::string id_no_nl = id_line;
        while (!id_no_nl.empty() && (id_no_nl.back() == '\n' || id_no_nl.back() == '\r')) id_no_nl.pop_back();

        if (id_no_nl.empty() || id_no_nl[0] != '@') {
            std::cerr << "Warning: malformed R1 record during handshake. Skipping." << std::endl;
            records_dropped++;
            continue;
        }

        std::string base_id = get_base_id_from_string(id_no_nl);
        if (r2_found_ids.count(base_id)) {
            if (gzputs(out_gz, id_line.c_str()) == -1 ||
                gzputs(out_gz, seq_line.c_str()) == -1 ||
                gzputs(out_gz, plus_line.c_str()) == -1 ||
                gzputs(out_gz, qual_line.c_str()) == -1)
            {
                int errnum = 0; const char* error_str = gzerror(out_gz, &errnum);
                std::cerr << "Error writing R1 record during handshake filtering. Error (" << errnum << "): " << (error_str ? error_str : "Unknown") << std::endl;
                error_occurred = true;
                break;
            }
            records_kept++;
        } else {
            records_dropped++;
        }
    }

    gzclose(in_gz);
    int close_status = gzclose(out_gz);
    if (close_status != Z_OK && close_status != Z_BUF_ERROR) {
        std::cerr << "Error: Failed to properly close R1 handshake output. Close status: " << close_status << std::endl;
        error_occurred = true;
    }

    std::cout << "Handshake summary: read=" << records_read
              << ", kept=" << records_kept
              << ", dropped (no R2 match)=" << records_dropped << std::endl;

    if (error_occurred) {
        std::error_code ec; fs::remove(r1_output_gz, ec);
        return false;
    }

    return true;
}

// ========================================================================
// == Main Function (MODIFIED for Chopping Arguments + Handshake)       ==
// ========================================================================
int main(int argc, char* argv[]) {
    // MODIFIED: Updated description to reflect simplified chopping
    cxxopts::Options options("fastq_filter_sort_extract_sortedR2",
        "Filters R1 FASTQ reads (Hamming/BaseComp), optionally chops passed R1 reads,\n"
        "outputs sorted, gzipped files for each R1 category (using parallel sort/write),\n"
        "extracts corresponding R2 reads PARALLELLY, and outputs them sorted and gzipped.\n"
        "NEW: Performs a two-way handshake to guarantee an exact 1:1 match between R1 and R2,\n"
        "even when R2 is truncated or missing records present in R1.");

    options.add_options()
        // --- Inputs / Outputs (original long names) ---
        ("r1_input",                       "Input R1 FASTQ file (.fastq.gz)",                                       cxxopts::value<std::string>())
        ("r2_input",                       "Input R2 FASTQ file (.fastq.gz)",                                       cxxopts::value<std::string>())
        ("r1_output_dir",                  "Output directory for sorted, filtered, and chopped R1 files",           cxxopts::value<std::string>())
        ("r2_output_dir",                  "Output directory for sorted, filtered R2 files",                        cxxopts::value<std::string>())
        ("hamming_filterout_dir",          "Output directory for sorted Hamming-rejected R1 files",                 cxxopts::value<std::string>())
        ("base_composition_filterout_dir", "Output directory for sorted Base Composition-rejected R1 files",        cxxopts::value<std::string>())

        // --- Hamming filter (original names) ---
        ("target_seq",                     "Target sequence for Hamming distance comparison",                       cxxopts::value<std::string>())
        ("seq_length",                     "Length of the sub-sequence region for Hamming comparison",              cxxopts::value<int>())
        ("start_pos",                      "Start position (0-based) of the Hamming sub-sequence in R1",            cxxopts::value<int>())
        ("threshold",                      "Maximum allowed Hamming distance",                                      cxxopts::value<int>())

        // --- Base composition filter (original names) ---
        ("base_composition_threshold",     "Minimum required base composition (ACGT fraction, 0.0-1.0)",            cxxopts::value<double>())
        ("start_pos_base_com",             "Start position (0-based) of the base composition region in R1",         cxxopts::value<int>())
        ("seq_length_base_com",            "Length of the base composition region in R1",                           cxxopts::value<int>())

        // --- Performance / runtime (original names; both - and _ accepted by cxxopts via aliases) ---
        ("threads",                        "Number of threads for filtering and R2 extraction stages",              cxxopts::value<int>()->default_value("1"))
        ("chunk_records",                  "Number of reads per chunk for processing",                              cxxopts::value<long long>()->default_value("100000"))
        ("temp_dir",                       "Base directory for temporary files (defaults to system temp or CWD)",   cxxopts::value<std::string>()->default_value(""))

        // NOTE: cxxopts treats option names with '-' literally; to support BOTH the original
        // hyphenated form (--sort-threads / --sort-mem-mb / --chop-reads / --chop-start /
        // --chop-length) AND the underscore form, we register the hyphenated form as the
        // primary key (matching your command line exactly).
        ("sort-mem-mb",                    "Approx memory limit (MB) per sort worker",                              cxxopts::value<size_t>()->default_value("1024"))
        ("sort-threads",                   "Number of parallel workers for the sort phase (Phase 1)",               cxxopts::value<int>()->default_value("0"))

        // --- Chopping (original hyphenated names) ---
        ("chop-reads",                     "Enable chopping of R1 reads that pass all filters",                     cxxopts::value<bool>()->default_value("false"))
        ("chop-start",                     "Start position (0-based) of the region to keep when chopping R1",       cxxopts::value<int>())
        ("chop-length",                    "Length of the region to keep when chopping R1",                         cxxopts::value<int>())

        ("h,help", "Print usage");

    Params params;
    bool enable_parallel_sort_alg_main = false;

    try {
        auto result = options.parse(argc, argv);

        if (result.count("help")) {
            std::cout << options.help() << std::endl;
            #ifdef USE_PARALLEL_SORT
                std::cout << "\nINFO: Compiled with USE_PARALLEL_SORT support." << std::endl;
            #else
                std::cout << "\nINFO: Compiled WITHOUT USE_PARALLEL_SORT support." << std::endl;
            #endif
            return 0;
        }

        // Required arguments check (using ORIGINAL argument names)
        std::vector<std::string> required_args = {
            "r1_input", "r2_input",
            "r1_output_dir", "r2_output_dir",
            "hamming_filterout_dir", "base_composition_filterout_dir",
            "target_seq", "seq_length", "start_pos", "threshold",
            "base_composition_threshold", "start_pos_base_com", "seq_length_base_com"
        };
        for (const auto& arg_name : required_args) {
            if (!result.count(arg_name)) {
                throw cxxopts::exceptions::exception("Required argument missing: --" + arg_name);
            }
        }

        params.r1_input                      = result["r1_input"].as<std::string>();
        params.r2_input                      = result["r2_input"].as<std::string>();
        params.r1_output_dir                 = result["r1_output_dir"].as<std::string>();
        params.r2_output_dir                 = result["r2_output_dir"].as<std::string>();
        params.hamming_filterout_dir         = result["hamming_filterout_dir"].as<std::string>();
        params.base_composition_filterout_dir= result["base_composition_filterout_dir"].as<std::string>();

        params.target_seq                    = result["target_seq"].as<std::string>();
        params.seq_length                    = result["seq_length"].as<int>();
        params.start_pos                     = result["start_pos"].as<int>();
        params.threshold                     = result["threshold"].as<int>();

        params.base_composition_threshold    = result["base_composition_threshold"].as<double>();
        params.base_comp_start_pos           = result["start_pos_base_com"].as<int>();
        params.base_comp_length              = result["seq_length_base_com"].as<int>();

        params.num_threads                   = result["threads"].as<int>();
        params.chunk_read_count              = result["chunk_records"].as<long long>();
        params.sort_memory_limit_mb          = result["sort-mem-mb"].as<size_t>();

        int sort_threads_arg                 = result["sort-threads"].as<int>();
        if (sort_threads_arg <= 0) {
            params.sort_threads = std::max(1, params.num_threads);
        } else {
            params.sort_threads = sort_threads_arg;
        }

        // Parse chopping options (original hyphenated names)
        params.chop_reads = result["chop-reads"].as<bool>();
        if (params.chop_reads) {
            if (!result.count("chop-start") || !result.count("chop-length")) {
                throw cxxopts::exceptions::exception("--chop-reads requires both --chop-start and --chop-length");
            }
            params.chop_start  = result["chop-start"].as<int>();
            params.chop_length = result["chop-length"].as<int>();
            if (params.chop_start < 0 || params.chop_length <= 0) {
                throw cxxopts::exceptions::exception("--chop-start must be >= 0 and --chop-length must be > 0");
            }
        }

        // Temp dir handling
        std::string temp_dir_arg = result["temp_dir"].as<std::string>();
        if (temp_dir_arg.empty()) {
            const char* sys_temp = std::getenv("TMPDIR");
            if (!sys_temp) sys_temp = std::getenv("TMP");
            if (!sys_temp) sys_temp = std::getenv("TEMP");
            if (sys_temp && fs::exists(sys_temp) && fs::is_directory(sys_temp)) {
                params.main_temp_dir_base = sys_temp;
            } else {
                params.main_temp_dir_base = fs::current_path();
            }
        } else {
            params.main_temp_dir_base = temp_dir_arg;
        }

        // Validate input files
        if (!fs::exists(params.r1_input)) throw std::runtime_error("R1 input file not found: " + params.r1_input.string());
        if (!fs::exists(params.r2_input)) throw std::runtime_error("R2 input file not found: " + params.r2_input.string());

        // Create output directories
        std::error_code ec;
        for (const auto& d : {params.r1_output_dir, params.r2_output_dir,
                              params.hamming_filterout_dir, params.base_composition_filterout_dir}) {
            fs::create_directories(d, ec);
            if (ec) throw std::runtime_error("Failed to create output directory: " + d.string() + " (" + ec.message() + ")");
        }

        // Validate parameters
        if (params.seq_length <= 0) throw std::runtime_error("--seq_length must be > 0");
        if (params.start_pos < 0)   throw std::runtime_error("--start_pos must be >= 0");
        if (params.threshold < 0)   throw std::runtime_error("--threshold must be >= 0");
        if ((int)params.target_seq.length() != params.seq_length) {
            throw std::runtime_error("Target sequence length (" + std::to_string(params.target_seq.length()) +
                                     ") does not match --seq_length (" + std::to_string(params.seq_length) + ")");
        }
        if (params.base_composition_threshold < 0.0 || params.base_composition_threshold > 1.0) {
            throw std::runtime_error("--base_composition_threshold must be in [0.0, 1.0]");
        }
        if (params.base_comp_length <= 0)   throw std::runtime_error("--seq_length_base_com must be > 0");
        if (params.base_comp_start_pos < 0) throw std::runtime_error("--start_pos_base_com must be >= 0");
        if (params.num_threads <= 0)        params.num_threads = 1;
        if (params.chunk_read_count <= 0)   params.chunk_read_count = 100000;
        if (params.sort_memory_limit_mb == 0) params.sort_memory_limit_mb = 1024;

        #ifdef USE_PARALLEL_SORT
            enable_parallel_sort_alg_main = (params.sort_threads > 0);
        #else
            enable_parallel_sort_alg_main = false;
        #endif

    } catch (const cxxopts::exceptions::exception& e) {
        std::cerr << "Argument parsing error: " << e.what() << std::endl;
        std::cerr << options.help() << std::endl;
        return 1;
    } catch (const std::exception& e) {
        std::cerr << "Configuration error: " << e.what() << std::endl;
        return 1;
    }

    // Print configuration
    std::cout << "==================== Configuration ====================" << std::endl;
    std::cout << "R1 input:                  " << params.r1_input << std::endl;
    std::cout << "R2 input:                  " << params.r2_input << std::endl;
    std::cout << "R1 output dir:             " << params.r1_output_dir << std::endl;
    std::cout << "R2 output dir:             " << params.r2_output_dir << std::endl;
    std::cout << "Hamming reject dir:        " << params.hamming_filterout_dir << std::endl;
    std::cout << "BaseComp reject dir:       " << params.base_composition_filterout_dir << std::endl;
    std::cout << "Target seq:                " << params.target_seq << std::endl;
    std::cout << "Hamming region:            start=" << params.start_pos << ", length=" << params.seq_length << ", threshold=" << params.threshold << std::endl;
    std::cout << "BaseComp region:           start=" << params.base_comp_start_pos << ", length=" << params.base_comp_length << ", threshold>=" << params.base_composition_threshold << std::endl;
    std::cout << "Chop reads:                " << (params.chop_reads ? "ENABLED" : "disabled") << std::endl;
    if (params.chop_reads) {
        std::cout << "  Chop region:             start=" << params.chop_start << ", length=" << params.chop_length << std::endl;
    }
    std::cout << "Threads (filter/extract):  " << params.num_threads << std::endl;
    std::cout << "Sort workers:              " << params.sort_threads << std::endl;
    std::cout << "Chunk size:                " << params.chunk_read_count << std::endl;
    std::cout << "Sort mem/worker:           " << params.sort_memory_limit_mb << " MB" << std::endl;
    std::cout << "Temp dir base:             " << params.main_temp_dir_base << std::endl;
    std::cout << "Two-way handshake:         ENABLED (guarantees R1/R2 1:1 match)" << std::endl;
    std::cout << "=======================================================" << std::endl;

    auto wall_start = std::chrono::steady_clock::now();

    // ===== Stage 1: Filtering =====
    std::unique_ptr<TempDir> main_temp_dir_ptr;
    fs::path main_temp_path;
    std::vector<fs::path> passed_chunks, hamming_chunks, basecomp_chunks;

    try {
        main_temp_dir_ptr = std::make_unique<TempDir>(params.main_temp_dir_base, "filter_sort_main_");
        main_temp_path = main_temp_dir_ptr->getPath();
        std::cout << "\nMain temp directory: " << main_temp_path << std::endl;

        std::cout << "\n--- Stage 1: Filtering (and optional chopping) R1 ---" << std::endl;
        igzstream r1_in(params.r1_input.c_str());
        if (!r1_in.good()) throw std::runtime_error("Failed to open R1 input: " + params.r1_input.string());

        std::vector<std::future<UnsortedChunkOutputPaths>> filter_futures;
        std::vector<std::string> current_chunk;
        current_chunk.reserve(params.chunk_read_count * 4);
        long long records_read = 0;
        int chunk_idx = 0;

        std::string line_id, line_seq, line_plus, line_qual;
        while (std::getline(r1_in, line_id)) {
            if (!std::getline(r1_in, line_seq)) handle_read_error(r1_in, params.r1_input, records_read, "sequence");
            if (!std::getline(r1_in, line_plus)) handle_read_error(r1_in, params.r1_input, records_read, "plus");
            if (!std::getline(r1_in, line_qual)) handle_read_error(r1_in, params.r1_input, records_read, "quality");

            current_chunk.push_back(std::move(line_id));
            current_chunk.push_back(std::move(line_seq));
            current_chunk.push_back(std::move(line_plus));
            current_chunk.push_back(std::move(line_qual));
            records_read++;

            if ((long long)(current_chunk.size() / 4) >= params.chunk_read_count) {
                filter_futures.push_back(std::async(std::launch::async, process_filter_chunk,
                                                    std::move(current_chunk), main_temp_path, chunk_idx++, std::cref(params)));
                current_chunk = {};
                current_chunk.reserve(params.chunk_read_count * 4);

                while (filter_futures.size() >= static_cast<size_t>(params.num_threads * 2 + 4)) {
                    bool found_ready = false;
                    for (auto it = filter_futures.begin(); it != filter_futures.end();) {
                        if (it->valid() && it->wait_for(std::chrono::milliseconds(0)) == std::future_status::ready) {
                            auto res = it->get();
                            if (res.success) {
                                passed_chunks.push_back(res.passed);
                                hamming_chunks.push_back(res.hamming_rejected);
                                basecomp_chunks.push_back(res.basecomp_rejected);
                            } else {
                                throw std::runtime_error("Filter chunk processing failed.");
                            }
                            it = filter_futures.erase(it);
                            found_ready = true;
                        } else {
                            ++it;
                        }
                    }
                    if (!found_ready) std::this_thread::sleep_for(std::chrono::milliseconds(20));
                }
            }
        }

        if (!current_chunk.empty()) {
            filter_futures.push_back(std::async(std::launch::async, process_filter_chunk,
                                                std::move(current_chunk), main_temp_path, chunk_idx++, std::cref(params)));
        }

        std::cout << "Total R1 records read: " << records_read << ". Waiting for filter tasks..." << std::flush;
        for (auto& fut : filter_futures) {
            auto res = fut.get();
            if (res.success) {
                passed_chunks.push_back(res.passed);
                hamming_chunks.push_back(res.hamming_rejected);
                basecomp_chunks.push_back(res.basecomp_rejected);
            } else {
                throw std::runtime_error("Filter chunk processing failed.");
            }
        }
        std::cout << " Done." << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Stage 1 (Filtering) error: " << e.what() << std::endl;
        return 2;
    }

    // ===== Stage 2: Sort R1 categories =====
    size_t sort_mem_bytes = params.sort_memory_limit_mb * 1024ULL * 1024ULL;

    auto get_out_path = [](const fs::path& dir, const fs::path& orig_input, const std::string& suffix) -> fs::path {
        std::string base = orig_input.filename().string();
        size_t gz = base.rfind(".gz"); if (gz != std::string::npos) base = base.substr(0, gz);
        size_t fq = base.rfind(".fastq"); if (fq != std::string::npos) base = base.substr(0, fq);
        else { size_t fq2 = base.rfind(".fq"); if (fq2 != std::string::npos) base = base.substr(0, fq2); }
        return dir / (base + suffix + ".fastq.gz");
    };

    fs::path r1_passed_out = get_out_path(params.r1_output_dir, params.r1_input, "_passed");
    fs::path r1_hamming_out = get_out_path(params.hamming_filterout_dir, params.r1_input, "_hamming_rejected");
    fs::path r1_basecomp_out = get_out_path(params.base_composition_filterout_dir, params.r1_input, "_basecomp_rejected");
    fs::path r2_sorted_out = get_out_path(params.r2_output_dir, params.r2_input, "_passed");

    try {
        if (!sort_and_merge_chunk_category(passed_chunks, r1_passed_out, main_temp_path, "r1_passed",
                                            sort_mem_bytes, params.sort_threads, enable_parallel_sort_alg_main))
            throw std::runtime_error("Failed to sort R1 passed category.");

        if (!sort_and_merge_chunk_category(hamming_chunks, r1_hamming_out, main_temp_path, "r1_hamming",
                                            sort_mem_bytes, params.sort_threads, enable_parallel_sort_alg_main))
            throw std::runtime_error("Failed to sort R1 hamming category.");

        if (!sort_and_merge_chunk_category(basecomp_chunks, r1_basecomp_out, main_temp_path, "r1_basecomp",
                                            sort_mem_bytes, params.sort_threads, enable_parallel_sort_alg_main))
            throw std::runtime_error("Failed to sort R1 basecomp category.");
    } catch (const std::exception& e) {
        std::cerr << "Stage 2 (Sort R1) error: " << e.what() << std::endl;
        return 3;
    }

    // ===== Stage 3: R2 extraction with two-way handshake ID collection =====
    std::unordered_set<std::string> r2_found_ids;
    fs::path r2_intermediate;

    try {
        r2_intermediate = extract_r2_reads_parallel(r1_passed_out, params.r2_input, r2_sorted_out,
                                                    main_temp_path, params, r2_found_ids);
    } catch (const std::exception& e) {
        std::cerr << "Stage 3 (R2 extraction) error: " << e.what() << std::endl;
        return 4;
    }

    // ===== Stage 4: Two-Way Handshake - filter R1 passed by R2 found IDs =====
    // This guarantees an exact 1:1 match between R1 and R2, even when R2 is truncated
    // or missing records that exist in R1.
    if (!r2_intermediate.empty() && !r2_found_ids.empty()) {
        std::cout << "\n========== Two-Way Handshake Stage ==========" << std::endl;

        // Check whether handshake will actually change anything. If R2 contained every R1 ID,
        // there's nothing to filter — but we still run it for safety/consistency.
        fs::path r1_handshake_temp = r1_passed_out;
        r1_handshake_temp += ".handshake.tmp.gz";

        if (!filter_r1_by_r2_ids(r1_passed_out, r1_handshake_temp, r2_found_ids)) {
            std::cerr << "Error: Two-way handshake R1 filtering failed. R1 passed file kept as-is, R1/R2 may not match 1:1." << std::endl;
            std::error_code ec; fs::remove(r1_handshake_temp, ec);
            return 5;
        }

        // Atomically replace R1 passed file with handshake-filtered version
        std::error_code ec;
        fs::remove(r1_passed_out, ec);
        fs::rename(r1_handshake_temp, r1_passed_out, ec);
        if (ec) {
            std::cerr << "Error: Failed to replace R1 passed file with handshake-filtered version: " << ec.message() << std::endl;
            return 5;
        }
        std::cout << "Two-way handshake complete. R1 passed file now contains only reads with matching R2." << std::endl;
    } else if (r2_intermediate.empty()) {
        std::cout << "\nNote: No R2 reads were extracted. Skipping two-way handshake stage." << std::endl;
        std::cout << "      R1 passed file remains untouched (may contain reads without R2 mates)." << std::endl;
    }

    // ===== Stage 5: Sort R2 extracted file =====
    if (!r2_intermediate.empty()) {
        try {
            std::vector<fs::path> r2_singleton = { r2_intermediate };
            if (!sort_and_merge_chunk_category(r2_singleton, r2_sorted_out, main_temp_path, "r2_extracted",
                                                sort_mem_bytes, params.sort_threads, enable_parallel_sort_alg_main))
                throw std::runtime_error("Failed to sort R2 extracted file.");

            // Cleanup intermediate
            std::error_code ec;
            if (fs::exists(r2_intermediate, ec)) fs::remove(r2_intermediate, ec);
        } catch (const std::exception& e) {
            std::cerr << "Stage 5 (Sort R2) error: " << e.what() << std::endl;
            return 6;
        }
    } else {
        // Create empty R2 output for consistency
        gzFile empty_out = gzopen(r2_sorted_out.c_str(), "wb");
        if (empty_out) gzclose(empty_out);
        std::cout << "Created empty R2 output file: " << r2_sorted_out << std::endl;
    }

    auto wall_end = std::chrono::steady_clock::now();
    auto elapsed_sec = std::chrono::duration_cast<std::chrono::seconds>(wall_end - wall_start).count();

    std::cout << "\n==================== Summary ====================" << std::endl;
    std::cout << "R1 passed (handshake-verified): " << r1_passed_out << std::endl;
    std::cout << "R1 Hamming-rejected:            " << r1_hamming_out << std::endl;
    std::cout << "R1 BaseComp-rejected:           " << r1_basecomp_out << std::endl;
    std::cout << "R2 passed (handshake-verified): " << r2_sorted_out << std::endl;
    std::cout << "Total wall time:                " << elapsed_sec << " seconds" << std::endl;
    std::cout << "=================================================" << std::endl;
    try {
        fs::path summary_path = params.r2_output_dir / "processing_summary.txt";
        std::ofstream summary(summary_path);
        if (!summary.is_open()) {
            std::cerr << "Warning: Failed to create processing_summary.txt at " << summary_path << std::endl;
        } else {
            // Timestamp
            std::time_t now_t = std::time(nullptr);
            char timebuf[64];
            std::strftime(timebuf, sizeof(timebuf), "%Y-%m-%d %H:%M:%S", std::localtime(&now_t));

            // Helper: count reads in a gzipped FASTQ (lines / 4)
            auto count_reads_gz = [](const fs::path& p) -> long long {
                if (!fs::exists(p)) return -1;
                igzstream in(p.c_str());
                if (!in.good()) return -1;
                long long lines = 0;
                std::string l;
                while (std::getline(in, l)) ++lines;
                return lines / 4;
            };

            long long n_r1_passed   = count_reads_gz(r1_passed_out);
            long long n_r1_hamming  = count_reads_gz(r1_hamming_out);
            long long n_r1_basecomp = count_reads_gz(r1_basecomp_out);
            long long n_r2_passed   = count_reads_gz(r2_sorted_out);

            summary << "============================================================\n";
            summary << "  FASTQ Filter / Sort / Pair Processing Summary\n";
            summary << "============================================================\n";
            summary << "Generated:                 " << timebuf << "\n";
            summary << "Total wall time (seconds): " << elapsed_sec << "\n";
            summary << "\n";
            summary << "------------------- Input Files -------------------\n";
            summary << "R1 input:                  " << params.r1_input  << "\n";
            summary << "R2 input:                  " << params.r2_input  << "\n";
            summary << "\n";
            summary << "------------------- Output Files ------------------\n";
            summary << "R1 passed (handshake):     " << r1_passed_out    << "\n";
            summary << "R1 Hamming-rejected:       " << r1_hamming_out   << "\n";
            summary << "R1 BaseComp-rejected:      " << r1_basecomp_out  << "\n";
            summary << "R2 passed (handshake):     " << r2_sorted_out    << "\n";
            summary << "\n";
            summary << "------------------- Filter Parameters -------------\n";
            summary << "Target sequence:           " << params.target_seq << "\n";
            summary << "Hamming start position:    " << params.start_pos  << "\n";
            summary << "Hamming region length:     " << params.seq_length << "\n";
            summary << "Hamming threshold (max):   " << params.threshold  << "\n";
            summary << "BaseComp start position:   " << params.base_comp_start_pos << "\n";
            summary << "BaseComp region length:    " << params.base_comp_length    << "\n";
            summary << "BaseComp threshold (min):  " << params.base_composition_threshold << "\n";
            summary << "\n";
            summary << "------------------- Chopping ----------------------\n";
            summary << "Chop reads enabled:        " << (params.chop_reads ? "YES" : "no") << "\n";
            if (params.chop_reads) {
                summary << "Chop start position:       " << params.chop_start  << "\n";
                summary << "Chop length:               " << params.chop_length << "\n";
            }
            summary << "\n";
            summary << "------------------- Runtime Settings --------------\n";
            summary << "Filter/extract threads:    " << params.num_threads        << "\n";
            summary << "Sort workers:              " << params.sort_threads       << "\n";
            summary << "Chunk size (reads):        " << params.chunk_read_count   << "\n";
            summary << "Sort memory/worker (MB):   " << params.sort_memory_limit_mb << "\n";
            summary << "Temp directory base:       " << params.main_temp_dir_base << "\n";
            summary << "Two-way handshake:         ENABLED (R1/R2 1:1 guaranteed)\n";
            summary << "\n";
            summary << "------------------- Read Counts -------------------\n";
            summary << "R1 passed reads:           " << n_r1_passed   << "\n";
            summary << "R1 Hamming-rejected reads: " << n_r1_hamming  << "\n";
            summary << "R1 BaseComp-rejected reads:" << n_r1_basecomp << "\n";
            summary << "R2 passed reads:           " << n_r2_passed   << "\n";
            if (n_r1_passed > 0 && n_r2_passed >= 0) {
                summary << "R1/R2 pair match:          "
                        << (n_r1_passed == n_r2_passed ? "PERFECT 1:1" : "MISMATCH")
                        << " (" << n_r2_passed << " / " << n_r1_passed << ")\n";
            }
            summary << "============================================================\n";
            summary.close();

            std::cout << "Processing summary written to: " << summary_path << std::endl;
        }
    } catch (const std::exception& e) {
        std::cerr << "Warning: Could not write processing_summary.txt: " << e.what() << std::endl;
    }

    return 0;
}