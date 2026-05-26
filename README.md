# Pattern-Filter: High-Performance Pre-processing of sc/snRNA-seq Data

Pattern-Filter is a high-performance C++17 command-line tool for filtering, sorting, and pre-processing paired-end FASTQ data from single-cell and single-nucleus RNA-sequencing experiments. It allows for precise quality control at the raw read level by validating expected molecular structures based on user-defined patterns, such as poly(T) tails, linkers, or template switching oligos (TSOs).

## Key Features

*   **Multi-criteria Filtering:** Filters R1 reads based on Hamming distance to a target sequence and/or base composition (e.g., % of A/C/G/T).
*   **Dynamic Read Trimming (--chop-reads):** Automatically crops validated R1 sequences and quality scores to a precise, user-defined length (e.g., exactly 28 bp for 10x v3) to ensure seamless compatibility with downstream aligners like Cell Ranger or STARsolo.
*   **Strict R1/R2 Synchronization:** Features a built-in two-way handshake validation to guarantee a perfect 1:1 match between output R1 and R2 reads, preventing downstream pipeline crashes caused by truncated or corrupted FASTQ files.
*   **Efficient & Parallelized:** Built for speed and memory efficiency, leveraging multi-threading for all major processing stages (filtering, sorting, and R2 extraction).
*   **Workflow-Agnostic:** Flexible enough to be configured for a wide variety of sc/snRNA-seq technologies (10x Genomics, Drop-seq, SPLiT-seq, etc.).
*   **Comprehensive Reporting:** Automatically generates a detailed processing_summary.txt logging exact read counts, filter outcomes, and parameters for complete reproducibility.

## Core Concepts

Pattern-Filter operates in several sequential phases:

1.  **Parallel Filtering (R1):** Reads from the R1 file are processed in chunks across multiple threads. Each read is checked against two main criteria:
    *   **Hamming Filter:** Compares a specific region of the read (defined by `--start_pos` and `--seq_length`) to a `--target_seq`. If the number of mismatches exceeds the `--threshold`, the read is rejected. This is ideal for finding fixed sequences like poly(T) tails or linkers.
    *   **Base Composition Filter:** Calculates the fraction of 'A', 'C', 'G', 'T' bases in a specified region (defined by `--start_pos_base_com` and `--seq_length_base_com`). If the fraction is below the `--base_composition_threshold`, the read is rejected. This is useful for validating the integrity of barcodes or UMIs.
2.  **In-Line Read Chopping (Optional):** If a read passes the filters and the --chop-reads flag is active, the tool precisely crops the read's sequence and quality scores based on --chop-start and --chop-length before saving it. Reads that are physically too short to be chopped are safely discarded to prevent errors.
3.  **Parallel Sorting (R1):** Reads from each category ('passed', 'Hamming rejected', 'base-comp rejected') are sorted by read ID in parallel. This is an external memory sort, capable of handling files larger than system RAM.

4.  **Parallel R2 Extraction & Handshake Validation:** For reads that passed the R1 filters, the tool finds their corresponding pairs in the original R2 file. A strict ID-matching handshake ensures no orphan reads slip through, outputting perfectly synchronized and sorted R1 and R2 files ready for alignment.

---

## 1. System Requirements & Installation

### System Requirements

*   **Operating System:** A modern Linux distribution (e.g., Ubuntu 18.04+, CentOS 7+) or macOS.
*   **Compiler:** A C++17 compliant compiler.
    *   **GCC:** version 7 or newer.
    *   **Clang:** version 5 or newer.
*   **Dependencies:**
    *   **Zlib (development library):** For handling `.gz` files.
    *   **pthreads:** Standard on Linux/macOS; required for `std::thread`.
    *   **Intel TBB (optional but recommended):** For the C++17 parallel sort algorithm (`std::execution::par`). Provides a significant speed boost.

### Installing Dependencies

#### On Debian/Ubuntu:
```bash
# Update package list
sudo apt update

# Install compiler and zlib
sudo apt install build-essential g++ zlib1g-dev

# (Optional, Recommended) Install Intel TBB
sudo apt install libtbb-dev

# For older systems like Ubuntu 16.04, you need a newer compiler:
# sudo add-apt-repository ppa:ubuntu-toolchain-r/test
# sudo apt update
# sudo apt install g++-7 zlib1g-dev libtbb-dev
# (You may need to use `g++-7` in the compilation command below)
```

#### On CentOS/Fedora/RHEL:
```bash
# Install development tools, compiler, and zlib
sudo dnf groupinstall "Development Tools"
sudo dnf install gcc-c++ zlib-devel

# (Optional, Recommended) Install Intel TBB
sudo dnf install tbb-devel
```

#### On macOS:
```bash
# Install Xcode Command Line Tools (includes Clang and zlib)
xcode-select --install

# Install dependencies with Homebrew
brew install tbb
```

---

## 2. Compilation

First, clone the repository:
```bash
git clone https://github.com/QiangSu/PatternFilter.git
cd PatternFilter
```

Compile the program using the following command. This command enables the high-performance parallel sorting algorithm (requires Intel TBB) and points to the included gzstream library:

```bash
g++ -std=c++17 -O3 -D USE_PARALLEL_SORT -o Pattern_Filter Pattern_Filter.cpp ./gzstream/gzstream.C -lz -pthread -ltbb -I ./gzstream
```

---

## 3. Usage

### Parameter Recipes for Common Workflows

The key to using Pattern-Filter is correctly specifying the regions to filter. Here are recipes for common scRNA-seq platforms. **Note:** Positions are 0-based.

| Workflow              | Read to Filter | Goal                               | `target_seq`        | `start_pos` | `seq_length` | `threshold` | Notes|
| --------------------- | -------------- | ---------------------------------- | ------------------- | ----------- | ------------ | ----------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| **10x Genomics 3'**   | R1             | Find poly(T) tail after UMI+BC     | `"TTTTTTTTTTTTTTT"` | `28`        | `15`         | `2`         | For v3 chemistry, the CB+UMI is 28 bp. We check for a poly(T) starting at position 28. Use `--base_composition_threshold 1 --start_pos_base_com 0 --seq_length_base_com 28` to ensure barcode/UMI integrity. |
| **10x Genomics 5'**   | R1             | Find TSO sequence after UMI+BC     | `"TTTCTTATATGGG"`   | `26`        | `13`         | `2`         | For v2 chemistry, the CB+UMI is 26 bp. This filter validates the Template Switching Oligo.                                                                                    |
| **Drop-seq**          | R1             | Find poly(T) tail after UMI+BC     | `"TTTTTTTTTTTTTTT"` | `20`        | `15`         | `2`         | The Cell Barcode (12 bp) + UMI (8 bp) is 20 bp. The filter looks for the poly(T) capture sequence immediately after.                                                         |
| **SPLiT-seq**         | R2             | Find Round 2 & 3 linkers           | *[See Note]*       | *[See Note]* | `30`         | `2`         | SPLiT-seq requires multi-step filtering on **R2**. First, filter on the R2 linker, then on the R3 linker. You must run the tool sequentially. See example command in the manuscript. |

### Example Command (10x Genomics 3' v3)
This command validates the 28 bp Barcode+UMI and filters for the poly(T) tail starting at position 28. It then chops the R1 read down to exactly 28 bp for seamless compatibility with downstream tools like Cell Ranger or STARsolo.

```bash
# First, create output directories
mkdir -p R1_passed R2_extracted Hamming_rejected Basecomp_rejected Temp_Dir

# Run Pattern_Filter
./Pattern_Filter \
  --r1_input /path/to/your_R1.fastq.gz \
  --r2_input /path/to/your_R2.fastq.gz \
  --r1_output_dir /path/to/output_R1_passed \
  --r2_output_dir /path/to/output_R2_extracted \
  --hamming_filterout_dir /path/to/R1_hamming_rejects \
  --base_composition_filterout_dir /path/to/R1_basecomp_rejects \
  --target_seq TTTTTTTTTTTTTTT \
  --threshold 2 \
  --start_pos 28 \
  --seq_length 15 \
  --base_composition_threshold 1.0 \
  --start_pos_base_com 0 \
  --seq_length_base_com 28 \
  --chop-reads \
  --chop-start 0 \
  --chop-length 28 \
  --chunk_records 10000000 \
  --threads 40 \
  --sort-threads 40 \
  --sort-mem-mb 204800
```

### 10x Genomics 5' (v2 Chemistry)
This command validates the 26 bp Barcode+UMI and filters for the Template Switching Oligo (TSO) sequence. R1 reads are chopped to 26 bp.

```bash
# First, create output directories
mkdir -p R1_passed R2_extracted Hamming_rejected Basecomp_rejected Temp_Dir

# Run Pattern_Filter
./Pattern_Filter \
  --r1_input /path/to/your_R1.fastq.gz \
  --r2_input /path/to/your_R2.fastq.gz \
  --r1_output_dir /path/to/output_R1_passed \
  --r2_output_dir /path/to/output_R2_extracted \
  --hamming_filterout_dir /path/to/R1_hamming_rejects \
  --base_composition_filterout_dir /path/to/R1_basecomp_rejects \
  --target_seq TTTCTTATATGGG \
  --threshold 2 \
  --start_pos 26 \
  --seq_length 13 \
  --base_composition_threshold 1.0 \
  --start_pos_base_com 0 \
  --seq_length_base_com 26 \
  --chop-reads \
  --chop-start 0 \
  --chop-length 26 \
  --chunk_records 10000000 \
  --threads 40 \
  --sort-threads 40 \
  --sort-mem-mb 204800
```

### Drop-seq
This command validates the 20 bp Barcode+UMI (12 bp cell barcode + 8 bp UMI) and filters for the poly(T) capture sequence immediately following it. R1 reads are chopped to 20 bp.

```bash
# First, create output directories
mkdir -p R1_passed R2_extracted Hamming_rejected Basecomp_rejected Temp_Dir

# Run Pattern_Filter
./Pattern_Filter \
  --r1_input /path/to/your_R1.fastq.gz \
  --r2_input /path/to/your_R2.fastq.gz \
  --r1_output_dir /path/to/output_R1_passed \
  --r2_output_dir /path/to/output_R2_extracted \
  --hamming_filterout_dir /path/to/R1_hamming_rejects \
  --base_composition_filterout_dir /path/to/R1_basecomp_rejects \
  --target_seq TTTTTTTTTTTTTTT \
  --threshold 2 \
  --start_pos 20 \
  --seq_length 15 \
  --base_composition_threshold 1.0 \
  --start_pos_base_com 0 \
  --seq_length_base_com 20 \
  --chop-reads \
  --chop-start 0 \
  --chop-length 20 \
  --chunk_records 10000000 \
  --threads 40 \
  --sort-threads 40 \
  --sort-mem-mb 204800
```

### SPLiT-seq (Multi-step filtering)
Because the structural linkers in SPLiT-seq are located on R2 (and differ by kit version), you must pass your R2 file to the --r1_input argument and your R1 file to the --r2_input argument to filter based on R2 structure. 

SPLiT-seq requires two sequential runs. First, filter for the Round 2 Linker. Then, take the passed files from Run 1 and filter them for the Round 3 Linker. (Note: Replace [YOUR_R2_LINKER] and [YOUR_R3_LINKER] with the actual linker sequences used in your specific SPLiT-seq kit, and adjust start_pos accordingly).

```bash
# First, create output directories
mkdir -p R1_passed R2_extracted Hamming_rejected Basecomp_rejected Temp_Dir
```

# Run 1: Filter Round 2 Linker
```bash
./Pattern_Filter \
  --r1_input /path/to/your_R2.fastq.gz \
  --r2_input /path/to/your_R1.fastq.gz \
  --r1_output_dir /path/to/Step1_Passed_R2 \
  --r2_output_dir /path/to/Step1_Passed_R1 \
  --hamming_filterout_dir /path/to/Step1_R2_hamming_rejects \
  --base_composition_filterout_dir /path/to/Step1_R2_basecomp_rejects \
  --target_seq [YOUR_R2_LINKER] \
  --threshold 2 \
  --start_pos 10 \
  --seq_length 30 \
  --chunk_records 10000000 \
  --threads 40 \
  --sort-threads 40 \
  --sort-mem-mb 204800
```

# Run 2: Filter Round 3 Linker
```bash
./Pattern_Filter \
  --r1_input /path/to/Step1_Passed_R2/filtered_R1.fastq.gz \
  --r2_input /path/to/Step1_Passed_R1/filtered_R2.fastq.gz \
  --r1_output_dir /path/to/Final_Passed_R2 \
  --r2_output_dir /path/to/Final_Passed_R1 \
  --hamming_filterout_dir /path/to/Step2_R2_hamming_rejects \
  --base_composition_filterout_dir /path/to/Step2_R2_basecomp_rejects \
  --target_seq [YOUR_R3_LINKER] \
  --threshold 2 \
  --start_pos 48 \
  --seq_length 30 \
  --chunk_records 10000000 \
  --threads 40 \
  --sort-threads 40 \
  --sort-mem-mb 204800
```

---

## 4. Troubleshooting

**Problem: Almost all my reads are rejected!**

This is the most common issue and usually has one of three causes:

1.  **Incorrect Filter Parameters:** Your `start_pos` and `seq_length` might be wrong for your data's read structure. Double-check the "Parameter Recipes" table above and verify the read structure of your specific chemistry version.

2.  **R1 Reads are Too Short for the Filter:** The program checks if a read is long enough to contain the filter region. For example, if you set --start_pos 28 --seq_length 15, the R1 read must be at least 43 bp long. If your reads are shorter, they will be rejected. Look for Reads Rejected (Too Short for Filter) in your processing_summary.txt.

3.  **Reads are Too Short to be Chopped:** If you enable --chop-reads --chop-length 28 but your raw sequencing read is only 26 bp long, the tool will safely discard it into the basecomp_rejected bin to avoid creating corrupted fastq files.

    *   **Solution:** If this happens, it is not a software failure, but a mismatch between your sequencing data structure and your parameters. You may need to disable the Hamming filter (--threshold 999) or disable read chopping if your data naturally fits downstream tools. Check the processing_summary.txt generated in your output directory for an exact breakdown of why reads failed.