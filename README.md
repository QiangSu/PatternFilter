# Fastq Pattern Filter and Sorter

This C++17 program filters paired-end FASTQ reads based on Hamming distance to a target sequence in R1 and base composition criteria in R1. It then sorts the R1 reads into 'passed', 'Hamming rejected', and 'base composition rejected' categories. Corresponding R2 reads for the 'passed' R1 set are extracted and sorted.

## Dependencies

*   A C++17 compliant compiler (e.g., g++ 7 or later, or Clang 5 or later)
*   Zlib development library (e.g., `zlib1g-dev` on Debian/Ubuntu, `zlib-devel` on CentOS/Fedora)
*   pthreads (usually included with the compiler toolchain on POSIX systems)

Included header-only/source libraries:
*   `cxxopts.hpp` (for command-line argument parsing) - in `include/`
*   `gzstream` (for reading/writing gzipped files) - in `gzstream/`

## Setting up the Build Environment (Installing Dependencies)

Before you can compile the program, you need to ensure the required development tools and libraries are installed on your system.

**1. C++17 Compiler:**
   You need a C++ compiler that supports the C++17 standard. GCC (g++) version 7 or newer, or Clang version 5 or newer, are common choices.

   *   **On Debian/Ubuntu based systems (like Ubuntu):**
       ```bash
       sudo apt update
       sudo apt install build-essential g++
       # Check version:
       g++ --version
       ```
       `build-essential` installs a lot of common development tools, including `g++` and `make`.

   *   **On Fedora/CentOS/RHEL based systems:**
       ```bash
       sudo dnf groupinstall "Development Tools" # For Fedora
       # or
       sudo yum groupinstall "Development Tools" # For CentOS/RHEL
       sudo dnf install gcc-c++ # For Fedora (if not covered by groupinstall)
       # or
       sudo yum install gcc-c++ # For CentOS/RHEL (if not covered by groupinstall)
       # Check version:
       g++ --version
       ```

   *   **On macOS:**
       Install Xcode Command Line Tools (which includes Clang):
       ```bash
       xcode-select --install
       # Check version:
       clang++ --version
       ```
       You might need to set your compiler to `clang++` in the compilation command if `g++` defaults to an older version or is an alias to `clang++`.

**2. Zlib Development Library:**
   This library is required for working with gzipped files (.gz). You need the development package, which includes header files.

   *   **On Debian/Ubuntu based systems:**
       ```bash
       sudo apt install zlib1g-dev
       ```

   *   **On Fedora/CentOS/RHEL based systems:**
       ```bash
       sudo dnf install zlib-devel # For Fedora
       # or
       sudo yum install zlib-devel # For CentOS/RHEL
       ```

   *   **On macOS:**
       Zlib is usually available if Xcode Command Line Tools are installed. If not, or for a specific version, you can use Homebrew:
       ```bash
       brew install zlib
       ```
       (Note: Homebrew might install it in a non-standard path, so you might need to adjust include/library paths during compilation if the system doesn't find it automatically.)

**3. pthreads:**
   The POSIX Threads (pthreads) library is used for `std::thread` and `std::future`.
   *   On **Linux and macOS**, pthreads is typically part of the standard C library and the compiler toolchain (e.g., GCC, Clang). You generally **do not need to install it separately**. The `-pthread` flag used during compilation tells the linker to link against the pthreads library.

**4. Included Libraries (`cxxopts.hpp` and `gzstream`):**
   The libraries `cxxopts.hpp` and `gzstream` (which consists of `gzstream.h` and `gzstream.C`) are **included directly in this repository** within the `include/` and `gzstream/` directories, respectively.
   *   **You do not need to download or install them separately** if you have cloned this repository. They are compiled/included along with `Pattern_Filter.cpp` as specified in the compilation command.
   *   Original sources (for reference):
       *   cxxopts: [https://github.com/jarro2783/cxxopts](https://github.com/jarro2783/cxxopts)
       *   gzstream: [http://www.cs.unc.edu/Research/compgeom/gzstream/](http://www.cs.unc.edu/Research/compgeom/gzstream/)

**5. Clone the Repository (if you haven't already):**
   First, get the code and the included libraries onto your local machine:
   ```bash
   git clone https://github.com/QiangSu/PatternFilter.git
   cd FastqPatternFilter
## Compilation

Once the dependencies are installed, navigate to the project root directory and run:

```bash
g++ -std=c++17 -O3 -Wall -I./include -I./gzstream ./Pattern_Filter.cpp ./gzstream/gzstream.C -o Pattern_Filter -pthread -lz -DUSE_PARALLEL_SORT.


##  Usage
Run the program with --help to see all available options:

./Pattern_Filter --help

An example command using the test data (make sure test_data/ exists from the clone):
mkdir -p R1_passed R2_extracted Hamming_rejected Basecomp_rejected Test_Temp_Dir

./Pattern_Filter \
  --r1_input test_data/extracted_100000_reads_musWT_STR_R1_001.fastq.gz \
  --r2_input test_data/extracted_100000_reads_musWT_STR_R2_001.fastq.gz \
  --r1_output_dir R1_passed \
  --r2_output_dir R2_extracted \
  --hamming_filterout_dir Hamming_rejected \
  --base_composition_filterout_dir Basecomp_rejected \
  --target_seq "AGCTAGCT" \
  --seq_length 8 \
  --start_pos 10 \
  --threshold 2 \
  --start_pos_base_com 20 \
  --seq_length_base_com 30 \
  --base_composition_threshold 0.8 \
  --threads 2 \
  --sort-threads 2 \
  --chunk_records 10000 \
  --sort-mem-mb 128 \
  --temp_dir Test_Temp_Dir \
  --parallel-sort-alg true






