# Pattern-Filter filtering sc/snRNA-seq data

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
   Get the code and the included libraries onto your local machine:
   ```bash
   git clone https://github.com/QiangSu/PatternFilter.git
   cd PatternFilter
   ```
   
   Compilation
   ```bash
   g++ -std=c++17 -O3 -D USE_PARALLEL_SORT -o Pattern_Filter Pattern_Filter.cpp -lz -pthread -ltbb -I ./include -I ./gzstream ./gzstream/gzstream.C
   ```
## Usage:
   Run the program with --help to see all available options:
   ```bash
   ./Pattern_Filter --help
   ```
   An example command using the test data (make sure test_data/ exists from the clone):
   mkdir -p R1_passed R2_extracted Hamming_rejected Basecomp_rejected Test_Temp_Dir
```bash
./Pattern_Filter \
  --r1_input test_data/extracted_100000_reads_musWT_STR_R1_001.fastq.gz \
  --r2_input test_data/extracted_100000_reads_musWT_STR_R2_001.fastq.gz \
  --r1_output_dir R1_passed \
  --r2_output_dir R2_extracted \
  --hamming_filterout_dir Hamming_rejected \
  --base_composition_filterout_dir Basecomp_rejected \
  --target_seq "TTTTTTTTTTTTTTT" \
  --seq_length 15 \
  --start_pos 28 \ # Positions are 0-based in this c++ program. If a desired position is 29 in a 1-based system, you should input 28.
  --threshold 2 \
  --start_pos_base_com 0 \
  --seq_length_base_com 28 \
  --base_composition_threshold 1 \
  --threads 40 \
  --sort-threads 40 \
  --chunk_records 10000000 \ # Suitable for processing very large datasets (e.g., sc/snRNA-seq)
  --sort-mem-mb 204800 \     # Approx. memory limit (MB) per sort worker - adjust based on system RAM
  --temp_dir Test_Temp_Dir \ # Base directory for temporary files
  --parallel-sort-alg true   # Enable parallel algorithm for in-memory sort step
```


# scRNA-Seq Simulation Pipeline (Splatter + Polyester + Custom Barcodes/UMIs)

This project provides a two-step pipeline to simulate single-cell RNA sequencing data, generating both ground truth counts and Cell Ranger-compatible FASTQ files. It combines the power of Splatter for realistic count simulation, Polyester for generating reads based on these counts, and a custom Python script to add realistic 10x Genomics barcodes and UMIs from a provided whitelist.

The final output (`.fastq.gz` files) serves as high-quality simulated "ground truth" data that can be directly processed by pipelines like Cell Ranger, enabling benchmarking and testing.

## Overview

The simulation process is split into two main steps:

1.  **Count Simulation (R - `step1_simulate_scRNA_seq_data.R`):** Uses the `splatter` package to simulate a matrix of gene expression counts for a specified number of cells and genes, incorporating realistic biological variation, differential expression between simulated cell types, and dropout. It selects genes from a provided reference transcriptome FASTA file.
2.  **FASTQ Generation and Formatting (Python - `step2_add_barcodes_umis_cellranger.py`):** Takes the counts simulated in Step 1, uses the `polyester` package (via its R interface, but processed in the Python script's logic) to generate read sequences, then crucially processes these sequences to:
    *   Assign a unique, valid 10x barcode (sampled from a whitelist) to each simulated cell.
    *   Generate random UMIs.
    *   Construct Cell Ranger-compatible paired-end reads (`_R1_` containing Barcode+UMI, `_R2_` containing the cDNA sequence).
    *   Combine reads from all simulated cells into a single pair of gzipped FASTQ files ready for Cell Ranger `count`.

## Simulation Pipeline

### Step 1: Simulate Counts (R)

*   **Script:** `step1_simulate_scRNA_seq_data.R`
*   **Tools:** `splatter`, `polyester`, `Biostrings`
*   **Input:** A reference transcriptome FASTA file (used to select realistic gene IDs).
*   **Process:**
    *   Sets a seed for reproducibility.
    *   Loads required R libraries.
    *   Defines base directories and input/output file paths.
    *   Reads gene/transcript IDs from the reference FASTA file.
    *   Initializes Splatter parameters, setting the number of cells, number of genes (selected from the FASTA), number of groups (cell types), and parameters controlling variability, differential expression, and dropout.
    *   Simulates the data using `splatSimulate`, generating a `SummarizedExperiment` object.
    *   Extracts the count matrix and cell group assignments (simulated cell types).
    *   Assigns the selected gene IDs from the FASTA as row names and "Cell_X" as column names.
    *   Saves the count matrix and cell type assignments as ground truth files.
    *   *Intermediate Step:* Uses Polyester (within this script) to generate initial *simulated* FASTQ files based on the simulated counts and the reference FASTA. These are initially named `sample_X_Y.fasta` and stored in a subdirectory. The script renames these files to remove leading zeros (`sample_01` becomes `sample_1`) but keeps the `.fasta` extension.
*   **Outputs:**
    *   `ground_truth_matrix.csv`: A CSV file containing the simulated count matrix (Genes as rows, Cells as columns).
    *   `ground_truth_celltypes.txt`: A tab-separated file mapping Cell IDs (`Cell_1`, etc.) to their simulated group/cell type.
    *   `simulated_fastq_.../sample_X_1.fasta`, `simulated_fastq_.../sample_X_2.fasta` (where X is 1 to `n_cells`): Intermediate FASTQ files generated by Polyester for each simulated cell. Note: Polyester's R1 contains the cDNA sequence in this typical configuration, and the names retain the `.fasta` extension after renaming leading zeros.
    *   `temp_filtered_transcriptome.fa`: A temporary FASTA file used internally by Polyester, which is removed upon successful completion.

### Step 2: Generate & Format FASTQ (Python)

*   **Script:** `step2_add_barcodes_umis_cellranger.py`
*   **Tools:** Standard Python libraries (`gzip`, `random`, `os`, `sys`)
*   **Input:**
    *   The intermediate Polyester FASTQ files generated by Step 1 (`simulated_fastq_.../sample_X_Y.fasta`).
    *   A gzipped 10x Genomics barcode whitelist file (`.txt.gz`).
*   **Process:**
    *   Sets base directories and input/output paths.
    *   Defines simulation parameters, importantly matching the number of cells (`n_cells`) from the R script.
    *   Reads and validates barcodes from the provided whitelist.
    *   Randomly samples `n_cells` *unique* barcodes from the whitelist and assigns one to each simulated cell (corresponding to `sample_1`, `sample_2`, etc.).
    *   Saves a mapping of input cell index to assigned barcode.
    *   Iterates through the intermediate Polyester FASTQ files (`sample_1_1.fasta`, `sample_2_1.fasta`, ..., `sample_N_1.fasta`).
    *   For each read in a `sample_X_1.fasta` file (which contains the simulated cDNA sequence):
        *   Generates a random UMI.
        *   Constructs the Cell Ranger R1 sequence: `assigned_barcode + UMI`.
        *   Uses the sequence from the Polyester file as the Cell Ranger R2 sequence.
        *   Creates Cell Ranger-compatible FASTQ headers embedding sample info and the barcode.
        *   Assigns placeholder quality scores (Phred 40 equivalent).
    *   Writes the newly constructed paired-end reads to a *single pair* of gzipped output FASTQ files (`_R1_...fastq.gz` and `_R2_...fastq.gz`) in the Cell Ranger standard format.
*   **Outputs:**
    *   `cellranger_simulated_data_.../simulated_S1_R1_001.fastq.gz`: The Cell Ranger compatible Read 1 file, containing Barcode and UMI sequences.
    *   `cellranger_simulated_data_.../simulated_S1_R2_001.fastq.gz`: The Cell Ranger compatible Read 2 file, containing the simulated cDNA sequences.
    *   `cellranger_simulated_data_.../simulated_barcode_mapping.txt`: A file showing which simulated cell (by input index `1` to `n_cells`) was assigned which specific 10x whitelist barcode.

## Prerequisites

*   R (>= 4.0)
    *   BiocManager
    *   splatter
    *   polyester
    *   Biostrings
    *   SummarizedExperiment
    *   Matrix
*   Python 3
*   A reference transcriptome FASTA file (e.g., `Mus_musculus.GRCm38.cdna.all.fa`). **Note:** The script expects transcript IDs in the header before the first space or dot (e.g., `ENSMUST00000123456` from `>ENSMUST00000123456.1 cdna...`). Ensure your FASTA headers are in a compatible format or modify the R script's parsing logic.
*   A 10x Genomics barcode whitelist file (`.txt.gz`, e.g., `3M-february-2018.txt.gz`).
*   Sufficient disk space (simulated FASTQ files, especially the intermediate ones, can be very large).

## Installation

1.  Clone this repository.
2.  Install the necessary R packages:

    ```R
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("splatter")
    BiocManager::install("polyester")
    BiocManager::install("Biostrings")
    BiocManager::install("SummarizedExperiment") # Often installed with others, but good to check
    # Matrix is usually base R, but explicitly loading is good practice
    ```
3.  Python 3 and its standard libraries (`gzip`, `random`, `os`, `sys`) are typically pre-installed or easy to obtain via package managers (apt, yum, brew) or Anaconda/Miniconda.

## Usage

1.  **Set Paths and Parameters:**
    *   Edit **`step1_simulate_scRNA_seq_data.R`**:
        *   Modify `base_dir` to the desired root output directory.
        *   Modify `reference_fasta` to the path of your reference transcriptome.
        *   Adjust `n_cells`, `n_genes`, `n_groups`, and other `splatSimulate` parameters as needed.
    *   Edit **`step2_add_barcodes_umis_cellranger.py`**:
        *   Modify `base_dir` to the *same* directory as set in the R script.
        *   Modify `barcode_whitelist_path` to the path of your 10x whitelist.
        *   Ensure `n_cells` matches the value set in the R script.
        *   Adjust `cr_sample_name` if you want a different prefix for the final Cell Ranger files.
        *   Confirm `barcode_length`, `umi_length`, and `cr_read2_length` match your intended Cell Ranger chemistry and R simulation parameters (`readlen`).

2.  **Run Step 1 (R):** Execute the R script from your terminal.

    ```bash
    Rscript step1_simulate_scRNA_seq_data.R
    ```
    This will generate the count matrix, cell types file, and the intermediate Polyester `.fasta` files in the specified output directory.

3.  **Run Step 2 (Python):** Once Step 1 has completed successfully, execute the Python script.

    ```bash
    python step2_add_barcodes_umis_cellranger.py
    ```
    This will read the intermediate `.fasta` files, add barcodes/UMIs, and write the final gzipped `.fastq.gz` files in the Cell Ranger output directory.

## Outputs for Cell Ranger

The final Cell Ranger compatible FASTQ files will be located in the directory specified by `cellranger_simulated_data_...` (derived from `base_dir` and `cr_sample_name`) and named like:

*   `simulated_S1_R1_001.fastq.gz`
*   `simulated_S1_R2_001.fastq.gz`

You can now point Cell Ranger `count` (or similar pipelines) to the directory containing these files using the `--fastqs` argument and specify the sample name using the `--sample` argument (which is `simulated` by default in the Python script, derived from `cr_sample_name`).

Example Cell Ranger command structure:

```bash
cellranger count --id=my_simulation_analysis \
                 --transcriptome=/path/to/your/cellranger/index \
                 --fastqs=/path/to/your/base_dir/cellranger_simulated_data_... \
                 --sample=simulated \
                 --localcores=16 --localmem=64

   
