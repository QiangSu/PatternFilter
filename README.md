# Fastq Pattern Filter and Sorter

This C++17 program filters paired-end FASTQ reads based on Hamming distance to a target sequence in R1 and base composition criteria in R1. It then sorts the R1 reads into 'passed', 'Hamming rejected', and 'base composition rejected' categories. Corresponding R2 reads for the 'passed' R1 set are extracted and sorted.

## Dependencies

*   A C++17 compliant compiler (e.g., g++ 7 or later)
*   Zlib development library (e.g., `zlib1g-dev` on Debian/Ubuntu, `zlib-devel` on CentOS/Fedora)
*   pthreads (usually included with the compiler toolchain)

Included header-only/source libraries:
*   `cxxopts.hpp` (for command-line argument parsing) - in `include/`
*   `gzstream` (for reading/writing gzipped files) - in `gzstream/`

## Compilation

To compile the program, navigate to the project root directory and run:

```bash
g++ -std=c++17 -O3 -Wall -I./include -I./gzstream ./Pattern_Filter.cpp ./gzstream/gzstream.C -o Pattern_Filter -pthread -lz -DUSE_PARALLEL_SORT

