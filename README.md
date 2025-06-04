# OpenCLigner

**OpenCLigner** is a sequence alignment tool written in C++ and OpenCL. It supports global (Needleman-Wunsch) and local (Smith-Waterman) alignment methods and is optimized for performance using parallel processing on compatible hardware. The tool can be compiled and run on Linux and Windows, and a Docker container is provided for easy deployment.

## Features

- Global alignment (Needleman-Wunsch algorithm)
- Local alignment (Smith-Waterman algorithm)
- OpenCL acceleration (wavefront parallelization)
- FASTA input file support
- Output in aligned FASTA format
- Docker-based build and execution
- Simple command-line interface

## Usage

### Installation

- Download **Release** OR `opencligner` and all available `.cl` cores
- Place them in the same folder
- To run the programm you should have OpenCL support on your system. It can be obtained by executing:
```bash
sudo apt install ocl-icd-opencl-dev
```
But GPU manufacturer driver is a recomended option.

### Command-line Options

```bash
./opencligner [options] -i <input1.fasta> -i <input2.fasta>
```
#### Options:

- -o <file>: Output file path (default: output.fasta)

- -m <method>: Alignment method (global or local) (default: global)

- -h, --help: Show usage instructions

#### Example:

```bash
./opencligner -m global -o result.fasta -i sequence1.fasta -i sequence2.fasta
```

## Building
### Requirements
- C++ compiler (e.g., g++)

- OpenCL SDK (e.g., Intel, AMD, NVIDIA)

- CMake (optional)

- Docker (optional, for Docker build)

### Docker Build
To build and run inside a Docker container:

```bash
docker build -t opencligner .
docker run --rm -v $(pwd):/app opencligner -m global -o /app/output.fasta /app/seq1.fasta /app/seq2.fasta
```
This builds the project inside an Ubuntu container with OpenCL headers and compilers installed.

## Possible Improvements and Future Work

- **Multiple Sequence Alignment (MSA)**  
  Extend support from pairwise alignment to multiple sequences using progressive or iterative methods.

- **Affine Gap Penalties**  
  Implement support for affine gap scoring (gap opening and gap extension penalties) to improve biological accuracy.

- **Performance Optimization**  
  - Explore tile-based or diagonal parallelism in OpenCL for better scalability  
  - Add device-specific optimizations (e.g., vectorization on GPUs)

- **FASTA Parsing Enhancements**  
  Improve support for multiline sequences and unusual formatting in FASTA headers.

- **Alignment Scoring Report**  
  Output alignment score, identity percentage, and mismatch/gap statistics.

- **Visual Output Format**  
  Export alignment result in:
  - Clustal format
  - SAM/BAM format for integration with bioinformatics tools
  - Visual HTML or SVG-based alignment view

- **Support for Substitution Matrices**  
  Allow the use of scoring matrices like BLOSUM or PAM for protein sequences.

- **Validation**  
  Add benchmarks comparing results with established tools (e.g., EMBOSS, NCBI BLAST).

- **Cross-Platform GUI**  
  Optionally provide a simple desktop interface for users unfamiliar with the command line.

- **Library API**  
  Expose alignment as a C++ library or C-compatible API for embedding in other software.

- **OpenCL Device Selection**  
  Allow users to choose which OpenCL platform and device to use for execution.

- **Parallel Input Processing**  
  Batch processing of multiple FASTA files with parallel alignment execution.


## Files
- `main.cpp`: Main entry point

- `needleman_wunsch_wavefront.cl`: OpenCL kernel for global alignment

- `smith_waterman_wavefront.cl`: OpenCL kernel for local alignment

- `Dockerfile`: Docker build configuration

- `LICENSE.txt`: License information

## License
This project is licensed under the MIT License - see the [LICENSE.txt](https://github.com/NosovMihail/OpenCLigner/blob/master/LICENSE.txt) file for details.

## Author
Created by Mikhail Nosov
