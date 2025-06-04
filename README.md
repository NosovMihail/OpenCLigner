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

Download `opencligner` and all available `.cl` cores. Place them in the same folder. To run the programm you should have OpenCL support on your system. It can be obtained by executing:
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

## Files
- main.cpp: Main entry point

- needleman_wunsch_wavefront.cl: OpenCL kernel for global alignment

- smith_waterman_wavefront.cl: OpenCL kernel for local alignment

- Dockerfile: Docker build configuration

- LICENSE.txt: License information

## License
This project is licensed under the MIT License - see the [LICENSE.txt](https://github.com/NosovMihail/OpenCLigner/blob/master/LICENSE.txt) file for details.

## Author
Created by Mikhail Nosov
