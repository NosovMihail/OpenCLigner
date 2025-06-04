#define CL_TARGET_OPENCL_VERSION 120
#include <CL/cl.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>
#include <algorithm>

#define CHECK_ERROR(err, msg) \
    if (err != CL_SUCCESS) { std::cerr << msg << ": " << err << std::endl; exit(1); }

#include <sstream>
#include <stdexcept>

std::pair<std::string, std::string> readFasta(const std::string & filename) {
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Cannot open input file: " + filename);
    }

    std::string header;
    std::stringstream seq;
    std::string line;

    while (std::getline(file, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!header.empty()) {
                throw std::runtime_error("Only one sequence per FASTA file is supported.");
            }
            header = line.substr(1);  // remove '>'
        }
        else {
            seq << line;
        }
    }

    if (header.empty()) {
        throw std::runtime_error("FASTA file does not contain a header line starting with '>'");
    }

    return { header, seq.str() };
}

void writeAlignmentFasta(
    const std::string& filename,
    const std::string& header1,
    const std::string& align1,
    const std::string& header2,
    const std::string& align2)
{
    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error("Cannot open output file: " + filename);
    }

    out << ">" << header1 << "\n" << align1 << "\n";
    out << ">" << header2 << "\n" << align2 << "\n";
}

void printHelp() {
    std::cout << "Usage:\n"
        << "  aligner -i <file1> -i <file2> [-i <fileN> ...] -o <output_file> [-m global|local]\n\n"
        << "Options:\n"
        << "  -i, --input <file>     Input file with one sequence (can be repeated)\n"
        << "  -o, --output <file>    Output file for alignment result\n"
        << "  -m, --method <mode>    Alignment method: global or local (default: global)\n"
        << "  -h, --help             Show this help message\n";
}

int main(int argc, char* argv[]) {
    std::vector<std::string> inputFiles;
    std::string outputFile = "output.fasta";
    std::string method = "global";
    const char* seq1 = "GAGCGTAGTCCCTTCAGAGAAAGCTCTCCGGATCAGAAAGTGTCGTAGAACGGTGTGTTGCCTTCTCTCCACAGTCTGTGGTAGGATACGCAACACGCTGGACCGATGCTAGHAGGGGAAGGA";
    const char* seq2 = "CCTCAAAAGAAGCTCTCTGGATCAGAACGGTCTTAGTACGGGGTGTTGCATGTATTCCCCACAGTTTGTGGTAGGAGACGCAACACGCTGGAAATGTCTCTAGTAGGAGAAGGAAAGCTTAGTCCCGCT";
    int len1 = strlen(seq1);
    int len2 = strlen(seq2);
    const int match = 2, mismatch = -1, gap = -1;

    // Обработка аргументов
    for (int i = 1; i < argc; ++i) {
        if ((strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0)) {
            printHelp();
            return 0;
        }
        else if ((strcmp(argv[i], "-i") == 0 || strcmp(argv[i], "--input") == 0) && i + 1 < argc) {
            inputFiles.push_back(argv[++i]);
        }
        else if ((strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--output") == 0) && i + 1 < argc) {
            outputFile = argv[++i];
        }
        else if ((strcmp(argv[i], "-m") == 0 || strcmp(argv[i], "--method") == 0) && i + 1 < argc) {
            method = argv[++i];
            if (method != "global" && method != "local") {
                std::cerr << "Error. Method shold be 'global' or 'local'\n";
                return 1;
            }
        }
        else {
            std::cerr << "Unknown option: " << argv[i] << "\n";
            printHelp();
            return 1;
        }
    }

    if (inputFiles.size() < 2) {
        std::cerr << "Error. At least two input files required.\n";
        printHelp();
        return 1;
    }

    std::vector<std::string> sequences;
    std::vector<std::string> headers;

    for (const auto& file : inputFiles) {
        try {
            std::pair<std::string, std::string> result = readFasta(file);
            headers.push_back(result.first);     // Заголовок (без '>')
            sequences.push_back(result.second);  // Последовательность
        }
        catch (const std::exception& e) {
            std::cerr << "Error reading file '" << file << "': " << e.what() << "\n";
            return 1;
        }
    }

    if (sequences.size() < 2) {
        std::cerr << "Error: at least two input files with FASTA sequences are required.\n";
        return 1;
    }

    std::string sequence1 = sequences[0];
    std::string sequence2 = sequences[1];
    std::string header1 = headers[0];
    std::string header2 = headers[1];

    len1 = sequence1.length();
    len2 = sequence2.length();

    std::ifstream kernelFile("smith_waterman_wavefront.cl");
    std::string kernelSource((std::istreambuf_iterator<char>(kernelFile)),
        std::istreambuf_iterator<char>());
    const char* kernelSourceStr = kernelSource.c_str();

    std::ifstream kernelFile2("needleman_wunsch_wavefront.cl");
    std::string kernelSource2((std::istreambuf_iterator<char>(kernelFile2)),
        std::istreambuf_iterator<char>());
    const char* kernelSourceStr2 = kernelSource2.c_str();

    std::vector<int> matrix((len1 + 1) * (len2 + 1), 0);
    int max_score = 0;
    int max_i_host = 0;
    int max_j_host = 0;

    cl_int err;
    cl_platform_id platform;
    cl_device_id device;
    cl_context context;
    cl_command_queue queue;

    err = clGetPlatformIDs(1, &platform, nullptr);
    CHECK_ERROR(err, "clGetPlatformIDs");

    err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_DEFAULT, 1, &device, nullptr);
    CHECK_ERROR(err, "clGetDeviceIDs");

    context = clCreateContext(nullptr, 1, &device, nullptr, nullptr, &err);
    CHECK_ERROR(err, "clCreateContext");

    queue = clCreateCommandQueue(context, device, 0, &err);
    CHECK_ERROR(err, "clCreateCommandQueue");

    cl_program program = (method == "global") ? clCreateProgramWithSource(context, 1, &kernelSourceStr2, nullptr, &err) : clCreateProgramWithSource(context, 1, &kernelSourceStr, nullptr, &err); //kernelSourceStr/2
    CHECK_ERROR(err, "clCreateProgramWithSource");

    err = clBuildProgram(program, 0, nullptr, nullptr, nullptr, nullptr);
    if (err != CL_SUCCESS) {
        size_t log_size;
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, nullptr, &log_size);
        std::vector<char> log(log_size);
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, log_size, log.data(), nullptr);
        std::cerr << "Build log:\n" << log.data() << std::endl;
        exit(1);
    }

    cl_kernel kernel = (method == "global") ? clCreateKernel(program, "needleman_wunsch_wavefront", &err) : clCreateKernel(program, "smith_waterman_wavefront", &err); //change method
    CHECK_ERROR(err, "clCreateKernel");

    cl_mem seq1_buf = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, len1, (void*)sequence1.data(), &err);
    cl_mem seq2_buf = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, len2, (void*)sequence2.data(), &err);
    cl_mem matrix_buf = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(int) * matrix.size(), matrix.data(), &err);
    cl_mem max_score_buf = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(int), &max_score, &err);
    cl_mem max_i_buf = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
        sizeof(int), &max_i_host, &err);
    cl_mem max_j_buf = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
        sizeof(int), &max_j_host, &err);

    clSetKernelArg(kernel, 0, sizeof(cl_mem), &seq1_buf);
    clSetKernelArg(kernel, 1, sizeof(cl_mem), &seq2_buf);
    clSetKernelArg(kernel, 2, sizeof(cl_mem), &matrix_buf);
    clSetKernelArg(kernel, 3, sizeof(cl_mem), &max_score_buf);
    clSetKernelArg(kernel, 4, sizeof(int), &len1);
    clSetKernelArg(kernel, 5, sizeof(int), &len2);
    clSetKernelArg(kernel, 7, sizeof(int), &match);
    clSetKernelArg(kernel, 8, sizeof(int), &mismatch);
    clSetKernelArg(kernel, 9, sizeof(int), &gap);
    clSetKernelArg(kernel, 10, sizeof(cl_mem), &max_i_buf);
    clSetKernelArg(kernel, 11, sizeof(cl_mem), &max_j_buf);

    for (int diag = 2; diag <= len1 + len2; ++diag) {
        int n = std::min(diag - 1, std::min(len1, len2));
        clSetKernelArg(kernel, 6, sizeof(int), &diag);

        size_t global_work_size = n;
        err = clEnqueueNDRangeKernel(queue, kernel, 1, nullptr, &global_work_size, nullptr, 0, nullptr, nullptr);
        CHECK_ERROR(err, "clEnqueueNDRangeKernel");
        clFinish(queue);
    }

    clEnqueueReadBuffer(queue, matrix_buf, CL_TRUE, 0, sizeof(int) * matrix.size(), matrix.data(), 0, nullptr, nullptr);
    clEnqueueReadBuffer(queue, max_score_buf, CL_TRUE, 0, sizeof(int), &max_score, 0, nullptr, nullptr);
    clEnqueueReadBuffer(queue, max_i_buf, CL_TRUE, 0, sizeof(int), &max_i_host, 0, nullptr, nullptr);
    clEnqueueReadBuffer(queue, max_j_buf, CL_TRUE, 0, sizeof(int), &max_j_host, 0, nullptr, nullptr);

    /*std::cout << "Matrix (" << len1 + 1 << " x " << len2 + 1 << "):\n";
    for (int i = 0; i <= len1; i++) {
        for (int j = 0; j <= len2; j++) {
            std::cout << matrix[i * (len2 + 1) + j] << "\t";
        }
        std::cout << "\n";
    }*/

    std::cout << "Max score: " << max_score << " at (" << max_i_host << ", " << max_j_host << ")\n";

    // Traceback
    std::string aligned_seq1 = "";
    std::string aligned_seq2 = "";

    int i = max_i_host;
    int j = max_j_host;

    while (i > 0 && j > 0 && matrix[i * (len2 + 1) + j] > 0) {
        int score = matrix[i * (len2 + 1) + j];
        int diag = matrix[(i - 1) * (len2 + 1) + (j - 1)];
        int up = matrix[(i - 1) * (len2 + 1) + j];
        int left = matrix[i * (len2 + 1) + (j - 1)];
        int match2 = match;
        if (seq1[i - 1] != seq2[j - 1]) {
            match2 = mismatch;
        }

        if (score == diag + match2) {
            aligned_seq1 = seq1[i - 1] + aligned_seq1;
            aligned_seq2 = seq2[j - 1] + aligned_seq2;
            i--;
            j--;
            //std::cout << seq1[i - 1] << " " << j << " " << i << "\n";
        }
        else if (score == up + gap) {
            aligned_seq1 = seq1[i - 1] + aligned_seq1;
            aligned_seq2 = "-" + aligned_seq2;
            i--;
            //std::cout << seq1[i - 1] << " " << j << " " << i << "\n";
        }
        else if (score == left + gap) {
            aligned_seq1 = "-" + aligned_seq1;
            aligned_seq2 = seq2[j - 1] + aligned_seq2;
            j--;
            //std::cout << seq2[j - 1] << " " << j << " " << i << "\n";
        }
        else {
            //std::cout << j << " " << i << "\n";
            break;
        }
    }

    //Debug lines
    std::cout << "Alignment result:\n";
    std::cout << aligned_seq1 << "\n";
    std::cout << aligned_seq2 << "\n";

    std::ofstream out(outputFile);
    if (!out.is_open()) {
        std::cerr << "Error, can't open file: " << outputFile << std::endl;
        return 1;
    }

    writeAlignmentFasta(outputFile, header1, aligned_seq1, header2, aligned_seq2);
    std::cout << "Result written to: " << outputFile << std::endl;

    // Cleanup
    clReleaseMemObject(seq1_buf);
    clReleaseMemObject(seq2_buf);
    clReleaseMemObject(matrix_buf);
    clReleaseMemObject(max_score_buf);
    clReleaseKernel(kernel);
    clReleaseProgram(program);
    clReleaseCommandQueue(queue);
    clReleaseContext(context);

    return 0;
}