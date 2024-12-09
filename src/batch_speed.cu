#include <cuda_runtime.h>
#include <cusparse_v2.h>
#include <iostream>
#include <vector>
#include <cmath>

#define CHECK_CUDA(call)                                                                 \
    {                                                                                   \
        cudaError_t err = call;                                                         \
        if (err != cudaSuccess) {                                                       \
            std::cerr << "CUDA error in " << __FILE__ << ":" << __LINE__ << ": "        \
                      << cudaGetErrorString(err) << std::endl;                          \
            exit(EXIT_FAILURE);                                                         \
        }                                                                               \
    }

#define CHECK_CUSPARSE(call)                                                            \
    {                                                                                   \
        cusparseStatus_t err = call;                                                    \
        if (err != CUSPARSE_STATUS_SUCCESS) {                                           \
            std::cerr << "cuSPARSE error in " << __FILE__ << ":" << __LINE__ << ": "    \
                      << err << std::endl;                                              \
            exit(EXIT_FAILURE);                                                         \
        }                                                                               \
    }






void test_cusparseDgtsv2StridedBatch_2(int m, int batchCount) {
    int batchStride = m; // Minimum batch stride

    // Allocate memory for diagonals and RHS
    size_t diag_size = batchCount * batchStride * sizeof(double);
    double *dl, *d, *du, *x, *x_2;
    CHECK_CUDA(cudaMalloc(&dl, diag_size));
    CHECK_CUDA(cudaMalloc(&d, diag_size));
    CHECK_CUDA(cudaMalloc(&du, diag_size));
    CHECK_CUDA(cudaMalloc(&x, diag_size));
    CHECK_CUDA(cudaMalloc(&x_2, diag_size));

    // Initialize diagonals and RHS with dummy data
    std::vector<double> h_dl(batchCount * batchStride, -1.0);
    std::vector<double> h_d(batchCount * batchStride, 4.0);
    std::vector<double> h_du(batchCount * batchStride, -1.0);
    std::vector<double> h_x(batchCount * batchStride, 1.0);
    std::vector<double> h_x_2(batchCount * batchStride, 2.0);


    CHECK_CUDA(cudaMemcpy(dl, h_dl.data(), diag_size, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d, h_d.data(), diag_size, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(du, h_du.data(), diag_size, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(x, h_x.data(), diag_size, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(x_2, h_x_2.data(), diag_size, cudaMemcpyHostToDevice));


    // Create cuSPARSE handle
    cusparseHandle_t handle;
    CHECK_CUSPARSE(cusparseCreate(&handle));

    // Query buffer size
    void *pBuffer = nullptr;
    size_t pBufferSizeInBytes;
    CHECK_CUSPARSE(cusparseDgtsv2StridedBatch_bufferSizeExt(handle, m, dl, d, du, nullptr, batchCount, batchStride, &pBufferSizeInBytes));
    CHECK_CUDA(cudaMalloc(&pBuffer, pBufferSizeInBytes));

    std::cout<<pBufferSizeInBytes<<std::endl;

    // Timing setup
    cudaEvent_t start, stop;
    CHECK_CUDA(cudaEventCreate(&start));
    CHECK_CUDA(cudaEventCreate(&stop));

    // Solve systems and measure time
    CHECK_CUDA(cudaEventRecord(start));
    CHECK_CUSPARSE(cusparseDgtsv2StridedBatch(handle, m, dl, d, du, x, batchCount, batchStride, pBuffer));
    CHECK_CUSPARSE(cusparseDgtsv2StridedBatch(handle, m, dl, d, du, x_2, batchCount, batchStride, pBuffer));
    CHECK_CUDA(cudaEventRecord(stop));
    CHECK_CUDA(cudaEventSynchronize(stop));

    float milliseconds = 0;
    CHECK_CUDA(cudaEventElapsedTime(&milliseconds, start, stop));
    std::cout << "Execution time: " << milliseconds << " ms" << std::endl;

    // Calculate and print norm of x
    std::vector<double> h_x_result(batchCount * batchStride);
    CHECK_CUDA(cudaMemcpy(h_x_result.data(), x, diag_size, cudaMemcpyDeviceToHost));
    double norm = 0.0;
    for (const double &val : h_x_result) {
        norm += val * val;
    }
    norm = std::sqrt(norm);
    std::cout << "Norm of x: " << norm << std::endl;

    std::vector<double> h_x_result_2(batchCount * batchStride);
    CHECK_CUDA(cudaMemcpy(h_x_result_2.data(), x_2, diag_size, cudaMemcpyDeviceToHost));
    double norm2 = 0.0;
    for (const double &val : h_x_result_2) {
        norm2 += val * val;
    }
    norm2 = std::sqrt(norm2);
    std::cout << "Norm of x: " << norm2 << std::endl;


    // Clean up
    cudaFree(dl);
    cudaFree(d);
    cudaFree(du);
    cudaFree(x);
    cudaFree(pBuffer);
    cusparseDestroy(handle);
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
}







void test_cusparseDgtsv2StridedBatch(int m, int batchCount) {
    int batchStride = m; // Minimum batch stride

    // Allocate memory for diagonals and RHS
    size_t diag_size = batchCount * batchStride * sizeof(double);
    double *dl, *d, *du, *x;
    CHECK_CUDA(cudaMalloc(&dl, diag_size));
    CHECK_CUDA(cudaMalloc(&d, diag_size));
    CHECK_CUDA(cudaMalloc(&du, diag_size));
    CHECK_CUDA(cudaMalloc(&x, diag_size));

    // Initialize diagonals and RHS with dummy data
    std::vector<double> h_dl(batchCount * batchStride, -1.0);
    std::vector<double> h_d(batchCount * batchStride, 4.0);
    std::vector<double> h_du(batchCount * batchStride, -1.0);
    std::vector<double> h_x(batchCount * batchStride, 1.0);


    CHECK_CUDA(cudaMemcpy(dl, h_dl.data(), diag_size, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d, h_d.data(), diag_size, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(du, h_du.data(), diag_size, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(x, h_x.data(), diag_size, cudaMemcpyHostToDevice));

    // Create cuSPARSE handle
    cusparseHandle_t handle;
    CHECK_CUSPARSE(cusparseCreate(&handle));

    // Query buffer size
    void *pBuffer = nullptr;
    size_t pBufferSizeInBytes;
    CHECK_CUSPARSE(cusparseDgtsv2StridedBatch_bufferSizeExt(handle, m/2, nullptr, nullptr, nullptr, nullptr, batchCount, 2*batchStride, &pBufferSizeInBytes));
    CHECK_CUDA(cudaMalloc(&pBuffer, pBufferSizeInBytes));

    std::cout<<pBufferSizeInBytes<<std::endl;
    std::cout<< 5 * batchCount * m * sizeof(double) <<std::endl;



    // Timing setup
    cudaEvent_t start, stop;
    CHECK_CUDA(cudaEventCreate(&start));
    CHECK_CUDA(cudaEventCreate(&stop));

    // Solve systems and measure time
    CHECK_CUDA(cudaEventRecord(start));
    CHECK_CUSPARSE(cusparseDgtsv2StridedBatch(handle, m, dl, d, du, x, batchCount, batchStride, pBuffer));
    CHECK_CUDA(cudaEventRecord(stop));
    CHECK_CUDA(cudaEventSynchronize(stop));

    float milliseconds = 0;
    CHECK_CUDA(cudaEventElapsedTime(&milliseconds, start, stop));
    std::cout << "Execution time: " << milliseconds << " ms" << std::endl;

    // Calculate and print norm of x
    std::vector<double> h_x_result(batchCount * batchStride);
    CHECK_CUDA(cudaMemcpy(h_x_result.data(), x, diag_size, cudaMemcpyDeviceToHost));
    double norm = 0.0;
    for (const double &val : h_x_result) {
        norm += val * val;
    }
    norm = std::sqrt(norm);
    std::cout << "Norm of x: " << norm << std::endl;

    // Clean up
    cudaFree(dl);
    cudaFree(d);
    cudaFree(du);
    cudaFree(x);
    cudaFree(pBuffer);
    cusparseDestroy(handle);
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
}

int main() {
    const int m = 4096;         // System size
    const int batchCount = 427; // Number of systems

    test_cusparseDgtsv2StridedBatch(m, batchCount * 2);

    // test_cusparseDgtsv2StridedBatch_2(m, batchCount);

    return 0;
}