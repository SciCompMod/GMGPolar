#pragma once

#include "gpu_vector.h"

#include <cuda_runtime.h>
#include <iostream>
#include <stdexcept>

// Kernel function to assign a value to all elements of the vector
template <typename T>
__global__ void assign_kernel(T* device_values, int size, T value)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < size) {
        device_values[idx] = value;
    }
}

// Host function for assigning a value to all elements of the vector
template <typename T>
void assign(GPU_Vector<T>& v, T value)
{
    int block_size = 256;
    int num_blocks = (v.size() + block_size - 1) / block_size;

    assign_kernel<T><<<num_blocks, block_size>>>(v.data(), v.size(), value);
    cudaDeviceSynchronize();
}

// Kernel function to add two vectors element-wise
template <typename T>
__global__ void add_kernel(T* result, const T* x, int size)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < size) {
        result[idx] += x[idx];
    }
}

// Host function for adding two vectors element-wise
template <typename T>
void add(GPU_Vector<T>& result, const GPU_Vector<T>& x)
{
    if (result.size() != x.size()) {
        throw std::invalid_argument("Vectors must be of the same size.");
    }

    int block_size = 256;
    int num_blocks = (result.size() + block_size - 1) / block_size;

    add_kernel<T><<<num_blocks, block_size>>>(result.data(), x.data(), result.size());
    cudaDeviceSynchronize();
}

// Kernel function to subtract two vectors element-wise
template <typename T>
__global__ void subtract_kernel(T* result, const T* x, int size)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < size) {
        result[idx] -= x[idx];
    }
}

// Host function for subtracting two vectors element-wise
template <typename T>
void subtract(GPU_Vector<T>& result, const GPU_Vector<T>& x)
{
    if (result.size() != x.size()) {
        throw std::invalid_argument("Vectors must be of the same size.");
    }

    int block_size = 256;
    int num_blocks = (result.size() + block_size - 1) / block_size;

    subtract_kernel<T><<<num_blocks, block_size>>>(result.data(), x.data(), result.size());
    cudaDeviceSynchronize();
}

// Kernel function for a linear combination of two vectors
template <typename T>
__global__ void linear_combination_kernel(T* x, T alpha, const T* y, T beta, int size)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < size) {
        x[idx] = alpha * x[idx] + beta * y[idx];
    }
}

// Host function for a linear combination of two vectors
template <typename T>
void linear_combination(GPU_Vector<T>& x, const T alpha, const GPU_Vector<T>& y, const T beta)
{
    if (x.size() != y.size()) {
        throw std::invalid_argument("Vectors must be of the same size.");
    }

    int block_size = 256;
    int num_blocks = (x.size() + block_size - 1) / block_size;

    linear_combination_kernel<T><<<num_blocks, block_size>>>(x.data(), alpha, y.data(), beta, x.size());
    cudaDeviceSynchronize();
}

// Kernel function to multiply a vector by a scalar
template <typename T>
__global__ void multiply_kernel(T* result, int size, T alpha)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < size) {
        result[idx] *= alpha;
    }
}

// Host function for multiplying a vector by a scalar
template <typename T>
void multiply(GPU_Vector<T>& v, T alpha)
{
    int block_size = 256;
    int num_blocks = (v.size() + block_size - 1) / block_size;

    multiply_kernel<T><<<num_blocks, block_size>>>(v.data(), v.size(), alpha);
    cudaDeviceSynchronize();
}

// Device function to perform warp-level reduction.
template <typename T, unsigned int blockSize>
__device__ void warpReduceAdd(volatile T* sdata, unsigned int tid)
{
    if (blockSize >= 64)
        sdata[tid] += sdata[tid + 32];
    if (blockSize >= 32)
        sdata[tid] += sdata[tid + 16];
    if (blockSize >= 16)
        sdata[tid] += sdata[tid + 8];
    if (blockSize >= 8)
        sdata[tid] += sdata[tid + 4];
    if (blockSize >= 4)
        sdata[tid] += sdata[tid + 2];
    if (blockSize >= 2)
        sdata[tid] += sdata[tid + 1];
}
// Device function to perform warp-level reduction.
template <typename T, unsigned int blockSize>
__device__ void warpReduceMax(volatile T* sdata, unsigned int tid)
{
    if (blockSize >= 64)
        sdata[tid] = max(sdata[tid], sdata[tid + 32]);
    if (blockSize >= 32)
        sdata[tid] = max(sdata[tid], sdata[tid + 16]);
    if (blockSize >= 16)
        sdata[tid] = max(sdata[tid], sdata[tid + 8]);
    if (blockSize >= 8)
        sdata[tid] = max(sdata[tid], sdata[tid + 4]);
    if (blockSize >= 4)
        sdata[tid] = max(sdata[tid], sdata[tid + 2]);
    if (blockSize >= 2)
        sdata[tid] = max(sdata[tid], sdata[tid + 1]);
}

template <typename T, unsigned int blockSize>
__global__ void dot_product_kernel(const T* lhs, const T* rhs, int size, T* result)
{
    // Shared memory allocation
    extern __shared__ __align__(sizeof(double)) unsigned char my_sdata[];
    T* sdata = reinterpret_cast<T*>(my_sdata);

    unsigned int tid      = threadIdx.x;
    unsigned int i        = blockIdx.x * (blockSize * 2) + tid;
    unsigned int gridSize = blockSize * 2 * gridDim.x;

    // Initialize thread-local sum
    T temp = 0;

    // Grid-stride loop to accumulate partial results
    while (i < size) {
        temp += lhs[i] * rhs[i];
        if (i + blockSize < size) {
            temp += lhs[i + blockSize] * rhs[i + blockSize];
        }
        i += gridSize;
    }

    // Store thread-local sum in shared memory
    sdata[tid] = temp;
    __syncthreads();

    // Reduction within the block
    if (blockSize >= 512) {
        if (tid < 256)
            sdata[tid] += sdata[tid + 256];
        __syncthreads();
    }
    if (blockSize >= 256) {
        if (tid < 128)
            sdata[tid] += sdata[tid + 128];
        __syncthreads();
    }
    if (blockSize >= 128) {
        if (tid < 64)
            sdata[tid] += sdata[tid + 64];
        __syncthreads();
    }

    // Final warp reduction
    if (tid < 32)
        warpReduceAdd<T, blockSize>(sdata, tid);

    // Accumulate the result of this block into the global result using atomicAdd
    if (tid == 0) {
        atomicAdd(result, sdata[0]);
    }
}

// Host function to compute the dot product of two vectors
template <typename T>
T dot_product(const GPU_Vector<T>& lhs, const GPU_Vector<T>& rhs)
{
    if (lhs.size() != rhs.size()) {
        throw std::invalid_argument("Vectors must be of the same size.");
    }

    int size      = lhs.size();
    T* device_lhs = lhs.data();
    T* device_rhs = rhs.data();

    const int threads_per_block = 256;
    int blocks                  = (size + threads_per_block * 2 - 1) / (threads_per_block * 2);

    // Allocate memory for the global result on the device
    T* device_result;
    cudaMalloc((void**)&device_result, sizeof(T));

    // Initialize the device result to 0
    T zero = 0;
    cudaMemcpy(device_result, &zero, sizeof(T), cudaMemcpyHostToDevice);

    // Launch the kernel
    size_t shared_memory_size = threads_per_block * sizeof(T);
    dot_product_kernel<T, threads_per_block>
        <<<blocks, threads_per_block, shared_memory_size>>>(device_lhs, device_rhs, size, device_result);
    cudaDeviceSynchronize();

    // Copy the final result back to the host
    T result;
    cudaMemcpy(&result, device_result, sizeof(T), cudaMemcpyDeviceToHost);

    // Cleanup
    cudaFree(device_result);

    return result;
}

// Kernel function for computing the L1 norm of a vector
template <typename T, unsigned int blockSize>
__global__ void l1_norm_kernel(const T* x, int size, T* result)
{
    // Shared memory allocation
    extern __shared__ __align__(sizeof(double)) unsigned char my_sdata[];
    T* sdata = reinterpret_cast<T*>(my_sdata);

    unsigned int tid      = threadIdx.x;
    unsigned int i        = blockIdx.x * (blockSize * 2) + tid;
    unsigned int gridSize = blockSize * 2 * gridDim.x;

    // Initialize thread-local sum
    T temp = 0;

    // Grid-stride loop to accumulate partial absolute values
    while (i < size) {
        temp += abs(x[i]);
        if (i + blockSize < size) {
            temp += abs(x[i + blockSize]);
        }
        i += gridSize;
    }

    // Store thread-local sum in shared memory
    sdata[tid] = temp;
    __syncthreads();

    // Reduction within the block
    if (blockSize >= 512) {
        if (tid < 256)
            sdata[tid] += sdata[tid + 256];
        __syncthreads();
    }
    if (blockSize >= 256) {
        if (tid < 128)
            sdata[tid] += sdata[tid + 128];
        __syncthreads();
    }
    if (blockSize >= 128) {
        if (tid < 64)
            sdata[tid] += sdata[tid + 64];
        __syncthreads();
    }

    // Final warp reduction
    if (tid < 32)
        warpReduceAdd<T, blockSize>(sdata, tid);

    // Accumulate the result of this block into the global result using atomicAdd
    if (tid == 0) {
        atomicAdd(result, sdata[0]);
    }
}

// Host function for computing the L1 norm of a vector
template <typename T>
T l1_norm(const GPU_Vector<T>& x)
{
    int size    = x.size();
    T* device_x = x.data();

    const int threads_per_block = 256;
    int blocks                  = (size + threads_per_block * 2 - 1) / (threads_per_block * 2);

    // Allocate memory for the global result on the device
    T* device_result;
    cudaMalloc((void**)&device_result, sizeof(T));

    // Initialize the device result to 0
    T zero = 0;
    cudaMemcpy(device_result, &zero, sizeof(T), cudaMemcpyHostToDevice);

    // Launch the kernel
    size_t shared_memory_size = threads_per_block * sizeof(T);
    l1_norm_kernel<T, threads_per_block>
        <<<blocks, threads_per_block, shared_memory_size>>>(device_x, size, device_result);
    cudaDeviceSynchronize();

    // Copy the final result back to the host
    T result;
    cudaMemcpy(&result, device_result, sizeof(T), cudaMemcpyDeviceToHost);

    // Cleanup
    cudaFree(device_result);

    return result;
}

// Kernel function for computing the L2 norm squared of a vector
template <typename T, unsigned int blockSize>
__global__ void l2_norm_squared_kernel(const T* x, int size, T* result)
{
    // Shared memory allocation
    extern __shared__ __align__(sizeof(double)) unsigned char my_sdata[];
    T* sdata = reinterpret_cast<T*>(my_sdata);

    unsigned int tid      = threadIdx.x;
    unsigned int i        = blockIdx.x * (blockSize * 2) + tid;
    unsigned int gridSize = blockSize * 2 * gridDim.x;

    // Initialize thread-local sum
    T temp = 0;

    // Grid-stride loop to accumulate partial squares
    while (i < size) {
        temp += x[i] * x[i];
        if (i + blockSize < size) {
            temp += x[i + blockSize] * x[i + blockSize];
        }
        i += gridSize;
    }

    // Store thread-local sum in shared memory
    sdata[tid] = temp;
    __syncthreads();

    // Reduction within the block
    if (blockSize >= 512) {
        if (tid < 256)
            sdata[tid] += sdata[tid + 256];
        __syncthreads();
    }
    if (blockSize >= 256) {
        if (tid < 128)
            sdata[tid] += sdata[tid + 128];
        __syncthreads();
    }
    if (blockSize >= 128) {
        if (tid < 64)
            sdata[tid] += sdata[tid + 64];
        __syncthreads();
    }

    // Final warp reduction
    if (tid < 32)
        warpReduceAdd<T, blockSize>(sdata, tid);

    // Accumulate the result of this block into the global result using atomicAdd
    if (tid == 0) {
        atomicAdd(result, sdata[0]);
    }
}

// Host function for computing the L2 norm squared of a vector
template <typename T>
T l2_norm_squared(const GPU_Vector<T>& x)
{
    int size    = x.size();
    T* device_x = x.data();

    const int threads_per_block = 256;
    int blocks                  = (size + threads_per_block * 2 - 1) / (threads_per_block * 2);

    // Allocate memory for the global result on the device
    T* device_result;
    cudaMalloc((void**)&device_result, sizeof(T));

    // Initialize the device result to 0
    T zero = 0;
    cudaMemcpy(device_result, &zero, sizeof(T), cudaMemcpyHostToDevice);

    // Launch the kernel
    size_t shared_memory_size = threads_per_block * sizeof(T);
    l2_norm_squared_kernel<T, threads_per_block>
        <<<blocks, threads_per_block, shared_memory_size>>>(device_x, size, device_result);
    cudaDeviceSynchronize();

    // Copy the final result back to the host
    T result;
    cudaMemcpy(&result, device_result, sizeof(T), cudaMemcpyDeviceToHost);

    // Cleanup
    cudaFree(device_result);

    return result;
}

// Host function for computing the L2 norm of a vector
template <typename T>
T l2_norm(const GPU_Vector<T>& x)
{
    T l2_squared = l2_norm_squared(x);
    return sqrt(l2_squared);
}

// Kernel function for computing the infinity norm of a vector
template <typename T, unsigned int blockSize>
__global__ void infinity_norm_kernel(const T* x, int size, T* result)
{
    // Shared memory allocation
    extern __shared__ __align__(sizeof(double)) unsigned char my_sdata[];
    T* sdata = reinterpret_cast<T*>(my_sdata);

    unsigned int tid      = threadIdx.x;
    unsigned int i        = blockIdx.x * (blockSize * 2) + tid;
    unsigned int gridSize = blockSize * 2 * gridDim.x;

    // Initialize thread-local maximum absolute value
    T temp = 0;

    // Grid-stride loop to accumulate the maximum absolute value
    while (i < size) {
        temp = max(temp, abs(x[i]));
        if (i + blockSize < size) {
            temp = max(temp, abs(x[i + blockSize]));
        }
        i += gridSize;
    }

    // Store thread-local maximum value in shared memory
    sdata[tid] = temp;
    __syncthreads();

    // Reduction within the block (find the maximum within the block)
    if (blockSize >= 512) {
        if (tid < 256)
            sdata[tid] = max(sdata[tid], sdata[tid + 256]);
        __syncthreads();
    }
    if (blockSize >= 256) {
        if (tid < 128)
            sdata[tid] = max(sdata[tid], sdata[tid + 128]);
        __syncthreads();
    }
    if (blockSize >= 128) {
        if (tid < 64)
            sdata[tid] = max(sdata[tid], sdata[tid + 64]);
        __syncthreads();
    }

    // Final warp reduction
    if (tid < 32)
        warpReduceMax<T, blockSize>(sdata, tid);

    // Store the result of the block reduction in shared memory
    if (tid == 0) {
        sdata[0] = sdata[0]; // block maximum
    }

    __syncthreads();

    // Final reduction step (after kernel execution)
    if (tid == 0) {
        // Here we will store the final maximum back to global memory using the result pointer
        T* global_result = result + blockIdx.x;
        *global_result   = sdata[0];
    }
}

// Host function for computing the infinity norm of a vector
template <typename T>
T infinity_norm(const GPU_Vector<T>& x)
{
    int size    = x.size();
    T* device_x = x.data();

    const int threads_per_block = 256;
    int blocks                  = (size + threads_per_block * 2 - 1) / (threads_per_block * 2);

    // Allocate memory for the intermediate block results on the device
    T* device_block_results;
    cudaMalloc((void**)&device_block_results, sizeof(T) * blocks);

    // Launch the kernel
    size_t shared_memory_size = threads_per_block * sizeof(T);
    infinity_norm_kernel<T, threads_per_block>
        <<<blocks, threads_per_block, shared_memory_size>>>(device_x, size, device_block_results);
    cudaDeviceSynchronize();

    // Perform the final reduction on the host
    T* host_block_results = new T[blocks];
    cudaMemcpy(host_block_results, device_block_results, sizeof(T) * blocks, cudaMemcpyDeviceToHost);

    T max_value = host_block_results[0];
    for (int i = 1; i < blocks; i++) {
        max_value = max(max_value, host_block_results[i]);
    }

    // Cleanup
    delete[] host_block_results;
    cudaFree(device_block_results);

    return max_value;
}
