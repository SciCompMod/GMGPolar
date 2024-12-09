#pragma once

#include "vector.h"
#include "gpu_vector.h"

#include <cuda_runtime.h>
#include <iostream>
#include <stdexcept>

template <typename T>
class GPU_Vector
{
public:
    GPU_Vector();
    GPU_Vector(const GPU_Vector& other);
    GPU_Vector(GPU_Vector&& other) noexcept;

    explicit GPU_Vector(int size);
    ~GPU_Vector();

    GPU_Vector& operator=(const GPU_Vector& other);
    GPU_Vector& operator=(GPU_Vector&& other) noexcept;

    __host__ __device__ int size() const noexcept;
    __host__ __device__ T* data() const noexcept;

    template <typename U>
    friend std::ostream& operator<<(std::ostream& stream, const GPU_Vector<U>& gpu_vector);

private:
    int size_;
    T* device_values_;

    void checkCudaError(cudaError_t error, const char* errorMessage) const;
};

// Copy data from host (Vector<T>) to device (GPU_Vector<T>)
template <typename T>
void copyHostToDevice(const Vector<T>& host_vector, GPU_Vector<T>& gpu_vector)
{
    if (host_vector.size() != gpu_vector.size()) {
        throw std::invalid_argument("Vector sizes must match for data copy.");
    }
    cudaMemcpy(gpu_vector.data(), host_vector.begin(), host_vector.size() * sizeof(T), cudaMemcpyHostToDevice);
}

// Copy data from device (GPU_Vector<T>) to host (Vector<T>)
template <typename T>
void copyDeviceToHost(const GPU_Vector<T>& gpu_vector, Vector<T>& host_vector)
{
    if (gpu_vector.size() != host_vector.size()) {
        throw std::invalid_argument("Vector sizes must match for data copy.");
    }
    cudaMemcpy(host_vector.begin(), gpu_vector.data(), gpu_vector.size() * sizeof(T), cudaMemcpyDeviceToHost);
}

// Default constructor
template <typename T>
GPU_Vector<T>::GPU_Vector()
    : size_(0)
    , device_values_(nullptr)
{
}

// Constructor with size
template <typename T>
GPU_Vector<T>::GPU_Vector(int size)
    : size_(size)
    , device_values_(nullptr)
{
    cudaError_t error = cudaMalloc(&device_values_, size_ * sizeof(T));
    checkCudaError(error, "Failed to allocate memory on the GPU.");
}

// Destructor
template <typename T>
GPU_Vector<T>::~GPU_Vector()
{
    if (device_values_) {
        cudaError_t error = cudaFree(device_values_);
        checkCudaError(error, "Failed to free memory on the GPU.");
    }
}

// Copy constructor
template <typename T>
GPU_Vector<T>::GPU_Vector(const GPU_Vector& other)
    : size_(other.size_)
{
    cudaError_t error = cudaMalloc(&device_values_, size_ * sizeof(T));
    checkCudaError(error, "Failed to allocate memory on the GPU.");

    error = cudaMemcpy(device_values_, other.device_values_, size_ * sizeof(T), cudaMemcpyDeviceToDevice);
    checkCudaError(error, "Failed to copy data to the GPU.");
}

// Move constructor
template <typename T>
GPU_Vector<T>::GPU_Vector(GPU_Vector&& other) noexcept
    : size_(other.size_)
    , device_values_(other.device_values_)
{
    other.size_          = 0;
    other.device_values_ = nullptr;
}

// Copy assignment operator
template <typename T>
GPU_Vector<T>& GPU_Vector<T>::operator=(const GPU_Vector& other)
{
    if (this == &other)
        return *this;

    if (size_ != other.size_) {
        size_ = other.size_;
        cudaFree(device_values_);
        cudaError_t error = cudaMalloc(&device_values_, size_ * sizeof(T));
        checkCudaError(error, "Failed to allocate memory on the GPU.");
    }

    cudaError_t error = cudaMemcpy(device_values_, other.device_values_, size_ * sizeof(T), cudaMemcpyDeviceToDevice);
    checkCudaError(error, "Failed to copy data to the GPU.");

    return *this;
}

// Move assignment operator
template <typename T>
GPU_Vector<T>& GPU_Vector<T>::operator=(GPU_Vector&& other) noexcept
{
    if (this != &other) {
        cudaFree(device_values_);
        size_                = other.size_;
        device_values_       = other.device_values_;
        other.size_          = 0;
        other.device_values_ = nullptr;
    }

    return *this;
}

// Get size of the vector
template <typename T>
__host__ __device__ int GPU_Vector<T>::size() const noexcept
{
    return size_;
}

// Get device pointer to the data
template <typename T>
__host__ __device__ T* GPU_Vector<T>::data() const noexcept
{
    return device_values_;
}

// Helper function to check CUDA errors
template <typename T>
void GPU_Vector<T>::checkCudaError(cudaError_t error, const char* errorMessage) const
{
    if (error != cudaSuccess) {
        std::cerr << errorMessage << ": " << cudaGetErrorString(error) << std::endl;
        throw std::runtime_error("CUDA error occurred.");
    }
}

template <typename U>
std::ostream& operator<<(std::ostream& stream, const GPU_Vector<U>& gpu_vector)
{
    Vector<U> host_vector(gpu_vector.size());
    copyDeviceToHost(gpu_vector, host_vector);
    return stream << host_vector;
}
