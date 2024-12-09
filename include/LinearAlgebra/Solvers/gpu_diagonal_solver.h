#pragma once

#include <cuda_runtime.h>
#include <iostream>
#include <stdexcept>
#include <cassert>

template <typename T>
class GPU_DiagonalSolver
{
public:
    GPU_DiagonalSolver()
        : matrix_dimension_(0), diagonal_values_(nullptr) {}

    GPU_DiagonalSolver(const GPU_DiagonalSolver& other)
        : matrix_dimension_(other.matrix_dimension_)
    {
        cudaMalloc(&diagonal_values_, matrix_dimension_ * sizeof(T));
        cudaMemcpy(diagonal_values_, other.diagonal_values_,
                   matrix_dimension_ * sizeof(T), cudaMemcpyDeviceToDevice);
    }

    GPU_DiagonalSolver(GPU_DiagonalSolver&& other) noexcept
        : matrix_dimension_(other.matrix_dimension_), diagonal_values_(other.diagonal_values_)
    {
        other.matrix_dimension_ = 0;
        other.diagonal_values_ = nullptr;
    }

    explicit GPU_DiagonalSolver(const int matrix_dimension)
        : matrix_dimension_(matrix_dimension)
    {
        if (matrix_dimension_ <= 0)
            throw std::invalid_argument("Matrix dimension must be positive.");
        cudaMalloc(&diagonal_values_, matrix_dimension_ * sizeof(T));
    }

    GPU_DiagonalSolver& operator=(const GPU_DiagonalSolver& other)
    {
        if (this == &other) return *this;

        if (diagonal_values_) cudaFree(diagonal_values_);

        matrix_dimension_ = other.matrix_dimension_;
        cudaMalloc(&diagonal_values_, matrix_dimension_ * sizeof(T));
        cudaMemcpy(diagonal_values_, other.diagonal_values_,
                   matrix_dimension_ * sizeof(T), cudaMemcpyDeviceToDevice);

        return *this;
    }

    GPU_DiagonalSolver& operator=(GPU_DiagonalSolver&& other) noexcept
    {
        if (this == &other) return *this;

        if (diagonal_values_) cudaFree(diagonal_values_);

        matrix_dimension_ = other.matrix_dimension_;
        diagonal_values_ = other.diagonal_values_;

        other.matrix_dimension_ = 0;
        other.diagonal_values_ = nullptr;

        return *this;
    }

    ~GPU_DiagonalSolver()
    {
        if (diagonal_values_)
            cudaFree(diagonal_values_);
    }

    __host__ __device__ __forceinline__ int rows() const { return matrix_dimension_; }
    __host__ __device__ __forceinline__ int columns() const { return matrix_dimension_; }

    __device__ __forceinline__ const T& diagonal(const int index) const
    {
        assert(0 <= index && index < matrix_dimension_);
        return diagonal_values_[index];
    }

    __device__ __forceinline__ T& diagonal(const int index)
    {
        assert(0 <= index && index < matrix_dimension_);
        return diagonal_values_[index];
    }

    __device__ void solveInPlace(T* sol_rhs) const
    {
        int idx = threadIdx.x + blockIdx.x * blockDim.x;
        if (idx < matrix_dimension_)
        {
            sol_rhs[idx] /= diagonal_values_[idx];
        }
    }

    template <typename U>
    friend std::ostream& operator<<(std::ostream& stream, const GPU_DiagonalSolver<U>& solver);

private:
    int matrix_dimension_;
    T* diagonal_values_;
};

template <typename T>
std::ostream& operator<<(std::ostream& stream, const GPU_DiagonalSolver<T>& solver)
{
    T* host_data = new T[solver.matrix_dimension_];
    cudaMemcpy(host_data, solver.diagonal_values_,
               solver.matrix_dimension_ * sizeof(T), cudaMemcpyDeviceToHost);

    stream << "DiagonalSolver (dimension: " << solver.matrix_dimension_ << "): [";
    for (int i = 0; i < solver.matrix_dimension_; ++i)
    {
        stream << host_data[i];
        if (i != solver.matrix_dimension_ - 1)
            stream << ", ";
    }
    stream << "]";

    delete[] host_data;
    return stream;
}