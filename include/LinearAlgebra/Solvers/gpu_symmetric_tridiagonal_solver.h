#pragma once

#include <cassert>
#include <iostream>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

template <typename T>
class GPU_SymmetricTridiagonalSolver
{
public:
    GPU_SymmetricTridiagonalSolver()
        : matrix_dimension_(0), main_diagonal_values_(nullptr), sub_diagonal_values_(nullptr), is_cyclic_(true), factorized_(false), gamma_(0.0) {}

    GPU_SymmetricTridiagonalSolver(const GPU_SymmetricTridiagonalSolver& other) = delete;
    GPU_SymmetricTridiagonalSolver& operator=(const GPU_SymmetricTridiagonalSolver& other) = delete;

    GPU_SymmetricTridiagonalSolver(GPU_SymmetricTridiagonalSolver&& other) noexcept
        : matrix_dimension_(other.matrix_dimension_),
          main_diagonal_values_(other.main_diagonal_values_),
          sub_diagonal_values_(other.sub_diagonal_values_),
          cyclic_corner_element_(other.cyclic_corner_element_),
          is_cyclic_(other.is_cyclic_),
          factorized_(other.factorized_),
          gamma_(other.gamma_) {
        other.matrix_dimension_ = 0;
        other.main_diagonal_values_ = nullptr;
        other.sub_diagonal_values_ = nullptr;
    }

    GPU_SymmetricTridiagonalSolver& operator=(GPU_SymmetricTridiagonalSolver&& other) noexcept {
        matrix_dimension_ = other.matrix_dimension_;
        main_diagonal_values_ = other.main_diagonal_values_;
        sub_diagonal_values_ = other.sub_diagonal_values_;
        cyclic_corner_element_ = other.cyclic_corner_element_;
        is_cyclic_ = other.is_cyclic_;
        factorized_ = other.factorized_;
        gamma_ = other.gamma_;

        other.matrix_dimension_ = 0;
        other.main_diagonal_values_ = nullptr;
        other.sub_diagonal_values_ = nullptr;
        return *this;
    }

    explicit GPU_SymmetricTridiagonalSolver(const int matrix_dimension)
        : matrix_dimension_(matrix_dimension),
          main_diagonal_values_(nullptr),
          sub_diagonal_values_(nullptr),
          is_cyclic_(true),
          factorized_(false),
          gamma_(0.0) {
        cudaMalloc(&main_diagonal_values_, matrix_dimension_ * sizeof(T));
        cudaMalloc(&sub_diagonal_values_, (matrix_dimension_ - 1) * sizeof(T));
    }

    ~GPU_SymmetricTridiagonalSolver() {
        if (main_diagonal_values_) cudaFree(main_diagonal_values_);
        if (sub_diagonal_values_) cudaFree(sub_diagonal_values_);
    }

    __host__ __device__ __forceinline__ void is_cyclic(bool value) { is_cyclic_ = value; }
    __host__ __device__ __forceinline__ bool is_cyclic() const { return is_cyclic_; }

    __host__ __device__ __forceinline__ int rows() const { return matrix_dimension_; }
    __host__ __device__ __forceinline__ int columns() const { return matrix_dimension_; }

    __device__ __forceinline__ const T& main_diagonal(const int index) const {
        assert(index >= 0 && index < matrix_dimension_);
        return main_diagonal_values_[index];
    }

    __device__ __forceinline__ T& main_diagonal(const int index) {
        assert(index >= 0 && index < matrix_dimension_);
        return main_diagonal_values_[index];
    }

    __device__ __forceinline__ const T& sub_diagonal(const int index) const {
        assert(index >= 0 && index < matrix_dimension_ - 1);
        return sub_diagonal_values_[index];
    }

    __device__ __forceinline__ T& sub_diagonal(const int index) {
        assert(index >= 0 && index < matrix_dimension_ - 1);
        return sub_diagonal_values_[index];
    }

    __host__ __device__ __forceinline__ const T& cyclic_corner_element() const {
        assert(is_cyclic_);
        return cyclic_corner_element_;
    }

    __device__ __forceinline__ T& cyclic_corner_element() {
        assert(is_cyclic_);
        return cyclic_corner_element_;
    }

    __device__ void solveInPlace(T* sol_rhs, T* temp = nullptr) {
        assert(matrix_dimension_ >= 2);
        if (is_cyclic_) {
            assert(temp != nullptr);
            solveSymmetricCyclicTridiagonal(sol_rhs, temp);
        } else {
            solveSymmetricTridiagonal(sol_rhs);
        }
    }

    __device__ void factorize(){
        if (!is_cyclic_) {
            if (!factorized_) {
                for (int i = 1; i < matrix_dimension_; i++) {
                    assert(main_diagonal(i - 1) != 0.0);
                    sub_diagonal(i - 1) /= main_diagonal(i - 1);
                    main_diagonal(i) -= sub_diagonal(i - 1) * sub_diagonal(i - 1) * main_diagonal(i - 1);
                }
                factorized_ = true;
            }
        } else {
            if (!factorized_) {
                gamma_ = -main_diagonal(0);
                main_diagonal(0) -= gamma_;
                main_diagonal(matrix_dimension_ - 1) -= cyclic_corner_element() * cyclic_corner_element() / gamma_;

                for (int i = 1; i < matrix_dimension_; i++) {
                    sub_diagonal(i - 1) /= main_diagonal(i - 1);
                    main_diagonal(i) -= sub_diagonal(i - 1) * sub_diagonal(i - 1) * main_diagonal(i - 1);
                }
                factorized_ = true;
            }
        }
    }


    template <typename U>
    friend std::ostream& operator<<(std::ostream& stream, const GPU_SymmetricTridiagonalSolver<U>& solver);

private:
    int matrix_dimension_;
    T* main_diagonal_values_;
    T* sub_diagonal_values_;
    T cyclic_corner_element_ = 0.0;
    bool is_cyclic_ = true;
    bool factorized_ = false;
    T gamma_ = 0.0;

    __device__ void solveSymmetricTridiagonal(T* x) {
        assert(factorized_);

        for (int i = 1; i < matrix_dimension_; i++) {
            x[i] -= sub_diagonal(i - 1) * x[i - 1];
        }

        for (int i = 0; i < matrix_dimension_; i++) {
            assert(main_diagonal(i) != 0.0);
            x[i] /= main_diagonal(i);
        }

        for (int i = matrix_dimension_ - 2; i >= 0; i--) {
            x[i] -= sub_diagonal(i) * x[i + 1];
        }
    }

    __device__ void solveSymmetricCyclicTridiagonal(T* x, T* u) {
        assert(factorized_);

        u[0] = gamma_;
        for (int i = 1; i < matrix_dimension_; i++) {
            x[i] -= sub_diagonal(i - 1) * x[i - 1];
            if (i < matrix_dimension_ - 1)
                u[i] = -sub_diagonal(i - 1) * u[i - 1];
            else
                u[i] = cyclic_corner_element() - sub_diagonal(i - 1) * u[i - 1];
        }

        for (int i = 0; i < matrix_dimension_; i++) {
            x[i] /= main_diagonal(i);
            u[i] /= main_diagonal(i);
        }

        for (int i = matrix_dimension_ - 2; i >= 0; i--) {
            x[i] -= sub_diagonal(i) * x[i + 1];
            u[i] -= sub_diagonal(i) * u[i + 1];
        }

        const double dot_x_u = x[0] + cyclic_corner_element() / gamma_ * x[matrix_dimension_ - 1];
        const double dot_u_u = u[0] + cyclic_corner_element() / gamma_ * u[matrix_dimension_ - 1];
        const double factor = dot_x_u / (1.0 + dot_u_u);

        for (int i = 0; i < matrix_dimension_; i++) {
            x[i] -= factor * u[i];
        }
    }
};





// template <typename U>
// std::ostream& operator<<(std::ostream& stream, const GPU_SymmetricTridiagonalSolver<U>& solver)
// {

//     std::cout<<solver.cyclic_corner_element_<<std::endl;





//     stream << "Symmetric Tridiagonal Matrix (Dimension: " << solver.matrix_dimension_ << ")\n";

//     // Allocate memory on the host for the main diagonal and sub diagonal
//     U* h_main_diagonal = new U[solver.matrix_dimension_];
//     U* h_sub_diagonal = new U[solver.matrix_dimension_ - 1];

//     // Copy data from device to host
//     cudaMemcpy(h_main_diagonal, solver.main_diagonal_values_, solver.matrix_dimension_ * sizeof(U), cudaMemcpyDeviceToHost);
//     cudaMemcpy(h_sub_diagonal, solver.sub_diagonal_values_, (solver.matrix_dimension_ - 1) * sizeof(U), cudaMemcpyDeviceToHost);

//     if (solver.factorized_) {
//         // Print the L, D decomposition if factorized
//         stream << "L Factor (Sub Diagonal Elements): [";
//         for (int i = 0; i < solver.matrix_dimension_ - 1; ++i) {
//             stream << h_sub_diagonal[i];
//             if (i != solver.matrix_dimension_ - 2)
//                 stream << ", ";
//         }
//         stream << "]\n";

//         stream << "D Factor (Diagonal Elements): [";
//         for (int i = 0; i < solver.matrix_dimension_; ++i) {
//             stream << h_main_diagonal[i];
//             if (i != solver.matrix_dimension_ - 1)
//                 stream << ", ";
//         }
//         stream << "]\n";
//     }
//     else {
//         // Print the matrix in its tridiagonal form if not factorized
//         stream << "Main Diagonal: [";
//         for (int i = 0; i < solver.matrix_dimension_; ++i) {
//             stream << h_main_diagonal[i];
//             if (i != solver.matrix_dimension_ - 1)
//                 stream << ", ";
//         }
//         stream << "]\n";

//         stream << "Sub Diagonal: [";
//         for (int i = 0; i < solver.matrix_dimension_ - 1; ++i) {
//             stream << h_sub_diagonal[i];
//             if (i != solver.matrix_dimension_ - 2)
//                 stream << ", ";
//         }
//         stream << "]\n";

//         if (solver.is_cyclic_) {
//             // U cyclic_corner_value;
//             // cudaMemcpy(&cyclic_corner_value, cyclic_corner_element, sizeof(U), cudaMemcpyDeviceToHost);

//             stream << "Cyclic Corner Element: " << solver.cyclic_corner_element() << "\n";
//         }
//         else {
//             stream << "Matrix is not cyclic.\n";
//         }
//     }

//     // Free allocated host memory
//     delete[] h_main_diagonal;
//     delete[] h_sub_diagonal;

//     return stream;
// }