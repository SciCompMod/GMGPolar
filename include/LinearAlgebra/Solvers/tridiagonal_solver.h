#pragma once

#include <Kokkos_Core.hpp>

#include <stdexcept>

#include "../../LinearAlgebra/Vector/vector.h"
#include "../../LinearAlgebra/Vector/vector_operations.h"

namespace gmgpolar
{

template <typename T>
class BatchedTridiagonalSolverBase
{
public:
    BatchedTridiagonalSolverBase(int matrix_dimension, int batch_count, bool is_cyclic = true)
        : matrix_dimension_(matrix_dimension)
        , batch_count_(batch_count)
        , main_diagonal_("BatchedTridiagonalSolver::main_diagonal", matrix_dimension * batch_count)
        , sub_diagonal_("BatchedTridiagonalSolver::sub_diagonal", matrix_dimension * batch_count)
        , is_cyclic_(is_cyclic)
    {
        if (matrix_dimension_ <= 0) {
            throw std::invalid_argument("matrix_dimension must be positive");
        }

        if (batch_count_ < 0) {
            throw std::invalid_argument("batch_count must be non-negative");
        }

        assign(main_diagonal_, T(0));
        assign(sub_diagonal_, T(0));
    }

    virtual ~BatchedTridiagonalSolverBase() = default;

    virtual void setup() = 0;

    virtual void solve(Vector<T> rhs, int batch_offset = 0, int batch_stride = 1) = 0;

    virtual void solve_diagonal(Vector<T> rhs, int batch_offset = 0, int batch_stride = 1) = 0;

    KOKKOS_INLINE_FUNCTION int matrixDimension() const
    {
        return matrix_dimension_;
    }

    KOKKOS_INLINE_FUNCTION int batchCount() const
    {
        return batch_count_;
    }

    KOKKOS_INLINE_FUNCTION
    const T& main_diagonal(int batch_idx, int index) const
    {
        return main_diagonal_(batch_idx * matrix_dimension_ + index);
    }

    KOKKOS_INLINE_FUNCTION
    void set_main_diagonal(int batch_idx, int index, const T& value) const
    {
        main_diagonal_(batch_idx * matrix_dimension_ + index) = value;
    }

    KOKKOS_INLINE_FUNCTION
    void increase_main_diagonal(int batch_idx, int index, const T& value) const
    {
        main_diagonal_(batch_idx * matrix_dimension_ + index) += value;
    }

    KOKKOS_INLINE_FUNCTION
    const T& sub_diagonal(int batch_idx, int index) const
    {
        return sub_diagonal_(batch_idx * matrix_dimension_ + index);
    }

    KOKKOS_INLINE_FUNCTION
    void set_sub_diagonal(int batch_idx, int index, const T& value) const
    {
        sub_diagonal_(batch_idx * matrix_dimension_ + index) = value;
    }

    KOKKOS_INLINE_FUNCTION
    void increase_sub_diagonal(int batch_idx, int index, const T& value) const
    {
        sub_diagonal_(batch_idx * matrix_dimension_ + index) += value;
    }

    KOKKOS_INLINE_FUNCTION
    const T& cyclic_corner(int batch_idx) const
    {
        return sub_diagonal_(batch_idx * matrix_dimension_ + matrix_dimension_ - 1);
    }

    KOKKOS_INLINE_FUNCTION
    void set_cyclic_corner(int batch_idx, const T& value) const
    {
        sub_diagonal_(batch_idx * matrix_dimension_ + matrix_dimension_ - 1) = value;
    }

    KOKKOS_INLINE_FUNCTION
    void increase_cyclic_corner(int batch_idx, const T& value) const
    {
        sub_diagonal_(batch_idx * matrix_dimension_ + matrix_dimension_ - 1) += value;
    }

protected:
    int matrix_dimension_;
    int batch_count_;

    Vector<T> main_diagonal_;
    Vector<T> sub_diagonal_;

    bool is_cyclic_;
};

} // namespace gmgpolar

#include "tridiagonal_solver_thomas.h"
#include "tridiagonal_solver_pcr.h"
#include "tridiagonal_solver_cr.h"
#include "tridiagonal_solver_crpcr.h"

namespace gmgpolar
{

template <typename T>
using BatchedTridiagonalSolver = BatchedTridiagonalSolverCRPCR<T>;

/*
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP) || defined(KOKKOS_ENABLE_SYCL)

template <typename T>
using BatchedTridiagonalSolver = BatchedTridiagonalSolverPCR<T>;

// template <typename T>
// using BatchedTridiagonalSolver = BatchedTridiagonalSolverCRPCR<T>;

#else

template <typename T>
using BatchedTridiagonalSolver = BatchedTridiagonalSolverThomas<T>;

#endif

*/

} // namespace gmgpolar
