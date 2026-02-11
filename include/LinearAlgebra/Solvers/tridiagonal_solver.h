#pragma once

#include <Kokkos_Core.hpp>

#include "../vector.h"
#include "../vector_operations.h"

template <typename T>
class BatchedTridiagonalSolver
{
public:
    BatchedTridiagonalSolver(int matrix_dimension, int batch_count, bool is_cyclic = true)
        : matrix_dimension_(matrix_dimension)
        , batch_count_(batch_count)
        , main_diagonal_("BatchedTridiagonalSolver::main_diagonal", matrix_dimension * batch_count)
        , sub_diagonal_("BatchedTridiagonalSolver::sub_diagonal", matrix_dimension * batch_count)
        , buffer_("BatchedTridiagonalSolver::buffer", is_cyclic ? matrix_dimension * batch_count : 0)
        , gamma_("BatchedTridiagonalSolver::gamma", is_cyclic ? batch_count : 0)
        , is_cyclic_(is_cyclic)
        , is_factorized_(false)
    {
        assign(main_diagonal_, T(0));
        assign(sub_diagonal_, T(0));
    }

    /* ------------------- */
    /* Accessors for sizes */
    /* ------------------- */

    int matrixDimension() const
    {
        return matrix_dimension_;
    }

    int batchCount() const
    {
        return batch_count_;
    }

    /* ---------------------------- */
    /* Accessors for matrix entries */
    /* ---------------------------- */

    KOKKOS_INLINE_FUNCTION
    const T& main_diagonal(const int batch_idx, const int index) const
    {
        return main_diagonal_(batch_idx * matrix_dimension_ + index);
    }
    KOKKOS_INLINE_FUNCTION
    T& main_diagonal(const int batch_idx, const int index)
    {
        return main_diagonal_(batch_idx * matrix_dimension_ + index);
    }

    KOKKOS_INLINE_FUNCTION
    const T& sub_diagonal(const int batch_idx, const int index) const
    {
        return sub_diagonal_(batch_idx * matrix_dimension_ + index);
    }
    KOKKOS_INLINE_FUNCTION
    T& sub_diagonal(const int batch_idx, const int index)
    {
        return sub_diagonal_(batch_idx * matrix_dimension_ + index);
    }

    KOKKOS_INLINE_FUNCTION
    const T& cyclic_corner(const int batch_idx) const
    {
        return sub_diagonal_(batch_idx * matrix_dimension_ + (matrix_dimension_ - 1));
    }

    KOKKOS_INLINE_FUNCTION
    T& cyclic_corner(const int batch_idx)
    {
        return sub_diagonal_(batch_idx * matrix_dimension_ + (matrix_dimension_ - 1));
    }

    /* ---------------------------------------------- */
    /* Setup: Cholesky Decomposition: A = L * D * L^T */
    /* ---------------------------------------------- */
    // This step factorizes the tridiagonal matrix into lower  triangular (L) and diagonal (D) matrices.
    // For cyclic systems, it also applies the Shermann-Morrison adjustment to account for the cyclic connection.

    void setup()
    {
        // Create local copies for lambda capture
        int matrix_dimension    = matrix_dimension_;
        Vector<T> main_diagonal = main_diagonal_;
        Vector<T> sub_diagonal  = sub_diagonal_;
        Vector<T> gamma         = gamma_;

        if (!is_cyclic_) {
            Kokkos::parallel_for(
                "SetupNonCyclic", batch_count_, KOKKOS_LAMBDA(const int batch_idx) {
                    // ----------------------------------- //
                    // Obtain offset for the current batch //
                    int offset = batch_idx * matrix_dimension;

                    // ---------------------- //
                    // Cholesky Decomposition //
                    for (int i = 1; i < matrix_dimension; i++) {
                        sub_diagonal(offset + i - 1) /= main_diagonal(offset + i - 1);
                        const T factor = sub_diagonal(offset + i - 1);
                        main_diagonal(offset + i) -= factor * factor * main_diagonal(offset + i - 1);
                    }
                });
        }
        else {
            Kokkos::parallel_for(
                "SetupCyclic", batch_count_, KOKKOS_LAMBDA(const int batch_idx) {
                    // ----------------------------------- //
                    // Obtain offset for the current batch //
                    int offset = batch_idx * matrix_dimension;

                    // ------------------------------------------------- //
                    // Shermann-Morrison Adjustment                      //
                    // - Modify the first and last main diagonal element //
                    // - Compute and store gamma for later use           //
                    // ------------------------------------------------- //
                    T cyclic_corner_element = sub_diagonal(offset + matrix_dimension - 1);
                    gamma(batch_idx)        = -main_diagonal(offset + 0);
                    main_diagonal(offset + 0) -= gamma(batch_idx);
                    main_diagonal(offset + matrix_dimension - 1) -=
                        cyclic_corner_element * cyclic_corner_element / gamma(batch_idx);

                    // ---------------------- //
                    // Cholesky Decomposition //
                    for (int i = 1; i < matrix_dimension; i++) {
                        sub_diagonal(offset + i - 1) /= main_diagonal(offset + i - 1);
                        const T factor = sub_diagonal(offset + i - 1);
                        main_diagonal(offset + i) -= factor * factor * main_diagonal(offset + i - 1);
                    }
                });
        }
        Kokkos::fence();
        is_factorized_ = true;
    }

    /* ---------------------------------------- */
    /* Solve: Forward and Backward Substitution */
    /* ---------------------------------------- */
    // This step solves the system Ax = b using the factorized form of A.
    // For cyclic systems, it also performs the Shermann-Morrison reconstruction to obtain the final solution.

    void solve(Vector<T> rhs, int batch_offset = 0, int batch_stride = 1)
    {
        if (!is_factorized_) {
            throw std::runtime_error("Error: Matrix must be factorized before solving.");
        }

        // Compute the effective number of batches to solve
        int effective_batch_count = (batch_count_ - batch_offset + batch_stride - 1) / batch_stride;

        // Create local copies for lambda capture
        int matrix_dimension    = matrix_dimension_;
        Vector<T> main_diagonal = main_diagonal_;
        Vector<T> sub_diagonal  = sub_diagonal_;
        Vector<T> buffer        = buffer_;
        Vector<T> gamma         = gamma_;

        if (!is_cyclic_) {
            Kokkos::parallel_for(
                "SolveNonCyclic", effective_batch_count, KOKKOS_LAMBDA(const int k) {
                    // ----------------------------------- //
                    // Obtain offset for the current batch //
                    int batch_idx = batch_stride * k + batch_offset;
                    int offset    = batch_idx * matrix_dimension;

                    // -------------------- //
                    // Forward Substitution //
                    for (int i = 1; i < matrix_dimension; i++) {
                        rhs(offset + i) -= sub_diagonal(offset + i - 1) * rhs(offset + i - 1);
                    }

                    // ---------------- //
                    // Diagonal Scaling //
                    for (int i = 0; i < matrix_dimension; i++) {
                        rhs(offset + i) /= main_diagonal(offset + i);
                    }

                    // --------------------- //
                    // Backward Substitution //
                    for (int i = matrix_dimension - 2; i >= 0; i--) {
                        rhs(offset + i) -= sub_diagonal(offset + i) * rhs(offset + i + 1);
                    }
                });
        }
        else {
            Kokkos::parallel_for(
                "SolveCyclic", effective_batch_count, KOKKOS_LAMBDA(const int k) {
                    // ----------------------------------- //
                    // Obtain offset for the current batch //
                    int batch_idx = batch_stride * k + batch_offset;
                    int offset    = batch_idx * matrix_dimension;

                    // -------------------- //
                    // Forward Substitution //
                    T cyclic_corner_element = sub_diagonal(offset + matrix_dimension - 1);
                    buffer(offset + 0)      = gamma(batch_idx);
                    for (int i = 1; i < matrix_dimension; i++) {
                        rhs(offset + i) -= sub_diagonal(offset + i - 1) * rhs(offset + i - 1);
                        if (i < matrix_dimension - 1)
                            buffer(offset + i) = 0.0 - sub_diagonal(offset + i - 1) * buffer(offset + i - 1);
                        else
                            buffer(offset + i) =
                                cyclic_corner_element - sub_diagonal(offset + i - 1) * buffer(offset + i - 1);
                    }

                    // ---------------- //
                    // Diagonal Scaling //
                    for (int i = 0; i < matrix_dimension; i++) {
                        rhs(offset + i) /= main_diagonal(offset + i);
                        buffer(offset + i) /= main_diagonal(offset + i);
                    }

                    // --------------------- //
                    // Backward Substitution //
                    for (int i = matrix_dimension - 2; i >= 0; i--) {
                        rhs(offset + i) -= sub_diagonal(offset + i) * rhs(offset + i + 1);
                        buffer(offset + i) -= sub_diagonal(offset + i) * buffer(offset + i + 1);
                    }

                    // ------------------------------- //
                    // Shermann-Morrison Reonstruction //
                    const T dot_product_x_v =
                        rhs(offset + 0) + cyclic_corner_element / gamma(batch_idx) * rhs(offset + matrix_dimension - 1);
                    const T dot_product_u_v = buffer(offset + 0) + cyclic_corner_element / gamma(batch_idx) *
                                                                       buffer(offset + matrix_dimension - 1);
                    const T factor = dot_product_x_v / (1.0 + dot_product_u_v);

                    for (int i = 0; i < matrix_dimension; i++) {
                        rhs(offset + i) -= factor * buffer(offset + i);
                    }
                });
        }
        Kokkos::fence();
    }

    /* ---------------------------- */
    /* Solve: Diagonal Scaling Only */
    /* ---------------------------- */
    // This step performs only the diagonal scaling part of the solve process.
    // It is useful when the matrix has a non-zero diagonal but zero off-diagonal entries.
    // Note that .setup() modifies main_diagonal(0) in the cyclic case.

    void solve_diagonal(Vector<T> rhs, int batch_offset = 0, int batch_stride = 1)
    {
        if (!is_factorized_) {
            throw std::runtime_error("Error: Matrix must be factorized before solving.");
        }

        // Compute the effective number of batches to solve
        int effective_batch_count = (batch_count_ - batch_offset + batch_stride - 1) / batch_stride;

        // Create local copies for lambda capture
        int matrix_dimension    = matrix_dimension_;
        Vector<T> main_diagonal = main_diagonal_;
        Vector<T> gamma         = gamma_;

        if (!is_cyclic_) {
            Kokkos::parallel_for(
                "SolveDiagonalNonCyclic", effective_batch_count, KOKKOS_LAMBDA(const int k) {
                    // ----------------------------------- //
                    // Obtain offset for the current batch //
                    int batch_idx = batch_stride * k + batch_offset;
                    int offset    = batch_idx * matrix_dimension;

                    // ---------------- //
                    // Diagonal Scaling //
                    for (int i = 0; i < matrix_dimension; i++) {
                        rhs(offset + i) /= main_diagonal(offset + i);
                    }
                });
        }
        else {
            Kokkos::parallel_for(
                "SolveDiagonalCyclic", effective_batch_count, KOKKOS_LAMBDA(const int k) {
                    // ----------------------------------- //
                    // Obtain offset for the current batch //
                    int batch_idx = batch_stride * k + batch_offset;
                    int offset    = batch_idx * matrix_dimension;

                    // ---------------- //
                    // Diagonal Scaling //
                    rhs(offset + 0) /= main_diagonal(offset + 0) + gamma(batch_idx);
                    for (int i = 1; i < matrix_dimension; i++) {
                        rhs(offset + i) /= main_diagonal(offset + i);
                    }
                });
        }
        Kokkos::fence();
    }

private:
    int matrix_dimension_;
    int batch_count_;

    Vector<T> main_diagonal_;
    Vector<T> sub_diagonal_;
    Vector<T> buffer_;
    Vector<T> gamma_;

    bool is_cyclic_;
    bool is_factorized_;
};