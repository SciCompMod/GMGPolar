#pragma once

#include <Kokkos_Core.hpp>

#include "../../LinearAlgebra/Vector/vector.h"
#include "../../LinearAlgebra/Vector/vector_operations.h"

namespace gmgpolar
{

template <typename T>
class BatchedTridiagonalSolverThomas : public BatchedTridiagonalSolverBase<T>
{
public:
    BatchedTridiagonalSolverThomas(int matrix_dimension, int batch_count, bool is_cyclic = true)
        : BatchedTridiagonalSolverBase<T>(matrix_dimension, batch_count, is_cyclic)
        , buffer_("BatchedTridiagonalSolverThomas::buffer", is_cyclic ? matrix_dimension * batch_count : 0)
        , gamma_("BatchedTridiagonalSolverThomas::gamma", is_cyclic ? batch_count : 0)
        , is_factorized_(false)
    {
        assign(buffer_, T(0));
        assign(gamma_, T(0));
    }

    /* ---------------------------------------------- */
    /* Setup: Cholesky Decomposition: A = L * D * L^T */
    /* ---------------------------------------------- */
    // This step factorizes the tridiagonal matrix into lower  triangular (L) and diagonal (D) matrices.
    // For cyclic systems, it also applies the Shermann-Morrison adjustment to account for the cyclic connection.

    void setup()
    {
        // Create local copies for lambda capture
        int matrix_dimension    = this->matrix_dimension_;
        Vector<T> main_diagonal = this->main_diagonal_;
        Vector<T> sub_diagonal  = this->sub_diagonal_;
        Vector<T> gamma         = gamma_;

        if (!this->is_cyclic_) {
            Kokkos::parallel_for(
                "SetupNonCyclic", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, this->batch_count_),
                KOKKOS_LAMBDA(const int batch_idx) {
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
                "SetupCyclic", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, this->batch_count_),
                KOKKOS_LAMBDA(const int batch_idx) {
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
        int effective_batch_count = (this->batch_count_ - batch_offset + batch_stride - 1) / batch_stride;

        // Create local copies for lambda capture
        int matrix_dimension    = this->matrix_dimension_;
        Vector<T> main_diagonal = this->main_diagonal_;
        Vector<T> sub_diagonal  = this->sub_diagonal_;
        Vector<T> buffer        = buffer_;
        Vector<T> gamma         = gamma_;

        if (!this->is_cyclic_) {
            Kokkos::parallel_for(
                "SolveNonCyclic", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, effective_batch_count),
                KOKKOS_LAMBDA(const int k) {
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
                "SolveCyclic", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, effective_batch_count),
                KOKKOS_LAMBDA(const int k) {
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
        int effective_batch_count = (this->batch_count_ - batch_offset + batch_stride - 1) / batch_stride;

        // Create local copies for lambda capture
        int matrix_dimension    = this->matrix_dimension_;
        Vector<T> main_diagonal = this->main_diagonal_;
        Vector<T> gamma         = gamma_;

        if (!this->is_cyclic_) {
            Kokkos::parallel_for(
                "SolveDiagonalNonCyclic", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, effective_batch_count),
                KOKKOS_LAMBDA(const int k) {
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
                "SolveDiagonalCyclic", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, effective_batch_count),
                KOKKOS_LAMBDA(const int k) {
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
    Vector<T> buffer_;
    Vector<T> gamma_;

    bool is_factorized_;
};
} // namespace gmgpolar
