#pragma once

#include <Kokkos_Core.hpp>

#include "../vector.h"

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
        Kokkos::deep_copy(main_diagonal_, T(0));
        Kokkos::deep_copy(sub_diagonal_, T(0));
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

    struct SetupNonCyclic {
        int m_matrix_dimension;
        Vector<T> m_main_diagonal;
        Vector<T> m_sub_diagonal;

        void operator()(const int batch_idx) const
        {
            // ----------------------------------- //
            // Obtain offset for the current batch //
            int offset = batch_idx * m_matrix_dimension;

            // ---------------------- //
            // Cholesky Decomposition //
            for (int i = 1; i < m_matrix_dimension; i++) {
                m_sub_diagonal(offset + i - 1) /= m_main_diagonal(offset + i - 1);
                const T factor = m_sub_diagonal(offset + i - 1);
                m_main_diagonal(offset + i) -= factor * factor * m_main_diagonal(offset + i - 1);
            }
        }
    };

    struct SetupCyclic {
        int m_matrix_dimension;
        Vector<T> m_main_diagonal;
        Vector<T> m_sub_diagonal;
        Vector<T> m_gamma;

        void operator()(const int batch_idx) const
        {
            // ----------------------------------- //
            // Obtain offset for the current batch //
            int offset = batch_idx * m_matrix_dimension;

            // ------------------------------------------------- //
            // Shermann-Morrison Adjustment                      //
            // - Modify the first and last main diagonal element //
            // - Compute and store gamma for later use           //
            // ------------------------------------------------- //
            T cyclic_corner_element = m_sub_diagonal(offset + m_matrix_dimension - 1);
            /* gamma_ = -main_diagonal(0);*/
            m_gamma(batch_idx) = -m_main_diagonal(offset + 0);
            /* main_diagonal(0) -= gamma_;*/
            m_main_diagonal(offset + 0) -= m_gamma(batch_idx);
            /* main_diagonal(matrix_dimension_ - 1) -= cyclic_corner_element()^2 / gamma_;*/
            m_main_diagonal(offset + m_matrix_dimension - 1) -=
                cyclic_corner_element * cyclic_corner_element / m_gamma(batch_idx);

            // ---------------------- //
            // Cholesky Decomposition //
            for (int i = 1; i < m_matrix_dimension; i++) {
                m_sub_diagonal(offset + i - 1) /= m_main_diagonal(offset + i - 1);
                const T factor = m_sub_diagonal(offset + i - 1);
                m_main_diagonal(offset + i) -= factor * factor * m_main_diagonal(offset + i - 1);
            }
        }
    };

    void setup()
    {
        if (!is_cyclic_) {
            SetupNonCyclic functor{matrix_dimension_, main_diagonal_, sub_diagonal_};
            Kokkos::parallel_for("SetupNonCyclic", batch_count_, functor);
        }
        else {
            SetupCyclic functor{matrix_dimension_, main_diagonal_, sub_diagonal_, gamma_};
            Kokkos::parallel_for("SetupCyclic", batch_count_, functor);
        }
        Kokkos::fence();
        is_factorized_ = true;
    }

    /* ---------------------------------------- */
    /* Solve: Forward and Backward Substitution */
    /* ---------------------------------------- */
    // This step solves the system Ax = b using the factorized form of A.
    // For cyclic systems, it also performs the Shermann-Morrison reconstruction to obtain the final solution.

    struct SolveNonCyclic {
        int m_matrix_dimension;
        Vector<T> m_main_diagonal;
        Vector<T> m_sub_diagonal;
        Vector<T> m_rhs;
        int m_batch_offset;
        int m_batch_stride;

        void operator()(const int k) const
        {
            // ----------------------------------- //
            // Obtain offset for the current batch //
            int batch_idx = m_batch_stride * k + m_batch_offset;
            int offset    = batch_idx * m_matrix_dimension;

            // -------------------- //
            // Forward Substitution //
            for (int i = 1; i < m_matrix_dimension; i++) {
                m_rhs(offset + i) -= m_sub_diagonal(offset + i - 1) * m_rhs(offset + i - 1);
            }
            // ---------------- //
            // Diagonal Scaling //
            for (int i = 0; i < m_matrix_dimension; i++) {
                m_rhs(offset + i) /= m_main_diagonal(offset + i);
            }
            // --------------------- //
            // Backward Substitution //
            for (int i = m_matrix_dimension - 2; i >= 0; i--) {
                m_rhs(offset + i) -= m_sub_diagonal(offset + i) * m_rhs(offset + i + 1);
            }
        }
    };

    struct SolveCyclic {
        int m_matrix_dimension;
        Vector<T> m_main_diagonal;
        Vector<T> m_sub_diagonal;
        Vector<T> m_buffer;
        Vector<T> m_gamma;
        Vector<T> m_rhs;
        int m_batch_offset;
        int m_batch_stride;

        void operator()(const int k) const
        {
            // ----------------------------------- //
            // Obtain offset for the current batch //
            int batch_idx = m_batch_stride * k + m_batch_offset;
            int offset    = batch_idx * m_matrix_dimension;

            // -------------------- //
            // Forward Substitution //
            T cyclic_corner_element = m_sub_diagonal(offset + m_matrix_dimension - 1);
            m_buffer(offset + 0)    = m_gamma(batch_idx);
            for (int i = 1; i < m_matrix_dimension; i++) {
                m_rhs(offset + i) -= m_sub_diagonal(offset + i - 1) * m_rhs(offset + i - 1);
                if (i < m_matrix_dimension - 1)
                    m_buffer(offset + i) = 0.0 - m_sub_diagonal(offset + i - 1) * m_buffer(offset + i - 1);
                else
                    m_buffer(offset + i) =
                        cyclic_corner_element - m_sub_diagonal(offset + i - 1) * m_buffer(offset + i - 1);
            }
            // ---------------- //
            // Diagonal Scaling //
            for (int i = 0; i < m_matrix_dimension; i++) {
                m_rhs(offset + i) /= m_main_diagonal(offset + i);
                m_buffer(offset + i) /= m_main_diagonal(offset + i);
            }
            // --------------------- //
            // Backward Substitution //
            for (int i = m_matrix_dimension - 2; i >= 0; i--) {
                m_rhs(offset + i) -= m_sub_diagonal(offset + i) * m_rhs(offset + i + 1);
                m_buffer(offset + i) -= m_sub_diagonal(offset + i) * m_buffer(offset + i + 1);
            }
            // ------------------------------- //
            // Shermann-Morrison Reonstruction //
            const T dot_product_x_v =
                m_rhs(offset + 0) + cyclic_corner_element / m_gamma(batch_idx) * m_rhs(offset + m_matrix_dimension - 1);
            const T dot_product_u_v = m_buffer(offset + 0) + cyclic_corner_element / m_gamma(batch_idx) *
                                                                 m_buffer(offset + m_matrix_dimension - 1);
            const T factor = dot_product_x_v / (1.0 + dot_product_u_v);

            for (int i = 0; i < m_matrix_dimension; i++) {
                m_rhs(offset + i) -= factor * m_buffer(offset + i);
            }
        }
    };

    void solve(Vector<T> rhs, int batch_offset = 0, int batch_stride = 1)
    {
        if (!is_factorized_) {
            throw std::runtime_error("Error: Matrix must be factorized before solving.");
        }

        // Compute the effective number of batches to solve
        int effective_batch_count = (batch_count_ - batch_offset + batch_stride - 1) / batch_stride;

        if (!is_cyclic_) {
            SolveNonCyclic functor{matrix_dimension_, main_diagonal_, sub_diagonal_, rhs, batch_offset, batch_stride};
            Kokkos::parallel_for("SolveNonCyclic", effective_batch_count, functor);
        }
        else {
            SolveCyclic functor{matrix_dimension_, main_diagonal_, sub_diagonal_, buffer_, gamma_, rhs,
                                batch_offset,      batch_stride};
            Kokkos::parallel_for("SolveCyclic", effective_batch_count, functor);
        }
        Kokkos::fence();
    }

    /* ---------------------------- */
    /* Solve: Diagonal Scaling Only */
    /* ---------------------------- */
    // This step performs only the diagonal scaling part of the solve process.
    // It is useful when the matrix has a non-zero diagonal but zero off-diagonal entries.
    // Note that .setup() doesn't modify the main diagonal in this case.

    struct SolveDiagonalNonCyclic {
        int m_matrix_dimension;
        Vector<T> m_main_diagonal;
        Vector<T> m_rhs;
        int m_batch_offset;
        int m_batch_stride;

        void operator()(const int k) const
        {
            // ----------------------------------- //
            // Obtain offset for the current batch //
            int batch_idx = m_batch_stride * k + m_batch_offset;
            int offset    = batch_idx * m_matrix_dimension;

            // ---------------- //
            // Diagonal Scaling //
            for (int i = 0; i < m_matrix_dimension; i++) {
                m_rhs(offset + i) /= m_main_diagonal(offset + i);
            }
        }
    };

    struct SolveDiagonalCyclic {
        int m_matrix_dimension;
        Vector<T> m_main_diagonal;
        Vector<T> m_gamma;
        Vector<T> m_rhs;
        int m_batch_offset;
        int m_batch_stride;

        void operator()(const int k) const
        {
            // ----------------------------------- //
            // Obtain offset for the current batch //
            int batch_idx = m_batch_stride * k + m_batch_offset;
            int offset    = batch_idx * m_matrix_dimension;

            // ---------------- //
            // Diagonal Scaling //
            m_rhs(offset + 0) /= m_main_diagonal(offset + 0) + m_gamma(batch_idx);
            for (int i = 1; i < m_matrix_dimension; i++) {
                m_rhs(offset + i) /= m_main_diagonal(offset + i);
            }
        }
    };

    void solve_diagonal(Vector<T> rhs, int batch_offset = 0, int batch_stride = 1)
    {
        if (!is_factorized_) {
            throw std::runtime_error("Error: Matrix must be factorized before solving.");
        }

        // Compute the effective number of batches to solve
        int effective_batch_count = (batch_count_ - batch_offset + batch_stride - 1) / batch_stride;

        if (!is_cyclic_) {
            SolveDiagonalNonCyclic functor{matrix_dimension_, main_diagonal_, rhs, batch_offset, batch_stride};
            Kokkos::parallel_for("SolveDiagonalNonCyclic", effective_batch_count, functor);
        }
        else {
            SolveDiagonalCyclic functor{matrix_dimension_, main_diagonal_, gamma_, rhs, batch_offset, batch_stride};
            Kokkos::parallel_for("SolveDiagonalCyclic", effective_batch_count, functor);
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
