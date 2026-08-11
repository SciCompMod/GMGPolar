#pragma once

#ifdef GMGPOLAR_USE_MUMPS

    #include "dmumps_c.h"
    #include "../../LinearAlgebra/Matrix/coo_matrix.h"
    #include "../../LinearAlgebra/Vector/vector.h"

/*
 * CooMumpsSolver - sparse direct solver backed by MUMPS (MUltifrontal Massively Parallel Solver).
 *
 * Accepts a matrix in COO (Coordinate) format and solves systems of the form Ax = b.
 * The solver takes ownership of the matrix, as MUMPS may reorder and modify it in place
 * during the analysis and factorization phases.
 *
 * Usage:
 *   1. Construct with a fully assembled COO matrix (all non-zeros must be present,
 *      including both triangles if the matrix is symmetric).
 *   2. Call solve(rhs) one or more times. The right-hand side vector is overwritten
 *      in place with the solution x on return.
 *
 * Symmetry:
 *   If the matrix is marked symmetric (is_symmetric = true), it is assumed to be
 *   symmetric positive definite and MUMPS is configured accordingly (SYM = 1).
 *   In that case, only the upper triangle is passed to MUMPS internally.
 *   This assumption is valid within GMGPolar because all system matrices arise from
 *   an invertible domain mapping, guaranteeing positive definiteness.
 *
 * Thread safety:
 *   Not thread-safe. Each instance manages its own MUMPS communicator and
 *   must not be shared across threads.
 */
namespace gmgpolar
{

class CooMumpsSolver
{
public:
    template <class MemorySpace>
    explicit CooMumpsSolver(SparseMatrixCOO<double, MemorySpace>&& matrix)
    {
        auto matrix_host = matrix.template mirror_view_and_copy<Kokkos::HostSpace>();
        if (matrix.is_symmetric()) {
            matrix_ = extractUpperTriangle(matrix_host);
        }
        else {
            matrix_ = std::move(matrix_host);
        }

        initialize();
    }

    ~CooMumpsSolver();

    // rhs is overwritten in-place with the solution on return.
    void solveInPlace(Vector<double>& rhs);

    /**
     * @brief Refactorize using a matrix with the same sparsity pattern as the
     * one originally used to construct this solver.
     *
     * Reuses the existing MUMPS analysis/ordering phase (job 1) and only
     * reruns the numeric factorization (job 2). The matrix must have the same
     * non-zero structure - same number of entries in the same order - as the
     * matrix originally passed to the constructor.
     *
     * Takes matrix by const reference (unlike the constructor, which takes
     * ownership) so the caller can keep refilling and reusing their own
     * matrix storage across repeated updates.
     */
    template <class MemorySpace>
    void updateValues(const SparseMatrixCOO<double, MemorySpace>& matrix)
    {
        auto matrix_host = matrix.template mirror_view_and_copy<Kokkos::HostSpace>();

        if (matrix_host.is_symmetric()) {
            matrix_ = extractUpperTriangle(matrix_host);
        }
        else {
            matrix_ = std::move(matrix_host);
        }

        assert(matrix_.non_zero_size() == mumps_solver_.nz &&
               "Matrix structure must match the originally factorized matrix");

        // MUMPS uses 1-based indexing.
        for (int i = 0; i < matrix_.non_zero_size(); i++) {
            matrix_.increment_row_index(i);
            matrix_.increment_col_index(i);
        }

        // irn/jcn must be repointed even though structure is unchanged, since
        // extractUpperTriangle (if symmetric) allocates a new matrix_ object
        // with a new underlying buffer each call.
        mumps_solver_.job = JOB_FACTORIZATION_PHASE; // reuse stored analysis/ordering
        mumps_solver_.irn = matrix_.row_indices_data();
        mumps_solver_.jcn = matrix_.column_indices_data();
        mumps_solver_.a   = matrix_.values_data();

        dmumps_c(&mumps_solver_);

        if (INFOG(1) != 0) {
            std::cerr << "MUMPS reported an error during factorization update "
                      << "(INFOG(1) = " << INFOG(1) << ").\n";
        }

        if (mumps_solver_.sym == SYM_POSITIVE_DEFINITE && INFOG(12) != 0) {
            std::cerr << "Matrix declared positive definite, "
                      << "but negative pivots were encountered during refactorization "
                      << "(INFOG(12) = " << INFOG(12) << ").\n";
        }
    }

private:
    void initialize();
    void finalize();
    void configureICNTL();
    void configureCNTL();

    /* ------------------------------------------------ */
    /* MUMPS uses 1-based indexing in the documentation */
    /* ------------------------------------------------ */
    int& ICNTL(int i)
    {
        return mumps_solver_.icntl[i - 1];
    }
    double& CNTL(int i)
    {
        return mumps_solver_.cntl[i - 1];
    }
    int& INFOG(int i)
    {
        return mumps_solver_.infog[i - 1];
    }

    SparseMatrixCOO<double, Kokkos::HostSpace>
    extractUpperTriangle(const SparseMatrixCOO<double, Kokkos::HostSpace>& matrix) const;

    /* ----------------------------------- */
    /* MUMPS jobs and constant definitions */
    /* ----------------------------------- */
    static constexpr int USE_COMM_WORLD   = -987654;
    static constexpr int PAR_NOT_PARALLEL = 0;
    static constexpr int PAR_PARALLEL     = 1;

    static constexpr int JOB_INIT               = -1;
    static constexpr int JOB_END                = -2;
    static constexpr int JOB_REMOVE_SAVED_DATA  = -3;
    static constexpr int JOB_FREE_INTERNAL_DATA = -4;
    static constexpr int JOB_SUPPRESS_OOC_FILES = -200;

    static constexpr int JOB_ANALYSIS_PHASE                  = 1;
    static constexpr int JOB_FACTORIZATION_PHASE             = 2;
    static constexpr int JOB_COMPUTE_SOLUTION                = 3;
    static constexpr int JOB_ANALYSIS_AND_FACTORIZATION      = 4;
    static constexpr int JOB_FACTORIZATION_AND_SOLUTION      = 5;
    static constexpr int JOB_ANALYSIS_FACTORIZATION_SOLUTION = 6;
    static constexpr int JOB_SAVE_INTERNAL_DATA              = 7;
    static constexpr int JOB_RESTORE_INTERNAL_DATA           = 8;
    static constexpr int JOB_DISTRIBUTE_RHS                  = 9;

    static constexpr int SYM_UNSYMMETRIC       = 0;
    static constexpr int SYM_POSITIVE_DEFINITE = 1;
    static constexpr int SYM_GENERAL_SYMMETRIC = 2;

    SparseMatrixCOO<double, Kokkos::HostSpace> matrix_;
    DMUMPS_STRUC_C mumps_solver_ = {};
};

} // namespace gmgpolar
#endif // GMGPOLAR_USE_MUMPS
