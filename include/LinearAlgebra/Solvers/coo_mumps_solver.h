#pragma once

#ifdef GMGPOLAR_USE_MUMPS

    #include "dmumps_c.h"
    #include "../../LinearAlgebra/Matrix/coo_matrix.h"
    #include "../../LinearAlgebra/Vector/vector.h"

/*
 * Wraps MUMPS for solving sparse linear systems given in COO format.
 * For general matrices, all non-zero entries must be provided.
 * For symmetric matrices (is_symmetric = true), only the upper or lower
 * triangular entries should be provided, and the matrix is assumed to be
 * positive definite. In GMGPolar this holds true since the domain mapping is invertible.
 */
class CooMumpsSolver
{
public:
    explicit CooMumpsSolver(SparseMatrixCOO<double> matrix);
    ~CooMumpsSolver();

    // rhs is overwritten in-place with the solution on return.
    void solve(Vector<double>& rhs);

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

    SparseMatrixCOO<double> matrix_;
    DMUMPS_STRUC_C mumps_solver_ = {};
};

#endif // GMGPOLAR_USE_MUMPS
