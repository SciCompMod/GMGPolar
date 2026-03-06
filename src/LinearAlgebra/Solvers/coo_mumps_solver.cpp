#include "../../../include/LinearAlgebra/Solvers/coo_mumps_solver.h"

#ifdef GMGPOLAR_USE_MUMPS

    #include <iostream>
    #include <cassert>
    #include <stdexcept>

CooMumpsSolver::CooMumpsSolver(SparseMatrixCOO<double> matrix)
    : matrix_(std::move(matrix))
{
    initialize();
}

CooMumpsSolver::~CooMumpsSolver()
{
    finalize();
}

void CooMumpsSolver::solve(Vector<double>& rhs)
{
    assert(std::ssize(rhs) == mumps_solver_.n);

    mumps_solver_.job  = JOB_COMPUTE_SOLUTION;
    mumps_solver_.nrhs = 1;
    mumps_solver_.lrhs = mumps_solver_.n; // leading dimension: must equal n for dense centralized RHS
    mumps_solver_.rhs  = rhs.data(); // in: RHS, out: solution (overwritten in-place)
    dmumps_c(&mumps_solver_);

    if (INFOG(1) != 0) {
        std::cerr << "MUMPS reported an error during solution phase " << "(INFOG(1) = " << INFOG(1) << ").\n";
    }
}

void CooMumpsSolver::initialize()
{
    assert(matrix_.rows() == matrix_.columns());

    /*
     * MUMPS uses 1-based indexing; our COO matrix uses 0-based indexing.
     * Adjust row and column indices to match MUMPS' requirements.
     */
    for (int i = 0; i < matrix_.non_zero_size(); i++) {
        matrix_.row_index(i) += 1;
        matrix_.col_index(i) += 1;
    }

    mumps_solver_.job          = JOB_INIT;
    mumps_solver_.par          = PAR_PARALLEL;
    mumps_solver_.sym          = matrix_.is_symmetric() ? SYM_POSITIVE_DEFINITE : SYM_UNSYMMETRIC;
    mumps_solver_.comm_fortran = USE_COMM_WORLD;
    dmumps_c(&mumps_solver_);

    configureICNTL();
    configureCNTL();

    mumps_solver_.job = JOB_ANALYSIS_AND_FACTORIZATION;
    mumps_solver_.n   = matrix_.rows();
    mumps_solver_.nz  = matrix_.non_zero_size();
    mumps_solver_.irn = matrix_.row_indices_data();
    mumps_solver_.jcn = matrix_.column_indices_data();
    mumps_solver_.a   = matrix_.values_data();
    dmumps_c(&mumps_solver_);

    if (INFOG(1) != 0) {
        std::cerr << "MUMPS reported an error during analysis/factorization " << "(INFOG(1) = " << INFOG(1) << ").\n";
        return;
    }

    if (mumps_solver_.sym == SYM_POSITIVE_DEFINITE && INFOG(12) != 0) {
        std::cerr << "Matrix declared positive definite, "
                  << "but negative pivots were encountered during factorization " << "(INFOG(12) = " << INFOG(12)
                  << ").\n";
    }
}

void CooMumpsSolver::finalize()
{
    mumps_solver_.job = JOB_END;
    dmumps_c(&mumps_solver_);
}

void CooMumpsSolver::configureICNTL()
{
    // All ICNTL values are left at their defaults unless noted below.
    // ICNTL(1) = 0: suppress error message output
    // ICNTL(3) = 0: suppress global information output
    // ICNTL(6) = 7: automatically choose permutation/scaling strategy
    // ICNTL(7) = 5: use METIS for fill-reducing ordering
    // ICNTL(48) = 0: disable tree parallelism (conflicts with OpenMP in newer MUMPS releases)

    ICNTL(1)  = 0; // Output stream for error messages
    ICNTL(2)  = 0; // Output stream for diagnostic printing and statistics local to each MPI process
    ICNTL(3)  = 0; // Output stream for global information, collected on the host
    ICNTL(4)  = 0; // Level of printing for error, warning, and diagnostic messages
    ICNTL(5)  = 0; // Controls the matrix input format
    ICNTL(6)  = 7; // Permutes the matrix to a zero-free diagonal and/or scales the matrix
    ICNTL(7)  = 5; // Symmetric permutation (ordering) to determine pivot order for sequential analysis
    ICNTL(8)  = 77; // Scaling strategy
    ICNTL(9)  = 1; // Computes the solution using A or A^T
    ICNTL(10) = 0; // Iterative refinement steps applied to the computed solution
    ICNTL(11) = 0; // Error analysis statistics
    ICNTL(12) = 0; // Ordering strategy for symmetric matrices
    ICNTL(13) = 0; // Controls the parallelism of the root node
    ICNTL(14) = matrix_.is_symmetric() ? 5 : 20; // Percentage increase in estimated working space
    ICNTL(15) = 0; // Exploits compression of the input matrix resulting from a block format
    ICNTL(16) = 0; // Controls the setting of the number of OpenMP threads
    // ICNTL(17) does not exist
    ICNTL(18) = 0; // Strategy for the distributed input matrix
    ICNTL(19) = 0; // Computes the Schur complement matrix
    ICNTL(20) = 0; // Format of the right-hand sides (dense, sparse, or distributed)
    ICNTL(21) = 0; // Distribution of the solution vectors (centralized or distributed)
    ICNTL(22) = 0; // In-core/out-of-core (OOC) factorization and solve
    ICNTL(23) = 0; // Maximum working memory in MegaBytes per working process
    ICNTL(24) = 0; // Detection of null pivot rows
    ICNTL(25) = 0; // Solution of a deficient matrix and null space basis computation
    ICNTL(26) = 0; // Solution phase when Schur complement has been computed
    ICNTL(27) = -32; // Blocking size for multiple right-hand sides
    ICNTL(28) = 0; // Sequential or parallel computation of the ordering
    ICNTL(29) = 0; // Parallel ordering tool when ICNTL(28)=1
    ICNTL(30) = 0; // User-specified entries in the inverse A^-1
    ICNTL(31) = 0; // Which factors may be discarded during factorization
    ICNTL(32) = 0; // Forward elimination of the right-hand sides during factorization
    ICNTL(33) = 0; // Computes the determinant of the input matrix
    ICNTL(34) = 0; // Conservation of OOC files during JOB=-3
    ICNTL(35) = 0; // Activation of the BLR feature
    ICNTL(36) = 0; // Choice of BLR factorization variant
    ICNTL(37) = 0; // BLR compression of the contribution blocks
    ICNTL(38) = 600; // Estimated compression rate of LU factors
    ICNTL(39) = 500; // Estimated compression rate of contribution blocks
    // ICNTL(40-47) do not exist
    ICNTL(48) = 0; // Multithreading with tree parallelism
    ICNTL(49) = 0; // Compact workarray id%S at end of factorization phase
    // ICNTL(50-55) do not exist
    ICNTL(56) = 0; // Detects pseudo-singularities; rank-revealing factorization of root node
    // ICNTL(57) does not exist
    ICNTL(58) = 2; // Options for symbolic factorization
    // ICNTL(59-60) do not exist
}

void CooMumpsSolver::configureCNTL()
{
    // All CNTL values are left at their defaults unless noted below.

    CNTL(1) = -1.0; // Relative threshold for numerical pivoting
    CNTL(2) = -1.0; // Stopping criterion for iterative refinement
    CNTL(3) = 0.0; // Threshold for null pivot row detection
    CNTL(4) = -1.0; // Threshold for static pivoting
    CNTL(5) = 0.0; // Fixation for null pivots (effective only when null pivot detection is active)
    // CNTL(6) does not exist
    CNTL(7) = 0.0; // Dropping parameter precision for BLR compression
    // CNTL(8-15) do not exist
}

#endif // GMGPOLAR_USE_MUMPS
