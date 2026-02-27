#include "../../../include/Smoother/SmootherGive/smootherGive.h"

#ifdef GMGPOLAR_USE_MUMPS

void SmootherGive::initializeMumpsSolver(DMUMPS_STRUC_C& mumps_solver, SparseMatrixCOO<double>& solver_matrix)
{
    /* 
     * MUMPS (a parallel direct solver) uses 1-based indexing, 
     * whereas the input matrix follows 0-based indexing. 
     * Adjust row and column indices to match MUMPS' requirements.
     */
    for (int i = 0; i < solver_matrix.non_zero_size(); i++) {
        solver_matrix.row_index(i) += 1;
        solver_matrix.col_index(i) += 1;
    }

    mumps_solver.job = JOB_INIT;
    mumps_solver.par = PAR_PARALLEL;
    /* The matrix is positive definite for invertible mappings. */
    /* Therefore we use SYM_POSITIVE_DEFINITE instead of SYM_GENERAL_SYMMETRIC. */
    mumps_solver.sym          = (solver_matrix.is_symmetric() ? SYM_POSITIVE_DEFINITE : SYM_UNSYMMETRIC);
    mumps_solver.comm_fortran = USE_COMM_WORLD;
    dmumps_c(&mumps_solver);

    ICNTL(mumps_solver, 1) = 0; // Output stream for error messages.
    ICNTL(mumps_solver, 2) = 0; // Output stream for diagnostic printing and statistics local to each MPI process.
    ICNTL(mumps_solver, 3) = 0; // Output stream for global information, collected on the host
    ICNTL(mumps_solver, 4) = 0; // Level of printing for error, warning, and diagnostic messages.
    ICNTL(mumps_solver, 5) = 0; // Controls the matrix input format
    ICNTL(mumps_solver, 6) = 7; // Permutes the matrix to a zero-free diagonal and/or scale the matrix
    ICNTL(mumps_solver, 7) =
        5; // Computes a symmetric permutation (ordering) to determine the pivot order to be used for the factorization in case of sequential analysis
    ICNTL(mumps_solver, 8)  = 77; // Describes the scaling strategy
    ICNTL(mumps_solver, 9)  = 1; // Computes the solution using A or A^T
    ICNTL(mumps_solver, 10) = 0; // Applies the iterative refinement to the computed solution
    ICNTL(mumps_solver, 11) = 0; // Computes statistics related to an error analysis of the linear system solved
    ICNTL(mumps_solver, 12) = 0; // Defines an ordering strategy for symmetric matrices and is used
    ICNTL(mumps_solver, 13) = 0; // Controls the parallelism of the root node
    ICNTL(mumps_solver, 14) = // Controls the percentage increase in the estimated working space
        (solver_matrix.is_symmetric() ? 5 : 20);
    ICNTL(mumps_solver, 15) = 0; // Exploits compression of the input matrix resulting from a block format
    ICNTL(mumps_solver, 16) = 0; // Controls the setting of the number of OpenMP threads
    // ICNTL(17) Doesn't exist
    ICNTL(mumps_solver, 18) = 0; // Defines the strategy for the distributed input matrix
    ICNTL(mumps_solver, 19) = 0; // Computes the Schur complement matrix
    ICNTL(mumps_solver, 20) = 0; // Determines the format (dense, sparse, or distributed) of the right-hand sides
    ICNTL(mumps_solver, 21) = 0; // Determines the distribution (centralized or distributed) of the solution vectors.
    ICNTL(mumps_solver, 22) = 0; // Controls the in-core/out-of-core (OOC) factorization and solve.
    ICNTL(mumps_solver, 23) = 0; // Corresponds to the maximum size of the working memory in MegaBytes that MUMPS can
    //                             allocate per working process
    ICNTL(mumps_solver, 24) = 0; // Controls the detection of “null pivot rows”.
    ICNTL(mumps_solver, 25) =
        0; // Allows the computation of a solution of a deficient matrix and also of a null space basis
    ICNTL(mumps_solver, 26) = 0; // Drives the solution phase if a Schur complement matrix has been computed
    ICNTL(mumps_solver, 27) = -32; // Controls the blocking size for multiple right-hand sides.
    ICNTL(mumps_solver, 28) = 0; // Determines whether a sequential or parallel computation of the ordering is performed
    ICNTL(mumps_solver, 29) =
        0; // Defines the parallel ordering tool (when ICNTL(28)=1) to be used to compute the fill-in reducing permutation.
    ICNTL(mumps_solver, 30) = 0; // Computes a user-specified set of entries in the inverse A^−1 of the original matrix
    ICNTL(mumps_solver, 31) = 0; // Indicates which factors may be discarded during the factorization.
    ICNTL(mumps_solver, 32) = 0; // Performs the forward elimination of the right-hand sides during the factorization
    ICNTL(mumps_solver, 33) = 0; // Computes the determinant of the input matrix.
    ICNTL(mumps_solver, 34) = 0; // Controls the conservation of the OOC files during JOB= –3
    ICNTL(mumps_solver, 35) = 0; // Controls the activation of the BLR feature
    ICNTL(mumps_solver, 36) = 0; // Controls the choice of BLR factorization variant
    ICNTL(mumps_solver, 37) = 0; // Controls the BLR compression of the contribution blocks
    ICNTL(mumps_solver, 38) = 600; // Estimates compression rate of LU factors
    ICNTL(mumps_solver, 39) = 500; // Estimates compression rate of contribution blocks
    // ICNTL(40-47) Don't exist
    ICNTL(mumps_solver, 48) = 0; // Multithreading with tree parallelism
    ICNTL(mumps_solver, 49) = 0; // Compact workarray id%S at the end of factorization phase
    // ICNTL(50-55) Don't exist
    ICNTL(mumps_solver, 56) =
        0; // Detects pseudo-singularities during factorization and factorizes the root node with a rankrevealing method
    // ICNTL(57) Doesn't exist
    ICNTL(mumps_solver, 58) = 2; // Defines options for symbolic factorization
    // ICNTL(59-60) Don't exist

    CNTL(mumps_solver, 1) = -1.0; // Relative threshold for numerical pivoting
    CNTL(mumps_solver, 2) = -1.0; // Stopping criterion for iterative refinement
    CNTL(mumps_solver, 3) = 0.0; // Determine null pivot rows
    CNTL(mumps_solver, 4) = -1.0; // Determines the threshold for static pivoting
    CNTL(mumps_solver, 5) =
        0.0; // Defines the fixation for null pivots and is effective only when null pivot row detection is active
    // CNTL(6) Doesn't exist
    CNTL(mumps_solver, 7) = 0.0; // Defines the precision of the dropping parameter used during BLR compression
    // CNTL(8-15) Don't exist

    mumps_solver.job = JOB_ANALYSIS_AND_FACTORIZATION;
    assert(solver_matrix.rows() == solver_matrix.columns());
    mumps_solver.n   = solver_matrix.rows();
    mumps_solver.nz  = solver_matrix.non_zero_size();
    mumps_solver.irn = solver_matrix.row_indices_data();
    mumps_solver.jcn = solver_matrix.column_indices_data();
    mumps_solver.a   = solver_matrix.values_data();
    dmumps_c(&mumps_solver);

    if (mumps_solver.sym == SYM_POSITIVE_DEFINITE && INFOG(mumps_solver, 12) != 0) {
        std::cout << "Warning: Smoother inner boundary matrix is not positive definite: Negative pivots in the "
                     "factorization phase."
                  << std::endl;
    }
}

void SmootherGive::finalizeMumpsSolver(DMUMPS_STRUC_C& mumps_solver)
{
    mumps_solver.job = JOB_END;
    dmumps_c(&mumps_solver);
}

#endif
