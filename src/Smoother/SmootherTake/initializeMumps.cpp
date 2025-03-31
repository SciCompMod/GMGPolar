#include "../../../include/Smoother/SmootherTake/smootherTake.h"

#ifdef GMGPOLAR_USE_MUMPS

void SmootherTake::initializeMumpsSolver(DMUMPS_STRUC_C& mumps_solver, SparseMatrixCOO<double>& solver_matrix)
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

    mumps_solver.ICNTL(1) = 0; // Output stream for error messages.
    mumps_solver.ICNTL(2) = 0; // Output stream for diagnostic printing and statistics local to each MPI process.
    mumps_solver.ICNTL(3) = 0; // Output stream for global information, collected on the host
    mumps_solver.ICNTL(4) = 0; // Level of printing for error, warning, and diagnostic messages.
    mumps_solver.ICNTL(5) = 0; // Controls the matrix input format
    mumps_solver.ICNTL(6) = 7; // Permutes the matrix to a zero-free diagonal and/or scale the matrix
    mumps_solver.ICNTL(7) =
        5; // Computes a symmetric permutation (ordering) to determine the pivot order to be used for the factorization in case of sequential analysis
    mumps_solver.ICNTL(8)  = 77; // Describes the scaling strategy
    mumps_solver.ICNTL(9)  = 1; // Computes the solution using A or A^T
    mumps_solver.ICNTL(10) = 0; // Applies the iterative refinement to the computed solution
    mumps_solver.ICNTL(11) = 0; // Computes statistics related to an error analysis of the linear system solved
    mumps_solver.ICNTL(12) = 0; // Defines an ordering strategy for symmetric matrices and is used
    mumps_solver.ICNTL(13) = 0; // Controls the parallelism of the root node
    mumps_solver.ICNTL(14) = // Controls the percentage increase in the estimated working space
        (solver_matrix.is_symmetric() ? 5 : 20);
    mumps_solver.ICNTL(15) = 0; // Exploits compression of the input matrix resulting from a block format
    mumps_solver.ICNTL(16) = 0; // Controls the setting of the number of OpenMP threads
    // ICNTL(17) Doesn't exist
    mumps_solver.ICNTL(18) = 0; // Defines the strategy for the distributed input matrix
    mumps_solver.ICNTL(19) = 0; // Computes the Schur complement matrix
    mumps_solver.ICNTL(20) = 0; // Determines the format (dense, sparse, or distributed) of the right-hand sides
    mumps_solver.ICNTL(21) = 0; // Determines the distribution (centralized or distributed) of the solution vectors.
    mumps_solver.ICNTL(22) = 0; // Controls the in-core/out-of-core (OOC) factorization and solve.
    mumps_solver.ICNTL(23) = 0; // Corresponds to the maximum size of the working memory in MegaBytes that MUMPS can
    //                             allocate per working process
    mumps_solver.ICNTL(24) = 0; // Controls the detection of “null pivot rows”.
    mumps_solver.ICNTL(25) =
        0; // Allows the computation of a solution of a deficient matrix and also of a null space basis
    mumps_solver.ICNTL(26) = 0; // Drives the solution phase if a Schur complement matrix has been computed
    mumps_solver.ICNTL(27) = -32; // Controls the blocking size for multiple right-hand sides.
    mumps_solver.ICNTL(28) = 0; // Determines whether a sequential or parallel computation of the ordering is performed
    mumps_solver.ICNTL(29) =
        0; // Defines the parallel ordering tool (when ICNTL(28)=1) to be used to compute the fill-in reducing permutation.
    mumps_solver.ICNTL(30) = 0; // Computes a user-specified set of entries in the inverse A^−1 of the original matrix
    mumps_solver.ICNTL(31) = 0; // Indicates which factors may be discarded during the factorization.
    mumps_solver.ICNTL(32) = 0; // Performs the forward elimination of the right-hand sides during the factorization
    mumps_solver.ICNTL(33) = 0; // Computes the determinant of the input matrix.
    mumps_solver.ICNTL(34) = 0; // Controls the conservation of the OOC files during JOB= –3
    mumps_solver.ICNTL(35) = 0; // Controls the activation of the BLR feature
    mumps_solver.ICNTL(36) = 0; // Controls the choice of BLR factorization variant
    mumps_solver.ICNTL(37) = 0; // Controls the BLR compression of the contribution blocks
    mumps_solver.ICNTL(38) = 600; // Estimates compression rate of LU factors
    mumps_solver.ICNTL(39) = 500; // Estimates compression rate of contribution blocks
    // ICNTL(40-47) Don't exist
    mumps_solver.ICNTL(48) = 1; // Multithreading with tree parallelism
    mumps_solver.ICNTL(49) = 0; // Compact workarray id%S at the end of factorization phase
    // ICNTL(50-55) Don't exist
    mumps_solver.ICNTL(56) =
        0; // Detects pseudo-singularities during factorization and factorizes the root node with a rankrevealing method
    // ICNTL(57) Doesn't exist
    mumps_solver.ICNTL(58) = 2; // Defines options for symbolic factorization
    // ICNTL(59-60) Don't exist

    mumps_solver.CNTL(1) = -1.0; // Relative threshold for numerical pivoting
    mumps_solver.CNTL(2) = -1.0; // Stopping criterion for iterative refinement
    mumps_solver.CNTL(3) = 0.0; // Determine null pivot rows
    mumps_solver.CNTL(4) = -1.0; // Determines the threshold for static pivoting
    mumps_solver.CNTL(5) =
        0.0; // Defines the fixation for null pivots and is effective only when null pivot row detection is active
    // CNTL(6) Doesn't exist
    mumps_solver.CNTL(7) = 0.0; // Defines the precision of the dropping parameter used during BLR compression
    // CNTL(8-15) Don't exist

    mumps_solver.job = JOB_ANALYSIS_AND_FACTORIZATION;
    assert(solver_matrix.rows() == solver_matrix.columns());
    mumps_solver.n   = solver_matrix.rows();
    mumps_solver.nz  = solver_matrix.non_zero_size();
    mumps_solver.irn = solver_matrix.row_indices_data();
    mumps_solver.jcn = solver_matrix.column_indices_data();
    mumps_solver.a   = solver_matrix.values_data();
    dmumps_c(&mumps_solver);

    if (mumps_solver.sym == SYM_POSITIVE_DEFINITE && mumps_solver.INFOG(12) != 0) {
        std::cout << "Warning: Smoother inner boundary matrix is not positive definite: Negative pivots in the "
                     "factorization phase."
                  << std::endl;
    }
}

void SmootherTake::finalizeMumpsSolver(DMUMPS_STRUC_C& mumps_solver)
{
    mumps_solver.job = JOB_END;
    dmumps_c(&mumps_solver);
}

#endif