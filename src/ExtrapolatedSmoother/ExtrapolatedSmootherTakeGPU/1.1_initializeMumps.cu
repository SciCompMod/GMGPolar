#include "../../../include/ExtrapolatedSmoother/ExtrapolatedSmootherTakeGPU/extrapolatedSmoother.h"

void ExtrapolatedSmootherTakeGPU::initializeMumps(){
    inner_boundary_mumps_solver_.job = JOB_INIT;
    inner_boundary_mumps_solver_.par = PAR_PARALLEL;
    /* The matrix is positive definite for invertible mappings. */
    /* Therefore we use SYM_POSITIVE_DEFINITE instead of SYM_GENERAL_SYMMETRIC. */
    inner_boundary_mumps_solver_.sym          = SYM_UNSYMMETRIC;
    inner_boundary_mumps_solver_.comm_fortran = USE_COMM_WORLD;
    dmumps_c(&inner_boundary_mumps_solver_);

    inner_boundary_mumps_solver_.ICNTL(1) = 0; // Output stream for error messages.
    inner_boundary_mumps_solver_.ICNTL(2) = 0; // Output stream for diagnostic printing and statistics local to each MPI process.
    inner_boundary_mumps_solver_.ICNTL(3) = 0; // Output stream for global information, collected on the host
    inner_boundary_mumps_solver_.ICNTL(4) = 0; // Level of printing for error, warning, and diagnostic messages.
    inner_boundary_mumps_solver_.ICNTL(5) = 0; // Controls the matrix input format
    inner_boundary_mumps_solver_.ICNTL(6) = 7; // Permutes the matrix to a zero-free diagonal and/or scale the matrix
    inner_boundary_mumps_solver_.ICNTL(7) =
        5; // Computes a symmetric permutation (ordering) to determine the pivot order to be used for the
    //                             factorization in case of sequential analysis
    inner_boundary_mumps_solver_.ICNTL(8)  = 77; // Describes the scaling strategy
    inner_boundary_mumps_solver_.ICNTL(9)  = 1; // Computes the solution using A or A^T
    inner_boundary_mumps_solver_.ICNTL(10) = 0; // Applies the iterative refinement to the computed solution
    inner_boundary_mumps_solver_.ICNTL(11) = 0; // Computes statistics related to an error analysis of the linear system solved
    inner_boundary_mumps_solver_.ICNTL(12) = 0; // Defines an ordering strategy for symmetric matrices and is used
    inner_boundary_mumps_solver_.ICNTL(13) = 0; // Controls the parallelism of the root node
    inner_boundary_mumps_solver_.ICNTL(14) = 20;// Controls the percentage increase in the estimated working space
    inner_boundary_mumps_solver_.ICNTL(15) = 0; // Exploits compression of the input matrix resulting from a block format
    inner_boundary_mumps_solver_.ICNTL(16) = 0; // Controls the setting of the number of OpenMP threads
    // ICNTL(17) Doesn't exist
    inner_boundary_mumps_solver_.ICNTL(18) = 0; // Defines the strategy for the distributed input matrix
    inner_boundary_mumps_solver_.ICNTL(19) = 0; // Computes the Schur complement matrix
    inner_boundary_mumps_solver_.ICNTL(20) = 0; // Determines the format (dense, sparse, or distributed) of the right-hand sides
    inner_boundary_mumps_solver_.ICNTL(21) = 0; // Determines the distribution (centralized or distributed) of the solution vectors.
    inner_boundary_mumps_solver_.ICNTL(22) = 0; // Controls the in-core/out-of-core (OOC) factorization and solve.
    inner_boundary_mumps_solver_.ICNTL(23) = 0; // Corresponds to the maximum size of the working memory in MegaBytes that MUMPS can
    //                             allocate per working process
    inner_boundary_mumps_solver_.ICNTL(24) = 0; // Controls the detection of “null pivot rows”.
    inner_boundary_mumps_solver_.ICNTL(25) =
        0; // Allows the computation of a solution of a deficient matrix and also of a null space basis
    inner_boundary_mumps_solver_.ICNTL(26) = 0; // Drives the solution phase if a Schur complement matrix has been computed
    inner_boundary_mumps_solver_.ICNTL(27) = -32; // Controls the blocking size for multiple right-hand sides.
    inner_boundary_mumps_solver_.ICNTL(28) = 0; // Determines whether a sequential or parallel computation of the ordering is performed
    inner_boundary_mumps_solver_.ICNTL(29) =
        0; // Defines the parallel ordering tool (when ICNTL(28)=1) to be used to compute the fill-in reducing permutation.
    inner_boundary_mumps_solver_.ICNTL(30) = 0; // Computes a user-specified set of entries in the inverse A^−1 of the original matrix
    inner_boundary_mumps_solver_.ICNTL(31) = 0; // Indicates which factors may be discarded during the factorization.
    inner_boundary_mumps_solver_.ICNTL(32) = 0; // Performs the forward elimination of the right-hand sides during the factorization
    inner_boundary_mumps_solver_.ICNTL(33) = 0; // Computes the determinant of the input matrix.
    inner_boundary_mumps_solver_.ICNTL(34) = 0; // Controls the conservation of the OOC files during JOB= –3
    inner_boundary_mumps_solver_.ICNTL(35) = 0; // Controls the activation of the BLR feature
    inner_boundary_mumps_solver_.ICNTL(36) = 0; // Controls the choice of BLR factorization variant
    inner_boundary_mumps_solver_.ICNTL(37) = 0; // Controls the BLR compression of the contribution blocks
    inner_boundary_mumps_solver_.ICNTL(38) = 600; // Estimates compression rate of LU factors
    inner_boundary_mumps_solver_.ICNTL(39) = 500; // Estimates compression rate of contribution blocks
    // ICNTL(40-47) Don't exist
    inner_boundary_mumps_solver_.ICNTL(48) = 1; // Multithreading with tree parallelism
    inner_boundary_mumps_solver_.ICNTL(49) = 0; // Compact workarray id%S at the end of factorization phase
    // ICNTL(50-55) Don't exist
    inner_boundary_mumps_solver_.ICNTL(56) =
        0; // Detects pseudo-singularities during factorization and factorizes the root node with a rankrevealing method
    // ICNTL(57) Doesn't exist
    inner_boundary_mumps_solver_.ICNTL(58) = 2; // Defines options for symbolic factorization
    // ICNTL(59-60) Don't exist

    inner_boundary_mumps_solver_.CNTL(1) = -1.0; // Relative threshold for numerical pivoting
    inner_boundary_mumps_solver_.CNTL(2) = -1.0; // Stopping criterion for iterative refinement
    inner_boundary_mumps_solver_.CNTL(3) = 0.0; // Determine null pivot rows
    inner_boundary_mumps_solver_.CNTL(4) = -1.0; // Determines the threshold for static pivoting
    inner_boundary_mumps_solver_.CNTL(5) =
        0.0; // Defines the fixation for null pivots and is effective only when null pivot row detection is active
    // CNTL(6) Doesn't exist
    inner_boundary_mumps_solver_.CNTL(7) = 0.0; // Defines the precision of the dropping parameter used during BLR compression
    // CNTL(8-15) Don't exist

    inner_boundary_mumps_solver_.job = JOB_ANALYSIS_AND_FACTORIZATION;
    inner_boundary_mumps_solver_.n   = level_.grid().ntheta();
    inner_boundary_mumps_solver_.nz  = inner_boundary_matrix_nnz_; 
    inner_boundary_mumps_solver_.irn = inner_boundary_matrix_row_indices_.get();
    inner_boundary_mumps_solver_.jcn = inner_boundary_matrix_column_indices_.get();
    inner_boundary_mumps_solver_.a   = inner_boundary_matrix_values_.get();

    dmumps_c(&inner_boundary_mumps_solver_);
}

void ExtrapolatedSmootherTakeGPU::finalizeMumpsSolver()
{
    inner_boundary_mumps_solver_.job = JOB_END;
    dmumps_c(&inner_boundary_mumps_solver_);
}