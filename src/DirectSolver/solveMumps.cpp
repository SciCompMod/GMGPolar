#include "../../include/DirectSolver/directSolver.h"

void DirectSolver::initializeMumps(DMUMPS_STRUC_C& mumps, const SparseMatrix<double>& matrixA){
    mumps.job = JOB_INIT;
    mumps.par = 1;
    mumps.sym = (matrixA.is_symmetric() ? 1 : 0);
    mumps.comm_fortran = USE_COMM_WORLD;
    dmumps_c(&mumps);

    mumps.ICNTL(1) =  0; // Output stream for error messages.
    mumps.ICNTL(2) =  0; // Output stream for diagnostic printing and statistics local to each MPI process.
    mumps.ICNTL(3) =  0; // Output stream for global information, collected on the host
    mumps.ICNTL(4) =  0; // Level of printing for error, warning, and diagnostic messages.
    mumps.ICNTL(5) =  0; // Controls the matrix input format
    mumps.ICNTL(6)  = 7; // Permutes the matrix to a zero-free diagonal and/or scale the matrix 
    mumps.ICNTL(7) =  5; // Computes a symmetric permutation (ordering) to determine the pivot order to be used for the 
    //                      factorization in case of sequential analysis
    mumps.ICNTL(8) = 77; // Describes the scaling strategy
    mumps.ICNTL(9) =  1; // Computes the solution using A or A^T
    mumps.ICNTL(10) = 0; // Applies the iterative refinement to the computed solution
    mumps.ICNTL(11) = 0; // Computes statistics related to an error analysis of the linear system solved
    mumps.ICNTL(12) = 0; // Defines an ordering strategy for symmetric matrices and is used
    mumps.ICNTL(13) = 0; // Controls the parallelism of the root node
    mumps.ICNTL(14) =    // Controls the percentage increase in the estimated working space
        (matrixA.is_symmetric() ? 5 : 20); 
    mumps.ICNTL(15) = 0; // Exploits compression of the input matrix resulting from a block format
    mumps.ICNTL(16) = 0; // Controls the setting of the number of OpenMP threads
    // ICNTL(17) Doesn't exist
    mumps.ICNTL(18) = 0; // Defines the strategy for the distributed input matrix
    mumps.ICNTL(19) = 0; // Computes the Schur complement matrix
    mumps.ICNTL(20) = 0; // Determines the format (dense, sparse, or distributed) of the right-hand sides
    mumps.ICNTL(21) = 0; // Determines the distribution (centralized or distributed) of the solution vectors.
    mumps.ICNTL(22) = 0; // Controls the in-core/out-of-core (OOC) factorization and solve.
    mumps.ICNTL(23) = 0; // Corresponds to the maximum size of the working memory in MegaBytes that MUMPS can
    //                      allocate per working process
    mumps.ICNTL(24) = 0; // Controls the detection of “null pivot rows”.
    mumps.ICNTL(25) = 0; // Allows the computation of a solution of a deficient matrix and also of a null space basis
    mumps.ICNTL(26) = 0; // Drives the solution phase if a Schur complement matrix has been computed
    mumps.ICNTL(27) = -32; // Controls the blocking size for multiple right-hand sides.
    mumps.ICNTL(28) = 0; // Determines whether a sequential or parallel computation of the ordering is performed
    mumps.ICNTL(29) = 0; // Defines the parallel ordering tool (when ICNTL(28)=1) to be used to compute the fill-in reducing permutation.
    mumps.ICNTL(30) = 0; // Computes a user-specified set of entries in the inverse A^−1 of the original matrix
    mumps.ICNTL(31) = 0; // Indicates which factors may be discarded during the factorization.
    mumps.ICNTL(32) = 0; // Performs the forward elimination of the right-hand sides during the factorization
    mumps.ICNTL(33) = 0; // Computes the determinant of the input matrix.
    mumps.ICNTL(34) = 0; // Controls the conservation of the OOC files during JOB= –3
    mumps.ICNTL(35) = 0; // Controls the activation of the BLR feature
    mumps.ICNTL(36) = 0; // Controls the choice of BLR factorization variant
    mumps.ICNTL(37) = 0; // Controls the BLR compression of the contribution blocks
    mumps.ICNTL(38) = 600; // Estimates compression rate of LU factors
    mumps.ICNTL(39) = 500; // Estimates compression rate of contribution blocks
    // ICNTL(40-47) Don't exist
    mumps.ICNTL(48) = 1; // Multithreading with tree parallelism
    mumps.ICNTL(49) = 0; // Compact workarray id%S at the end of factorization phase
    // ICNTL(50-55) Don't exist
    mumps.ICNTL(56) = 0; // Detects pseudo-singularities during factorization and factorizes the root node with a rankrevealing method
    // ICNTL(57) Doesn't exist
    mumps.ICNTL(58) = 2; // Defines options for symbolic factorization
    // ICNTL(59-60) Don't exist

    mumps.CNTL(1) = -1.0 ; // Relative threshold for numerical pivoting
    mumps.CNTL(2) = -1.0 ; // Stopping criterion for iterative refinement
    mumps.CNTL(3) = 0.0 ; // Determine null pivot rows
    mumps.CNTL(4) = -1.0 ; // Determines the threshold for static pivoting
    mumps.CNTL(5) = 0.0 ; // Defines the fixation for null pivots and is effective only when null pivot row detection is active
    // CNTL(6) Doesn't exist 
    mumps.CNTL(7) = 0.0 ; // Defines the precision of the dropping parameter used during BLR compression
    // CNTL(8-15) Don't exist

    mumps.job = JOB_ANALYSIS_AND_FACTORIZATION; 
    assert(matrixA.rows() == matrixA.columns());
    mumps.n = matrixA.rows();
    mumps.nz = matrixA.non_zero_size();
    mumps.irn = matrixA.row_indices_data();
    mumps.jcn = matrixA.column_indices_data();
    mumps.a = matrixA.values_data();
    dmumps_c(&mumps);
}

void DirectSolver::solveMumps(Vector<double>& result_rhs){
    mumps_.job = JOB_COMPUTE_SOLUTION;
    mumps_.nrhs = 1;
    mumps_.nz_rhs = result_rhs.size();
    mumps_.rhs = result_rhs.begin();
    mumps_.lrhs = result_rhs.size();
    dmumps_c(&mumps_);
    if (mumps_.info[0] != 0) {
        std::cerr << "Error solving the direct system: " << mumps_.info[0] << std::endl;
    }
}

void DirectSolver::deleteMumps(DMUMPS_STRUC_C& mumps){
    mumps.job = JOB_END;
    dmumps_c(&mumps);
}