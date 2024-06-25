#include "../../include/Smoother/smoother.h"

void Smoother::initializeMumps(DMUMPS_STRUC_C& Asc_mumps, const SparseMatrix<double>& Asc_matrix){
    Asc_mumps.job = JOB_INIT;
    Asc_mumps.par = 1;
    Asc_mumps.sym = (Asc_matrix.is_symmetric() ? 1 : 0);
    Asc_mumps.comm_fortran = USE_COMM_WORLD;
    dmumps_c(&Asc_mumps);
    
    Asc_mumps.ICNTL(1) =  0; // Output stream for error messages.
    Asc_mumps.ICNTL(2) =  0; // Output stream for diagnostic printing and statistics local to each MPI process.
    Asc_mumps.ICNTL(3) =  0; // Output stream for global information, collected on the host
    Asc_mumps.ICNTL(4) =  0; // Level of printing for error, warning, and diagnostic messages.
    Asc_mumps.ICNTL(5) =  0; // Controls the matrix input format
    Asc_mumps.ICNTL(6)  = 7; // Permutes the matrix to a zero-free diagonal and/or scale the matrix 
    Asc_mumps.ICNTL(7) =  5; // Computes a symmetric permutation (ordering) to determine the pivot order to be used for the 
    //                      factorization in case of sequential analysis
    Asc_mumps.ICNTL(8) = 77; // Describes the scaling strategy
    Asc_mumps.ICNTL(9) =  1; // Computes the solution using A or A^T
    Asc_mumps.ICNTL(10) = 0; // Applies the iterative refinement to the computed solution
    Asc_mumps.ICNTL(11) = 0; // Computes statistics related to an error analysis of the linear system solved
    Asc_mumps.ICNTL(12) = 0; // Defines an ordering strategy for symmetric matrices and is used
    Asc_mumps.ICNTL(13) = 0; // Controls the parallelism of the root node
    Asc_mumps.ICNTL(14) =    // Controls the percentage increase in the estimated working space
        (Asc_matrix.is_symmetric() ? 5 : 20); 
    Asc_mumps.ICNTL(15) = 0; // Exploits compression of the input matrix resulting from a block format
    Asc_mumps.ICNTL(16) = 0; // Controls the setting of the number of OpenMP threads
    // ICNTL(17) Doesn't exist
    Asc_mumps.ICNTL(18) = 0; // Defines the strategy for the distributed input matrix
    Asc_mumps.ICNTL(19) = 0; // Computes the Schur complement matrix
    Asc_mumps.ICNTL(20) = 0; // Determines the format (dense, sparse, or distributed) of the right-hand sides
    Asc_mumps.ICNTL(21) = 0; // Determines the distribution (centralized or distributed) of the solution vectors.
    Asc_mumps.ICNTL(22) = 0; // Controls the in-core/out-of-core (OOC) factorization and solve.
    Asc_mumps.ICNTL(23) = 0; // Corresponds to the maximum size of the working memory in MegaBytes that MUMPS can
    //                      allocate per working process
    Asc_mumps.ICNTL(24) = 0; // Controls the detection of “null pivot rows”.
    Asc_mumps.ICNTL(25) = 0; // Allows the computation of a solution of a deficient matrix and also of a null space basis
    Asc_mumps.ICNTL(26) = 0; // Drives the solution phase if a Schur complement matrix has been computed
    Asc_mumps.ICNTL(27) = -32; // Controls the blocking size for multiple right-hand sides.
    Asc_mumps.ICNTL(28) = 1; // Determines whether a sequential or parallel computation of the ordering is performed
    // 0: automatic, 1: sequential, 2: parallel
    Asc_mumps.ICNTL(29) = 0; // Defines the parallel ordering tool (when ICNTL(28)=1) to be used to compute the fill-in reducing permutation.
    Asc_mumps.ICNTL(30) = 0; // Computes a user-specified set of entries in the inverse A^−1 of the original matrix
    Asc_mumps.ICNTL(31) = 0; // Indicates which factors may be discarded during the factorization.
    Asc_mumps.ICNTL(32) = 0; // Performs the forward elimination of the right-hand sides during the factorization
    Asc_mumps.ICNTL(33) = 0; // Computes the determinant of the input matrix.
    Asc_mumps.ICNTL(34) = 0; // Controls the conservation of the OOC files during JOB= –3
    Asc_mumps.ICNTL(35) = 0; // Controls the activation of the BLR feature
    Asc_mumps.ICNTL(36) = 0; // Controls the choice of BLR factorization variant
    Asc_mumps.ICNTL(37) = 0; // Controls the BLR compression of the contribution blocks
    Asc_mumps.ICNTL(38) = 600; // Estimates compression rate of LU factors
    Asc_mumps.ICNTL(39) = 500; // Estimates compression rate of contribution blocks
    // ICNTL(40-47) Don't exist
    Asc_mumps.ICNTL(48) = 1; // Multithreading with tree parallelism
    Asc_mumps.ICNTL(49) = 0; // Compact workarray id%S at the end of factorization phase
    // ICNTL(50-55) Don't exist
    Asc_mumps.ICNTL(56) = 0; // Detects pseudo-singularities during factorization and factorizes the root node with a rankrevealing method
    // ICNTL(57) Doesn't exist
    Asc_mumps.ICNTL(58) = 2; // Defines options for symbolic factorization
    // ICNTL(59-60) Don't exist

    Asc_mumps.CNTL(1) = -1.0 ; // Relative threshold for numerical pivoting
    Asc_mumps.CNTL(2) = -1.0 ; // Stopping criterion for iterative refinement
    Asc_mumps.CNTL(3) = 0.0 ; // Determine null pivot rows
    Asc_mumps.CNTL(4) = -1.0 ; // Determines the threshold for static pivoting
    Asc_mumps.CNTL(5) = 0.0 ; // Defines the fixation for null pivots and is effective only when null pivot row detection is active
    // CNTL(6) Doesn't exist 
    Asc_mumps.CNTL(7) = 0.0 ; // Defines the precision of the dropping parameter used during BLR compression
    // CNTL(8-15) Don't exist

    Asc_mumps.job = JOB_ANALYSIS_AND_FACTORIZATION; 
    assert(Asc_matrix.rows() == Asc_matrix.columns());
    Asc_mumps.n = Asc_matrix.rows();
    Asc_mumps.nz = Asc_matrix.non_zero_size();
    Asc_mumps.irn = Asc_matrix.row_indices_data();
    Asc_mumps.jcn = Asc_matrix.column_indices_data();
    Asc_mumps.a = Asc_matrix.values_data();
    dmumps_c(&Asc_mumps);   
}

void Smoother::deleteMumps(DMUMPS_STRUC_C& Asc_mumps){
    Asc_mumps.job = JOB_END;
    dmumps_c(&Asc_mumps);
}