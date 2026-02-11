#pragma once

#include "../smoother.h"

#include "../../LinearAlgebra/Solvers/tridiagonal_solver.h"

#ifdef GMGPOLAR_USE_MUMPS
    #include "dmumps_c.h"
    #include "mpi.h"
#endif

// The smoother implements the coupled circle-radial smoothing procedure.
// It performs iterative updates on different parts of the grid based
// on the circle/radial section of the grid and black/white line coloring.
//
// The smoother solves linear systems of the form:
//   A_sc * u_sc = f_sc − A_sc^ortho * u_sc^ortho
// where:
//   - s in {Circle, Radial} denotes the smoother section type,
//   - c in {Black, White} denotes the coloring (even/odd line sub-system).
//
// The update sequence is as follows:
//   1. Black-Circle update (u_bc):
//      A_bc * u_bc = f_bc − A_bc^ortho * u_bc^ortho
//   2. White-Circle update (u_wc):
//      A_wc * u_wc = f_wc − A_wc^ortho * u_wc^ortho
//   3. Black-Radial update (u_br):
//      A_br * u_br = f_br − A_br^ortho * u_br^ortho
//   4. White-Radial update (u_wr):
//      A_wr * u_wr = f_wr − A_wr^ortho * u_wr^ortho
//
// Algorithm details:
//   - 'rhs' corresponds to the f vector, 'x' stores the final solution,
//     and 'temp' is used for temporary storage during updates.
//   - First, temp is updated with f_sc − A_sc^ortho * u_sc^ortho.
//   - The system is then solved in-place in temp, and the results
//     are copied back to x.
//   - Using 'temp' isn't strictly necessary as all updates could be performed in place in 'x'.
//   - The stencil is applied using the A-Give method.
//
// Solver and matrix structure:
//   - The matrix A_sc is block tridiagonal due to the smoother-based
//     grid indexing, which allows efficient line-wise factorization.
//   - The inner boundary requires special handling because it
//     contains an additional across-origin coupling, making it
//     non-tridiagonal; therefore, a more general solver is used there.
//     When using the MUMPS solver, the matrix is assembled in COO format.
//     When using the in-house solver, the matrix is stored in CSR format.
//   - Circular line matrices are cyclic tridiagonal due to angular
//     periodicity, whereas radial line matrices are strictly tridiagonal.
//   - Dirichlet boundary contributions in radial matrices are shifted
//     into the right-hand side to make A_sc symmetric.

class SmootherGive : public Smoother
{
public:
    // Constructs the coupled circle-radial smoother.
    // Builds the A_sc smoother matrices and prepares the solvers.
    explicit SmootherGive(const PolarGrid& grid, const LevelCache& level_cache, const DomainGeometry& domain_geometry,
                          const DensityProfileCoefficients& density_profile_coefficients, bool DirBC_Interior,
                          int num_omp_threads);

    // If MUMPS is enabled, this cleans up the inner boundary solver.
    ~SmootherGive() override;

    // Performs one full coupled smoothing sweep:
    //   BC -> WC -> BR -> WR
    // using temp as RHS workspace.
    void smoothing(Vector<double> x, ConstVector<double> rhs, Vector<double> temp) override;

private:
    // Sequential implementation to showcase simplified version without complex parallelization patterns.
    void smoothingSequential(Vector<double> x, ConstVector<double> rhs, Vector<double> temp);
    // Parallel implementation using OpenMP for loops.
    // Scedule every 2nd/4th line in parallel to avoid race conditions arising from the A-Give distribution.
    void smoothingForLoop(Vector<double> x, ConstVector<double> rhs, Vector<double> temp);

    /* ------------------- */
    /* Tridiagonal solvers */
    /* ------------------- */

    // Batched solver for cyclic-tridiagonal circle line A_sc matrices.
    BatchedTridiagonalSolver<double> circle_tridiagonal_solver_;

    // Batched solver for tridiagonal radial circle line A_sc matrices.
    BatchedTridiagonalSolver<double> radial_tridiagonal_solver_;

    // The A_sc matrix on i_r = 0 (inner circle) is NOT tridiagonal because
    // it potentially includes across-origin coupling. Therefore, it is assembled
    // into a sparse matrix and solved using a general-purpose sparse solver.
    // When using the MUMPS solver, the matrix is assembled in COO format.
    // When using the in-house solver, the matrix is stored in CSR format.
#ifdef GMGPOLAR_USE_MUMPS
    using MatrixType = SparseMatrixCOO<double>;
    DMUMPS_STRUC_C inner_boundary_mumps_solver_;
#else
    using MatrixType = SparseMatrixCSR<double>;
    SparseLUSolver<double> inner_boundary_lu_solver_;
#endif
    // Sparse matrix for the non-tridiagonal inner boundary circle block.
    MatrixType inner_boundary_circle_matrix_;

    // Note:
    //   - circle_tridiagonal_solver_[batch=0] is unused. Use the COO/CSR matrix instead.
    //   - circle_tridiagonal_solver_[batch=i_r] solves circle line i_r.
    //   - radial_tridiagonal_solver_[batch=i_theta] solves radial line i_theta.

    /* ------------------- */
    /* Stencil definitions */
    /* ------------------- */

    // Stencils encode neighborhood connectivity for A_sc matrix assembly.
    // It is only used in the construction of COO/CSR matrices.
    // Thus it is only used for the interior boundary matrix and not needed for the tridiagonal matrices.
    // The Stencil class stores the offset for each position.
    // - Non-zero matrix indicesare obtained via `ptr + offset`
    // - A offset value of `-1` means the position is not included in the stencil pattern.
    // - Other values (0, 1, 2, ..., stencil_size - 1) correspond to valid stencil indices.

    // clang-format off
    const Stencil stencil_DB_ = {
        -1, -1, -1,
        -1,  0, -1,
        -1, -1, -1
    };
    const Stencil circle_stencil_across_origin_ = {
        -1,  3, -1,
        1,  0, -1,
        -1,  2, -1
    };

    // clang-format on

    // Select correct stencil depending on the grid position.
    const Stencil& getStencil(int i_r) const; /* Only i_r = 0 implemented */
    // Number of nonzero A_sc entries.
    int getNonZeroCountCircleAsc(int i_r) const; /* Only i_r = 0 implemented */
    // Obtain a ptr to index into COO matrices.
    // It accumulates all stencil sizes within a line up to, but excluding the current node.
    int getCircleAscIndex(int i_r, int i_theta) const; /* Only i_r = 0 implemented */

    /* --------------- */
    /* Matrix assembly */
    /* --------------- */

    // Build all A_sc matrices for circle and radial smoothers.
    void buildAscMatrices();
    // Build A_sc matrix block for a single circular line.
    void buildAscCircleSection(int i_r);
    // Build A_sc matrix block for a single radial line.
    void buildAscRadialSection(int i_theta);
    // Build A_sc for a specific node (i_r, i_theta)
    void nodeBuildAscGive(int i_r, int i_theta, const PolarGrid& grid, bool DirBC_Interior,
                          MatrixType& inner_boundary_circle_matrix,
                          BatchedTridiagonalSolver<double>& circle_tridiagonal_solver,
                          BatchedTridiagonalSolver<double>& radial_tridiagonal_solver, double arr, double att,
                          double art, double detDF, double coeff_beta);

    /* ---------------------- */
    /* Orthogonal application */
    /* ---------------------- */

    // Compute temp = f_sc − A_sc^ortho * u_sc^ortho   (precomputed right-hand side)
    // where x = u_sc and rhs = f_sc
    void applyAscOrthoCircleSection(int i_r, const SmootherColor smoother_color, ConstVector<double> x,
                                    ConstVector<double> rhs, Vector<double> temp);
    void applyAscOrthoRadialSection(int i_theta, const SmootherColor smoother_color, ConstVector<double> x,
                                    ConstVector<double> rhs, Vector<double> temp);

    /* ----------------- */
    /* Line-wise solvers */
    /* ----------------- */

    // Solve the linear system:
    //     A_sc * u_sc = f_sc − A_sc^ortho * u_sc^ortho
    // Parameter mapping:
    //   x    = u_sc   (solution vector for section s and color c)
    //   temp = f_sc − A_sc^ortho * u_sc^ortho   (precomputed right-hand side)
    // where:
    //   s in {Circle, Radial}  denotes the smoother section type,
    //   c in {Black, White}    denotes the line coloring.
    void solveBlackCircleSection(Vector<double> x, Vector<double> temp);
    void solveWhiteCircleSection(Vector<double> x, Vector<double> temp);
    void solveBlackRadialSection(Vector<double> x, Vector<double> temp);
    void solveWhiteRadialSection(Vector<double> x, Vector<double> temp);

    /* ----------------------------------- */
    /* Initialize and destroy MUMPS solver */
    /* ----------------------------------- */
#ifdef GMGPOLAR_USE_MUMPS
    // Initialize sparse MUMPS solver with assembled COO matrix.
    void initializeMumpsSolver(DMUMPS_STRUC_C& mumps_solver, SparseMatrixCOO<double>& solver_matrix);
    // Release MUMPS internal memory and MPI structures.
    void finalizeMumpsSolver(DMUMPS_STRUC_C& mumps_solver);
#endif
};
