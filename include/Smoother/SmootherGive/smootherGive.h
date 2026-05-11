#pragma once

#include "../smoother.h"

namespace gmgpolar
{

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

template <class LevelCacheType>
class SmootherGive : public Smoother<LevelCacheType>
{
public:
    // Constructs the coupled circle-radial smoother.
    // Builds the A_sc smoother matrices and prepares the solvers.
    explicit SmootherGive(const PolarGrid& grid, const LevelCacheType& level_cache, bool DirBC_Interior,
                          int num_omp_threads);

    // Performs one full coupled smoothing sweep:
    //   BC -> WC -> BR -> WR
    // Parallel implementation using OpenMP:
    // Scedule every 2nd/4th line in parallel to avoid race conditions arising from the A-Give distribution.
    // Sceduling every 3rd line in parallel would also be possible, but is less natural for the 2 coloring.
    void smoothing(Vector<double> x, ConstVector<double> rhs, Vector<double> temp) override;

private:
    /* ------------------- */
    /* Stencil definitions */
    /* ------------------- */

    // The stencil definitions must be defined before the declaration of the inner_boundary_mumps_solver_,
    // since the mumps solver will be built in the member initializer of the Smoother class.

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

    /* ------------------- */
    /* Tridiagonal solvers */
    /* ------------------- */

    // Batched solver for cyclic-tridiagonal circle line A_sc matrices.
    BatchedTridiagonalSolver<double> circle_tridiagonal_solver_;

    // Batched solver for tridiagonal radial line A_sc matrices.
    BatchedTridiagonalSolver<double> radial_tridiagonal_solver_;

    // Note:
    //   - circle_tridiagonal_solver_[batch=0] is unused. Use the COO/CSR matrix instead.
    //   - circle_tridiagonal_solver_[batch=i_r] solves circle line i_r.
    //   - radial_tridiagonal_solver_[batch=i_theta] solves radial line i_theta.

    /* ------------------------ */
    /* Interior boundary solver */
    /* ------------------------ */

    // The inner circle matrix (i_r = 0) is NOT tridiagonal due to across-origin coupling.
    // It is solved using a general-purpose sparse solver.
    // - MUMPS: matrix assembled in COO format; solver owns the matrix internally.
    // - In-house: matrix stored in CSR; solver does not own the matrix.

#ifdef GMGPOLAR_USE_MUMPS
    using InnerBoundaryMatrix = SparseMatrixCOO<double>;
    using InnerBoundarySolver = CooMumpsSolver;
#else
    using InnerBoundaryMatrix = SparseMatrixCSR<double>;
    using InnerBoundarySolver = SparseLUSolver<double>;

    // Stored only for the in-house solver (CSR).
    InnerBoundaryMatrix inner_boundary_circle_matrix_;
#endif

    // Solver object (owns matrix if MUMPS, references if in-house solver).
    InnerBoundarySolver inner_boundary_solver_;

    /* -------------- */
    /* Stencil access */
    /* -------------- */

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
    void buildTridiagonalSolverMatrices();
    void buildTridiagonalCircleSection(int i_r);
    void buildTridiagonalRadialSection(int i_theta);
    // Build the tridiagonal solver matrices for a specific node (i_r, i_theta)
    void nodeBuildTridiagonalSolverMatrices(int i_r, int i_theta, const PolarGrid& grid, bool DirBC_Interior,
                                            BatchedTridiagonalSolver<double>& circle_tridiagonal_solver,
                                            BatchedTridiagonalSolver<double>& radial_tridiagonal_solver, double arr,
                                            double att, double art, double detDF, double coeff_beta);

    // Build the solver matrix for the interior boundary (i_r = 0) which is non-tridiagonal due to across-origin coupling.
    InnerBoundaryMatrix buildInteriorBoundarySolverMatrix();
    // Build the solver matrix for a specific node (i_r = 0, i_theta) on the interior boundary.
    void nodeBuildInteriorBoundarySolverMatrix_i_r_0(int i_theta, const PolarGrid& grid, bool DirBC_Interior,
                                                     InnerBoundaryMatrix& matrix, double arr, double att, double art,
                                                     double detDF, double coeff_beta);
    void nodeBuildInteriorBoundarySolverMatrix_i_r_1(int i_theta, const PolarGrid& grid, bool DirBC_Interior,
                                                     InnerBoundaryMatrix& matrix, double arr, double att, double art,
                                                     double detDF, double coeff_beta);

    /* ---------------------- */
    /* Orthogonal application */
    /* ---------------------- */

    // Compute temp = f_sc − A_sc^ortho * u_sc^ortho   (precomputed right-hand side)
    // where x = u_sc and rhs = f_sc
    void applyAscOrthoCircleSection(const int i_r, const SmootherColor smoother_color, ConstVector<double> x,
                                    ConstVector<double> rhs, Vector<double> temp);
    void applyAscOrthoRadialSection(const int i_theta, const SmootherColor smoother_color, ConstVector<double> x,
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
};

#include "smootherGive.inl"
#include "buildInnerBoundaryAsc.inl"
#include "buildTridiagonalAsc.inl"
#include "applyAscOrtho.inl"
#include "solveAscSystem.inl"
#include "matrixStencil.inl"

} // namespace gmgpolar