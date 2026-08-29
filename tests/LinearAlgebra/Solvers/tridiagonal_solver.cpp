#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <vector>
#include <cmath>

#include <LinearAlgebra/Solvers/tridiagonal_solver.h>
#include <LinearAlgebra/Vector/vector.h>
using namespace gmgpolar;

// clang-format off
void test_non_cyclic_tridiagonal_n_4()
{
    int batch_count = 4;
    int matrix_dimension = 4;
    bool is_cyclic = false;

    BatchedTridiagonalSolver<double> solver(matrix_dimension, batch_count, is_cyclic);

    // System 1: {{2, 1, 0,0},{1,4,2,0},{0,2,6,3},{0,0,3,8}} * {{a},{b},{c},{d}} = {{1},{2},{3},{4}}
    // a = 70/209, b = 69/209, c = 36/209, d = 91/209
    // System 2: {{3,1,0,0},{1,5,2,0},{0,2,7,3},{0,0,3,9}} * {{a},{b},{c},{d}} = {{2},{3},{4},{5}}
    // a = 29/54, b = 7/18, c = 7/27, d = 38/81
    // System 3: {{4,1,0,0},{1,6,2,0},{0,2,8,3},{0,0,3,10}} * {{a},{b},{c},{d}} = {{3},{4},{5},{6}}
    // a = 938/1473, b = 667/1473, c = 476/1473, d = 247/491
    // System 4: {{5,1,0,0},{1,7,2,0},{0,2,9,3},{0,0,3,11}} * {{a},{b},{c},{d}} = {{4},{5},{6},{7}}
    // a = 248/355, b = 36/71, c = 267/710, d = 379/710

	Kokkos::parallel_for(
					"Test",
					1,
					KOKKOS_LAMBDA(const int) {
    solver.set_main_diagonal(0,0, 2.0); solver.set_sub_diagonal(0,0, 1.0);
    solver.set_main_diagonal(0,1, 4.0); solver.set_sub_diagonal(0,1, 2.0);
    solver.set_main_diagonal(0,2, 6.0); solver.set_sub_diagonal(0,2, 3.0);
    solver.set_main_diagonal(0,3, 8.0);

    solver.set_main_diagonal(1,0, 3.0); solver.set_sub_diagonal(1,0, 1.0);
    solver.set_main_diagonal(1,1, 5.0); solver.set_sub_diagonal(1,1, 2.0);
    solver.set_main_diagonal(1,2, 7.0); solver.set_sub_diagonal(1,2, 3.0);
    solver.set_main_diagonal(1,3, 9.0);

    solver.set_main_diagonal(2,0, 4.0); solver.set_sub_diagonal(2,0, 1.0);
    solver.set_main_diagonal(2,1, 6.0); solver.set_sub_diagonal(2,1, 2.0);
    solver.set_main_diagonal(2,2, 8.0); solver.set_sub_diagonal(2,2, 3.0);
    solver.set_main_diagonal(2,3, 10.0);

    solver.set_main_diagonal(3,0, 5.0); solver.set_sub_diagonal(3,0, 1.0);
    solver.set_main_diagonal(3,1, 7.0); solver.set_sub_diagonal(3,1, 2.0);
    solver.set_main_diagonal(3,2, 9.0); solver.set_sub_diagonal(3,2, 3.0);
    solver.set_main_diagonal(3,3, 11.0);
					});

    HostVector<double> h_rhs("h_rhs", matrix_dimension * batch_count);

    // Initialize RHS for each system
    h_rhs(0) = 1.0; h_rhs(1) = 2.0; h_rhs(2) = 3.0; h_rhs(3) = 4.0;
    h_rhs(4) = 2.0; h_rhs(5) = 3.0; h_rhs(6) = 4.0; h_rhs(7) = 5.0;
    h_rhs(8) = 3.0; h_rhs(9) = 4.0; h_rhs(10) = 5.0; h_rhs(11) = 6.0;
    h_rhs(12) = 4.0; h_rhs(13) = 5.0; h_rhs(14) = 6.0; h_rhs(15) = 7.0; 

    solver.setup();

	auto rhs = Kokkos::create_mirror_view_and_copy(DefaultMemorySpace(), h_rhs);

    int offset, stride;
    // Solve each even system separately
    offset = 0; stride = 2;
    solver.solve(rhs, offset, stride);
    // Solve each odd system separately
    offset = 1; stride = 2;
    solver.solve(rhs, offset, stride);

	Kokkos::deep_copy(h_rhs, rhs);

    // Verify solutions
    double tol = 1e-12;

    EXPECT_NEAR(h_rhs(0), 70.0/209.0, tol);
    EXPECT_NEAR(h_rhs(1), 69.0/209.0, tol);
    EXPECT_NEAR(h_rhs(2), 36.0/209.0, tol);
    EXPECT_NEAR(h_rhs(3), 91.0/209.0, tol);

    EXPECT_NEAR(h_rhs(4), 29.0/54.0, tol);
    EXPECT_NEAR(h_rhs(5), 7.0/18.0, tol);
    EXPECT_NEAR(h_rhs(6), 7.0/27.0, tol);
    EXPECT_NEAR(h_rhs(7), 38.0/81.0, tol);

    EXPECT_NEAR(h_rhs(8), 938.0/1473.0, tol);
    EXPECT_NEAR(h_rhs(9), 667.0/1473.0, tol);
    EXPECT_NEAR(h_rhs(10), 476.0/1473.0, tol);
    EXPECT_NEAR(h_rhs(11), 247.0/491.0, tol);

    EXPECT_NEAR(h_rhs(12), 248.0/355.0, tol);
    EXPECT_NEAR(h_rhs(13), 36.0/71.0, tol);
    EXPECT_NEAR(h_rhs(14), 267.0/710.0, tol);
    EXPECT_NEAR(h_rhs(15), 379.0/710.0, tol);
}
TEST(BatchedTridiagonalSolvers, non_cyclic_tridiagonal_n_4)
{
		// Call a function named function due to cuda restriction
		test_non_cyclic_tridiagonal_n_4();
}

void test_cyclic_tridiagonal_n_4()
{
    int batch_count = 4;
    int matrix_dimension = 4;
    bool is_cyclic = true;

    BatchedTridiagonalSolver<double> solver(matrix_dimension, batch_count, is_cyclic);

    // System 1: {{2, 1, 0,-1},{1,4,2,0},{0,2,6,3},{-1,0,3,8}} * {{a},{b},{c},{d}} = {{1},{2},{3},{4}}
    // a = 42/67, b = 18/67, c = 10/67, d = 35/67
    // System 2: {{3,1,0,-2},{1,5,2,0},{0,2,7,3},{-2,0,3,9}} * {{a},{b},{c},{d}} = {{2},{3},{4},{5}}
    // a = 287/274, b = 89/274, c = 45/274, d = 201/274
    // System 3: {{4,1,0,-3},{1,6,2,0},{0,2,8,3},{-3,0,3,10}} * {{a},{b},{c},{d}} = {{3},{4},{5},{6}}
    // a = 1532/1113, b = 8/21, c = 188/1113, d = 51/53
    // System 4: {{5,1,0,-4},{1,7,2,0},{0,2,9,3},{-4,0,3,11}} * {{a},{b},{c},{d}} = {{4},{5},{6},{7}}
    // a = 271/162, b = 23/54, c = 14/81, d = 97/81

	Kokkos::parallel_for(
					"Test",
					1,
					KOKKOS_LAMBDA(const int) {
    solver.set_main_diagonal(0,0, 2.0); solver.set_sub_diagonal(0,0, 1.0);
    solver.set_main_diagonal(0,1, 4.0); solver.set_sub_diagonal(0,1, 2.0);
    solver.set_main_diagonal(0,2, 6.0); solver.set_sub_diagonal(0,2, 3.0);
    solver.set_main_diagonal(0,3, 8.0); solver.set_cyclic_corner(0, -1.0);

    solver.set_main_diagonal(1,0, 3.0); solver.set_sub_diagonal(1,0, 1.0);
    solver.set_main_diagonal(1,1, 5.0); solver.set_sub_diagonal(1,1, 2.0);
    solver.set_main_diagonal(1,2, 7.0); solver.set_sub_diagonal(1,2, 3.0);
    solver.set_main_diagonal(1,3, 9.0); solver.set_cyclic_corner(1, -2.0);

    solver.set_main_diagonal(2,0, 4.0); solver.set_sub_diagonal(2,0, 1.0);
    solver.set_main_diagonal(2,1, 6.0); solver.set_sub_diagonal(2,1, 2.0);
    solver.set_main_diagonal(2,2, 8.0); solver.set_sub_diagonal(2,2, 3.0);
    solver.set_main_diagonal(2,3, 10.0); solver.set_cyclic_corner(2, -3.0);

    solver.set_main_diagonal(3,0, 5.0); solver.set_sub_diagonal(3,0, 1.0);
    solver.set_main_diagonal(3,1, 7.0); solver.set_sub_diagonal(3,1, 2.0);
    solver.set_main_diagonal(3,2, 9.0); solver.set_sub_diagonal(3,2, 3.0);
    solver.set_main_diagonal(3,3, 11.0); solver.set_cyclic_corner(3, -4.0);
					});

    HostVector<double> h_rhs("h_rhs", matrix_dimension * batch_count);

    // Initialize RHS for each system
    h_rhs(0) = 1.0; h_rhs(1) = 2.0; h_rhs(2) = 3.0; h_rhs(3) = 4.0;
    h_rhs(4) = 2.0; h_rhs(5) = 3.0; h_rhs(6) = 4.0; h_rhs(7) = 5.0;
    h_rhs(8) = 3.0; h_rhs(9) = 4.0; h_rhs(10) = 5.0; h_rhs(11) = 6.0;
    h_rhs(12) = 4.0; h_rhs(13) = 5.0; h_rhs(14) = 6.0; h_rhs(15) = 7.0; 

	auto rhs = Kokkos::create_mirror_view_and_copy(DefaultMemorySpace(), h_rhs);

    solver.setup();

    int offset, stride;
    // Solve each even system separately
    offset = 0; stride = 2;
    solver.solve(rhs, offset, stride);
    // Solve each odd system separately
    offset = 1; stride = 2;
    solver.solve(rhs, offset, stride);

	Kokkos::deep_copy(h_rhs, rhs);

    // Verify solutions
    double tol = 1e-12;

    EXPECT_NEAR(h_rhs(0), 42.0/67.0, tol);
    EXPECT_NEAR(h_rhs(1), 18.0/67.0, tol);
    EXPECT_NEAR(h_rhs(2), 10.0/67.0, tol);
    EXPECT_NEAR(h_rhs(3), 35.0/67.0, tol);

    EXPECT_NEAR(h_rhs(4), 287.0/274.0, tol);
    EXPECT_NEAR(h_rhs(5), 89.0/274.0, tol);
    EXPECT_NEAR(h_rhs(6), 45.0/274.0, tol);
    EXPECT_NEAR(h_rhs(7), 201.0/274.0, tol);

    EXPECT_NEAR(h_rhs(8), 1532.0/1113.0, tol);
    EXPECT_NEAR(h_rhs(9), 8.0/21.0, tol);
    EXPECT_NEAR(h_rhs(10), 188.0/1113.0, tol);
    EXPECT_NEAR(h_rhs(11), 51.0/53.0, tol);

    EXPECT_NEAR(h_rhs(12), 271.0/162.0, tol);
    EXPECT_NEAR(h_rhs(13), 23.0/54.0, tol);
    EXPECT_NEAR(h_rhs(14), 14.0/81.0, tol);
    EXPECT_NEAR(h_rhs(15), 97.0/81.0, tol);
}
TEST(BatchedTridiagonalSolvers, cyclic_tridiagonal_n_4) {
		// Call a function named function due to cuda restriction
		test_cyclic_tridiagonal_n_4();
}

void test_non_cyclic_diagonal_n_4() {
    int batch_count = 4;
    int matrix_dimension = 4;
    bool is_cyclic = false;

    BatchedTridiagonalSolver<double> solver(matrix_dimension, batch_count, is_cyclic);

	Kokkos::parallel_for(
					"Test",
					1,
					KOKKOS_LAMBDA(const int) {
    solver.set_main_diagonal(0,0, 2.0);
    solver.set_main_diagonal(0,1, 4.0);
    solver.set_main_diagonal(0,2, 6.0);
    solver.set_main_diagonal(0,3, 8.0);

    solver.set_main_diagonal(1,0, 3.0);
    solver.set_main_diagonal(1,1, 5.0);
    solver.set_main_diagonal(1,2, 7.0);
    solver.set_main_diagonal(1,3, 9.0);

    solver.set_main_diagonal(2,0, 4.0);
    solver.set_main_diagonal(2,1, 6.0);
    solver.set_main_diagonal(2,2, 8.0);
    solver.set_main_diagonal(2,3, 10.0);

    solver.set_main_diagonal(3,0, 5.0);
    solver.set_main_diagonal(3,1, 7.0);
    solver.set_main_diagonal(3,2, 9.0);
    solver.set_main_diagonal(3,3, 11.0);
					});

    HostVector<double> h_rhs("h_rhs", matrix_dimension * batch_count);

    // Initialize RHS for each system
    h_rhs(0) = 1.0; h_rhs(1) = 2.0; h_rhs(2) = 3.0; h_rhs(3) = 4.0;
    h_rhs(4) = 2.0; h_rhs(5) = 3.0; h_rhs(6) = 4.0; h_rhs(7) = 5.0;
    h_rhs(8) = 3.0; h_rhs(9) = 4.0; h_rhs(10) = 5.0; h_rhs(11) = 6.0;
    h_rhs(12) = 4.0; h_rhs(13) = 5.0; h_rhs(14) = 6.0; h_rhs(15) = 7.0; 

	auto rhs = Kokkos::create_mirror_view_and_copy(DefaultMemorySpace(), h_rhs);

    solver.setup();

    int offset, stride;
    // Solve each even system separately
    offset = 0; stride = 2;
    solver.solve_diagonal(rhs, offset, stride);
    // Solve each odd system separately
    offset = 1; stride = 2;
    solver.solve_diagonal(rhs, offset, stride);

	Kokkos::deep_copy(h_rhs, rhs);

    // Verify solutions
    double tol = 1e-12;

    EXPECT_NEAR(h_rhs(0), 1.0/2.0, tol);
    EXPECT_NEAR(h_rhs(1), 2.0/4.0, tol);
    EXPECT_NEAR(h_rhs(2), 3.0/6.0, tol);
    EXPECT_NEAR(h_rhs(3), 4.0/8.0, tol);

    EXPECT_NEAR(h_rhs(4), 2.0/3.0, tol);
    EXPECT_NEAR(h_rhs(5), 3.0/5.0, tol);
    EXPECT_NEAR(h_rhs(6), 4.0/7.0, tol);
    EXPECT_NEAR(h_rhs(7), 5.0/9.0, tol);

    EXPECT_NEAR(h_rhs(8), 3.0/4.0, tol);
    EXPECT_NEAR(h_rhs(9), 4.0/6.0, tol);
    EXPECT_NEAR(h_rhs(10), 5.0/8.0, tol);
    EXPECT_NEAR(h_rhs(11), 6.0/10.0, tol);

    EXPECT_NEAR(h_rhs(12), 4.0/5.0, tol);
    EXPECT_NEAR(h_rhs(13), 5.0/7.0, tol);
    EXPECT_NEAR(h_rhs(14), 6.0/9.0, tol);
    EXPECT_NEAR(h_rhs(15), 7.0/11.0, tol);
}
TEST(BatchedTridiagonalSolvers, non_cyclic_diagonal_n_4) {
		// Call a function named function due to cuda restriction
		test_non_cyclic_diagonal_n_4();
}

void test_cyclic_diagonal_n_4() {
    int batch_count = 4;
    int matrix_dimension = 4;
    bool is_cyclic = true;

    BatchedTridiagonalSolver<double> solver(matrix_dimension, batch_count, is_cyclic);

	Kokkos::parallel_for(
					"Test",
					1,
					KOKKOS_LAMBDA(const int) {
    solver.set_main_diagonal(0,0, 2.0);
    solver.set_main_diagonal(0,1, 4.0);
    solver.set_main_diagonal(0,2, 6.0);
    solver.set_main_diagonal(0,3, 8.0);

    solver.set_main_diagonal(1,0, 3.0);
    solver.set_main_diagonal(1,1, 5.0);
    solver.set_main_diagonal(1,2, 7.0);
    solver.set_main_diagonal(1,3, 9.0);

    solver.set_main_diagonal(2,0, 4.0);
    solver.set_main_diagonal(2,1, 6.0);
    solver.set_main_diagonal(2,2, 8.0);
    solver.set_main_diagonal(2,3, 10.0);

    solver.set_main_diagonal(3,0, 5.0);
    solver.set_main_diagonal(3,1, 7.0);
    solver.set_main_diagonal(3,2, 9.0);
    solver.set_main_diagonal(3,3, 11.0);
					});

    HostVector<double> h_rhs("h_rhs", matrix_dimension * batch_count);

    // Initialize RHS for each system
    h_rhs(0) = 1.0; h_rhs(1) = 2.0; h_rhs(2) = 3.0; h_rhs(3) = 4.0;
    h_rhs(4) = 2.0; h_rhs(5) = 3.0; h_rhs(6) = 4.0; h_rhs(7) = 5.0;
    h_rhs(8) = 3.0; h_rhs(9) = 4.0; h_rhs(10) = 5.0; h_rhs(11) = 6.0;
    h_rhs(12) = 4.0; h_rhs(13) = 5.0; h_rhs(14) = 6.0; h_rhs(15) = 7.0; 

	auto rhs = Kokkos::create_mirror_view_and_copy(DefaultMemorySpace(), h_rhs);

    solver.setup();

    int offset, stride;
    // Solve each even system separately
    offset = 0; stride = 2;
    solver.solve_diagonal(rhs, offset, stride);
    // Solve each odd system separately
    offset = 1; stride = 2;
    solver.solve_diagonal(rhs, offset, stride);

	Kokkos::deep_copy(h_rhs, rhs);

    // Verify solutions
    double tol = 1e-12;

    EXPECT_NEAR(h_rhs(0), 1.0/2.0, tol);
    EXPECT_NEAR(h_rhs(1), 2.0/4.0, tol);
    EXPECT_NEAR(h_rhs(2), 3.0/6.0, tol);
    EXPECT_NEAR(h_rhs(3), 4.0/8.0, tol);

    EXPECT_NEAR(h_rhs(4), 2.0/3.0, tol);
    EXPECT_NEAR(h_rhs(5), 3.0/5.0, tol);
    EXPECT_NEAR(h_rhs(6), 4.0/7.0, tol);
    EXPECT_NEAR(h_rhs(7), 5.0/9.0, tol);

    EXPECT_NEAR(h_rhs(8), 3.0/4.0, tol);
    EXPECT_NEAR(h_rhs(9), 4.0/6.0, tol);
    EXPECT_NEAR(h_rhs(10), 5.0/8.0, tol);
    EXPECT_NEAR(h_rhs(11), 6.0/10.0, tol);

    EXPECT_NEAR(h_rhs(12), 4.0/5.0, tol);
    EXPECT_NEAR(h_rhs(13), 5.0/7.0, tol);
    EXPECT_NEAR(h_rhs(14), 6.0/9.0, tol);
    EXPECT_NEAR(h_rhs(15), 7.0/11.0, tol);
}
TEST(BatchedTridiagonalSolvers, cyclic_diagonal_n_4) {
		// Call a function named function due to cuda restriction
		test_cyclic_diagonal_n_4();
}


// -----------------------------------------------------------------------------------------------
// Additional edge-case coverage.
//
// The four tests above only ever exercise matrix_dimension_ == 4 (a clean power of two, well
// above every boundary-clamping/degenerate special case) with batch_count_ == matrix_dimension_
// and a perfectly even 50/50 stride split. None of them touch:
//   - matrix_dimension_ == 1 (the early-return branch in setup(), and the RangePolicy-only branch
//     in solve() — literally never entered by any existing test)
//   - matrix_dimension_ == 2 (every PCR step clamps BOTH neighbors to the opposite boundary)
//   - non-power-of-two matrix_dimension_ (num_steps_ = ceil(log2(n)) taking a non-exact value)
//   - batch_count_ == 1 (TeamPolicy league_size == 1)
//   - batch_count_ > matrix_dimension_ (the inverse of the "few, long lines" regime the file's own
//     header comments assume)
//   - matrix_dimension_ large enough that Kokkos::AUTO's chosen team_size is very likely < n,
//     forcing the `for (i = rank; i < n; i += team_size)` strided loops to actually iterate more
//     than once per thread on at least some backends
//   - batch_offset_/batch_stride_ combinations that don't split batch_count_ evenly
//   - calling solve()/solve_diagonal() with their default arguments at all
//
// Since hand-deriving exact fractions stops being practical past n=4, these tests check against a
// dense Gaussian-elimination reference solver instead. It shares no code or algorithmic structure
// with either the old Thomas solver or this PCR solver — it validates the actual linear system,
// not just that one tridiagonal-specific algorithm agrees with another.
// -----------------------------------------------------------------------------------------------

#include <algorithm>
#include <utility>

// Deterministic (not random) but distinct-per-(batch,index) diagonally dominant system generator,
// shared between "build the solver's input" and "build the reference's input" below, so both
// sides are guaranteed to describe the same matrix without needing a host<->device data transfer
// of arbitrary-sized arrays into a KOKKOS_LAMBDA.
KOKKOS_INLINE_FUNCTION double sysDiag(int b, int i)
{
    return 8.0 + 2.0 * i + 0.5 * b; // baseline large enough to keep the corner + off-diagonals
                                    // diagonally dominant up to the largest n/batch used below
}
KOKKOS_INLINE_FUNCTION double sysSub(int b, int i)
{
    return 1.0 + 0.1 * i;
}
KOKKOS_INLINE_FUNCTION double sysCorner(int b)
{
    return 0.7 + 0.05 * b;
}
KOKKOS_INLINE_FUNCTION double sysRhs(int b, int i)
{
    return 1.0 + 0.3 * i + 0.2 * b;
}

// Fully independent reference: build the dense n x n matrix explicitly (tridiagonal, plus the
// wraparound corner entries when is_cyclic) and solve it with plain partial-pivot Gaussian
// elimination. O(n^3) is irrelevant at the sizes used in these tests.
static void denseReferenceSolve(int n, const std::vector<double>& diag, const std::vector<double>& subdiag,
                                 double corner, bool is_cyclic, std::vector<double>& rhs)
{
    std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) {
        A[i][i] = diag[i];
    }
    for (int i = 0; i < n - 1; i++) {
        A[i][i + 1] = subdiag[i];
        A[i + 1][i] = subdiag[i];
    }
    if (is_cyclic && n > 1) {
        A[0][n - 1] = corner;
        A[n - 1][0] = corner;
    }

    for (int col = 0; col < n; col++) {
        int pivot = col;
        for (int row = col + 1; row < n; row++) {
            if (std::fabs(A[row][col]) > std::fabs(A[pivot][col]))
                pivot = row;
        }
        std::swap(A[col], A[pivot]);
        std::swap(rhs[col], rhs[pivot]);

        for (int row = col + 1; row < n; row++) {
            double factor = A[row][col] / A[col][col];
            for (int c = col; c < n; c++) {
                A[row][c] -= factor * A[col][c];
            }
            rhs[row] -= factor * rhs[col];
        }
    }

    for (int row = n - 1; row >= 0; row--) {
        double sum = rhs[row];
        for (int c = row + 1; c < n; c++) {
            sum -= A[row][c] * rhs[c];
        }
        rhs[row] = sum / A[row][row];
    }
}

// Builds a solver instance sized (matrix_dimension, batch_count, is_cyclic), fills it with the
// deterministic sysDiag/sysSub/sysCorner/sysRhs data, calls setup() once then solve() with the
// given (batch_offset, batch_stride), and checks every system that offset/stride combination
// actually covers against denseReferenceSolve. Does NOT assume batch_offset/stride cover the
// whole batch — callers needing full coverage make multiple calls against the same solver (see
// test_uneven_batch_stride_split below) or pass offset=0/stride=1.
static void runCorrectnessCheck(int matrix_dimension, int batch_count, bool is_cyclic, int batch_offset,
                                 int batch_stride)
{
    BatchedTridiagonalSolver<double> solver(matrix_dimension, batch_count, is_cyclic);

    Kokkos::parallel_for(
        "Fill", 1, KOKKOS_LAMBDA(const int) {
            for (int b = 0; b < batch_count; b++) {
                for (int i = 0; i < matrix_dimension; i++) {
                    solver.set_main_diagonal(b, i, sysDiag(b, i));
                }
                for (int i = 0; i < matrix_dimension - 1; i++) {
                    solver.set_sub_diagonal(b, i, sysSub(b, i));
                }
                if (is_cyclic && matrix_dimension > 1) {
                    solver.set_cyclic_corner(b, sysCorner(b));
                }
            }
        });

    HostVector<double> h_rhs("h_rhs", matrix_dimension * batch_count);
    for (int b = 0; b < batch_count; b++) {
        for (int i = 0; i < matrix_dimension; i++) {
            h_rhs(b * matrix_dimension + i) = sysRhs(b, i);
        }
    }

    auto rhs = Kokkos::create_mirror_view_and_copy(DefaultMemorySpace(), h_rhs);
    solver.setup();
    solver.solve(rhs, batch_offset, batch_stride);
    Kokkos::deep_copy(h_rhs, rhs);

    const double tol = 1e-9;
    for (int k = 0;; k++) {
        int batch_idx = batch_stride * k + batch_offset;
        if (batch_idx >= batch_count)
            break;

        std::vector<double> diag(matrix_dimension);
        std::vector<double> subdiag(matrix_dimension > 1 ? matrix_dimension - 1 : 0);
        std::vector<double> rhs_ref(matrix_dimension);
        for (int i = 0; i < matrix_dimension; i++) {
            diag[i]    = sysDiag(batch_idx, i);
            rhs_ref[i] = sysRhs(batch_idx, i);
        }
        for (int i = 0; i < matrix_dimension - 1; i++) {
            subdiag[i] = sysSub(batch_idx, i);
        }
        double corner = sysCorner(batch_idx);

        denseReferenceSolve(matrix_dimension, diag, subdiag, corner, is_cyclic, rhs_ref);

        for (int i = 0; i < matrix_dimension; i++) {
            EXPECT_NEAR(h_rhs(batch_idx * matrix_dimension + i), rhs_ref[i], tol)
                << "n=" << matrix_dimension << " batch=" << batch_idx << " index=" << i;
        }
    }
}

// Diagonal-only counterpart: solve_diagonal() deliberately ignores off-diagonal coupling, so its
// only well-defined expected answer is plain rhs[i] / diag[i] against the ORIGINAL diagonal — the
// solver is filled with zero off-diagonal entries (including a zero corner even when is_cyclic) so
// that expectation actually holds. See setup()'s Sherman-Morrison-Woodbury b[0] adjustment: for a
// genuinely diagonal cyclic system (corner == 0) the +gamma/-gamma correction solve_diagonal()
// applies cancels back to the original diagonal exactly — a nonzero corner here would not be a
// "diagonal matrix" in the first place, so it's intentionally not tested through this path.
static void runDiagonalCorrectnessCheck(int matrix_dimension, int batch_count, bool is_cyclic, int batch_offset,
                                         int batch_stride)
{
    BatchedTridiagonalSolver<double> solver(matrix_dimension, batch_count, is_cyclic);

    Kokkos::parallel_for(
        "FillDiag", 1, KOKKOS_LAMBDA(const int) {
            for (int b = 0; b < batch_count; b++) {
                for (int i = 0; i < matrix_dimension; i++) {
                    solver.set_main_diagonal(b, i, sysDiag(b, i));
                }
            }
        });

    HostVector<double> h_rhs("h_rhs", matrix_dimension * batch_count);
    for (int b = 0; b < batch_count; b++) {
        for (int i = 0; i < matrix_dimension; i++) {
            h_rhs(b * matrix_dimension + i) = sysRhs(b, i);
        }
    }

    auto rhs = Kokkos::create_mirror_view_and_copy(DefaultMemorySpace(), h_rhs);
    solver.setup();
    solver.solve_diagonal(rhs, batch_offset, batch_stride);
    Kokkos::deep_copy(h_rhs, rhs);

    const double tol = 1e-12;
    for (int k = 0;; k++) {
        int batch_idx = batch_stride * k + batch_offset;
        if (batch_idx >= batch_count)
            break;
        for (int i = 0; i < matrix_dimension; i++) {
            double expected = sysRhs(batch_idx, i) / sysDiag(batch_idx, i);
            EXPECT_NEAR(h_rhs(batch_idx * matrix_dimension + i), expected, tol)
                << "n=" << matrix_dimension << " batch=" << batch_idx << " index=" << i;
        }
    }
}

// --- matrix_dimension_ == 1: the degenerate early-return in setup() and the RangePolicy-only
// branch in solve() — not entered by any existing test. Also checks solve() and solve_diagonal()
// agree with each other here, since for n=1 there's no off-diagonal coupling to distinguish them.
void test_trivial_matrix_dimension_1()
{
    for (bool is_cyclic : {false, true}) {
        runCorrectnessCheck(/*n=*/1, /*batch_count=*/3, is_cyclic, /*offset=*/0, /*stride=*/1);
        runDiagonalCorrectnessCheck(/*n=*/1, /*batch_count=*/3, is_cyclic, /*offset=*/0, /*stride=*/1);
    }
}
TEST(BatchedTridiagonalSolvers, trivial_matrix_dimension_1)
{
    test_trivial_matrix_dimension_1();
}

// --- matrix_dimension_ == 2: smallest nontrivial size — every PCR step's iLeft/iRight clamps to
// the opposite boundary immediately (delta=1 already reaches both ends of a length-2 system).
void test_matrix_dimension_2()
{
    for (bool is_cyclic : {false, true}) {
        runCorrectnessCheck(/*n=*/2, /*batch_count=*/3, is_cyclic, /*offset=*/0, /*stride=*/1);
    }
}
TEST(BatchedTridiagonalSolvers, matrix_dimension_2)
{
    test_matrix_dimension_2();
}

// --- Non-power-of-two matrix_dimension_: num_steps_ = ceil(log2(n)) is not exact, so the last PCR
// step operates on a system that's "overshot" past what n strictly requires. Covers odd, even, and
// values straddling a power-of-two boundary on both sides (3, 5, 6, 7 bracket 4 and 8).
void test_non_power_of_two_dimensions()
{
    for (int n : {3, 5, 6, 7}) {
        for (bool is_cyclic : {false, true}) {
            runCorrectnessCheck(n, /*batch_count=*/3, is_cyclic, /*offset=*/0, /*stride=*/1);
        }
    }
}
TEST(BatchedTridiagonalSolvers, non_power_of_two_dimensions)
{
    test_non_power_of_two_dimensions();
}

// --- batch_count_ == 1: TeamPolicy league_size == 1, and batch_count_ > matrix_dimension_: the
// inverse of the "few, long lines" regime the file's parallelization-model comments assume.
void test_batch_count_extremes()
{
    for (bool is_cyclic : {false, true}) {
        runCorrectnessCheck(/*n=*/5, /*batch_count=*/1, is_cyclic, /*offset=*/0, /*stride=*/1);
        runCorrectnessCheck(/*n=*/2, /*batch_count=*/10, is_cyclic, /*offset=*/0, /*stride=*/1);
    }
}
TEST(BatchedTridiagonalSolvers, batch_count_extremes)
{
    test_batch_count_extremes();
}

// --- matrix_dimension_ large enough that Kokkos::AUTO's chosen team_size is plausibly < n on at
// least some backends, forcing the `for (i = rank; i < n; i += team_size)` strided loops to
// actually run more than one iteration per thread — every existing test uses n=4, small enough
// that most backends would just grant team_size >= n and never exercise that loop's stride > 0
// case at all.
void test_large_dimension_strided_team_loop()
{
    for (bool is_cyclic : {false, true}) {
        runCorrectnessCheck(/*n=*/33, /*batch_count=*/2, is_cyclic, /*offset=*/0, /*stride=*/1);
    }
    runDiagonalCorrectnessCheck(/*n=*/33, /*batch_count=*/2, /*is_cyclic=*/true, /*offset=*/0, /*stride=*/1);
}
TEST(BatchedTridiagonalSolvers, large_dimension_strided_team_loop)
{
    test_large_dimension_strided_team_loop();
}

// --- Uneven batch_offset_/batch_stride_ split: every existing test splits batch_count_ evenly in
// half (offset 0/1, stride 2, batch_count_ divisible by 2). Here batch_count_ = 5 is not divisible
// by stride = 3, so effective_batch_count_'s ceiling-division arithmetic produces a genuinely
// uneven split (2, 2, 1) across three solve() calls sharing ONE setup() — also checks that calling
// solve() more than twice against the same stored trajectory still works, not just twice.
void test_uneven_batch_stride_split()
{
    const int matrix_dimension = 4;
    const int batch_count      = 5;

    for (bool is_cyclic : {false, true}) {
        BatchedTridiagonalSolver<double> solver(matrix_dimension, batch_count, is_cyclic);

        Kokkos::parallel_for(
            "Fill", 1, KOKKOS_LAMBDA(const int) {
                for (int b = 0; b < batch_count; b++) {
                    for (int i = 0; i < matrix_dimension; i++) {
                        solver.set_main_diagonal(b, i, sysDiag(b, i));
                    }
                    for (int i = 0; i < matrix_dimension - 1; i++) {
                        solver.set_sub_diagonal(b, i, sysSub(b, i));
                    }
                    if (is_cyclic) {
                        solver.set_cyclic_corner(b, sysCorner(b));
                    }
                }
            });

        HostVector<double> h_rhs("h_rhs", matrix_dimension * batch_count);
        for (int b = 0; b < batch_count; b++) {
            for (int i = 0; i < matrix_dimension; i++) {
                h_rhs(b * matrix_dimension + i) = sysRhs(b, i);
            }
        }

        auto rhs = Kokkos::create_mirror_view_and_copy(DefaultMemorySpace(), h_rhs);
        solver.setup();

        const int stride = 3;
        solver.solve(rhs, /*offset=*/0, stride); // covers batch 0, 3
        solver.solve(rhs, /*offset=*/1, stride); // covers batch 1, 4
        solver.solve(rhs, /*offset=*/2, stride); // covers batch 2

        Kokkos::deep_copy(h_rhs, rhs);

        const double tol = 1e-9;
        for (int b = 0; b < batch_count; b++) {
            std::vector<double> diag(matrix_dimension), subdiag(matrix_dimension - 1), rhs_ref(matrix_dimension);
            for (int i = 0; i < matrix_dimension; i++) {
                diag[i]    = sysDiag(b, i);
                rhs_ref[i] = sysRhs(b, i);
            }
            for (int i = 0; i < matrix_dimension - 1; i++) {
                subdiag[i] = sysSub(b, i);
            }
            denseReferenceSolve(matrix_dimension, diag, subdiag, sysCorner(b), is_cyclic, rhs_ref);

            for (int i = 0; i < matrix_dimension; i++) {
                EXPECT_NEAR(h_rhs(b * matrix_dimension + i), rhs_ref[i], tol) << "batch " << b << " index " << i;
            }
        }
    }
}
TEST(BatchedTridiagonalSolvers, uneven_batch_stride_split)
{
    test_uneven_batch_stride_split();
}

// --- Default batch_offset_/batch_stride_ arguments: every existing test passes offset/stride
// explicitly. This checks solve()/solve_diagonal() called with no arguments at all (full batch in
// a single call) still compile and behave as offset=0, stride=1.
void test_default_batch_arguments()
{
    const int matrix_dimension = 4;
    const int batch_count      = 3;

    for (bool is_cyclic : {false, true}) {
        BatchedTridiagonalSolver<double> solver(matrix_dimension, batch_count, is_cyclic);

        Kokkos::parallel_for(
            "Fill", 1, KOKKOS_LAMBDA(const int) {
                for (int b = 0; b < batch_count; b++) {
                    for (int i = 0; i < matrix_dimension; i++) {
                        solver.set_main_diagonal(b, i, sysDiag(b, i));
                    }
                    for (int i = 0; i < matrix_dimension - 1; i++) {
                        solver.set_sub_diagonal(b, i, sysSub(b, i));
                    }
                    if (is_cyclic) {
                        solver.set_cyclic_corner(b, sysCorner(b));
                    }
                }
            });

        HostVector<double> h_rhs("h_rhs", matrix_dimension * batch_count);
        for (int b = 0; b < batch_count; b++) {
            for (int i = 0; i < matrix_dimension; i++) {
                h_rhs(b * matrix_dimension + i) = sysRhs(b, i);
            }
        }

        auto rhs = Kokkos::create_mirror_view_and_copy(DefaultMemorySpace(), h_rhs);
        solver.setup();
        solver.solve(rhs); // no offset/stride passed

        Kokkos::deep_copy(h_rhs, rhs);

        const double tol = 1e-9;
        for (int b = 0; b < batch_count; b++) {
            std::vector<double> diag(matrix_dimension), subdiag(matrix_dimension - 1), rhs_ref(matrix_dimension);
            for (int i = 0; i < matrix_dimension; i++) {
                diag[i]    = sysDiag(b, i);
                rhs_ref[i] = sysRhs(b, i);
            }
            for (int i = 0; i < matrix_dimension - 1; i++) {
                subdiag[i] = sysSub(b, i);
            }
            denseReferenceSolve(matrix_dimension, diag, subdiag, sysCorner(b), is_cyclic, rhs_ref);

            for (int i = 0; i < matrix_dimension; i++) {
                EXPECT_NEAR(h_rhs(b * matrix_dimension + i), rhs_ref[i], tol) << "batch " << b << " index " << i;
            }
        }
    }
}
TEST(BatchedTridiagonalSolvers, default_batch_arguments)
{
    test_default_batch_arguments();
}