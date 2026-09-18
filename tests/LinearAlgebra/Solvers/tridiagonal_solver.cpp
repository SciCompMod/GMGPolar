#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <vector>
#include <cmath>
#include <algorithm>
#include <utility>

#include <LinearAlgebra/Solvers/tridiagonal_solver.h>

#include <LinearAlgebra/Vector/vector.h>
using namespace gmgpolar;

// clang-format off

// -----------------------------------------------------------------------------------------------
// Every test in this file is written as a function template on the concrete solver type and is
// instantiated for all four batched tridiagonal solver backends: Thomas, CR, PCR and CRPCR.
//
// Previously these tests only used `BatchedTridiagonalSolver<double>`, a type alias that picks a
// SINGLE one of the four implementations depending on the active Kokkos backend (Thomas on
// host-only builds, CRPCR when CUDA/HIP/SYCL is enabled). That means CR and PCR were never
// exercised by this file at all, and depending on the build, either Thomas or CRPCR was silently
// skipped too.
//
// Test bodies are still free functions (not TEST/TYPED_TEST member functions) for the same reason
// as before: CUDA extended device lambdas cannot be defined inside a function with internal
// linkage, so the Kokkos::parallel_for bodies need to live in ordinary, externally-linked
// functions. The INSTANTIATE_FOR_ALL_SOLVERS macro below generates the four TEST(...) cases (one
// per solver backend) that each call the templated test function with the corresponding solver
// type.
// -----------------------------------------------------------------------------------------------

#define INSTANTIATE_FOR_ALL_SOLVERS(TESTNAME)                                                                 \
    TEST(BatchedTridiagonalSolvers_Thomas, TESTNAME) { test_##TESTNAME<BatchedTridiagonalSolverThomas<double>>(); } \
    TEST(BatchedTridiagonalSolvers_CR,     TESTNAME) { test_##TESTNAME<BatchedTridiagonalSolverCR<double>>(); }     \
    TEST(BatchedTridiagonalSolvers_PCR,    TESTNAME) { test_##TESTNAME<BatchedTridiagonalSolverPCR<double>>(); }    \
    TEST(BatchedTridiagonalSolvers_CRPCR,  TESTNAME) { test_##TESTNAME<BatchedTridiagonalSolverCRPCR<double>>(); }

// =================================================================================================
// Hand-derived 4x4 systems (unchanged from before, just templated on SolverType)
// =================================================================================================

template <typename SolverType>
void test_non_cyclic_tridiagonal_n_4()
{
    int batch_count = 4;
    int matrix_dimension = 4;
    bool is_cyclic = false;

    SolverType solver(matrix_dimension, batch_count, is_cyclic);

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
INSTANTIATE_FOR_ALL_SOLVERS(non_cyclic_tridiagonal_n_4)

template <typename SolverType>
void test_cyclic_tridiagonal_n_4()
{
    int batch_count = 4;
    int matrix_dimension = 4;
    bool is_cyclic = true;

    SolverType solver(matrix_dimension, batch_count, is_cyclic);

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
INSTANTIATE_FOR_ALL_SOLVERS(cyclic_tridiagonal_n_4)

template <typename SolverType>
void test_non_cyclic_diagonal_n_4() {
    int batch_count = 4;
    int matrix_dimension = 4;
    bool is_cyclic = false;

    SolverType solver(matrix_dimension, batch_count, is_cyclic);

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
INSTANTIATE_FOR_ALL_SOLVERS(non_cyclic_diagonal_n_4)

template <typename SolverType>
void test_cyclic_diagonal_n_4() {
    int batch_count = 4;
    int matrix_dimension = 4;
    bool is_cyclic = true;

    SolverType solver(matrix_dimension, batch_count, is_cyclic);

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
INSTANTIATE_FOR_ALL_SOLVERS(cyclic_diagonal_n_4)


// -----------------------------------------------------------------------------------------------
// Additional edge-case coverage (templated on SolverType, unchanged in spirit from before).
//
// The four tests above only ever exercise matrix_dimension_ == 4 (a clean power of two, well
// above every boundary-clamping/degenerate special case) with batch_count_ == matrix_dimension_
// and a perfectly even 50/50 stride split. The tests below additionally touch:
//   - matrix_dimension_ == 1 (the early-return branch in setup(), and the RangePolicy-only branch
//     in solve())
//   - matrix_dimension_ == 2 (every PCR/CR step clamps neighbors to the opposite boundary)
//   - non-power-of-two matrix_dimension_ (num_steps_ = ceil(log2(n)) taking a non-exact value)
//   - batch_count_ == 1 (TeamPolicy league_size == 1)
//   - batch_count_ > matrix_dimension_ (the inverse of the "few, long lines" regime the file's own
//     header comments assume)
//   - matrix_dimension_ large enough that Kokkos::AUTO's chosen team_size is very likely < n,
//     forcing strided team loops to actually iterate more than once per thread
//   - batch_offset_/batch_stride_ combinations that don't split batch_count_ evenly
//   - calling solve()/solve_diagonal() with their default arguments at all
//   - matrix_dimension_ in the thousands, well beyond anything reachable via dense O(n^3)
//     Gaussian-elimination verification
//
// Since hand-deriving exact fractions stops being practical past n=4, these tests check against
// independently-implemented reference solvers instead. They share no code with any of the four
// solver classes under test - only the mathematics of solving the actual linear system.
// -----------------------------------------------------------------------------------------------

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

// Fully independent reference for small/moderate n: build the dense n x n matrix explicitly
// (tridiagonal, plus the wraparound corner entries when is_cyclic) and solve it with plain
// partial-pivot Gaussian elimination. O(n^3) is irrelevant at the sizes used for these tests
// (up to ~1024) but far too slow for the n=4000-scale tests further below.
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
        A[0][n - 1] += corner;
        A[n - 1][0] += corner;
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

// Independent O(n) banded reference solver, used only for the very-large-matrix_dimension_ tests
// where the O(n^3) dense Gaussian elimination above would be impractically slow (n=4000 => ~6.4e10
// flops per system). Any exact direct solve of a tridiagonal system is necessarily some variant of
// forward-elimination/back-substitution - there is no way to be algorithmically "independent" of
// that shape at O(n) - but this implementation is written from scratch here and does not call,
// share code with, or share any intermediate representation with any of the four solver classes
// under test (Thomas/CR/PCR/CRPCR), so it still independently catches bugs in each of them.
static void bandedReferenceSolveNonCyclic(int n, const std::vector<double>& diag, const std::vector<double>& subdiag,
                                           std::vector<double>& rhs)
{
    std::vector<double> c_prime(n, 0.0), d_prime(n, 0.0);
    c_prime[0] = (n > 1) ? subdiag[0] / diag[0] : 0.0;
    d_prime[0] = rhs[0] / diag[0];
    for (int i = 1; i < n; i++) {
        double denom = diag[i] - subdiag[i - 1] * c_prime[i - 1];
        c_prime[i]   = (i < n - 1) ? subdiag[i] / denom : 0.0;
        d_prime[i]   = (rhs[i] - subdiag[i - 1] * d_prime[i - 1]) / denom;
    }
    rhs[n - 1] = d_prime[n - 1];
    for (int i = n - 2; i >= 0; i--) {
        rhs[i] = d_prime[i] - c_prime[i] * rhs[i + 1];
    }
}

// Cyclic counterpart via the standard Sherman-Morrison technique for bordered tridiagonal systems
// (turns one cyclic solve into two plain tridiagonal solves plus a rank-1 correction). n == 2 is
// special-cased because there the "corner" position (0,1)/(1,0) coincides with the regular
// off-diagonal position rather than being a separate matrix entry (matching how denseReferenceSolve
// and the solver classes under test build their n==2 cyclic matrix).
static void bandedReferenceSolveCyclic(int n, const std::vector<double>& diag, const std::vector<double>& subdiag,
                                        double corner, std::vector<double>& rhs)
{
    if (n == 1) {
        rhs[0] = rhs[0] / (diag[0] + 2.0 * corner);
        return;
    }
    if (n == 2) {
        double a00  = diag[0];
        double aoff = subdiag[0] + corner;
        double a11  = diag[1];
        double det  = a00 * a11 - aoff * aoff;
        double r0 = rhs[0], r1 = rhs[1];
        rhs[0] = (r0 * a11 - aoff * r1) / det;
        rhs[1] = (a00 * r1 - aoff * r0) / det;
        return;
    }

    double alpha = corner;
    double beta  = corner;
    double gamma = -diag[0];

    std::vector<double> bb = diag;
    bb[0] -= gamma;
    bb[n - 1] -= alpha * beta / gamma;

    std::vector<double> u(n, 0.0);
    u[0]     = gamma;
    u[n - 1] = alpha;

    std::vector<double> x = rhs;
    bandedReferenceSolveNonCyclic(n, bb, subdiag, x);

    std::vector<double> z = u;
    bandedReferenceSolveNonCyclic(n, bb, subdiag, z);

    double fact = (x[0] + beta * x[n - 1] / gamma) / (1.0 + z[0] + beta * z[n - 1] / gamma);

    for (int i = 0; i < n; i++) {
        rhs[i] = x[i] - fact * z[i];
    }
}

// Builds a solver instance sized (matrix_dimension, batch_count, is_cyclic), fills it with the
// deterministic sysDiag/sysSub/sysCorner/sysRhs data, calls setup() once then solve() with the
// given (batch_offset, batch_stride), and checks every system that offset/stride combination
// actually covers against denseReferenceSolve. Does NOT assume batch_offset/stride cover the
// whole batch - callers needing full coverage make multiple calls against the same solver (see
// test_uneven_batch_stride_split below) or pass offset=0/stride=1. Intended for small/moderate n
// (dense O(n^3) reference) - use runHugeCorrectnessCheck for n in the thousands.
template <typename SolverType>
static void runCorrectnessCheck(int matrix_dimension, int batch_count, bool is_cyclic, int batch_offset,
                                 int batch_stride)
{
    SolverType solver(matrix_dimension, batch_count, is_cyclic);

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

// Same idea as runCorrectnessCheck but for matrix_dimension_ in the thousands: uses the O(n)
// bandedReferenceSolve{Non}Cyclic reference instead of the O(n^3) dense one.
template <typename SolverType>
static void runHugeCorrectnessCheck(int matrix_dimension, int batch_count, bool is_cyclic, int batch_offset,
                                     int batch_stride)
{
    SolverType solver(matrix_dimension, batch_count, is_cyclic);

    Kokkos::parallel_for(
        "FillHuge", 1, KOKKOS_LAMBDA(const int) {
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

    const double tol = 1e-6;
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

        if (is_cyclic) {
            bandedReferenceSolveCyclic(matrix_dimension, diag, subdiag, sysCorner(batch_idx), rhs_ref);
        } else {
            bandedReferenceSolveNonCyclic(matrix_dimension, diag, subdiag, rhs_ref);
        }

        for (int i = 0; i < matrix_dimension; i++) {
            EXPECT_NEAR(h_rhs(batch_idx * matrix_dimension + i), rhs_ref[i], tol)
                << "n=" << matrix_dimension << " batch=" << batch_idx << " index=" << i;
        }
    }
}

// Diagonal-only counterpart: solve_diagonal() deliberately ignores off-diagonal coupling, so its
// only well-defined expected answer is plain rhs[i] / diag[i] against the ORIGINAL diagonal - the
// solver is filled with zero off-diagonal entries (including a zero corner even when is_cyclic) so
// that expectation actually holds. See setup()'s Sherman-Morrison-Woodbury b[0] adjustment: for a
// genuinely diagonal cyclic system (corner == 0) the +gamma/-gamma correction solve_diagonal()
// applies cancels back to the original diagonal exactly - a nonzero corner here would not be a
// "diagonal matrix" in the first place, so it's intentionally not tested through this path.
template <typename SolverType>
static void runDiagonalCorrectnessCheck(int matrix_dimension, int batch_count, bool is_cyclic, int batch_offset,
                                         int batch_stride)
{
    SolverType solver(matrix_dimension, batch_count, is_cyclic);

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
// branch in solve() - not entered by the hand-derived n=4 tests. Also checks solve() and
// solve_diagonal() agree with each other here, since for n=1 there's no off-diagonal coupling to
// distinguish them.
template <typename SolverType>
void test_trivial_matrix_dimension_1()
{
    for (bool is_cyclic : {false, true}) {
        runCorrectnessCheck<SolverType>(/*n=*/1, /*batch_count=*/3, is_cyclic, /*offset=*/0, /*stride=*/1);
        runDiagonalCorrectnessCheck<SolverType>(/*n=*/1, /*batch_count=*/3, is_cyclic, /*offset=*/0, /*stride=*/1);
    }
}
INSTANTIATE_FOR_ALL_SOLVERS(trivial_matrix_dimension_1)

// --- matrix_dimension_ == 2: smallest nontrivial size - every PCR/CR step's iLeft/iRight clamps
// to the opposite boundary immediately (delta=1 already reaches both ends of a length-2 system).
template <typename SolverType>
void test_matrix_dimension_2()
{
    for (bool is_cyclic : {false, true}) {
        runCorrectnessCheck<SolverType>(/*n=*/2, /*batch_count=*/3, is_cyclic, /*offset=*/0, /*stride=*/1);
    }
}
INSTANTIATE_FOR_ALL_SOLVERS(matrix_dimension_2)

// --- Non-power-of-two matrix_dimension_: num_steps_ = ceil(log2(n)) is not exact, so the last
// reduction step operates on a system that's "overshot" past what n strictly requires. Covers odd,
// even, and values straddling a power-of-two boundary on both sides (3, 5, 6, 7 bracket 4 and 8).
template <typename SolverType>
void test_non_power_of_two_dimensions()
{
    for (int n : {3, 5, 6, 7}) {
        for (bool is_cyclic : {false, true}) {
            runCorrectnessCheck<SolverType>(n, /*batch_count=*/3, is_cyclic, /*offset=*/0, /*stride=*/1);
        }
    }
}
INSTANTIATE_FOR_ALL_SOLVERS(non_power_of_two_dimensions)

// --- batch_count_ == 1: TeamPolicy league_size == 1, and batch_count_ > matrix_dimension_: the
// inverse of the "few, long lines" regime the file's parallelization-model comments assume.
template <typename SolverType>
void test_batch_count_extremes()
{
    for (bool is_cyclic : {false, true}) {
        runCorrectnessCheck<SolverType>(/*n=*/5, /*batch_count=*/1, is_cyclic, /*offset=*/0, /*stride=*/1);
        runCorrectnessCheck<SolverType>(/*n=*/2, /*batch_count=*/10, is_cyclic, /*offset=*/0, /*stride=*/1);
    }
}
INSTANTIATE_FOR_ALL_SOLVERS(batch_count_extremes)

// --- matrix_dimension_ large enough that Kokkos::AUTO's chosen team_size is plausibly < n on at
// least some backends, forcing strided team loops to actually run more than one iteration per
// thread - the hand-derived tests all use n=4, small enough that most backends would just grant
// team_size >= n and never exercise that loop's stride > 0 case at all.
template <typename SolverType>
void test_large_dimension_strided_team_loop()
{
    for (bool is_cyclic : {false, true}) {
        runCorrectnessCheck<SolverType>(/*n=*/33, /*batch_count=*/2, is_cyclic, /*offset=*/0, /*stride=*/1);
    }
    runDiagonalCorrectnessCheck<SolverType>(/*n=*/33, /*batch_count=*/2, /*is_cyclic=*/true, /*offset=*/0,
                                             /*stride=*/1);
}
INSTANTIATE_FOR_ALL_SOLVERS(large_dimension_strided_team_loop)

// --- Uneven batch_offset_/batch_stride_ split: the hand-derived tests split batch_count_ evenly
// in half (offset 0/1, stride 2, batch_count_ divisible by 2). Here batch_count_ = 5 is not
// divisible by stride = 3, so effective_batch_count_'s ceiling-division arithmetic produces a
// genuinely uneven split (2, 2, 1) across three solve() calls sharing ONE setup() - also checks
// that calling solve() more than twice against the same stored trajectory still works.
template <typename SolverType>
void test_uneven_batch_stride_split()
{
    const int matrix_dimension = 4;
    const int batch_count      = 5;

    for (bool is_cyclic : {false, true}) {
        SolverType solver(matrix_dimension, batch_count, is_cyclic);

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
INSTANTIATE_FOR_ALL_SOLVERS(uneven_batch_stride_split)

// --- Default batch_offset_/batch_stride_ arguments: checks solve()/solve_diagonal() called with
// no arguments at all (full batch in a single call) still compile and behave as offset=0, stride=1.
template <typename SolverType>
void test_default_batch_arguments()
{
    const int matrix_dimension = 4;
    const int batch_count      = 3;

    for (bool is_cyclic : {false, true}) {
        SolverType solver(matrix_dimension, batch_count, is_cyclic);

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
INSTANTIATE_FOR_ALL_SOLVERS(default_batch_arguments)

// --- Large matrix_dimension_ (still within reach of the O(n^3) dense reference): exercises the
// reduction recursion over many more steps than the hand-derived tests, and stresses the strided
// team loops with a dimension unlikely to fit inside a single team's thread count on any backend.
// Covers both a non-power-of-two size (513, just past 2^9) and an exact power-of-two size (1024).
template <typename SolverType>
void test_very_large_matrix_dimension_non_cyclic()
{
    runCorrectnessCheck<SolverType>(/*n=*/513, /*batch_count=*/2, /*is_cyclic=*/false, /*offset=*/0, /*stride=*/1);
    runCorrectnessCheck<SolverType>(/*n=*/1024, /*batch_count=*/2, /*is_cyclic=*/false, /*offset=*/0, /*stride=*/1);
}
INSTANTIATE_FOR_ALL_SOLVERS(very_large_matrix_dimension_non_cyclic)

template <typename SolverType>
void test_very_large_matrix_dimension_cyclic()
{
    runCorrectnessCheck<SolverType>(/*n=*/513, /*batch_count=*/2, /*is_cyclic=*/true, /*offset=*/0, /*stride=*/1);
    runCorrectnessCheck<SolverType>(/*n=*/1024, /*batch_count=*/2, /*is_cyclic=*/true, /*offset=*/0, /*stride=*/1);
}
INSTANTIATE_FOR_ALL_SOLVERS(very_large_matrix_dimension_cyclic)

// --- Large matrix_dimension_ combined with a non-trivial batch_offset_/batch_stride_ split:
// checks that the batch-selection arithmetic and the large-n recursion compose correctly together,
// rather than each only being tested in isolation. batch_count_ = 5 with stride = 2 covers batches
// {0, 2, 4} on one call and {1, 3} on the other, sharing a single setup() at large n.
template <typename SolverType>
void test_very_large_matrix_dimension_with_batch_stride_non_cyclic()
{
    const int matrix_dimension = 777;
    const int batch_count      = 5;
    const bool is_cyclic       = false;

    SolverType solver(matrix_dimension, batch_count, is_cyclic);

    Kokkos::parallel_for(
        "Fill", 1, KOKKOS_LAMBDA(const int) {
            for (int b = 0; b < batch_count; b++) {
                for (int i = 0; i < matrix_dimension; i++) {
                    solver.set_main_diagonal(b, i, sysDiag(b, i));
                }
                for (int i = 0; i < matrix_dimension - 1; i++) {
                    solver.set_sub_diagonal(b, i, sysSub(b, i));
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

    const int stride = 2;
    solver.solve(rhs, /*offset=*/0, stride); // covers batch 0, 2, 4
    solver.solve(rhs, /*offset=*/1, stride); // covers batch 1, 3

    Kokkos::deep_copy(h_rhs, rhs);

    const double tol = 1e-8;
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
INSTANTIATE_FOR_ALL_SOLVERS(very_large_matrix_dimension_with_batch_stride_non_cyclic)

template <typename SolverType>
void test_very_large_matrix_dimension_with_batch_stride_cyclic()
{
    const int matrix_dimension = 777;
    const int batch_count      = 5;
    const bool is_cyclic       = true;

    SolverType solver(matrix_dimension, batch_count, is_cyclic);

    Kokkos::parallel_for(
        "Fill", 1, KOKKOS_LAMBDA(const int) {
            for (int b = 0; b < batch_count; b++) {
                for (int i = 0; i < matrix_dimension; i++) {
                    solver.set_main_diagonal(b, i, sysDiag(b, i));
                }
                for (int i = 0; i < matrix_dimension - 1; i++) {
                    solver.set_sub_diagonal(b, i, sysSub(b, i));
                }
                solver.set_cyclic_corner(b, sysCorner(b));
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

    const int stride = 2;
    solver.solve(rhs, /*offset=*/0, stride); // covers batch 0, 2, 4
    solver.solve(rhs, /*offset=*/1, stride); // covers batch 1, 3

    Kokkos::deep_copy(h_rhs, rhs);

    const double tol = 1e-8;
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
INSTANTIATE_FOR_ALL_SOLVERS(very_large_matrix_dimension_with_batch_stride_cyclic)
