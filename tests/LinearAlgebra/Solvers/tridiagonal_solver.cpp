#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <vector>
#include <cmath>

#include "../../../include/LinearAlgebra/Solvers/tridiagonal_solver.h"
#include "../../../include/LinearAlgebra/vector.h"

// Helper function to initialize tridiagonal system
template <typename T>
void initialize_tridiagonal_system(BatchedTridiagonalSolver<T>& solver, Vector<T>& rhs, int batch_idx, int matrix_dim,
                                   bool is_cyclic)
{
    // Create a simple test system: A*x = b
    // Main diagonal = 4, Sub diagonal = -1
    for (int i = 0; i < matrix_dim; i++) {
        solver.main_diagonal(batch_idx, i) = 4.0;
        if (i < matrix_dim - 1) {
            solver.sub_diagonal(batch_idx, i) = -1.0;
        }
    }

    if (is_cyclic) {
        // For cyclic, set the corner element
        solver.sub_diagonal(batch_idx, matrix_dim - 1) = -1.0;
    }

    // Set RHS to 1.0 for simple test
    for (int i = 0; i < matrix_dim; i++) {
        rhs(batch_idx * matrix_dim + i) = 1.0;
    }
}

// Helper function to verify solution
// NOTE: This function uses the ORIGINAL matrix values (before setup() modifies them)
// The test matrix has main_diagonal = 4.0 and sub_diagonal = -1.0
template <typename T>
bool verify_solution(const Vector<T>& x, int batch_idx, int matrix_dim, bool is_cyclic, T tolerance = 1e-6)
{
    // Reconstruct the original test matrix (main=4, sub=-1) and compute A*x
    // NOTE: We use hardcoded values because setup() has already modified the solver's internal state
    std::vector<T> Ax(matrix_dim, T(0));

    for (int i = 0; i < matrix_dim; ++i) {
        Ax[i] += T(4.0) * x(batch_idx * matrix_dim + i);
        if (i > 0) {
            Ax[i] += T(-1.0) * x(batch_idx * matrix_dim + i - 1);
        }
        if (i < matrix_dim - 1) {
            Ax[i] += T(-1.0) * x(batch_idx * matrix_dim + i + 1);
        }
    }

    if (is_cyclic) {
        Ax[0] += T(-1.0) * x(batch_idx * matrix_dim + matrix_dim - 1);
        Ax[matrix_dim - 1] += T(-1.0) * x(batch_idx * matrix_dim + 0);
    }

    // Check if Ax ≈ b (original RHS = 1.0)
    T max_error = T(0);
    for (int i = 0; i < matrix_dim; ++i) {
        T error = std::abs(Ax[i] - T(1.0));
        if (error > max_error)
            max_error = error;
    }

    return max_error < tolerance;
}

// Helper function to compute actual error for more detailed testing
// NOTE: Uses hardcoded original matrix values (main=4.0, sub=-1.0)
template <typename T>
T compute_solution_error(const Vector<T>& x, int batch_idx, int matrix_dim, bool is_cyclic)
{
    std::vector<T> Ax(matrix_dim, T(0));

    for (int i = 0; i < matrix_dim; ++i) {
        Ax[i] += T(4.0) * x(batch_idx * matrix_dim + i);
        if (i > 0) {
            Ax[i] += T(-1.0) * x(batch_idx * matrix_dim + i - 1);
        }
        if (i < matrix_dim - 1) {
            Ax[i] += T(-1.0) * x(batch_idx * matrix_dim + i + 1);
        }
    }

    if (is_cyclic) {
        Ax[0] += T(-1.0) * x(batch_idx * matrix_dim + matrix_dim - 1);
        Ax[matrix_dim - 1] += T(-1.0) * x(batch_idx * matrix_dim + 0);
    }

    T max_error = T(0);
    for (int i = 0; i < matrix_dim; ++i) {
        T error = std::abs(Ax[i] - T(1.0));
        if (error > max_error)
            max_error = error;
    }

    return max_error;
}

// Test fixture for BatchedTridiagonalSolver tests
class BatchedTridiagonalSolverTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        // Kokkos initialization is handled in main()
    }

    void TearDown() override
    {
        // Cleanup if needed
    }

    // Common test parameters
    static constexpr int default_matrix_dim   = 10;
    static constexpr int default_batch_count  = 8;
    static constexpr double default_tolerance = 1e-6;
};

// ============================================================================
// Basic Functionality Tests
// ============================================================================

TEST_F(BatchedTridiagonalSolverTest, NonCyclicAllBatches)
{
    const int matrix_dim  = default_matrix_dim;
    const int batch_count = default_batch_count;

    BatchedTridiagonalSolver<double> solver(matrix_dim, batch_count, false);
    Vector<double> rhs("rhs", matrix_dim * batch_count);

    // Initialize all batches
    for (int b = 0; b < batch_count; b++) {
        initialize_tridiagonal_system(solver, rhs, b, matrix_dim, false);
    }

    // Factorize
    solver.setup();

    // Solve all batches
    solver.solve(rhs, 0, 1);

    // Verify all batches
    for (int b = 0; b < batch_count; b++) {
        EXPECT_TRUE(verify_solution(rhs, b, matrix_dim, false, default_tolerance))
            << "Batch " << b << " failed verification";
    }
}

TEST_F(BatchedTridiagonalSolverTest, CyclicAllBatches)
{
    const int matrix_dim  = default_matrix_dim;
    const int batch_count = default_batch_count;

    BatchedTridiagonalSolver<double> solver(matrix_dim, batch_count, true);
    Vector<double> rhs("rhs", matrix_dim * batch_count);

    // Initialize all batches
    for (int b = 0; b < batch_count; b++) {
        initialize_tridiagonal_system(solver, rhs, b, matrix_dim, true);
    }

    // Factorize
    solver.setup();

    // Solve all batches
    solver.solve(rhs, 0, 1);

    // Verify all batches
    for (int b = 0; b < batch_count; b++) {
        EXPECT_TRUE(verify_solution(rhs, b, matrix_dim, true, default_tolerance))
            << "Batch " << b << " failed verification";
    }
}

// ============================================================================
// Stride and Offset Tests
// ============================================================================

TEST_F(BatchedTridiagonalSolverTest, NonCyclicEvenBatchesStride2Offset0)
{
    const int matrix_dim  = default_matrix_dim;
    const int batch_count = default_batch_count;

    BatchedTridiagonalSolver<double> solver(matrix_dim, batch_count, false);
    Vector<double> rhs("rhs", matrix_dim * batch_count);

    // Initialize all batches
    for (int b = 0; b < batch_count; b++) {
        initialize_tridiagonal_system(solver, rhs, b, matrix_dim, false);
    }

    solver.setup();

    // Solve only even batches (0, 2, 4, 6)
    solver.solve(rhs, 0, 2);

    // Verify only even batches
    std::vector<int> even_batches = {0, 2, 4, 6};
    for (int b : even_batches) {
        EXPECT_TRUE(verify_solution(rhs, b, matrix_dim, false, default_tolerance))
            << "Even batch " << b << " failed verification";
    }
}

TEST_F(BatchedTridiagonalSolverTest, NonCyclicOddBatchesStride2Offset1)
{
    const int matrix_dim  = default_matrix_dim;
    const int batch_count = default_batch_count;

    BatchedTridiagonalSolver<double> solver(matrix_dim, batch_count, false);
    Vector<double> rhs("rhs", matrix_dim * batch_count);

    // Initialize all batches
    for (int b = 0; b < batch_count; b++) {
        initialize_tridiagonal_system(solver, rhs, b, matrix_dim, false);
    }

    solver.setup();

    // Solve only odd batches (1, 3, 5, 7)
    solver.solve(rhs, 1, 2);

    // Verify only odd batches
    std::vector<int> odd_batches = {1, 3, 5, 7};
    for (int b : odd_batches) {
        EXPECT_TRUE(verify_solution(rhs, b, matrix_dim, false, default_tolerance))
            << "Odd batch " << b << " failed verification";
    }
}

TEST_F(BatchedTridiagonalSolverTest, CyclicStride3Offset1)
{
    const int matrix_dim  = default_matrix_dim;
    const int batch_count = default_batch_count;

    BatchedTridiagonalSolver<double> solver(matrix_dim, batch_count, true);
    Vector<double> rhs("rhs", matrix_dim * batch_count);

    // Initialize all batches
    for (int b = 0; b < batch_count; b++) {
        initialize_tridiagonal_system(solver, rhs, b, matrix_dim, true);
    }

    solver.setup();

    // Solve batches with stride 3, offset 1 (1, 4, 7)
    solver.solve(rhs, 1, 3);

    // Verify
    std::vector<int> batch_indices = {1, 4, 7};
    for (int b : batch_indices) {
        EXPECT_TRUE(verify_solution(rhs, b, matrix_dim, true, default_tolerance))
            << "Batch " << b << " failed verification";
    }
}

// ============================================================================
// Edge Cases
// ============================================================================

TEST_F(BatchedTridiagonalSolverTest, SingleBatchWithOffset)
{
    const int matrix_dim   = default_matrix_dim;
    const int batch_count  = default_batch_count;
    const int target_batch = 5;

    BatchedTridiagonalSolver<double> solver(matrix_dim, batch_count, false);
    Vector<double> rhs("rhs", matrix_dim * batch_count);

    // Initialize all batches
    for (int b = 0; b < batch_count; b++) {
        initialize_tridiagonal_system(solver, rhs, b, matrix_dim, false);
    }

    solver.setup();

    // Solve only batch 5
    solver.solve(rhs, target_batch, batch_count);

    // Verify only batch 5
    EXPECT_TRUE(verify_solution(rhs, target_batch, matrix_dim, false, default_tolerance))
        << "Single batch " << target_batch << " failed verification";
}

TEST_F(BatchedTridiagonalSolverTest, SmallMatrixSize)
{
    const int matrix_dim  = 3;
    const int batch_count = 4;

    BatchedTridiagonalSolver<double> solver(matrix_dim, batch_count, false);
    Vector<double> rhs("rhs", matrix_dim * batch_count);

    for (int b = 0; b < batch_count; b++) {
        initialize_tridiagonal_system(solver, rhs, b, matrix_dim, false);
    }

    solver.setup();
    solver.solve(rhs, 0, 1);

    for (int b = 0; b < batch_count; b++) {
        EXPECT_TRUE(verify_solution(rhs, b, matrix_dim, false, default_tolerance))
            << "Small matrix batch " << b << " failed verification";
    }
}

TEST_F(BatchedTridiagonalSolverTest, SingleBatchSingleElement)
{
    const int matrix_dim  = 1;
    const int batch_count = 1;

    BatchedTridiagonalSolver<double> solver(matrix_dim, batch_count, false);
    Vector<double> rhs("rhs", matrix_dim * batch_count);

    // For 1x1 matrix, just set main diagonal
    solver.main_diagonal(0, 0) = 4.0;
    rhs(0)                     = 1.0;

    solver.setup();
    solver.solve(rhs, 0, 1);

    // Solution should be 1.0 / 4.0 = 0.25
    EXPECT_NEAR(rhs(0), 0.25, default_tolerance);
}

// ============================================================================
// Error Handling Tests
// ============================================================================

TEST_F(BatchedTridiagonalSolverTest, SolveBeforeSetupThrows)
{
    const int matrix_dim  = default_matrix_dim;
    const int batch_count = default_batch_count;

    BatchedTridiagonalSolver<double> solver(matrix_dim, batch_count, false);
    Vector<double> rhs("rhs", matrix_dim * batch_count);

    // Try to solve without calling setup() first
    EXPECT_THROW({ solver.solve(rhs, 0, 1); }, std::runtime_error);
}

TEST_F(BatchedTridiagonalSolverTest, DiagonalSolveBeforeSetupThrows)
{
    const int matrix_dim  = default_matrix_dim;
    const int batch_count = default_batch_count;

    BatchedTridiagonalSolver<double> solver(matrix_dim, batch_count, false);
    Vector<double> rhs("rhs", matrix_dim * batch_count);

    // Try to solve_diagonal without calling setup() first
    EXPECT_THROW({ solver.solve_diagonal(rhs, 0, 1); }, std::runtime_error);
}

// ============================================================================
// Diagonal Solve Tests
// ============================================================================

TEST_F(BatchedTridiagonalSolverTest, DiagonalSolveAllBatches)
{
    const int matrix_dim  = default_matrix_dim;
    const int batch_count = default_batch_count;

    BatchedTridiagonalSolver<double> solver(matrix_dim, batch_count, false);
    Vector<double> rhs("rhs", matrix_dim * batch_count);

    // Initialize with diagonal-only system
    for (int b = 0; b < batch_count; b++) {
        for (int i = 0; i < matrix_dim; i++) {
            solver.main_diagonal(b, i) = 2.0;
            rhs(b * matrix_dim + i)    = 4.0;
        }
    }

    solver.setup();
    solver.solve_diagonal(rhs, 0, 1);

    // Each solution element should be 4.0 / 2.0 = 2.0
    for (int b = 0; b < batch_count; b++) {
        for (int i = 0; i < matrix_dim; i++) {
            EXPECT_NEAR(rhs(b * matrix_dim + i), 2.0, default_tolerance)
                << "Diagonal solve failed at batch " << b << ", index " << i;
        }
    }
}

TEST_F(BatchedTridiagonalSolverTest, DiagonalSolveWithStride)
{
    const int matrix_dim  = default_matrix_dim;
    const int batch_count = default_batch_count;

    BatchedTridiagonalSolver<double> solver(matrix_dim, batch_count, false);
    Vector<double> rhs("rhs", matrix_dim * batch_count);

    // Initialize with diagonal-only system
    for (int b = 0; b < batch_count; b++) {
        for (int i = 0; i < matrix_dim; i++) {
            solver.main_diagonal(b, i) = 3.0;
            rhs(b * matrix_dim + i)    = 6.0;
        }
    }

    solver.setup();
    solver.solve_diagonal(rhs, 1, 2); // Solve odd batches only

    // Check odd batches: should be 6.0 / 3.0 = 2.0
    std::vector<int> odd_batches = {1, 3, 5, 7};
    for (int b : odd_batches) {
        for (int i = 0; i < matrix_dim; i++) {
            EXPECT_NEAR(rhs(b * matrix_dim + i), 2.0, default_tolerance)
                << "Diagonal solve with stride failed at batch " << b << ", index " << i;
        }
    }
}

// ============================================================================
// Numerical Accuracy Tests
// ============================================================================

TEST_F(BatchedTridiagonalSolverTest, AccuracyNonCyclic)
{
    const int matrix_dim  = default_matrix_dim;
    const int batch_count = default_batch_count;

    BatchedTridiagonalSolver<double> solver(matrix_dim, batch_count, false);
    Vector<double> rhs("rhs", matrix_dim * batch_count);

    for (int b = 0; b < batch_count; b++) {
        initialize_tridiagonal_system(solver, rhs, b, matrix_dim, false);
    }

    solver.setup();
    solver.solve(rhs, 0, 1);

    // Check that error is well below tolerance
    for (int b = 0; b < batch_count; b++) {
        double error = compute_solution_error(rhs, b, matrix_dim, false);
        EXPECT_LT(error, 1e-10) << "Error too large for batch " << b;
    }
}

TEST_F(BatchedTridiagonalSolverTest, AccuracyCyclic)
{
    const int matrix_dim  = default_matrix_dim;
    const int batch_count = default_batch_count;

    BatchedTridiagonalSolver<double> solver(matrix_dim, batch_count, true);
    Vector<double> rhs("rhs", matrix_dim * batch_count);

    for (int b = 0; b < batch_count; b++) {
        initialize_tridiagonal_system(solver, rhs, b, matrix_dim, true);
    }

    solver.setup();
    solver.solve(rhs, 0, 1);

    // Check that error is well below tolerance
    for (int b = 0; b < batch_count; b++) {
        double error = compute_solution_error(rhs, b, matrix_dim, true);
        EXPECT_LT(error, 1e-10) << "Error too large for cyclic batch " << b;
    }
}

// ============================================================================
// Performance/Stress Tests
// ============================================================================

TEST_F(BatchedTridiagonalSolverTest, LargeBatchCount)
{
    const int matrix_dim  = default_matrix_dim;
    const int large_batch = 1000;

    BatchedTridiagonalSolver<double> solver(matrix_dim, large_batch, false);
    Vector<double> rhs("rhs", matrix_dim * large_batch);

    // Initialize
    for (int b = 0; b < large_batch; b++) {
        initialize_tridiagonal_system(solver, rhs, b, matrix_dim, false);
    }

    // Factorize and solve
    solver.setup();
    solver.solve(rhs, 0, 1);

    // Verify a few random batches
    std::vector<int> test_batches = {0, 250, 500, 750, 999};
    for (int b : test_batches) {
        EXPECT_TRUE(verify_solution(rhs, b, matrix_dim, false, default_tolerance))
            << "Large batch test failed at batch " << b;
    }
}

TEST_F(BatchedTridiagonalSolverTest, LargeMatrixDimension)
{
    const int matrix_dim  = 100;
    const int batch_count = 4;

    BatchedTridiagonalSolver<double> solver(matrix_dim, batch_count, false);
    Vector<double> rhs("rhs", matrix_dim * batch_count);

    for (int b = 0; b < batch_count; b++) {
        initialize_tridiagonal_system(solver, rhs, b, matrix_dim, false);
    }

    solver.setup();
    solver.solve(rhs, 0, 1);

    for (int b = 0; b < batch_count; b++) {
        EXPECT_TRUE(verify_solution(rhs, b, matrix_dim, false, default_tolerance))
            << "Large matrix dimension test failed at batch " << b;
    }
}

// ============================================================================
// Accessor Tests
// ============================================================================

TEST_F(BatchedTridiagonalSolverTest, MainDiagonalAccessors)
{
    const int matrix_dim  = 5;
    const int batch_count = 2;

    BatchedTridiagonalSolver<double> solver(matrix_dim, batch_count, false);

    // Test write access
    for (int b = 0; b < batch_count; b++) {
        for (int i = 0; i < matrix_dim; i++) {
            solver.main_diagonal(b, i) = static_cast<double>(b * matrix_dim + i);
        }
    }

    // Test read access BEFORE setup()
    for (int b = 0; b < batch_count; b++) {
        for (int i = 0; i < matrix_dim; i++) {
            EXPECT_DOUBLE_EQ(solver.main_diagonal(b, i), static_cast<double>(b * matrix_dim + i))
                << "Main diagonal accessor failed at batch " << b << ", index " << i;
        }
    }
}

TEST_F(BatchedTridiagonalSolverTest, SubDiagonalAccessors)
{
    const int matrix_dim  = 5;
    const int batch_count = 2;

    BatchedTridiagonalSolver<double> solver(matrix_dim, batch_count, false);

    // Test write access
    for (int b = 0; b < batch_count; b++) {
        for (int i = 0; i < matrix_dim; i++) {
            solver.sub_diagonal(b, i) = static_cast<double>(100 + b * matrix_dim + i);
        }
    }

    // Test read access BEFORE setup()
    for (int b = 0; b < batch_count; b++) {
        for (int i = 0; i < matrix_dim; i++) {
            EXPECT_DOUBLE_EQ(solver.sub_diagonal(b, i), static_cast<double>(100 + b * matrix_dim + i))
                << "Sub diagonal accessor failed at batch " << b << ", index " << i;
        }
    }
}

TEST_F(BatchedTridiagonalSolverTest, CyclicCornerAccessors)
{
    const int matrix_dim  = 5;
    const int batch_count = 2;

    BatchedTridiagonalSolver<double> solver(matrix_dim, batch_count, true);

    // Test write access via cyclic_corner
    for (int b = 0; b < batch_count; b++) {
        solver.cyclic_corner(b) = static_cast<double>(200 + b);
    }

    // Test read access via cyclic_corner BEFORE setup()
    for (int b = 0; b < batch_count; b++) {
        EXPECT_DOUBLE_EQ(solver.cyclic_corner(b), static_cast<double>(200 + b))
            << "Cyclic corner accessor failed at batch " << b;
    }

    // Verify that cyclic_corner actually accesses sub_diagonal at the right location
    for (int b = 0; b < batch_count; b++) {
        EXPECT_DOUBLE_EQ(solver.cyclic_corner(b), solver.sub_diagonal(b, matrix_dim - 1))
            << "Cyclic corner should access sub_diagonal at index matrix_dim-1";
    }
}

TEST_F(BatchedTridiagonalSolverTest, SetupModifiesInternalState)
{
    // This test verifies that setup() modifies the internal diagonal values
    const int matrix_dim  = 5;
    const int batch_count = 2;

    BatchedTridiagonalSolver<double> solver(matrix_dim, batch_count, false);

    // Set up initial values
    for (int b = 0; b < batch_count; b++) {
        for (int i = 0; i < matrix_dim; i++) {
            solver.main_diagonal(b, i) = 4.0;
            if (i < matrix_dim - 1) {
                solver.sub_diagonal(b, i) = -1.0;
            }
        }
    }

    // Store original values for comparison
    std::vector<double> original_main(matrix_dim * batch_count);
    std::vector<double> original_sub(matrix_dim * batch_count);
    for (int b = 0; b < batch_count; b++) {
        for (int i = 0; i < matrix_dim; i++) {
            original_main[b * matrix_dim + i] = solver.main_diagonal(b, i);
            original_sub[b * matrix_dim + i]  = solver.sub_diagonal(b, i);
        }
    }

    // Call setup() - this performs Cholesky factorization
    solver.setup();

    // Verify that values have been modified
    bool main_changed = false;
    bool sub_changed  = false;
    for (int b = 0; b < batch_count; b++) {
        for (int i = 0; i < matrix_dim; i++) {
            if (std::abs(solver.main_diagonal(b, i) - original_main[b * matrix_dim + i]) > 1e-10) {
                main_changed = true;
            }
            if (std::abs(solver.sub_diagonal(b, i) - original_sub[b * matrix_dim + i]) > 1e-10) {
                sub_changed = true;
            }
        }
    }

    EXPECT_TRUE(main_changed) << "setup() should modify main diagonal values (Cholesky factorization)";
    EXPECT_TRUE(sub_changed) << "setup() should modify sub diagonal values (Cholesky factorization)";
}
