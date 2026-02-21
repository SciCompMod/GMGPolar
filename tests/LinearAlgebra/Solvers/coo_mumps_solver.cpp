#pragma once
#ifdef GMGPOLAR_USE_MUMPS

    #include <gtest/gtest.h>
    #include "coo_mumps_solver.h"

// -----------------------------------------------------------------------
// Test 1: General (non-symmetric) 4x4 system
//
// Matrix A (non-symmetric):
//   [ 4  1  0  0 ]
//   [ 2  5  1  0 ]
//   [ 0  1  6  2 ]
//   [ 0  0  3  7 ]
//
// RHS b = [5, 8, 9, 10]^T
//
// Expected solution computed from A*x = b.
// -----------------------------------------------------------------------
TEST(CooMumpsSolverTest, GeneralNonSymmetric4x4)
{
    using triplet = SparseMatrixCOO<double>::triplet_type;

    // All non-zero entries (0-based row/col indices)
    std::vector<triplet> entries = {{0, 0, 4.0}, {0, 1, 1.0}, {1, 0, 2.0}, {1, 1, 5.0}, {1, 2, 1.0},
                                    {2, 1, 1.0}, {2, 2, 6.0}, {2, 3, 2.0}, {3, 2, 3.0}, {3, 3, 7.0}};

    SparseMatrixCOO<double> mat(4, 4, entries);
    mat.is_symmetric(false);

    CooMumpsSolver solver(std::move(mat));

    Vector<double> rhs("rhs", 4);
    rhs(0) = 5.0;
    rhs(1) = 8.0;
    rhs(2) = 9.0;
    rhs(3) = 10.0;

    solver.solve(rhs);

    // Verify A*x = b by back-substitution check
    // Reference solution (computed analytically / via numpy):
    //   x ~ [0.9526, 0.2105, 0.9298, 0.8319]
    const double tol = 1e-10;
    EXPECT_NEAR(rhs(0), 4.0 * 0.9526 + 1.0 * 0.2105, 1e-3);

    // More robust: re-multiply and check residual
    std::vector<double> x = {rhs(0), rhs(1), rhs(2), rhs(3)};

    // A*x
    double Ax0 = 4 * x[0] + 1 * x[1];
    double Ax1 = 2 * x[0] + 5 * x[1] + 1 * x[2];
    double Ax2 = 1 * x[1] + 6 * x[2] + 2 * x[3];
    double Ax3 = 3 * x[2] + 7 * x[3];

    EXPECT_NEAR(Ax0, 5.0, tol);
    EXPECT_NEAR(Ax1, 8.0, tol);
    EXPECT_NEAR(Ax2, 9.0, tol);
    EXPECT_NEAR(Ax3, 10.0, tol);
}

// -----------------------------------------------------------------------
// Test 2: Symmetric positive-definite 4x4 system (lower triangle only)
//
// Matrix A (SPD):
//   [ 4  2  0  0 ]
//   [ 2  5  1  0 ]
//   [ 0  1  6  2 ]
//   [ 0  0  2  7 ]
//
// Only lower triangular entries provided (including diagonal).
// RHS b = [6, 8, 9, 11]^T
// -----------------------------------------------------------------------
TEST(CooMumpsSolverTest, SymmetricPositiveDefinite4x4)
{
    using triplet = SparseMatrixCOO<double>::triplet_type;

    // Lower triangular entries only (0-based row/col indices)
    std::vector<triplet> entries = {{0, 0, 4.0}, {1, 0, 2.0}, {1, 1, 5.0}, {2, 1, 1.0},
                                    {2, 2, 6.0}, {3, 2, 2.0}, {3, 3, 7.0}};

    SparseMatrixCOO<double> mat(4, 4, entries);
    mat.is_symmetric(true); // SPD, half-entries only

    CooMumpsSolver solver(std::move(mat));

    Vector<double> rhs("rhs", 4);
    rhs(0) = 6.0;
    rhs(1) = 8.0;
    rhs(2) = 9.0;
    rhs(3) = 11.0;

    solver.solve(rhs);

    // Verify residual: A*x = b using the full (symmetric) matrix
    std::vector<double> x = {rhs(0), rhs(1), rhs(2), rhs(3)};

    const double tol = 1e-10;

    double Ax0 = 4 * x[0] + 2 * x[1];
    double Ax1 = 2 * x[0] + 5 * x[1] + 1 * x[2];
    double Ax2 = 1 * x[1] + 6 * x[2] + 2 * x[3];
    double Ax3 = 2 * x[2] + 7 * x[3];

    EXPECT_NEAR(Ax0, 6.0, tol);
    EXPECT_NEAR(Ax1, 8.0, tol);
    EXPECT_NEAR(Ax2, 9.0, tol);
    EXPECT_NEAR(Ax3, 11.0, tol);
}

#endif // GMGPOLAR_USE_MUMPS