#ifdef GMGPOLAR_USE_MUMPS

    #include <gtest/gtest.h>

    #include <LinearAlgebra/Vector/vector.h>
    #include <LinearAlgebra/Matrix/coo_matrix.h>
    #include <LinearAlgebra/Solvers/coo_mumps_solver.h>
using namespace gmgpolar;

void fill_rhs(Vector<double>& rhs)
{
    Kokkos::parallel_for("fill rhs", rhs.size(), KOKKOS_LAMBDA(int i) { rhs(i) = (i + 1) * 2.0; });
}

// -----------------------------------------------------------------------
// Test 1: General (non-symmetric) 4x4 system
//
// Matrix A (non-symmetric):
//   [ 1  0  2  0 ]
//   [ 3  0  4  5 ]
//   [ 0  6  7  0 ]
//   [ 0  8  0  9 ]
//
// RHS b = [2, 4, 6, 8]^T
// -----------------------------------------------------------------------
TEST(CooMumpsSolverTest, GeneralNonSymmetric4x4)
{
    using triplet = SparseMatrixCOO<double>::triplet_type;

    // All non-zero entries (0-based row/col indices)
    std::vector<triplet> entries = {{0, 0, 1.0}, {0, 2, 2.0}, {1, 0, 3.0}, {1, 2, 4.0}, {1, 3, 5.0},
                                    {2, 1, 6.0}, {2, 2, 7.0}, {3, 1, 8.0}, {3, 3, 9.0}};

    SparseMatrixCOO<double> mat(4, 4, entries);
    mat.is_symmetric(false);

    CooMumpsSolver solver(std::move(mat));

    Vector<double> rhs("rhs", 4);
    fill_rhs(rhs);

    solver.solveInPlace(rhs);

    HostVector<double> solution("solution", 4);
    solution(0) = 140.0 / 43.0;
    solution(1) = 149.0 / 86.0;
    solution(2) = -27.0 / 43.0;
    solution(3) = -28.0 / 43.0;

    auto rhs_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), rhs);

    const double tol = 1e-10;
    EXPECT_NEAR(rhs_host(0), solution(0), tol);
    EXPECT_NEAR(rhs_host(1), solution(1), tol);
    EXPECT_NEAR(rhs_host(2), solution(2), tol);
    EXPECT_NEAR(rhs_host(3), solution(3), tol);
}

// -----------------------------------------------------------------------
// Test 2: Symmetric positive-definite 4x4 system
//
// Matrix A (SPD):
//   [ 4  0  2  0 ]
//   [ 0  5  1  3 ]
//   [ 2  1  6  2 ]
//   [ 0  3  2  7 ]
//
// RHS b = [2, 4, 6, 8]^T
// -----------------------------------------------------------------------
TEST(CooMumpsSolverTest, SymmetricPositiveDefinite4x4)
{
    using triplet = SparseMatrixCOO<double>::triplet_type;

    std::vector<triplet> entries = {{0, 0, 4.0}, {1, 1, 5.0}, {2, 0, 2.0}, {2, 1, 1.0}, {2, 2, 6.0}, {3, 1, 3.0},
                                    {3, 2, 2.0}, {3, 3, 7.0}, {0, 2, 2.0}, {1, 2, 1.0}, {1, 3, 3.0}, {2, 3, 2.0}};

    SparseMatrixCOO<double> mat(4, 4, entries);
    mat.is_symmetric(true);

    CooMumpsSolver solver(std::move(mat));

    Vector<double> rhs("rhs", 4);
    fill_rhs(rhs);

    solver.solveInPlace(rhs);

    HostVector<double> solution("solution", 4);
    solution(0) = 9.0 / 46.0;
    solution(1) = 3.0 / 23.0;
    solution(2) = 14.0 / 23.0;
    solution(3) = 21.0 / 23.0;

    auto rhs_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), rhs);

    const double tol = 1e-10;
    EXPECT_NEAR(rhs_host(0), solution(0), tol);
    EXPECT_NEAR(rhs_host(1), solution(1), tol);
    EXPECT_NEAR(rhs_host(2), solution(2), tol);
    EXPECT_NEAR(rhs_host(3), solution(3), tol);
}

#endif // GMGPOLAR_USE_MUMPS
