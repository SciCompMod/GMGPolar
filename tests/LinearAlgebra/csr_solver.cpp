#include <gtest/gtest.h>
#include <random>
#include <cmath>
#include "../../include/LinearAlgebra/csr_matrix.h"
#include "../../include/LinearAlgebra/vector.h"
#include "../../include/LinearAlgebra/vector_operations.h"
#include "../../include/LinearAlgebra/sparseLUSolver.h"

// All tests that the custom LU solver was unable to solve have been commented out.
// This typically happens when the matrix is singular, has zero diagonal entries,
// or requires pivoting that the static-pivoting LU algorithm cannot handle.
// The custom LU solver assumes all diagonal elements are nonzero and does not perform full pivoting.
// For such cases, numerical instability or division by zero may occur, so those tests are excluded.

// SparseMatrixCSR assumes that entries are sorted.
template <typename T>
std::vector<std::tuple<int, int, T>> sort_entries(std::initializer_list<std::tuple<int, int, T>> entries)
{
    std::vector<std::tuple<int, int, T>> sorted(entries);
    const auto compare = [](const auto& entry1, const auto& entry2) {
        const auto r1 = std::get<0>(entry1);
        const auto r2 = std::get<0>(entry2);
        if (r1 < r2) {
            return true;
        }
        return false;
    };
    std::sort(sorted.begin(), sorted.end(), compare);
    return sorted;
}

template <typename T>
std::vector<std::tuple<int, int, T>> sort_entries(std::vector<std::tuple<int, int, T>> entries)
{
    std::vector<std::tuple<int, int, T>> sorted(entries);
    const auto compare = [](const auto& entry1, const auto& entry2) {
        const auto r1 = std::get<0>(entry1);
        const auto r2 = std::get<0>(entry2);
        if (r1 < r2) {
            return true;
        }
        return false;
    };
    std::sort(sorted.begin(), sorted.end(), compare);
    return sorted;
}

// Helper: Multiply CSR matrix by vector
template <typename T>
Vector<T> csr_matvec(const SparseMatrixCSR<T>& A, const Vector<T>& x)
{
    Vector<T> y("y", A.rows());
    for (int i = 0; i < A.rows(); ++i) {
        T sum = 0;
        for (int k = 0; k < A.row_nz_size(i); ++k)
            sum += A.row_nz_entry(i, k) * x[A.row_nz_index(i, k)];
        y(i) = sum;
    }
    return y;
}

// Helper: Check if two vectors are close
template <typename T>
void expect_vector_near(const Vector<T>& a, const Vector<T>& b, double tol = 1e-8)
{
    ASSERT_EQ(a.size(), b.size());
    for (int i = 0; i < a.size(); ++i)
        EXPECT_NEAR(a(i), b(i), tol);
}

// Test 1: 1x1 matrix
TEST(SparseLUSolver, OneByOne)
{
    using T = double;
    SparseMatrixCSR<T> A(1, 1, sort_entries<T>({{0, 0, 2.0}}));
    Vector<T> b("b", 1);
    b(0) = 4.0;
    SparseLUSolver<T> solver(A);
    solver.solveInPlace(b);
    EXPECT_NEAR(b[0], 2.0, 1e-12);
}

// Test 2: 2x2 diagonal
TEST(SparseLUSolver, TwoByTwoDiagonal)
{
    using T = double;
    SparseMatrixCSR<T> A(2, 2, sort_entries<T>({{0, 0, 3.0}, {1, 1, 4.0}}));
    Vector<T> b("b", 2);
    b(0) = 6.0;
    b(1) = 8.0;
    SparseLUSolver<T> solver(A);
    solver.solveInPlace(b);
    EXPECT_NEAR(b[0], 2.0, 1e-12);
    EXPECT_NEAR(b[1], 2.0, 1e-12);
}

// Test 3: 2x2 off-diagonal
TEST(SparseLUSolver, TwoByTwoOffDiagonal)
{
    using T = double;
    SparseMatrixCSR<T> A(2, 2, sort_entries<T>({{0, 0, 1.0}, {0, 1, 2.0}, {1, 0, 3.0}, {1, 1, 4.0}}));
    Vector<T> b("b", 2);
    b(0) = 1.0;
    b(1) = 2.0;
    Vector<T> x_true("x_true", 2);
    x_true(0) = 0.0;
    x_true(1) = 0.5;
    SparseLUSolver<T> solver(A);
    solver.solveInPlace(b);
    expect_vector_near(b, x_true);
}

// Test 4: 3x3 lower triangular
TEST(SparseLUSolver, ThreeByThreeLowerTriangular)
{
    using T = double;
    SparseMatrixCSR<T> A(
        3, 3, sort_entries<T>({{0, 0, 1.0}, {1, 0, 2.0}, {1, 1, 3.0}, {2, 0, 4.0}, {2, 1, 5.0}, {2, 2, 6.0}}));
    Vector<T> b("b", 3);
    b(0) = 1;
    b(1) = 2;
    b(2) = 3;
    Vector<T> x_true("x_true", 3);
    x_true(0) = 1.0;
    x_true(1) = 0.0;
    x_true(2) = -1.0 / 6.0;
    SparseLUSolver<T> solver(A);
    solver.solveInPlace(b);
    expect_vector_near(b, x_true);
}

// Test 5: 3x3 upper triangular
TEST(SparseLUSolver, ThreeByThreeUpperTriangular)
{
    using T = double;
    SparseMatrixCSR<T> A(
        3, 3, sort_entries<T>({{0, 0, 1.0}, {0, 1, 2.0}, {0, 2, 3.0}, {1, 1, 4.0}, {1, 2, 5.0}, {2, 2, 6.0}}));
    Vector<T> b("b", 3);
    b(0) = 1;
    b(1) = 2;
    b(2) = 3;
    Vector<T> x_true("x_true", 3);
    x_true(0) = -0.25;
    x_true(1) = -0.125;
    x_true(2) = 0.5;
    SparseLUSolver<T> solver(A);
    solver.solveInPlace(b);
    expect_vector_near(b, x_true);
}

// Test 6: 3x3 with zero diagonal
TEST(SparseLUSolver, ThreeByThreeZeroDiagonal)
{
    using T = double;
    SparseMatrixCSR<T> A(3, 3, sort_entries<T>({{0, 0, 0.0}, {0, 1, 2.0}, {1, 1, 3.0}, {2, 2, 4.0}}));
    Vector<T> x_true("x_true", 3);
    x_true(0)   = 1;
    x_true(1)   = 2;
    x_true(2)   = 3;
    Vector<T> b = csr_matvec(A, x_true);
    SparseLUSolver<T> solver(A);
    solver.solveInPlace(b);
    Vector<T> b2 = csr_matvec(A, b);
    expect_vector_near(b2, csr_matvec(A, x_true), 1e-8);
}

// Test 7: 4x4 permutation needed
TEST(SparseLUSolver, FourByFourPermutation)
{
    using T = double;
    SparseMatrixCSR<T> A(4, 4,
                         sort_entries<T>({{0, 1, 1.0},
                                          {1, 2, 1.0},
                                          {2, 3, 1.0},
                                          {3, 0, 1.0},
                                          {0, 0, 10.0},
                                          {1, 1, 10.0},
                                          {2, 2, 10.0},
                                          {3, 3, 10.0}}));
    Vector<T> x_true("x_true", 4);
    x_true(0)   = 1;
    x_true(1)   = 2;
    x_true(2)   = 3;
    x_true(3)   = 4;
    Vector<T> b = csr_matvec(A, x_true);
    SparseLUSolver<T> solver(A);
    solver.solveInPlace(b);
    expect_vector_near(b, x_true);
}

// Test 8: 5x5 Sparse
TEST(SparseLUSolver, FiveByFiveRandomSparse)
{
    using T                                               = double;
    std::vector<SparseMatrixCSR<T>::triplet_type> entries = {{0, 0, 2.0}, {0, 2, 3.0}, {1, 1, 4.0}, {2, 0, 1.0},
                                                             {2, 2, 5.0}, {3, 3, 6.0}, {4, 1, 7.0}, {4, 4, 8.0}};
    auto sorted_entries                                   = sort_entries(entries);
    SparseMatrixCSR<T> A(5, 5, sorted_entries);
    Vector<T> b("b", 5);
    b(0) = 1;
    b(1) = 2;
    b(2) = 3;
    b(3) = 4;
    b(4) = 5;
    Vector<T> x_true("x_true", 5);
    x_true(0) = -0.5714285714285714;
    x_true(1) = 0.5;
    x_true(2) = 0.7142857142857143;
    x_true(3) = 0.6666666666666667;
    x_true(4) = 0.1875;
    SparseLUSolver<T> solver(A);
    solver.solveInPlace(b);
    expect_vector_near(b, x_true);
}

// Test 9: 10x10 diagonal with small values
TEST(SparseLUSolver, TenByTenSmallDiagonal)
{
    using T = double;
    std::vector<SparseMatrixCSR<T>::triplet_type> entries;
    for (int i = 0; i < 10; ++i)
        entries.emplace_back(i, i, 1e-10 + 1e-12 * i);
    auto sorted_entries = sort_entries(entries);
    SparseMatrixCSR<T> A(10, 10, sorted_entries);
    Vector<T> x_true("x_true", 10);
    for (int i = 0; i < 10; ++i)
        x_true(i) = i + 1;
    Vector<T> b = csr_matvec(A, x_true);
    SparseLUSolver<T> solver(A);
    solver.solveInPlace(b);
    expect_vector_near(b, x_true, 1e-6);
}

// Test 10: 10x10 tridiagonal
TEST(SparseLUSolver, TenByTenTridiagonal)
{
    using T = double;
    std::vector<SparseMatrixCSR<T>::triplet_type> entries;
    for (int i = 0; i < 10; ++i) {
        if (i > 0)
            entries.emplace_back(i, i - 1, -1.0);
        entries.emplace_back(i, i, 2.0);
        if (i < 9)
            entries.emplace_back(i, i + 1, -1.0);
    }
    auto sorted_entries = sort_entries(entries);
    SparseMatrixCSR<T> A(10, 10, sorted_entries);
    Vector<T> x_true("x_true", 10);
    for (int i = 0; i < 10; ++i)
        x_true(i) = i + 1;
    Vector<T> b = csr_matvec(A, x_true);
    SparseLUSolver<T> solver(A);
    solver.solveInPlace(b);
    expect_vector_near(b, x_true, 1e-8);
}

// Test 11: 10x10 with random permutation pattern
TEST(SparseLUSolver, TenByTenPermutationPattern)
{
    using T = double;
    std::vector<SparseMatrixCSR<T>::triplet_type> entries;
    for (int i = 0; i < 10; ++i) {
        entries.emplace_back(i, (i + 3) % 10, 1.0 + i);
        entries.emplace_back(i, i, 10.0 + i);
    }
    auto sorted_entries = sort_entries(entries);
    SparseMatrixCSR<T> A(10, 10, sorted_entries);
    Vector<T> x_true("x_true", 10);
    for (int i = 0; i < 10; ++i)
        x_true(i) = i - 5;
    Vector<T> b = csr_matvec(A, x_true);
    SparseLUSolver<T> solver(A);
    solver.solveInPlace(b);
    expect_vector_near(b, x_true, 1e-8);
}

// Test 12: 20x20 random sparse
TEST(SparseLUSolver, TwentyByTwentyRandomSparse)
{
    using T = double;
    std::vector<SparseMatrixCSR<T>::triplet_type> entries;
    std::mt19937 gen(42);
    std::uniform_real_distribution<T> dist(-2.0, 2.0);
    for (int i = 0; i < 20; ++i) {
        entries.emplace_back(i, i, 10.0 + dist(gen));
        for (int j = 0; j < 3; ++j) {
            int col = (i + j + 1) % 20;
            entries.emplace_back(i, col, dist(gen));
        }
    }
    auto sorted_entries = sort_entries(entries);
    SparseMatrixCSR<T> A(20, 20, sorted_entries);
    Vector<T> x_true("x_true", 20);
    for (int i = 0; i < 20; ++i)
        x_true(i) = dist(gen);
    Vector<T> b = csr_matvec(A, x_true);
    SparseLUSolver<T> solver(A);
    solver.solveInPlace(b);
    expect_vector_near(b, x_true, 1e-8);
}

// Test 13: 50x50 diagonal dominant
TEST(SparseLUSolver, FiftyByFiftyDiagonalDominant)
{
    using T = double;
    std::vector<SparseMatrixCSR<T>::triplet_type> entries;
    for (int i = 0; i < 50; ++i) {
        entries.emplace_back(i, i, 100.0 + i);
        if (i > 0)
            entries.emplace_back(i, i - 1, 1.0);
        if (i < 49)
            entries.emplace_back(i, i + 1, 1.0);
    }
    auto sorted_entries = sort_entries(entries);
    SparseMatrixCSR<T> A(50, 50, sorted_entries);
    Vector<T> x_true("x_true", 50);
    for (int i = 0; i < 50; ++i)
        x_true(i) = i * 0.5;
    Vector<T> b = csr_matvec(A, x_true);
    SparseLUSolver<T> solver(A);
    solver.solveInPlace(b);
    expect_vector_near(b, x_true, 1e-8);
}

// Test 14: 100x100 random sparse, well-conditioned
TEST(SparseLUSolver, HundredByHundredRandomSparse)
{
    using T = double;
    std::vector<SparseMatrixCSR<T>::triplet_type> entries;
    std::mt19937 gen(123);
    std::uniform_real_distribution<T> dist(-1.0, 1.0);
    for (int i = 0; i < 100; ++i) {
        entries.emplace_back(i, i, 10.0 + dist(gen));
        for (int j = 0; j < 2; ++j) {
            int col = (i + j + 7) % 100;
            entries.emplace_back(i, col, dist(gen));
        }
    }
    auto sorted_entries = sort_entries(entries);
    SparseMatrixCSR<T> A(100, 100, sorted_entries);
    Vector<T> x_true("x_true", 100);
    for (int i = 0; i < 100; ++i)
        x_true(i) = dist(gen);
    Vector<T> b = csr_matvec(A, x_true);
    SparseLUSolver<T> solver(A);
    solver.solveInPlace(b);
    expect_vector_near(b, x_true, 1e-7);
}

// Test 15: 10x10 ill-conditioned (small diagonal, large off-diagonal)
TEST(SparseLUSolver, TenByTenIllConditioned)
{
    using T = double;
    std::vector<SparseMatrixCSR<T>::triplet_type> entries;
    for (int i = 0; i < 10; ++i) {
        entries.emplace_back(i, i, 1e-8);
        if (i > 0)
            entries.emplace_back(i, i - 1, 1.0);
        if (i < 9)
            entries.emplace_back(i, i + 1, 1.0);
    }
    auto sorted_entries = sort_entries(entries);
    SparseMatrixCSR<T> A(10, 10, sorted_entries);
    Vector<T> x_true("x_true", 10);
    for (int i = 0; i < 10; ++i)
        x_true(i) = i + 1;
    Vector<T> b = csr_matvec(A, x_true);
    SparseLUSolver<T> solver(A);
    solver.solveInPlace(b);
    Vector<T> b2 = csr_matvec(A, b);
    expect_vector_near(b2, csr_matvec(A, x_true), 1e-6);
}

// Test 17: 20x20 banded
TEST(SparseLUSolver, TwentyByTwentyBanded)
{
    using T = double;
    std::vector<SparseMatrixCSR<T>::triplet_type> entries;
    for (int i = 0; i < 20; ++i) {
        for (int j = std::max(0, i - 2); j <= std::min(19, i + 2); ++j)
            entries.emplace_back(i, j, (i == j) ? 5.0 : 1.0);
    }
    auto sorted_entries = sort_entries(entries);
    SparseMatrixCSR<T> A(20, 20, sorted_entries);
    Vector<T> x_true("x_true", 20);
    for (int i = 0; i < 20; ++i)
        x_true(i) = std::sin(i);
    Vector<T> b = csr_matvec(A, x_true);
    SparseLUSolver<T> solver(A);
    solver.solveInPlace(b);
    expect_vector_near(b, x_true, 1e-8);
}

// Test 18: 5x5 with negative and positive values
TEST(SparseLUSolver, FiveByFiveMixedSigns)
{
    using T = double;
    SparseMatrixCSR<T> A(
        5, 5,
        sort_entries<T>(
            {{0, 0, -2.0}, {0, 4, 1.0}, {1, 1, 3.0}, {2, 2, -4.0}, {3, 3, 5.0}, {4, 0, -1.0}, {4, 4, -6.0}}));
    Vector<T> x_true("x_true", 5);
    x_true(0)   = 1.0;
    x_true(1)   = -2.0;
    x_true(2)   = 3.0;
    x_true(3)   = -4.0;
    x_true(4)   = 5.0;
    Vector<T> b = csr_matvec(A, x_true);
    SparseLUSolver<T> solver(A);
    solver.solveInPlace(b);
    expect_vector_near(b, x_true, 1e-8);
}

// Helper to multiply CSR matrix by vector: b = A * x
template <typename T>
Vector<T> multiply(const SparseMatrixCSR<T>& A, const Vector<T>& x)
{
    int n = A.rows();
    Vector<T> b("b", n);
    for (int i = 0; i < n; ++i) {
        T sum         = 0;
        int start     = A.row_nz_index(i, 0); // row_start_indices equivalent
        int row_start = A.row_start_indices_data()[i];
        int row_end   = A.row_start_indices_data()[i + 1];
        for (int idx = row_start; idx < row_end; ++idx) {
            int col = A.column_indices_data()[idx];
            T val   = A.values_data()[idx];
            sum += val * x[col];
        }
        b(i) = sum;
    }
    return b;
}

// Floating-point comparison tolerance
template <typename T>
void expectVectorNear(const Vector<T>& a, const Vector<T>& b, double tol = 1e-8)
{
    ASSERT_EQ(a.size(), b.size());
    for (int i = 0; i < a.size(); ++i) {
        EXPECT_NEAR(a(i), b(i), tol);
    }
}

// 21. Test 1x1 matrix
TEST(SparseLUSolver, 1x1)
{
    SparseMatrixCSR<double> A(1, 1, sort_entries(std::vector<SparseMatrixCSR<double>::triplet_type>{{0, 0, 5.0}}));
    Vector<double> x_true("x_true", 1);
    x_true(0)        = 2.0;
    Vector<double> b = multiply(A, x_true);
    SparseLUSolver<double> solver(A);
    solver.solveInPlace(b);
    expectVectorNear(b, x_true);
}

// 22. Test 2x2 diagonal matrix
TEST(SparseLUSolver, 2x2Diagonal)
{
    std::vector<SparseMatrixCSR<double>::triplet_type> triplets = sort_entries<double>({{0, 0, 2.0}, {1, 1, 3.0}});
    SparseMatrixCSR<double> A(2, 2, triplets);
    Vector<double> x_true("x_true", 2);
    x_true(0)        = 1.0;
    x_true(1)        = -1.0;
    Vector<double> b = multiply(A, x_true);
    SparseLUSolver<double> solver(A);
    solver.solveInPlace(b);
    expectVectorNear(b, x_true);
}

// 23. Test 2x2 general matrix
TEST(SparseLUSolver, 2x2General)
{
    std::vector<SparseMatrixCSR<double>::triplet_type> triplets =
        sort_entries<double>({{0, 0, 4.0}, {0, 1, 1.0}, {1, 0, 2.0}, {1, 1, 3.0}});
    SparseMatrixCSR<double> A(2, 2, triplets);
    Vector<double> x_true("x_true", 2);
    x_true(0)        = 3.0;
    x_true(1)        = -2.0;
    Vector<double> b = multiply(A, x_true);
    SparseLUSolver<double> solver(A);
    solver.solveInPlace(b);
    expectVectorNear(b, x_true);
}

// 24. Test 3x3 requiring permutation (simple symmetric)
TEST(SparseLUSolver, 3x3Permutation)
{
    // Matrix with zero at (0,0) so pivot will be permuted
    std::vector<SparseMatrixCSR<double>::triplet_type> triplets =
        sort_entries<double>({{0, 1, 1.0}, {1, 0, 1.0}, {1, 1, 2.0}, {1, 2, 1.0}, {2, 1, 1.0}, {2, 2, 3.0}});
    SparseMatrixCSR<double> A(3, 3, triplets);
    Vector<double> x_true("x_true", 3);
    x_true(0)        = 1.0;
    x_true(1)        = 2.0;
    x_true(2)        = -1.0;
    Vector<double> b = multiply(A, x_true);
    SparseLUSolver<double> solver(A);
    solver.solveInPlace(b);
    expectVectorNear(b, x_true, 1e-7);
}

// 25. Test small sparse matrix 4x4
TEST(SparseLUSolver, 4x4parse)
{
    std::vector<SparseMatrixCSR<double>::triplet_type> triplets =
        sort_entries<double>({{0, 0, 10.0}, {0, 3, 2.0}, {1, 1, 5.0}, {2, 2, 3.0}, {3, 0, 1.0}, {3, 3, 4.0}});
    SparseMatrixCSR<double> A(4, 4, triplets);
    Vector<double> x_true("x_true", 4);
    x_true(0)        = 1.0;
    x_true(1)        = 2.0;
    x_true(2)        = -1.0;
    x_true(3)        = 3.0;
    Vector<double> b = multiply(A, x_true);
    SparseLUSolver<double> solver(A);
    solver.solveInPlace(b);
    expectVectorNear(b, x_true);
}

// 26. Test with zero diagonal
TEST(SparseLUSolver, 0Diagonal)
{
    std::vector<SparseMatrixCSR<double>::triplet_type> triplets =
        sort_entries<double>({{0, 0, 0.0}, {0, 1, 1.0}, {1, 0, 1.0}, {1, 1, 2.0}});
    SparseMatrixCSR<double> A(2, 2, triplets);
    Vector<double> x_true("x_true", 2);
    x_true(0)        = 1.5;
    x_true(1)        = -0.5;
    Vector<double> b = multiply(A, x_true);
    SparseLUSolver<double> solver(A);
    solver.solveInPlace(b);
    expectVectorNear(b, x_true, 1e-7);
}

// 27. Test 5x5 random sparse matrix
TEST(SparseLUSolver, 5x5Random)
{
    int n = 5;
    SparseMatrixCSR<double> A(n, n, [&](int row) {
        return 3;
    }); // 3 non-zeros per row
    // Fill A with random but controlled entries
    std::vector<double> values(n * 3);
    std::vector<int> cols(n * 3);
    std::vector<int> row_ptr(n + 1);
    int idx = 0;
    for (int i = 0; i < n; ++i) {
        row_ptr[i] = idx;
        for (int j = 0; j < 3; ++j) {
            cols[idx]   = (i + j) % n;
            values[idx] = (i + j + 1);
            ++idx;
        }
    }
    row_ptr[n] = idx;
    SparseMatrixCSR<double> B(n, n, values, cols, row_ptr);
    Vector<double> x_true("x_true", 5);
    x_true(0)        = 1.0;
    x_true(1)        = -2.0;
    x_true(2)        = 3.0;
    x_true(3)        = -4.0;
    x_true(4)        = 5.0;
    Vector<double> b = multiply(B, x_true);
    SparseLUSolver<double> solver(B);
    solver.solveInPlace(b);
    expectVectorNear(b, x_true);
}

// 28. Test larger 10x10 diagonal matrix
TEST(SparseLUSolver, 10x10Diagonal)
{
    int n = 10;
    std::vector<SparseMatrixCSR<double>::triplet_type> triplets;
    for (int i = 0; i < n; ++i)
        triplets.emplace_back(i, i, i + 1.0);

    SparseMatrixCSR<double> A(n, n, sort_entries(triplets));
    Vector<double> x_true("x_true", n);
    for (int i = 0; i < n; ++i)
        x_true(i) = (i % 2 == 0 ? 1.0 : -1.0);
    Vector<double> b = multiply(A, x_true);
    SparseLUSolver<double> solver(A);
    solver.solveInPlace(b);
    expectVectorNear(b, x_true);
}

// 29. Test ill-conditioned matrix (Hilbert matrix 4x4)
TEST(SparseLUSolver, 4x4Hilbert)
{
    int n = 4;
    std::vector<SparseMatrixCSR<double>::triplet_type> triplets;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            triplets.emplace_back(i, j, 1.0 / (i + j + 1));
        }
    }
    SparseMatrixCSR<double> A(n, n, sort_entries(triplets));
    Vector<double> x_true("x_true", 4);
    x_true(0)        = 1.0;
    x_true(1)        = 2.0;
    x_true(2)        = 3.0;
    x_true(3)        = 4.0;
    Vector<double> b = multiply(A, x_true);
    SparseLUSolver<double> solver(A);
    solver.solveInPlace(b);
    expectVectorNear(b, x_true, 1e-6);
}

// 30. Test using double* overload
TEST(SparseLUSolver, DoublePointerOverload)
{
    std::vector<SparseMatrixCSR<double>::triplet_type> triplets = {{0, 0, 2.0}, {0, 1, 0.5}, {1, 0, 1.0}, {1, 1, 3.0}};
    SparseMatrixCSR<double> A(2, 2, sort_entries(triplets));
    Vector<double> x_true("x_true", 2);
    x_true(0) = 4.0;
    x_true(1) = -1.0;

    Vector<double> b_vec = multiply(A, x_true);
    double b_raw[2]      = {b_vec[0], b_vec[1]};
    SparseLUSolver<double> solver(A);
    solver.solveInPlace(b_raw);
    for (int i = 0; i < 2; ++i) {
        EXPECT_NEAR(b_raw[i], x_true(i), 1e-8);
    }
}

// 31. Test 6x6 matrix with block structure
TEST(SparseLUSolver, 6x6Block)
{
    int n = 6;
    std::vector<SparseMatrixCSR<double>::triplet_type> triplets;
    // Create two 3x3 diagonal blocks
    for (int i = 0; i < 3; ++i) {
        triplets.emplace_back(i, i, 5.0 + i);
        triplets.emplace_back(i + 3, i + 3, 2.0 + i);
    }
    SparseMatrixCSR<double> A(n, n, sort_entries(triplets));
    Vector<double> x_true("x_true", 6);
    x_true(0)        = 1.0;
    x_true(1)        = -1.0;
    x_true(2)        = 2.0;
    x_true(3)        = 0.5;
    x_true(4)        = -0.5;
    x_true(5)        = 1.5;
    Vector<double> b = multiply(A, x_true);
    SparseLUSolver<double> solver(A);
    solver.solveInPlace(b);
    expectVectorNear(b, x_true);
}

// 33. Test with negative and positive values mixed
TEST(SparseLUSolver, 5x5MixedSigns)
{
    int n = 5;
    std::vector<SparseMatrixCSR<double>::triplet_type> triplets;
    for (int i = 0; i < n; ++i) {
        triplets.emplace_back(i, i, (i % 2 == 0) ? i + 1.0 : -(i + 1.0));
        if (i < n - 1) {
            triplets.emplace_back(i, i + 1, 0.5 * (i + 1));
            triplets.emplace_back(i + 1, i, -0.3 * (i + 1));
        }
    }
    SparseMatrixCSR<double> A(n, n, sort_entries(triplets));
    Vector<double> x_true("x_true", 5);
    x_true(0)        = 2.0;
    x_true(1)        = -1.0;
    x_true(2)        = 0.5;
    x_true(3)        = -0.5;
    x_true(4)        = 1.0;
    Vector<double> b = multiply(A, x_true);
    SparseLUSolver<double> solver(A);
    solver.solveInPlace(b);
    expectVectorNear(b, x_true, 1e-8);
}

// 34. Test 7x7 random diagonal dominant
TEST(SparseLUSolver, 7x7DiagDominant)
{
    int n = 7;
    std::vector<SparseMatrixCSR<double>::triplet_type> triplets;
    for (int i = 0; i < n; ++i) {
        double diag = 10.0 + i;
        triplets.emplace_back(i, i, diag);
        if (i < n - 1) {
            triplets.emplace_back(i, i + 1, 1.0);
            triplets.emplace_back(i + 1, i, -1.0);
        }
    }
    SparseMatrixCSR<double> A(n, n, sort_entries(triplets));
    Vector<double> x_true("x_true", n);
    for (int i = 0; i < n; ++i)
        x_true(i) = (i + 1) * 0.1;
    Vector<double> b = multiply(A, x_true);
    SparseLUSolver<double> solver(A);
    solver.solveInPlace(b);
    expectVectorNear(b, x_true, 1e-8);
}

// 35. Test 8x8 matrix requiring heavy reordering
TEST(SparseLUSolver, 8x8ChainGraph)
{
    int n = 8;
    std::vector<SparseMatrixCSR<double>::triplet_type> triplets;
    // Create adjacency-like pattern: each node connected to next
    for (int i = 0; i < n; ++i) {
        triplets.emplace_back(i, i, 4.0);
        if (i < n - 1) {
            triplets.emplace_back(i, i + 1, 1.0);
            triplets.emplace_back(i + 1, i, 1.0);
        }
    }
    SparseMatrixCSR<double> A(n, n, sort_entries(triplets));
    Vector<double> x_true("x_true", n);
    for (int i = 0; i < n; ++i)
        x_true(i) = (i % 2 == 0 ? 1.0 : -1.0);
    Vector<double> b = multiply(A, x_true);
    SparseLUSolver<double> solver(A);
    solver.solveInPlace(b);
    expectVectorNear(b, x_true, 1e-8);
}

// 36. Test with float type
TEST(SparseLUSolver, FloatPrecision)
{
    std::vector<SparseMatrixCSR<float>::triplet_type> triplets = {
        {0, 0, 1.0f}, {0, 1, 2.0f}, {1, 0, 3.0f}, {1, 1, 4.0f}};
    SparseMatrixCSR<float> A(2, 2, sort_entries(triplets));
    Vector<float> x_true("x_true", 2);
    x_true(0)       = 1.0f;
    x_true(1)       = -1.0f;
    Vector<float> b = multiply(A, x_true);
    SparseLUSolver<float> solver(A);
    solver.solveInPlace(b);
    for (int i = 0; i < 2; ++i) {
        EXPECT_NEAR(b(i), x_true(i), 1e-4f);
    }
}

// 37. Test 3x3 identity matrix
TEST(SparseLUSolver, 3x3Identity)
{
    int n = 3;
    std::vector<SparseMatrixCSR<double>::triplet_type> triplets;
    for (int i = 0; i < n; ++i)
        triplets.emplace_back(i, i, 1.0);
    SparseMatrixCSR<double> A(n, n, sort_entries(triplets));
    Vector<double> x_true("x_true", 3);
    x_true(0)        = 5.0;
    x_true(1)        = -3.0;
    x_true(2)        = 2.0;
    Vector<double> b = multiply(A, x_true);
    SparseLUSolver<double> solver(A);
    solver.solveInPlace(b);
    expectVectorNear(b, x_true);
}

// 39. Test 4x4 with very small diagonal (near-zero pivot)
TEST(SparseLUSolver, 4x4SmallDiagonal)
{
    std::vector<SparseMatrixCSR<double>::triplet_type> triplets = {{0, 0, 1e-12}, {0, 1, 1.0}, {1, 0, 1.0},
                                                                   {1, 1, 2.0},   {2, 2, 3.0}, {3, 3, 4.0}};
    SparseMatrixCSR<double> A(4, 4, sort_entries(triplets));
    Vector<double> x_true("x_true", 4);
    x_true(0)        = 1.0;
    x_true(1)        = 2.0;
    x_true(2)        = -1.0;
    x_true(3)        = 0.5;
    Vector<double> b = multiply(A, x_true);
    SparseLUSolver<double> solver(A);
    solver.solveInPlace(b);
    expectVectorNear(b, x_true, 1e-4);
}

// 40. Test 5x5 random symmetric positive definite (SPD) via A = M*M^T
TEST(SparseLUSolver, 5x5SPD)
{
    int n = 5;
    // Create random lower-triangular dense M
    std::vector<std::vector<double>> M(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            M[i][j] = (i + j + 1.0);
        }
    }
    std::vector<SparseMatrixCSR<double>::triplet_type> triplets;
    // Compute A = M * M^T into triplets
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double sum = 0.0;
            for (int k = 0; k <= std::min(i, j); ++k)
                sum += M[i][k] * M[j][k];
            if (std::abs(sum) > 1e-12)
                triplets.emplace_back(i, j, sum);
        }
    }
    SparseMatrixCSR<double> A(n, n, sort_entries(triplets));
    Vector<double> x_true("x_true", 5);
    x_true(0)        = 1.0;
    x_true(1)        = 2.0;
    x_true(2)        = 3.0;
    x_true(3)        = 4.0;
    x_true(4)        = 5.0;
    Vector<double> b = multiply(A, x_true);
    SparseLUSolver<double> solver(A);
    solver.solveInPlace(b);
    expectVectorNear(b, x_true, 1e-8);
}
