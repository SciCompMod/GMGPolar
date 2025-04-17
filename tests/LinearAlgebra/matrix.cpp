#include <gtest/gtest.h>

#include "../../include/LinearAlgebra/coo_matrix.h"

// Alias for readability
using triplet = std::tuple<int, int, double>;

/* Default Construct (double) */

TEST(SparseMatrixDouble, default_construct)
{
    const SparseMatrixCOO<double> m;
    (void)m;
}

/* Size Construct (double) */

TEST(SparseMatrixDouble, size_construct)
{
    const SparseMatrixCOO<double> m(5, 5, 0);
    ASSERT_EQ(m.rows(), 5);
    ASSERT_EQ(m.columns(), 5);
    ASSERT_EQ(m.non_zero_size(), 0);
}

/* Entries Construct (double) */

TEST(SparseMatrixDouble, entries_construct)
{
    std::vector<triplet> entries = {{0, 0, 1.0}, {1, 2, 2.0}, {2, 1, 3.0}};
    const SparseMatrixCOO<double> m(3, 3, entries);

    ASSERT_EQ(m.rows(), 3);
    ASSERT_EQ(m.columns(), 3);
    ASSERT_EQ(m.non_zero_size(), 3);

    EXPECT_EQ(m.row_index(0), 0);
    EXPECT_EQ(m.col_index(0), 0);
    EXPECT_DOUBLE_EQ(m.value(0), 1.0);

    EXPECT_EQ(m.row_index(1), 1);
    EXPECT_EQ(m.col_index(1), 2);
    EXPECT_DOUBLE_EQ(m.value(1), 2.0);

    EXPECT_EQ(m.row_index(2), 2);
    EXPECT_EQ(m.col_index(2), 1);
    EXPECT_DOUBLE_EQ(m.value(2), 3.0);
}

/* Copy Construct (double) */

TEST(SparseMatrixDouble, copy_construct)
{
    std::vector<triplet> entries = {{0, 0, 1.0}, {1, 2, 2.0}, {2, 1, 3.0}};
    const SparseMatrixCOO<double> m1(3, 3, entries);
    const SparseMatrixCOO<double> m2 = m1;

    ASSERT_EQ(m2.rows(), 3);
    ASSERT_EQ(m2.columns(), 3);
    ASSERT_EQ(m2.non_zero_size(), 3);

    EXPECT_EQ(m2.row_index(0), 0);
    EXPECT_EQ(m2.col_index(0), 0);
    EXPECT_DOUBLE_EQ(m2.value(0), 1.0);

    EXPECT_EQ(m2.row_index(1), 1);
    EXPECT_EQ(m2.col_index(1), 2);
    EXPECT_DOUBLE_EQ(m2.value(1), 2.0);

    EXPECT_EQ(m2.row_index(2), 2);
    EXPECT_EQ(m2.col_index(2), 1);
    EXPECT_DOUBLE_EQ(m2.value(2), 3.0);
}

/* Move Construct (double) */

TEST(SparseMatrixDouble, move_construct)
{
    std::vector<triplet> entries = {{0, 0, 1.0}, {1, 2, 2.0}, {2, 1, 3.0}};
    SparseMatrixCOO<double> m1(3, 3, entries);

    ASSERT_EQ(m1.rows(), 3);
    ASSERT_EQ(m1.columns(), 3);
    ASSERT_EQ(m1.non_zero_size(), 3);

    SparseMatrixCOO<double> m2 = std::move(m1);

    ASSERT_EQ(m2.rows(), 3);
    ASSERT_EQ(m2.columns(), 3);
    ASSERT_EQ(m2.non_zero_size(), 3);

    EXPECT_EQ(m2.row_index(0), 0);
    EXPECT_EQ(m2.col_index(0), 0);
    EXPECT_DOUBLE_EQ(m2.value(0), 1.0);

    EXPECT_EQ(m2.row_index(1), 1);
    EXPECT_EQ(m2.col_index(1), 2);
    EXPECT_DOUBLE_EQ(m2.value(1), 2.0);

    EXPECT_EQ(m2.row_index(2), 2);
    EXPECT_EQ(m2.col_index(2), 1);
    EXPECT_DOUBLE_EQ(m2.value(2), 3.0);
}

/* Copy Assign (double) */

TEST(SparseMatrixDouble, copy_assign)
{
    std::vector<triplet> entries = {{0, 0, 1.0}, {1, 2, 2.0}, {2, 1, 3.0}};
    const SparseMatrixCOO<double> m1(3, 3, entries);
    SparseMatrixCOO<double> m2;

    m2 = m1;

    ASSERT_EQ(m2.rows(), 3);
    ASSERT_EQ(m2.columns(), 3);
    ASSERT_EQ(m2.non_zero_size(), 3);

    EXPECT_EQ(m2.row_index(0), 0);
    EXPECT_EQ(m2.col_index(0), 0);
    EXPECT_DOUBLE_EQ(m2.value(0), 1.0);

    EXPECT_EQ(m2.row_index(1), 1);
    EXPECT_EQ(m2.col_index(1), 2);
    EXPECT_DOUBLE_EQ(m2.value(1), 2.0);

    EXPECT_EQ(m2.row_index(2), 2);
    EXPECT_EQ(m2.col_index(2), 1);
    EXPECT_DOUBLE_EQ(m2.value(2), 3.0);
}

/* Move Assign (double) */

TEST(SparseMatrixDouble, move_assign)
{
    std::vector<triplet> entries = {{0, 0, 1.0}, {1, 2, 2.0}, {2, 1, 3.0}};
    SparseMatrixCOO<double> m1(3, 3, entries);
    SparseMatrixCOO<double> m2;

    m2 = std::move(m1);

    ASSERT_EQ(m2.rows(), 3);
    ASSERT_EQ(m2.columns(), 3);
    ASSERT_EQ(m2.non_zero_size(), 3);

    EXPECT_EQ(m2.row_index(0), 0);
    EXPECT_EQ(m2.col_index(0), 0);
    EXPECT_DOUBLE_EQ(m2.value(0), 1.0);

    EXPECT_EQ(m2.row_index(1), 1);
    EXPECT_EQ(m2.col_index(1), 2);
    EXPECT_DOUBLE_EQ(m2.value(1), 2.0);

    EXPECT_EQ(m2.row_index(2), 2);
    EXPECT_EQ(m2.col_index(2), 1);
    EXPECT_DOUBLE_EQ(m2.value(2), 3.0);
}

/* Value Modify (double) */

TEST(SparseMatrixDouble, value_modify)
{
    std::vector<triplet> entries = {{0, 0, 1.0}, {1, 2, 2.0}, {2, 1, 3.0}};
    SparseMatrixCOO<double> m(3, 3, entries);

    ASSERT_EQ(m.rows(), 3);
    ASSERT_EQ(m.columns(), 3);
    ASSERT_EQ(m.non_zero_size(), 3);

    EXPECT_DOUBLE_EQ(m.value(2), 3.0);

    m.value(2) = 5.0;
    EXPECT_DOUBLE_EQ(m.value(2), 5.0);
}

/* Modify Construction (double) */

TEST(SparseMatrixDouble, modify_construction)
{
    const int rows    = 4;
    const int columns = 5;
    const int nnz     = 6;
    SparseMatrixCOO<double> matrixA(rows, columns, nnz);

    matrixA.row_index(0) = 0;
    matrixA.col_index(0) = 1;
    matrixA.value(0)     = 2.0;

    matrixA.row_index(1) = 1;
    matrixA.col_index(1) = 0;
    matrixA.value(1)     = 8.0;

    matrixA.row_index(2) = 1;
    matrixA.col_index(2) = 3;
    matrixA.value(2)     = -6.0;

    matrixA.row_index(3) = 2;
    matrixA.col_index(3) = 0;
    matrixA.value(3)     = -5.0;

    matrixA.row_index(4) = 3;
    matrixA.col_index(4) = 1;
    matrixA.value(4)     = 7.0;

    matrixA.row_index(5) = 3;
    matrixA.col_index(5) = 3;
    matrixA.value(5)     = 3.0;

    ASSERT_EQ(matrixA.rows(), 4);
    ASSERT_EQ(matrixA.columns(), 5);
    ASSERT_EQ(matrixA.non_zero_size(), 6);

    std::vector<triplet> entries = {{0, 1, 2.0}, {1, 0, 8.0}, {1, 3, -6.0}, {2, 0, -5.0}, {3, 1, 7.0}, {3, 3, 3.0}};

    SparseMatrixCOO<double> matrixB(rows, columns, entries);

    for (int i = 0; i < nnz; i++) {
        ASSERT_EQ(matrixA.row_index(i), matrixB.row_index(i));
        ASSERT_EQ(matrixA.col_index(i), matrixB.col_index(i));
        EXPECT_DOUBLE_EQ(matrixA.value(i), matrixB.value(i));
    }
}

/* Symmetry Check (double) */

TEST(SparseMatrixDouble, symmetry_check)
{
    std::vector<triplet> entries = {{0, 1, 2.0}, {1, 1, -3.0}, {1, 2, 3.0}};
    SparseMatrixCOO<double> m(3, 3, entries);

    ASSERT_FALSE(m.is_symmetric());
    m.is_symmetric(true);
    ASSERT_TRUE(m.is_symmetric());
}

