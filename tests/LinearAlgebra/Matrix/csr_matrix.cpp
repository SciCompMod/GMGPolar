#include <gtest/gtest.h>

#include "../../../include/LinearAlgebra/Matrix/csr_matrix.h"
#include "../../../include/LinearAlgebra/Solvers/csr_lu_solver.h"
#include "../../../include/LinearAlgebra/Vector/vector.h"

TEST(SparseMatrixCSR, default_construct)
{
    const SparseMatrixCSR<double> v;
    (void)v;
}

TEST(SparseMatrixCSR, value_construct)
{
    using triplet = SparseMatrixCSR<double>::triplet_type;

    const SparseMatrixCSR<double> m(4, 3,
                                    {triplet{0, 0, 1.0}, triplet{0, 2, 2.0}, triplet{1, 0, 3.0}, triplet{2, 0, 4.0},
                                     triplet{2, 1, 5.0}, triplet{2, 2, 6.0}, triplet{3, 2, 7.0}});

    EXPECT_EQ(m.rows(), 4);
    EXPECT_EQ(m.columns(), 3);
    EXPECT_EQ(m.non_zero_size(), 7);

    EXPECT_EQ(m.row_nz_size(0), 2);
    EXPECT_EQ(m.row_nz_size(1), 1);
    EXPECT_EQ(m.row_nz_size(2), 3);
    EXPECT_EQ(m.row_nz_size(3), 1);

    const double dense_values[4][3] = {{1.0, 0.0, 2.0}, {3.0, 0.0, 0.0}, {4.0, 5.0, 6.0}, {0.0, 0.0, 7.0}};

    EXPECT_DOUBLE_EQ(m.row_nz_entry(0, 0), dense_values[0][m.row_nz_index(0, 0)]);
    EXPECT_DOUBLE_EQ(m.row_nz_entry(0, 1), dense_values[0][m.row_nz_index(0, 1)]);
    EXPECT_DOUBLE_EQ(m.row_nz_entry(1, 0), dense_values[1][m.row_nz_index(1, 0)]);
    EXPECT_DOUBLE_EQ(m.row_nz_entry(2, 0), dense_values[2][m.row_nz_index(2, 0)]);
    EXPECT_DOUBLE_EQ(m.row_nz_entry(2, 1), dense_values[2][m.row_nz_index(2, 1)]);
    EXPECT_DOUBLE_EQ(m.row_nz_entry(2, 2), dense_values[2][m.row_nz_index(2, 2)]);
    EXPECT_DOUBLE_EQ(m.row_nz_entry(3, 0), dense_values[3][m.row_nz_index(3, 0)]);
}

TEST(SparseMatrixCSR, copy_construct)
{
    using triplet = SparseMatrixCSR<double>::triplet_type;

    const SparseMatrixCSR<double> m1(4, 3,
                                     {triplet{0, 0, 1.0}, triplet{0, 2, 2.0}, triplet{1, 0, 3.0}, triplet{2, 0, 4.0},
                                      triplet{2, 1, 5.0}, triplet{2, 2, 6.0}, triplet{3, 2, 7.0}});
    const SparseMatrixCSR<double> m2 = m1;

    const double dense_values[4][3] = {{1.0, 0.0, 2.0}, {3.0, 0.0, 0.0}, {4.0, 5.0, 6.0}, {0.0, 0.0, 7.0}};

    EXPECT_EQ(m1.rows(), 4);
    EXPECT_EQ(m1.columns(), 3);
    EXPECT_EQ(m1.non_zero_size(), 7);

    EXPECT_EQ(m1.row_nz_size(0), 2);
    EXPECT_EQ(m1.row_nz_size(1), 1);
    EXPECT_EQ(m1.row_nz_size(2), 3);
    EXPECT_EQ(m1.row_nz_size(3), 1);

    EXPECT_DOUBLE_EQ(m1.row_nz_entry(0, 0), dense_values[0][m1.row_nz_index(0, 0)]);
    EXPECT_DOUBLE_EQ(m1.row_nz_entry(0, 1), dense_values[0][m1.row_nz_index(0, 1)]);
    EXPECT_DOUBLE_EQ(m1.row_nz_entry(1, 0), dense_values[1][m1.row_nz_index(1, 0)]);
    EXPECT_DOUBLE_EQ(m1.row_nz_entry(2, 0), dense_values[2][m1.row_nz_index(2, 0)]);
    EXPECT_DOUBLE_EQ(m1.row_nz_entry(2, 1), dense_values[2][m1.row_nz_index(2, 1)]);
    EXPECT_DOUBLE_EQ(m1.row_nz_entry(2, 2), dense_values[2][m1.row_nz_index(2, 2)]);
    EXPECT_DOUBLE_EQ(m1.row_nz_entry(3, 0), dense_values[3][m1.row_nz_index(3, 0)]);

    EXPECT_EQ(m2.rows(), 4);
    EXPECT_EQ(m2.columns(), 3);
    EXPECT_EQ(m2.non_zero_size(), 7);

    EXPECT_EQ(m2.row_nz_size(0), 2);
    EXPECT_EQ(m2.row_nz_size(1), 1);
    EXPECT_EQ(m2.row_nz_size(2), 3);
    EXPECT_EQ(m2.row_nz_size(3), 1);

    EXPECT_DOUBLE_EQ(m2.row_nz_entry(0, 0), dense_values[0][m2.row_nz_index(0, 0)]);
    EXPECT_DOUBLE_EQ(m2.row_nz_entry(0, 1), dense_values[0][m2.row_nz_index(0, 1)]);
    EXPECT_DOUBLE_EQ(m2.row_nz_entry(1, 0), dense_values[1][m2.row_nz_index(1, 0)]);
    EXPECT_DOUBLE_EQ(m2.row_nz_entry(2, 0), dense_values[2][m2.row_nz_index(2, 0)]);
    EXPECT_DOUBLE_EQ(m2.row_nz_entry(2, 1), dense_values[2][m2.row_nz_index(2, 1)]);
    EXPECT_DOUBLE_EQ(m2.row_nz_entry(2, 2), dense_values[2][m2.row_nz_index(2, 2)]);
    EXPECT_DOUBLE_EQ(m2.row_nz_entry(3, 0), dense_values[3][m2.row_nz_index(3, 0)]);
}

TEST(SparseMatrixCSR, copy_assign)
{
    const double dense_values[4][3] = {{1.0, 0.0, 2.0}, {3.0, 0.0, 0.0}, {4.0, 5.0, 6.0}, {0.0, 0.0, 7.0}};

    SparseMatrixCSR<double> m2;

    {
        using triplet = SparseMatrixCSR<double>::triplet_type;

        const SparseMatrixCSR<double> m1(4, 3,
                                         {triplet{0, 0, 1.0}, triplet{0, 2, 2.0}, triplet{1, 0, 3.0},
                                          triplet{2, 0, 4.0}, triplet{2, 1, 5.0}, triplet{2, 2, 6.0},
                                          triplet{3, 2, 7.0}});

        m2 = m1;

        EXPECT_EQ(m1.rows(), 4);
        EXPECT_EQ(m1.columns(), 3);
        EXPECT_EQ(m1.non_zero_size(), 7);

        EXPECT_EQ(m1.row_nz_size(0), 2);
        EXPECT_EQ(m1.row_nz_size(1), 1);
        EXPECT_EQ(m1.row_nz_size(2), 3);
        EXPECT_EQ(m1.row_nz_size(3), 1);

        EXPECT_DOUBLE_EQ(m1.row_nz_entry(0, 0), dense_values[0][m1.row_nz_index(0, 0)]);
        EXPECT_DOUBLE_EQ(m1.row_nz_entry(0, 1), dense_values[0][m1.row_nz_index(0, 1)]);
        EXPECT_DOUBLE_EQ(m1.row_nz_entry(1, 0), dense_values[1][m1.row_nz_index(1, 0)]);
        EXPECT_DOUBLE_EQ(m1.row_nz_entry(2, 0), dense_values[2][m1.row_nz_index(2, 0)]);
        EXPECT_DOUBLE_EQ(m1.row_nz_entry(2, 1), dense_values[2][m1.row_nz_index(2, 1)]);
        EXPECT_DOUBLE_EQ(m1.row_nz_entry(2, 2), dense_values[2][m1.row_nz_index(2, 2)]);
        EXPECT_DOUBLE_EQ(m1.row_nz_entry(3, 0), dense_values[3][m1.row_nz_index(3, 0)]);

        EXPECT_EQ(m2.rows(), 4);
        EXPECT_EQ(m2.columns(), 3);
        EXPECT_EQ(m2.non_zero_size(), 7);

        EXPECT_EQ(m2.row_nz_size(0), 2);
        EXPECT_EQ(m2.row_nz_size(1), 1);
        EXPECT_EQ(m2.row_nz_size(2), 3);
        EXPECT_EQ(m2.row_nz_size(3), 1);

        EXPECT_DOUBLE_EQ(m2.row_nz_entry(0, 0), dense_values[0][m2.row_nz_index(0, 0)]);
        EXPECT_DOUBLE_EQ(m2.row_nz_entry(0, 1), dense_values[0][m2.row_nz_index(0, 1)]);
        EXPECT_DOUBLE_EQ(m2.row_nz_entry(1, 0), dense_values[1][m2.row_nz_index(1, 0)]);
        EXPECT_DOUBLE_EQ(m2.row_nz_entry(2, 0), dense_values[2][m2.row_nz_index(2, 0)]);
        EXPECT_DOUBLE_EQ(m2.row_nz_entry(2, 1), dense_values[2][m2.row_nz_index(2, 1)]);
        EXPECT_DOUBLE_EQ(m2.row_nz_entry(2, 2), dense_values[2][m2.row_nz_index(2, 2)]);
        EXPECT_DOUBLE_EQ(m2.row_nz_entry(3, 0), dense_values[3][m2.row_nz_index(3, 0)]);
    }

    EXPECT_EQ(m2.rows(), 4);
    EXPECT_EQ(m2.columns(), 3);
    EXPECT_EQ(m2.non_zero_size(), 7);

    EXPECT_EQ(m2.row_nz_size(0), 2);
    EXPECT_EQ(m2.row_nz_size(1), 1);
    EXPECT_EQ(m2.row_nz_size(2), 3);
    EXPECT_EQ(m2.row_nz_size(3), 1);

    EXPECT_DOUBLE_EQ(m2.row_nz_entry(0, 0), dense_values[0][m2.row_nz_index(0, 0)]);
    EXPECT_DOUBLE_EQ(m2.row_nz_entry(0, 1), dense_values[0][m2.row_nz_index(0, 1)]);
    EXPECT_DOUBLE_EQ(m2.row_nz_entry(1, 0), dense_values[1][m2.row_nz_index(1, 0)]);
    EXPECT_DOUBLE_EQ(m2.row_nz_entry(2, 0), dense_values[2][m2.row_nz_index(2, 0)]);
    EXPECT_DOUBLE_EQ(m2.row_nz_entry(2, 1), dense_values[2][m2.row_nz_index(2, 1)]);
    EXPECT_DOUBLE_EQ(m2.row_nz_entry(2, 2), dense_values[2][m2.row_nz_index(2, 2)]);
    EXPECT_DOUBLE_EQ(m2.row_nz_entry(3, 0), dense_values[3][m2.row_nz_index(3, 0)]);
}

TEST(SparseMatrixCSR, move_construct)
{
    const double dense_values[4][3] = {{1.0, 0.0, 2.0}, {3.0, 0.0, 0.0}, {4.0, 5.0, 6.0}, {0.0, 0.0, 7.0}};

    using triplet = SparseMatrixCSR<double>::triplet_type;

    const SparseMatrixCSR<double> m1(4, 3,
                                     {triplet{0, 0, 1.0}, triplet{0, 2, 2.0}, triplet{1, 0, 3.0}, triplet{2, 0, 4.0},
                                      triplet{2, 1, 5.0}, triplet{2, 2, 6.0}, triplet{3, 2, 7.0}});

    EXPECT_EQ(m1.rows(), 4);
    EXPECT_EQ(m1.columns(), 3);
    EXPECT_EQ(m1.non_zero_size(), 7);

    EXPECT_EQ(m1.row_nz_size(0), 2);
    EXPECT_EQ(m1.row_nz_size(1), 1);
    EXPECT_EQ(m1.row_nz_size(2), 3);
    EXPECT_EQ(m1.row_nz_size(3), 1);

    EXPECT_DOUBLE_EQ(m1.row_nz_entry(0, 0), dense_values[0][m1.row_nz_index(0, 0)]);
    EXPECT_DOUBLE_EQ(m1.row_nz_entry(0, 1), dense_values[0][m1.row_nz_index(0, 1)]);
    EXPECT_DOUBLE_EQ(m1.row_nz_entry(1, 0), dense_values[1][m1.row_nz_index(1, 0)]);
    EXPECT_DOUBLE_EQ(m1.row_nz_entry(2, 0), dense_values[2][m1.row_nz_index(2, 0)]);
    EXPECT_DOUBLE_EQ(m1.row_nz_entry(2, 1), dense_values[2][m1.row_nz_index(2, 1)]);
    EXPECT_DOUBLE_EQ(m1.row_nz_entry(2, 2), dense_values[2][m1.row_nz_index(2, 2)]);
    EXPECT_DOUBLE_EQ(m1.row_nz_entry(3, 0), dense_values[3][m1.row_nz_index(3, 0)]);

    const SparseMatrixCSR<double> m2 = std::move(m1);

    EXPECT_EQ(m2.rows(), 4);
    EXPECT_EQ(m2.columns(), 3);
    EXPECT_EQ(m2.non_zero_size(), 7);

    EXPECT_EQ(m2.row_nz_size(0), 2);
    EXPECT_EQ(m2.row_nz_size(1), 1);
    EXPECT_EQ(m2.row_nz_size(2), 3);
    EXPECT_EQ(m2.row_nz_size(3), 1);

    EXPECT_DOUBLE_EQ(m2.row_nz_entry(0, 0), dense_values[0][m2.row_nz_index(0, 0)]);
    EXPECT_DOUBLE_EQ(m2.row_nz_entry(0, 1), dense_values[0][m2.row_nz_index(0, 1)]);
    EXPECT_DOUBLE_EQ(m2.row_nz_entry(1, 0), dense_values[1][m2.row_nz_index(1, 0)]);
    EXPECT_DOUBLE_EQ(m2.row_nz_entry(2, 0), dense_values[2][m2.row_nz_index(2, 0)]);
    EXPECT_DOUBLE_EQ(m2.row_nz_entry(2, 1), dense_values[2][m2.row_nz_index(2, 1)]);
    EXPECT_DOUBLE_EQ(m2.row_nz_entry(2, 2), dense_values[2][m2.row_nz_index(2, 2)]);
    EXPECT_DOUBLE_EQ(m2.row_nz_entry(3, 0), dense_values[3][m2.row_nz_index(3, 0)]);
}

TEST(SparseMatrixCSR, move_assign)
{
    const double dense_values[4][3] = {{1.0, 0.0, 2.0}, {3.0, 0.0, 0.0}, {4.0, 5.0, 6.0}, {0.0, 0.0, 7.0}};

    using triplet = SparseMatrixCSR<double>::triplet_type;

    SparseMatrixCSR<double> m2;

    {
        const SparseMatrixCSR<double> m1(4, 3,
                                         {triplet{0, 0, 1.0}, triplet{0, 2, 2.0}, triplet{1, 0, 3.0},
                                          triplet{2, 0, 4.0}, triplet{2, 1, 5.0}, triplet{2, 2, 6.0},
                                          triplet{3, 2, 7.0}});

        EXPECT_EQ(m1.rows(), 4);
        EXPECT_EQ(m1.columns(), 3);
        EXPECT_EQ(m1.non_zero_size(), 7);

        EXPECT_EQ(m1.row_nz_size(0), 2);
        EXPECT_EQ(m1.row_nz_size(1), 1);
        EXPECT_EQ(m1.row_nz_size(2), 3);
        EXPECT_EQ(m1.row_nz_size(3), 1);

        EXPECT_DOUBLE_EQ(m1.row_nz_entry(0, 0), dense_values[0][m1.row_nz_index(0, 0)]);
        EXPECT_DOUBLE_EQ(m1.row_nz_entry(0, 1), dense_values[0][m1.row_nz_index(0, 1)]);
        EXPECT_DOUBLE_EQ(m1.row_nz_entry(1, 0), dense_values[1][m1.row_nz_index(1, 0)]);
        EXPECT_DOUBLE_EQ(m1.row_nz_entry(2, 0), dense_values[2][m1.row_nz_index(2, 0)]);
        EXPECT_DOUBLE_EQ(m1.row_nz_entry(2, 1), dense_values[2][m1.row_nz_index(2, 1)]);
        EXPECT_DOUBLE_EQ(m1.row_nz_entry(2, 2), dense_values[2][m1.row_nz_index(2, 2)]);
        EXPECT_DOUBLE_EQ(m1.row_nz_entry(3, 0), dense_values[3][m1.row_nz_index(3, 0)]);

        m2 = std::move(m1);

        EXPECT_EQ(m2.rows(), 4);
        EXPECT_EQ(m2.columns(), 3);
        EXPECT_EQ(m2.non_zero_size(), 7);

        EXPECT_EQ(m2.row_nz_size(0), 2);
        EXPECT_EQ(m2.row_nz_size(1), 1);
        EXPECT_EQ(m2.row_nz_size(2), 3);
        EXPECT_EQ(m2.row_nz_size(3), 1);

        EXPECT_DOUBLE_EQ(m2.row_nz_entry(0, 0), dense_values[0][m2.row_nz_index(0, 0)]);
        EXPECT_DOUBLE_EQ(m2.row_nz_entry(0, 1), dense_values[0][m2.row_nz_index(0, 1)]);
        EXPECT_DOUBLE_EQ(m2.row_nz_entry(1, 0), dense_values[1][m2.row_nz_index(1, 0)]);
        EXPECT_DOUBLE_EQ(m2.row_nz_entry(2, 0), dense_values[2][m2.row_nz_index(2, 0)]);
        EXPECT_DOUBLE_EQ(m2.row_nz_entry(2, 1), dense_values[2][m2.row_nz_index(2, 1)]);
        EXPECT_DOUBLE_EQ(m2.row_nz_entry(2, 2), dense_values[2][m2.row_nz_index(2, 2)]);
        EXPECT_DOUBLE_EQ(m2.row_nz_entry(3, 0), dense_values[3][m2.row_nz_index(3, 0)]);
    }

    EXPECT_EQ(m2.rows(), 4);
    EXPECT_EQ(m2.columns(), 3);
    EXPECT_EQ(m2.non_zero_size(), 7);

    EXPECT_EQ(m2.row_nz_size(0), 2);
    EXPECT_EQ(m2.row_nz_size(1), 1);
    EXPECT_EQ(m2.row_nz_size(2), 3);
    EXPECT_EQ(m2.row_nz_size(3), 1);

    EXPECT_DOUBLE_EQ(m2.row_nz_entry(0, 0), dense_values[0][m2.row_nz_index(0, 0)]);
    EXPECT_DOUBLE_EQ(m2.row_nz_entry(0, 1), dense_values[0][m2.row_nz_index(0, 1)]);
    EXPECT_DOUBLE_EQ(m2.row_nz_entry(1, 0), dense_values[1][m2.row_nz_index(1, 0)]);
    EXPECT_DOUBLE_EQ(m2.row_nz_entry(2, 0), dense_values[2][m2.row_nz_index(2, 0)]);
    EXPECT_DOUBLE_EQ(m2.row_nz_entry(2, 1), dense_values[2][m2.row_nz_index(2, 1)]);
    EXPECT_DOUBLE_EQ(m2.row_nz_entry(2, 2), dense_values[2][m2.row_nz_index(2, 2)]);
    EXPECT_DOUBLE_EQ(m2.row_nz_entry(3, 0), dense_values[3][m2.row_nz_index(3, 0)]);
}

TEST(SparseMatrixCSR, value_construct_modify)
{
    using triplet = SparseMatrixCSR<double>::triplet_type;

    SparseMatrixCSR<double> m(4, 3,
                              {triplet{0, 0, 1.0}, triplet{0, 2, 2.0}, triplet{1, 0, 3.0}, triplet{2, 0, 4.0},
                               triplet{2, 1, 5.0}, triplet{2, 2, 6.0}, triplet{3, 2, 7.0}});

    EXPECT_EQ(m.rows(), 4);
    EXPECT_EQ(m.columns(), 3);
    EXPECT_EQ(m.non_zero_size(), 7);

    EXPECT_EQ(m.row_nz_size(0), 2);
    EXPECT_EQ(m.row_nz_size(1), 1);
    EXPECT_EQ(m.row_nz_size(2), 3);
    EXPECT_EQ(m.row_nz_size(3), 1);

    const double dense_values[4][3] = {{1.0, 0.0, 2.0}, {3.0, 0.0, 0.0}, {4.0, 5.0, 6.0}, {0.0, 0.0, 7.0}};

    EXPECT_DOUBLE_EQ(m.row_nz_entry(0, 0), dense_values[0][m.row_nz_index(0, 0)]);
    EXPECT_DOUBLE_EQ(m.row_nz_entry(0, 1), dense_values[0][m.row_nz_index(0, 1)]);
    EXPECT_DOUBLE_EQ(m.row_nz_entry(1, 0), dense_values[1][m.row_nz_index(1, 0)]);
    EXPECT_DOUBLE_EQ(m.row_nz_entry(2, 0), dense_values[2][m.row_nz_index(2, 0)]);
    EXPECT_DOUBLE_EQ(m.row_nz_entry(2, 1), dense_values[2][m.row_nz_index(2, 1)]);
    EXPECT_DOUBLE_EQ(m.row_nz_entry(2, 2), dense_values[2][m.row_nz_index(2, 2)]);
    EXPECT_DOUBLE_EQ(m.row_nz_entry(3, 0), dense_values[3][m.row_nz_index(3, 0)]);

    m.row_nz_entry(1, 0) = 8.0;
    m.row_nz_entry(2, 2) = 9.0;

    EXPECT_EQ(m.rows(), 4);
    EXPECT_EQ(m.columns(), 3);
    EXPECT_EQ(m.non_zero_size(), 7);

    EXPECT_EQ(m.row_nz_size(0), 2);
    EXPECT_EQ(m.row_nz_size(1), 1);
    EXPECT_EQ(m.row_nz_size(2), 3);
    EXPECT_EQ(m.row_nz_size(3), 1);

    EXPECT_DOUBLE_EQ(m.row_nz_entry(0, 0), dense_values[0][m.row_nz_index(0, 0)]);
    EXPECT_DOUBLE_EQ(m.row_nz_entry(0, 1), dense_values[0][m.row_nz_index(0, 1)]);
    EXPECT_DOUBLE_EQ(m.row_nz_entry(1, 0), 8.0);
    EXPECT_DOUBLE_EQ(m.row_nz_entry(2, 0), dense_values[2][m.row_nz_index(2, 0)]);
    EXPECT_DOUBLE_EQ(m.row_nz_entry(2, 1), dense_values[2][m.row_nz_index(2, 1)]);
    EXPECT_DOUBLE_EQ(m.row_nz_entry(2, 2), 9.0);
    EXPECT_DOUBLE_EQ(m.row_nz_entry(3, 0), dense_values[3][m.row_nz_index(3, 0)]);
}

TEST(SparseMatrixCSR, lu_solver_3x3)
{
    using triplet = SparseMatrixCSR<double>::triplet_type;

    SparseMatrixCSR<double> matrix(3, 3,
                                   {triplet{0, 0, 2.0}, triplet{0, 2, -10.0}, triplet{1, 0, -1.0}, triplet{1, 1, 2.0},
                                    triplet{1, 2, -1.0}, triplet{2, 0, -7.0}, triplet{2, 2, 2.0}});

    Vector<double> rhs("rhs", 3);
    rhs(0) = 1.0;
    rhs(1) = -5;
    rhs(2) = 3.0;
    Vector<double> exact_solution("exact solution", 3);
    exact_solution(0) = -16.0 / 33.0;
    exact_solution(1) = -125.0 / 44.0;
    exact_solution(2) = -13.0 / 66.0;

    SparseLUSolver solver(matrix);
    solver.solveInPlace(rhs);

    EXPECT_DOUBLE_EQ(rhs(0), exact_solution(0));
    EXPECT_DOUBLE_EQ(rhs(1), exact_solution(1));
    EXPECT_DOUBLE_EQ(rhs(2), exact_solution(2));
}

TEST(SparseMatrixCSR, lu_solver_5x5)
{
    using triplet = SparseMatrixCSR<double>::triplet_type;

    SparseMatrixCSR<double> matrix(5, 5,
                                   {triplet{0, 0, 2.0}, triplet{0, 3, 7.0}, triplet{0, 4, 3.0}, triplet{1, 0, -1.0},
                                    triplet{1, 1, -5.0}, triplet{1, 2, -1.0}, triplet{1, 3, 5.0}, triplet{1, 4, 6.0},
                                    triplet{2, 1, -7.0}, triplet{2, 2, 1.0}, triplet{2, 3, 2.0}, triplet{2, 4, -4.0},
                                    triplet{3, 2, -4.0}, triplet{3, 3, 6.0}, triplet{3, 4, 9.0}, triplet{4, 0, 2.0},
                                    triplet{4, 1, 4.0}, triplet{4, 3, -4.0}, triplet{4, 4, 9.0}});

    Vector<double> rhs("rhs", 5);
    rhs(0) = 1.0;
    rhs(1) = -5;
    rhs(2) = 3.0;
    rhs(3) = 7.0;
    rhs(4) = 2.0;
    Vector<double> exact_solution("exact_solution", 5);
    exact_solution(0) = 2792.0 / 567.0;
    exact_solution(1) = -589.0 / 648.0;
    exact_solution(2) = -7615.0 / 1512.0;
    exact_solution(3) = -1013.0 / 1134.0;
    exact_solution(4) = -109.0 / 126.0;

    SparseLUSolver solver(matrix);
    solver.solveInPlace(rhs);

    EXPECT_DOUBLE_EQ(rhs(0), exact_solution(0));
    EXPECT_DOUBLE_EQ(rhs(1), exact_solution(1));
    EXPECT_DOUBLE_EQ(rhs(2), exact_solution(2));
    EXPECT_DOUBLE_EQ(rhs(3), exact_solution(3));
    EXPECT_DOUBLE_EQ(rhs(4), exact_solution(4));
}
