#include <random>

#include <gtest/gtest.h>

#include "../../include/LinearAlgebra/Vector/vector.h"
#include "../../include/LinearAlgebra/Vector/vector_operations.h"
#include "../../include/LinearAlgebra/Solvers/symmetricTridiagonalSolver.h"

TEST(SymmetricTridiagonalSolver, diagonal_dominant_n_2)
{
    const int n = 2;

    SymmetricTridiagonalSolver<double> solver(n);
    solver.is_cyclic(false);
    solver.main_diagonal(0) = -4.0;
    solver.main_diagonal(1) = 8.0;
    solver.sub_diagonal(0)  = 1.0;

    Vector<double> rhs("rhs", n);
    rhs[0] = 3.0;
    rhs[1] = -2.0;

    Vector<double> exact_solution("exact_solution", n);
    exact_solution[0] = -26.0 / 33.0;
    exact_solution[1] = -5.0 / 33.0;

    Vector<double> temp("temp", n);
    solver.solveInPlace(rhs.data(), temp.data());

    EXPECT_NEAR(rhs[0], exact_solution[0], 1e-12);
    EXPECT_NEAR(rhs[1], exact_solution[1], 1e-12);

    EXPECT_TRUE(equals(ConstVector<double>(rhs), ConstVector<double>(exact_solution)));
}

TEST(SymmetricTridiagonalSolver, not_diagonal_dominant_n_2)
{
    const int n = 2;

    SymmetricTridiagonalSolver<double> solver(n);
    solver.is_cyclic(false);
    solver.main_diagonal(0) = -4.0;
    solver.main_diagonal(1) = 8.0;
    solver.sub_diagonal(0)  = 1000.0;

    Vector<double> rhs("rhs", n);
    rhs[0] = 3.0;
    rhs[1] = -2.0;

    Vector<double> exact_solution("exact_solution", n);
    exact_solution[0] = -23.0 / 11364.0;
    exact_solution[1] = 17.0 / 5682.0;

    Vector<double> temp("temp", n);
    solver.solveInPlace(rhs.data(), temp.data());

    EXPECT_NEAR(rhs[0], exact_solution[0], 1e-12);
    EXPECT_NEAR(rhs[1], exact_solution[1], 1e-12);

    EXPECT_TRUE(equals(ConstVector<double>(rhs), ConstVector<double>(exact_solution)));
}

TEST(SymmetricTridiagonalSolver, diagonal_dominant_n_3)
{
    const int n = 3;

    SymmetricTridiagonalSolver<double> solver(n);
    solver.is_cyclic(false);
    solver.main_diagonal(0) = -5.0;
    solver.main_diagonal(1) = 10.0;
    solver.main_diagonal(2) = 12.0;
    solver.sub_diagonal(0)  = 1.0;
    solver.sub_diagonal(1)  = -9.0;

    Vector<double> rhs("rhs", n);
    rhs[0] = 65.0;
    rhs[1] = -10.5;
    rhs[2] = 13.6;

    Vector<double> exact_solution("exact_solution", n);
    exact_solution[0] = -4231.0 / 345.0;
    exact_solution[1] = 254.0 / 69.0;
    exact_solution[2] = 2687.0 / 690.0;

    Vector<double> temp("temp", n);
    solver.solveInPlace(rhs.data(), temp.data());

    EXPECT_NEAR(rhs[0], exact_solution[0], 1e-12);
    EXPECT_NEAR(rhs[1], exact_solution[1], 1e-12);
    EXPECT_NEAR(rhs[2], exact_solution[2], 1e-12);

    EXPECT_TRUE(equals(ConstVector<double>(rhs), ConstVector<double>(exact_solution)));
}

TEST(SymmetricTridiagonalSolver, not_diagonal_dominant_n_3)
{
    const int n = 3;

    SymmetricTridiagonalSolver<double> solver(n);
    solver.is_cyclic(false);
    solver.main_diagonal(0) = -1.0;
    solver.main_diagonal(1) = 0.001;
    solver.main_diagonal(2) = 0.6;
    solver.sub_diagonal(0)  = -20.0;
    solver.sub_diagonal(1)  = -9.0;

    Vector<double> rhs("rhs", n);
    rhs[0] = 65.0;
    rhs[1] = -10.5;
    rhs[2] = 13.6;

    Vector<double> exact_solution("exact_solution", n);
    exact_solution[0] = 4904935.0 / 265001.0;
    exact_solution[1] = -1106500.0 / 265001.0;
    exact_solution[2] = -31772432.0 / 795003.0;

    Vector<double> temp("temp", n);
    solver.solveInPlace(rhs.data(), temp.data());

    EXPECT_NEAR(rhs[0], exact_solution[0], 1e-12);
    EXPECT_NEAR(rhs[1], exact_solution[1], 1e-12);
    EXPECT_NEAR(rhs[2], exact_solution[2], 1e-12);

    EXPECT_DOUBLE_EQ(rhs[0], exact_solution[0]);

    EXPECT_TRUE(equals(ConstVector<double>(rhs), ConstVector<double>(exact_solution)));
}

TEST(SymmetricTridiagonalSolver, random_tridiagonal_n_10)
{
    const int n            = 10;
    const double precision = 1e-9;

    SymmetricTridiagonalSolver<double> solver(n);
    solver.is_cyclic(false);

    // Create random number generator
    std::random_device rd;
    std::mt19937 gen(42);
    std::uniform_real_distribution<> dis(-100.0, 100.0);

    // Fill main_diagonal with random values
    for (int i = 0; i < n; ++i) {
        solver.main_diagonal(i) = dis(gen);
    }

    // Fill sub_diagonal with random values
    for (int i = 0; i < n - 1; ++i) {
        solver.sub_diagonal(i) = dis(gen);
    }

    const SymmetricTridiagonalSolver<double> copy_solver = solver;

    // Fill rhs with random values
    Vector<double> rhs("rhs", n);
    for (int i = 0; i < n; ++i) {
        rhs[i] = dis(gen);
    }

    Vector<double> copy_rhs("copy_rhs", rhs.size());
    Kokkos::deep_copy(copy_rhs, rhs);

    Vector<double> temp("temp", n);
    solver.solveInPlace(rhs.data(), temp.data());

    EXPECT_NEAR(copy_solver.main_diagonal(0) * rhs[0] + copy_solver.sub_diagonal(0) * rhs[1], copy_rhs[0], precision);
    for (int i = 1; i < n - 1; ++i) {
        EXPECT_NEAR(copy_solver.sub_diagonal(i - 1) * rhs[i - 1] + copy_solver.main_diagonal(i) * rhs[i] +
                        copy_solver.sub_diagonal(i) * rhs[i + 1],
                    copy_rhs[i], precision);
    }
    EXPECT_NEAR(copy_solver.sub_diagonal(n - 2) * rhs[n - 2] + copy_solver.main_diagonal(n - 1) * rhs[n - 1],
                copy_rhs[n - 1], precision);
}

TEST(SymmetricTridiagonalSolver, random_tridiagonal_n_100)
{
    const int n = 100;

    const double precision = 1e-8;

    SymmetricTridiagonalSolver<double> solver(n);
    solver.is_cyclic(false);

    // Create random number generator
    std::random_device rd;
    std::mt19937 gen(42);
    std::uniform_real_distribution<> dis(-100.0, 100.0);

    // Fill main_diagonal with random values
    for (int i = 0; i < n; ++i) {
        solver.main_diagonal(i) = dis(gen);
    }

    // Fill sub_diagonal with random values
    for (int i = 0; i < n - 1; ++i) {
        solver.sub_diagonal(i) = dis(gen);
    }

    const SymmetricTridiagonalSolver<double> copy_solver = solver;

    // Fill rhs with random values
    Vector<double> rhs("rhs", n);
    for (int i = 0; i < n; ++i) {
        rhs[i] = dis(gen);
    }

    Vector<double> copy_rhs("copy_rhs", rhs.size());
    Kokkos::deep_copy(copy_rhs, rhs);

    Vector<double> temp("temp", n);
    solver.solveInPlace(rhs.data(), temp.data());

    EXPECT_NEAR(copy_solver.main_diagonal(0) * rhs[0] + copy_solver.sub_diagonal(0) * rhs[1], copy_rhs[0], precision);
    for (int i = 1; i < n - 1; ++i) {
        EXPECT_NEAR(copy_solver.sub_diagonal(i - 1) * rhs[i - 1] + copy_solver.main_diagonal(i) * rhs[i] +
                        copy_solver.sub_diagonal(i) * rhs[i + 1],
                    copy_rhs[i], precision);
    }
    EXPECT_NEAR(copy_solver.sub_diagonal(n - 2) * rhs[n - 2] + copy_solver.main_diagonal(n - 1) * rhs[n - 1],
                copy_rhs[n - 1], precision);
}

TEST(SymmetricTridiagonalSolver, random_tridiagonal_n_1000)
{
    const int n = 1000;

    const double precision = 1e-8;

    SymmetricTridiagonalSolver<double> solver(n);
    solver.is_cyclic(false);

    // Create random number generator
    std::random_device rd;
    std::mt19937 gen(42);
    std::uniform_real_distribution<> dis(-100.0, 100.0);

    // Fill main_diagonal with random values
    for (int i = 0; i < n; ++i) {
        solver.main_diagonal(i) = dis(gen);
    }

    // Fill sub_diagonal with random values
    for (int i = 0; i < n - 1; ++i) {
        solver.sub_diagonal(i) = dis(gen);
    }

    const SymmetricTridiagonalSolver<double> copy_solver = solver;

    // Fill rhs with random values
    Vector<double> rhs("rhs", n);
    for (int i = 0; i < n; ++i) {
        rhs[i] = dis(gen);
    }

    Vector<double> copy_rhs("copy_rhs", rhs.size());
    Kokkos::deep_copy(copy_rhs, rhs);

    Vector<double> temp("temp", n);
    solver.solveInPlace(rhs.data(), temp.data());

    EXPECT_NEAR(copy_solver.main_diagonal(0) * rhs[0] + copy_solver.sub_diagonal(0) * rhs[1], copy_rhs[0], precision);
    for (int i = 1; i < n - 1; ++i) {
        EXPECT_NEAR(copy_solver.sub_diagonal(i - 1) * rhs[i - 1] + copy_solver.main_diagonal(i) * rhs[i] +
                        copy_solver.sub_diagonal(i) * rhs[i + 1],
                    copy_rhs[i], precision);
    }
    EXPECT_NEAR(copy_solver.sub_diagonal(n - 2) * rhs[n - 2] + copy_solver.main_diagonal(n - 1) * rhs[n - 1],
                copy_rhs[n - 1], precision);
}

TEST(SymmetricTridiagonalSolver, random_tridiagonal_n_10000)
{
    const int n = 10000;

    const double precision = 1e-7;

    SymmetricTridiagonalSolver<double> solver(n);
    solver.is_cyclic(false);

    // Create random number generator
    std::random_device rd;
    std::mt19937 gen(42);
    std::uniform_real_distribution<> dis(-100.0, 100.0);

    // Fill main_diagonal with random values
    for (int i = 0; i < n; ++i) {
        solver.main_diagonal(i) = dis(gen);
    }

    // Fill sub_diagonal with random values
    for (int i = 0; i < n - 1; ++i) {
        solver.sub_diagonal(i) = dis(gen);
    }

    const SymmetricTridiagonalSolver<double> copy_solver = solver;

    // Fill rhs with random values
    Vector<double> rhs("rhs", n);
    for (int i = 0; i < n; ++i) {
        rhs[i] = dis(gen);
    }

    Vector<double> copy_rhs("copy_rhs", rhs.size());
    Kokkos::deep_copy(copy_rhs, rhs);

    Vector<double> temp("temp", n);
    solver.solveInPlace(rhs.data(), temp.data());

    EXPECT_NEAR(copy_solver.main_diagonal(0) * rhs[0] + copy_solver.sub_diagonal(0) * rhs[1], copy_rhs[0], precision);
    for (int i = 1; i < n - 1; ++i) {
        EXPECT_NEAR(copy_solver.sub_diagonal(i - 1) * rhs[i - 1] + copy_solver.main_diagonal(i) * rhs[i] +
                        copy_solver.sub_diagonal(i) * rhs[i + 1],
                    copy_rhs[i], precision);
    }
    EXPECT_NEAR(copy_solver.sub_diagonal(n - 2) * rhs[n - 2] + copy_solver.main_diagonal(n - 1) * rhs[n - 1],
                copy_rhs[n - 1], precision);
}

TEST(SymmetricTridiagonalSolver, random_tridiagonal_boosted_subdiagonal_LOW_PRECISION_n_10000)
{
    const int n = 10000;

    const double precision = 1e-3;

    const double subdiagonal_boost = 10000.0;

    SymmetricTridiagonalSolver<double> solver(n);
    solver.is_cyclic(false);

    // Create random number generator
    std::random_device rd;
    std::mt19937 gen(79);
    std::uniform_real_distribution<> dis(-100.0, 100.0);

    // Fill main_diagonal with random values
    for (int i = 0; i < n; ++i) {
        solver.main_diagonal(i) = dis(gen);
    }

    // Fill sub_diagonal with random values
    for (int i = 0; i < n - 1; ++i) {
        solver.sub_diagonal(i) = subdiagonal_boost * dis(gen);
    }

    const SymmetricTridiagonalSolver<double> copy_solver = solver;

    // Fill rhs with random values
    Vector<double> rhs("rhs", n);
    for (int i = 0; i < n; ++i) {
        rhs[i] = dis(gen);
    }

    Vector<double> copy_rhs("copy_rhs", rhs.size());
    Kokkos::deep_copy(copy_rhs, rhs);

    Vector<double> temp("temp", n);
    solver.solveInPlace(rhs.data(), temp.data());

    EXPECT_NEAR(copy_solver.main_diagonal(0) * rhs[0] + copy_solver.sub_diagonal(0) * rhs[1], copy_rhs[0], precision);
    for (int i = 1; i < n - 1; ++i) {
        EXPECT_NEAR(copy_solver.sub_diagonal(i - 1) * rhs[i - 1] + copy_solver.main_diagonal(i) * rhs[i] +
                        copy_solver.sub_diagonal(i) * rhs[i + 1],
                    copy_rhs[i], precision);
    }
    EXPECT_NEAR(copy_solver.sub_diagonal(n - 2) * rhs[n - 2] + copy_solver.main_diagonal(n - 1) * rhs[n - 1],
                copy_rhs[n - 1], precision);
}
