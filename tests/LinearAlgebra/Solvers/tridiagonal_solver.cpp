#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <vector>
#include <cmath>

#include "../../../include/LinearAlgebra/Solvers/tridiagonal_solver.h"
#include "../../../include/LinearAlgebra/Vector/vector.h"

// clang-format off
TEST(BatchedTridiagonalSolvers, non_cyclic_tridiagonal_n_4)
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

    solver.main_diagonal(0,0) = 2.0; solver.sub_diagonal(0,0) = 1.0;
    solver.main_diagonal(0,1) = 4.0; solver.sub_diagonal(0,1) = 2.0;
    solver.main_diagonal(0,2) = 6.0; solver.sub_diagonal(0,2) = 3.0;
    solver.main_diagonal(0,3) = 8.0;

    solver.main_diagonal(1,0) = 3.0; solver.sub_diagonal(1,0) = 1.0;
    solver.main_diagonal(1,1) = 5.0; solver.sub_diagonal(1,1) = 2.0;
    solver.main_diagonal(1,2) = 7.0; solver.sub_diagonal(1,2) = 3.0;
    solver.main_diagonal(1,3) = 9.0;    

    solver.main_diagonal(2,0) = 4.0; solver.sub_diagonal(2,0) = 1.0;
    solver.main_diagonal(2,1) = 6.0; solver.sub_diagonal(2,1) = 2.0;
    solver.main_diagonal(2,2) = 8.0; solver.sub_diagonal(2,2) = 3.0;
    solver.main_diagonal(2,3) = 10.0;

    solver.main_diagonal(3,0) = 5.0; solver.sub_diagonal(3,0) = 1.0;
    solver.main_diagonal(3,1) = 7.0; solver.sub_diagonal(3,1) = 2.0;
    solver.main_diagonal(3,2) = 9.0; solver.sub_diagonal(3,2) = 3.0;
    solver.main_diagonal(3,3) = 11.0;

    Vector<double> rhs("rhs", matrix_dimension * batch_count);

    // Initialize RHS for each system
    rhs(0) = 1.0; rhs(1) = 2.0; rhs(2) = 3.0; rhs(3) = 4.0;
    rhs(4) = 2.0; rhs(5) = 3.0; rhs(6) = 4.0; rhs(7) = 5.0;
    rhs(8) = 3.0; rhs(9) = 4.0; rhs(10) = 5.0; rhs(11) = 6.0;
    rhs(12) = 4.0; rhs(13) = 5.0; rhs(14) = 6.0; rhs(15) = 7.0; 

    solver.setup();

    int offset, stride;
    // Solve each even system separately
    offset = 0; stride = 2;
    solver.solve(rhs, offset, stride);
    // Solve each odd system separately
    offset = 1; stride = 2;
    solver.solve(rhs, offset, stride);

    // Verify solutions
    double tol = 1e-12;

    EXPECT_NEAR(rhs(0), 70.0/209.0, tol);
    EXPECT_NEAR(rhs(1), 69.0/209.0, tol);
    EXPECT_NEAR(rhs(2), 36.0/209.0, tol);
    EXPECT_NEAR(rhs(3), 91.0/209.0, tol);

    EXPECT_NEAR(rhs(4), 29.0/54.0, tol);
    EXPECT_NEAR(rhs(5), 7.0/18.0, tol);
    EXPECT_NEAR(rhs(6), 7.0/27.0, tol);
    EXPECT_NEAR(rhs(7), 38.0/81.0, tol);

    EXPECT_NEAR(rhs(8), 938.0/1473.0, tol);
    EXPECT_NEAR(rhs(9), 667.0/1473.0, tol);
    EXPECT_NEAR(rhs(10), 476.0/1473.0, tol);
    EXPECT_NEAR(rhs(11), 247.0/491.0, tol);

    EXPECT_NEAR(rhs(12), 248.0/355.0, tol);
    EXPECT_NEAR(rhs(13), 36.0/71.0, tol);
    EXPECT_NEAR(rhs(14), 267.0/710.0, tol);
    EXPECT_NEAR(rhs(15), 379.0/710.0, tol);
}

TEST(BatchedTridiagonalSolvers, cyclic_tridiagonal_n_4)
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

    solver.main_diagonal(0,0) = 2.0; solver.sub_diagonal(0,0) = 1.0;
    solver.main_diagonal(0,1) = 4.0; solver.sub_diagonal(0,1) = 2.0;
    solver.main_diagonal(0,2) = 6.0; solver.sub_diagonal(0,2) = 3.0;
    solver.main_diagonal(0,3) = 8.0; solver.cyclic_corner(0) = -1.0;

    solver.main_diagonal(1,0) = 3.0; solver.sub_diagonal(1,0) = 1.0;
    solver.main_diagonal(1,1) = 5.0; solver.sub_diagonal(1,1) = 2.0;
    solver.main_diagonal(1,2) = 7.0; solver.sub_diagonal(1,2) = 3.0;
    solver.main_diagonal(1,3) = 9.0; solver.cyclic_corner(1) = -2.0;

    solver.main_diagonal(2,0) = 4.0; solver.sub_diagonal(2,0) = 1.0;
    solver.main_diagonal(2,1) = 6.0; solver.sub_diagonal(2,1) = 2.0;
    solver.main_diagonal(2,2) = 8.0; solver.sub_diagonal(2,2) = 3.0;
    solver.main_diagonal(2,3) = 10.0; solver.cyclic_corner(2) = -3.0;

    solver.main_diagonal(3,0) = 5.0; solver.sub_diagonal(3,0) = 1.0;
    solver.main_diagonal(3,1) = 7.0; solver.sub_diagonal(3,1) = 2.0;
    solver.main_diagonal(3,2) = 9.0; solver.sub_diagonal(3,2) = 3.0;
    solver.main_diagonal(3,3) = 11.0; solver.cyclic_corner(3) = -4.0;

    Vector<double> rhs("rhs", matrix_dimension * batch_count);

    // Initialize RHS for each system
    rhs(0) = 1.0; rhs(1) = 2.0; rhs(2) = 3.0; rhs(3) = 4.0;
    rhs(4) = 2.0; rhs(5) = 3.0; rhs(6) = 4.0; rhs(7) = 5.0;
    rhs(8) = 3.0; rhs(9) = 4.0; rhs(10) = 5.0; rhs(11) = 6.0;
    rhs(12) = 4.0; rhs(13) = 5.0; rhs(14) = 6.0; rhs(15) = 7.0; 

    solver.setup();

    int offset, stride;
    // Solve each even system separately
    offset = 0; stride = 2;
    solver.solve(rhs, offset, stride);
    // Solve each odd system separately
    offset = 1; stride = 2;
    solver.solve(rhs, offset, stride);

    // Verify solutions
    double tol = 1e-12;

    EXPECT_NEAR(rhs(0), 42.0/67.0, tol);
    EXPECT_NEAR(rhs(1), 18.0/67.0, tol);
    EXPECT_NEAR(rhs(2), 10.0/67.0, tol);
    EXPECT_NEAR(rhs(3), 35.0/67.0, tol);

    EXPECT_NEAR(rhs(4), 287.0/274.0, tol);
    EXPECT_NEAR(rhs(5), 89.0/274.0, tol);
    EXPECT_NEAR(rhs(6), 45.0/274.0, tol);
    EXPECT_NEAR(rhs(7), 201.0/274.0, tol);

    EXPECT_NEAR(rhs(8), 1532.0/1113.0, tol);
    EXPECT_NEAR(rhs(9), 8.0/21.0, tol);
    EXPECT_NEAR(rhs(10), 188.0/1113.0, tol);
    EXPECT_NEAR(rhs(11), 51.0/53.0, tol);

    EXPECT_NEAR(rhs(12), 271.0/162.0, tol);
    EXPECT_NEAR(rhs(13), 23.0/54.0, tol);
    EXPECT_NEAR(rhs(14), 14.0/81.0, tol);
    EXPECT_NEAR(rhs(15), 97.0/81.0, tol);
}

TEST(BatchedTridiagonalSolvers, non_cyclic_diagonal_n_4)
{
    int batch_count = 4;
    int matrix_dimension = 4;
    bool is_cyclic = false;

    BatchedTridiagonalSolver<double> solver(matrix_dimension, batch_count, is_cyclic);

    solver.main_diagonal(0,0) = 2.0;
    solver.main_diagonal(0,1) = 4.0;
    solver.main_diagonal(0,2) = 6.0;
    solver.main_diagonal(0,3) = 8.0;

    solver.main_diagonal(1,0) = 3.0;
    solver.main_diagonal(1,1) = 5.0;
    solver.main_diagonal(1,2) = 7.0;
    solver.main_diagonal(1,3) = 9.0;

    solver.main_diagonal(2,0) = 4.0;
    solver.main_diagonal(2,1) = 6.0;
    solver.main_diagonal(2,2) = 8.0;
    solver.main_diagonal(2,3) = 10.0;

    solver.main_diagonal(3,0) = 5.0;
    solver.main_diagonal(3,1) = 7.0;
    solver.main_diagonal(3,2) = 9.0;
    solver.main_diagonal(3,3) = 11.0;

    Vector<double> rhs("rhs", matrix_dimension * batch_count);

    // Initialize RHS for each system
    rhs(0) = 1.0; rhs(1) = 2.0; rhs(2) = 3.0; rhs(3) = 4.0;
    rhs(4) = 2.0; rhs(5) = 3.0; rhs(6) = 4.0; rhs(7) = 5.0;
    rhs(8) = 3.0; rhs(9) = 4.0; rhs(10) = 5.0; rhs(11) = 6.0;
    rhs(12) = 4.0; rhs(13) = 5.0; rhs(14) = 6.0; rhs(15) = 7.0; 

    solver.setup();

    int offset, stride;
    // Solve each even system separately
    offset = 0; stride = 2;
    solver.solve_diagonal(rhs, offset, stride);
    // Solve each odd system separately
    offset = 1; stride = 2;
    solver.solve_diagonal(rhs, offset, stride);

    // Verify solutions
    double tol = 1e-12;

    EXPECT_NEAR(rhs(0), 1.0/2.0, tol);
    EXPECT_NEAR(rhs(1), 2.0/4.0, tol);
    EXPECT_NEAR(rhs(2), 3.0/6.0, tol);
    EXPECT_NEAR(rhs(3), 4.0/8.0, tol);

    EXPECT_NEAR(rhs(4), 2.0/3.0, tol);
    EXPECT_NEAR(rhs(5), 3.0/5.0, tol);
    EXPECT_NEAR(rhs(6), 4.0/7.0, tol);
    EXPECT_NEAR(rhs(7), 5.0/9.0, tol);

    EXPECT_NEAR(rhs(8), 3.0/4.0, tol);
    EXPECT_NEAR(rhs(9), 4.0/6.0, tol);
    EXPECT_NEAR(rhs(10), 5.0/8.0, tol);
    EXPECT_NEAR(rhs(11), 6.0/10.0, tol);

    EXPECT_NEAR(rhs(12), 4.0/5.0, tol);
    EXPECT_NEAR(rhs(13), 5.0/7.0, tol);
    EXPECT_NEAR(rhs(14), 6.0/9.0, tol);
    EXPECT_NEAR(rhs(15), 7.0/11.0, tol);
}

TEST(BatchedTridiagonalSolvers, cyclic_diagonal_n_4)
{
    int batch_count = 4;
    int matrix_dimension = 4;
    bool is_cyclic = true;

    BatchedTridiagonalSolver<double> solver(matrix_dimension, batch_count, is_cyclic);

    solver.main_diagonal(0,0) = 2.0;
    solver.main_diagonal(0,1) = 4.0;
    solver.main_diagonal(0,2) = 6.0;
    solver.main_diagonal(0,3) = 8.0;

    solver.main_diagonal(1,0) = 3.0;
    solver.main_diagonal(1,1) = 5.0;
    solver.main_diagonal(1,2) = 7.0;
    solver.main_diagonal(1,3) = 9.0;

    solver.main_diagonal(2,0) = 4.0;
    solver.main_diagonal(2,1) = 6.0;
    solver.main_diagonal(2,2) = 8.0;
    solver.main_diagonal(2,3) = 10.0;

    solver.main_diagonal(3,0) = 5.0;
    solver.main_diagonal(3,1) = 7.0;
    solver.main_diagonal(3,2) = 9.0;
    solver.main_diagonal(3,3) = 11.0;

    Vector<double> rhs("rhs", matrix_dimension * batch_count);

    // Initialize RHS for each system
    rhs(0) = 1.0; rhs(1) = 2.0; rhs(2) = 3.0; rhs(3) = 4.0;
    rhs(4) = 2.0; rhs(5) = 3.0; rhs(6) = 4.0; rhs(7) = 5.0;
    rhs(8) = 3.0; rhs(9) = 4.0; rhs(10) = 5.0; rhs(11) = 6.0;
    rhs(12) = 4.0; rhs(13) = 5.0; rhs(14) = 6.0; rhs(15) = 7.0; 

    solver.setup();

    int offset, stride;
    // Solve each even system separately
    offset = 0; stride = 2;
    solver.solve_diagonal(rhs, offset, stride);
    // Solve each odd system separately
    offset = 1; stride = 2;
    solver.solve_diagonal(rhs, offset, stride);

    // Verify solutions
    double tol = 1e-12;

    EXPECT_NEAR(rhs(0), 1.0/2.0, tol);
    EXPECT_NEAR(rhs(1), 2.0/4.0, tol);
    EXPECT_NEAR(rhs(2), 3.0/6.0, tol);
    EXPECT_NEAR(rhs(3), 4.0/8.0, tol);

    EXPECT_NEAR(rhs(4), 2.0/3.0, tol);
    EXPECT_NEAR(rhs(5), 3.0/5.0, tol);
    EXPECT_NEAR(rhs(6), 4.0/7.0, tol);
    EXPECT_NEAR(rhs(7), 5.0/9.0, tol);

    EXPECT_NEAR(rhs(8), 3.0/4.0, tol);
    EXPECT_NEAR(rhs(9), 4.0/6.0, tol);
    EXPECT_NEAR(rhs(10), 5.0/8.0, tol);
    EXPECT_NEAR(rhs(11), 6.0/10.0, tol);

    EXPECT_NEAR(rhs(12), 4.0/5.0, tol);
    EXPECT_NEAR(rhs(13), 5.0/7.0, tol);
    EXPECT_NEAR(rhs(14), 6.0/9.0, tol);
    EXPECT_NEAR(rhs(15), 7.0/11.0, tol);
}
