#include <gtest/gtest.h>

#include "../../include/LinearAlgebra/vector.h"
#include "../../include/LinearAlgebra/vector_operations.h"

/* bool equals(const Vector<T>& lhs, const Vector<T>& rhs); */

TEST(VectorOperations, equals_vector_vector)
{
    const Vector<double> v  = {1, 2, 3};
    const Vector<double> v1 = {1, 2, 3};
    const Vector<double> v2 = {1, 2, 3, 4};
    const Vector<double> v3 = {1, 2};
    const Vector<double> v4 = {1, 3, 3};

    EXPECT_TRUE(equals(v, v));
    EXPECT_TRUE(equals(v, v1));

    EXPECT_FALSE(equals(v, v2));
    EXPECT_FALSE(equals(v, v3));
    EXPECT_FALSE(equals(v, v4));
}

/* void assign(Vector<T>& lhs, const T& value); */

TEST(VectorOperations, assign_vector_scalar)
{
    Vector<double> v = {1, 2, 3};
    assign(v, 5.0);
    const Vector<double> expected_result = {5, 5, 5};
    EXPECT_TRUE(equals(v, expected_result));
}

/* void add(Vector<T>& result, const Vector<T>& x); */

TEST(VectorOperations, add_vector_vector)
{
    Vector<double> v1       = {1, 2, 3};
    const Vector<double> v2 = {-1, -5, 2};
    add(v1, v2);
    const Vector<double> expected_result = {0, -3, 5};
    EXPECT_TRUE(equals(v1, expected_result));
}

/* void subtract(Vector<T>& result, const Vector<T>& x) */

TEST(VectorOperations, subtract_vector_vector)
{
    Vector<double> v1       = {1, 2, 3};
    const Vector<double> v2 = {-1, -5, 2};
    subtract(v1, v2);
    const Vector<double> expected_result = {2, 7, 1};
    EXPECT_TRUE(equals(v1, expected_result));
}

/* void linear_combination(Vector<T>& x, const T& alpha, const Vector<T>& y, const T& beta); */

TEST(VectorOperations, linear_combination)
{
    Vector<double> v1       = {1, 2, 3};
    const Vector<double> v2 = {-1, -5, 2};
    const double alpha      = -3.0;
    const double beta       = 2.0;
    linear_combination(v1, alpha, v2, beta);
    const Vector<double> expected_result = {-5.0, -16.0, -5.0};
    EXPECT_TRUE(equals(v1, expected_result));
}

/* void multiply(Vector<T>& x, const T& alpha); */

TEST(VectorOperations, multiply_vector_scalar)
{
    Vector<double> v1  = {1, -2, 3};
    const double alpha = -3.0;
    multiply(v1, alpha);
    const Vector<double> expected_result = {-3.0, 6.0, -9.0};
    EXPECT_TRUE(equals(v1, expected_result));
}

/* T dot_product(const Vector<T>& lhs, const Vector<T>& rhs); */

TEST(VectorOperations, dot_product)
{
    const Vector<double> v1 = {1, 2, 3};
    const Vector<double> v2 = {1, 5, 2};
    EXPECT_DOUBLE_EQ(dot_product(v1, v2), 17.0);
}

/* T l1_norm(const Vector<T>& x); */

TEST(VectorOperations, l1_vector_norm)
{
    const Vector<double> v = {1, -5, 2};
    EXPECT_DOUBLE_EQ(l1_norm(v), 8.0);
}

/* T l2_norm_squared(const Vector<T>& x); */

TEST(VectorOperations, l2_vector_norm_squared)
{
    const Vector<double> v = {1, -5, 2};
    EXPECT_DOUBLE_EQ(l2_norm_squared(v), 30.0);
}

/* T infinity_norm(const Vector<T>& x); */

TEST(VectorOperations, infinity_vector_norm)
{
    const Vector<double> v = {1, -5, 2};
    EXPECT_DOUBLE_EQ(infinity_norm(v), 5.0);
}
