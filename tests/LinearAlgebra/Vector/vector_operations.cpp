#include <gtest/gtest.h>

#include "../../../include/LinearAlgebra/Vector/vector.h"
#include "../../../include/LinearAlgebra/Vector/vector_operations.h"
using namespace gmgpolar;

/* bool equals(HostConstVector<T>& lhs, HostConstVector<T>& rhs); */

TEST(VectorOperations, equals_vector_vector)
{
    HostVector<double> v("v", 3);
    HostVector<double> v1("v1", 3);
    HostVector<double> v2("v2", 4);
    HostVector<double> v3("v3", 2);
    HostVector<double> v4("v4", 3);
    v(0) = 1;
    v(1) = 2;
    v(2) = 3;

    v1(0) = 1;
    v1(1) = 2;
    v1(2) = 3;

    v2(0) = 1;
    v2(1) = 2;
    v2(2) = 3;
    v2(3) = 4;

    v3(0) = 1;
    v3(1) = 2;

    v4(0) = 1;
    v4(1) = 3;
    v4(2) = 3;

    EXPECT_TRUE(equals(HostConstVector<double>(v), HostConstVector<double>(v)));
    EXPECT_TRUE(equals(HostConstVector<double>(v), HostConstVector<double>(v1)));
    EXPECT_FALSE(equals(HostConstVector<double>(v), HostConstVector<double>(v2)));
    EXPECT_FALSE(equals(HostConstVector<double>(v), HostConstVector<double>(v3)));
    EXPECT_FALSE(equals(HostConstVector<double>(v), HostConstVector<double>(v4)));
}

/* void assign(HostVector<T>& lhs, const T& value); */

TEST(VectorOperations, assign_vector_scalar)
{
    HostVector<double> v("v", 3);
    v(0) = 1;
    v(1) = 2;
    v(2) = 3;
    assign(v, 5.0);
    HostVector<double> expected_result("expected_result", 3);
    Kokkos::deep_copy(expected_result, 5);
    EXPECT_TRUE(equals(HostConstVector<double>(v), HostConstVector<double>(expected_result)));
}

/* void add(HostVector<T> result, HostConstVector<T> x); */

TEST(VectorOperations, add_vector_vector)
{
    HostVector<double> v1("v1", 3);
    HostVector<double> v2("v2", 3);
    v1(0) = 1;
    v1(1) = 2;
    v1(2) = 3;

    v2(0) = -1;
    v2(1) = -5;
    v2(2) = 2;
    add(v1, HostConstVector<double>(v2));

    HostVector<double> expected_result("expected_result", 3);
    expected_result(0) = 0;
    expected_result(1) = -3;
    expected_result(2) = 5;
    EXPECT_TRUE(equals(HostConstVector<double>(v1), HostConstVector<double>(expected_result)));
}

/* void subtract(HostVector<T>& result, HostConstVector<T>& x) */

TEST(VectorOperations, subtract_vector_vector)
{
    HostVector<double> v1("v1", 3);
    HostVector<double> v2("v2", 3);
    v1(0) = 1;
    v1(1) = 2;
    v1(2) = 3;

    v2(0) = -1;
    v2(1) = -5;
    v2(2) = 2;
    subtract(v1, HostConstVector<double>(v2));

    HostVector<double> expected_result("expected_result", 3);
    expected_result(0) = 2;
    expected_result(1) = 7;
    expected_result(2) = 1;
    EXPECT_TRUE(equals(HostConstVector<double>(v1), HostConstVector<double>(expected_result)));
}

/* void linear_combination(HostVector<T>& x, const T& alpha, HostConstVector<T>& y, const T& beta); */

TEST(VectorOperations, linear_combination)
{
    HostVector<double> v1("v1", 3);
    HostVector<double> v2("v2", 3);
    v1(0) = 1;
    v1(1) = 2;
    v1(2) = 3;

    v2(0)              = -1;
    v2(1)              = -5;
    v2(2)              = 2;
    const double alpha = -3.0;
    const double beta  = 2.0;
    linear_combination(v1, alpha, HostConstVector<double>(v2), beta);
    HostVector<double> expected_result("expected_result", 3);
    expected_result(0) = -5.;
    expected_result(1) = -16.;
    expected_result(2) = -5.;
    EXPECT_TRUE(equals(HostConstVector<double>(v1), HostConstVector<double>(expected_result)));
}

/* void multiply(HostVector<T>& x, const T& alpha); */

TEST(VectorOperations, multiply_vector_scalar)
{
    HostVector<double> v1("v1", 3);
    v1(0)              = 1;
    v1(1)              = -2;
    v1(2)              = 3;
    const double alpha = -3.0;
    multiply(v1, alpha);
    HostVector<double> expected_result("expected_result", 3);
    expected_result(0) = -3.;
    expected_result(1) = 6.;
    expected_result(2) = -9.;
    EXPECT_TRUE(equals(HostConstVector<double>(v1), HostConstVector<double>(expected_result)));
}

/* T dot_product(HostConstVector<T>& lhs, HostConstVector<T>& rhs); */

TEST(VectorOperations, dot_product)
{
    HostVector<double> v1("v1", 3);
    v1(0) = 1;
    v1(1) = 2;
    v1(2) = 3;
    HostVector<double> v2("v2", 3);
    v2(0) = 1;
    v2(1) = 5;
    v2(2) = 2;
    EXPECT_DOUBLE_EQ(dot_product(HostConstVector<double>(v1), HostConstVector<double>(v2)), 17.0);
}

/* T l1_norm(HostConstVector<T>& x); */

TEST(VectorOperations, l1_vector_norm)
{
    HostVector<double> v("v", 3);
    v(0) = 1;
    v(1) = -5;
    v(2) = 2;
    EXPECT_DOUBLE_EQ(l1_norm(HostConstVector<double>(v)), 8.0);
}

/* T l2_norm(HostConstVector<T>& x); */

TEST(VectorOperations, l2_vector_norm)
{
    HostVector<double> v("v", 3);
    v(0) = 1;
    v(1) = -5;
    v(2) = 2;
    HostConstVector<double> const_v(v);
    EXPECT_DOUBLE_EQ(l2_norm(HostConstVector<double>(const_v)), std::sqrt(30.0));
}

TEST(VectorOperations, zero_l2_vector_norm)
{
    HostVector<double> v("v", 3);
    v(0) = 0;
    v(1) = 0;
    v(2) = 0;
    HostConstVector<double> const_v(v);
    EXPECT_DOUBLE_EQ(l2_norm(HostConstVector<double>(const_v)), 0.0);
}

TEST(VectorOperations, underflow_l2_vector_norm)
{
    HostVector<double> v("v", 3);
    v(0) = 1e-300;
    v(1) = 1e-300;
    v(2) = 1e-300;
    HostConstVector<double> const_v(v);
    EXPECT_DOUBLE_EQ(l2_norm(HostConstVector<double>(const_v)), std::sqrt(3.0) * 1e-300);
}

TEST(VectorOperations, overflow_l2_vector_norm)
{
    HostVector<double> v("v", 3);
    v(0) = 1e+300;
    v(1) = 1e+300;
    v(2) = 1e+300;
    HostConstVector<double> const_v(v);
    EXPECT_DOUBLE_EQ(l2_norm(HostConstVector<double>(const_v)), std::sqrt(3.0) * 1e+300);
}

/* T infinity_norm(HostConstVector<T>& x); */

TEST(VectorOperations, infinity_vector_norm)
{
    HostVector<double> v("v", 3);
    v(0) = 1;
    v(1) = -5;
    v(2) = 2;
    HostConstVector<double> const_v(v);
    EXPECT_DOUBLE_EQ(infinity_norm(HostConstVector<double>(const_v)), 5.0);
}
