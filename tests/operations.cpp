#include <gtest/gtest.h>

#include "../include/linear_algebra/vector.h"
#include "../include/linear_algebra/operations.h"

TEST(Operations, equals_vector_vector)
{
    const Vector<double> v1 = {1, 2, 3};
    const Vector<double> v2 = {1, 2, 3};
    const Vector<double> v3 = {1, 3, 3};

    EXPECT_TRUE(equals(v1, v1));
    EXPECT_TRUE(equals(v2, v2));
    EXPECT_TRUE(equals(v3, v3));

    EXPECT_TRUE(equals(v1, v2));
    EXPECT_TRUE(equals(v2, v1));

    EXPECT_FALSE(equals(v1, v3));
    EXPECT_FALSE(equals(v3, v1));

    EXPECT_FALSE(equals(v2, v3));
    EXPECT_FALSE(equals(v3, v2));
}

TEST(Operations, assign_vector_scalar)
{
    Vector<double> v = {1, 2, 3};
    assign(v, 5.0);
    const Vector<double> expected_result = {5, 5, 5};
    EXPECT_TRUE(equals(v, expected_result));
}

TEST(Operations, subtract_vector_vector)
{
    const Vector<double> v1 = {1, 2, 3};
    const Vector<double> v2 = {1, 5, 2};
    Vector<double> r(3);
    subtract(r, v1, v2);
    const Vector<double> expected_result = {0, -3, 1};
    EXPECT_TRUE(equals(r, expected_result));
}

TEST(Operations, dot_product)
{
    const Vector<double> v1 = {1, 2, 3};
    const Vector<double> v2 = {1, 5, 2};
    EXPECT_DOUBLE_EQ(dot_product(v1, v2), 17.0);
}

int main(int argc, char* argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}