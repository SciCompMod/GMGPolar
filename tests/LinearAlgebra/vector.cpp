#include <complex>

#include <gtest/gtest.h>

#include "../../include/LinearAlgebra/vector.h"

#define EXPECT_DOUBLE_COMPLEX_EQ(a, b)                                                                                 \
    EXPECT_DOUBLE_EQ((a).real(), (b).real());                                                                          \
    EXPECT_DOUBLE_EQ((a).imag(), (b).imag())

// -------------- //
// Vector<double> //
// -------------- //

/* Default Construct (double) */

TEST(VectorDouble, default_construct)
{
    const Vector<double> v;
    (void)v;
}

/* Value Construct (double) */

TEST(VectorDouble, value_construct)
{
    const Vector<double> v = {1.0, 2.0, 3.0};
    ASSERT_EQ(v.size(), 3);
    EXPECT_DOUBLE_EQ(v[0], 1.0);
    EXPECT_DOUBLE_EQ(v[1], 2.0);
    EXPECT_DOUBLE_EQ(v[2], 3.0);
}

/* Size Construct (double) */

TEST(VectorDouble, size_construct)
{
    const Vector<double> v(5);
    ASSERT_EQ(v.size(), 5);
}

/* Copy Construct (double) */

TEST(VectorDouble, copy_construct)
{
    const Vector<double> v1 = {1.0, 2.0, 3.0};
    const Vector<double> v2 = v1;

    ASSERT_EQ(v1.size(), 3);
    EXPECT_DOUBLE_EQ(v1[0], 1.0);
    EXPECT_DOUBLE_EQ(v1[1], 2.0);
    EXPECT_DOUBLE_EQ(v1[2], 3.0);

    ASSERT_EQ(v2.size(), 3);
    EXPECT_DOUBLE_EQ(v2[0], 1.0);
    EXPECT_DOUBLE_EQ(v2[1], 2.0);
    EXPECT_DOUBLE_EQ(v2[2], 3.0);
}

/* Copy Assign (double) */

TEST(VectorDouble, copy_assign)
{
    Vector<double> v2;

    {
        const Vector<double> v1 = {1.0, 2.0, 3.0};
        v2                      = v1;

        ASSERT_EQ(v1.size(), 3);
        EXPECT_DOUBLE_EQ(v1[0], 1.0);
        EXPECT_DOUBLE_EQ(v1[1], 2.0);
        EXPECT_DOUBLE_EQ(v1[2], 3.0);

        ASSERT_EQ(v2.size(), 3);
        EXPECT_DOUBLE_EQ(v2[0], 1.0);
        EXPECT_DOUBLE_EQ(v2[1], 2.0);
        EXPECT_DOUBLE_EQ(v2[2], 3.0);
    }

    ASSERT_EQ(v2.size(), 3);
    EXPECT_DOUBLE_EQ(v2[0], 1.0);
    EXPECT_DOUBLE_EQ(v2[1], 2.0);
    EXPECT_DOUBLE_EQ(v2[2], 3.0);
}

/* Move Construct (double) */

TEST(VectorDouble, move_construct)
{
    Vector<double> v1 = {1.0, 2.0, 3.0};

    ASSERT_EQ(v1.size(), 3);
    EXPECT_DOUBLE_EQ(v1[0], 1.0);
    EXPECT_DOUBLE_EQ(v1[1], 2.0);
    EXPECT_DOUBLE_EQ(v1[2], 3.0);

    const Vector<double> v2 = std::move(v1);

    ASSERT_EQ(v2.size(), 3);
    EXPECT_DOUBLE_EQ(v2[0], 1.0);
    EXPECT_DOUBLE_EQ(v2[1], 2.0);
    EXPECT_DOUBLE_EQ(v2[2], 3.0);
}

/* Move Assign (double) */

TEST(VectorDouble, move_assign)
{
    Vector<double> v2;

    {
        const Vector<double> v1 = {1.0, 2.0, 3.0};
        ASSERT_EQ(v1.size(), 3);
        EXPECT_DOUBLE_EQ(v1[0], 1.0);
        EXPECT_DOUBLE_EQ(v1[1], 2.0);
        EXPECT_DOUBLE_EQ(v1[2], 3.0);

        v2 = std::move(v1);

        ASSERT_EQ(v2.size(), 3);
        EXPECT_DOUBLE_EQ(v2[0], 1.0);
        EXPECT_DOUBLE_EQ(v2[1], 2.0);
        EXPECT_DOUBLE_EQ(v2[2], 3.0);
    }

    ASSERT_EQ(v2.size(), 3);
    EXPECT_DOUBLE_EQ(v2[0], 1.0);
    EXPECT_DOUBLE_EQ(v2[1], 2.0);
    EXPECT_DOUBLE_EQ(v2[2], 3.0);
}

/* Value Construct Modify (double) */

TEST(VectorDouble, value_construct_modify)
{
    Vector<double> v = {1.0, 2.0, 3.0};
    ASSERT_EQ(v.size(), 3);
    EXPECT_DOUBLE_EQ(v[0], 1.0);
    EXPECT_DOUBLE_EQ(v[1], 2.0);
    EXPECT_DOUBLE_EQ(v[2], 3.0);

    v[2] = 5.0;
    ASSERT_EQ(v.size(), 3);
    EXPECT_DOUBLE_EQ(v[0], 1.0);
    EXPECT_DOUBLE_EQ(v[1], 2.0);
    EXPECT_DOUBLE_EQ(v[2], 5.0);
}

/* Size Construct Modify (double) */

TEST(VectorDouble, size_construct_modify)
{
    Vector<double> v(5);
    ASSERT_EQ(v.size(), 5);

    v[0] = 1.0;
    v[1] = 2.0;
    v[2] = 3.0;
    v[3] = 4.0;
    v[4] = 5.0;

    EXPECT_DOUBLE_EQ(v[0], 1.0);
    EXPECT_DOUBLE_EQ(v[1], 2.0);
    EXPECT_DOUBLE_EQ(v[2], 3.0);
    EXPECT_DOUBLE_EQ(v[3], 4.0);
    EXPECT_DOUBLE_EQ(v[4], 5.0);
}

// ------------- //
// Vector<float> //
// ------------- //

TEST(VectorFloat, default_construct)
{
    const Vector<float> v;
    (void)v;
}

TEST(VectorFloat, value_construct)
{
    const Vector<float> v = {1.0, 2.0, 3.0};
    ASSERT_EQ(v.size(), 3);
    EXPECT_FLOAT_EQ(v[0], 1.0);
    EXPECT_FLOAT_EQ(v[1], 2.0);
    EXPECT_FLOAT_EQ(v[2], 3.0);
}

TEST(VectorFloat, size_construct)
{
    const Vector<float> v(5);
    ASSERT_EQ(v.size(), 5);
}

TEST(VectorFloat, copy_construct)
{
    const Vector<float> v1 = {1.0, 2.0, 3.0};
    const Vector<float> v2 = v1;

    ASSERT_EQ(v1.size(), 3);
    EXPECT_FLOAT_EQ(v1[0], 1.0);
    EXPECT_FLOAT_EQ(v1[1], 2.0);
    EXPECT_FLOAT_EQ(v1[2], 3.0);

    ASSERT_EQ(v2.size(), 3);
    EXPECT_FLOAT_EQ(v2[0], 1.0);
    EXPECT_FLOAT_EQ(v2[1], 2.0);
    EXPECT_FLOAT_EQ(v2[2], 3.0);
}

TEST(VectorFloat, copy_assign)
{
    Vector<float> v2;

    {
        const Vector<float> v1 = {1.0, 2.0, 3.0};
        v2                     = v1;

        ASSERT_EQ(v1.size(), 3);
        EXPECT_FLOAT_EQ(v1[0], 1.0);
        EXPECT_FLOAT_EQ(v1[1], 2.0);
        EXPECT_FLOAT_EQ(v1[2], 3.0);

        ASSERT_EQ(v2.size(), 3);
        EXPECT_FLOAT_EQ(v2[0], 1.0);
        EXPECT_FLOAT_EQ(v2[1], 2.0);
        EXPECT_FLOAT_EQ(v2[2], 3.0);
    }

    ASSERT_EQ(v2.size(), 3);
    EXPECT_FLOAT_EQ(v2[0], 1.0);
    EXPECT_FLOAT_EQ(v2[1], 2.0);
    EXPECT_FLOAT_EQ(v2[2], 3.0);
}

TEST(VectorFloat, move_construct)
{
    Vector<float> v1 = {1.0, 2.0, 3.0};

    ASSERT_EQ(v1.size(), 3);
    EXPECT_FLOAT_EQ(v1[0], 1.0);
    EXPECT_FLOAT_EQ(v1[1], 2.0);
    EXPECT_FLOAT_EQ(v1[2], 3.0);

    const Vector<float> v2 = std::move(v1);

    ASSERT_EQ(v2.size(), 3);
    EXPECT_FLOAT_EQ(v2[0], 1.0);
    EXPECT_FLOAT_EQ(v2[1], 2.0);
    EXPECT_FLOAT_EQ(v2[2], 3.0);
}

TEST(VectorFloat, move_assign)
{
    Vector<float> v2;

    {
        const Vector<float> v1 = {1.0, 2.0, 3.0};
        ASSERT_EQ(v1.size(), 3);
        EXPECT_FLOAT_EQ(v1[0], 1.0);
        EXPECT_FLOAT_EQ(v1[1], 2.0);
        EXPECT_FLOAT_EQ(v1[2], 3.0);

        v2 = std::move(v1);

        ASSERT_EQ(v2.size(), 3);
        EXPECT_FLOAT_EQ(v2[0], 1.0);
        EXPECT_FLOAT_EQ(v2[1], 2.0);
        EXPECT_FLOAT_EQ(v2[2], 3.0);
    }

    ASSERT_EQ(v2.size(), 3);
    EXPECT_FLOAT_EQ(v2[0], 1.0);
    EXPECT_FLOAT_EQ(v2[1], 2.0);
    EXPECT_FLOAT_EQ(v2[2], 3.0);
}

TEST(VectorFloat, value_construct_modify)
{
    Vector<float> v = {1.0, 2.0, 3.0};
    ASSERT_EQ(v.size(), 3);
    EXPECT_FLOAT_EQ(v[0], 1.0);
    EXPECT_FLOAT_EQ(v[1], 2.0);
    EXPECT_FLOAT_EQ(v[2], 3.0);

    v[2] = 5.0;
    ASSERT_EQ(v.size(), 3);
    EXPECT_FLOAT_EQ(v[0], 1.0);
    EXPECT_FLOAT_EQ(v[1], 2.0);
    EXPECT_FLOAT_EQ(v[2], 5.0);
}

TEST(VectorFloat, size_construct_modify)
{
    Vector<float> v(5);
    ASSERT_EQ(v.size(), 5);

    v[0] = 1.0;
    v[1] = 2.0;
    v[2] = 3.0;
    v[3] = 4.0;
    v[4] = 5.0;

    EXPECT_FLOAT_EQ(v[0], 1.0);
    EXPECT_FLOAT_EQ(v[1], 2.0);
    EXPECT_FLOAT_EQ(v[2], 3.0);
    EXPECT_FLOAT_EQ(v[3], 4.0);
    EXPECT_FLOAT_EQ(v[4], 5.0);
}

// ---------------------------- //
// Vector<std::complex<double>> //
// ---------------------------- //

TEST(VectorDoubleComplex, default_construct)
{
    const Vector<std::complex<double>> v;
    (void)v;
}

TEST(VectorDoubleComplex, value_construct)
{
    using namespace std::complex_literals;
    const Vector<std::complex<double>> v = {1.0 + 4.0i, 2.0 + 5.0i, 3.0 + 6.0i};
    ASSERT_EQ(v.size(), 3);
    EXPECT_DOUBLE_COMPLEX_EQ(v[0], 1.0 + 4.0i);
    EXPECT_DOUBLE_COMPLEX_EQ(v[1], 2.0 + 5.0i);
    EXPECT_DOUBLE_COMPLEX_EQ(v[2], 3.0 + 6.0i);
}

TEST(VectorDoubleComplex, size_construct)
{
    const Vector<std::complex<double>> v(5);
    ASSERT_EQ(v.size(), 5);
}

TEST(VectorDoubleComplex, copy_construct)
{
    using namespace std::complex_literals;
    const Vector<std::complex<double>> v1 = {1.0 + 4.0i, 2.0 + 5.0i, 3.0 + 6.0i};
    const Vector<std::complex<double>> v2 = v1;

    ASSERT_EQ(v1.size(), 3);
    EXPECT_DOUBLE_COMPLEX_EQ(v1[0], 1.0 + 4.0i);
    EXPECT_DOUBLE_COMPLEX_EQ(v1[1], 2.0 + 5.0i);
    EXPECT_DOUBLE_COMPLEX_EQ(v1[2], 3.0 + 6.0i);

    ASSERT_EQ(v2.size(), 3);
    EXPECT_DOUBLE_COMPLEX_EQ(v2[0], 1.0 + 4.0i);
    EXPECT_DOUBLE_COMPLEX_EQ(v2[1], 2.0 + 5.0i);
    EXPECT_DOUBLE_COMPLEX_EQ(v2[2], 3.0 + 6.0i);
}

TEST(VectorDoubleComplex, copy_assign)
{
    using namespace std::complex_literals;
    Vector<std::complex<double>> v2;

    {
        const Vector<std::complex<double>> v1 = {1.0 + 4.0i, 2.0 + 5.0i, 3.0 + 6.0i};
        v2                                    = v1;

        ASSERT_EQ(v1.size(), 3);
        EXPECT_DOUBLE_COMPLEX_EQ(v1[0], 1.0 + 4.0i);
        EXPECT_DOUBLE_COMPLEX_EQ(v1[1], 2.0 + 5.0i);
        EXPECT_DOUBLE_COMPLEX_EQ(v1[2], 3.0 + 6.0i);

        ASSERT_EQ(v2.size(), 3);
        EXPECT_DOUBLE_COMPLEX_EQ(v2[0], 1.0 + 4.0i);
        EXPECT_DOUBLE_COMPLEX_EQ(v2[1], 2.0 + 5.0i);
        EXPECT_DOUBLE_COMPLEX_EQ(v2[2], 3.0 + 6.0i);
    }

    ASSERT_EQ(v2.size(), 3);
    EXPECT_DOUBLE_COMPLEX_EQ(v2[0], 1.0 + 4.0i);
    EXPECT_DOUBLE_COMPLEX_EQ(v2[1], 2.0 + 5.0i);
    EXPECT_DOUBLE_COMPLEX_EQ(v2[2], 3.0 + 6.0i);
}

TEST(VectorDoubleComplex, move_construct)
{
    using namespace std::complex_literals;
    Vector<std::complex<double>> v1 = {1.0 + 4.0i, 2.0 + 5.0i, 3.0 + 6.0i};

    ASSERT_EQ(v1.size(), 3);
    EXPECT_DOUBLE_COMPLEX_EQ(v1[0], 1.0 + 4.0i);
    EXPECT_DOUBLE_COMPLEX_EQ(v1[1], 2.0 + 5.0i);
    EXPECT_DOUBLE_COMPLEX_EQ(v1[2], 3.0 + 6.0i);

    const Vector<std::complex<double>> v2 = std::move(v1);

    ASSERT_EQ(v2.size(), 3);
    EXPECT_DOUBLE_COMPLEX_EQ(v2[0], 1.0 + 4.0i);
    EXPECT_DOUBLE_COMPLEX_EQ(v2[1], 2.0 + 5.0i);
    EXPECT_DOUBLE_COMPLEX_EQ(v2[2], 3.0 + 6.0i);
}

TEST(VectorDoubleComplex, move_assign)
{
    using namespace std::complex_literals;
    Vector<std::complex<double>> v2;

    {
        const Vector<std::complex<double>> v1 = {1.0 + 4.0i, 2.0 + 5.0i, 3.0 + 6.0i};
        ASSERT_EQ(v1.size(), 3);
        EXPECT_DOUBLE_COMPLEX_EQ(v1[0], 1.0 + 4.0i);
        EXPECT_DOUBLE_COMPLEX_EQ(v1[1], 2.0 + 5.0i);
        EXPECT_DOUBLE_COMPLEX_EQ(v1[2], 3.0 + 6.0i);

        v2 = std::move(v1);

        ASSERT_EQ(v2.size(), 3);
        EXPECT_DOUBLE_COMPLEX_EQ(v2[0], 1.0 + 4.0i);
        EXPECT_DOUBLE_COMPLEX_EQ(v2[1], 2.0 + 5.0i);
        EXPECT_DOUBLE_COMPLEX_EQ(v2[2], 3.0 + 6.0i);
    }

    ASSERT_EQ(v2.size(), 3);
    EXPECT_DOUBLE_COMPLEX_EQ(v2[0], 1.0 + 4.0i);
    EXPECT_DOUBLE_COMPLEX_EQ(v2[1], 2.0 + 5.0i);
    EXPECT_DOUBLE_COMPLEX_EQ(v2[2], 3.0 + 6.0i);
}

TEST(VectorDoubleComplex, value_construct_modify)
{
    using namespace std::complex_literals;
    Vector<std::complex<double>> v = {1.0 + 4.0i, 2.0 + 5.0i, 3.0 + 6.0i};
    ASSERT_EQ(v.size(), 3);
    EXPECT_DOUBLE_COMPLEX_EQ(v[0], 1.0 + 4.0i);
    EXPECT_DOUBLE_COMPLEX_EQ(v[1], 2.0 + 5.0i);
    EXPECT_DOUBLE_COMPLEX_EQ(v[2], 3.0 + 6.0i);

    v[2] = 5.0;
    ASSERT_EQ(v.size(), 3);
    EXPECT_DOUBLE_COMPLEX_EQ(v[0], 1.0 + 4.0i);
    EXPECT_DOUBLE_COMPLEX_EQ(v[1], 2.0 + 5.0i);
    EXPECT_DOUBLE_COMPLEX_EQ(v[2], 5.0 + 0.0i);
}

TEST(VectorDoubleComplex, size_construct_modify)
{
    using namespace std::complex_literals;
    Vector<std::complex<double>> v(5);
    ASSERT_EQ(v.size(), 5);

    v[0] = 1.0 + 6.0i;
    v[1] = 2.0 + 7.0i;
    v[2] = 3.0 + 8.0i;
    v[3] = 4.0 + 9.0i;
    v[4] = 5.0 + 10.0i;

    EXPECT_DOUBLE_COMPLEX_EQ(v[0], 1.0 + 6.0i);
    EXPECT_DOUBLE_COMPLEX_EQ(v[1], 2.0 + 7.0i);
    EXPECT_DOUBLE_COMPLEX_EQ(v[2], 3.0 + 8.0i);
    EXPECT_DOUBLE_COMPLEX_EQ(v[3], 4.0 + 9.0i);
    EXPECT_DOUBLE_COMPLEX_EQ(v[4], 5.0 + 10.0i);
}

int main(int argc, char* argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}