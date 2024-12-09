#include <gtest/gtest.h>
#include <complex>

#include "../../include/LinearAlgebra/Vector/gpu_vector.h"

#define EXPECT_DOUBLE_COMPLEX_EQ(a, b) \
    do { \
        EXPECT_DOUBLE_EQ((a).real(), (b).real()); \
        EXPECT_DOUBLE_EQ((a).imag(), (b).imag()); \
    } while (0)


/* Default Construct (double) */
TEST(GPU_VectorDouble, default_construct)
{
    const GPU_Vector<double> v;
    (void)v;
}

/* Size Construct (double) */
TEST(GPU_VectorDouble, size_construct)
{
    const Vector<double> v = {1.0, 2.0, 3.0};

    GPU_Vector<double> gpu_v(v.size());
    ASSERT_EQ(gpu_v.size(), 3);
    copyHostToDevice(v, gpu_v);

    Vector<double> host_v(gpu_v.size());
    copyDeviceToHost(gpu_v, host_v);

    ASSERT_EQ(host_v.size(), 3);
    EXPECT_DOUBLE_EQ(host_v[0], 1.0);
    EXPECT_DOUBLE_EQ(host_v[1], 2.0);
    EXPECT_DOUBLE_EQ(host_v[2], 3.0);
}


/* Copy Construct (double) */
TEST(GPU_VectorDouble, copy_construct)
{
    const Vector<double> v = {1.0, 2.0, 3.0};

    GPU_Vector<double> gpu_v(v.size());
    ASSERT_EQ(gpu_v.size(), 3);
    copyHostToDevice(v, gpu_v);

    GPU_Vector<double> gpu_w = gpu_v;

    Vector<double> host_v(gpu_w.size());
    copyDeviceToHost(gpu_w, host_v);

    ASSERT_EQ(host_v.size(), 3);
    EXPECT_DOUBLE_EQ(host_v[0], 1.0);
    EXPECT_DOUBLE_EQ(host_v[1], 2.0);
    EXPECT_DOUBLE_EQ(host_v[2], 3.0);
}

/* Copy Assign (double) */
TEST(GPU_VectorDouble, copy_assign)
{
    const Vector<double> v = {1.0, 2.0, 3.0};

    GPU_Vector<double> gpu_v(v.size());
    ASSERT_EQ(gpu_v.size(), 3);
    copyHostToDevice(v, gpu_v);

    GPU_Vector<double> gpu_w;
    gpu_w = gpu_v;

    Vector<double> host_v(gpu_w.size());
    copyDeviceToHost(gpu_w, host_v);

    ASSERT_EQ(host_v.size(), 3);
    EXPECT_DOUBLE_EQ(host_v[0], 1.0);
    EXPECT_DOUBLE_EQ(host_v[1], 2.0);
    EXPECT_DOUBLE_EQ(host_v[2], 3.0);
}

/* Move Construct (double) */
TEST(GPU_VectorDouble, move_construct)
{
    const Vector<double> v = {1.0, 2.0, 3.0};

    GPU_Vector<double> gpu_v(v.size());
    ASSERT_EQ(gpu_v.size(), 3);
    copyHostToDevice(v, gpu_v);

    GPU_Vector<double> gpu_w = std::move(gpu_v);

    Vector<double> host_v(gpu_w.size());
    copyDeviceToHost(gpu_w, host_v);

    ASSERT_EQ(host_v.size(), 3);
    EXPECT_DOUBLE_EQ(host_v[0], 1.0);
    EXPECT_DOUBLE_EQ(host_v[1], 2.0);
    EXPECT_DOUBLE_EQ(host_v[2], 3.0);
}

/* Move Assign (double) */
TEST(GPU_VectorDouble, move_assign)
{
    const Vector<double> v = {1.0, 2.0, 3.0};

    GPU_Vector<double> gpu_v(v.size());
    ASSERT_EQ(gpu_v.size(), 3);
    copyHostToDevice(v, gpu_v);

    GPU_Vector<double> gpu_w;
    gpu_w = std::move(gpu_v);

    Vector<double> host_v(gpu_w.size());
    copyDeviceToHost(gpu_w, host_v);

    ASSERT_EQ(host_v.size(), 3);
    EXPECT_DOUBLE_EQ(host_v[0], 1.0);
    EXPECT_DOUBLE_EQ(host_v[1], 2.0);
    EXPECT_DOUBLE_EQ(host_v[2], 3.0);
}


/* Default Construct (std::complex<double>) */
TEST(GPU_VectorComplex, default_construct)
{
    using namespace std::complex_literals;
    const GPU_Vector<std::complex<double>> v;
    (void)v;
}

/* Size Construct (std::complex<double>) */
TEST(GPU_VectorComplex, size_construct)
{
    using namespace std::complex_literals;
    const Vector<std::complex<double>> v = {1.0 + 4.0i, 2.0 + 5.0i, 3.0 + 6.0i};

    GPU_Vector<std::complex<double>> gpu_v(v.size());
    ASSERT_EQ(gpu_v.size(), 3);
    copyHostToDevice(v, gpu_v);

    Vector<std::complex<double>> host_v(gpu_v.size());
    copyDeviceToHost(gpu_v, host_v);

    ASSERT_EQ(host_v.size(), 3);
    EXPECT_DOUBLE_COMPLEX_EQ(host_v[0], 1.0 + 4.0i);
    EXPECT_DOUBLE_COMPLEX_EQ(host_v[1], 2.0 + 5.0i);
    EXPECT_DOUBLE_COMPLEX_EQ(host_v[2], 3.0 + 6.0i);
}

/* Copy Construct (std::complex<double>) */
TEST(GPU_VectorComplex, copy_construct)
{
    using namespace std::complex_literals;
    const Vector<std::complex<double>> v = {1.0 + 4.0i, 2.0 + 5.0i, 3.0 + 6.0i};

    GPU_Vector<std::complex<double>> gpu_v(v.size());
    ASSERT_EQ(gpu_v.size(), 3);
    copyHostToDevice(v, gpu_v);

    GPU_Vector<std::complex<double>> gpu_w = gpu_v;

    Vector<std::complex<double>> host_v(gpu_w.size());
    copyDeviceToHost(gpu_w, host_v);

    ASSERT_EQ(host_v.size(), 3);
    EXPECT_DOUBLE_COMPLEX_EQ(host_v[0], 1.0 + 4.0i);
    EXPECT_DOUBLE_COMPLEX_EQ(host_v[1], 2.0 + 5.0i);
    EXPECT_DOUBLE_COMPLEX_EQ(host_v[2], 3.0 + 6.0i);
}

/* Copy Assign (std::complex<double>) */
TEST(GPU_VectorComplex, copy_assign)
{
    using namespace std::complex_literals;
    const Vector<std::complex<double>> v = {1.0 + 4.0i, 2.0 + 5.0i, 3.0 + 6.0i};

    GPU_Vector<std::complex<double>> gpu_v(v.size());
    ASSERT_EQ(gpu_v.size(), 3);
    copyHostToDevice(v, gpu_v);

    GPU_Vector<std::complex<double>> gpu_w;
    gpu_w = gpu_v;

    Vector<std::complex<double>> host_v(gpu_w.size());
    copyDeviceToHost(gpu_w, host_v);

    ASSERT_EQ(host_v.size(), 3);
    EXPECT_DOUBLE_COMPLEX_EQ(host_v[0], 1.0 + 4.0i);
    EXPECT_DOUBLE_COMPLEX_EQ(host_v[1], 2.0 + 5.0i);
    EXPECT_DOUBLE_COMPLEX_EQ(host_v[2], 3.0 + 6.0i);
}

/* Move Construct (std::complex<double>) */
TEST(GPU_VectorComplex, move_construct)
{
    using namespace std::complex_literals;
    const Vector<std::complex<double>> v = {1.0 + 4.0i, 2.0 + 5.0i, 3.0 + 6.0i};

    GPU_Vector<std::complex<double>> gpu_v(v.size());
    ASSERT_EQ(gpu_v.size(), 3);
    copyHostToDevice(v, gpu_v);

    GPU_Vector<std::complex<double>> gpu_w = std::move(gpu_v);

    Vector<std::complex<double>> host_v(gpu_w.size());
    copyDeviceToHost(gpu_w, host_v);

    ASSERT_EQ(host_v.size(), 3);
    EXPECT_DOUBLE_COMPLEX_EQ(host_v[0], 1.0 + 4.0i);
    EXPECT_DOUBLE_COMPLEX_EQ(host_v[1], 2.0 + 5.0i);
    EXPECT_DOUBLE_COMPLEX_EQ(host_v[2], 3.0 + 6.0i);
}

/* Move Assign (std::complex<double>) */
TEST(GPU_VectorComplex, move_assign)
{
    using namespace std::complex_literals;
    const Vector<std::complex<double>> v = {1.0 + 4.0i, 2.0 + 5.0i, 3.0 + 6.0i};

    GPU_Vector<std::complex<double>> gpu_v(v.size());
    ASSERT_EQ(gpu_v.size(), 3);
    copyHostToDevice(v, gpu_v);

    GPU_Vector<std::complex<double>> gpu_w;
    gpu_w = std::move(gpu_v);

    Vector<std::complex<double>> host_v(gpu_w.size());
    copyDeviceToHost(gpu_w, host_v);

    ASSERT_EQ(host_v.size(), 3);
    EXPECT_DOUBLE_COMPLEX_EQ(host_v[0], 1.0 + 4.0i);
    EXPECT_DOUBLE_COMPLEX_EQ(host_v[1], 2.0 + 5.0i);
    EXPECT_DOUBLE_COMPLEX_EQ(host_v[2], 3.0 + 6.0i);
}

int main(int argc, char* argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}