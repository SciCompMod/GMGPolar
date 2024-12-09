#include <gtest/gtest.h>

#include "../../include/LinearAlgebra/Vector/vector.h"
#include "../../include/LinearAlgebra/Vector/vector_operations.h"
#include "../../include/LinearAlgebra/Vector/gpu_vector.h"
#include "../../include/LinearAlgebra/Vector/gpu_vector_operations.h"

/* void assign(GPU_Vector<T>& lhs, const T& value); */

TEST(GPU_VectorOperations, assign_vector_scalar)
{
    Vector<double> v = {1, 2, 3};

    GPU_Vector<double> gpu_v(v.size());
    copyHostToDevice(v, gpu_v);

    assign(gpu_v, 5.0);

    Vector<double> v_host(gpu_v.size());
    copyDeviceToHost(gpu_v, v_host);

    const Vector<double> expected_result = {5, 5, 5};
    EXPECT_EQ(v_host.size(), expected_result.size()) << "Vector sizes do not match.";
    for (size_t i = 0; i < v_host.size(); ++i) {
        EXPECT_DOUBLE_EQ(v_host[i], expected_result[i]) << "Element mismatch at index " << i;
    }
}

/* void add(GPU_Vector<T>& result, const GPU_Vector<T>& x); */

TEST(GPU_VectorOperations, add_vector_vector)
{
    const Vector<double> v1 = {1, 2, 3};
    const Vector<double> v2 = {-1, -5, 2};

    GPU_Vector<double> gpu_v1(v1.size());
    copyHostToDevice(v1, gpu_v1);

    GPU_Vector<double> gpu_v2(v2.size());
    copyHostToDevice(v2, gpu_v2);

    add(gpu_v1, gpu_v2);

    Vector<double> v_host(gpu_v1.size());
    copyDeviceToHost(gpu_v1, v_host);

    const Vector<double> expected_result = {0, -3, 5};
    EXPECT_EQ(v_host.size(), expected_result.size()) << "Vector sizes do not match.";
    for (size_t i = 0; i < v_host.size(); ++i) {
        EXPECT_DOUBLE_EQ(v_host[i], expected_result[i]) << "Element mismatch at index " << i;
    }
}

/* void subtract(GPU_Vector<T>& result, const GPU_Vector<T>& x) */

TEST(GPU_VectorOperations, subtract_vector_vector)
{
    const Vector<double> v1 = {1, 2, 3};
    const Vector<double> v2 = {-1, -5, 2};

    GPU_Vector<double> gpu_v1(v1.size());
    copyHostToDevice(v1, gpu_v1);

    GPU_Vector<double> gpu_v2(v2.size());
    copyHostToDevice(v2, gpu_v2);

    subtract(gpu_v1, gpu_v2);

    Vector<double> v_host(gpu_v1.size());
    copyDeviceToHost(gpu_v1, v_host);

    const Vector<double> expected_result = {2, 7, 1};
    EXPECT_EQ(v_host.size(), expected_result.size()) << "Vector sizes do not match.";
    for (size_t i = 0; i < v_host.size(); ++i) {
        EXPECT_DOUBLE_EQ(v_host[i], expected_result[i]) << "Element mismatch at index " << i;
    }
}

/* void linear_combination(GPU_Vector<T>& x, const T& alpha, const GPU_Vector<T>& y, const T& beta); */

TEST(GPU_VectorOperations, linear_combination)
{
    const Vector<double> v1 = {1, 2, 3};
    const Vector<double> v2 = {-1, -5, 2};

    const double alpha    = -3.0;
    const double beta       = 2.0;

    GPU_Vector<double> gpu_v1(v1.size());
    copyHostToDevice(v1, gpu_v1);

    GPU_Vector<double> gpu_v2(v2.size());
    copyHostToDevice(v2, gpu_v2);

    linear_combination(gpu_v1, alpha, gpu_v2, beta);

    Vector<double> v_host(gpu_v1.size());
    copyDeviceToHost(gpu_v1, v_host);

    const Vector<double> expected_result = {-5.0, -16.0, -5.0};
    EXPECT_EQ(v_host.size(), expected_result.size()) << "Vector sizes do not match.";
    for (size_t i = 0; i < v_host.size(); ++i) {
        EXPECT_DOUBLE_EQ(v_host[i], expected_result[i]) << "Element mismatch at index " << i;
    }
}


/* void multiply(GPU_Vector<T>& x, const T& alpha); */

TEST(GPU_VectorOperations, multiply_vector_scalar)
{
    const Vector<double> v1 = {1, -2, 3};
    const double alpha    = -3.0;

    GPU_Vector<double> gpu_v1(v1.size());
    copyHostToDevice(v1, gpu_v1);

    multiply(gpu_v1, alpha);

    Vector<double> v_host(gpu_v1.size());
    copyDeviceToHost(gpu_v1, v_host);

    const Vector<double> expected_result = {-3.0, 6.0, -9.0};
    EXPECT_EQ(v_host.size(), expected_result.size()) << "Vector sizes do not match.";
    for (size_t i = 0; i < v_host.size(); ++i) {
        EXPECT_DOUBLE_EQ(v_host[i], expected_result[i]) << "Element mismatch at index " << i;
    }
}

/* T dot_product(const GPU_Vector<T>& lhs, const GPU_Vector<T>& rhs); */

TEST(GPU_VectorOperations, dot_product)
{
    const Vector<double> v1 = {1, 2, 3};
    const Vector<double> v2 = {1, 5, 2};

    GPU_Vector<double> gpu_v1(v1.size());
    copyHostToDevice(v1, gpu_v1);

    GPU_Vector<double> gpu_v2(v2.size());
    copyHostToDevice(v2, gpu_v2);

    EXPECT_DOUBLE_EQ(dot_product(gpu_v1, gpu_v2), 17.0);
}

/* T l1_norm(const GPU_Vector<T>& x); */

TEST(GPU_VectorOperations, l1_vector_norm)
{
    const Vector<double> v = {1, -5, 2};

    GPU_Vector<double> gpu_v(v.size());
    copyHostToDevice(v, gpu_v);

    EXPECT_DOUBLE_EQ(l1_norm(gpu_v), 8.0);
}

/* T l2_norm(const GPU_Vector<T>& x); */

TEST(GPU_VectorOperations, l2_vector_norm)
{
    const Vector<double> v = {1, -5, 2};
    GPU_Vector<double> gpu_v(v.size());
    copyHostToDevice(v, gpu_v);

    EXPECT_DOUBLE_EQ(l2_norm(gpu_v), sqrt(30.0));
}

/* T infinity_norm(const GPU_Vector<T>& x); */

TEST(GPU_VectorOperations, infinity_vector_norm)
{
    const Vector<double> v = {1, -5, 2};

    GPU_Vector<double> gpu_v(v.size());
    copyHostToDevice(v, gpu_v);

    EXPECT_DOUBLE_EQ(infinity_norm(gpu_v), 5.0);
}

/* Test for larger GPU_Vectors */

/* T dot_product(const GPU_Vector<T>& lhs, const GPU_Vector<T>& rhs); */

TEST(GPU_VectorOperations, large_dot_product)
{
    const int n = 100000;
    Vector<double> v1(n);
    Vector<double> v2(n);
    assign(v1, 2.0);
    assign(v2, -3.0);

    GPU_Vector<double> gpu_v1(v1.size());
    copyHostToDevice(v1, gpu_v1);

    GPU_Vector<double> gpu_v2(v2.size());
    copyHostToDevice(v2, gpu_v2);

    EXPECT_DOUBLE_EQ(dot_product(v1, v2), n * 2.0 * (-3.0));
    EXPECT_DOUBLE_EQ(dot_product(gpu_v1, gpu_v2), n * 2.0 * (-3.0));
}

/* T l1_norm(const GPU_Vector<T>& x); */

TEST(GPU_VectorOperations, large_l1_vector_norm)
{
    const int n = 100000;
    Vector<double> v(n);
    assign(v, -3.0);

    GPU_Vector<double> gpu_v(v.size());
    copyHostToDevice(v, gpu_v);

    EXPECT_DOUBLE_EQ(l1_norm(v), n * 3.0);
    EXPECT_DOUBLE_EQ(l1_norm(gpu_v), n * 3.0);
}

/* T l2_norm(const GPU_Vector<T>& x); */

TEST(GPU_VectorOperations, large_l2_vector_norm)
{
    const int n = 100000;
    Vector<double> v(n);
    assign(v, -3.0);

    GPU_Vector<double> gpu_v(v.size());
    copyHostToDevice(v, gpu_v);

    EXPECT_DOUBLE_EQ(l2_norm(v), sqrt(n * (-3.0) * (-3.0)));
    EXPECT_DOUBLE_EQ(l2_norm(gpu_v), sqrt(n * (-3.0) * (-3.0)));
}

/* T infinity_norm(const GPU_Vector<T>& x); */

TEST(GPU_VectorOperations, large_infinity_vector_norm)
{
    const int n = 100000;
    Vector<double> v(n);
    assign(v, -3.0);

    GPU_Vector<double> gpu_v(v.size());
    copyHostToDevice(v, gpu_v);

    EXPECT_DOUBLE_EQ(infinity_norm(v), 3.0);
    EXPECT_DOUBLE_EQ(infinity_norm(gpu_v), 3.0);
}


int main(int argc, char* argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}