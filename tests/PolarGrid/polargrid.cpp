#include <gtest/gtest.h>
#include "../../include/PolarGrid/polargrid.h"
using namespace gmgpolar;

double compare_coords(PolarGrid<DefaultMemorySpace> grid, std::vector<double> radii, std::vector<double> angles)
{
    Kokkos::View<double*, Kokkos::HostSpace> host_vector_radi(radii.data(), radii.size());
    Kokkos::View<double*, Kokkos::HostSpace> host_vector_angles(angles.data(), angles.size());

    auto expected_r     = Kokkos::create_mirror_view_and_copy(DefaultMemorySpace(), host_vector_radi);
    auto expected_theta = Kokkos::create_mirror_view_and_copy(DefaultMemorySpace(), host_vector_angles);

    double r_err     = 0;
    double theta_err = 0;

    Kokkos::parallel_reduce(
        "r_test", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, grid.nr()),
        KOKKOS_LAMBDA(int i_r, double& error) { error += expected_r(i_r) - grid.radius(i_r); }, r_err);

    Kokkos::parallel_reduce(
        "tetha_test", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, grid.ntheta()),
        KOKKOS_LAMBDA(int i_theta, double& t_error) { t_error += expected_theta(i_theta) - grid.theta(i_theta); },
        theta_err);

    return std::max(r_err, theta_err);
}

double expected_indices(PolarGrid<DefaultMemorySpace> grid)
{

    double idx_err;
    Kokkos::parallel_reduce(
        "indices",
        Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>({0, 0}, {grid.nr(), grid.ntheta()}),
        KOKKOS_LAMBDA(const int i, const int j, double& error) {
            int node_index = grid.index(i, j);
            int r_out, theta_out;
            grid.multiIndex(node_index, r_out, theta_out);
            error += (i - r_out) + (j - theta_out);
        },
        idx_err);
    return idx_err;
}

double expected_spacing(PolarGrid<DefaultMemorySpace> grid, Vector<double> rad_expected, Vector<double> ang_expected)
{
    double radial_spacing_err  = 0;
    double angular_spacing_err = 0;
    Kokkos::parallel_reduce(
        "radial_spacing", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, rad_expected.size()),
        KOKKOS_LAMBDA(int i_r, double& error) { error += rad_expected(i_r) - grid.radialSpacing(i_r); },
        radial_spacing_err);

    Kokkos::parallel_reduce(
        "angular_spacing", Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, ang_expected.size()),
        KOKKOS_LAMBDA(int i_theta, double& error) { error += ang_expected(i_theta) - grid.angularSpacing(i_theta); },
        angular_spacing_err);
    return std::max(radial_spacing_err, angular_spacing_err);
}

TEST(PolarGridTest, DefaultConstructor)
{
    PolarGrid<DefaultMemorySpace> grid;
}

TEST(PolarGridTest, VectorConstructor)
{
    std::vector<double> radii  = {0.1, 0.2, 0.5, 0.9, 1.3};
    std::vector<double> angles = {0, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};
    PolarGrid<DefaultMemorySpace> grid(radii, angles);
}

TEST(PolarGridTest, NumberOfNodes)
{
    std::vector<double> radii  = {0.1, 0.2, 0.5, 0.9, 1.3};
    std::vector<double> angles = {0, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};
    PolarGrid<DefaultMemorySpace> grid(radii, angles);
    ASSERT_EQ(grid.numberOfNodes(), radii.size() * (angles.size() - 1));
}

TEST(PolarGridTest, AccessorsTest)
{
    std::vector<double> radii  = {0.1, 0.2, 0.5, 0.9, 1.3};
    std::vector<double> angles = {0, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};
    PolarGrid<DefaultMemorySpace> grid(radii, angles);

    ASSERT_DOUBLE_EQ(compare_coords(grid, radii, angles), 0.);
}

TEST(PolarGridTest, GridJumpTest)
{
    std::vector<double> radii  = {0.1, 0.2, 0.5, 0.9, 1.3};
    std::vector<double> angles = {0, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};
    double splitting_radius    = 0.4;
    PolarGrid<DefaultMemorySpace> grid(radii, angles, splitting_radius);

    ASSERT_DOUBLE_EQ(compare_coords(grid, radii, angles), 0.);
}

TEST(PolarGridTest, IndexingTest)
{
    std::vector<double> radii  = {0.1, 0.2, 0.25, 0.5, 0.8, 0.9, 1.3, 1.4, 2.0};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};
    double splitting_radius = 0.6;
    PolarGrid<DefaultMemorySpace> grid(radii, angles, splitting_radius);

    ASSERT_EQ(expected_indices(grid), 0.);

    for (int i = 0; i < grid.numberOfNodes(); i++) {
        int r_out, theta_out;
        grid.multiIndex(i, r_out, theta_out);
        int node_index = grid.index(r_out, theta_out);
        ASSERT_EQ(node_index, i);
    }
}

TEST(PolarGridTest, IndexingValuesTest)
{
    std::vector<double> radii  = {0.1, 0.2, 0.25, 0.5, 0.8, 0.9, 1.3, 1.4, 2.0};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};
    double splitting_radius = 0.6;
    PolarGrid<DefaultMemorySpace> grid(radii, angles, splitting_radius);

    {
        int node_index = grid.index(2, 6);
        ASSERT_EQ(node_index, 22);
        int r_out, theta_out;
        grid.multiIndex(node_index, r_out, theta_out);
        ASSERT_EQ(2, r_out);
        ASSERT_EQ(6, theta_out);
    }

    {
        int node_index = grid.index(3, 2);
        ASSERT_EQ(node_index, 26);
        int r_out, theta_out;
        grid.multiIndex(node_index, r_out, theta_out);
        ASSERT_EQ(3, r_out);
        ASSERT_EQ(2, theta_out);
    }

    {
        int node_index = grid.index(6, 4);
        ASSERT_EQ(node_index, 54);
        int r_out, theta_out;
        grid.multiIndex(node_index, r_out, theta_out);
        ASSERT_EQ(6, r_out);
        ASSERT_EQ(4, theta_out);
    }

    {
        int node_index = grid.index(4, 7);
        ASSERT_EQ(node_index, 67);
        int r_out, theta_out;
        grid.multiIndex(node_index, r_out, theta_out);
        ASSERT_EQ(4, r_out);
        ASSERT_EQ(7, theta_out);
    }
}

TEST(PolarGridTest, CoordinatesTest)
{
    std::vector<double> radii  = {0.1, 0.2, 0.25, 0.5, 0.8, 0.9, 1.3, 1.4, 2.0};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};
    double splitting_radius = 0.6;
    PolarGrid<DefaultMemorySpace> grid(radii, angles, splitting_radius);
    // Check that coordinates are correct
    ASSERT_DOUBLE_EQ(compare_coords(grid, radii, angles), 0.);
}

TEST(PolarGridTest, SpacingTest)
{
    std::vector<double> radii  = {0.1, 0.2, 0.25, 0.5, 0.8, 0.9, 1.3, 1.4, 2.0};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};
    PolarGrid<DefaultMemorySpace> grid(radii, angles);
    HostVector<double> h_exp_rad_spacing("h_rad_spacing", grid.nr() - 1);
    HostVector<double> h_exp_theta_spacing("h_theta_spacing", grid.ntheta() - 1);
    for (int i_r = 0; i_r < grid.nr() - 1; i_r++) {
        h_exp_rad_spacing(i_r) = radii[i_r + 1] - radii[i_r];
    }
    for (int i_theta = 0; i_theta < grid.ntheta() - 1; i_theta++) {
        h_exp_theta_spacing(i_theta) = angles[i_theta + 1] - angles[i_theta];
    }
    auto exp_rad_spacing   = Kokkos::create_mirror_view_and_copy(DefaultMemorySpace(), h_exp_rad_spacing);
    auto exp_theta_spacing = Kokkos::create_mirror_view_and_copy(DefaultMemorySpace(), h_exp_theta_spacing);

    // Test radial and angular spacings
    ASSERT_EQ(expected_spacing(grid, exp_rad_spacing, exp_theta_spacing), 0.);
}
