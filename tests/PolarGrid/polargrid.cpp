#include <gtest/gtest.h>
#include "../../include/PolarGrid/polargrid.h"

TEST(PolarGridTest, DefaultConstructor)
{
    PolarGrid grid;
}

TEST(PolarGridTest, VectorConstructor)
{
    std::vector<double> radii  = {0.1, 0.2, 0.5, 0.9, 1.3};
    std::vector<double> angles = {0, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};
    PolarGrid grid(radii, angles);
}

TEST(PolarGridTest, NumberOfNodes)
{
    std::vector<double> radii  = {0.1, 0.2, 0.5, 0.9, 1.3};
    std::vector<double> angles = {0, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};
    PolarGrid grid(radii, angles);
    ASSERT_EQ(grid.numberOfNodes(), radii.size() * (angles.size() - 1));
}

TEST(PolarGridTest, AccessorsTest)
{
    std::vector<double> radii  = {0.1, 0.2, 0.5, 0.9, 1.3};
    std::vector<double> angles = {0, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};
    PolarGrid grid(radii, angles);
    ASSERT_DOUBLE_EQ(grid.radius(0), 0.1);
    ASSERT_DOUBLE_EQ(grid.radius(1), 0.2);
    ASSERT_DOUBLE_EQ(grid.radius(4), 1.3);
    ASSERT_DOUBLE_EQ(grid.theta(0), 0);
    ASSERT_DOUBLE_EQ(grid.theta(1), M_PI / 8);
    ASSERT_DOUBLE_EQ(grid.theta(4), M_PI + M_PI / 8);
}

TEST(PolarGridTest, GridJumpTest)
{
    std::vector<double> radii  = {0.1, 0.2, 0.5, 0.9, 1.3};
    std::vector<double> angles = {0, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};
    double splitting_radius    = 0.4;
    PolarGrid grid(radii, angles, splitting_radius);
    ASSERT_DOUBLE_EQ(grid.radius(0), 0.1);
    ASSERT_DOUBLE_EQ(grid.radius(1), 0.2);
    ASSERT_DOUBLE_EQ(grid.radius(4), 1.3);
    ASSERT_DOUBLE_EQ(grid.theta(0), 0);
    ASSERT_DOUBLE_EQ(grid.theta(1), M_PI / 8);
    ASSERT_DOUBLE_EQ(grid.theta(4), M_PI + M_PI / 8);
}

TEST(PolarGridTest, IndexingTest)
{
    std::vector<double> radii  = {0.1, 0.2, 0.25, 0.5, 0.8, 0.9, 1.3, 1.4, 2.0};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};
    double splitting_radius = 0.6;
    PolarGrid grid(radii, angles, splitting_radius);

    for (int i = 0; i < grid.nr(); i++) {
        for (int j = 0; j < grid.ntheta(); j++) {
            int node_index = grid.index(i, j);
            int r_out, theta_out;
            grid.multiIndex(node_index, r_out, theta_out);
            ASSERT_EQ(i, r_out);
            ASSERT_EQ(j, theta_out);
        }
    }

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
    PolarGrid grid(radii, angles, splitting_radius);

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
    PolarGrid grid(radii, angles, splitting_radius);

    ASSERT_DOUBLE_EQ(grid.radius(3), 0.5);
    ASSERT_DOUBLE_EQ(grid.theta(2), M_PI / 8);

    ASSERT_DOUBLE_EQ(grid.radius(7), 1.4);
    ASSERT_DOUBLE_EQ(grid.theta(1), M_PI / 16);
}

TEST(PolarGridTest, SpacingTest)
{
    std::vector<double> radii  = {0.1, 0.2, 0.25, 0.5, 0.8, 0.9, 1.3, 1.4, 2.0};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};
    PolarGrid grid(radii, angles);

    // Test radial spacings
    ASSERT_DOUBLE_EQ(grid.radialSpacing(0), 0.2 - 0.1);
    ASSERT_DOUBLE_EQ(grid.radialSpacing(1), 0.25 - 0.2);
    ASSERT_DOUBLE_EQ(grid.radialSpacing(2), 0.5 - 0.25);

    // Test angular spacings
    ASSERT_DOUBLE_EQ(grid.angularSpacing(0), M_PI / 16 - 0);
    ASSERT_DOUBLE_EQ(grid.angularSpacing(1), M_PI / 8 - M_PI / 16);
    ASSERT_DOUBLE_EQ(grid.angularSpacing(2), M_PI / 2 - M_PI / 8);
}
