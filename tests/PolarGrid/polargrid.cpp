#include <gtest/gtest.h>
#include "../../include/PolarGrid/polargrid.h"

#include <cmath>
#include <vector>

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

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}