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
            MultiIndex alpha(i, j);
            int node_index  = grid.index(alpha);
            MultiIndex beta = grid.multiIndex(node_index);
            ASSERT_EQ(alpha[0], beta[0]);
            ASSERT_EQ(alpha[1], beta[1]);
        }
    }

    for (int i = 0; i < grid.numberOfNodes(); i++) {
        MultiIndex alpha = grid.multiIndex(i);
        int node_index   = grid.index(alpha);
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
        MultiIndex alpha(2, 6);
        int node_index = grid.index(alpha);
        ASSERT_EQ(node_index, 22);
        MultiIndex beta = grid.multiIndex(node_index);
        ASSERT_EQ(alpha[0], beta[0]);
        ASSERT_EQ(alpha[1], beta[1]);
    }

    {
        MultiIndex alpha(3, 2);
        int node_index = grid.index(alpha);
        ASSERT_EQ(node_index, 26);
        MultiIndex beta = grid.multiIndex(node_index);
        ASSERT_EQ(alpha[0], beta[0]);
        ASSERT_EQ(alpha[1], beta[1]);
    }

    {
        MultiIndex alpha(6, 4);
        int node_index = grid.index(alpha);
        ASSERT_EQ(node_index, 54);
        MultiIndex beta = grid.multiIndex(node_index);
        ASSERT_EQ(alpha[0], beta[0]);
        ASSERT_EQ(alpha[1], beta[1]);
    }

    {
        MultiIndex alpha(4, 7);
        int node_index = grid.index(alpha);
        ASSERT_EQ(node_index, 67);
        MultiIndex beta = grid.multiIndex(node_index);
        ASSERT_EQ(alpha[0], beta[0]);
        ASSERT_EQ(alpha[1], beta[1]);
    }
}

TEST(PolarGridTest, CoordinatesTest)
{
    std::vector<double> radii  = {0.1, 0.2, 0.25, 0.5, 0.8, 0.9, 1.3, 1.4, 2.0};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};
    double splitting_radius = 0.6;
    PolarGrid grid(radii, angles, splitting_radius);

    MultiIndex alpha(3, 2);
    Point P1 = grid.polarCoordinates(alpha);
    ASSERT_DOUBLE_EQ(P1[0], 0.5);
    ASSERT_DOUBLE_EQ(P1[1], M_PI / 8);

    alpha[0] += 4;
    alpha[1] -= 1;
    P1 = grid.polarCoordinates(alpha);
    ASSERT_DOUBLE_EQ(P1[0], 1.4);
    ASSERT_DOUBLE_EQ(P1[1], M_PI / 16);
}

TEST(PolarGridTest, NeighborDistanceTest)
{
    std::vector<double> radii  = {0.1, 0.2, 0.25, 0.5, 0.8, 0.9, 1.3, 1.4, 2.0};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};
    double splitting_radius = 0.6;
    PolarGrid grid(radii, angles, splitting_radius);

    std::array<std::pair<double, double>, space_dimension> neighbor_distance;
    {
        MultiIndex alpha(3, 2);
        grid.adjacentNeighborDistances(alpha, neighbor_distance);
        ASSERT_DOUBLE_EQ(neighbor_distance[0].first, 0.5 - 0.25);
        ASSERT_DOUBLE_EQ(neighbor_distance[0].second, 0.8 - 0.5);
        ASSERT_DOUBLE_EQ(neighbor_distance[1].first, M_PI / 8 - M_PI / 16);
        ASSERT_DOUBLE_EQ(neighbor_distance[1].second, M_PI / 2 - M_PI / 8);
    }

    {
        MultiIndex alpha(8, 7);
        grid.adjacentNeighborDistances(alpha, neighbor_distance);
        ASSERT_DOUBLE_EQ(neighbor_distance[0].first, 2.0 - 1.4);
        ASSERT_DOUBLE_EQ(neighbor_distance[0].second, 0.0);
        ASSERT_DOUBLE_EQ(neighbor_distance[1].first, M_PI / 2 - M_PI / 8);
        ASSERT_DOUBLE_EQ(neighbor_distance[1].second, M_PI / 2);
    }

    {
        MultiIndex alpha(4, 0);
        grid.adjacentNeighborDistances(alpha, neighbor_distance);
        ASSERT_DOUBLE_EQ(neighbor_distance[0].first, 0.8 - 0.5);
        ASSERT_DOUBLE_EQ(neighbor_distance[0].second, 0.9 - 0.8);
        ASSERT_DOUBLE_EQ(neighbor_distance[1].first, M_PI / 2);
        ASSERT_DOUBLE_EQ(neighbor_distance[1].second, M_PI / 16);
    }

    {
        MultiIndex alpha(0, 5);
        grid.adjacentNeighborDistances(alpha, neighbor_distance);
        ASSERT_DOUBLE_EQ(neighbor_distance[0].first, 0.0);
        ASSERT_DOUBLE_EQ(neighbor_distance[0].second, 0.2 - 0.1);
        ASSERT_NEAR(neighbor_distance[1].first, M_PI / 16, 1e-15);
        ASSERT_NEAR(neighbor_distance[1].second, M_PI / 8 - M_PI / 16, 1e-15);
        ASSERT_EQ(equals(neighbor_distance[1].first, M_PI / 16), true);
        ASSERT_EQ(equals(neighbor_distance[1].second, M_PI / 8 - M_PI / 16), true);
    }
}

TEST(PolarGridTest, NeighborsTest)
{
    std::vector<double> radii  = {0.1, 0.2, 0.25, 0.5, 0.8, 0.9, 1.3, 1.4, 2.0};
    std::vector<double> angles = {
        0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};
    double splitting_radius = 0.6;
    PolarGrid grid(radii, angles, splitting_radius);

    std::array<std::pair<int, int>, space_dimension> neighbors;

    {
        MultiIndex alpha(3, 2);
        grid.adjacentNeighborsOf(alpha, neighbors);
        ASSERT_EQ(neighbors[0].first, 18);
        ASSERT_EQ(neighbors[0].second, 42);
        ASSERT_EQ(neighbors[1].first, 25);
        ASSERT_EQ(neighbors[1].second, 27);
        grid.diagonalNeighborsOf(alpha, neighbors);
        ASSERT_EQ(neighbors[0].first, 17);
        ASSERT_EQ(neighbors[0].second, 37);
        ASSERT_EQ(neighbors[1].first, 19);
        ASSERT_EQ(neighbors[1].second, 47);
    }

    {
        MultiIndex alpha(8, 7);
        grid.adjacentNeighborsOf(alpha, neighbors);
        ASSERT_EQ(neighbors[0].first, 70);
        ASSERT_EQ(neighbors[0].second, -1);
        ASSERT_EQ(neighbors[1].first, 66);
        ASSERT_EQ(neighbors[1].second, 36);
        grid.diagonalNeighborsOf(alpha, neighbors);
        ASSERT_EQ(neighbors[0].first, 65);
        ASSERT_EQ(neighbors[0].second, -1);
        ASSERT_EQ(neighbors[1].first, 35);
        ASSERT_EQ(neighbors[1].second, -1);
    }

    {
        MultiIndex alpha(4, 0);
        grid.adjacentNeighborsOf(alpha, neighbors);
        ASSERT_EQ(neighbors[0].first, 24);
        ASSERT_EQ(neighbors[0].second, 33);
        ASSERT_EQ(neighbors[1].first, 67);
        ASSERT_EQ(neighbors[1].second, 37);
        grid.diagonalNeighborsOf(alpha, neighbors);
        ASSERT_EQ(neighbors[0].first, 31);
        ASSERT_EQ(neighbors[0].second, 68);
        ASSERT_EQ(neighbors[1].first, 25);
        ASSERT_EQ(neighbors[1].second, 38);
    }

    {
        MultiIndex alpha(0, 5);
        grid.adjacentNeighborsOf(alpha, neighbors);
        ASSERT_EQ(neighbors[0].first, -1);
        ASSERT_EQ(neighbors[0].second, 13);
        ASSERT_EQ(neighbors[1].first, 4);
        ASSERT_EQ(neighbors[1].second, 6);
        grid.diagonalNeighborsOf(alpha, neighbors);
        ASSERT_EQ(neighbors[0].first, -1);
        ASSERT_EQ(neighbors[0].second, 12);
        ASSERT_EQ(neighbors[1].first, -1);
        ASSERT_EQ(neighbors[1].second, 14);
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}