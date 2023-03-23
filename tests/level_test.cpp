#include <gtest/gtest.h>
#include "gmgpolar.h"

class Level_test : public ::testing::Test
{
protected:
    Level_test()
        : test_level(0)
    {
    }
    void SetUp() override
    {
        gyro::init_params();
    }

    level test_level;
};

TEST_F(Level_test, Nodes_Build_r)
{

    gyro::icntl[Param::fac_ani] = 0;
    test_level.build_r();
    EXPECT_EQ(pow(2, gyro::icntl[Param::nr_exp]) + 1, test_level.r.size());
    gyro::icntl[Param::fac_ani] = 2;
    EXPECT_EQ(pow(2, gyro::icntl[Param::nr_exp]) + 1, test_level.r.size());
}

TEST_F(Level_test, Anisotropy_r)
{
    //testing anisotropy. refinement around 2/3 r
    ASSERT_NE(gyro::icntl[Param::alpha_coeff], 1);
    gyro::icntl[Param::fac_ani] = 0;
    while (gyro::icntl[Param::fac_ani] < gyro::icntl[Param::nr_exp]) {
        test_level.build_r();
        double i_fine   = floor(0.66 * test_level.nr);
        double h_fine   = test_level.r[i_fine + 1] - test_level.r[i_fine];
        double h_coarse = test_level.r[1] - test_level.r[0];

        EXPECT_NEAR(pow(2, gyro::icntl[Param::fac_ani]) * h_fine, h_coarse, h_fine / 10); //relative error
        std::cout << "fine mesh size " + std::to_string(h_fine) + " coarse mesh size " + std::to_string(h_coarse)
                  << std::endl;
        gyro::icntl[Param::fac_ani]++;
    }
}

TEST_F(Level_test, Nodes_Build_theta) //first learn implementation and how anisotropy is included
{

    gyro::icntl[Param::periodic] = 0;
    test_level.build_r();
    test_level.build_theta();
    int ntheta = pow(2, ceil(log2(test_level.nr))) + 1;
    EXPECT_EQ(ntheta, test_level.ntheta);

    gyro::icntl[Param::periodic] = 1;
    test_level.build_r();
    test_level.build_theta();

    ntheta = pow(2, ceil(log2(test_level.nr)));
    EXPECT_EQ(ntheta, test_level.ntheta);
}

TEST_F(Level_test, Anisotropy_theta)
{

    gyro::icntl[Param::theta_aniso] = 1;
    test_level.build_r();
    test_level.build_theta();

    double h_fine   = test_level.theta[14] - test_level.theta[13]; //based on the construction.
    double h_coarse = test_level.theta[1] - test_level.theta[0];

    EXPECT_LT(h_fine, h_coarse);
}

/*
TEST(level_test, Nodes_on_boundary) //build bound doesnt work
{
    gyro::init_params();
    level test_level(0);
    gyro::icntl[Param::verbose] = 4;

    test_level.build_r();
    test_level.build_theta();

    test_level.m = test_level.nr * test_level.ntheta;
    test_level.build_bound();

    for (int k = 0; k < test_level.m; k++) {
        std::cout << test_level.is_bound[k];
        std::cout << " ";
    }
    std::cout << "" << std::endl;
}

*/

TEST_F(Level_test, Test_Nonzeros)
{

    gyro::icntl[Param::DirBC_Interior] = 1; //no discretization across the interior
    gyro::icntl[Param::mod_pk]         = 0;
    gyro::icntl[Param::periodic]       = 1;
    gyro::icntl[Param::nr_exp]         = 6;
    gyro::icntl[Param::ntheta_exp]     = 3;
    test_level.build_r();
    test_level.build_theta();
    test_level.store_theta_n_co();
    test_level.define_nz();

    int nz_test = (test_level.nr - 2) * test_level.ntheta_int;

    EXPECT_EQ(nz_test, test_level.nz);
}

//TODO: working on testing if the matrix is build correctly. Need to understand how the values are stored.
//Especially if we discretize across the origin I am confused how we jump somehow 5*index. For tomorrow