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

    int& nr_exp = gyro::icntl[Param::nr_exp];

    for (nr_exp = 3; nr_exp < 6; nr_exp++) {
        gyro::icntl[Param::fac_ani] = 0;
        test_level.build_r();
        EXPECT_EQ(pow(2, nr_exp) + 1, test_level.r.size());
        gyro::icntl[Param::fac_ani] = 2;
        test_level.build_r();
        EXPECT_EQ(pow(2, nr_exp + 1) + 1, test_level.r.size());
    }
}

TEST_F(Level_test, Anisotropy_r)
{
    //testing anisotropy. refinement around 2/3 r
    ASSERT_NE(gyro::icntl[Param::alpha_coeff], 1);
    int& ani_r = gyro::icntl[Param::fac_ani];
    ani_r      = 0;
    while (ani_r < gyro::icntl[Param::nr_exp]) {
        test_level.build_r();
        double i_fine   = floor(0.66 * test_level.nr);
        double h_fine   = test_level.r[i_fine + 1] - test_level.r[i_fine];
        double h_coarse = test_level.r[1] - test_level.r[0];

        EXPECT_NEAR(pow(2, ani_r) * h_fine, h_coarse, h_fine / 10); //relative error
        std::cout << "fine mesh size " + std::to_string(h_fine) + " coarse mesh size " + std::to_string(h_coarse)
                  << std::endl;
        ani_r++;
    }
}

TEST_F(Level_test, Nodes_Build_theta) //first learn implementation and how anisotropy is included
{
    int& ntheta_exp = gyro::icntl[Param::ntheta_exp];
    for (ntheta_exp = 3; ntheta_exp < 6; ntheta_exp++) {
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
}
//todo
TEST_F(Level_test, Anisotropy_theta) //unfinished
{

    gyro::icntl[Param::theta_aniso] = 1;
    test_level.build_r();
    test_level.build_theta();

    double h_fine   = test_level.theta[14] - test_level.theta[13]; //based on the construction.
    double h_coarse = test_level.theta[1] - test_level.theta[0];

    EXPECT_LT(h_fine, h_coarse);
}

//build bound function is unfinished

TEST_F(Level_test, Test_Nonzero_size)
{

    gyro::icntl[Param::periodic] = 1;

    int& geom    = gyro::icntl[Param::mod_pk];
    int& dir_int = gyro::icntl[Param::DirBC_Interior];
    for (int z = 0; z < 16; z++) { //test for different mesh sizes
        gyro::icntl[Param::nr_exp]     = (int)(z / 4) + 3;
        gyro::icntl[Param::ntheta_exp] = z % 4 + 3;

        geom = 0; //5-point stencil

        int dof = 0;
        for (dir_int = 0; dir_int < 2; dir_int++) {

            test_level.define_nz();

            int bd = dir_int + 1;
            dof    = (test_level.nr - bd) * test_level.ntheta_int;
            EXPECT_EQ(5 * dof, test_level.nz);

            /*Every dof is found five times in the entry.
             boundary nodes once and neighbors to DB nodes 4 times*/
        }

        geom = 1; //9 point stencil

        for (dir_int = 0; dir_int < 2; dir_int++) {

            test_level.define_nz();

            int bd      = dir_int + 1;
            dof         = (test_level.nr - bd) * test_level.ntheta_int;
            int entries = 9 * dof - 2 * (bd + (1 - dir_int)) * test_level.ntheta_int;

            EXPECT_EQ(entries, test_level.nz);

            /*DB_boundary neighbors 6 stencil points +1 for Db nodes. and 7 point stencil for disc across 
                the interior*/
        }
    }
}

TEST_F(Level_test, Test_get_ptr)
{
    int& dir_int = gyro::icntl[Param::DirBC_Interior];

    dir_int = 1;

    test_level.build_r();
    test_level.build_theta();
    test_level.store_theta_n_co();

    //z = r*ntheta_int + theta
    EXPECT_EQ(0, test_level.get_ptr(0, 0));

    EXPECT_EQ(test_level.ntheta_int - 1, test_level.get_ptr(test_level.ntheta_int - 1, 0));
    for (int r = 0; r < test_level.nr; r++) {
        EXPECT_EQ(test_level.get_ptr(0, r), test_level.get_ptr(test_level.ntheta_int, r)); //periodic
    }
}

TEST_F(Level_test, Test_get_stencil) //unfinished
{
    test_level.build_r();
    int& geom = gyro::icntl[Param::mod_pk];
    for (geom = 0; geom < 3; geom++) {
        for (int r = 0; r < test_level.nr; r++) {
            std::vector<int> s = test_level.get_stencil(r);
            int val            = *std::max_element(s.begin(), s.end());
            EXPECT_LE(val, (geom == 0) ? 4 : 8);
            EXPECT_LE(0, val);
        }
    }
}