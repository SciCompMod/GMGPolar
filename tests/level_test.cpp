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
    int& geom    = gyro::icntl[Param::mod_pk];
    int& dir_int = gyro::icntl[Param::DirBC_Interior];

    dir_int = 1;

    test_level.build_r();
    test_level.build_theta();
    test_level.store_theta_n_co();
    int m = test_level.nr * test_level.ntheta;

    //z = r*ntheta_int + theta
    EXPECT_EQ(0, test_level.get_ptr(0, 0));

    EXPECT_EQ(test_level.ntheta_int - 1, test_level.get_ptr(test_level.ntheta_int - 1, 0));
    for (int r = 0; r < test_level.nr; r++) {
        EXPECT_EQ(test_level.get_ptr(0, r), test_level.get_ptr(test_level.ntheta_int, r)); //periodic
    }
}

TEST_F(Level_test, Test_get_stencil)
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

TEST_F(Level_test, Test_prolongation)
{

    gmgpolar test_prolon;
    test_prolon.create_grid_polar();
    test_prolon.check_geom();
    test_prolon.define_coarse_nodes();

    level& prolon_level = *(test_prolon.v_level[0]);
    int ctheta_int      = test_prolon.v_level[1]->ntheta_int;
    int cr_int          = test_prolon.v_level[1]->nr_int;

    prolon_level.m  = test_prolon.v_level[0]->nr * test_prolon.v_level[0]->ntheta;
    prolon_level.mc = test_prolon.v_level[1]->nr * test_prolon.v_level[1]->ntheta;

    std::vector<double> u_test(prolon_level.mc);
    for (int z = 0; z < u_test.size(); z++) {
        u_test[z] = z;
    }

    std::vector<double> sol =
        prolon_level.apply_prolongation_bi0(u_test, prolon_level.mc, prolon_level.m, prolon_level.coarse_nodes_list_r,
                                            prolon_level.coarse_nodes_list_theta, 0);

    for (int j = 0; j < prolon_level.nr_int + 1; j++) {
        for (int i = 0; i < prolon_level.ntheta_int; i++) {
            if (j % 2 == 0 && i % 2 == 0) {
                EXPECT_EQ(u_test[(j / 2) * ctheta_int + (i / 2)], sol[j * prolon_level.ntheta_int + i])
                    << "coarse node injection is failing";
            }
            else if (i % 2 != 0 && j % 2 == 0) { // theta_i is fine node
                double k_qm1 = prolon_level.theta[i] - prolon_level.theta[i - 1];
                double k_q =
                    (i < prolon_level.ntheta_int - 1) ? prolon_level.theta[i + 1] - prolon_level.theta[i] : k_qm1;

                int i1     = (j / 2) * ctheta_int + (i - 1) / 2;
                int i2     = i1 + 1;
                double val = 1 / (k_q + k_qm1) * (k_q * u_test[i1] + k_qm1 * u_test[i2]);

                EXPECT_NEAR(val, sol[j * prolon_level.ntheta_int + i], 1e-8);
            }
        }
    }
}

//TODO: working on testing if the matrix is build correctly. Need to understand how the values are stored.
//Especially if we discretize across the origin I am confused how we jump somehow 5*index. For tomorrow

//TODO: Last test fails for i== ntheta_int -1. why ?