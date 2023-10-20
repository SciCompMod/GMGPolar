#include <gtest/gtest.h>
#include "gmgpolar.h"
#include "mockgrid.h"

class test_grid : public ::testing::TestWithParam<int>
{
protected:
    test_grid()
        : test_level(0)
    {
    }
    void SetUp() override
    {
        gyro::init_params();
        gyro::icntl[Param::verbose]  = 0;
        gyro::dcntl[Param::R0]       = 1e-5;
        gyro::icntl[Param::periodic] = 1;
        gyro::f_grid_r               = "";
        gyro::f_grid_theta           = "";
        gyro::f_sol_in               = "";
        gyro::f_sol_out              = "";
        gyro::icntl[Param::fac_ani]  = 3;
        gyro::select_functions_class(gyro::icntl[Param::alpha_coeff], gyro::icntl[Param::beta_coeff],
                                     gyro::icntl[Param::mod_pk], gyro::icntl[Param::prob]);

        /*
        for (int smoother = 0; smoother < 4; smoother++) { //so that the constructor ~level() works;
            int *array_temp, *array_temp2;
            array_temp     = new int[1];
            array_temp[0]  = 0;
            array_temp2    = new int[1];
            array_temp2[0] = 0;
            test_level.dep_Asc_ortho.push_back(array_temp);
            test_level.dep_Asc.push_back(array_temp2);
            test_level.size_Asc_ortho.push_back(1);
            test_level.size_Asc.push_back(1);
        }*/
    }

    level test_level;
};

TEST_P(test_grid, Nodes_Build_r)
{
    const int& val_size            = GetParam();
    gyro::icntl[Param::nr_exp]     = (int)(val_size / 3) + 3;
    gyro::icntl[Param::ntheta_exp] = (val_size % 3) + 3;

    if (gyro::icntl[Param::nr_exp] == 3)
        gyro::icntl[Param::fac_ani] = 2;

    int& nr_exp = gyro::icntl[Param::nr_exp];

    gyro::icntl[Param::fac_ani] = 0; //uniform grid
    test_level.build_r();
    EXPECT_EQ(pow(2, nr_exp) + 1, test_level.r.size());
    for (int i = 0; i < test_level.nr; i++) {
        EXPECT_NEAR(gyro::dcntl[Param::R0] + i * (gyro::dcntl[Param::R] - gyro::dcntl[Param::R0]) / (test_level.nr - 1),
                    test_level.r[i], 1e6); //uniform grid
    }
    gyro::icntl[Param::fac_ani] = 2;
    test_level.build_r();
    EXPECT_EQ(pow(2, nr_exp + 1) + 1, test_level.r.size());
}

TEST_P(test_grid, Anisotropy_r)
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

TEST_P(test_grid, coarse_grid_one_level)
{

    const int& val_size            = GetParam();
    gyro::icntl[Param::nr_exp]     = (int)(val_size / 3) + 3;
    gyro::icntl[Param::ntheta_exp] = (val_size % 3) + 3;

    if (gyro::icntl[Param::nr_exp] == 3)
        gyro::icntl[Param::fac_ani] = 2;

    int& nr_exp = gyro::icntl[Param::nr_exp];

    gmgpolar coarse_test;
    create_grid(coarse_test);

    coarse_test.v_level[0]->define_coarse_nodes_onelevel(coarse_test.v_level[1]);

    level& fine   = *(coarse_test.v_level[0]);
    level& coarse = *(coarse_test.v_level[1]);

    EXPECT_NE(coarse.nr, fine.nr);
    EXPECT_NE(coarse.ntheta, fine.ntheta);
    EXPECT_EQ(coarse.nr, coarse.r.size());
    EXPECT_EQ(coarse.ntheta, coarse.theta.size());

    for (int j = 0; j < fine.nr; j++) {
        if (j % 2 == 0) {
            EXPECT_EQ(coarse.r[j / 2], fine.r[j]);
        }
        for (int i = 0; i < fine.ntheta_int; i++) {
            if (i % 2 == 0) {
                EXPECT_EQ(coarse.theta[i / 2], fine.theta[i]);
            }
            if (j % 2 == 0 && i % 2 == 0) {
                EXPECT_EQ(fine.coarse_nodes_list_r[j * fine.ntheta_int + i], j);
                EXPECT_EQ(fine.coarse_nodes_list_theta[j * fine.ntheta_int + i], i);
            }
            else {
                EXPECT_EQ(fine.coarse_nodes_list_r[j * fine.ntheta_int + i], -1);
                EXPECT_EQ(fine.coarse_nodes_list_theta[j * fine.ntheta_int + i], -1);
            }
        }
    }
}

TEST_P(test_grid, define_coarse_nodes)
{
    gyro::icntl[Param::nr_exp] =
        8; //test two (levels=9 or 10) instances where we can not achieve desired number of levels

    gmgpolar coarse_test;
    create_grid(coarse_test, 0);
    coarse_test.levels = GetParam() + 2;
    coarse_test.define_coarse_nodes();

    if (gyro::icntl[Param::nr_exp] < coarse_test.levels_orig) { //can not reach desired number of levels
        EXPECT_EQ(coarse_test.levels, gyro::icntl[Param::nr_exp]);
        EXPECT_EQ(coarse_test.v_level.size(), coarse_test.levels_orig);
    }
    else {
        EXPECT_EQ(coarse_test.v_level.size(), coarse_test.levels);
        EXPECT_EQ(coarse_test.levels, coarse_test.levels_orig);
    }

    for (int k = 0; k < coarse_test.levels; k++) {
        EXPECT_EQ(coarse_test.v_level[k]->nr, pow(2, gyro::icntl[Param::nr_exp] - k) + 1);
    }
}

TEST_P(test_grid, store_theta_n_co)
{

    gyro::icntl[Param::periodic] = 0;

    gmgpolar reference;
    create_grid(reference, 0);

    test_level.nr     = reference.v_level[0]->nr;
    test_level.r      = reference.v_level[0]->r; //deep copy
    test_level.ntheta = reference.v_level[0]->ntheta;
    test_level.theta  = reference.v_level[0]->theta;

    test_level.store_theta_n_co();

    EXPECT_EQ(test_level.ntheta_int, reference.v_level[0]->ntheta_int);
    EXPECT_EQ(test_level.ntheta - 1, test_level.ntheta_int);

    std::vector<level*>().swap(reference.v_level);
    gyro::icntl[Param::periodic] = 1;
    create_grid(reference, 0);

    test_level.nr     = reference.v_level[0]->nr;
    test_level.r      = reference.v_level[0]->r; //deep copy
    test_level.ntheta = reference.v_level[0]->ntheta;
    test_level.theta  = reference.v_level[0]->theta;

    test_level.store_theta_n_co();
    std::cout << reference.v_level[0]->ntheta << std::endl;
    EXPECT_EQ(test_level.ntheta_int, reference.v_level[0]->ntheta_int);
    EXPECT_EQ(test_level.ntheta, test_level.ntheta_int);
}

TEST_P(test_grid, create_grid_polar_divide)
{

    gmgpolar test_case;
    gyro::icntl[Param::periodic] = GetParam() % 2;
    create_grid(test_case, 0);
    int& divide                   = gyro::icntl[Param::divideBy2];
    int prev_nr                   = test_case.v_level[0]->nr;
    int prev_ntheta               = test_case.v_level[0]->ntheta;
    std::vector<double> r_tmp     = test_case.v_level[0]->r;
    std::vector<double> theta_tmp = test_case.v_level[0]->theta;
    divide                        = GetParam();
    test_case.create_grid_polar_divide();

    int compare_theta =
        pow(2, GetParam()) * prev_ntheta - (1 - gyro::icntl[Param::periodic]) * (pow(2, GetParam()) - 1);
    int compare_r = pow(2, GetParam()) * prev_nr - (pow(2, GetParam()) - 1);
    EXPECT_EQ(compare_r, test_case.v_level[0]->nr);
    EXPECT_EQ(compare_theta, test_case.v_level[0]->ntheta);

    int reference = (int)pow(2, divide);
    for (int k = 0; k < prev_nr - 1; k++) {
        for (int j = 0; j <= reference; j++) {
            double comp = ((double)(reference - j) / (double)reference) * r_tmp[k] +
                          ((double)j / (double)reference) * r_tmp[k + 1];
            EXPECT_NEAR(test_case.v_level[0]->r[k * reference + j], comp, 1e-6);
        }
    }
    int& per = gyro::icntl[Param::periodic];
    for (int k = 0; k < (per ? prev_ntheta : prev_ntheta - 1); k++) {
        for (int j = 0; j <= reference; j++) {
            if (k < prev_ntheta - 1) {
                double comp = ((double)(reference - j) / (double)reference) * theta_tmp[k] +
                              ((double)j / (double)reference) * theta_tmp[k + 1];
                EXPECT_NEAR(test_case.v_level[0]->theta[k * reference + j], comp, 1e-6);
            }
            else {
                if (j < reference) {
                    double comp = ((double)(reference - j) / (double)reference) * theta_tmp[k] +
                                  ((double)j / (double)reference) * 2 * PI;
                    EXPECT_NEAR(test_case.v_level[0]->theta[k * reference + j], comp, 1e-6);
                }
            }
        }
    }
}
INSTANTIATE_TEST_SUITE_P(Problem_size, test_grid, ::testing::Values(0, 1, 2, 3, 4, 5, 6, 7, 8));