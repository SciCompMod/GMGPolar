#include <gtest/gtest.h>
#include "gmgpolar.h"

class Level_test : public ::testing::TestWithParam<int>
{
protected:
    Level_test()
        : test_level(0)
    {
    }
    void SetUp() override
    {
        gyro::init_params();
        gyro::icntl[Param::verbose]        = 0;
        gyro::icntl[Param::debug]          = 0;
        gyro::icntl[Param::extrapolation]  = 0;
        gyro::icntl[Param::DirBC_Interior] = 1;
        gyro::icntl[Param::check_error]    = 1;
        gyro::dcntl[Param::R0]             = 1e-5;
        gyro::f_grid_r                     = "";
        gyro::f_grid_theta                 = "";
        gyro::f_sol_in                     = "";
        gyro::f_sol_out                    = "";
        gyro::icntl[Param::nr_exp]         = 4;
        gyro::icntl[Param::ntheta_exp]     = 4;
        gyro::icntl[Param::fac_ani]        = 3;
        gyro::select_functions_class(gyro::icntl[Param::alpha_coeff], gyro::icntl[Param::beta_coeff],
                                     gyro::icntl[Param::mod_pk], gyro::icntl[Param::prob]);

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
        }
    }

    level test_level;
};

TEST_P(Level_test, Nodes_Build_r)
{
    const int& val_size            = GetParam();
    gyro::icntl[Param::nr_exp]     = (int)(val_size / 3) + 3;
    gyro::icntl[Param::ntheta_exp] = (val_size % 3) + 3;

    if (gyro::icntl[Param::nr_exp] == 3)
        gyro::icntl[Param::fac_ani] = 2;

    int& nr_exp = gyro::icntl[Param::nr_exp];

    gyro::icntl[Param::fac_ani] = 0;
    test_level.build_r();
    EXPECT_EQ(pow(2, nr_exp) + 1, test_level.r.size());
    gyro::icntl[Param::fac_ani] = 2;
    test_level.build_r();
    EXPECT_EQ(pow(2, nr_exp + 1) + 1, test_level.r.size());
}

TEST_P(Level_test, Anisotropy_r)
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

TEST_P(Level_test, Nodes_Build_theta) //first learn implementation and how anisotropy is included
{
    const int& val_size            = GetParam();
    gyro::icntl[Param::nr_exp]     = (int)(val_size / 3) + 3;
    gyro::icntl[Param::ntheta_exp] = (val_size % 3) + 3;

    if (gyro::icntl[Param::nr_exp] == 3)
        gyro::icntl[Param::fac_ani] = 2;

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

//build bound function is unfinished

TEST_P(Level_test, Test_Nonzero_size)
{
    const int& val_size            = GetParam();
    gyro::icntl[Param::nr_exp]     = (int)(val_size / 3) + 3;
    gyro::icntl[Param::ntheta_exp] = (val_size % 3) + 3;

    if (gyro::icntl[Param::nr_exp] == 3)
        gyro::icntl[Param::fac_ani] = 2;

    gyro::icntl[Param::periodic] = 1;

    int& geom    = gyro::icntl[Param::mod_pk];
    int& dir_int = gyro::icntl[Param::DirBC_Interior];

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

TEST_P(Level_test, Test_get_ptr)
{
    const int& val_size            = GetParam();
    gyro::icntl[Param::nr_exp]     = (int)(val_size / 3) + 3;
    gyro::icntl[Param::ntheta_exp] = (val_size % 3) + 3;

    if (gyro::icntl[Param::nr_exp] == 3)
        gyro::icntl[Param::fac_ani] = 2;

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

TEST_P(Level_test, Test_get_stencil) //unfinished
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

TEST_P(Level_test, mapping_usc_to_u_extrapol)
{

    const int& val_size            = GetParam();
    gyro::icntl[Param::nr_exp]     = (int)(val_size / 3) + 3;
    gyro::icntl[Param::ntheta_exp] = (val_size % 3) + 3;

    if (gyro::icntl[Param::nr_exp] == 3)
        gyro::icntl[Param::fac_ani] = 2;

    int& extrap   = gyro::icntl[Param::extrapolation];
    extrap        = 1;
    int& periodic = gyro::icntl[Param::periodic];
    periodic      = 1;

    test_level.build_r();
    test_level.build_theta();
    test_level.store_theta_n_co();

    EXPECT_EQ(test_level.nr % 2 == 1, 1);
    //test for circle smoother first

    test_level.define_line_splitting();

    std::vector<int> nblocks(4);
    nblocks[0] = ceil(test_level.delete_circles * 0.5);
    nblocks[1] = floor(test_level.delete_circles * 0.5);
    nblocks[2] = test_level.ntheta_int * 0.5;
    nblocks[3] = nblocks[2];

    for (int smoother = 0; smoother < 4; ++smoother) {
        int ntheta_smoothed;
        int only_fine = extrap && smoother == 0;
        if (only_fine) {
            ntheta_smoothed = test_level.ntheta_int / 2;
        }
        else {
            ntheta_smoothed = test_level.ntheta_int;
        }

        for (int k = 0; k < nblocks[smoother] * ntheta_smoothed; k++) {
            if (smoother < 2) {
                int glob_row = 2 * (k / ntheta_smoothed) + (smoother == 1); //white cirle uneven lines
                int glob_col = only_fine ? (2 * (k % ntheta_smoothed) + 1) : k % ntheta_smoothed;

                int glob_ind = glob_row * test_level.ntheta_int + glob_col;

                EXPECT_EQ(test_level.mapping_usc_to_u(k, smoother), glob_ind)
                    << "fails for smoother " + std::to_string(smoother) + " index " + std::to_string(k) +
                           " ntheta_int" + std::to_string(test_level.ntheta_int);
            }
            else {
                int nr_smoothed =
                    (extrap && smoother == 2)
                        ? (test_level.nr - test_level.delete_circles + (test_level.delete_circles % 2 == 0)) / 2
                        : test_level.nr - test_level.delete_circles;
                int glob_row = (extrap && smoother == 2)
                                   ? (2 * (k % nr_smoothed) + (test_level.delete_circles % 2 == 0))
                                   : k % nr_smoothed;
                int glob_col = 2 * (k / nr_smoothed) + (smoother == 3);

                int glob_ind = (glob_row + test_level.delete_circles) * test_level.ntheta_int + glob_col;

                EXPECT_EQ(test_level.mapping_usc_to_u(k, smoother), glob_ind)
                    << "fails for smoother " + std::to_string(smoother) + " index " + std::to_string(k) +
                           " ntheta_int " + std::to_string(test_level.ntheta_int) + "delete_circ " +
                           std::to_string(test_level.delete_circles);
            }
        }
        //TODO: This test just rewrites the function.. maybe test some szenarios better
    }
}

TEST_P(Level_test, mapping_usc_to_u_noextrapol)
{

    const int& val_size            = GetParam();
    gyro::icntl[Param::nr_exp]     = (int)(val_size / 3) + 3;
    gyro::icntl[Param::ntheta_exp] = (val_size % 3) + 3;

    if (gyro::icntl[Param::nr_exp] == 3)
        gyro::icntl[Param::fac_ani] = 2;

    int& extrap   = gyro::icntl[Param::extrapolation];
    extrap        = 0;
    int& periodic = gyro::icntl[Param::periodic];
    periodic      = 1;

    test_level.build_r();
    test_level.build_theta();
    test_level.store_theta_n_co();

    EXPECT_EQ(test_level.nr % 2 == 1, 1);
    //test for circle smoother first

    test_level.define_line_splitting();

    std::vector<int> nblocks(4);
    nblocks[0] = ceil(test_level.delete_circles * 0.5);
    nblocks[1] = floor(test_level.delete_circles * 0.5);
    nblocks[2] = test_level.ntheta_int * 0.5;
    nblocks[3] = nblocks[2];

    for (int smoother = 0; smoother < 4; ++smoother) {
        int ntheta_smoothed;
        int only_fine = extrap && smoother == 0;
        if (only_fine) {
            ntheta_smoothed = test_level.ntheta_int / 2;
        }
        else {
            ntheta_smoothed = test_level.ntheta_int;
        }

        for (int k = 0; k < nblocks[smoother] * ntheta_smoothed; k++) {
            if (smoother < 2) {
                int glob_row = 2 * (k / ntheta_smoothed) + (smoother == 1); //white cirle uneven lines
                int glob_col = only_fine ? (2 * (k % ntheta_smoothed) + 1) : k % ntheta_smoothed;

                int glob_ind = glob_row * test_level.ntheta_int + glob_col;

                EXPECT_EQ(test_level.mapping_usc_to_u(k, smoother), glob_ind)
                    << "fails for smoother " + std::to_string(smoother) + " index " + std::to_string(k) +
                           " ntheta_int" + std::to_string(test_level.ntheta_int);
            }
            else {
                int nr_smoothed =
                    (extrap && smoother == 2)
                        ? (test_level.nr - test_level.delete_circles + (test_level.delete_circles % 2 == 0)) / 2
                        : test_level.nr - test_level.delete_circles;
                int glob_row = (extrap && smoother == 2)
                                   ? (2 * (k % nr_smoothed) + (test_level.delete_circles % 2 == 0))
                                   : k % nr_smoothed;
                int glob_col = 2 * (k / nr_smoothed) + (smoother == 3);

                int glob_ind = (glob_row + test_level.delete_circles) * test_level.ntheta_int + glob_col;

                EXPECT_EQ(test_level.mapping_usc_to_u(k, smoother), glob_ind)
                    << "fails for smoother " + std::to_string(smoother) + " index " + std::to_string(k) +
                           " ntheta_int " + std::to_string(test_level.ntheta_int) + "delete_circ " +
                           std::to_string(test_level.delete_circles);
            }
        }
        //TODO: This test just rewrites the function.. maybe test some szenarios better
    }
}

TEST_P(Level_test, mapping_usc_to_u_array)
{
}

INSTANTIATE_TEST_SUITE_P(Level_size, Level_test, ::testing::Values(0, 1, 2, 3, 4, 5, 6, 7, 8));