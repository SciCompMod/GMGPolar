#include <gtest/gtest.h>
#include "gmgpolar.h"

#include <gtest/gtest.h>
#include "gmgpolar.h"

#include "mockgrid.h"
class test_smoother : public ::testing::TestWithParam<int>
{
protected:
    test_smoother()
        : test_level(0)
    {
    }
    void SetUp() override
    {
        //initialize default parameters.
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

TEST_P(test_smoother, test_line_splitting)
{

    const int& val_size            = GetParam();
    gyro::icntl[Param::nr_exp]     = (int)(val_size / 3) + 3;
    gyro::icntl[Param::ntheta_exp] = (val_size % 3) + 3;

    gyro::icntl[Param::fac_ani] = (gyro::icntl[Param::nr_exp] == 3) ? 2 : 3;

    EXPECT_EQ(gyro::icntl[Param::theta_aniso], 0) << "We expect constant grid size in direction theta";

    gmgpolar instance;
    create_grid(instance, 0);
    test_level = *(instance.v_level[0]);
    test_level.define_line_splitting();

    int i = test_level.delete_circles;
    std::cout << i << std::endl;
    double k_j = test_level.theta[1] - test_level.theta[0];
    double h_i = test_level.r[i] - test_level.r[i - 1];
    EXPECT_GT(k_j * test_level.r[i], h_i) << "case 1" << std::flush;
}

TEST_P(
    test_smoother,
    mapping_usc_to_u_extrapol) //maps local indices (k,smoother) i.e the k.th coordinate from Black row etc to global one
{
    //
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

    EXPECT_EQ(test_level.nr % 2, 1);
    //test for circle smoother first
    EXPECT_EQ(test_level.ntheta % 2, 0);
    test_level.define_line_splitting();

    std::vector<int> nblocks(4);
    nblocks[0] = ceil(test_level.delete_circles * 0.5);
    nblocks[1] = floor(test_level.delete_circles * 0.5);
    nblocks[2] = test_level.ntheta_int * 0.5;
    nblocks[3] = nblocks[2];

    for (int smoother = 0; smoother < 4; ++smoother) {
        int nodesperblock;
        int only_fine_circ = extrap && (smoother == 0);
        int only_fine_rad  = extrap && (smoother == 2);
        if (only_fine_circ) {
            nodesperblock = test_level.ntheta_int / 2;
        }
        else {
            if (only_fine_rad) {
                nodesperblock = (test_level.nr - test_level.delete_circles - (test_level.delete_circles % 2 == 0)) / 2;
            }
            else {
                if (smoother < 2) {
                    nodesperblock = test_level.ntheta_int; //nodes in one circle
                }
                else {
                    nodesperblock = test_level.nr - test_level.delete_circles; //nodes in one radial
                }
            }
        }

        for (int k = 0; k < nblocks[smoother] * nodesperblock; k++) {
            if (smoother < 2) {
                int glob_row = 2 * (k / nodesperblock) +
                               smoother; //white (smoother=1) cirle uneven lines. Times two since zebra coloring
                int glob_col = only_fine_circ ? (2 * (k % nodesperblock) + 1) : k % nodesperblock;

                int glob_ind = glob_row * test_level.ntheta_int + glob_col;

                EXPECT_EQ(test_level.mapping_usc_to_u(k, smoother), glob_ind)
                    << "fails for smoother " + std::to_string(smoother) + " index " + std::to_string(k) +
                           " ntheta_int" + std::to_string(test_level.ntheta_int);
            }
            else {
                int glob_row = only_fine_rad ? (2 * (k % nodesperblock) + (test_level.delete_circles % 2 == 0))
                                             : k % nodesperblock;
                int glob_col = 2 * (k / nodesperblock) + (smoother == 3);

                int glob_ind = (glob_row + test_level.delete_circles) * test_level.ntheta_int + glob_col;

                EXPECT_EQ(test_level.mapping_usc_to_u(k, smoother), glob_ind)
                    << "fails for smoother " + std::to_string(smoother) + " index " + std::to_string(k) +
                           " ntheta_int " + std::to_string(test_level.ntheta_int) + "delete_circ " +
                           std::to_string(test_level.delete_circles);
            }
        }
    }
}

TEST_P(test_smoother, mapping_usc_to_u_no_extrapol)
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

    EXPECT_EQ(test_level.nr % 2, 1);
    //test for circle smoother first
    EXPECT_EQ(test_level.ntheta % 2, 0);
    test_level.define_line_splitting();

    std::vector<int> nblocks(4);
    nblocks[0] = ceil(test_level.delete_circles * 0.5);
    nblocks[1] = floor(test_level.delete_circles * 0.5);
    nblocks[2] = test_level.ntheta_int * 0.5;
    nblocks[3] = nblocks[2];

    for (int smoother = 0; smoother < 4; ++smoother) {
        int nodesperblock;
        if (smoother < 2) {
            nodesperblock = test_level.ntheta_int;
        }
        else {
            nodesperblock = test_level.nr - test_level.delete_circles;
        }

        for (int k = 0; k < nblocks[smoother] * nodesperblock; k++) {
            if (smoother < 2) {
                int glob_row = 2 * (k / nodesperblock) +
                               smoother; //white (smoother=1) cirle uneven lines. Times two since zebra coloring
                int glob_col = k % nodesperblock;

                int glob_ind = glob_row * test_level.ntheta_int + glob_col;

                EXPECT_EQ(test_level.mapping_usc_to_u(k, smoother), glob_ind)
                    << "fails for smoother " + std::to_string(smoother) + " index " + std::to_string(k) +
                           " ntheta_int" + std::to_string(test_level.ntheta_int);
            }
            else {
                int glob_row = k % nodesperblock;
                int glob_col = 2 * (k / nodesperblock) + (smoother == 3);

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

/*
TEST_P(test_smoother, test_get_row)
{
    test_level.build_r();
    test_level.build_theta();
    test_level.store_theta_n_co();

    test_level.define_line_splitting();

    int& extrap = gyro::icntl[Param::extrapolation]; //level is 0 already

    extrap = 1;

    std::vector<int> nblocks(4);
    nblocks[0] = ceil(test_level.delete_circles * 0.5);
    nblocks[1] = floor(test_level.delete_circles * 0.5);
    nblocks[2] = test_level.ntheta_int * 0.5;
    nblocks[3] = nblocks[2];

    for (int j = 0; j < test_level.nr; j++) {
        int smoother;
        if (!(j % 2) && j < test_level.delete_circles) {
            smoother = 0;
        }
        else {
            if (j < test_level.delete_circles) {
                smoother = 1;
            }
            else {
                smoother = 3;
            }
        }
        std::vector<int> res = test_level.get_row(j, smoother, extrap, 1, 1);
        if (smoother < 2) {
            if (extrap == 0) {
                for (int k = 0; k < test_level.ntheta_int; k++) {
                    EXPECT_EQ(res[k], k);
                }
            }
            else {
                if (smoother == 0) {
                    for (int k = 1; k < test_level.ntheta_int; k += 2) {
                        EXPECT_EQ(res[k], (k - 1) / 2);
                    }
                }
                else {
                    for (int k = 0; k < test_level.ntheta_int; k++) {
                        EXPECT_EQ(res[k], k);
                    }
                }
            }
        }
        else {
            if (extrap == 0) {
                for (int k = 0; k < test_level.ntheta_int; k++) {
                    EXPECT_EQ(res[k], j - test_level.delete_circles);
                }
            }
            else {
                if (test_level.delete_circles % 2) {
                    int jump = !((j - test_level.delete_circles) % 2);
                    for (int k = jump; k < test_level.ntheta_int; k += (jump + 1)) {
                        int eq = (jump == 0 && !(k % 2)) ? (j - test_level.delete_circles) / 2
                                                         : j - test_level.delete_circles;
                        EXPECT_EQ(res[k], eq);
                    }
                }
                else {
                    int jump = (j - test_level.delete_circles) % 2;
                    for (int k = jump; k < test_level.ntheta_int; k += (jump + 1)) {
                        int eq = (jump == 0 && !(k % 2)) ? (j - test_level.delete_circles) / 2
                                                         : j - test_level.delete_circles;
                        EXPECT_EQ(res[k], eq);
                    }
                }
            }
        }
    }
}
*/

/*
TEST_F(Test_smoother, Test_build_Asc)
{

   
    level test_level = new level(0);

    test_level.define_line_splitting() ;
    test_level.nblocks = std::vector<int>(5);
    test_level.nblocks[0] = ceil(v_level[l]->delete_circles * 0.5);
    test_level.nblocks[1] = floor(v_level[l]->delete_circles * 0.5);
    test_level.nblocks[2] = v_level[l]->ntheta_int * 0.5;
    test_level.nblocks[3] = v_level[l]->nblocks[2];
    test_level.nblocks[4] = std::max(v_level[l]->nblocks[0], v_level[l]->nblocks[2]);

    for (int smoother = 0; smoother < 4; smoother++) {
                int size = 0, size_ortho = 0, *array_temp, *array_temp2;
                if (smoother < 2) {
                    size_ortho = test_level.delete_circles + 1;
                    size       = test_level.delete_circles;
                }
                else if (smoother > 1) {
                    size_ortho = test_level.ntheta_int;
                    size       = test_level.ntheta_int;
                }
                array_temp = new int[size_ortho];
                for (int i = 0; i < size_ortho; i++)
                    array_temp[i] = 0;
                array_temp2 = new int[size];
                for (int i = 0; i < size; i++)
                    array_temp2[i] = 0;
                test_level.dep_Asc_ortho.push_back(array_temp);
                test_level.dep_Asc.push_back(array_temp2);
                test_level.size_Asc_ortho.push_back(size_ortho);
                test_level.size_Asc.push_back(size);
            }
}
*/

INSTANTIATE_TEST_SUITE_P(Problem_size, test_smoother, ::testing::Values(0, 1, 2, 3, 4, 5, 6, 7, 8));