#include <gtest/gtest.h>
#include <fstream>
#include <sched.h>
#include <time.h>
#include "cmdline.h"
#include "gmgpolar.h"

class Test_Parameters : public ::testing::Test
{
protected:
    void SetUp() override
    {
        //Results of the run with default parameters.
        input.assign({5, 0, 0, 4, 4, 3, 0, 1, 0, 3136, 49, 64, 4, 12});

        initparam              = 9;
        int initarr[initparam] = {Param::prob,   Param::alpha_coeff,    Param::beta_coeff,
                                  Param::nr_exp, Param::ntheta_exp,     Param::fac_ani,
                                  Param::mod_pk, Param::DirBC_Interior, Param::divideBy2};
        data.resize(initparam);
        std::copy(initarr, initarr + initparam, data.begin());

        //Initializing default parameters
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
    }

    int initparam;
    std::vector<int> input;
    std::vector<int> data;
};

TEST_F(Test_Parameters, Initialize_Parameters)
{
    for (int z = 0; z < initparam; z++) {
        ASSERT_EQ(gyro::icntl[data[z]], input[z]) << "Initial data is not assigned properly";
    }
}

TEST_F(Test_Parameters, DOF_on_finest_grid)
{
    gmgpolar gmgtest;
    gmgtest.create_grid_polar(); //only the finest grid is now created
    EXPECT_EQ(gmgtest.v_level.size(), 1);
    int nodes_r     = gmgtest.v_level[0]->nr;
    int nodes_theta = gmgtest.v_level[0]->ntheta;
    std::cout << "nr ntheta no prob" << std::endl;
    int finest_nodes = nodes_r * nodes_theta;
    EXPECT_EQ(finest_nodes, input[initparam]);
    std::cout << "input not assigned properly" << std::endl;
    EXPECT_EQ(nodes_r, input[initparam + 1]);
    std::cout << "initparam+1" << std::endl;
    EXPECT_EQ(nodes_theta, input[initparam + 2]);
    std::cout << "initparam+2" << std::endl;
}

TEST_F(Test_Parameters, Test_multigrid_Iterations)
{
    gmgpolar gmgtest2;
    gmgtest2.create_grid_polar();
    gmgtest2.polar_multigrid();

    EXPECT_EQ(gmgtest2.levels, input[initparam + 3]);
    int iterations =
        gyro::icntl[Param::extrapolation] < 2 ? gmgtest2.nrm_2_res.size() - 1 : gmgtest2.nrm_2_err.size() - 1;
    EXPECT_EQ(iterations, input[input.size() - 1]) << "Multigrid iterations do not match";
}