#include <gtest/gtest.h>
#include "gmgpolar.h"
#include "mockgrid.h"

class test_operator : public testing::Test
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

        create_grid(test_solver, 0);

        test_level.nr         = reference.v_level[0]->nr;
        test_level.r          = reference.v_level[0]->r;
        test_level.ntheta     = reference.v_level[0]->ntheta;
        test_level.theta      = reference.v_level[0]->theta;
        test_level.ntheta_int = reference.v_level[0]->ntheta_int;
        test_level.nr_int     = reference.v_level[0]->nr_int;
        test_level.thetaplus  = reference.v_level[0]->thetaplus;
        test_level.hplus      = reference.v_level[0]->hplus;
    }

    gmgpolar test_solver;
    level test_level;
};

TEST_G(test_operator, test_stencil)
{
}
TEST_F(test_operator, build_A)
{
}
