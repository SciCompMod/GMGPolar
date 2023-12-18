#include <gtest/gtest.h>
#include "gmgpolar.h"
#include "mockgrid.h"

class test_operator : public testing::Test
{
protected:
    test_operator()
        : test_level(0)
    {
    }
    void SetUp() override
    {
        gyro::init_params();
        gyro::icntl[Param::nr_exp]   = 7;
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
        std::cout << "1" << std::endl;
        create_grid(test_solver, 0);
        std::cout << "1" << std::endl;
        test_level.nr         = test_solver.v_level[0]->nr;
        test_level.r          = test_solver.v_level[0]->r;
        test_level.ntheta     = test_solver.v_level[0]->ntheta;
        test_level.theta      = test_solver.v_level[0]->theta;
        test_level.ntheta_int = test_solver.v_level[0]->ntheta_int;
        test_level.nr_int     = test_solver.v_level[0]->nr_int;
        test_level.thetaplus  = test_solver.v_level[0]->thetaplus;
        test_level.hplus      = test_solver.v_level[0]->hplus;
        std::cout << "2" << std::endl;
    }

    gmgpolar test_solver;
    level test_level;
};

int set_nonzero(int mod_pk, int nr, int ntheta_int)
{
    int nz             = 0;
    int stencil_points = mod_pk == 0 ? 5 : 9;
    //todo: differentiate between dirichlet, discretization across interior, no dirichlet
    int dirichlets = (gyro::icntl[Param::DirBC_Interior]) ? 2 : 1;
    nz += dirichlets * ntheta_int;
    if ((dirichlets == 2 && nr >= 4) || (dirichlets == 1 && nr >= 3)) {
        int linked_nodes = 2 * (mod_pk != 0) + 1; //1 in case circular, 3 else
        nz += (stencil_points - linked_nodes) * ntheta_int * dirichlets; //nodes linked to dirichlet
    }
    else {
        if ((dirichlets == 2 && nr == 3)) {
            nz += 3 * ntheta_int;
        }
    }
    int interior_unlinked_inr = std::max(nr - (2 * dirichlets), 0);
    nz += interior_unlinked_inr * ntheta_int * stencil_points;

    return nz;
    //if(!gyro::icntl[Param::DirBC_Interior])
}

TEST_F(test_operator, test_nonzero)
{

    gyro::icntl[Param::DirBC_Interior] = 1;
    int& ntheta_int                    = test_level.ntheta_int;

    test_level.define_nz();

    EXPECT_EQ(test_level.nz, set_nonzero(gyro::icntl[Param::mod_pk], test_level.nr, ntheta_int));
}

TEST_F(test_operator, build_A)
{
    gyro::icntl[Param::DirBC_Interior] = 1;
    int& ntheta_int                    = test_level.ntheta_int;
    int& nr_int                        = test_level.nr_int;
    int nz                             = set_nonzero(gyro::icntl[Param::mod_pk], test_level.nr, ntheta_int);
    test_level.build_A();

    int size = test_level.vals.size();
    std::cout << size << std::endl;
    std::cout << test_level.row_indices.size() << std::endl;
    std::cout << test_level.col_indices.size() << std ::endl;
    int gridp = (nr_int + 1) * ntheta_int;
    if (gyro::icntl[Param::DirBC_Interior]) {
        for (int k = 0; k < ntheta_int; ++k) {
            EXPECT_EQ(test_level.row_indices[k], k);
            EXPECT_EQ(test_level.col_indices[k], k);
            EXPECT_EQ(test_level.vals[k], 1.0);

            EXPECT_EQ(test_level.vals[nz - (k + 1)], 1.0);
            EXPECT_EQ(test_level.row_indices[nz - (k + 1)], nz - (k + 1));
            EXPECT_EQ(test_level.col_indices[nz - (k + 1)], nz - (k + 1));
        }
    }
    for (int j = 0; j < nr_int + 1; ++j) {
        for (int i = 0; i < ntheta_int; ++i) {
            if ((j == 0 || j == nr_int)) {
            }
        }
    }
}
