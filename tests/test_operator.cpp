#include <gtest/gtest.h>
#include "gmgpolar.h"
#include "mockgrid.h"
#include <math.h>

class test_operator : public testing::TestWithParam<std::tuple<int, bool>>
{
protected:
    test_operator()
        : test_level(0)
    {
    }
    void SetUp() override
    {
        gyro::init_params();
        gyro::icntl[Param::mod_pk]         = 0;
        gyro::icntl[Param::beta_coeff]     = 1;
        gyro::icntl[Param::nr_exp]         = 3;
        gyro::icntl[Param::verbose]        = 0;
        gyro::dcntl[Param::R0]             = 1e-5;
        gyro::icntl[Param::periodic]       = 1;
        gyro::icntl[Param::DirBC_Interior] = 1;
        gyro::f_grid_r                     = "";
        gyro::f_grid_theta                 = "";
        gyro::f_sol_in                     = "";
        gyro::f_sol_out                    = "";
        gyro::icntl[Param::fac_ani]        = 3;
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
        test_level.m          = test_level.nr * test_level.ntheta;
        std::cout << "2" << std::endl;
    }

    gmgpolar test_solver;
    level test_level;
};

int set_nonzero(int mod_pk, int nr, int ntheta_int)
{
    int nz             = 0;
    int stencil_points = (mod_pk == 0) ? 5 : 9;
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
    nz += (2 - dirichlets) * (5 + 2 * (mod_pk != 0)) * ntheta_int; //discretization across the interior
    int interior_unlinked_inr =
        std::max(nr - (2 * dirichlets) - (2 - dirichlets), 0); //2 times since for each side we consider linked aswell
    nz += interior_unlinked_inr * ntheta_int * stencil_points;

    return nz;
    //if(!gyro::icntl[Param::DirBC_Interior])
}

std::vector<double> calculate_DFinvTranspose(int mod_pk, double r_j, double theta_i)
{
    std::vector<double> res(4);
    switch (mod_pk) {
    case 0:
        res = {cos(theta_i), -(1 / r_j) * sin(theta_i), sin(theta_i), (1 / r_j) * cos(theta_i)}; //a11, a12, a21, a22
        break;
    case 1:
        //bla
        break;
    default:
        //blabla
        break;
    }
    return res;
}
double calculate_detDF(int mod_pk, double r_j, double theta_i)
{
    double res;
    switch (mod_pk) {
    case 0:
        res = r_j / pow(gyro::dcntl[Param::R], 2); //why divide everything with max_R ???
        break;
    case 1:
        res =
            abs(0.91 * r_j - 0.4 * 1.3 * pow(r_j, 2) * cos(theta_i)); //Cartesian sonnendrucker shafranox makes 0 sense
        break;
    default:
        res = 0;
        break;
    }
    return res;
}

int emulate_get_ptr(int i, int j, level& test_level)
{
    int index = j * test_level.ntheta_int + i;

    if (j == 0) {
        if (gyro::icntl[Param::DirBC_Interior]) {
            return i;
        }
        else {
            int res = (5 + 2 * (gyro::icntl[Param::mod_pk] > 0)) * i;
            return res;
        }
    }
    if (j == 1 && gyro::icntl[Param::DirBC_Interior]) { //nodes linked to dirichlet-bc
        int res = test_level.ntheta_int + (4 + 2 * (gyro::icntl[Param::mod_pk] > 0)) * i;
        return res;
    }
    if (j > gyro::icntl[Param::DirBC_Interior]) {

        int ptr = (5 + 2 * (gyro::icntl[Param::mod_pk] > 0) +
                   (std::min(j, test_level.nr_int - 1) - (gyro::icntl[Param::DirBC_Interior] + 1)) *
                       (5 + 4 * (gyro::icntl[Param::mod_pk] > 0))) *
                  test_level.ntheta_int;

        /*
                    in dirichlet case (1+4+2*(mod_pk>0))*ntheta_int. 1 for Dirichlet nodes, 4+2*(mod_pk>0) for linked to
                    Dirichlets. In this case the interior nodes get counted from j>1 upwards until j<test_level.nr_int-1 
                    
                    in across the origin case we have (5+2*(mod_pk))*ntheta_int nodes which are discretized across the origin
                    and the interior nodes get counted from j>0 upwards until j<test_level.nr_int-1
                */

        if (j < test_level.nr_int - 1) {
            ptr += i * (5 + 4 * (gyro::icntl[Param::mod_pk] > 0));
            return ptr;
        }
        if (j == test_level.nr_int - 1) { //nodes linked to outer dirichlet BC
            ptr += i * (4 + 2 * (gyro::icntl[Param::mod_pk] > 0));
            return ptr;
        }
        if (j == test_level.nr_int) {
            ptr += (4 + 2 * (gyro::icntl[Param::mod_pk] > 0)) * test_level.ntheta_int + i;
            return ptr;
        }
    }
}

TEST_P(test_operator, test_nonzero)
{

    gyro::icntl[Param::DirBC_Interior] = std::get<1>(GetParam());
    gyro::icntl[Param::mod_pk]         = std::get<0>(GetParam());
    int& ntheta_int                    = test_level.ntheta_int;

    test_level.define_nz();
    EXPECT_EQ(test_level.nz, set_nonzero(gyro::icntl[Param::mod_pk], test_level.nr, test_level.ntheta_int))
        << "mod_pk: " + std::to_string(gyro::icntl[Param::mod_pk]) +
               " Dir_BC: " + std::to_string(gyro::icntl[Param::DirBC_Interior]);
}

TEST_P(test_operator, build_beta)
{
    test_level.store_theta_n_co(); //below one r[j], whole line of theta
    test_level.betaVec.resize(test_level.nr * test_level.ntheta_int);
    test_level.build_betaVec();
    for (int j = 0; j < test_level.nr_int; ++j) {
        double hjmin1 = j > 0 ? (test_level.r[j] - test_level.r[j - 1]) : 2 * test_level.r[0];
        double hj     = test_level.r[j + 1] - test_level.r[j];
        for (int i = 0; i < test_level.ntheta_int; ++i) {
            double det    = calculate_detDF(gyro::icntl[Param::mod_pk], test_level.r[j], test_level.theta[i]);
            double kimin1 = i > 0 ? (test_level.theta[i] - test_level.theta[i - 1])
                                  : 2 * PI - test_level.theta[test_level.ntheta_int - 1];
            double ki = i < test_level.ntheta_int - 1 ? (test_level.theta[i + 1] - test_level.theta[i])
                                                      : 2 * PI - test_level.theta[i];
            int ind = j * test_level.ntheta_int + i;

            double res = 0.25 * (hj + hjmin1) * (ki + kimin1) * gyro::coeff_beta(test_level.r[j], 0) * fabs(det);
            EXPECT_NEAR(res, test_level.betaVec[ind], 1e-6) << std::to_string(j) + " " + std::to_string(i);
        }
    }
}

TEST_P(test_operator, get_ptr)
{
    gyro::icntl[Param::DirBC_Interior] = std::get<1>(GetParam());
    gyro::icntl[Param::mod_pk]         = std::get<0>(GetParam());
    for (int i = 0; i < test_level.ntheta_int; ++i) {
        for (int j = 0; j < test_level.nr; ++j) {
            int index = j * test_level.ntheta_int + i;

            EXPECT_EQ(test_level.get_ptr(i, j), emulate_get_ptr(i, j, test_level));
        }
    }
}
/*
TEST_F(test_operator, build_A)
{
    int& ntheta_int = test_level.ntheta_int;
    int& nr_int     = test_level.nr_int;
    int nz          = set_nonzero(gyro::icntl[Param::mod_pk], test_level.nr, ntheta_int);
    test_level.vals.resize(nz);
    test_level.row_indices.resize(nz);
    test_level.col_indices.resize(nz);
    test_level.store_theta_n_co();
    test_level.betaVec.resize(test_level.nr * test_level.ntheta_int);
    test_level.build_betaVec();
    test_level.build_A();

    int size = test_level.vals.size();
    std::cout << test_level.nr << std::endl;
    std::cout << test_level.ntheta << std::endl;

    EXPECT_EQ(size, nz);
    EXPECT_EQ(test_level.row_indices.size(), nz);
    EXPECT_EQ(test_level.col_indices.size(), nz);
    int gridp = (nr_int + 1) * ntheta_int;
    if (gyro::icntl[Param::DirBC_Interior]) {
        for (int k = 0; k < ntheta_int; ++k) {
            std::cout << "test?" << std::endl;
            EXPECT_EQ(test_level.row_indices[k], k);
            EXPECT_EQ(test_level.col_indices[k], k);
            EXPECT_EQ(test_level.vals[k], 1.0);
            std::cout << "bis hier kein seg" << std::endl;
            EXPECT_EQ(test_level.vals[nz - (k + 1)], 1.0);
            EXPECT_EQ(test_level.row_indices[nz - (k + 1)], gridp - (k + 1));
            EXPECT_EQ(test_level.col_indices[nz - (k + 1)], gridp - (k + 1));

            //Linked Nodes 5p stencil
            int idx = ntheta_int + k;
            test_level.row_indices[]
        }
    }
}
*/
INSTANTIATE_TEST_SUITE_P(Param, test_operator,
                         ::testing::Values(std::make_tuple(0, false), std::make_tuple(1, false),
                                           std::make_tuple(2, false), std::make_tuple(0, true),
                                           std::make_tuple(1, true), std::make_tuple(2, true)));