#include <gtest/gtest.h>
#include <tuple>
#include "gmgpolar.h"
#include "mockgrid.h"
class test_prolongation
    : public ::testing::TestWithParam<
          std::tuple<int, bool>> //tuple includes data on grid size and Dirichlet boundary conditions
{
protected:
    void SetUp() override
    {
        //initialize default parameters.
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
    }
};

/*!
 *  \brief Test the bilinear prolongation operator used in the multigrid cycle coarse-grid correction. 
 *  
 *  The Test creates an arbitrary grid-function on the coarser level and prolongates it. 
 *  On the fine level we iterate over all nodes and test the result based on whether the node is fine in theta, r or in both.
 *
 *  Parametrized tests are used to test for different grid sizes and with or without Dirichlet boundary conditions.
 */

TEST_P(test_prolongation, test_bilinear_prolongation)
{
    //we vary the grid size to guarantee that the problem works for all sizes
    gyro::icntl[Param::nr_exp]         = (int)(std::get<0>(GetParam()) / 3) + 3;
    gyro::icntl[Param::ntheta_exp]     = (std::get<0>(GetParam()) % 3) + 3;
    gyro::icntl[Param::DirBC_Interior] = std::get<1>(GetParam());
    if (gyro::icntl[Param::nr_exp] == 3)
        gyro::icntl[Param::fac_ani] = 2; //anisotropy should not exceed grid size

    gmgpolar test_p;
    create_grid(test_p);

    level& p_level = *(test_p.v_level[0]);
    int ctheta_int = test_p.v_level[1]->ntheta_int; //number of coarse nodes in theta direction

    p_level.m  = test_p.v_level[0]->nr * test_p.v_level[0]->ntheta; //fine grid size
    p_level.mc = test_p.v_level[1]->nr * test_p.v_level[1]->ntheta; //coarser grid size

    std::vector<double> u_test(p_level.mc);
    for (int z = 0; z < p_level.mc; z++) {
        u_test[z] =
            1 - z +
            pow(PI,
                -z * z); //constructing arbitrary grid-function on coarse level to test our prolongation operator with.
    }

    std::vector<double> sol =
        p_level.apply_prolongation_bi(u_test); //apply prolongation operator on arbitrary grid-function

    for (int j = 0; j < p_level.nr_int + 1; j++) {
        for (int i = 0; i < p_level.ntheta_int; i++) {
            if (j % 2 == 0 && i % 2 == 0) { //testing coarse node injection
                EXPECT_EQ(u_test[(j / 2) * ctheta_int + (i / 2)], sol[j * p_level.ntheta_int + i])
                    << "coarse node injection is failing";
            }
            else if (i % 2 != 0 && j % 2 == 0) {
            // as numbering in angle (theta_i) is odd, we have a fine node with
            // coarse neighbors at (r_j, theta_i - k_{i-1}) and (r_j, theta_i + k_i)

                double k_qm1 = p_level.theta[i] - p_level.theta[i - 1]; //calculate k_{q-1}
                double k_q   = (i < p_level.ntheta_int - 1) ? p_level.theta[i + 1] - p_level.theta[i]
                                                          : 2 * PI - p_level.theta[i]; //k_q

                int i1 = (j / 2) * ctheta_int + (i - 1) / 2; //bottom coarse node in theta
                int i2 = (i < p_level.ntheta_int - 1) ? i1 + 1 : (j / 2) * ctheta_int; //top coarse node in theta

                double val = (k_q * u_test[i1] + k_qm1 * u_test[i2]);

                EXPECT_NEAR(val, (k_q + k_qm1) * sol[j * p_level.ntheta_int + i], 1e-8) //compare values as in the paper
                    << "Bilinear Prolongation fails for Index (r,theta) : (" + std::to_string(j) + "," +
                           std::to_string(i) + ")";
                ;
            }
            else if (i % 2 == 0 && j % 2 != 0) {
            // as numbering in radius (r_j) is odd, we have a fine node with
            // coarse neighbors at (r_j - h_{j-1}, theta_i) and (r_j+h_j, theta_i )
                double h_pm1 = p_level.r[j] - p_level.r[j - 1]; //h_{p-1}
                double h_p   = p_level.r[j + 1] - p_level.r[j]; //h_p

                double v1 = u_test[(j - 1) / 2 * ctheta_int + (i / 2)]; // left coarse node in r
                double v2 = u_test[(j + 1) / 2 * ctheta_int + (i / 2)]; // right coarse node in r

                double val = (h_p * v1 + h_pm1 * v2);

                EXPECT_NEAR(val, (h_p + h_pm1) * sol[j * p_level.ntheta_int + i], 1e-8) //compare values
                    << "Bilinear Prolongation fails for Index (r,theta) : (" + std::to_string(j) + "," +
                           std::to_string(i) + ")";
                ;
            }
            else // both are fine nodes
            {
                double h_pm1 = p_level.r[j] - p_level.r[j - 1];
                double h_p   = p_level.r[j + 1] - p_level.r[j];

                double k_qm1 = p_level.theta[i] - p_level.theta[i - 1];
                double k_q =
                    (i < p_level.ntheta_int - 1) ? p_level.theta[i + 1] - p_level.theta[i] : 2 * PI - p_level.theta[i];

                /*the stencil-corners corresponding to the values (r_j,theta_i)*/
                /*bottom corresponds to lower indices in theta direction. left to lower indices in radius direction*/

                double bottom_left;
                double top_left;
                double bottom_right;
                double top_right;
                ASSERT_NE(p_level.nr_int - 1, 1);

                if (i < p_level.ntheta - 1) {
                    bottom_left  = u_test[((j - 1) / 2) * ctheta_int + (i - 1) / 2];
                    top_left     = u_test[((j - 1) / 2) * ctheta_int + (i + 1) / 2];
                    bottom_right = u_test[((j + 1) / 2) * ctheta_int + (i - 1) / 2];
                    top_right    = u_test[((j + 1) / 2) * ctheta_int + (i + 1) / 2];
                }
                else {
                    bottom_left  = u_test[((j - 1) / 2) * ctheta_int + (i - 1) / 2];
                    top_left     = u_test[((j - 1) / 2) * ctheta_int];
                    bottom_right = u_test[((j + 1) / 2) * ctheta_int + (i - 1) / 2];
                    top_right    = u_test[((j + 1) / 2) * ctheta_int];
                }

                double val = ((h_p * k_q * bottom_left) + (h_p * k_qm1 * top_left) + (h_pm1 * k_q * bottom_right) +
                              (h_pm1 * k_qm1 * top_right)); //calculate value as in the paper

                EXPECT_NEAR(val, (h_p + h_pm1) * (k_q + k_qm1) * sol[j * p_level.ntheta_int + i], 1e-8)
                    << "Bilinear Prolongation fails for Index (r,theta) : (" + std::to_string(j) + "," +
                           std::to_string(i) + ")";
            }
        }
    }
}

/*!
 *  \brief Test the injection prolongation operator used in the implicit extrapolation step of the multigrid-cycle. 
 *  
 *  The Test creates an arbitrary grid-function on the coarser level and injects it. 
 *  On the fine level we iterate over all nodes and test the result based on whether the node is fine or not.
 *
 *  Parametrized tests are used to test for different grid sizes and with or without Dirichlet boundary conditions.
 */

TEST_P(test_prolongation, test_injection_prolongation)
{
    gyro::icntl[Param::nr_exp]         = (int)(std::get<0>(GetParam()) / 3) + 3;
    gyro::icntl[Param::ntheta_exp]     = (std::get<0>(GetParam()) % 3) + 3;
    gyro::icntl[Param::DirBC_Interior] = std::get<1>(GetParam());

    if (gyro::icntl[Param::nr_exp] == 3)
        gyro::icntl[Param::fac_ani] = 2;

    gmgpolar test_p;

    create_grid(test_p);

    level& p_level = *(test_p.v_level[0]);
    int ctheta_int = test_p.v_level[1]->ntheta_int;

    p_level.m  = test_p.v_level[0]->nr * test_p.v_level[0]->ntheta;
    p_level.mc = test_p.v_level[1]->nr * test_p.v_level[1]->ntheta;

    std::vector<double> u_test(p_level.mc);
    for (int z = 0; z < p_level.mc; z++) {
        u_test[z] = z; //arbitrary grid-function that we prolongate defined on the coarse nodes
    }

    std::vector<double> sol = p_level.apply_prolongation_inj(u_test);

    EXPECT_EQ((int)sol.size(), p_level.m);
    //testing the embedding. since we embed the coarse gird function, the values of the fine node should be zero
    for (int j = 0; j < p_level.nr_int + 1; j++) {
        for (int i = 0; i < p_level.ntheta_int; i++) {
            if (j % 2 == 0 && i % 2 == 0) {
                EXPECT_EQ(u_test[(j / 2) * ctheta_int + (i / 2)], sol[j * p_level.ntheta_int + i])
                    << "The Injection value fails at Index (r,theta): (" + std::to_string(j) + "," + std::to_string(i) +
                           ")";
            }
            else {
                EXPECT_EQ(sol[j * p_level.ntheta_int + i], 0); //fine nodes set to zero.
            }
        }
    }
}
/*!
 *  \brief Test the extrapolation prolongation operator used in the implicit extrapolation step of the multigrid-cycle. 
 *  
 *  The Test creates an arbitrary grid-function on the coarser level and prolongates it. 
 *  On the fine level we iterate over all nodes and test the result based on whether the node is fine in theta, r or in both.
 *  The value is determined by a Triangulation of the grid on the coarser level.
 *  
 *  Parametrized tests are used to test for different grid sizes and with or without Dirichlet boundary conditions.
 */

TEST_P(test_prolongation, test_extrapolation_prolongation)
{
    gyro::icntl[Param::nr_exp]         = (int)(std::get<0>(GetParam()) / 3) + 3;
    gyro::icntl[Param::ntheta_exp]     = (std::get<0>(GetParam()) % 3) + 3;
    gyro::icntl[Param::DirBC_Interior] = std::get<1>(GetParam());

    if (gyro::icntl[Param::nr_exp] == 3)
        gyro::icntl[Param::fac_ani] = 2;

    gmgpolar test_p;

    create_grid(test_p);

    level& p_level = *(test_p.v_level[0]);
    int ctheta_int = test_p.v_level[1]->ntheta_int;

    p_level.m  = test_p.v_level[0]->nr * test_p.v_level[0]->ntheta;
    p_level.mc = test_p.v_level[1]->nr * test_p.v_level[1]->ntheta;

    std::vector<double> u_test(p_level.mc);
    for (int z = 0; z < p_level.mc; z++) {
        u_test[z] = z; //arbitrary grid-function to test the prolongation
    }
    std::vector<double> sol = p_level.apply_prolongation_ex(u_test);

    for (int j = 0; j < p_level.nr_int + 1; j++) {
        for (int i = 0; i < p_level.ntheta_int; i++) {
            if (j % 2 == 0 && i % 2 == 0) { //coarse node injection
                EXPECT_EQ(u_test[(j / 2) * ctheta_int + (i / 2)], sol[j * p_level.ntheta_int + i])
                    << "coarse node injection is failing";
            }
            else if (i % 2 != 0 && j % 2 == 0) { // theta_i is fine node. upper edge of the triangle

                int i1 = (j / 2) * ctheta_int + (i - 1) / 2;
                int i2 = (i < p_level.ntheta_int - 1) ? i1 + 1 : (j / 2) * ctheta_int; //cyclic coordinates

                double val = 0.5 * (u_test[i1] + u_test[i2]);

                EXPECT_NEAR(val, sol[j * p_level.ntheta_int + i], 1e-8)
                    << "Extrapolated Prolongation fails for Index (r,theta): (" + std::to_string(j) + "," +
                           std::to_string(i) + ")";
            }
            else if (i % 2 == 0) { // r_j fine node, theta_i coarse. lower edge of the triangle

                double v1 = u_test[(j - 1) / 2 * ctheta_int + (i / 2)];
                double v2 = u_test[(j + 1) / 2 * ctheta_int + (i / 2)];

                double val = 0.5 * (v1 + v2);

                EXPECT_NEAR(val, sol[j * p_level.ntheta_int + i], 1e-8)
                    << "Extrapolated Prolongation fails for Index (r,theta): (" + std::to_string(j) + "," +
                           std::to_string(i) + ")";
            }
            else // both are fine nodes
            {

                /*in the triangulation we now consider that the fine node is on the hypothenuse of the triangle.
                 The corresponding coarse nodes are located on the triangles nodes that span the hypothenuse */

                double top_left;
                double bottom_right;

                if (i < p_level.ntheta_int - 1) {
                    top_left     = u_test[((j - 1) / 2) * ctheta_int + (i + 1) / 2];
                    bottom_right = u_test[((j + 1) / 2) * ctheta_int + (i - 1) / 2];
                }
                else {
                    top_left     = u_test[((j - 1) / 2) * ctheta_int];
                    bottom_right = u_test[((j + 1) / 2) * ctheta_int + (i - 1) / 2];
                }

                double val = 0.5 * (top_left + bottom_right);

                EXPECT_NEAR(val, sol[j * p_level.ntheta_int + i], 1e-8)
                    << "Extrapolated Prolongation fails for Index (r,theta): (" + std::to_string(j) + "," +
                           std::to_string(i) + ")";
                ;
            }
        }
    }
}

INSTANTIATE_TEST_SUITE_P(Prolongation_size, test_prolongation,
                         ::testing::Values(std::make_tuple(0, false), std::make_tuple(1, false),
                                           std::make_tuple(2, false), std::make_tuple(3, false),
                                           std::make_tuple(4, false), std::make_tuple(5, false),
                                           std::make_tuple(6, false), std::make_tuple(7, false),
                                           std::make_tuple(8, false), std::make_tuple(0, true),
                                           std::make_tuple(1, true), std::make_tuple(2, true), std::make_tuple(3, true),
                                           std::make_tuple(4, true), std::make_tuple(5, true), std::make_tuple(6, true),
                                           std::make_tuple(7, true), std::make_tuple(8, true)));