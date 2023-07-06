#include <gtest/gtest.h>
#include "gmgpolar.h"
class test_prolongation : public ::testing::TestWithParam<int>
{
protected:
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
    }
};

TEST_P(test_prolongation, test_bilinear_prolongation)
{
    const int& val_size            = GetParam();
    gyro::icntl[Param::nr_exp]     = (int)(val_size / 3) + 3;
    gyro::icntl[Param::ntheta_exp] = (val_size % 3) + 3;

    if (gyro::icntl[Param::nr_exp] == 3)
        gyro::icntl[Param::fac_ani] = 2;

    gmgpolar test_p;
    test_p.create_grid_polar();
    test_p.check_geom();
    test_p.define_coarse_nodes();

    level& p_level = *(test_p.v_level[0]);
    int ctheta_int = test_p.v_level[1]->ntheta_int;

    p_level.m  = test_p.v_level[0]->nr * test_p.v_level[0]->ntheta;
    p_level.mc = test_p.v_level[1]->nr * test_p.v_level[1]->ntheta;

    std::vector<double> u_test(p_level.mc);
    for (int z = 0; z < p_level.mc; z++) {
        u_test[z] = 1 - z + pow(PI, -z * z);
    }

    std::vector<double> sol = p_level.apply_prolongation_bi(u_test);

    for (int j = 0; j < p_level.nr_int + 1; j++) {
        for (int i = 0; i < p_level.ntheta_int; i++) {
            if (j % 2 == 0 && i % 2 == 0) {
                EXPECT_EQ(u_test[(j / 2) * ctheta_int + (i / 2)], sol[j * p_level.ntheta_int + i])
                    << "coarse node injection is failing";
            }
            else if (i % 2 != 0 && j % 2 == 0) { // theta_i is fine node

                double k_qm1 = p_level.theta[i] - p_level.theta[i - 1];
                double k_q =
                    (i < p_level.ntheta_int - 1) ? p_level.theta[i + 1] - p_level.theta[i] : 2 * PI - p_level.theta[i];

                int i1 = (j / 2) * ctheta_int + (i - 1) / 2;
                int i2 = (i < p_level.ntheta_int - 1) ? i1 + 1 : (j / 2) * ctheta_int; //cyclic coordinates

                double val = (k_q * u_test[i1] + k_qm1 * u_test[i2]);

                EXPECT_NEAR(val, (k_q + k_qm1) * sol[j * p_level.ntheta_int + i], 1e-8)
                    << "Bilinear Prolongation fails for Index (r,theta) : (" + std::to_string(j) + "," +
                           std::to_string(i) + ")";
                ;
            }
            else if (i % 2 == 0) { // r_j fine node, theta_i coarse
                double h_pm1 = p_level.r[j] - p_level.r[j - 1];
                double h_p   = p_level.r[j + 1] - p_level.r[j];

                double v1 = u_test[(j - 1) / 2 * ctheta_int + (i / 2)];
                double v2 = u_test[(j + 1) / 2 * ctheta_int + (i / 2)];

                double val = (h_p * v1 + h_pm1 * v2);

                EXPECT_NEAR(val, (h_p + h_pm1) * sol[j * p_level.ntheta_int + i], 1e-8)
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
                              (h_pm1 * k_qm1 * top_right));

                EXPECT_NEAR(val, (h_p + h_pm1) * (k_q + k_qm1) * sol[j * p_level.ntheta_int + i], 1e-8)
                    << "Bilinear Prolongation fails for Index (r,theta) : (" + std::to_string(j) + "," +
                           std::to_string(i) + ")";
            }
        }
    }
}

TEST_P(test_prolongation, test_injection_prolongation)
{
    const int& val_size            = GetParam();
    gyro::icntl[Param::nr_exp]     = (int)(val_size / 3) + 3;
    gyro::icntl[Param::ntheta_exp] = (val_size % 3) + 3;

    if (gyro::icntl[Param::nr_exp] == 3)
        gyro::icntl[Param::fac_ani] = 2;

    gmgpolar test_p;
    test_p.create_grid_polar();
    test_p.check_geom();
    test_p.define_coarse_nodes();

    level& p_level = *(test_p.v_level[0]);
    int ctheta_int = test_p.v_level[1]->ntheta_int;

    p_level.m  = test_p.v_level[0]->nr * test_p.v_level[0]->ntheta;
    p_level.mc = test_p.v_level[1]->nr * test_p.v_level[1]->ntheta;

    std::vector<double> u_test(p_level.mc);
    for (int z = 0; z < p_level.mc; z++) {
        u_test[z] = z;
    }

    std::vector<double> sol = p_level.apply_prolongation_inj(u_test);

    EXPECT_EQ((int)sol.size(), p_level.m);

    for (int j = 0; j < p_level.nr_int + 1; j++) {
        for (int i = 0; i < p_level.ntheta_int; i++) {
            if (j % 2 == 0 && i % 2 == 0) {
                EXPECT_EQ(u_test[(j / 2) * ctheta_int + (i / 2)], sol[j * p_level.ntheta_int + i])
                    << "The Injection value fails at Index (r,theta): (" + std::to_string(j) + "," + std::to_string(i) +
                           ")";
            }
            else {
                EXPECT_EQ(sol[j * p_level.ntheta_int + i], 0);
            }
        }
    }
}

TEST_P(test_prolongation, test_extrapolation_prolongation)
{
    const int& val_size            = GetParam();
    gyro::icntl[Param::nr_exp]     = (int)(val_size / 3) + 3;
    gyro::icntl[Param::ntheta_exp] = (val_size % 3) + 3;

    if (gyro::icntl[Param::nr_exp] == 3)
        gyro::icntl[Param::fac_ani] = 2;

    gmgpolar test_p;
    test_p.create_grid_polar();
    test_p.check_geom();
    test_p.define_coarse_nodes();

    level& p_level = *(test_p.v_level[0]);
    int ctheta_int = test_p.v_level[1]->ntheta_int;

    p_level.m  = test_p.v_level[0]->nr * test_p.v_level[0]->ntheta;
    p_level.mc = test_p.v_level[1]->nr * test_p.v_level[1]->ntheta;

    std::vector<double> u_test(p_level.mc);
    for (int z = 0; z < p_level.mc; z++) {
        u_test[z] = z;
    }

    std::vector<double> sol = p_level.apply_prolongation_ex(u_test);

    for (int j = 0; j < p_level.nr_int + 1; j++) {
        for (int i = 0; i < p_level.ntheta_int; i++) {
            if (j % 2 == 0 && i % 2 == 0) {
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

                /*in the triangulation we now consider that the fine node is on the hypothenuse of the triangle*/

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

INSTANTIATE_TEST_SUITE_P(Prolongation_size, test_prolongation, ::testing::Values(0, 1, 2, 3, 4, 5, 6, 7, 8));