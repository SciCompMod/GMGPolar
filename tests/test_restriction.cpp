#include <gtest/gtest.h>
#include "gmgpolar.h"
class test_restriction : public ::testing::TestWithParam<int>
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

TEST_P(test_restriction, test_bilinear_restriction)
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
    int cr_int     = test_p.v_level[1]->nr_int;

    p_level.m  = test_p.v_level[0]->nr * test_p.v_level[0]->ntheta;
    p_level.mc = test_p.v_level[1]->nr * test_p.v_level[1]->ntheta;

    std::vector<double> u_test(p_level.m);
    for (int z = 0; z < p_level.m; z++) {
        u_test[z] = 1 - z;
    }

    std::vector<double> sol = p_level.apply_restriction_bi(u_test);

    /*
     * The restriction operator accumulates the values in the direct vicinity and stores them in the corresponding node.
     * We hence calculate 8 values for every coarse node. These 8 values are the fine nodes in the vicinity that are to be accumulated
     * in the coarse node
     *
     * we treat values in r as the x-axis. values in theta as the y-axis. Hence the vector 'adjacent' stores the values in the order:

     * **bottom_left, left, top_left, bottom, top, bottom_right, right, top_right;**
     */

    for (int j = 0; j < cr_int + 1; j++) {
        for (int i = 0; i < ctheta_int; i++) {
            std::vector<double> adjacent(8, -1.0);
            if (j == 0) { //interior circle- nodes to the left disappear
                adjacent[0] = 0.0;
                adjacent[1] = 0.0;
                adjacent[2] = 0.0;
            }
            if (j == cr_int) { //exterior circle- nodes to the right disappear
                adjacent[5] = 0.0;
                adjacent[6] = 0.0;
                adjacent[7] = 0.0;
            }

            int k = 0;
            std::vector<std::tuple<double, double>> h_p(8, {-1, -1}); // (h_p, h_pm1)
            std::vector<std::tuple<double, double>> k_q(8, {-1, -1}); // (k_q, k_qm1)

            for (int z = -1; z < 2; z++) {
                for (int y = -1; y < 2; y++) {
                    if (z != 0 || y != 0) {
                        if (adjacent[k] != 0.0) {

                            adjacent[k] =
                                u_test[(2 * j + z) * p_level.ntheta_int +
                                       ((i != 0 || y > -1) ? (2 * i + y) : p_level.ntheta_int - 1)]; //adjacent value

                            h_p[k] = {p_level.r[2 * j + z + 1] - p_level.r[2 * j + z],
                                      p_level.r[2 * j + z] -
                                          p_level.r[2 * j + z - 1]}; //distance in r coordinate for the adjacent node

                            if (i > 0) {
                                double indx =
                                    (2 * i + y < p_level.ntheta_int - 1) ? p_level.theta[2 * i + y + 1] : 2 * PI;

                                k_q[k] = {indx - p_level.theta[2 * i + y],
                                          p_level.theta[2 * i + y] - p_level.theta[2 * i + y - 1]};
                            }
                            else {
                                switch (y) {
                                case 1:
                                    k_q[k] = {p_level.theta[2] - p_level.theta[1], p_level.theta[1] - p_level.theta[0]};
                                case 0:
                                    k_q[k] = {p_level.theta[1] - p_level.theta[0],
                                              2 * PI - p_level.theta[p_level.ntheta_int - 1]};
                                default:
                                    k_q[k] = {2 * PI - p_level.theta[p_level.ntheta_int - 1],
                                              p_level.theta[p_level.ntheta_int - 1] -
                                                  p_level.theta[p_level.ntheta_int - 2]};
                                }
                            }
                        }
                        k += 1;
                    }
                }
            }
            //values given by the operator
            std::vector<double> vals{
                std::get<1>(h_p[0]) * std::get<1>(k_q[0]) /
                    ((std::get<0>(h_p[0]) + std::get<1>(h_p[0])) * (std::get<0>(k_q[0]) + std::get<1>(k_q[0]))),
                std::get<1>(h_p[1]) / (std::get<0>(h_p[1]) + std::get<1>(h_p[1])),
                std::get<1>(h_p[2]) * std::get<0>(k_q[2]) /
                    ((std::get<0>(h_p[2]) + std::get<1>(h_p[2])) * (std::get<0>(k_q[2]) + std::get<1>(k_q[2]))),
                std::get<1>(k_q[3]) / (std::get<0>(k_q[3]) + std::get<1>(k_q[3])),
                std::get<0>(k_q[4]) / (std::get<0>(k_q[4]) + std::get<1>(k_q[4])),
                std::get<0>(h_p[5]) * std::get<1>(k_q[5]) /
                    ((std::get<0>(h_p[5]) + std::get<1>(h_p[5])) * (std::get<0>(k_q[5]) + std::get<1>(k_q[5]))),
                std::get<0>(h_p[6]) / (std::get<0>(h_p[6]) + std::get<1>(h_p[6])),
                std::get<0>(h_p[7]) * std::get<0>(k_q[7]) /
                    ((std::get<0>(h_p[7]) + std::get<1>(h_p[7])) * (std::get<0>(k_q[7]) + std::get<1>(k_q[7])))};

            double finval = u_test[(2 * j) * p_level.ntheta_int + (2 * i)];
            for (int z = 0; z < 8; z++) {
                finval += vals[z] * adjacent[z];
            }

            EXPECT_NEAR(finval, sol[j * ctheta_int + i], 1e-6)
                << "The test fails at Index for (r,theta): (" + std::to_string(j) + "," + std::to_string(i) + ")";
        }
    }
}

TEST_P(test_restriction, test_injection_restriction)
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
    int cr_int     = test_p.v_level[1]->nr_int;

    p_level.m  = test_p.v_level[0]->nr * test_p.v_level[0]->ntheta;
    p_level.mc = test_p.v_level[1]->nr * test_p.v_level[1]->ntheta;

    std::vector<double> u_test(p_level.m);
    for (int z = 0; z < p_level.m; z++) {
        u_test[z] = z;
    }

    std::vector<double> sol = p_level.apply_restriction_inj(u_test);

    EXPECT_EQ((int)sol.size(), p_level.mc);

    for (int j = 0; j < cr_int + 1; j++) {
        for (int i = 0; i < ctheta_int; i++) {
            EXPECT_EQ(sol[j * ctheta_int + i], u_test[(2 * j) * p_level.ntheta_int + (2 * i)])
                << "The injection restriction fails at Index (r,theta): (" + std::to_string(j) + "," +
                       std::to_string(i) + ")";
        }
    }
}

TEST_P(test_restriction, test_extrapolation_restriction)
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
    int cr_int     = test_p.v_level[1]->nr_int;

    p_level.m  = test_p.v_level[0]->nr * test_p.v_level[0]->ntheta;
    p_level.mc = test_p.v_level[1]->nr * test_p.v_level[1]->ntheta;

    std::vector<double> u_test(p_level.m);
    for (int z = 0; z < p_level.m; z++) {
        u_test[z] = 1 - z;
    }

    std::vector<double> sol = p_level.apply_restriction_ex(u_test);

    /* based on the triangulation, at most 6 adjacent fine nodes are considered:
     * **left,top_left,bottom,top,bottom_right,right**
     */
    for (int j = 0; j < cr_int + 1; j++) {
        for (int i = 0; i < ctheta_int; i++) {
            std::vector<double> adjacent(6, -1.0);
            if (j == 0) { //interior circle
                adjacent[0] = 0.0;
                adjacent[1] = 0.0;
            }
            if (j == cr_int) { //exterior circle
                adjacent[4] = 0.0;
                adjacent[5] = 0.0;
            }

            int k = 0;

            for (int z = -1; z < 2; z++) {
                for (int y = -1; y < 2; y++) {
                    if ((z != 0 || y != 0) && (z * y != 1)) {
                        if (adjacent[k] != 0.0) {

                            adjacent[k] = u_test[(2 * j + z) * p_level.ntheta_int +
                                                 ((i != 0 || y > -1) ? (2 * i + y) : p_level.ntheta_int - 1)];
                        }
                        k += 1;
                    }
                }
            }

            double finval = u_test[(2 * j) * p_level.ntheta_int + (2 * i)];
            for (int z = 0; z < 8; z++) {
                finval += 0.5 * adjacent[z];
            }

            EXPECT_NEAR(finval, sol[j * ctheta_int + i], 1e-6)
                << "The test fails at Index for (r,theta): (" + std::to_string(j) + "," + std::to_string(i) + ")";
        }
    }
}

INSTANTIATE_TEST_SUITE_P(Restriction_size, test_restriction, ::testing::Values(0, 1, 2, 3, 4, 5, 6, 7, 8));