#include <gtest/gtest.h>
#include "gmgpolar.h"

TEST(Test_prolongation, Test_bilinear_prolongation)
{

    gyro::icntl[Param::DirBC_Interior] = 1; //test should fail if this is 0 no?
    gyro::icntl[Param::optimized]      = 1;

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

    std::vector<double> sol;
    if (gyro::icntl[Param::optimized] == 0) {
        sol = p_level.apply_prolongation_bi0(u_test, p_level.mc, p_level.m, p_level.coarse_nodes_list_r,
                                             p_level.coarse_nodes_list_theta, 0);
    }
    else {
        p_level.define_nz_P();
        p_level.ri_prol = std::vector<int>(p_level.nz_P);
        p_level.ci_prol = std::vector<int>(p_level.nz_P);
        p_level.v_prol  = std::vector<double>(p_level.nz_P);
        p_level.build_prolongation_bi();

        sol = p_level.apply_prolongation_bi(u_test);
    }
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

                EXPECT_NEAR(val, (k_q + k_qm1) * sol[j * p_level.ntheta_int + i], 1e-8);
            }
            else if (i % 2 == 0) { // r_j fine node, theta_i coarse
                double h_pm1 = p_level.r[j] - p_level.r[j - 1];
                double h_p   = p_level.r[j + 1] - p_level.r[j];

                double v1 = (j > 1) ? u_test[(j - 1) / 2 * ctheta_int + (i / 2)] : 0; //Dirichlet BC
                double v2 = (j < p_level.nr_int - 1) ? u_test[(j + 1) / 2 * ctheta_int + (i / 2)] : 0;

                double val = (h_p * v1 + h_pm1 * v2);

                EXPECT_NEAR(val, (h_p + h_pm1) * sol[j * p_level.ntheta_int + i], 1e-8);
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
                int boundary_terms = (j == 1) * ((i == p_level.ntheta_int - 1) + 2 * (i < p_level.ntheta_int - 1)) +
                                     (j == p_level.nr_int - 1) * (3 + (i == p_level.ntheta_int - 1)) +
                                     (j > 1 && j < p_level.nr_int - 1) * (i == p_level.ntheta_int - 1) * 5;

                switch (boundary_terms) {
                default: //(j>1 && j<nr_int-1 && i<ntheta_int-1)
                    bottom_left  = u_test[((j - 1) / 2) * ctheta_int + (i - 1) / 2];
                    top_left     = u_test[((j - 1) / 2) * ctheta_int + (i + 1) / 2];
                    bottom_right = u_test[((j + 1) / 2) * ctheta_int + (i - 1) / 2];
                    top_right    = u_test[((j + 1) / 2) * ctheta_int + (i + 1) / 2];
                    break;
                case 1: // (j==1 && i==ntheta_int-1)
                    bottom_left  = 0;
                    top_left     = 0;
                    bottom_right = u_test[((j + 1) / 2) * ctheta_int + (i - 1) / 2];
                    top_right    = u_test[((j + 1) / 2) * ctheta_int];
                    break;
                case 2: //(j==1 && i<ntheta_int -1)
                    bottom_left  = 0;
                    top_left     = 0;
                    bottom_right = u_test[((j + 1) / 2) * ctheta_int + (i - 1) / 2];
                    top_right    = u_test[((j + 1) / 2) * ctheta_int + (i + 1) / 2];
                    break;
                case 3: //(j==nr_int-1 && i<ntheta_int-1)
                    bottom_left  = u_test[((j - 1) / 2) * ctheta_int + (i - 1) / 2];
                    top_left     = u_test[((j - 1) / 2) * ctheta_int + (i + 1) / 2];
                    bottom_right = 0;
                    top_right    = 0;
                    break;
                case 4: //(j==nr_int-1 && i== ntheta_int-1)
                    bottom_left  = u_test[((j - 1) / 2) * ctheta_int + (i - 1) / 2];
                    top_left     = u_test[((j - 1) / 2) * ctheta_int];
                    bottom_right = 0;
                    top_right    = 0;
                    break;
                case 5: //(j interior, i==ntheta_int-1)
                    bottom_left  = u_test[((j - 1) / 2) * ctheta_int + (i - 1) / 2];
                    top_left     = u_test[((j - 1) / 2) * ctheta_int];
                    bottom_right = u_test[((j + 1) / 2) * ctheta_int + (i - 1) / 2];
                    top_right    = u_test[((j + 1) / 2) * ctheta_int];
                    break;
                }

                double val = ((h_p * k_q * bottom_left) + (h_p * k_qm1 * top_left) + (h_pm1 * k_q * bottom_right) +
                              (h_pm1 * k_qm1 * top_right));

                EXPECT_NEAR(val, (h_p + h_pm1) * (k_q + k_qm1) * sol[j * p_level.ntheta_int + i], 1e-8);
            }
        }
    }
}

TEST(Test_prolongation, Test_injection_prolongation)
{
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

    p_level.ri_prol_inj = std::vector<int>(p_level.nz_P_inj);
    p_level.ci_prol_inj = std::vector<int>(p_level.nz_P_inj);
    p_level.v_prol_inj  = std::vector<double>(p_level.nz_P_inj);
    p_level.build_prolongation_inj();

    std::vector<double> sol = p_level.apply_prolongation_inj(u_test);

    EXPECT_EQ((int)sol.size(), p_level.m);

    for (int j = 0; j < p_level.nr_int + 1; j++) {
        for (int i = 0; i < p_level.ntheta_int; i++) {
            if (j % 2 == 0 && i % 2 == 0) {
                EXPECT_EQ(u_test[(j / 2) * ctheta_int + (i / 2)], sol[j * p_level.ntheta_int + i])
                    << "the injection value is wrong at j=" + std::to_string(j) + " and i=" + std::to_string(i);
            }
            else {
                EXPECT_EQ(sol[j * p_level.ntheta_int + i], 0);
            }
        }
    }
}

TEST(Test_prolongation, Test_extrapolation_prolongation)
{
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

    p_level.define_nz_P();
    p_level.ri_prol_ex = std::vector<int>(p_level.nz_P_ex);
    p_level.ci_prol_ex = std::vector<int>(p_level.nz_P_ex);
    p_level.v_prol_ex  = std::vector<double>(p_level.nz_P_ex);
    p_level.build_prolongation_ex();

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

                EXPECT_NEAR(val, sol[j * p_level.ntheta_int + i], 1e-8);
            }
            else if (i % 2 == 0) { // r_j fine node, theta_i coarse. lower edge of the triangle

                double v1 = (j > 1) ? u_test[(j - 1) / 2 * ctheta_int + (i / 2)] : 0; //zero Dirichlet BC
                double v2 = (j < p_level.nr_int - 1) ? u_test[(j + 1) / 2 * ctheta_int + (i / 2)] : 0;

                double val = 0.5 * (v1 + v2);

                EXPECT_NEAR(val, sol[j * p_level.ntheta_int + i], 1e-8);
            }
            else // both are fine nodes
            {

                /*in the triangulation we now consider that the fine node is on the hypothenuse of the triangle*/

                double top_left;
                double bottom_right;

                ASSERT_NE(p_level.nr_int - 1, 1);
                int boundary_terms = (j == 1) + (j == p_level.nr_int - 1) * (2 + (i == p_level.ntheta_int - 1)) +
                                     (j > 1 && j < p_level.nr_int - 1) * (i == p_level.ntheta_int - 1) * 4;

                switch (boundary_terms) {
                default: //(j>1 && j<nr_int-1 && i<ntheta_int-1)
                    top_left     = u_test[((j - 1) / 2) * ctheta_int + (i + 1) / 2];
                    bottom_right = u_test[((j + 1) / 2) * ctheta_int + (i - 1) / 2];
                    break;
                case 1: // (j==1)
                    top_left     = 0;
                    bottom_right = u_test[((j + 1) / 2) * ctheta_int + (i - 1) / 2];
                    break;
                case 2: //(j==nr_int-1 && i<ntheta_int-1)
                    top_left     = u_test[((j - 1) / 2) * ctheta_int + (i + 1) / 2];
                    bottom_right = 0;
                    break;
                case 3: //(j==nr_int-1 && i== ntheta_int-1)
                    top_left     = u_test[((j - 1) / 2) * ctheta_int];
                    bottom_right = 0;
                    break;
                case 4: //(j interior, i==ntheta_int-1)
                    top_left     = u_test[((j - 1) / 2) * ctheta_int];
                    bottom_right = u_test[((j + 1) / 2) * ctheta_int + (i - 1) / 2];
                    break;
                }

                double val = 0.5 * (top_left + bottom_right);

                EXPECT_NEAR(val, sol[j * p_level.ntheta_int + i], 1e-8);
            }
        }
    }
}

//TODO:  Test restriction operator. Maybe fix interior dirichlet boundary ?
