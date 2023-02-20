/* 
* Copyright (C) 2019-2023 The GMGPolar Development Team
*
* Authors: Philippe Leleux, Christina Schwarz, Martin J. Kühn, Carola Kruse, Ulrich Rüde
*
* Contact: 
*    Carola Kruse <Carola.Kruse@CERFACS.fr>
*    Martin J. Kuehn <Martin.Kuehn@DLR.de>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/

/*!
 * \file debug.cpp
 * \brief Comparison of function results with the matlab code (A, f, P, smoothers)
 * \author M. Kuehn, C. Kruse, P. Leleux, C. Schwarz
 * \version 0.0
 */
#include "gmgpolar.h"

inline void unit_test(std::stringstream& fname, int m, std::vector<double>& uC, std::string vect_name, int test)
{
    std::fstream f;
    std::vector<double> data;
    if (gyro::icntl[Param::debug] == 1) {
        f.open(fname.str().c_str(), std::ios::in);
        if (!f) {
            std::cout << "No such file: " << fname.str() << "\n";
        }
        else {
            while (1) {
                double p;
                f >> p;
                if (f.eof())
                    break;
                data.push_back(p);
            }
        }
        f.close();
        fname.clear();
        fname.str(std::string());
        if (data.size() == 0) {
            std::cout << "Warning: array to test is empty.";
        }
        else {
            if (test && data.size() > 0) {
                for (int i = 0; i < m; i++) {
                    if (fabs(uC[i] - data[i]) > 1e-12)
                        std::cout << "i: " << i << ", uC[i]: " << uC[i] << ", " << vect_name << ": " << data[i] << "\n";
                    assert(fabs(uC[i] - data[i]) < 1e-12);
                }
            }
            if (!test) {
                uC = data;
            }
        }
    }
    else {
        if (test && uC.size() == 0) {
            std::cout << "Warning: array to test is empty.\n";
            return;
        }
        if (!test)
            std::srand(static_cast<unsigned int>(std::time(nullptr)));
        f.open(fname.str().c_str(), std::ios::out);
        fname.clear();
        fname.str(std::string());
        if (!f) {
            std::cout << "No such file\n";
        }
        else {
            for (int j = 0; j < m; j++) {
                if (!test) {
                    double random_nb = (double)std::rand() / (double)RAND_MAX;
                    uC.push_back(random_nb);
                }
                f << std::setprecision(17) << uC[j] << "\n";
            }
        }
        f.close();
        fname.clear();
        fname.str(std::string());
    }
}

inline void unit_test(std::stringstream& fname, int m, std::vector<int>& uC, std::string vect_name, int test)
{
    std::fstream f;
    std::vector<int> data;
    if (gyro::icntl[Param::debug] == 1) {
        f.open(fname.str().c_str(), std::ios::in);
        if (!f) {
            std::cout << "No such file: " << fname.str() << "\n";
        }
        else {
            while (1) {
                int p;
                f >> p;
                if (f.eof())
                    break;
                data.push_back(p);
            }
        }
        f.close();
        fname.clear();
        fname.str(std::string());
        if (data.size() == 0) {
            std::cout << "Warning: array to test is empty.";
        }
        else {
            for (int i = 0; i < m; i++) {
                if (uC[i] != data[i])
                    std::cout << "i: " << i << ", uC[i]: " << uC[i] << ", " << vect_name << ": " << data[i] << "\n";
                assert(uC[i] == data[i]);
            }
        }
    }
    else {
        if (uC.size() == 0) { // Sometimes Asc is empty
            std::cout << "Warning: array to test is empty.\n";
            return;
        }
        f.open(fname.str().c_str(), std::ios::out);
        fname.clear();
        fname.str(std::string());
        if (!f) {
            std::cout << "No such file\n";
        }
        else {
            for (int j = 0; j < m; j++) {
                f << std::setprecision(17) << uC[j] << "\n";
            }
        }
        f.close();
        fname.clear();
        fname.str(std::string());
    }
}

/*!
 *  \brief Check the implementation
 *
 * Compare the results of the current implementation with
 * structures extracted from the MATLAB implementation.
 *
 */
void gmgpolar::debug()
{
    std::stringstream ss, folder, fname;
    std::fstream f, fr, fc;
    int optimized      = gyro::icntl[Param::optimized];
    int DirBC_Interior = gyro::icntl[Param::DirBC_Interior];
    int smoother       = gyro::icntl[Param::smoother];
    int extrapolation  = gyro::icntl[Param::extrapolation];
    double R           = gyro::dcntl[Param::R];
    int prob           = gyro::icntl[Param::prob];
    int mod_pk         = gyro::icntl[Param::mod_pk];
    int alpha          = gyro::icntl[Param::alpha_coeff];
    int beta           = gyro::icntl[Param::beta_coeff];

    ss << "debug_data/";
    if (gyro::icntl[Param::debug] != 1) {
        const int dir_err0 = mkdir(ss.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        // if (-1 == dir_err0) {
        //     printf("Error creating directory!\n");
        // }
    }

    ss << "OPTI-" << optimized << "_BC-" << DirBC_Interior << "_SMOOTH-" << smoother << "_EXTR-" << extrapolation
       << "_R-" << R << "_PROB-" << prob << "_GEOM-" << mod_pk << "_A-" << alpha << "_B-" << beta << "/";

    if (gyro::icntl[Param::debug] != 1) {
        const int dir_err1 = mkdir(ss.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        // if (-1 == dir_err1) {
        //     printf("Error creating directory!\n");
        // }
    }

    ///////////////////////////////////////////////////////////
    // Parameters: params.txt (compute_rho, cycle, delta_e, fac_ani, kappa_eps, levels, maxiter, mod_pk, nr, ntheta, periodic, plotit, R0, r0_DB, rDBC(1), rDBC(2), Rmax, solveit, v1, v2)
    ///////////////////////////////////////////////////////////
    std::cout << "Checking parameters...\n";
    fname << ss.str() << "params.txt";
    std::vector<double> parray;
    if (gyro::icntl[Param::debug] == 1) {
        f.open(fname.str().c_str(), std::ios::in);
        if (!f) {
            std::cout << "No such file: " << fname.str() << "\n";
        }
        else {
            while (1) {
                double p;
                f >> p;
                parray.push_back(p);
                // std::cout << p << "\n";
                if (f.eof())
                    break;
            }

            assert(gyro::icntl[Param::optimized] == parray[0]);
            // assert(gyro::icntl[Param::verbose] == parray[1]);
            // assert(gyro::icntl[Param::openmp] == parray[2]);
            // assert(gyro::icntl[Param::matrix_free] == parray[3]);
            // assert(gyro::icntl[Param::debug] == parray[1]);
            assert(gyro::icntl[Param::nr_exp] == (int)parray[1]);
            assert(gyro::icntl[Param::ntheta_exp] == (int)parray[2]);
            assert(gyro::icntl[Param::fac_ani] == (int)parray[3]);
            assert(gyro::icntl[Param::theta_aniso] == (int)parray[4]);
            assert(gyro::icntl[Param::v1] == (int)parray[5]);
            assert(gyro::icntl[Param::v2] == (int)parray[6]);
            assert(gyro::icntl[Param::cycle] == (int)parray[7]);
            assert(gyro::icntl[Param::mod_pk] == (int)parray[8]);
            assert(gyro::icntl[Param::compute_rho] == (int)parray[9]);
            assert(gyro::icntl[Param::level] == (int)parray[10]);
            assert(gyro::icntl[Param::plotit] == (int)parray[11]);
            assert(gyro::icntl[Param::solveit] == (int)parray[12]);
            assert(gyro::icntl[Param::maxiter] == (int)parray[13]);
            assert(gyro::icntl[Param::periodic] == (int)parray[14]);
            assert(gyro::icntl[Param::origin_NOT_coarse] == (int)parray[15]);
            assert(gyro::icntl[Param::smoother] == (int)parray[16]);
            assert(gyro::icntl[Param::discr] == (int)parray[17]);
            assert(gyro::icntl[Param::extrapolation] == (int)parray[18]);
            assert(gyro::icntl[Param::DirBC_Interior] == (int)parray[19]);
            assert(gyro::icntl[Param::paraview] == (int)parray[20]);
            assert(gyro::icntl[Param::divideBy2] == (int)parray[21]);
            assert(gyro::icntl[Param::prob] == (int)parray[22]);
            assert(gyro::icntl[Param::alpha_coeff] == (int)parray[23]);
            assert(gyro::icntl[Param::beta_coeff] == (int)parray[24]);
            assert(gyro::icntl[Param::res_norm] == (int)parray[25]);
            assert(fabs(gyro::dcntl[Param::r0_DB] - parray[26]) < 1e-12);
            assert(fabs(gyro::dcntl[Param::R0] - parray[27]) < 1e-12);
            assert(fabs(gyro::dcntl[Param::R] - parray[28]) < 1e-12);
            assert(fabs(gyro::dcntl[Param::THETA0] - parray[29]) < 1e-12);
            assert(fabs(gyro::dcntl[Param::THETA] - parray[30]) < 1e-12);
            assert(fabs(gyro::dcntl[Param::kappa_eps] - parray[31]) < 1e-12);
            assert(fabs(gyro::dcntl[Param::delta_e] - parray[32]) < 1e-12);
            assert(fabs(gyro::dcntl[Param::tol_bound_check] - parray[33]) < 1e-12);
            assert(fabs(gyro::dcntl[Param::rel_red_conv] - parray[34]) < 1e-12);
        }
        f.close();
        fname.clear(); //clear any bits set
        fname.str(std::string());
    }
    else {
        f.open(fname.str().c_str(), std::ios::out);
        if (!f) {
            std::cout << "No such file: " << fname.str() << "\n";
        }
        else {
            f << gyro::icntl[Param::optimized] << "\n";
            // f << gyro::icntl[Param::verbose] << "\n";
            // f << gyro::icntl[Param::openmp] << "\n";
            // f << gyro::icntl[Param::matrix_free] << "\n";
            // f << gyro::icntl[Param::debug] << "\n";
            f << gyro::icntl[Param::nr_exp] << "\n";
            f << gyro::icntl[Param::ntheta_exp] << "\n";
            f << gyro::icntl[Param::fac_ani] << "\n";
            f << gyro::icntl[Param::theta_aniso] << "\n";
            f << gyro::icntl[Param::v1] << "\n";
            f << gyro::icntl[Param::v2] << "\n";
            f << gyro::icntl[Param::cycle] << "\n";
            f << gyro::icntl[Param::mod_pk] << "\n";
            f << gyro::icntl[Param::compute_rho] << "\n";
            f << gyro::icntl[Param::level] << "\n";
            f << gyro::icntl[Param::plotit] << "\n";
            f << gyro::icntl[Param::solveit] << "\n";
            f << gyro::icntl[Param::maxiter] << "\n";
            f << gyro::icntl[Param::periodic] << "\n";
            f << gyro::icntl[Param::origin_NOT_coarse] << "\n";
            f << gyro::icntl[Param::smoother] << "\n";
            f << gyro::icntl[Param::discr] << "\n";
            f << gyro::icntl[Param::extrapolation] << "\n";
            f << gyro::icntl[Param::DirBC_Interior] << "\n";
            f << gyro::icntl[Param::paraview] << "\n";
            f << gyro::icntl[Param::divideBy2] << "\n";
            f << gyro::icntl[Param::prob] << "\n";
            f << gyro::icntl[Param::alpha_coeff] << "\n";
            f << gyro::icntl[Param::beta_coeff] << "\n";
            f << gyro::icntl[Param::res_norm] << "\n";
            f << std::setprecision(17) << gyro::dcntl[Param::r0_DB] << "\n";
            f << std::setprecision(17) << gyro::dcntl[Param::R0] << "\n";
            f << std::setprecision(17) << gyro::dcntl[Param::R] << "\n";
            f << std::setprecision(17) << gyro::dcntl[Param::THETA0] << "\n";
            f << std::setprecision(17) << gyro::dcntl[Param::THETA] << "\n";
            f << std::setprecision(17) << gyro::dcntl[Param::kappa_eps] << "\n";
            f << std::setprecision(17) << gyro::dcntl[Param::delta_e] << "\n";
            f << std::setprecision(17) << gyro::dcntl[Param::tol_bound_check] << "\n";
            f << std::setprecision(17) << gyro::dcntl[Param::rel_red_conv] << "\n";
        }
        f.close();
        fname.clear(); //clear any bits set
        fname.str(std::string());
    }

    // std::cout << "No test on the size of the i,j,val vectors: explicit zeros removed in matlab.\n";

    std::vector<double> u_rand0;
    for (int l = 0; l < levels; l++) {
        std::cout << "--------- debug on level " << l << "--------- \n";
        folder << ss.str() << l << "/";
        int m = v_level[l]->m;

        if (gyro::icntl[Param::debug] == 1) {
            const int dir_err2 = mkdir(folder.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            // if (-1 == dir_err2) {
            //     printf("Error creating directory!\n");
            // }
        }

        // Random test vector
        std::vector<double> u_rand;
        fname << folder.str() << "rand.txt";
        unit_test(fname, m, u_rand, "u_rand", 0);
        if (l == 0)
            u_rand0 = u_rand;

        ///////////////////////////////////////////////////////////
        // Grids: r.txt, theta.txt
        ///////////////////////////////////////////////////////////
        std::cout << "Checking grid (r, theta)...\n";
        std::vector<double> r_tmp, theta_tmp;

        fname << folder.str() << "r.txt";
        unit_test(fname, v_level[l]->nr, v_level[l]->r, "v_level[l]->r", 1);

        fname << folder.str() << "theta.txt";
        unit_test(fname, v_level[l]->ntheta, v_level[l]->theta, "v_level[l]->theta", 1);

        ///////////////////////////////////////////////////////////
        // Linear problem:
        ///////////////////////////////////////////////////////////
        // - row_indices.txt, col_indices.txt, vals.txt
        std::cout << "Checking Operator A...\n";
        std::vector<double> uC(m, 0);
        if (gyro::icntl[Param::optimized] == 0)
            v_level[l]->apply_A0(u_rand, uC);
        else
            v_level[l]->apply_A(u_rand, uC);
        fname << folder.str() << "Au.txt";
        unit_test(fname, m, uC, "Au", 1);

        // - row_Ac_LU.txt, col_Ac_LU.txt, vals_Ac_LU.txt
        if (l == levels - 1) {
            std::cout << "Checking Coarsest grid Ac...\n";
            uC = std::vector<double>(m, 0);
            for (std::size_t j = 0; j < u_rand.size(); j++) {
                uC[v_level[l]->row_Ac_LU[j]] += v_level[l]->vals_Ac_LU[j] * u_rand[v_level[l]->col_Ac_LU[j]];
            }
            fname << folder.str() << "Acu.txt";
            unit_test(fname, m, uC, "Acu", 1);

            uC = std::vector<double>(m, 0);
            std::cout << "Checking Operator Ac^-1...\n";
            fname << folder.str() << "Acinvu.txt";
            std::vector<double> Ainv_coarse = v_level[l]->solve_gaussian_elimination(
                v_level[l]->row_Ac_LU, v_level[l]->col_Ac_LU, v_level[l]->vals_Ac_LU, u_rand);
            unit_test(fname, m, uC, "Acinvu", 1);
        }

        // - f.txt
        if (l == 0 || (l == 1 && gyro::icntl[Param::extrapolation] > 0)) {
            std::cout << "Checking RHS...\n";
            fname << folder.str() << "f.txt";
            unit_test(fname, m, v_level[l]->fVec, "v_level[l]->fVec", 1);
        }

        // - betaVec.txt
        if (gyro::icntl[Param::beta_coeff] > 0) {
            std::cout << "Checking Beta vector...\n";
            fname << folder.str() << "betaVec.txt";
            unit_test(fname, v_level[l]->betaVec.size(), v_level[l]->betaVec, "v_level[l]->betaVec", 1);
        }

        ///////////////////////////////////////////////////////////
        // Multigrid (l-1):
        if (l == levels - 1)
            break;
        ///////////////////////////////////////////////////////////
        int mc = v_level[l]->mc;

        // - coarse_nodes_r.txt, coarse_nodes_theta.txt
        std::cout << "Checking coarse grids r...\n";
        fname << folder.str() << "coarse_nodes_r.txt";
        unit_test(fname, v_level[l]->coarse_nodes_list_r.size(), v_level[l]->coarse_nodes_list_r,
                  "v_level[l]->coarse_nodes_list_r", 1);

        std::cout << "Checking coarse grids theta...\n";
        fname << folder.str() << "coarse_nodes_theta.txt";
        unit_test(fname, v_level[l]->coarse_nodes_list_theta.size(), v_level[l]->coarse_nodes_list_theta,
                  "v_level[l]->coarse_nodes_list_theta", 1);

        // - row_indices_prol.txt, col_indices_prol.txt, vals_prol.txt
        std::cout << "Checking Prolongation P...\n";
        fname << folder.str() << "Pu.txt";
        if (gyro::icntl[Param::optimized] == 0)
            uC = v_level[l]->apply_prolongation_bi0(u_rand, mc, m, v_level[l]->coarse_nodes_list_r,
                                                    v_level[l]->coarse_nodes_list_theta, 0);
        else
            uC = v_level[l]->apply_prolongation_bi(u_rand);
        unit_test(fname, m, uC, "Pu", 1);

        if (gyro::icntl[Param::extrapolation] > 0 && l == 0) {
            // - row_indices_prol_ex.txt, col_indices_prol_ex.txt, vals_prol_ex.txt
            std::cout << "Checking Prolongation Pex...\n";
            fname << folder.str() << "Pexu.txt";
            if (gyro::icntl[Param::optimized] == 0)
                uC = v_level[l]->apply_prolongation_ex0(u_rand, mc, m, v_level[l]->coarse_nodes_list_r,
                                                        v_level[l]->coarse_nodes_list_theta, 0);
            else
                uC = v_level[l]->apply_prolongation_ex(u_rand);
            unit_test(fname, m, uC, "Pexu", 1);

            // - row_indices_prol_inj.txt, col_indices_prol_inj.txt, vals_prol_inj.txt
            std::cout << "Checking Prolongation Pinj...\n";
            fname << folder.str() << "Pinju.txt";
            if (gyro::icntl[Param::optimized] == 0)
                uC = v_level[l]->apply_prolongation_inj0(u_rand, mc, m, v_level[l]->coarse_nodes_list_r,
                                                         v_level[l]->coarse_nodes_list_theta, 0);
            else
                uC = v_level[l]->apply_prolongation_inj(u_rand);
            unit_test(fname, m, uC, "Pinju", 1);
        }

        // - row_indices_Asc_s_c.txt, col_indices_Asc_s_c.txt, vals_Asc_s_c.txt
        std::cout << "Checking Smoother Matrices Asc...\n";
        uC.clear();
        int start = 0, end = 0;
        std::vector<double> Ascu;
        for (int smoother = 0; smoother < 2; smoother++) {
            for (int color = 0; color < 2; color++) {
                int smoother_ind = smoother * 2 + color;
                int m_sc         = v_level[l]->m_sc[smoother_ind];
                for (int ij = 0; ij < v_level[l]->nblocks[smoother_ind]; ij++) {
                    start = end;
                    end += m_sc;
                    std::vector<double>::const_iterator first = u_rand.begin() + start;
                    std::vector<double>::const_iterator last  = u_rand.begin() + end;
                    std::vector<double> u_rand_sc(first, last);
                    Ascu = std::vector<double>(m_sc, 0);
                    for (int j = 0; j < v_level[l]->A_Zebra_r_LU_row[smoother_ind][ij].size(); j++) {
                        Ascu[v_level[l]->A_Zebra_r_LU_row[smoother_ind][ij][j]] +=
                            v_level[l]->A_Zebra_v_LU_row[smoother_ind][ij][j] *
                            u_rand_sc[v_level[l]->A_Zebra_c_LU_row[smoother_ind][ij][j]];
                    }
                    uC.insert(uC.end(), Ascu.begin(), Ascu.end());
                }
            }
        }
        fname << folder.str() << "Ascu.txt";
        unit_test(fname, m, uC, "Ascu", 1);

        std::cout << "Checking Asc^-1...\n";
        uC.clear();
        start = 0;
        end   = 0;
        std::vector<double> Ascinvu;
        for (int smoother = 0; smoother < 2; smoother++) {
            for (int color = 0; color < 2; color++) {
                int smoother_ind = smoother * 2 + color;
                int m_sc         = v_level[l]->m_sc[smoother_ind];
                for (int ij = 0; ij < v_level[l]->nblocks[smoother_ind]; ij++) {
                    start = end;
                    end += m_sc;
                    std::vector<double>::const_iterator first = u_rand.begin() + start;
                    std::vector<double>::const_iterator last  = u_rand.begin() + end;
                    std::vector<double> u_rand_sc(first, last);
                    if (v_level[l]->A_Zebra_r_row[smoother_ind][ij].size() == 0) {
                        continue;
                    }
                    // Dirichlet or size of system is 1 (diagonal solve)
                    if ((smoother_ind == 0 && ij == 0 && gyro::icntl[Param::DirBC_Interior]) || m_sc == 1 ||
                        (gyro::icntl[Param::extrapolation] && l == 0 && smoother_ind % 2 == 0 &&
                         !(smoother_ind == 0 && !gyro::icntl[Param::DirBC_Interior] && ij == 0))) {
                        Ascinvu = v_level[l]->solve_diag(v_level[l]->A_Zebra_v_row[smoother_ind][ij], u_rand_sc);
                    }
                    // Circle (not across)
                    else if (smoother_ind < 2 &&
                             !(smoother_ind == 0 && !gyro::icntl[Param::DirBC_Interior] && ij == 0)) {
                        Ascinvu = v_level[l]->solve_circle(v_level[l]->A_Zebra_r_LU_row[smoother_ind][ij],
                                                           v_level[l]->A_Zebra_c_LU_row[smoother_ind][ij],
                                                           v_level[l]->A_Zebra_v_LU_row[smoother_ind][ij], u_rand_sc);
                    }
                    // Radial (not across)
                    else if (smoother_ind > 1) {
                        Ascinvu = v_level[l]->solve_radial(v_level[l]->A_Zebra_r_LU_row[smoother_ind][ij],
                                                           v_level[l]->A_Zebra_c_LU_row[smoother_ind][ij],
                                                           v_level[l]->A_Zebra_v_LU_row[smoother_ind][ij], u_rand_sc);
                    }
                    // Across (direct solver)
                    else {
                        Ascinvu = v_level[l]->solve_gaussian_elimination(v_level[l]->A_Zebra_r_LU_row[smoother_ind][ij],
                                                                         v_level[l]->A_Zebra_c_LU_row[smoother_ind][ij],
                                                                         v_level[l]->A_Zebra_v_LU_row[smoother_ind][ij],
                                                                         u_rand_sc);
                    }
                    uC.insert(uC.end(), Ascinvu.begin(), Ascinvu.end());
                }
            }
        }
        fname << folder.str() << "Ascinvu.txt";
        unit_test(fname, m, uC, "Ascinvu", 1);

        // - row_indices_Asc_ortho_s_c.txt, col_indices_Asc_ortho_s_c.txt, vals_Asc_ortho_s_c.txt
        std::cout << "Checking Smoother Complement Matrices Asc_ortho...\n";
        uC.clear();
        for (int smoother = 0; smoother < 2; smoother++) {
            for (int color = 0; color < 2; color++) {
                int smoother_ind = smoother * 2 + color;
                int m_sc         = v_level[l]->nblocks[smoother_ind] * v_level[l]->m_sc[smoother_ind];
                std::vector<double> Asc_orthou(m_sc, 0);
                v_level[l]->u = u_rand;
                if (gyro::icntl[Param::optimized] == 0)
                    v_level[l]->apply_Asc_ortho0(Asc_orthou, smoother_ind);
                else {
                    int sm = (smoother_ind < 2) ? 1 - smoother_ind : 5 - smoother_ind;
                    if (gyro::icntl[Param::smoother] == 13) {
                        if (smoother < 2) { //circle
                            v_level[l]->apply_Asc_ortho(Asc_orthou, v_level[l]->u, smoother_ind, 0, 0,
                                                        v_level[l]->dep_Asc[smoother_ind], v_level[l]->dep_Asc[sm],
                                                        v_level[l]->dep_Asc[1], v_level[l]->dep_Asc[smoother_ind]);
                        }
                        else { //radial
                            v_level[l]->apply_Asc_ortho(Asc_orthou, v_level[l]->u, smoother_ind, 0, 0,
                                                        v_level[l]->dep_Asc[smoother_ind], v_level[l]->dep_Asc[sm],
                                                        v_level[l]->dep_Asc[1], v_level[l]->dep_Asc[smoother_ind]);
                        }
                    }
                    else {
                        v_level[l]->apply_Asc_ortho(Asc_orthou, v_level[l]->u, smoother_ind, 0, 0,
                                                    v_level[l]->dep_Asc[smoother_ind], v_level[l]->dep_Asc[sm],
                                                    v_level[l]->dep_Asc[1], v_level[l]->dep_Asc_ortho[smoother_ind]);
                    }
                }
                uC.insert(uC.end(), Asc_orthou.begin(), Asc_orthou.end());
            }
        }
        fname << folder.str() << "Asc_orthou.txt";
        unit_test(fname, m, uC, "Asc_orthou", 1);

        folder.clear(); //clear any bits set
        folder.str(std::string());
    }

    ///////////////////////////////////////////////////////////
    // Checking whole multigrid cycle
    ///////////////////////////////////////////////////////////
    std::cout << "Checking multigrid iteration...\n";
    v_level[0]->u = u_rand0;
    if (gyro::icntl[Param::extrapolation] > 0)
        v_level[1]->fVec_initial = v_level[1]->fVec; //Extrapolation: store the initial fVec on level 1
    multigrid_cycle_extrapol(0);
    fname << ss.str() << "0/MG_cycle_rand.txt";
    unit_test(fname, v_level[0]->m, v_level[0]->u, "MG_cycle_rand", 1);

    std::cout << "Checking residual...\n";
    std::vector<double> res_M;
    compute_residual(0, gyro::icntl[Param::extrapolation]);
    fname << ss.str() << "0/error_residual.txt";
    unit_test(fname, v_level[0]->m, v_level[0]->res, "residual", 1);

    std::cout << "Checking error...\n";
    std::vector<double> error_M;
    std::vector<double> error = compute_error();
    fname << ss.str() << "0/error_rand.txt";
    unit_test(fname, v_level[0]->m, error, "error", 1);
} /* ----- end of gmgpolar::debug ----- */
