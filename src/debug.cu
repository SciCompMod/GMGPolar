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

#ifdef CUDA

    /*!
 * \file debug.cpp
 * \brief Comparison of function results with the matlab code (A, f, P, smoothers)
 * \author M. Kuehn, C. Kruse, P. Leleux, C. Schwarz
 * \version 0.0
 */
    #include "gmgpolar.h"
    #include <sstream>
    #include <fstream>
    #include <assert.h>
    #include "cuda.h"

    #define NTHREADS 256

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
    std::fstream f;
    int fac_ani        = gyro::icntl[Param::fac_ani];
    int mod_pk         = gyro::icntl[Param::mod_pk];
    int DirBC_Interior = gyro::icntl[Param::DirBC_Interior];
    ss << "outputs_MATLAB/" << fac_ani << "_" << mod_pk << "_" << DirBC_Interior << "_0/";

    double *dev_uc, *dev_uc2, *dev_usc, *dev_Pu, *dev_Ru;

    ///////////////////////////////////////////////////////////
    // Parameters: params.txt (compute_rho, cycle, delta_e, fac_ani, kappa_eps, levels, maxiter, mod_pk, nr, ntheta, periodic, plotit, R0, r0_DB, rDBC(1), rDBC(2), Rmax, solveit, v1, v2)
    ///////////////////////////////////////////////////////////
    std::cout << "Checking parameters...\n";
    fname << ss.str() << "params.txt";
    f.open(fname.str().c_str(), std::ios::in);
    if (!f) {
        std::cout << "No such file\n";
    }
    else {
        std::vector<double> parray;
        while (1) {
            double p;
            f >> p;
            parray.push_back(p);
            if (f.eof())
                break;
        }

        assert(gyro::icntl[Param::compute_rho] == parray[0]);
        assert(gyro::icntl[Param::cycle] == parray[1]);
        assert(gyro::dcntl[Param::delta_e] == parray[2]);
        assert(gyro::icntl[Param::fac_ani] == parray[3]);
        assert(gyro::dcntl[Param::kappa_eps] == parray[4]);
        assert(levels == parray[5]);
        assert(gyro::icntl[Param::maxiter] == parray[6]);
        assert(gyro::icntl[Param::mod_pk] == parray[7]);
        assert(v_level[0]->nr - 1 == parray[8]);
        assert(v_level[0]->ntheta - 1 == parray[9]);
        assert(gyro::icntl[Param::periodic] == parray[10]);
        assert(gyro::icntl[Param::plotit] == parray[11]);
        assert(gyro::dcntl[Param::R0] == parray[12]);
        assert(gyro::dcntl[Param::r0_DB] == parray[13]);
        assert(gyro::dcntl[Param::R] == parray[14]);
        assert(gyro::icntl[Param::solveit] == parray[15]);
        assert(gyro::icntl[Param::v1] == parray[16]);
        assert(gyro::icntl[Param::v2] == parray[17]);
    }
    f.close();
    fname.clear(); //clear any bits set
    fname.str(std::string());

    // std::cout << "No test on the size of the i,j,val vectors: explicit zeros removed in matlab.\n";

    std::vector<double> u_rand0;
    for (int l = 0; l < levels; l++) {
        v_level[l]->res = std::vector<double>(v_level[l]->m);
        int n_blocks    = ceil((double)v_level[l]->m / (double)NTHREADS);

        if (l == 0)
            cudaMemcpy(v_level[l]->dev_fVec, &(v_level[l]->fVec[0]), v_level[l]->m * sizeof(double),
                       cudaMemcpyHostToDevice);
        cudaMemcpy(v_level[l]->dev_res, &(v_level[l]->res[0]), v_level[l]->m * sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(v_level[l]->dev_r, &(v_level[l]->r[0]), v_level[l]->r.size() * sizeof(double),
                   cudaMemcpyHostToDevice);
        cudaMemcpy(v_level[l]->dev_hplus, &(v_level[l]->hplus[0]), v_level[l]->hplus.size() * sizeof(double),
                   cudaMemcpyHostToDevice);
        cudaMemcpy(v_level[l]->dev_theta, &(v_level[l]->theta[0]), v_level[l]->theta.size() * sizeof(double),
                   cudaMemcpyHostToDevice);
        cudaMemcpy(v_level[l]->dev_thetaplus, &(v_level[l]->thetaplus[0]),
                   v_level[l]->thetaplus.size() * sizeof(double), cudaMemcpyHostToDevice);

        std::cout << "--------- debug on level " << l << "--------- \n";
        folder << ss.str() << l + 1 << "/";

        // Random test vector
        std::vector<double> u_rand;
        fname << folder.str() << "rand.txt";
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
                u_rand.push_back(p);
            }
        }
        f.close();
        fname.clear();
        fname.str(std::string());
        if (l == 0)
            u_rand0 = u_rand;

        ///////////////////////////////////////////////////////////
        // Grids: r.txt, theta.txt
        ///////////////////////////////////////////////////////////
        std::cout << "Checking grid (r, theta)...\n";
        std::vector<double> r_tmp, theta_tmp;

        fname << folder.str() << "r.txt";
        f.open(fname.str().c_str(), std::ios::in);
        if (!f) {
            std::cout << "No such file\n";
        }
        else {
            while (1) {
                double p;
                f >> p;
                if (f.eof())
                    break;
                r_tmp.push_back(p);
            }
            assert(r_tmp.size() == v_level[l]->nr);
            for (int i = 0; i < r_tmp.size(); i++) {
                if (r_tmp[i] - v_level[l]->r[i] > 1e-4)
                    assert(r_tmp[i] - v_level[l]->r[i] < 1e-12);
            }
        }
        f.close();
        fname.clear();
        fname.str(std::string());

        fname << folder.str() << "theta.txt";
        f.open(fname.str().c_str(), std::ios::in);
        if (!f) {
            std::cout << "No such file\n";
        }
        else {
            while (1) {
                double p;
                f >> p;
                if (f.eof())
                    break;
                theta_tmp.push_back(p);
            }
            assert(theta_tmp.size() == v_level[l]->ntheta);
            for (int i = 0; i < theta_tmp.size(); i++) {
                if (theta_tmp[i] - v_level[l]->theta[i] > 1e-4)
                    std::cout << "theta_tmp[i]: " << theta_tmp[i] << ", v_level[l]->theta[i]: " << v_level[l]->theta[i]
                              << "\n";
                assert(theta_tmp[i] - v_level[l]->theta[i] < 1e-4);
            }
        }
        f.close();
        fname.clear();
        fname.str(std::string());

        ///////////////////////////////////////////////////////////
        // Linear problem :
        ///////////////////////////////////////////////////////////
        std::cout << "Checking Operator A...\n";
        // -row_indices.txt, col_indices.txt, vals.txt
        std::vector<int> ri, ci;
        std::vector<double> v;

        fname << folder.str() << "row_indices.txt";
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
                ri.push_back(p - 1);
            }
        }
        f.close();
        fname.clear();
        fname.str(std::string());

        fname << folder.str() << "col_indices.txt";
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
                ci.push_back(p - 1);
            }
        }
        f.close();
        fname.clear();
        fname.str(std::string());

        fname << folder.str() << "vals.txt";
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
                v.push_back(p);
            }
        }
        f.close();
        fname.clear();
        fname.str(std::string());

        int m = v_level[l]->m;

        for (int i = 0; i < m; i++) {
            std::vector<double> ei(m, 0), uC(m, 0), uM(m, 0);
            ei[i] = 1;
            uC    = std::vector<double>(m, 0);

            cudaMemcpy(v_level[l]->dev_u, &(ei[0]), v_level[l]->m * sizeof(double), cudaMemcpyHostToDevice);
            cudaMemcpy(v_level[l]->dev_Au, &uC[0], v_level[l]->m * sizeof(double), cudaMemcpyHostToDevice);

            apply_A_SA<<<n_blocks, NTHREADS>>>(
                v_level[l]->m, v_level[l]->dev_u, v_level[l]->dev_Au, v_level[l]->nr_int, v_level[l]->dev_r,
                v_level[l]->dev_hplus, v_level[l]->ntheta_int, v_level[l]->dev_theta, v_level[l]->dev_thetaplus,
                gyro::icntl[Param::DirBC_Interior], gyro::icntl[Param::mod_pk], gyro::icntl[Param::prob],
                gyro::dcntl[Param::kappa_eps], gyro::dcntl[Param::delta_e], gyro::dcntl[Param::R]);

            cudaMemcpy(&(uC[0]), (void*)v_level[l]->dev_Au, v_level[l]->m * sizeof(double), cudaMemcpyDeviceToHost);

            for (int j = 0; j < ri.size(); j++) {
                uM[ri[j]] += v[j] * ei[ci[j]];
            }
            for (int j = 0; j < m; j++) {
                if (uC[j] - uM[j] > 1e-3)
                    std::cout << "j: " << j << ", uC[j]: " << uC[j] << ", uM[j]: " << uM[j] << "\n";
                assert(uC[j] - uM[j] < 1e-3);
            }
        }

        if (l == levels) {
            std::cout << "Checking Operator Ac^-1...\n";
            std::vector<double> Ainv_coarse_M;
            fname << folder.str() << "invA_coarse_rand.txt";
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
                    Ainv_coarse_M.push_back(p);
                }
            }
            f.close();
            fname.clear();
            fname.str(std::string());
            std::vector<double> Ainv_coarse = v_level[l]->solve_gaussian_elimination(
                v_level[l]->row_Ac_LU, v_level[l]->col_Ac_LU, v_level[l]->vals_Ac_LU, u_rand);
            for (int j = 0; j < m; j++) {
                if (Ainv_coarse[j] - Ainv_coarse_M[j] > 1e-3)
                    std::cout << "j: " << j << ", Ainv_coarse[j]: " << Ainv_coarse[j]
                              << ", Ainv_coarse_M[j]: " << Ainv_coarse_M[j] << "\n";
                assert(Ainv_coarse[j] - Ainv_coarse_M[j] < 1e-3);
            }
        }

        // - f.txt
        if (l == 0) {
            std::cout << "Checking RHS...\n";
            std::vector<double> rhs;
            fname << folder.str() << "f.txt";
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
                    rhs.push_back(p);
                }
                for (int i = 0; i < v_level[l]->fVec.size(); i++) {
                    if (rhs[i] - v_level[l]->fVec[i] > 1e-2)
                        std::cout << "i: " << i << ", rhs[i]: " << rhs[i]
                                  << ", v_level[l]->fVec[i]: " << v_level[l]->fVec[i] << "\n";
                    assert(rhs[i] - v_level[l]->fVec[i] < 1e-2);
                }
            }
            f.close();
            fname.clear();
            fname.str(std::string());
        }

        ///////////////////////////////////////////////////////////
        // Multigrid (l-1):
        if (l == levels - 1) {
            break;
        }
        ///////////////////////////////////////////////////////////
        int mc = v_level[l]->mc;

        cudaMalloc((void**)&dev_uc, v_level[l]->mc * sizeof(double));
        cudaMalloc((void**)&dev_uc2, v_level[l]->m * sizeof(double));
        cudaMalloc((void**)&dev_Pu, v_level[l]->m * sizeof(double));
        cudaMalloc((void**)&dev_Ru, v_level[l]->mc * sizeof(double));

        cudaMemcpy(v_level[l]->dev_coarse_nodes_list_r, &(v_level[l]->coarse_nodes_list_r[0]),
                   v_level[l]->coarse_nodes_list_r.size() * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(v_level[l]->dev_coarse_nodes_list_theta, &(v_level[l]->coarse_nodes_list_theta[0]),
                   v_level[l]->coarse_nodes_list_theta.size() * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(v_level[l]->dev_coarse_nodes_list_type_shift, &(v_level[l]->coarse_nodes_list_type_shift[0]),
                   sizeof(int) * v_level[l]->coarse_nodes_list_type_shift.size(), cudaMemcpyHostToDevice);
        cudaMemcpy(v_level[l]->dev_coarse_nodes_list_type_start, &(v_level[l]->coarse_nodes_list_type_start[0]),
                   sizeof(int) * v_level[l]->coarse_nodes_list_type_start.size(), cudaMemcpyHostToDevice);
        cudaMemcpy(v_level[l]->dev_coarse_nodes_list_type0, &(v_level[l]->coarse_nodes_list_type[0][0]),
                   sizeof(int) * v_level[l]->coarse_nodes_list_type[0].size(), cudaMemcpyHostToDevice);
        cudaMemcpy(v_level[l]->dev_coarse_nodes_list_type1, &(v_level[l]->coarse_nodes_list_type[1][0]),
                   sizeof(int) * v_level[l]->coarse_nodes_list_type[1].size(), cudaMemcpyHostToDevice);
        cudaMemcpy(v_level[l]->dev_coarse_nodes_list_type2, &(v_level[l]->coarse_nodes_list_type[2][0]),
                   sizeof(int) * v_level[l]->coarse_nodes_list_type[2].size(), cudaMemcpyHostToDevice);
        cudaMemcpy(v_level[l]->dev_coarse_nodes_list_type3, &(v_level[l]->coarse_nodes_list_type[3][0]),
                   sizeof(int) * v_level[l]->coarse_nodes_list_type[3].size(), cudaMemcpyHostToDevice);
        cudaMemcpy(v_level[l]->dev_coarse_nodes_list_type4, &(v_level[l]->coarse_nodes_list_type[4][0]),
                   sizeof(int) * v_level[l]->coarse_nodes_list_type[4].size(), cudaMemcpyHostToDevice);
        cudaMemcpy(v_level[l]->dev_coarse_nodes_list_type5, &(v_level[l]->coarse_nodes_list_type[5][0]),
                   sizeof(int) * v_level[l]->coarse_nodes_list_type[5].size(), cudaMemcpyHostToDevice);
        cudaMemcpy(v_level[l]->dev_coarse_nodes_list_type6, &(v_level[l]->coarse_nodes_list_type[6][0]),
                   sizeof(int) * v_level[l]->coarse_nodes_list_type[6].size(), cudaMemcpyHostToDevice);
        cudaMemcpy(v_level[l]->dev_coarse_nodes_list_type7, &(v_level[l]->coarse_nodes_list_type[7][0]),
                   sizeof(int) * v_level[l]->coarse_nodes_list_type[7].size(), cudaMemcpyHostToDevice);
        cudaMemcpy(v_level[l]->dev_coarse_nodes_list_i0, &(v_level[l]->coarse_nodes_list_i[0][0]),
                   sizeof(int) * v_level[l]->coarse_nodes_list_i[0].size(), cudaMemcpyHostToDevice);
        cudaMemcpy(v_level[l]->dev_coarse_nodes_list_i1, &(v_level[l]->coarse_nodes_list_i[1][0]),
                   sizeof(int) * v_level[l]->coarse_nodes_list_i[1].size(), cudaMemcpyHostToDevice);
        cudaMemcpy(v_level[l]->dev_coarse_nodes_list_i2, &(v_level[l]->coarse_nodes_list_i[2][0]),
                   sizeof(int) * v_level[l]->coarse_nodes_list_i[2].size(), cudaMemcpyHostToDevice);
        cudaMemcpy(v_level[l]->dev_coarse_nodes_list_i3, &(v_level[l]->coarse_nodes_list_i[3][0]),
                   sizeof(int) * v_level[l]->coarse_nodes_list_i[3].size(), cudaMemcpyHostToDevice);
        cudaMemcpy(v_level[l]->dev_coarse_nodes_list_i4, &(v_level[l]->coarse_nodes_list_i[4][0]),
                   sizeof(int) * v_level[l]->coarse_nodes_list_i[4].size(), cudaMemcpyHostToDevice);
        cudaMemcpy(v_level[l]->dev_coarse_nodes_list_i5, &(v_level[l]->coarse_nodes_list_i[5][0]),
                   sizeof(int) * v_level[l]->coarse_nodes_list_i[5].size(), cudaMemcpyHostToDevice);
        cudaMemcpy(v_level[l]->dev_coarse_nodes_list_i6, &(v_level[l]->coarse_nodes_list_i[6][0]),
                   sizeof(int) * v_level[l]->coarse_nodes_list_i[6].size(), cudaMemcpyHostToDevice);
        cudaMemcpy(v_level[l]->dev_coarse_nodes_list_i7, &(v_level[l]->coarse_nodes_list_i[7][0]),
                   sizeof(int) * v_level[l]->coarse_nodes_list_i[7].size(), cudaMemcpyHostToDevice);
        cudaMemcpy(v_level[l]->dev_coarse_nodes_list_j0, &(v_level[l]->coarse_nodes_list_j[0][0]),
                   sizeof(int) * v_level[l]->coarse_nodes_list_j[0].size(), cudaMemcpyHostToDevice);
        cudaMemcpy(v_level[l]->dev_coarse_nodes_list_j1, &(v_level[l]->coarse_nodes_list_j[1][0]),
                   sizeof(int) * v_level[l]->coarse_nodes_list_j[1].size(), cudaMemcpyHostToDevice);
        cudaMemcpy(v_level[l]->dev_coarse_nodes_list_j2, &(v_level[l]->coarse_nodes_list_j[2][0]),
                   sizeof(int) * v_level[l]->coarse_nodes_list_j[2].size(), cudaMemcpyHostToDevice);
        cudaMemcpy(v_level[l]->dev_coarse_nodes_list_j3, &(v_level[l]->coarse_nodes_list_j[3][0]),
                   sizeof(int) * v_level[l]->coarse_nodes_list_j[3].size(), cudaMemcpyHostToDevice);
        cudaMemcpy(v_level[l]->dev_coarse_nodes_list_j4, &(v_level[l]->coarse_nodes_list_j[4][0]),
                   sizeof(int) * v_level[l]->coarse_nodes_list_j[4].size(), cudaMemcpyHostToDevice);
        cudaMemcpy(v_level[l]->dev_coarse_nodes_list_j5, &(v_level[l]->coarse_nodes_list_j[5][0]),
                   sizeof(int) * v_level[l]->coarse_nodes_list_j[5].size(), cudaMemcpyHostToDevice);
        cudaMemcpy(v_level[l]->dev_coarse_nodes_list_j6, &(v_level[l]->coarse_nodes_list_j[6][0]),
                   sizeof(int) * v_level[l]->coarse_nodes_list_j[6].size(), cudaMemcpyHostToDevice);
        cudaMemcpy(v_level[l]->dev_coarse_nodes_list_j7, &(v_level[l]->coarse_nodes_list_j[7][0]),
                   sizeof(int) * v_level[l]->coarse_nodes_list_j[7].size(), cudaMemcpyHostToDevice);

        // - row_indices_prol.txt, col_indices_prol.txt, vals_prol.txt
        std::cout << "Checking Prolongation P...\n";
        ri.clear();
        ci.clear();
        v.clear();
        fname << folder.str() << "row_indices_prol.txt";
        f.open(fname.str().c_str(), std::ios::in);
        if (!f) {
            std::cout << "No such file\n";
        }
        else {
            while (1) {
                int p;
                f >> p;
                if (f.eof())
                    break;
                ri.push_back(p - 1);
            }
        }
        f.close();
        fname.clear();
        fname.str(std::string());

        fname << folder.str() << "col_indices_prol.txt";
        f.open(fname.str().c_str(), std::ios::in);
        if (!f) {
            std::cout << "No such file\n";
        }
        else {
            while (1) {
                int p;
                f >> p;
                if (f.eof())
                    break;
                ci.push_back(p - 1);
            }
        }
        f.close();
        fname.clear();
        fname.str(std::string());

        fname << folder.str() << "vals_prol.txt";
        f.open(fname.str().c_str(), std::ios::in);
        if (!f) {
            std::cout << "No such file\n";
        }
        else {
            while (1) {
                double p;
                f >> p;
                if (f.eof())
                    break;
                v.push_back(p);
            }
        }
        f.close();
        fname.clear();
        fname.str(std::string());

        for (int i = 0; i < mc; i++) {
            std::vector<double> ei(mc, 0), uC(m, 0), uM(m, 0);
            ei[i] = 1;
            uC    = std::vector<double>(m, 0);

            cudaMemcpy(dev_uc, &ei[0], v_level[l]->mc * sizeof(double), cudaMemcpyHostToDevice);
            cudaMemcpy(dev_Pu, &uC[0], v_level[l]->m * sizeof(double), cudaMemcpyHostToDevice);

            apply_prolongation_bi_gpu<<<n_blocks, NTHREADS>>>(
                dev_Pu, dev_uc, v_level[l]->m, v_level[l]->mc, v_level[l]->dev_thetaplus, v_level[l]->ntheta_int,
                v_level[l]->dev_hplus, v_level[l]->dev_coarse_nodes_list_type_shift,
                v_level[l]->dev_coarse_nodes_list_type_start, v_level[l]->coarse_nodes_list_i[0].size(),
                v_level[l]->coarse_nodes_list_i[1].size(), v_level[l]->coarse_nodes_list_i[2].size(),
                v_level[l]->coarse_nodes_list_i[3].size(), v_level[l]->coarse_nodes_list_i[4].size(),
                v_level[l]->coarse_nodes_list_i[5].size(), v_level[l]->coarse_nodes_list_i[6].size(),
                v_level[l]->coarse_nodes_list_i[7].size(), v_level[l]->dev_coarse_nodes_list_type0,
                v_level[l]->dev_coarse_nodes_list_type1, v_level[l]->dev_coarse_nodes_list_type2,
                v_level[l]->dev_coarse_nodes_list_type3, v_level[l]->dev_coarse_nodes_list_type4,
                v_level[l]->dev_coarse_nodes_list_type5, v_level[l]->dev_coarse_nodes_list_type6,
                v_level[l]->dev_coarse_nodes_list_type7, v_level[l]->dev_coarse_nodes_list_i0,
                v_level[l]->dev_coarse_nodes_list_i1, v_level[l]->dev_coarse_nodes_list_i2,
                v_level[l]->dev_coarse_nodes_list_i3, v_level[l]->dev_coarse_nodes_list_i4,
                v_level[l]->dev_coarse_nodes_list_i5, v_level[l]->dev_coarse_nodes_list_i6,
                v_level[l]->dev_coarse_nodes_list_i7, v_level[l]->dev_coarse_nodes_list_j0,
                v_level[l]->dev_coarse_nodes_list_j1, v_level[l]->dev_coarse_nodes_list_j2,
                v_level[l]->dev_coarse_nodes_list_j3, v_level[l]->dev_coarse_nodes_list_j4,
                v_level[l]->dev_coarse_nodes_list_j5, v_level[l]->dev_coarse_nodes_list_j6,
                v_level[l]->dev_coarse_nodes_list_j7);

            cudaMemcpy(&uC[0], dev_Pu, v_level[l]->m * sizeof(double), cudaMemcpyDeviceToHost);

            // uC    = v_level[l]->apply_prolongation_bi(ei, 0);
            for (int j = 0; j < ri.size(); j++) {
                uM[ri[j]] += v[j] * ei[ci[j]];
            }
            for (int j = 0; j < m; j++) {
                if (uC[j] - uM[j] > 1e-3)
                    std::cout << "j: " << j << ", uC[j]: " << uC[j] << ", uM[j]: " << uM[j] << "\n";
                assert(uC[j] - uM[j] < 1e-3);
            }
        }
        std::cout << "Checking Restriction R...\n";
        // for (int i = 0; i < m; i++) {
        std::vector<double> ei(m, 0), uC(mc, 0), uM(mc, 0);
        // ei[i] = 1;
        uC = std::vector<double>(mc, 0);

        // cudaMemcpy(dev_uc2, &ei[0], v_level[l]->m * sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(dev_uc2, &u_rand[0], v_level[l]->m * sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(dev_Ru, &uC[0], v_level[l]->mc * sizeof(double), cudaMemcpyHostToDevice);

        apply_restriction_bi_gpu<<<n_blocks, NTHREADS>>>(
            dev_Ru, dev_uc2, v_level[l]->m, v_level[l]->mc, v_level[l]->dev_thetaplus, v_level[l]->ntheta_int,
            v_level[l]->dev_hplus, v_level[l]->dev_coarse_nodes_list_type_shift,
            v_level[l]->dev_coarse_nodes_list_type_start, v_level[l]->coarse_nodes_list_type[0].size(),
            v_level[l]->coarse_nodes_list_type[1].size(), v_level[l]->coarse_nodes_list_type[2].size(),
            v_level[l]->coarse_nodes_list_type[3].size(), v_level[l]->coarse_nodes_list_type[4].size(),
            v_level[l]->coarse_nodes_list_type[5].size(), v_level[l]->coarse_nodes_list_type[6].size(),
            v_level[l]->coarse_nodes_list_type[7].size(), v_level[l]->dev_coarse_nodes_list_type0,
            v_level[l]->dev_coarse_nodes_list_type1, v_level[l]->dev_coarse_nodes_list_type2,
            v_level[l]->dev_coarse_nodes_list_type3, v_level[l]->dev_coarse_nodes_list_type4,
            v_level[l]->dev_coarse_nodes_list_type5, v_level[l]->dev_coarse_nodes_list_type6,
            v_level[l]->dev_coarse_nodes_list_type7, v_level[l]->dev_coarse_nodes_list_i0,
            v_level[l]->dev_coarse_nodes_list_i1, v_level[l]->dev_coarse_nodes_list_i2,
            v_level[l]->dev_coarse_nodes_list_i3, v_level[l]->dev_coarse_nodes_list_i4,
            v_level[l]->dev_coarse_nodes_list_i5, v_level[l]->dev_coarse_nodes_list_i6,
            v_level[l]->dev_coarse_nodes_list_i7, v_level[l]->dev_coarse_nodes_list_j0,
            v_level[l]->dev_coarse_nodes_list_j1, v_level[l]->dev_coarse_nodes_list_j2,
            v_level[l]->dev_coarse_nodes_list_j3, v_level[l]->dev_coarse_nodes_list_j4,
            v_level[l]->dev_coarse_nodes_list_j5, v_level[l]->dev_coarse_nodes_list_j6,
            v_level[l]->dev_coarse_nodes_list_j7);

        cudaMemcpy(&uC[0], dev_Ru, v_level[l]->mc * sizeof(double), cudaMemcpyDeviceToHost);

        // uC    = v_level[l]->apply_prolongation_bi(ei, 0);
        for (int j = 0; j < ri.size(); j++) {
            // uM[ci[j]] += v[j] * ei[ri[j]];
            uM[ci[j]] += v[j] * u_rand[ri[j]];
        }
        for (int j = 0; j < mc; j++) {
            if (uC[j] - uM[j] > 1e-3)
                std::cout << "j: " << j << ", uC[j]: " << uC[j] << ", uM[j]: " << uM[j] << "\n";
            assert(uC[j] - uM[j] < 1e-3);
        }
        // }

        std::cout << "Checking Smoother Matrices Asc...\n";
        // - row_indices_Asc_s_c.txt, col_indices_Asc_s_c.txt, vals_Asc_s_c.txt
        for (int smoother = 0; smoother < 2; smoother++) {
            for (int color = 0; color < 2; color++) {
                std::cout << "\t...Asc\n";
                std::vector<int> ri, ci;
                std::vector<double> v;

                fname << folder.str() << "row_indices_Asc_" << smoother + 1 << "_" << color + 1 << ".txt";
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
                        ri.push_back(p - 1);
                    }
                }
                f.close();
                fname.clear();
                fname.str(std::string());

                fname << folder.str() << "col_indices_Asc_" << smoother + 1 << "_" << color + 1 << ".txt";
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
                        ci.push_back(p - 1);
                    }
                }
                f.close();
                fname.clear();
                fname.str(std::string());

                fname << folder.str() << "vals_Asc_" << smoother + 1 << "_" << color + 1 << ".txt";
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
                        v.push_back(p);
                    }
                }
                f.close();
                fname.clear();
                fname.str(std::string());

                int smoother_ind = smoother * 2 + color;
                int m_sc         = v_level[l]->m_sc[smoother_ind];
                for (int i = 0; i < m_sc; i++) {
                    std::vector<double> ei(m_sc, 0), uC(m_sc, 0), uM(m_sc, 0);
                    ei[i] = 1;
                    for (int j = 0; j < v_level[l]->A_Zebra_r[smoother_ind].size(); j++) {
                        uC[v_level[l]->A_Zebra_r[smoother_ind][j]] +=
                            v_level[l]->A_Zebra_v[smoother_ind][j] * ei[v_level[l]->A_Zebra_c[smoother_ind][j]];
                    }
                    for (int j = 0; j < ri.size(); j++) {
                        uM[ri[j]] += v[j] * ei[ci[j]];
                    }
                    for (int j = 0; j < m_sc; j++) {
                        if (uC[j] - uM[j] > 1e-3)
                            std::cout << "j: " << j << ", uC[j]: " << uC[j] << ", uM[j]: " << uM[j] << "\n";
                        assert(uC[j] - uM[j] < 1e-3);
                    }
                }

                std::cout << "\t...fsc\n";
                std::vector<double> fsc_M;
                fname << folder.str() << "fsc_rand_" << smoother + 1 << "_" << color + 1 << ".txt";
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
                        fsc_M.push_back(p);
                    }
                }
                f.close();
                fname.clear();
                fname.str(std::string());
                std::vector<double> fVec_orig = v_level[l]->fVec;
                v_level[l]->fVec              = u_rand;
                std::vector<double> f_sc;
                v_level[l]->build_fsc(f_sc, smoother_ind);
                for (int i = 0; i < m_sc; i++) {
                    assert(f_sc[i] - fsc_M[i] < 1e-3);
                }
                v_level[l]->fVec = fVec_orig;

                std::cout << "\t...Asc^-1\n";
                std::vector<double> u_sc_M;
                fname << folder.str() << "usmooth_rand_" << smoother + 1 << "_" << color + 1 << ".txt";
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
                        u_sc_M.push_back(p);
                    }
                }
                f.close();
                fname.clear();
                fname.str(std::string());
                std::vector<double> u_sc = v_level[l]->solve_gaussian_elimination(
                    v_level[l]->A_Zebra_r_LU[smoother_ind], v_level[l]->A_Zebra_c_LU[smoother_ind],
                    v_level[l]->A_Zebra_v_LU[smoother_ind], f_sc);
                for (int i = 0; i < m_sc; i++) {
                    assert(u_sc[i] - u_sc_M[i] < 1e-3);
                }
            }
        }

        std::cout << "Checking Smoother Complement Matrices Asc_ortho...\n";

        // - row_indices_Asc_ortho_s_c.txt, col_indices_Asc_ortho_s_c.txt, vals_Asc_ortho_s_c.txt
        for (int smoother = 0; smoother < 2; smoother++) {
            for (int color = 0; color < 2; color++) {
                std::vector<int> ri, ci;
                std::vector<double> v;

                fname << folder.str() << "row_indices_Asc_Mix_" << smoother + 1 << "_" << color + 1 << ".txt";
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
                        ri.push_back(p - 1);
                    }
                }
                f.close();
                fname.clear();
                fname.str(std::string());

                fname << folder.str() << "col_indices_Asc_Mix_" << smoother + 1 << "_" << color + 1 << ".txt";
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
                        ci.push_back(p - 1);
                    }
                }
                f.close();
                fname.clear();
                fname.str(std::string());

                fname << folder.str() << "vals_Asc_Mix_" << smoother + 1 << "_" << color + 1 << ".txt";
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
                        v.push_back(p);
                    }
                }
                f.close();
                fname.clear();
                fname.str(std::string());
                // gyro::disp(ri, "ri");
                // gyro::disp(ci, "ci");
                // gyro::disp(v, "v");

                int smoother_ind = smoother * 2 + color;
                int m_sc         = v_level[l]->m_sc[smoother_ind];

                cudaMalloc((void**)&dev_usc, m_sc * sizeof(double));

                for (int i = 0; i < m; i++) {
                    std::vector<double> ei(m, 0), uC(m_sc, 0), uM(m_sc, 0);
                    ei[i]         = 1;
                    v_level[l]->u = ei;

                    cudaMemcpy(dev_usc, &(uC[0]), m_sc * sizeof(double), cudaMemcpyHostToDevice);
                    cudaMemcpy(v_level[l]->dev_u, &(v_level[l]->u[0]), m * sizeof(double), cudaMemcpyHostToDevice);
                    apply_Asc_ortho_gpu<<<n_blocks, NTHREADS>>>(
                        m, dev_usc, v_level[l]->dev_u, smoother_ind, v_level[l]->nr_int, v_level[l]->dev_r,
                        v_level[l]->dev_hplus, v_level[l]->dev_coarse_nodes_list_r, v_level[l]->ntheta_int,
                        v_level[l]->dev_theta, v_level[l]->dev_thetaplus, v_level[l]->dev_coarse_nodes_list_theta,
                        v_level[l]->delete_circles, gyro::icntl[Param::DirBC_Interior], gyro::icntl[Param::mod_pk],
                        gyro::dcntl[Param::r0_DB], gyro::icntl[Param::prob], gyro::dcntl[Param::kappa_eps],
                        gyro::dcntl[Param::delta_e], gyro::dcntl[Param::R]);

                    cudaMemcpy(&uC[0], dev_usc, m_sc * sizeof(double), cudaMemcpyDeviceToHost);

                    for (int j = 0; j < ri.size(); j++) {
                        uM[ri[j]] += v[j] * ei[ci[j]];
                    }
                    // gyro::disp(uC, "uC");
                    // gyro::disp(uM, "uM");
                    for (int j = 0; j < m_sc; j++) {
                        // printf("j: %d, uC[j]: %f, uM[j]: %f\n", j, uC[j], uM[j]);
                        if (uC[j] + uM[j] > 1e-3)
                            std::cout << "j: " << j << ", uC[j]: " << uC[j] << ", uM[j]: " << uM[j] << "\n";
                        assert(uC[j] + uM[j] < 1e-3);
                    }
                }

                cudaFree(dev_usc);
            }
        }

        folder.clear(); //clear any bits set
        folder.str(std::string());

        cudaFree(dev_Pu);
        cudaFree(dev_Ru);
        cudaFree(dev_uc);
        cudaFree(dev_uc2);
    }

    ///////////////////////////////////////////////////////////
    // Checking whole multigrid cycle
    ///////////////////////////////////////////////////////////
    std::cout << "Checking parameters...\n";
    std::vector<double> u_cycle;
    fname << ss.str() << "1/MG_cycle_rand.txt";
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
            u_cycle.push_back(p);
        }
    }
    f.close();
    fname.clear(); //clear any bits set
    fname.str(std::string());
    v_level[0]->u = u_rand0;
    multigrid_cycle_extrapol(0);
    for (int i = 0; i < v_level[0]->m; i++) {
        assert(v_level[0]->u[i] - u_cycle[i] < 1e-3);
    }

    std::vector<double> error_M;
    fname << ss.str() << "1/error_rand.txt";
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
            error_M.push_back(p);
        }
    }
    f.close();
    fname.clear(); //clear any bits set
    fname.str(std::string());
    std::vector<double> error = compute_error();
    for (int i = 0; i < v_level[0]->m; i++) {
        assert(error[i] - error_M[i] < 1e-3);
    }
} /* ----- end of gmgpolar::debug ----- */

#endif
