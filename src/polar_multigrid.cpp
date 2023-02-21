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
 * \file polar_multigrid.cpp
 * \brief Implementation of the whole gmgpolar (grid + multigrid)
 * \author M. Kuehn, C. Kruse, P. Leleux, C. Schwarz
 * \version 0.0
 */
#include "gmgpolar.h"
#include <unistd.h>
#include <random>
#include <algorithm>
#include <functional>

template <typename Generator>
void getrands(std::vector<double>& x, Generator& gen, unsigned num)
{
    generate_n(std::back_inserter(x), num, ref(gen));
}

/*!
 *  \brief Solves the problem gyro with multigrid on levels
 *
 *  Solves the problem gyro with multigrid on levels
 *
 */
void gmgpolar::polar_multigrid()
{
    if (gyro::icntl[Param::verbose] > 1)
        std::cout << "Define the coarse grids...\n";
    check_geom();
    define_coarse_nodes();

    if (gyro::icntl[Param::verbose] > 2) {
        std::cout << "\n";
        if (gyro::icntl[Param::prob] == 1)
            std::cout << "Solving the Poisson equation ";
        else if (gyro::icntl[Param::prob] == 5)
            std::cout << "Solving the Cartesian problem (from Zoni2019) ";
        else if (gyro::icntl[Param::prob] == 6)
            std::cout << "Solving the Poloidal problem ";
        if (gyro::icntl[Param::alpha_coeff] == 0)
            std::cout << "with coefficient alpha in atan (Sonnendrucker) ";
        else if (gyro::icntl[Param::alpha_coeff] == 1)
            std::cout << "with coefficient alpha in exp(tanh) and steep variation at 0.5 ";
        else if (gyro::icntl[Param::alpha_coeff] == 2)
            std::cout << "with coefficient alpha in exp(tanh) (Zoni2019) and variation at 0.7 ";
        if (gyro::icntl[Param::beta_coeff] == 0)
            std::cout << "and coefficient beta=0\n";
        else if (gyro::icntl[Param::beta_coeff] == 1)
            std::cout << "and coefficient beta=1/alpha\n";
        if (!gyro::icntl[Param::mod_pk])
            std::cout << "Considering POLAR coordinates with ";
        else if (gyro::icntl[Param::mod_pk] == 1)
            std::cout << "Considering the Shafranov geometry with kappa=" << gyro::dcntl[Param::kappa_eps]
                      << ", delta=" << gyro::dcntl[Param::delta_e] << ", and ";
        else if (gyro::icntl[Param::mod_pk] == 2)
            std::cout << "Considering the Czarny geometry with kappa=" << gyro::dcntl[Param::kappa_eps]
                      << ", delta=" << gyro::dcntl[Param::delta_e] << ", and ";
        std::cout << gyro::dcntl[Param::R0] << " <= r <= " << gyro::dcntl[Param::R] << "\n";
        std::cout << "Using 9 point star FINITE DIFFERENCES.\n";
        std::cout << "Using FULL EXTRAPOLATION\n\n";
    }

    gyro::dcntl[Param::r0_DB] = -1e6;
    if (gyro::icntl[Param::DirBC_Interior])
        gyro::dcntl[Param::r0_DB] = gyro::dcntl[Param::R0];

    int testA = 0;
    if (testA) {
        int testPerf        = 0;
        int display_vectors = 0;
        int nb_tests        = 10;
        int l               = 0;
        v_level[l]->m       = v_level[l]->nr * v_level[l]->ntheta;
        v_level[l]->define_nz();
        int m = v_level[l]->m;
        std::cout << "\n***** Problem size " << m << " (" << v_level[0]->nr << ", " << v_level[0]->ntheta << ")\n";

        double t, t_applyA = 0, t_applyA0 = 0;

        for (int i = 0; i < nb_tests; i++) {
            std::cout << "Test " << i << "\n";
            std::vector<double> ei, uC(m, 0);

            std::uniform_real_distribution<double> unif(0.0, 1.0);
            std::mt19937 re(std::random_device{}());
            auto generator = std::bind(unif, std::ref(re));
            getrands(ei, generator, m);
            if (display_vectors)
                gyro::disp(ei, "ei");

            // Apply A0
            if (!testPerf) {
                TIC;
                uC = std::vector<double>(m, 0);
                v_level[l]->apply_A0(ei, uC);
                t_applyA0 += TOC;
                if (display_vectors)
                    gyro::disp(uC, "uC (A0)");
            }
            // Apply A
            TIC;
            uC = std::vector<double>(m, 0);
            v_level[l]->apply_A(ei, uC);
            t_applyA += TOC;
            if (display_vectors)
                gyro::disp(uC, "uC (A)");
        }
        if (!testPerf)
            std::cout << "t_applyA0: " << t_applyA0 << std::endl;
        std::cout << "t_applyA: " << t_applyA << std::endl;
    }
    else {
        if (gyro::icntl[Param::verbose] > 1)
            std::cout
                << "Building discretized system, restriction and interpolation operators, and defining splittings...\n";
        prepare_op_levels();

        int m = v_level[0]->m;

        std::string cycle_str = "?";
        if (gyro::icntl[Param::cycle] == 1)
            cycle_str = "V";
        else if (gyro::icntl[Param::cycle] == 2)
            cycle_str = "W";
        // if (gyro::icntl[Param::verbose] > 1)
        std::cout << "\nProb: " << gyro::icntl[Param::prob] << ", alpha_coeff: " << gyro::icntl[Param::alpha_coeff]
                  << ", beta_coeff: " << gyro::icntl[Param::beta_coeff] << " ***** Problem size " << m << " ("
                  << v_level[0]->nr << ", " << v_level[0]->ntheta << "), " << levels
                  << " grids, nr_exp=" << gyro::icntl[Param::nr_exp] << ", aniso=" << gyro::icntl[Param::fac_ani]
                  << " ***** smoother=" << gyro::icntl[Param::smoother] << ", r0=" << gyro::dcntl[Param::R0]
                  << ", extrapolation=" << gyro::icntl[Param::extrapolation]
                  << ", mod_pk=" << gyro::icntl[Param::mod_pk] << ", DirBC=" << gyro::icntl[Param::DirBC_Interior]
                  << ", divide=" << gyro::icntl[Param::divideBy2] << " *****\n";

        double scaling = 1.0;
        if (gyro::icntl[Param::compute_rho])
            scaling = scaling / (sqrt(m));

        v_level[0]->u.assign(m, scaling); //create an empty vector u

        if (gyro::icntl[Param::debug] > 0) {
            debug();
        }
        else if (gyro::icntl[Param::write_radii_angles] == 0) {
            if (gyro::icntl[Param::verbose] > 1)
                std::cout << "Multigrid iteration....!\n\n";
            multigrid_iter();

            //        if (gyro::icntl[Param::verbose] > 1)
            for (int l = 0; l < levels; l++) {
                std::cout << "LEVEL " << l << "\n";
                std::cout << "\nt_smoothing: " << v_level[l]->t_smoothing << ", t_f_sc: " << v_level[l]->t_f_sc
                          << ", t_Asc_ortho: " << v_level[l]->t_Asc_ortho << ", t_Asc: " << v_level[l]->t_Asc << "\n";
                std::cout << "\nt_get_ptr: " << v_level[l]->t_get_ptr
                          << ", t_get_stencil: " << v_level[l]->t_get_stencil
                          << ", t_get_smoother: " << v_level[l]->t_get_smoother
                          << ", t_get_row: " << v_level[l]->t_get_row << "\n";
                std::cout << "\n";
            }

            //        if (gyro::icntl[Param::verbose] > 0) {
            std::cout << "\nt_setup: " << t_setup << ", t_build: " << t_build << ", t_facto_Ac: " << t_facto_Ac
                      << ", t_build_P: " << t_build_P << ", t_build_Asc: " << t_build_Asc
                      << ", t_facto_Asc: " << t_facto_Asc << "\n";
            std::cout << "t_total_(fine): " << t_total << ", t_smoothing: " << t_smoothing
                      << ", t_residual: " << t_residual << ", t_restriction: " << t_restriction << ", t_Ac: " << t_Ac
                      << ", t_prolongation: " << t_prolongation << ", t_fine_residual: " << t_fine_residual
                      << ", t_error: " << t_error << "\n";
            std::cout << "t_applyA: " << t_applyA << std::endl;
            //        }

            //        if (gyro::icntl[Param::verbose] > 0) {
            std::cout << "\nt_coeff: " << gyro::dcntl[Param::t_coeff]
                      << ", t_arr_art_att: " << gyro::dcntl[Param::t_arr_art_att]
                      << ", t_sol: " << gyro::dcntl[Param::t_sol] << ", t_detDFinv: " << gyro::dcntl[Param::t_detDFinv]
                      << ", t_trafo: " << gyro::dcntl[Param::t_trafo] << "\n";
            //        }
        }
    }
} /* ----- end of gmgpolar::polar_multigrid ----- */

/*!
 *  \brief Check the geometry construction
 *
 *  Check the geometry construction
 *
 */
void gmgpolar::check_geom()
{
    // intervals(!) not nodes in r direction (corresponds to nodes-1,
    // i.e., without considering the origin which is here always defined
    // as a node of the finest mesh only)
    int ntheta = v_level[0]->ntheta - 1;
    int nr     = v_level[0]->nr - 1;

    if (gyro::icntl[Param::level] == -1 || gyro::icntl[Param::level] < 2) {
        // at least two levels/grids.. and if more than two grids, at least 3 nodes
        // per circle and per row on coarsest mesh!!
        levels = std::max(2, (int)floor(std::min(log2(nr + 1) - 1, log2(ntheta + 2) - 1)));
    }
    else {
        levels = gyro::icntl[Param::level];
    }
    std::cout << "Desired number of levels: " << levels << "\n";

    double Rmax               = v_level[0]->r[nr];
    double R0                 = v_level[0]->r[0];
    gyro::dcntl[Param::r0_DB] = -1;
    if (R0 > 0)
        gyro::dcntl[Param::r0_DB] = R0;
    if (gyro::icntl[Param::verbose] > 2)
        std::cout << "levels: " << levels << ", Rmax: " << Rmax << "\n";
    if (Rmax != gyro::dcntl[Param::R])
        throw std::runtime_error("Program stopped... r[end] != R...");
    if (R0 != gyro::dcntl[Param::R0])
        throw std::runtime_error("Program stopped... r[0] != R0...");
    if (gyro::dcntl[Param::extrapolation] > 0 && levels < 3)
        throw std::runtime_error("ATTENTION: extrapolation technique needs at least three levels...");
} /* ----- end of gmgpolar::polar_multigrid ----- */

/*!
 *  \brief Prepare the levels for operator construction
 *
 *  Prepare the levels for operator construction
 *
 */
void gmgpolar::prepare_op_levels()
{
    double t;
    TIC;
    t_setup = t;

    for (int l = levels - 1; l >= 0; l--) { //define m on the coarsest level
        v_level[l]->m = v_level[l]->nr * v_level[l]->ntheta;
        v_level[l]->define_nz();

        v_level[l]->betaVec = std::vector<double>(v_level[l]->m, 0);
        v_level[l]->build_betaVec();

        if (l < levels - 1)
            v_level[l]->mc = v_level[l + 1]->m;
        if (gyro::icntl[Param::verbose] > 2)
            std::cout << "Create boundary array\n";
        v_level[l]->build_bound();

        if (l == 0 && !gyro::f_sol_in.empty()) {
            v_level[0]->read_sol();
        }

        if (gyro::icntl[Param::verbose] > 1)
            std::cout << "Create operator on level " << l << "\n";
        if (l == levels - 1 || gyro::icntl[Param::matrix_free] == 0) {
            TIC;
            if (gyro::icntl[Param::optimized] == 0 || v_level[l]->nr_int < 2 + gyro::icntl[Param::DirBC_Interior]) {
                std::cout << "Using the original construction of A and the RHS since nr_int is very "
                             "small.\n\n\n\n";
                v_level[l]->build_A0();
                v_level[l]->build_rhs0();
            }
            else {
                v_level[l]->row_indices = std::vector<int>(v_level[l]->nz);
                v_level[l]->col_indices = std::vector<int>(v_level[l]->nz);
                v_level[l]->vals        = std::vector<double>(v_level[l]->nz, 0);
                v_level[l]->fVec        = std::vector<double>(v_level[l]->m);
                v_level[l]->build_A();
                if (l == 0 || (l == 1 && gyro::icntl[Param::extrapolation] > 0)) {
                    v_level[l]->build_rhs();
                }
            }
            t_build += TOC;
            TIC;

            if (l == levels - 1) {
                if (gyro::icntl[Param::verbose] > 1)
                    std::cout << "Factorizing coarse operator...\n";
                TIC;
#ifdef USE_MUMPS
                if (gyro::icntl[Param::optimized] == 0) {
#endif
                    v_level[l]->row_Ac_LU  = std::vector<int>(v_level[l]->row_indices);
                    v_level[l]->col_Ac_LU  = std::vector<int>(v_level[l]->col_indices);
                    v_level[l]->vals_Ac_LU = std::vector<double>(v_level[l]->vals);
                    v_level[l]->row_indices.clear();
                    v_level[l]->col_indices.clear();
                    v_level[l]->vals.clear();
                    v_level[l]->row_indices.shrink_to_fit();
                    v_level[l]->col_indices.shrink_to_fit();
                    v_level[l]->vals.shrink_to_fit();
                    v_level[l]->facto_gaussian_elimination(v_level[l]->row_Ac_LU, v_level[l]->col_Ac_LU,
                                                           v_level[l]->vals_Ac_LU, v_level[l]->m);
#ifdef USE_MUMPS
                }
                else
                    v_level[l]->facto_mumps(v_level[l]->mumps_Ac, v_level[l]->row_indices, v_level[l]->col_indices,
                                            v_level[l]->vals, v_level[l]->m);
#endif
                t_facto_Ac += TOC;
            }
            TIC;
        }
        else {
            TIC;
            if (gyro::icntl[Param::optimized] == 0)
                v_level[l]->build_rhs0();
            else {
                v_level[l]->fVec = std::vector<double>(v_level[l]->m);
                if (l == 0 || (l == 1 && gyro::icntl[Param::extrapolation] > 0)) {
                    v_level[l]->build_rhs();
                }
            }
            t_build += TOC;
            TIC;
        }

        // Prolongation defined on all except the coarsest level
        if (l < levels - 1) {
            v_level[l]->define_line_splitting();
            if (gyro::icntl[Param::verbose] > 1)
                std::cout << "delete_circles: " << v_level[l]->delete_circles << "\n";

            // Number of blocks per smoother
            v_level[l]->nblocks    = std::vector<int>(5);
            v_level[l]->nblocks[0] = ceil(v_level[l]->delete_circles * 0.5);
            v_level[l]->nblocks[1] = floor(v_level[l]->delete_circles * 0.5);
            v_level[l]->nblocks[2] = v_level[l]->ntheta_int * 0.5;
            v_level[l]->nblocks[3] = v_level[l]->nblocks[2];
            v_level[l]->nblocks[4] = std::max(v_level[l]->nblocks[0], v_level[l]->nblocks[2]);

            for (int smoother = 0; smoother < 4; smoother++) {
                int size = 0, size_ortho = 0, *array_temp, *array_temp2;
                if (smoother < 2) {
                    size_ortho = v_level[l]->delete_circles + 1;
                    size       = v_level[l]->delete_circles;
                }
                else if (smoother > 1) {
                    size_ortho = v_level[l]->ntheta_int;
                    size       = v_level[l]->ntheta_int;
                }
                array_temp = new int[size_ortho];
                for (int i = 0; i < size_ortho; i++)
                    array_temp[i] = 0;
                array_temp2 = new int[size];
                for (int i = 0; i < size; i++)
                    array_temp2[i] = 0;
                v_level[l]->dep_Asc_ortho.push_back(array_temp);
                v_level[l]->dep_Asc.push_back(array_temp2);
                v_level[l]->size_Asc_ortho.push_back(size_ortho);
                v_level[l]->size_Asc.push_back(size);
            }

            // Smoother matrices Asc
            TIC;
            //build matrices A_sc
            if (gyro::icntl[Param::verbose] > 1)
                std::cout << "build Asc on level " << l << ", m: " << v_level[l]->m << ", mc: " << v_level[l]->mc
                          << "\n";
            v_level[l]->define_m_nz_Asc();
            if (gyro::icntl[Param::optimized] == 0) {
                v_level[l]->build_Asc0();
            }
            else {
                // 1 block matrix per row (column) for the circle (radial) smoother
                for (int smoother = 0; smoother < 4; smoother++) {
                    v_level[l]->A_Zebra_r_row[smoother].assign(v_level[l]->nblocks[smoother], std::vector<int>());
                    v_level[l]->A_Zebra_c_row[smoother].assign(v_level[l]->nblocks[smoother], std::vector<int>());
                    v_level[l]->A_Zebra_v_row[smoother].assign(v_level[l]->nblocks[smoother], std::vector<double>());
                    for (int ij = 0; ij < v_level[l]->nblocks[smoother]; ij++) {
                        int nsc                                 = v_level[l]->define_nz_Asc_ij(smoother, ij, 0);
                        v_level[l]->A_Zebra_r_row[smoother][ij] = std::vector<int>(nsc);
                        v_level[l]->A_Zebra_c_row[smoother][ij] = std::vector<int>(nsc);
                        v_level[l]->A_Zebra_v_row[smoother][ij] = std::vector<double>(nsc, 0);
                    }
                }
                // define Asc blocks
                v_level[l]->build_Asc();

                // 1 block matrix per row (column) for the circle (radial) smoother
                v_level[l]->A_Zebra_Mix_r.assign(4, std::vector<int>());
                v_level[l]->A_Zebra_Mix_c.assign(4, std::vector<int>());
                v_level[l]->A_Zebra_Mix_v.assign(4, std::vector<double>());
                for (int smoother = 0; smoother < 4; smoother++) {
                    int nsc                             = v_level[l]->nz_sc_ortho[smoother];
                    v_level[l]->A_Zebra_Mix_r[smoother] = std::vector<int>(nsc);
                    v_level[l]->A_Zebra_Mix_c[smoother] = std::vector<int>(nsc);
                    v_level[l]->A_Zebra_Mix_v[smoother] = std::vector<double>(nsc, 0);

                    // define Asc_ortho block
                    v_level[l]->build_Asc_ortho(smoother);

                    // Build vectors necessary for the parallel application of Asc_ortho:
                    // - ptr contains the nz entry for the points in the first radial line
                    // - shift contains the number of entries per node
                    int size_radial_line = v_level[l]->nr_int - v_level[l]->delete_circles;
                    v_level[l]->shift_vect_s2 =
                        std::vector<int>(size_radial_line); // shift between 2 radial lines for smoother 2
                    v_level[l]->shift_vect_s3 = std::vector<int>(size_radial_line); // idem for smoother 3
                    v_level[l]->ptr_vect_s2 = std::vector<int>(size_radial_line); // ptr to a radial line for smoother 2
                    v_level[l]->ptr_vect_s3 = std::vector<int>(size_radial_line); // idem for smoother 3
                    std::vector<int> ptr_vect;
                    for (int j = v_level[l]->delete_circles; j < v_level[l]->nr_int; j++) {
                        ptr_vect                                                  = v_level[l]->get_ptr_sc(j, 2, 1);
                        v_level[l]->ptr_vect_s2[j - v_level[l]->delete_circles]   = ptr_vect[0];
                        v_level[l]->ptr_vect_s3[j - v_level[l]->delete_circles]   = ptr_vect[1];
                        v_level[l]->shift_vect_s2[j - v_level[l]->delete_circles] = ptr_vect[2] - ptr_vect[0];
                        v_level[l]->shift_vect_s3[j - v_level[l]->delete_circles] = ptr_vect[3] - ptr_vect[1];
                    }
                }
            }

            t_build_Asc += TOC;
            TIC;

            if (gyro::icntl[Param::optimized] == 0) {
#ifdef USE_MUMPS
                if (gyro::icntl[Param::optimized] == 0) {
#endif
                    v_level[l]->A_Zebra_r_LU.assign(4, std::vector<int>());
                    v_level[l]->A_Zebra_c_LU.assign(4, std::vector<int>());
                    v_level[l]->A_Zebra_v_LU.assign(4, std::vector<double>());
                    for (int smoother = 0; smoother < 4; smoother++) {
                        v_level[l]->A_Zebra_r_LU[smoother] = std::vector<int>(v_level[l]->A_Zebra_r[smoother]);
                        v_level[l]->A_Zebra_c_LU[smoother] = std::vector<int>(v_level[l]->A_Zebra_c[smoother]);
                        v_level[l]->A_Zebra_v_LU[smoother] = std::vector<double>(v_level[l]->A_Zebra_v[smoother]);
                        v_level[l]->A_Zebra_r[smoother].clear();
                        v_level[l]->A_Zebra_c[smoother].clear();
                        v_level[l]->A_Zebra_v[smoother].clear();
                        v_level[l]->A_Zebra_r[smoother].shrink_to_fit();
                        v_level[l]->A_Zebra_c[smoother].shrink_to_fit();
                        v_level[l]->A_Zebra_v[smoother].shrink_to_fit();
                        if (gyro::icntl[Param::verbose] > 1)
                            std::cout << "Factorizing smoother " << smoother << "...\n";
                        v_level[l]->facto_gaussian_elimination(
                            v_level[l]->A_Zebra_r_LU[smoother], v_level[l]->A_Zebra_c_LU[smoother],
                            v_level[l]->A_Zebra_v_LU[smoother], v_level[l]->m_sc[smoother]);
                    }
                    v_level[l]->A_Zebra_r.clear();
                    v_level[l]->A_Zebra_c.clear();
                    v_level[l]->A_Zebra_v.clear();
                    v_level[l]->A_Zebra_r.shrink_to_fit();
                    v_level[l]->A_Zebra_c.shrink_to_fit();
                    v_level[l]->A_Zebra_v.shrink_to_fit();
#ifdef USE_MUMPS
                }
                else
                    for (int smoother = 0; smoother < 4; smoother++)
                        v_level[l]->facto_mumps(v_level[l]->mumps_A_Zebra[smoother], v_level[l]->A_Zebra_r[smoother],
                                                v_level[l]->A_Zebra_c[smoother], v_level[l]->A_Zebra_v[smoother],
                                                v_level[l]->m_sc[smoother]);
#endif
                t_facto_Asc += TOC;
                TIC;
            }
            else {
                for (int smoother = 0; smoother < 4; smoother++) {
                    TIC;
                    v_level[l]->A_Zebra_r_LU_row[smoother].assign(v_level[l]->nblocks[smoother], std::vector<int>());
                    v_level[l]->A_Zebra_c_LU_row[smoother].assign(v_level[l]->nblocks[smoother], std::vector<int>());
                    v_level[l]->A_Zebra_v_LU_row[smoother].assign(v_level[l]->nblocks[smoother], std::vector<double>());
                    t_facto_Asc += TOC;
                    TIC;
                    for (int ij = 0; ij < v_level[l]->nblocks[smoother]; ij++) {
                        // Diagonal smoother if:
                        // - DB and first radius
                        // - only 1 element (should not happen)
                        // - extrapolation on level 0 for smoothers 0 and 2 (except across the origin stencil)
                        if ((smoother == 0 && ij == 0 && gyro::icntl[Param::DirBC_Interior]) ||
                            v_level[l]->m_sc[smoother] == 1 ||
                            (gyro::icntl[Param::extrapolation] && l == 0 && smoother % 2 == 0 &&
                             !(smoother == 0 && !gyro::icntl[Param::DirBC_Interior] && ij == 0))) {
                            if (gyro::icntl[Param::verbose] > 2)
                                std::cout << "No diagonal factorization\n";
                            v_level[l]->A_Zebra_v_LU_row[smoother][ij] =
                                std::vector<double>(v_level[l]->A_Zebra_v_row[smoother][ij]);
                        }
                        // Circle smoother if:
                        // - smoother 0 or 1
                        // - and not across
                        else if (smoother < 2 && !(smoother == 0 && !gyro::icntl[Param::DirBC_Interior] && ij == 0)) {
                            TIC;
                            if (gyro::icntl[Param::verbose] > 2)
                                std::cout << "Circle factorization\n";
                            // Initialize the LU factors
                            v_level[l]->fill_in_circle(ij, smoother);
                            // Factorization
                            v_level[l]->facto_circle(
                                v_level[l]->A_Zebra_r_LU_row[smoother][ij], v_level[l]->A_Zebra_c_LU_row[smoother][ij],
                                v_level[l]->A_Zebra_v_LU_row[smoother][ij], v_level[l]->m_sc[smoother]);
#pragma omp atomic
                            t_facto_Asc += TOC;
                            TIC;
                        }
                        else if (smoother > 1) {
                            TIC;
                            // Initialize the LU factors
                            v_level[l]->A_Zebra_r_LU_row[smoother][ij] =
                                std::vector<int>(v_level[l]->A_Zebra_r_row[smoother][ij]);
                            v_level[l]->A_Zebra_c_LU_row[smoother][ij] =
                                std::vector<int>(v_level[l]->A_Zebra_c_row[smoother][ij]);
                            v_level[l]->A_Zebra_v_LU_row[smoother][ij] =
                                std::vector<double>(v_level[l]->A_Zebra_v_row[smoother][ij]);
                            if (gyro::icntl[Param::verbose] > 2)
                                std::cout << "Radial factorization\n";
                            // Factorization
                            v_level[l]->facto_radial(
                                v_level[l]->A_Zebra_r_LU_row[smoother][ij], v_level[l]->A_Zebra_c_LU_row[smoother][ij],
                                v_level[l]->A_Zebra_v_LU_row[smoother][ij], v_level[l]->m_sc[smoother]);
#pragma omp atomic
                            t_facto_Asc += TOC;
                            TIC;
                        }
                        if (smoother == 0 && !gyro::icntl[Param::DirBC_Interior] && ij == 0) {
                            TIC;
                            v_level[l]->A_Zebra_r_LU_row[smoother][ij] =
                                std::vector<int>(v_level[l]->A_Zebra_r_row[smoother][ij]);
                            v_level[l]->A_Zebra_c_LU_row[smoother][ij] =
                                std::vector<int>(v_level[l]->A_Zebra_c_row[smoother][ij]);
                            v_level[l]->A_Zebra_v_LU_row[smoother][ij] =
                                std::vector<double>(v_level[l]->A_Zebra_v_row[smoother][ij]);
                            if (gyro::icntl[Param::verbose] > 2)
                                std::cout << "Across factorization...";
#ifdef USE_MUMPS
                            if (gyro::icntl[Param::optimized] == 0) {
#endif
                                if (gyro::icntl[Param::verbose] > 2)
                                    std::cout << "using in-house direct solver.\n";
                                // Factorization
                                v_level[l]->facto_gaussian_elimination(v_level[l]->A_Zebra_r_LU_row[smoother][ij],
                                                                       v_level[l]->A_Zebra_c_LU_row[smoother][ij],
                                                                       v_level[l]->A_Zebra_v_LU_row[smoother][ij],
                                                                       v_level[l]->m_sc[smoother]);
#ifdef USE_MUMPS
                            }
                            else {
                                if (gyro::icntl[Param::verbose] > 2)
                                    std::cout << "using MUMPS.\n";
                                v_level[l]->facto_mumps(
                                    v_level[l]->mumps_across, v_level[l]->A_Zebra_r_row[smoother][ij],
                                    v_level[l]->A_Zebra_c_row[smoother][ij], v_level[l]->A_Zebra_v_row[smoother][ij],
                                    v_level[l]->m_sc[smoother]);
                            }
#endif
#pragma omp atomic
                            t_facto_Asc += TOC;
                            TIC;
                        }
                        v_level[l]->A_Zebra_r_row[smoother][ij].clear();
                        v_level[l]->A_Zebra_c_row[smoother][ij].clear();
                        v_level[l]->A_Zebra_v_row[smoother][ij].clear();
                        v_level[l]->A_Zebra_r_row[smoother][ij].shrink_to_fit();
                        v_level[l]->A_Zebra_c_row[smoother][ij].shrink_to_fit();
                        v_level[l]->A_Zebra_v_row[smoother][ij].shrink_to_fit();
                    }
                    v_level[l]->A_Zebra_r_row[smoother].clear();
                    v_level[l]->A_Zebra_c_row[smoother].clear();
                    v_level[l]->A_Zebra_v_row[smoother].clear();
                    v_level[l]->A_Zebra_r_row[smoother].shrink_to_fit();
                    v_level[l]->A_Zebra_c_row[smoother].shrink_to_fit();
                    v_level[l]->A_Zebra_v_row[smoother].shrink_to_fit();
                }
            }

            // Prolongation defined on all except the coarsest level
            TIC;
            if (gyro::icntl[Param::matrix_free] == 0) {
                v_level[l]->define_nz_P();

                // // Number of nonzeros in the matrix
                if (gyro::icntl[Param::verbose] > 2)
                    std::cout << "nz_P: " << v_level[l]->nz_P << "\n";
                v_level[l]->ri_prol = std::vector<int>(v_level[l]->nz_P);
                v_level[l]->ci_prol = std::vector<int>(v_level[l]->nz_P);
                v_level[l]->v_prol  = std::vector<double>(v_level[l]->nz_P);
                v_level[l]->build_prolongation_bi();

                // // Number of nonzeros in the matrix
                if (gyro::icntl[Param::verbose] > 2)
                    std::cout << "nz_P_inj: " << v_level[l]->nz_P_inj << "\n";
                v_level[l]->ri_prol_inj = std::vector<int>(v_level[l]->nz_P_inj);
                v_level[l]->ci_prol_inj = std::vector<int>(v_level[l]->nz_P_inj);
                v_level[l]->v_prol_inj  = std::vector<double>(v_level[l]->nz_P_inj);
                v_level[l]->build_prolongation_inj();

                // // Number of nonzeros in the matrix
                if (gyro::icntl[Param::verbose] > 2)
                    std::cout << "nz_P_ex: " << v_level[l]->nz_P_ex << "\n";
                v_level[l]->ri_prol_ex = std::vector<int>(v_level[l]->nz_P_ex);
                v_level[l]->ci_prol_ex = std::vector<int>(v_level[l]->nz_P_ex);
                v_level[l]->v_prol_ex  = std::vector<double>(v_level[l]->nz_P_ex);
                v_level[l]->build_prolongation_ex();

                t_build_P += TOC;
                TIC;
            }
        }
    }

    if (gyro::icntl[Param::verbose] > 2)
        std::cout << "Operators created.\n";

    t       = t_setup;
    t_setup = TOC;
} /* ----- end of level::prepare_op_levels ----- */