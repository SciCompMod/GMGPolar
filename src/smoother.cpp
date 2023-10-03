/* 
* Copyright (C) 2019-2023 The GMGPolar Development Team
*
* Authors: Philippe Leleux, Christina Schwarz, Martin J. Kühn, Carola Kruse, Ulrich Rüde, Ulrich Rüde
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
 * \file multigrid_iter.cpp
 * \brief Implementation of the smoother
 * \author M. Kuehn, C. Kruse, P. Leleux, C. Schwarz
 * \version 0.0
 */

#include "level.h"
#include <sstream>

/*!
 *  \brief Returns the index i+1 in (0, ntheta)
 *
 * Returns (i + 1) % ntheta_int
 * 
 * \param i: theta index
 * \param ntheta_int: size of theta
 * 
 * \return (i + 1) % ntheta_int
 *
 */
inline int iplus(int i, int ntheta_int)
{
    return (i + 1 > ntheta_int - 1) ? 0 : i + 1;
} /* ----- end of inline iplus ----- */

/*!
 *  \brief Returns the index i-1 in (0, ntheta)
 *
 * Returns (i + ntheta_int - 1) % ntheta_int
 * 
 * \param i: theta index
 * \param ntheta_int: size of theta
 * 
 * \return (i + ntheta_int - 1) % ntheta_int
 *
 */
inline int imoins(int i, int ntheta_int)
{
    return (i - 1 < 0) ? ntheta_int - 1 : i - 1;
} /* ----- end of inline imoins ----- */

/*!
 *  \brief Defines the smoother
 *
 * Defines the smoother
 * so far only different types of zebra splittings are implemented
 * (basically, two colors, although more can be defined but which might not be used)
 *
 */
void level::define_line_splitting()
{
    // useless array q
    int i;
    for (i = 2; i < nr_int - 2; i++) {
        // assume that k_ij=theta_{i,j+1}-theta_{i,j} is constant for fixed i and variable j (i.e., constant per circle)
        double q = (theta[1] - theta[0]) / (r[i] - r[i - 1]);
        if (r[i] > 1 / q)
            break;
    }
    delete_circles = i; // delete_circles = delete_circles(end);

    if (gyro::icntl[Param::verbose] > 2)
        std::cout << "Shifting from circle to radial at radius " << delete_circles << "\n";
} /* ----- end of destructor level::define_line_splitting ----- */

/*! \brief Applies smoothing
 *
 * For all lines of the smoother s and the colour c.
 * The matrix A_sc is for all lines of points of the smoother s/c.
 * The result is in u (or u_previous_cu_previous_r).
 * 
 * if gyro::icntl[Param::smoother]=
 * - 3: block Gauss Seidel for all 4 smoothers.
 * - 13: block Jacobi between Circle and Radial, and block Gauss Seidel for B/W inside each.
 * 
 * \param smoother: the smoother
*/
void level::multigrid_smoothing(int smoother, int v, std::vector<double>& f_Asc_u, int nblocks, int c, int* dep_Asc_cur,
                                int* dep_Asc_prev, int* dep_Asc1, int* dep_Asc_ortho_cur)
{
    double t, t_smoothing_tmp;
    int extrapol = gyro::icntl[Param::extrapolation] == 1 && l == 0;

    TIC;
    t_smoothing_tmp = t;

    if (gyro::icntl[Param::matrix_free] == 1) {
        apply_Asc_ortho(f_Asc_u, u, smoother, v, c, dep_Asc_cur, dep_Asc_prev, dep_Asc1, dep_Asc_ortho_cur);
    }
    else {
        /*        if (gyro::icntl[Param::openmp] == 1) {
            for (auto k = 0; k < A_Zebra_Mix_r[smoother].size(); k++) {
                f_Asc_u[A_Zebra_Mix_r[smoother][k]] -= A_Zebra_Mix_v[smoother][k] * u[A_Zebra_Mix_c[smoother][k]];
            }
        }
        else {
*/
        if (smoother < 2) {
            int start_j = 0;
            if (gyro::icntl[Param::DirBC_Interior])
                start_j = 1;

            int start;
            if (smoother == 0)
                start = (gyro::icntl[Param::DirBC_Interior]) ? 2 : 0;
            else if (smoother == 1)
                start = 1;
            int shift = 2;
            for (int j = start; j < delete_circles; j += shift) {
                int odd_j = j % 2;
#pragma omp task shared(u, f_Asc_u) firstprivate(smoother, j, odd_j) depend(in                                         \
                                                                            : dep_Asc_prev[j - odd_j]),                \
    depend(in                                                                                                          \
           : dep_Asc_prev[j + odd_j]),                                                                                 \
    depend(out                                                                                                         \
           : dep_Asc_ortho_cur[j])
                {
                    int smoother_tmp          = get_smoother(0, j);
                    std::vector<int> ptr_vect = get_ptr_sc(j, smoother_tmp, 1);
                    int ptr_start             = ptr_vect[0];
                    int ptr_end               = ptr_start + (ptr_vect[2] - ptr_vect[0]) * ntheta_int / 2;
                    for (int k = ptr_start; k < ptr_end; k++) {
                        f_Asc_u[A_Zebra_Mix_r[smoother][k]] -=
                            A_Zebra_Mix_v[smoother][k] * u[A_Zebra_Mix_c[smoother][k]];
                    }
                    dep_Asc_ortho_cur[j] = 1;
                }
            }
        }
        else {
            std::vector<int> shift_vect, ptr_vect;
            if (smoother == 2) {
                ptr_vect   = ptr_vect_s2;
                shift_vect = shift_vect_s2;
            }
            else {
                ptr_vect   = ptr_vect_s3;
                shift_vect = shift_vect_s3;
            }
            for (int i = smoother - 2; i < ntheta_int; i += 2) {
                int im        = imoins(i, ntheta_int);
                int ip        = iplus(i, ntheta_int);
                int odd_d     = delete_circles % 2;
                int ind_m_odd = i, ind_p_odd = i;
                int odd_i = i % 2;
                if (odd_i) {
                    ind_m_odd = im;
                    ind_p_odd = ip;
                }
#pragma omp task shared(u, f_Asc_u) firstprivate(smoother, i, odd_d, ind_m_odd, ind_p_odd)                             \
    depend(in                                                                                                          \
           : dep_Asc_prev[ind_m_odd]) depend(in                                                                        \
                                             : dep_Asc_prev[ind_p_odd], dep_Asc1[delete_circles - 1 - odd_d])          \
        depend(out                                                                                                     \
               : dep_Asc_ortho_cur[i])
                {
                    int start_s2_extr = (smoother == 2 && extrapol && delete_circles % 2 == 0) ? 1 : 0;
                    int shift_s2_extr = (smoother == 2 && extrapol) ? 2 : 1;
                    for (int j = delete_circles + start_s2_extr; j < nr_int; j += shift_s2_extr) {
                        for (int k =
                                 ptr_vect[j - delete_circles] + (i - smoother + 2) / 2 * shift_vect[j - delete_circles];
                             k < ptr_vect[j - delete_circles] +
                                     ((i - smoother + 2) / 2 + 1) * shift_vect[j - delete_circles];
                             k++) {
                            f_Asc_u[A_Zebra_Mix_r[smoother][k]] -=
                                A_Zebra_Mix_v[smoother][k] * u[A_Zebra_Mix_c[smoother][k]];
                        }
                    }
                    dep_Asc_ortho_cur[i] = 1;
                }
            }
        }
        //        }
    }

    t_Asc_ortho += TOC;
    TIC;

    /*    if (gyro::icntl[Param::openmp] == 1) {
        for (int k = 0; k < nblocks; k++) {
            int ij_glob = k * 2 + smoother % 2;
            int ind = ij_glob, indm = ij_glob - 1, indp = ij_glob + 1;
            if (smoother > 1 && gyro::icntl[Param::matrix_free] == 1) {
                indm = imoins(ij_glob, ntheta_int);
                indp = iplus(ij_glob, ntheta_int);
            }
            TIC;
            // create u_sc, f_sc
            std::vector<double> u_sc;
            std::vector<double> f_sc(m_sc[smoother]);
            std::vector<double> f_total(m_sc[smoother]), f_total2(m_sc[smoother]);

            std::vector<double>::const_iterator first = f_Asc_u.begin() + k * m_sc[smoother];
            std::vector<double>::const_iterator last  = f_Asc_u.begin() + (k + 1) * m_sc[smoother];
            f_total                                   = std::vector<double>(first, last);
            for (int i = 0; i < m_sc[smoother]; i++)
                f_Asc_u[k * m_sc[smoother] + i] = 0;

            //compute f_total = f_sc - A_sc_ortho * u
            build_fsc(f_sc, fVec, smoother, 0, k * m_sc[smoother], (k + 1) * m_sc[smoother]);
            for (int i = 0; i < m_sc[smoother]; i++) {
                f_sc[i] += f_total[i];
            }

            // #pragma omp atomic
            t_f_sc += TOC;
            TIC;

            // Dirichlet or size of system is 1 (diagonal solve)
            if ((smoother == 0 && k == 0 && gyro::icntl[Param::DirBC_Interior]) || m_sc[smoother] == 1 ||
                (gyro::icntl[Param::extrapolation] && l == 0 && smoother % 2 == 0 &&
                 !(smoother == 0 && !gyro::icntl[Param::DirBC_Interior] && k == 0))) {
                u_sc = solve_diag(A_Zebra_v_LU_row[smoother][k], f_sc);
            }
            // Circle (not across)
            else if (smoother < 2 && !(smoother == 0 && !gyro::icntl[Param::DirBC_Interior] && k == 0)) {
                u_sc = solve_circle(A_Zebra_r_LU_row[smoother][k], A_Zebra_c_LU_row[smoother][k],
                                    A_Zebra_v_LU_row[smoother][k], f_sc);
            }
            // Radial (not across)
            else if (smoother > 1) {
                u_sc = solve_radial(A_Zebra_r_LU_row[smoother][k], A_Zebra_c_LU_row[smoother][k],
                                    A_Zebra_v_LU_row[smoother][k], f_sc);
            }
            // Across (direct solver)
            else {
#ifdef GMGPOLAR_USE_MUMPS
                if (gyro::icntl[Param::optimized] == 0) {
#endif
                    u_sc = solve_gaussian_elimination(A_Zebra_r_LU_row[smoother][k], A_Zebra_c_LU_row[smoother][k],
                                                      A_Zebra_v_LU_row[smoother][k], f_sc);
#ifdef GMGPOLAR_USE_MUMPS
                }
                else
                    u_sc = solve_mumps(mumps_across, f_sc);
#endif
            }

            t_Asc += TOC;
            TIC;

            build_fsc(u_sc, u, smoother, 1, k * m_sc[smoother], (k + 1) * m_sc[smoother]);

            t_f_sc += TOC;
            TIC;
        }
    }
    else {
*/
    for (int k = 0; k < nblocks; k++) {
        int ij_glob = k * 2 + smoother % 2;
        int ind = ij_glob, indm = ij_glob - 1, indp = ij_glob + 1;
        if (smoother > 1 && gyro::icntl[Param::matrix_free] == 1) {
            indm = imoins(ij_glob, ntheta_int);
            indp = iplus(ij_glob, ntheta_int);
        }
#pragma omp task firstprivate(k, ij_glob, ind, indm, indp, t) shared(f_Asc_u)                                          \
    depend(in                                                                                                          \
           : dep_Asc_ortho_cur[indm], dep_Asc_ortho_cur[ind], dep_Asc_ortho_cur[indp]) depend(out                      \
                                                                                              : dep_Asc_cur[ind])
        {
            TIC;
            // create u_sc, f_sc
            std::vector<double> u_sc;
            std::vector<double> f_sc(m_sc[smoother]);
            std::vector<double> f_total(m_sc[smoother]), f_total2(m_sc[smoother]);

            std::vector<double>::const_iterator first = f_Asc_u.begin() + k * m_sc[smoother];
            std::vector<double>::const_iterator last  = f_Asc_u.begin() + (k + 1) * m_sc[smoother];
            f_total                                   = std::vector<double>(first, last);
            for (int i = 0; i < m_sc[smoother]; i++)
                f_Asc_u[k * m_sc[smoother] + i] = 0;

            //compute f_total = f_sc - A_sc_ortho * u
            build_fsc(f_sc, fVec, smoother, 0, k * m_sc[smoother], (k + 1) * m_sc[smoother]);
            for (int i = 0; i < m_sc[smoother]; i++) {
                f_sc[i] += f_total[i];
            }

            // #pragma omp atomic
            t_f_sc += TOC;
            TIC;

            // Dirichlet or size of system is 1 (diagonal solve)
            if ((smoother == 0 && k == 0 && gyro::icntl[Param::DirBC_Interior]) || m_sc[smoother] == 1 ||
                (gyro::icntl[Param::extrapolation] && l == 0 && smoother % 2 == 0 &&
                 !(smoother == 0 && !gyro::icntl[Param::DirBC_Interior] && k == 0))) {
                u_sc = solve_diag(A_Zebra_v_LU_row[smoother][k], f_sc);
            }
            // Circle (not across)
            else if (smoother < 2 && !(smoother == 0 && !gyro::icntl[Param::DirBC_Interior] && k == 0)) {
                u_sc = solve_circle(A_Zebra_r_LU_row[smoother][k], A_Zebra_c_LU_row[smoother][k],
                                    A_Zebra_v_LU_row[smoother][k], f_sc);
            }
            // Radial (not across)
            else if (smoother > 1) {
                u_sc = solve_radial(A_Zebra_r_LU_row[smoother][k], A_Zebra_c_LU_row[smoother][k],
                                    A_Zebra_v_LU_row[smoother][k], f_sc);
            }
            // Across (direct solver)
            else {
#ifdef GMGPOLAR_USE_MUMPS
                if (gyro::icntl[Param::optimized] == 0) {
#endif
                    u_sc = solve_gaussian_elimination(A_Zebra_r_LU_row[smoother][k], A_Zebra_c_LU_row[smoother][k],
                                                      A_Zebra_v_LU_row[smoother][k], f_sc);
#ifdef GMGPOLAR_USE_MUMPS
                }
                else
                    u_sc = solve_mumps(mumps_across, f_sc);
#endif
            }

#pragma omp atomic
            t_Asc += TOC;
            TIC;

            build_fsc(u_sc, u, smoother, 1, k * m_sc[smoother], (k + 1) * m_sc[smoother]);

#pragma omp atomic
            t_f_sc += TOC;
            TIC;

            dep_Asc_cur[ind] = 1;
        }
        //        }
    }

    t = t_smoothing_tmp;
    t_smoothing += TOC;
} /* ----- end of level::multigrid_smoothing ----- */

/*! \brief Create the RHS part corresponding to Asc_ortho for a smoother (on level l)
 *
 * Create f_sc, i.e. the RHS part corresponding to Asc_ortho for a smoother (on level l)
 *
 * \param f_sc: the RHS part (out)
 * \param smoother: the smoother
*/
void level::build_fsc(std::vector<double>& f_sc, std::vector<double>& f, int smoother, int loc_to_glob, int start,
                      int end)
{
    std::vector<int> ind = mapping_usc_to_u(start, end, smoother);
    if (loc_to_glob) {
        //computation of indices in the total vector u corresponding to the indices in u_sc
        for (int i = start; i < end; ++i) {
            int index = ind[i - start];
            f[index]  = f_sc[i % m_sc[smoother]];
        }
    }
    else {
        //computation of indices in the total vector u corresponding to the indices in u_sc
        for (int i = start; i < end; ++i) {
            int index                = ind[i - start];
            f_sc[i % m_sc[smoother]] = f[index];
        }
    }
} /* ----- end of level::build_fsc ----- */

/*! \brief Build the matrix A_sc explicitely on the level l
 *
 * Asc corresponds to all lines of a smoother s with colour c
 * so we need to build 4 matrices Asc for the 4 smoothers on the level l, all stored in A_Zebra
 * 0,1: circular white/black
 * 2,3: radial white/black
 *
 * this function is called only once at the beginning to build the matrices Asc on all levels
 * the matrix should contain only the parts corresponding to (s/c),
 * thus it is not of size m, but smaller (number of points corresponding to the smoother s/c),
 * so we need to treat the indices (!)
*/
void level::build_Asc()
{
    int extrapol = gyro::icntl[Param::extrapolation] == 1 && l == 0;
    int shift_ij = extrapol ? 1 : 0;
    int smoother_prev, smoother_cur, smoother_next, smoother;
    int start_j, start_loop, shift_loop, ip, im, ptr, row, col;
    double coeff, coeff2, coeff3, kt, ktmin1, hs, hsmin1, beta_val;

    std::vector<int> smoother23(ntheta_int);
    std::vector<int> row_vect_prev, row_vect, row_vect_next;
    std::vector<int> ptr_vect_prev, ptr_vect, ptr_vect_next;
    std::vector<int> stencil_prev, stencil_cur, stencil_next, stencil;
    std::vector<double> arr_vect, arr_vect2, att_vect, art_vect;
    std::vector<std::vector<int>> stencil23_prev(2), stencil23_cur(2), stencil23_next(2);
    for (int i = 0; i < ntheta_int; i += 2) {
        smoother23[i]     = 2;
        smoother23[i + 1] = 3;
    }

    /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!      Circle Smoother     !!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     */
    // Take boundary condition into account: Dirichlet-RB
    smoother = 0;
    if (gyro::icntl[Param::DirBC_Interior]) { // (r[0],0) is on Dirichlet boundary
        for (int i = 0; i < ntheta_int / (shift_ij + 1); i++) {
            A_Zebra_r_row[smoother][0][i] = i;
            A_Zebra_c_row[smoother][0][i] = i;
            A_Zebra_v_row[smoother][0][i] = 1.0;
        }
        start_j = 1;
    }
    else { // (r[0],0) is not on Dirichlet boundary
        start_j = 0;
    }

    for (int j = start_j; j < delete_circles; j++) {
        /* Initialize moother, row_vect, ptr, and stencil */
        if (j > start_j) {
            smoother_prev = smoother_cur;
            smoother_cur  = smoother_next;
            row_vect_prev = row_vect;
            row_vect      = row_vect_next;
            ptr_vect_prev = ptr_vect;
            ptr_vect      = ptr_vect_next;
            stencil_prev  = stencil_cur;
            stencil_cur   = stencil_next;
        }
        else {
            smoother_cur = j;
            row_vect     = get_row(j, smoother_cur, extrapol, 1, 1);
            ptr_vect     = get_ptr_sc(j, smoother_cur, 0);
            stencil_cur  = get_stencil_sc(j, smoother_cur, 0);
        }
        if (j < delete_circles - 1)
            smoother_next = get_smoother(0, j + 1);
        else
            smoother_next = 2;
        row_vect_next = get_row(j + 1, smoother_next, extrapol, 1, 1);
        ptr_vect_next = get_ptr_sc(j + 1, smoother_next, 0);
        if (j < delete_circles - 1)
            stencil_next = get_stencil_sc(j + 1, smoother_next, 0);
        else {
            stencil23_next[0] = get_stencil_sc(j + 1, 2, 0);
            stencil23_next[1] = get_stencil_sc(j + 1, 3, 0);
        }
        hs = hplus[j];
        gyro::arr_att_art(r[j], theta, sin_theta, cos_theta, ntheta_int, arr_vect, att_vect, art_vect, 0);

        /* Update (j-1, i) */
        if (j > start_j) {
            hsmin1   = hplus[j - 1];
            smoother = smoother_prev;
            if (extrapol && smoother == 0) {
                start_loop = 1;
                shift_loop = 2;
            }
            else {
                start_loop = 0;
                shift_loop = 1;
            }
            for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                kt     = thetaplus_per[i + 1];
                ktmin1 = thetaplus_per[i];

                // Update (j-1, i)
                ptr     = ptr_vect_prev[i];
                stencil = stencil_prev;
                coeff   = 0.5 * (kt + ktmin1) * arr_vect[i] / hsmin1;
                col     = row_vect_prev[i];

                A_Zebra_v_row[smoother][(j - 1) / 2][ptr + stencil[Param::middle]] += coeff;
            }
        }
        else {
            smoother = smoother_cur;
            stencil  = stencil_cur;
            // Across: bottom update
            if (!gyro::icntl[Param::DirBC_Interior]) {
                hsmin1    = 2 * r[0];
                arr_vect2 = gyro::arr(r[j], theta_PI, sin_theta_PI, cos_theta_PI, ntheta_int, 0);
                for (int i = shift_ij; i < ntheta_int; i += shift_ij + 1) {
                    kt     = thetaplus_per[i + 1];
                    ktmin1 = thetaplus_per[i];

                    ptr = ptr_vect[i];
                    row = row_vect[i];
                    col = row_vect[(i + ntheta_int / 2) % ntheta_int];
                    A_Zebra_r_row[smoother][j / 2][ptr + stencil[Param::bottom]] = row;
                    A_Zebra_c_row[smoother][j / 2][ptr + stencil[Param::bottom]] = col;

                    coeff  = 0.5 * (kt + ktmin1) * arr_vect2[i] / hsmin1;
                    coeff2 = 0.5 * (kt + ktmin1) * arr_vect[i] / hsmin1;

                    A_Zebra_v_row[smoother][j / 2][ptr + stencil[Param::bottom]] += -coeff - coeff2;

                    A_Zebra_v_row[smoother][j / 2][ptr + stencil[Param::middle]] += coeff;
                }
            }
            // Dirichlet
            else {
                hsmin1 = hplus[j - 1];
                // DB contribution arr (r(0))
                arr_vect2 = gyro::arr(r[j - 1], theta, sin_theta, cos_theta, ntheta_int, 0);

                for (int i = 0; i < ntheta_int; i++) {
                    kt     = thetaplus_per[i + 1];
                    ktmin1 = thetaplus_per[i];

                    ptr   = ptr_vect[i];
                    row   = row_vect[i];
                    coeff = 0.5 * (kt + ktmin1) * arr_vect2[i] / hsmin1;

                    A_Zebra_v_row[smoother][j / 2][ptr + stencil[Param::middle]] += coeff;
                }
            }
        }

        /* Update (j, x) */
        smoother = smoother_cur;
        if (!extrapol || smoother == 1) {
            for (int i = 0; i < ntheta_int; i++) {
                kt      = thetaplus_per[i + 1];
                ktmin1  = thetaplus_per[i];
                stencil = stencil_cur;
                ip      = iplus(i, ntheta_int);
                im      = imoins(i, ntheta_int);

                // Update (j, i)
                ptr                                                          = ptr_vect[i];
                row                                                          = row_vect[i];
                A_Zebra_r_row[smoother][j / 2][ptr + stencil[Param::middle]] = row;
                A_Zebra_c_row[smoother][j / 2][ptr + stencil[Param::middle]] = row;

                coeff  = 0.5 * (kt + ktmin1) * arr_vect[i];
                coeff2 = 0.5 * (hs + hsmin1) * att_vect[i];
                col    = row_vect[im];

                A_Zebra_v_row[smoother][j / 2][ptr + stencil[Param::left]] += -coeff2 / ktmin1;
                col = row_vect[ip];

                A_Zebra_v_row[smoother][j / 2][ptr + stencil[Param::right]] += -coeff2 / kt;

                A_Zebra_v_row[smoother][j / 2][ptr + stencil[Param::middle]] +=
                    coeff / hs + coeff / hsmin1 + coeff2 / ktmin1 + coeff2 / kt;

                // Beta coefficient
                beta_val = betaVec[j * ntheta_int + i];

                A_Zebra_v_row[smoother][j / 2][ptr + stencil[Param::middle]] += beta_val;

                // Update (j, i+1)
                ptr                                                        = ptr_vect[ip];
                row                                                        = row_vect[ip];
                col                                                        = row_vect[i];
                coeff                                                      = 0.5 * (hs + hsmin1) * att_vect[i] / kt;
                A_Zebra_r_row[smoother][j / 2][ptr + stencil[Param::left]] = row;
                A_Zebra_c_row[smoother][j / 2][ptr + stencil[Param::left]] = col;

                A_Zebra_v_row[smoother][j / 2][ptr + stencil[Param::left]] += -coeff;

                A_Zebra_v_row[smoother][j / 2][ptr + stencil[Param::middle]] += coeff;

                // Update (j, i-1)
                ptr   = ptr_vect[im];
                row   = row_vect[im];
                col   = row_vect[i];
                coeff = 0.5 * (hs + hsmin1) * att_vect[i] / ktmin1;
                A_Zebra_r_row[smoother][j / 2][ptr + stencil[Param::right]] = row;
                A_Zebra_c_row[smoother][j / 2][ptr + stencil[Param::right]] = col;

                A_Zebra_v_row[smoother][j / 2][ptr + stencil[Param::right]] += -coeff;

                A_Zebra_v_row[smoother][j / 2][ptr + stencil[Param::middle]] += coeff;
            }
        }
        if (extrapol && smoother == 0) {
            start_loop = 0;
            shift_loop = 2;
            for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                kt      = thetaplus_per[i + 1];
                ktmin1  = thetaplus_per[i];
                stencil = stencil_cur;
                ip      = iplus(i, ntheta_int);
                im      = imoins(i, ntheta_int);
                col     = row_vect[i];

                // Update (j, i+1)
                ptr   = ptr_vect[iplus(i, ntheta_int)];
                coeff = 0.5 * (hs + hsmin1) * att_vect[i] / kt;
                row   = row_vect[ip];

                A_Zebra_v_row[smoother][j / 2][ptr + stencil[Param::middle]] += coeff;

                // Update (j, i-1)
                ptr   = ptr_vect[imoins(i, ntheta_int)];
                coeff = 0.5 * (hs + hsmin1) * att_vect[i] / ktmin1;
                row   = row_vect[im];

                A_Zebra_v_row[smoother][j / 2][ptr + stencil[Param::middle]] += coeff;

                // Update (j, i) [for i+1]
                kt                                                           = thetaplus_per[i + 2];
                ktmin1                                                       = thetaplus_per[i + 1];
                ptr                                                          = ptr_vect[i + 1];
                coeff                                                        = 0.5 * (kt + ktmin1) * arr_vect[i + 1];
                coeff2                                                       = 0.5 * (hs + hsmin1) * att_vect[i + 1];
                row                                                          = row_vect[i + 1];
                A_Zebra_r_row[smoother][j / 2][ptr + stencil[Param::middle]] = row;
                A_Zebra_c_row[smoother][j / 2][ptr + stencil[Param::middle]] = row;

                A_Zebra_v_row[smoother][j / 2][ptr + stencil[Param::middle]] +=
                    coeff / hs + coeff / hsmin1 + coeff2 / ktmin1 + coeff2 / kt;

                // Beta coefficient
                beta_val = betaVec[j * ntheta_int + ip];

                A_Zebra_v_row[smoother][j / 2][ptr + stencil[Param::middle]] += beta_val;
            }
        }

        /* Update (j+1, i) */
        if (j < delete_circles - 1) {
            smoother = smoother_next;
            if (extrapol && smoother == 0) {
                start_loop = 1;
                shift_loop = 2;
            }
            else {
                start_loop = 0;
                shift_loop = 1;
            }
            for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                kt     = thetaplus_per[i + 1];
                ktmin1 = thetaplus_per[i];

                ptr     = ptr_vect_next[i];
                stencil = stencil_next;
                coeff   = 0.5 * (kt + ktmin1) * arr_vect[i] / hs;
                row     = row_vect_next[i];

                A_Zebra_v_row[smoother][(j + 1) / 2][ptr + stencil[Param::middle]] += coeff;
            }
        }
        else if (j == delete_circles - 1) {
            if (extrapol && (j + 1) % 2 == 0) {
                start_loop = 1;
                shift_loop = 2;
            }
            else {
                start_loop = 0;
                shift_loop = 1;
            }
            for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                kt     = thetaplus_per[i + 1];
                ktmin1 = thetaplus_per[i];

                ptr      = ptr_vect_next[i];
                smoother = smoother23[i];
                stencil  = stencil23_next[smoother - 2];
                coeff    = 0.5 * (kt + ktmin1) * arr_vect[i] / hs;
                row      = row_vect_next[i];

                A_Zebra_v_row[smoother][i / 2][ptr + stencil[Param::middle]] += coeff;
            }
        }
    }

    /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!      Radial Smoother     !!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     */
    for (int j = delete_circles; j < nr_int; j++) {
        /* Initialize moother, row_vect, ptr, and stencil */
        if (j == delete_circles) {
            smoother_prev = get_smoother(0, j - 1);
            stencil_prev  = get_stencil_sc(j - 1, smoother_prev, 0);
            row_vect_prev = get_row(j - 1, smoother_prev, extrapol, 1, 1);
            ptr_vect_prev = get_ptr_sc(j - 1, smoother_prev, 0);
        }
        else {
            stencil23_prev[0] = get_stencil_sc(j - 1, 2, 0);
            stencil23_prev[1] = get_stencil_sc(j - 1, 3, 0);
            row_vect_prev     = get_row(j - 1, 2, extrapol, 1, 1);
            ptr_vect_prev     = get_ptr_sc(j - 1, 2, 0);
        }
        row_vect         = row_vect_next;
        ptr_vect         = ptr_vect_next;
        stencil23_cur[0] = get_stencil_sc(j, 2, 0);
        stencil23_cur[1] = get_stencil_sc(j, 3, 0);
        if (j < nr_int - 1) {
            row_vect_next     = get_row(j + 1, 2, extrapol, 1, 1);
            ptr_vect_next     = get_ptr_sc(j + 1, 2, 0);
            stencil23_next[0] = get_stencil_sc(j + 1, 2, 0);
            stencil23_next[1] = get_stencil_sc(j + 1, 3, 0);
        }
        hs     = hplus[j];
        hsmin1 = hplus[j - 1];
        gyro::arr_att_art(r[j], theta, sin_theta, cos_theta, ntheta_int, arr_vect, att_vect, art_vect, 0);
        if (j == nr_int - 1)
            arr_vect2 = gyro::arr(r[j + 1], theta, sin_theta, cos_theta, ntheta_int, 0);
        /* Update (j-1, i) */
        if (j == delete_circles) {
            smoother = smoother_prev;
            stencil  = stencil_prev;
            if (extrapol && smoother == 0) {
                start_loop = 1;
                shift_loop = 2;
            }
            else {
                start_loop = 0;
                shift_loop = 1;
            }
            for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                kt     = thetaplus_per[i + 1];
                ktmin1 = thetaplus_per[i];

                ptr   = ptr_vect_prev[i];
                coeff = 0.5 * (kt + ktmin1) * arr_vect[i] / hsmin1;
                row   = row_vect_prev[i];

                A_Zebra_v_row[smoother][(j - 1) / 2][ptr + stencil[Param::middle]] += coeff;
            }
        }
        else {
            if (extrapol) {
                start_loop = 1;
                shift_loop = 2;
            }
            else {
                start_loop = 0;
                shift_loop = 1;
            }
            for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                kt       = thetaplus_per[i + 1];
                ktmin1   = thetaplus_per[i];
                smoother = smoother23[i];

                ptr                                                       = ptr_vect_prev[i];
                stencil                                                   = stencil23_prev[smoother - 2];
                coeff                                                     = 0.5 * (kt + ktmin1) * arr_vect[i] / hsmin1;
                row                                                       = row_vect_prev[i];
                col                                                       = row_vect[i];
                A_Zebra_r_row[smoother][i / 2][ptr + stencil[Param::top]] = row;
                A_Zebra_c_row[smoother][i / 2][ptr + stencil[Param::top]] = col;

                A_Zebra_v_row[smoother][i / 2][ptr + stencil[Param::top]] += -coeff;

                A_Zebra_v_row[smoother][i / 2][ptr + stencil[Param::middle]] += coeff;
            }
            start_loop = 0;
            shift_loop = 2;
            if (extrapol && j % 2 == 0) {
                for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                    kt       = thetaplus_per[i + 1];
                    ktmin1   = thetaplus_per[i];
                    smoother = smoother23[i];

                    ptr     = ptr_vect_prev[i];
                    stencil = stencil23_prev[smoother - 2];
                    coeff   = 0.5 * (kt + ktmin1) * arr_vect[i] / hsmin1;
                    row     = row_vect_prev[i];

                    A_Zebra_v_row[smoother][i / 2][ptr + stencil[Param::middle]] += coeff;
                }
            }
        }

        /* Update (j, x) */
        if (extrapol) {
            start_loop = 1;
            shift_loop = 2;
        }
        else {
            start_loop = 0;
            shift_loop = 1;
        }
        if (j == delete_circles) {
            for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                kt       = thetaplus_per[i + 1];
                ktmin1   = thetaplus_per[i];
                smoother = smoother23[i];
                // Update (j, i)
                ptr                                                          = ptr_vect[i];
                stencil                                                      = stencil23_cur[smoother - 2];
                coeff                                                        = 0.5 * (kt + ktmin1) * arr_vect[i];
                coeff2                                                       = 0.5 * (hs + hsmin1) * att_vect[i];
                row                                                          = row_vect[i];
                A_Zebra_r_row[smoother][i / 2][ptr + stencil[Param::middle]] = row;
                A_Zebra_c_row[smoother][i / 2][ptr + stencil[Param::middle]] = row;
                col                                                          = row_vect_prev[i];

                A_Zebra_v_row[smoother][i / 2][ptr + stencil[Param::top]] += -coeff / hs;

                A_Zebra_v_row[smoother][i / 2][ptr + stencil[Param::middle]] +=
                    coeff / hs + coeff / hsmin1 + coeff2 / ktmin1 + coeff2 / kt;

                // Beta coefficient
                beta_val = betaVec[j * ntheta_int + i];

                A_Zebra_v_row[smoother][i / 2][ptr + stencil[Param::middle]] += beta_val;
            }
        }
        else if (j == nr_int - 1) {
            for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                kt       = thetaplus_per[i + 1];
                ktmin1   = thetaplus_per[i];
                smoother = smoother23[i];

                // Update (j, i)
                ptr     = ptr_vect[i];
                stencil = stencil23_cur[smoother - 2];
                coeff   = 0.5 * (kt + ktmin1) * arr_vect[i];
                // Contribution to middle (top) from DB
                coeff3                                                       = 0.5 * (kt + ktmin1) * arr_vect2[i] / hs;
                coeff2                                                       = 0.5 * (hs + hsmin1) * att_vect[i];
                row                                                          = row_vect[i];
                A_Zebra_r_row[smoother][i / 2][ptr + stencil[Param::middle]] = row;
                A_Zebra_c_row[smoother][i / 2][ptr + stencil[Param::middle]] = row;
                col                                                          = row_vect_next[i];

                A_Zebra_v_row[smoother][i / 2][ptr + stencil[Param::bottom]] += -coeff / hsmin1;

                A_Zebra_v_row[smoother][i / 2][ptr + stencil[Param::middle]] +=
                    coeff / hsmin1 + coeff / hs + coeff3 + coeff2 / ktmin1 + coeff2 / kt;

                // Beta coefficient
                beta_val = betaVec[j * ntheta_int + i];

                A_Zebra_v_row[smoother][i / 2][ptr + stencil[Param::middle]] += beta_val;
            }
        }
        else {
            for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                kt       = thetaplus_per[i + 1];
                ktmin1   = thetaplus_per[i];
                smoother = smoother23[i];

                // Update (j, i)
                ptr                                                          = ptr_vect[i];
                stencil                                                      = stencil23_cur[smoother - 2];
                coeff                                                        = 0.5 * (kt + ktmin1) * arr_vect[i];
                coeff2                                                       = 0.5 * (hs + hsmin1) * att_vect[i];
                row                                                          = row_vect[i];
                A_Zebra_r_row[smoother][i / 2][ptr + stencil[Param::middle]] = row;
                A_Zebra_c_row[smoother][i / 2][ptr + stencil[Param::middle]] = row;
                col                                                          = row_vect_next[i];

                A_Zebra_v_row[smoother][i / 2][ptr + stencil[Param::top]] += -coeff / hs;

                A_Zebra_v_row[smoother][i / 2][ptr + stencil[Param::middle]] +=
                    coeff / hs + coeff / hsmin1 + coeff2 / ktmin1 + coeff2 / kt;
                col = row_vect_prev[i];

                A_Zebra_v_row[smoother][i / 2][ptr + stencil[Param::bottom]] += -coeff / hsmin1;

                // Beta coefficient
                beta_val = betaVec[j * ntheta_int + i];

                A_Zebra_v_row[smoother][i / 2][ptr + stencil[Param::middle]] += beta_val;
            }
        }
        start_loop = 0;
        shift_loop = 2;
        if (extrapol && j % 2 == 1) {
            if (j == nr_int - 1) {
                for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                    kt     = thetaplus_per[i + 1];
                    ktmin1 = thetaplus_per[i];

                    // Update (j, i)
                    smoother = smoother23[i];
                    stencil  = stencil23_cur[smoother - 2];
                    ptr      = ptr_vect[i];
                    coeff    = 0.5 * (kt + ktmin1) * arr_vect[i];
                    // Contribution to middle (top) from DB
                    coeff3 = 0.5 * (kt + ktmin1) * arr_vect2[i] / hs;
                    coeff2 = 0.5 * (hs + hsmin1) * att_vect[i];
                    row    = row_vect[i];
                    A_Zebra_r_row[smoother][i / 2][ptr + stencil[Param::middle]] = row;
                    A_Zebra_c_row[smoother][i / 2][ptr + stencil[Param::middle]] = row;

                    A_Zebra_v_row[smoother][i / 2][ptr + stencil[Param::middle]] +=
                        coeff / hsmin1 + coeff / hs + coeff3 + coeff2 / ktmin1 + coeff2 / kt;

                    // Beta coefficient
                    beta_val = betaVec[j * ntheta_int + i];

                    A_Zebra_v_row[smoother][i / 2][ptr + stencil[Param::middle]] += beta_val;
                }
            }
            else {
                for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                    kt     = thetaplus_per[i + 1];
                    ktmin1 = thetaplus_per[i];

                    // Update (j, i)
                    smoother                                                     = smoother23[i];
                    stencil                                                      = stencil23_cur[smoother - 2];
                    ptr                                                          = ptr_vect[i];
                    coeff                                                        = 0.5 * (kt + ktmin1) * arr_vect[i];
                    coeff2                                                       = 0.5 * (hs + hsmin1) * att_vect[i];
                    row                                                          = row_vect[i];
                    A_Zebra_r_row[smoother][i / 2][ptr + stencil[Param::middle]] = row;
                    A_Zebra_c_row[smoother][i / 2][ptr + stencil[Param::middle]] = row;

                    A_Zebra_v_row[smoother][i / 2][ptr + stencil[Param::middle]] +=
                        coeff / hs + coeff / hsmin1 + coeff2 / ktmin1 + coeff2 / kt;

                    // Beta coefficient
                    beta_val = betaVec[j * ntheta_int + i];

                    A_Zebra_v_row[smoother][i / 2][ptr + stencil[Param::middle]] += beta_val;
                }
            }
        }
        if (extrapol && j % 2 == 0) {
            start_loop = 0;
            shift_loop = 2;
        }
        else {
            start_loop = 0;
            shift_loop = 1;
        }
        for (int i = start_loop; i < ntheta_int; i += shift_loop) {
            kt       = thetaplus_per[i + 1];
            ktmin1   = thetaplus_per[i];
            ip       = iplus(i, ntheta_int);
            im       = imoins(i, ntheta_int);
            smoother = smoother23[ip];
            stencil  = stencil23_cur[smoother - 2];

            // Update (j, i+1)
            ptr   = ptr_vect[ip];
            coeff = 0.5 * (hs + hsmin1) * att_vect[i] / kt;
            row   = row_vect[ip];

            A_Zebra_v_row[smoother][ip / 2][ptr + stencil[Param::middle]] += coeff;

            // Update (j, i-1)
            smoother = smoother23[im];
            ptr      = ptr_vect[im];
            coeff    = 0.5 * (hs + hsmin1) * att_vect[i] / ktmin1;
            row      = row_vect[im];

            A_Zebra_v_row[smoother][im / 2][ptr + stencil[Param::middle]] += coeff;
        }

        /* Update (j+1, i) */
        if (j < nr_int - 1) {
            if (extrapol) {
                start_loop = 1;
                shift_loop = 2;
            }
            else {
                start_loop = 0;
                shift_loop = 1;
            }
            for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                kt       = thetaplus_per[i + 1];
                ktmin1   = thetaplus_per[i];
                smoother = smoother23[i];

                ptr                                                          = ptr_vect_next[i];
                stencil                                                      = stencil23_next[smoother - 2];
                coeff                                                        = 0.5 * (kt + ktmin1) * arr_vect[i] / hs;
                row                                                          = row_vect_next[i];
                col                                                          = row_vect[i];
                A_Zebra_r_row[smoother][i / 2][ptr + stencil[Param::bottom]] = row;
                A_Zebra_c_row[smoother][i / 2][ptr + stencil[Param::bottom]] = col;

                A_Zebra_v_row[smoother][i / 2][ptr + stencil[Param::bottom]] += -coeff;

                A_Zebra_v_row[smoother][i / 2][ptr + stencil[Param::middle]] += coeff;
            }
            start_loop = 0;
            shift_loop = 2;
            if (extrapol && j % 2 == 0) {
                for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                    kt       = thetaplus_per[i + 1];
                    ktmin1   = thetaplus_per[i];
                    smoother = smoother23[i];

                    ptr     = ptr_vect_next[i];
                    stencil = stencil23_next[smoother - 2];
                    coeff   = 0.5 * (kt + ktmin1) * arr_vect[i] / hs;
                    row     = row_vect_next[i];

                    A_Zebra_v_row[smoother][i / 2][ptr + stencil[Param::middle]] += coeff;
                }
            }
        }
    }

    // Take boundary condition into account: Dirichlet-RB
    int j    = nr_int;
    ptr_vect = get_ptr_sc(j, 2, 0);
    row_vect = get_row(j, 2, extrapol, 1, 1);
    for (int i = shift_ij; i < ntheta_int; i += shift_ij + 1) {
        ptr      = ptr_vect[i];
        smoother = smoother23[i];
        // row                      = get_row_radial(i, nr_int, smoother, extrapol, ntheta_int, delete_circles);
        row                                 = row_vect[i];
        A_Zebra_r_row[smoother][i / 2][ptr] = row;
        A_Zebra_c_row[smoother][i / 2][ptr] = row;
        A_Zebra_v_row[smoother][i / 2][ptr] += 1.0;
    }
} /* ----- end of level::build_Asc ----- */

/*! \brief Applies the matrix A_sc_ortho for a smoother explicitely on the level l
 *
 * Asc_ortho corresponds to all lines of a smoother s with colour c and the columns not in Asc
 *
 * \param smoother_todo: the smoother of this Asc_ortho matrix
*/
void level::build_Asc_ortho(int smoother_todo)
{
    int start_j, j, odd_j;
    int extrapol = gyro::icntl[Param::extrapolation] == 1 && l == 0;

    int smoother = smoother_todo;
    int smoother_prev, smoother_cur, smoother_next;
    int smoother_prev_circle, smoother_cur_left, smoother_cur_circle, smoother_cur_right, smoother_next_circle;
    int start_loop, shift_loop;
    int row, row2, col, col2, ptr;
    double coeff, coeff2, coeff3, kt, ktmin1, hs, hsmin1, val, val2;
    int base_prec;
    if (smoother == 1)
        base_prec = delete_circles + 1;
    else if (smoother == 0)
        base_prec = -(delete_circles + 1);
    else if (smoother == 3)
        base_prec = ntheta_int;
    else if (smoother == 2)
        base_prec = -ntheta_int;

    std::vector<int> row_vect_prev, row_vect, row_vect_next;
    std::vector<int> smoother_vect_prev, smoother_vect, smoother_vect_next;
    std::vector<int> ptr_vect_prev, ptr_vect, ptr_vect_next;
    std::vector<int> stencil_prev, stencil_cur, stencil_cur_left, stencil_cur_right, stencil_next, stencil;
    std::vector<double> arr_vect, arr_vect2, att_vect, art_vect;

    /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!      Circle Smoother     !!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     */
    // Take boundary condition into account: Dirichlet-RB
    if (gyro::icntl[Param::DirBC_Interior]) { // (r[0],0) is on Dirichlet boundary
        dep_Asc_ortho[0][0] = 1;
        dep_Asc_ortho[1][0] = 1;
        start_j             = 1;
    }
    else { // (r[0],0) is not on Dirichlet boundary
        start_j = 0;
    }

    if (smoother < 2) {
        /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!! Link Radial-Circle (circle) !!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     */
        j     = delete_circles;
        odd_j = j % 2;

        smoother_prev      = get_smoother(0, j - 1);
        smoother_vect      = get_smoother_radial(j);
        smoother_vect_next = get_smoother_radial(j + 1);
        ptr_vect_prev      = get_ptr_sc(j - 1, smoother_prev, 1);
        ptr_vect           = get_ptr_sc(j, smoother_vect[0], 1);
        ptr_vect_next      = get_ptr_sc(j + 1, smoother_vect_next[0], 1);
        stencil_prev       = get_stencil_sc(j - 1, smoother_prev, 1);
        stencil_cur        = get_stencil_sc(j, smoother_vect[0], 1);
        stencil_next       = get_stencil_sc(j + 1, smoother_vect_next[0], 1);
        row_vect_prev      = get_row(j - 1, smoother_prev, extrapol, 0, 1);
        row_vect           = get_row(j, 2, extrapol, 0, 1);
        row_vect_next      = get_row(j + 1, 2, extrapol, 0, 1);
        hs                 = hplus[j];
        hsmin1             = hplus[j - 1];
        gyro::arr_att_art(r[j], theta, sin_theta, cos_theta, ntheta_int, arr_vect, att_vect, art_vect, 0);
        if (smoother_prev == smoother) {
            if (extrapol && smoother == 0) {
                start_loop = 1;
                shift_loop = 2;
            }
            else {
                start_loop = 0;
                shift_loop = 1;
            }
            stencil = stencil_prev;
            for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                kt     = thetaplus_per[i + 1];
                ktmin1 = thetaplus_per[i];
                ptr    = ptr_vect_prev[i];

                // Update (j-1, i)
                row                                                = row_vect_prev[i];
                col                                                = j * ntheta_int + i;
                coeff                                              = 0.5 * (kt + ktmin1) * arr_vect[i] / hsmin1;
                val                                                = -coeff;
                A_Zebra_Mix_r[smoother][ptr + stencil[Param::top]] = row;
                A_Zebra_Mix_c[smoother][ptr + stencil[Param::top]] = col;
                A_Zebra_Mix_v[smoother][ptr + stencil[Param::top]] += val;
            }
            if (gyro::icntl[Param::mod_pk] > 0) {
                if (extrapol && smoother == 0) {
                    start_loop = 1;
                    shift_loop = 2;
                }
                else {
                    start_loop = 0;
                    shift_loop = 1;
                }
                for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                    kt     = thetaplus_per[i + 1];
                    ktmin1 = thetaplus_per[i];
                    ptr    = ptr_vect_prev[i];

                    // Update (j-1, i)
                    row                                                     = row_vect_prev[i];
                    col                                                     = j * ntheta_int + imoins(i, ntheta_int);
                    val                                                     = art_vect[i];
                    A_Zebra_Mix_r[smoother][ptr + stencil[Param::top_left]] = row;
                    A_Zebra_Mix_c[smoother][ptr + stencil[Param::top_left]] = col;
                    A_Zebra_Mix_v[smoother][ptr + stencil[Param::top_left]] += val;
                    col2                                                     = j * ntheta_int + iplus(i, ntheta_int);
                    val2                                                     = -art_vect[i];
                    A_Zebra_Mix_r[smoother][ptr + stencil[Param::top_right]] = row;
                    A_Zebra_Mix_c[smoother][ptr + stencil[Param::top_right]] = col2;
                    A_Zebra_Mix_v[smoother][ptr + stencil[Param::top_right]] += val2;
                }
            }
        }

        /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!! Link Circle-Radial !!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     */
        j     = delete_circles - 1;
        odd_j = j % 2;

        if (delete_circles > 2 || !gyro::icntl[Param::DirBC_Interior]) {
            smoother_prev      = get_smoother(0, j - 1);
            smoother_cur       = get_smoother(0, j);
            smoother_vect_next = get_smoother_radial(j + 1);
            ptr_vect_prev      = get_ptr_sc(j - 1, smoother_prev, 1);
            ptr_vect           = get_ptr_sc(j, smoother_cur, 1);
            ptr_vect_next      = get_ptr_sc(j + 1, smoother_vect_next[0], 1);
            stencil_prev       = get_stencil_sc(j - 1, smoother_prev, 1);
            stencil_cur        = get_stencil_sc(j, smoother_cur, 1);
            stencil_next       = get_stencil_sc(j + 1, smoother_vect_next[0], 1);
            row_vect_prev      = get_row(j - 1, smoother_prev, extrapol, 0, 1);
            row_vect           = get_row(j, smoother_cur, extrapol, 0, 1);
            row_vect_next      = get_row(j + 1, 2, extrapol, 0, 1);
            hs                 = hplus[j];
            hsmin1             = hplus[j - 1];
            gyro::arr_att_art(r[j], theta, sin_theta, cos_theta, ntheta_int, arr_vect, att_vect, art_vect, 0);
            if (smoother_cur == smoother) {
                if (extrapol && smoother == 0) {
                    start_loop = 1;
                    shift_loop = 2;
                }
                else {
                    start_loop = 0;
                    shift_loop = 1;
                }
                stencil = stencil_cur;
                for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                    kt     = thetaplus_per[i + 1];
                    ktmin1 = thetaplus_per[i];
                    ptr    = ptr_vect[i];

                    // Update (j, i)
                    row   = row_vect[i];
                    coeff = 0.5 * (kt + ktmin1) * arr_vect[i];
                    col   = (j + 1) * ntheta_int + i;
                    val   = -coeff / hs;
                    A_Zebra_Mix_v[smoother][ptr + stencil[Param::top]] += val;
                    col2 = (j - 1) * ntheta_int + i;
                    val2 = -coeff / hsmin1;
                    A_Zebra_Mix_v[smoother][ptr + stencil[Param::bottom]] += val2;
                }
                if (extrapol && smoother == 0) {
                    start_loop = 0;
                    shift_loop = 2;
                    for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                        kt     = thetaplus_per[i + 1];
                        ktmin1 = thetaplus_per[i];

                        // Update (j, i+1)
                        ptr                                                 = ptr_vect[iplus(i, ntheta_int)];
                        row                                                 = row_vect[iplus(i, ntheta_int)];
                        col                                                 = j * ntheta_int + i;
                        coeff                                               = 0.5 * (hs + hsmin1) * att_vect[i];
                        val                                                 = -coeff / kt;
                        A_Zebra_Mix_r[smoother][ptr + stencil[Param::left]] = row;
                        A_Zebra_Mix_c[smoother][ptr + stencil[Param::left]] = col;
                        A_Zebra_Mix_v[smoother][ptr + stencil[Param::left]] += val;

                        // Update (j, i-1)
                        ptr                                                  = ptr_vect[imoins(i, ntheta_int)];
                        row                                                  = row_vect[imoins(i, ntheta_int)];
                        val                                                  = -coeff / ktmin1;
                        A_Zebra_Mix_r[smoother][ptr + stencil[Param::right]] = row;
                        A_Zebra_Mix_c[smoother][ptr + stencil[Param::right]] = col;
                        A_Zebra_Mix_v[smoother][ptr + stencil[Param::right]] += val;

                        ptr    = ptr_vect[i + 1];
                        kt     = thetaplus_per[i + 2];
                        ktmin1 = thetaplus_per[i + 1];
                        row    = row_vect[i + 1];
                        coeff2 = 0.5 * (hs + hsmin1) * att_vect[i + 1];
                        col    = j * ntheta_int + imoins(i + 1, ntheta_int);
                        val    = -coeff2 / ktmin1;
                        A_Zebra_Mix_v[smoother][ptr + stencil[Param::left]] += val;
                        col2 = j * ntheta_int + iplus(i + 1, ntheta_int);
                        val2 = -coeff2 / kt;
                        A_Zebra_Mix_v[smoother][ptr + stencil[Param::right]] += val2;
                    }
                }
            }
            else if (smoother_prev == smoother) {
                if (extrapol && smoother == 0) {
                    start_loop = 1;
                    shift_loop = 2;
                }
                else {
                    start_loop = 0;
                    shift_loop = 1;
                }
                stencil = stencil_prev;
                for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                    kt     = thetaplus_per[i + 1];
                    ktmin1 = thetaplus_per[i];
                    ptr    = ptr_vect_prev[i];

                    // Update (j-1, i)
                    row                                                = row_vect_prev[i];
                    col                                                = j * ntheta_int + i;
                    coeff                                              = 0.5 * (kt + ktmin1) * arr_vect[i] / hsmin1;
                    val                                                = -coeff;
                    A_Zebra_Mix_r[smoother][ptr + stencil[Param::top]] = row;
                    A_Zebra_Mix_c[smoother][ptr + stencil[Param::top]] = col;
                    A_Zebra_Mix_v[smoother][ptr + stencil[Param::top]] += val;
                }
            }
            if (gyro::icntl[Param::mod_pk] > 0) {
                if (smoother_cur == smoother) {
                    start_loop = 0;
                    if (extrapol && smoother == 0)
                        shift_loop = 2;
                    else
                        shift_loop = 1;
                    stencil = stencil_cur;
                    for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                        kt     = thetaplus_per[i + 1];
                        ktmin1 = thetaplus_per[i];

                        // Update (j, i+1)
                        ptr                                                        = ptr_vect[iplus(i, ntheta_int)];
                        row                                                        = row_vect[iplus(i, ntheta_int)];
                        col                                                        = (j - 1) * ntheta_int + i;
                        val                                                        = -art_vect[i];
                        A_Zebra_Mix_r[smoother][ptr + stencil[Param::bottom_left]] = row;
                        A_Zebra_Mix_c[smoother][ptr + stencil[Param::bottom_left]] = col;
                        A_Zebra_Mix_v[smoother][ptr + stencil[Param::bottom_left]] += val;
                        col2                                                    = (j + 1) * ntheta_int + i;
                        val2                                                    = art_vect[i];
                        A_Zebra_Mix_r[smoother][ptr + stencil[Param::top_left]] = row;
                        A_Zebra_Mix_c[smoother][ptr + stencil[Param::top_left]] = col2;
                        A_Zebra_Mix_v[smoother][ptr + stencil[Param::top_left]] += val2;

                        // Update (j, i-1)
                        ptr                                                         = ptr_vect[imoins(i, ntheta_int)];
                        row2                                                        = row_vect[imoins(i, ntheta_int)];
                        A_Zebra_Mix_r[smoother][ptr + stencil[Param::bottom_right]] = row2;
                        A_Zebra_Mix_c[smoother][ptr + stencil[Param::bottom_right]] = col;
                        A_Zebra_Mix_v[smoother][ptr + stencil[Param::bottom_right]] -= val;

                        A_Zebra_Mix_r[smoother][ptr + stencil[Param::top_right]] = row2;
                        A_Zebra_Mix_c[smoother][ptr + stencil[Param::top_right]] = col2;
                        A_Zebra_Mix_v[smoother][ptr + stencil[Param::top_right]] -= val2;
                    }
                }
                else if (smoother_prev == smoother) {
                    stencil = stencil_prev;
                    if (extrapol && smoother == 0) {
                        start_loop = 1;
                        shift_loop = 2;
                    }
                    else {
                        start_loop = 0;
                        shift_loop = 1;
                    }
                    for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                        kt     = thetaplus_per[i + 1];
                        ktmin1 = thetaplus_per[i];
                        ptr    = ptr_vect_prev[i];

                        // Update (j-1, i)
                        row = row_vect_prev[i];
                        col = j * ntheta_int + imoins(i, ntheta_int);
                        val = art_vect[i];
                        A_Zebra_Mix_r[smoother][ptr + stencil[Param::top_left]] = row;
                        A_Zebra_Mix_c[smoother][ptr + stencil[Param::top_left]] = col;
                        A_Zebra_Mix_v[smoother][ptr + stencil[Param::top_left]] += val;
                        col2 = j * ntheta_int + iplus(i, ntheta_int);
                        val2 = -art_vect[i];
                        A_Zebra_Mix_r[smoother][ptr + stencil[Param::top_right]] = row;
                        A_Zebra_Mix_c[smoother][ptr + stencil[Param::top_right]] = col2;
                        A_Zebra_Mix_v[smoother][ptr + stencil[Param::top_right]] += val2;
                    }
                }
            }
        }

        /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!  Interior nodes   !!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     */
        for (int j = delete_circles - 2; j > start_j; j--) {
            odd_j = j % 2;

            smoother_prev = get_smoother(0, j - 1);
            smoother_cur  = get_smoother(0, j);
            smoother_next = get_smoother(0, j + 1);
            ptr_vect_prev = get_ptr_sc(j - 1, smoother_prev, 1);
            ptr_vect      = get_ptr_sc(j, smoother_cur, 1);
            ptr_vect_next = get_ptr_sc(j + 1, smoother_next, 1);
            stencil_prev  = get_stencil_sc(j - 1, smoother_prev, 1);
            stencil_cur   = get_stencil_sc(j, smoother_cur, 1);
            stencil_next  = get_stencil_sc(j + 1, smoother_next, 1);
            row_vect_prev = get_row(j - 1, smoother_prev, extrapol, 0, 1);
            row_vect      = get_row(j, smoother_cur, extrapol, 0, 1);
            row_vect_next = get_row(j + 1, smoother_next, extrapol, 0, 1);
            hs            = hplus[j];
            hsmin1        = hplus[j - 1];
            gyro::arr_att_art(r[j], theta, sin_theta, cos_theta, ntheta_int, arr_vect, att_vect, art_vect, 0);
            if (smoother_cur == smoother) {
                stencil = stencil_cur;
                if (extrapol && smoother == 0) {
                    start_loop = 1;
                    shift_loop = 2;
                }
                else {
                    start_loop = 0;
                    shift_loop = 1;
                }
                for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                    kt     = thetaplus_per[i + 1];
                    ktmin1 = thetaplus_per[i];
                    ptr    = ptr_vect[i];

                    // Update (j, i)
                    row   = row_vect[i];
                    coeff = 0.5 * (kt + ktmin1) * arr_vect[i];
                    col   = (j + 1) * ntheta_int + i;
                    val   = -coeff / hs;
                    A_Zebra_Mix_v[smoother][ptr + stencil[Param::top]] += val;
                    col2 = (j - 1) * ntheta_int + i;
                    val2 = -coeff / hsmin1;
                    A_Zebra_Mix_v[smoother][ptr + stencil[Param::bottom]] += val2;
                }
                if (extrapol && smoother == 0) {
                    start_loop = 0;
                    shift_loop = 2;
                    for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                        kt     = thetaplus_per[i + 1];
                        ktmin1 = thetaplus_per[i];

                        // Update (j, i+1)
                        ptr                                                 = ptr_vect[iplus(i, ntheta_int)];
                        row                                                 = row_vect[iplus(i, ntheta_int)];
                        col                                                 = j * ntheta_int + i;
                        coeff                                               = 0.5 * (hs + hsmin1) * att_vect[i];
                        val                                                 = -coeff / kt;
                        A_Zebra_Mix_r[smoother][ptr + stencil[Param::left]] = row;
                        A_Zebra_Mix_c[smoother][ptr + stencil[Param::left]] = col;
                        A_Zebra_Mix_v[smoother][ptr + stencil[Param::left]] += val;

                        // Update (j, i-1)
                        ptr                                                  = ptr_vect[imoins(i, ntheta_int)];
                        row                                                  = row_vect[imoins(i, ntheta_int)];
                        val                                                  = -coeff / ktmin1;
                        A_Zebra_Mix_r[smoother][ptr + stencil[Param::right]] = row;
                        A_Zebra_Mix_c[smoother][ptr + stencil[Param::right]] = col;
                        A_Zebra_Mix_v[smoother][ptr + stencil[Param::right]] += val;

                        // Update (j, i)
                        kt     = thetaplus_per[i + 2];
                        ktmin1 = thetaplus_per[i + 1];
                        ptr    = ptr_vect[iplus(i, ntheta_int)];

                        row    = row_vect[i + 1];
                        coeff2 = 0.5 * (hs + hsmin1) * att_vect[i + 1];
                        col    = j * ntheta_int + imoins(i + 1, ntheta_int);
                        val    = -coeff2 / ktmin1;
                        A_Zebra_Mix_v[smoother][ptr + stencil[Param::left]] += val;
                        col2 = j * ntheta_int + iplus(i + 1, ntheta_int);
                        val2 = -coeff2 / kt;
                        A_Zebra_Mix_v[smoother][ptr + stencil[Param::right]] += val2;
                        // ???
                    }
                }
            }
            else {
                if (extrapol && smoother == 0) {
                    start_loop = 1;
                    shift_loop = 2;
                }
                else {
                    start_loop = 0;
                    shift_loop = 1;
                }
                for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                    kt     = thetaplus_per[i + 1];
                    ktmin1 = thetaplus_per[i];

                    // Update (j-1, i)
                    stencil = stencil_prev;
                    ptr     = ptr_vect_prev[i];

                    row   = row_vect_prev[i];
                    col   = j * ntheta_int + i;
                    coeff = 0.5 * (kt + ktmin1) * arr_vect[i] / hsmin1;
                    val   = -coeff;

                    A_Zebra_Mix_r[smoother][ptr + stencil[Param::top]] = row;
                    A_Zebra_Mix_c[smoother][ptr + stencil[Param::top]] = col;
                    A_Zebra_Mix_v[smoother][ptr + stencil[Param::top]] += val;

                    // Update (j+1, i)
                    stencil = stencil_next;
                    ptr     = ptr_vect_next[i];

                    row                                                   = row_vect_next[i];
                    col                                                   = j * ntheta_int + i;
                    coeff                                                 = 0.5 * (kt + ktmin1) * arr_vect[i] / hs;
                    val                                                   = -coeff;
                    A_Zebra_Mix_r[smoother][ptr + stencil[Param::bottom]] = row;
                    A_Zebra_Mix_c[smoother][ptr + stencil[Param::bottom]] = col;
                    A_Zebra_Mix_v[smoother][ptr + stencil[Param::bottom]] += val;
                }
            }
            if (gyro::icntl[Param::mod_pk] > 0) {
                if (smoother_cur == smoother) {
                    stencil    = stencil_cur;
                    start_loop = 0;
                    if (extrapol && smoother == 0)
                        shift_loop = 2;
                    else
                        shift_loop = 1;
                    for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                        kt     = thetaplus_per[i + 1];
                        ktmin1 = thetaplus_per[i];

                        // Update (j, i+1)
                        ptr                                                        = ptr_vect[iplus(i, ntheta_int)];
                        row                                                        = row_vect[iplus(i, ntheta_int)];
                        col                                                        = (j - 1) * ntheta_int + i;
                        val                                                        = -art_vect[i];
                        A_Zebra_Mix_r[smoother][ptr + stencil[Param::bottom_left]] = row;
                        A_Zebra_Mix_c[smoother][ptr + stencil[Param::bottom_left]] = col;
                        A_Zebra_Mix_v[smoother][ptr + stencil[Param::bottom_left]] += val;
                        col2                                                    = (j + 1) * ntheta_int + i;
                        val2                                                    = art_vect[i];
                        A_Zebra_Mix_r[smoother][ptr + stencil[Param::top_left]] = row;
                        A_Zebra_Mix_c[smoother][ptr + stencil[Param::top_left]] = col2;
                        A_Zebra_Mix_v[smoother][ptr + stencil[Param::top_left]] += val2;

                        // Update (j, i-1)
                        ptr                                                         = ptr_vect[imoins(i, ntheta_int)];
                        row2                                                        = row_vect[imoins(i, ntheta_int)];
                        A_Zebra_Mix_r[smoother][ptr + stencil[Param::bottom_right]] = row2;
                        A_Zebra_Mix_c[smoother][ptr + stencil[Param::bottom_right]] = col;
                        A_Zebra_Mix_v[smoother][ptr + stencil[Param::bottom_right]] -= val;
                        A_Zebra_Mix_r[smoother][ptr + stencil[Param::top_right]] = row2;
                        A_Zebra_Mix_c[smoother][ptr + stencil[Param::top_right]] = col2;
                        A_Zebra_Mix_v[smoother][ptr + stencil[Param::top_right]] -= val2;
                    }
                }
                else {
                    if (extrapol && smoother == 0) {
                        start_loop = 1;
                        shift_loop = 2;
                    }
                    else {
                        start_loop = 0;
                        shift_loop = 1;
                    }
                    for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                        kt     = thetaplus_per[i + 1];
                        ktmin1 = thetaplus_per[i];

                        // Update (j-1, i)
                        stencil = stencil_prev;
                        ptr     = ptr_vect_prev[i];

                        row = row_vect_prev[i];
                        col = j * ntheta_int + imoins(i, ntheta_int);
                        val = art_vect[i];
                        A_Zebra_Mix_r[smoother][ptr + stencil[Param::top_left]] = row;
                        A_Zebra_Mix_c[smoother][ptr + stencil[Param::top_left]] = col;
                        A_Zebra_Mix_v[smoother][ptr + stencil[Param::top_left]] += val;
                        col2 = j * ntheta_int + iplus(i, ntheta_int);
                        val2 = -art_vect[i];
                        A_Zebra_Mix_r[smoother][ptr + stencil[Param::top_right]] = row;
                        A_Zebra_Mix_c[smoother][ptr + stencil[Param::top_right]] = col2;
                        A_Zebra_Mix_v[smoother][ptr + stencil[Param::top_right]] += val2;

                        // Update (j+1, i)
                        stencil = stencil_next;
                        ptr     = ptr_vect_next[i];

                        row                                                        = row_vect_next[i];
                        val                                                        = -art_vect[i];
                        val2                                                       = art_vect[i];
                        A_Zebra_Mix_r[smoother][ptr + stencil[Param::bottom_left]] = row;
                        A_Zebra_Mix_c[smoother][ptr + stencil[Param::bottom_left]] = col;
                        A_Zebra_Mix_v[smoother][ptr + stencil[Param::bottom_left]] += val;

                        A_Zebra_Mix_r[smoother][ptr + stencil[Param::bottom_right]] = row;
                        A_Zebra_Mix_c[smoother][ptr + stencil[Param::bottom_right]] = col2;
                        A_Zebra_Mix_v[smoother][ptr + stencil[Param::bottom_right]] += val2;
                    }
                }
            }
        }

        /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 * !!!!!!!!!!!    First lines    !!!!!!!!!!!
                 * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 */
        j     = start_j;
        odd_j = j % 2;

        smoother_cur = j;
        ptr_vect     = get_ptr_sc(j, smoother_cur, 1);
        row_vect     = get_row(j, smoother_cur, extrapol, 0, 1);
        if (j < delete_circles - 1) {
            smoother_next = (j + 1) % 2;
            row_vect_next = get_row(j + 1, smoother_next, extrapol, 0, 1);
        }
        else if (j == delete_circles - 1) {
            smoother_vect_next = get_smoother_radial(j + 1);
            smoother_next      = smoother_vect_next[0];
            row_vect_next      = get_row(j + 1, 2, extrapol, 0, 1);
        }
        ptr_vect_next = get_ptr_sc(j + 1, smoother_next, 1);
        stencil_cur   = get_stencil_sc(j, smoother_cur, 1);
        stencil_next  = get_stencil_sc(j + 1, smoother_next, 1);

        hs = hplus[j];
        if (j > 0) {
            hsmin1 = hplus[j - 1];
        }
        else {
            hsmin1 = 2 * r[0]; // across the origin
        }
        // Across and DB_int updates
        gyro::arr_att_art(r[j], theta, sin_theta, cos_theta, ntheta_int, arr_vect, att_vect, art_vect, 0);
        if (smoother_cur == smoother) {
            stencil = stencil_cur;
            if (extrapol && smoother == 0) {
                start_loop = 1;
                shift_loop = 2;
            }
            else {
                start_loop = 0;
                shift_loop = 1;
            }
            for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                kt     = thetaplus_per[i + 1];
                ktmin1 = thetaplus_per[i];
                ptr    = ptr_vect[i];

                // Update (j, i)
                row   = row_vect[i];
                coeff = 0.5 * (kt + ktmin1) * arr_vect[i];
                col   = (j + 1) * ntheta_int + i;
                val   = -coeff / hs;

                A_Zebra_Mix_v[smoother][ptr + stencil[Param::top]] += val;
            }
            if (extrapol && smoother == 0) {
                start_loop = 0;
                shift_loop = 2;
                for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                    kt     = thetaplus_per[i + 1];
                    ktmin1 = thetaplus_per[i];

                    // Update (j, i)
                    ptr   = ptr_vect[iplus(i, ntheta_int)];
                    row   = row_vect[iplus(i, ntheta_int)];
                    col   = j * ntheta_int + i;
                    coeff = 0.5 * (hs + hsmin1) * att_vect[i];
                    val   = -coeff / kt;
                    A_Zebra_Mix_v[smoother][ptr + stencil[Param::right]] += val;

                    // Update (j, i-1)
                    ptr = ptr_vect[imoins(i, ntheta_int)];
                    row = row_vect[imoins(i, ntheta_int)];
                    val = -coeff / ktmin1;
                    A_Zebra_Mix_v[smoother][ptr + stencil[Param::left]] += val;

                    // Update (j, i)
                    kt     = thetaplus_per[i + 2];
                    ktmin1 = thetaplus_per[i + 1];
                    ptr    = ptr_vect[i + 1];

                    row                                                  = row_vect[i + 1];
                    coeff2                                               = 0.5 * (hs + hsmin1) * att_vect[i + 1];
                    col                                                  = j * ntheta_int + imoins(i + 1, ntheta_int);
                    val                                                  = -coeff2 / ktmin1;
                    A_Zebra_Mix_r[smoother][ptr + stencil[Param::right]] = row;
                    A_Zebra_Mix_c[smoother][ptr + stencil[Param::right]] = col;
                    A_Zebra_Mix_v[smoother][ptr + stencil[Param::right]] += val;
                    col2                                                = j * ntheta_int + iplus(i + 1, ntheta_int);
                    val2                                                = -coeff2 / kt;
                    A_Zebra_Mix_r[smoother][ptr + stencil[Param::left]] = row;
                    A_Zebra_Mix_c[smoother][ptr + stencil[Param::left]] = col2;
                    A_Zebra_Mix_v[smoother][ptr + stencil[Param::left]] += val2;
                }
            }
            if (gyro::icntl[Param::mod_pk] > 0) {
                start_loop = 0;
                if (extrapol && smoother == 0)
                    shift_loop = 2;
                else
                    shift_loop = 1;
                for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                    kt     = thetaplus_per[i + 1];
                    ktmin1 = thetaplus_per[i];

                    // Update (j, i+1)
                    ptr = ptr_vect[iplus(i, ntheta_int)];
                    row = row_vect[iplus(i, ntheta_int)];
                    col = (j + 1) * ntheta_int + i;
                    val = art_vect[i];

                    A_Zebra_Mix_r[smoother][ptr + stencil[Param::top_left]] = row;
                    A_Zebra_Mix_c[smoother][ptr + stencil[Param::top_left]] = col;
                    A_Zebra_Mix_v[smoother][ptr + stencil[Param::top_left]] += val;

                    // Update (j, i-1)
                    ptr = ptr_vect[imoins(i, ntheta_int)];
                    row = row_vect[imoins(i, ntheta_int)];

                    A_Zebra_Mix_r[smoother][ptr + stencil[Param::top_right]] = row;
                    A_Zebra_Mix_c[smoother][ptr + stencil[Param::top_right]] = col;
                    A_Zebra_Mix_v[smoother][ptr + stencil[Param::top_right]] -= val;
                }
            }
        }

        // Update (j+1, i)
        if (j < delete_circles - 1 && smoother_next == smoother) {
            stencil = get_stencil_sc(j + 1, smoother, 1);
            if (extrapol && smoother == 0) {
                start_loop = 1;
                shift_loop = 2;
            }
            else {
                start_loop = 0;
                shift_loop = 1;
            }
            for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                kt     = thetaplus_per[i + 1];
                ktmin1 = thetaplus_per[i];
                ptr    = ptr_vect_next[i];

                row   = row_vect_next[i];
                col   = j * ntheta_int + i;
                coeff = 0.5 * (kt + ktmin1) * arr_vect[i] / hs;

                val = -coeff;

                A_Zebra_Mix_r[smoother][ptr + stencil[Param::bottom]] = row;
                A_Zebra_Mix_c[smoother][ptr + stencil[Param::bottom]] = col;
                A_Zebra_Mix_v[smoother][ptr + stencil[Param::bottom]] += val;
            }
            if (gyro::icntl[Param::mod_pk] > 0)
                for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                    ptr = ptr_vect_next[i];

                    row                                                        = row_vect_next[i];
                    col                                                        = j * ntheta_int + imoins(i, ntheta_int);
                    val                                                        = -art_vect[i];
                    A_Zebra_Mix_r[smoother][ptr + stencil[Param::bottom_left]] = row;
                    A_Zebra_Mix_c[smoother][ptr + stencil[Param::bottom_left]] = col;
                    A_Zebra_Mix_v[smoother][ptr + stencil[Param::bottom_left]] += val;

                    col2                                                        = j * ntheta_int + iplus(i, ntheta_int);
                    val2                                                        = art_vect[i];
                    A_Zebra_Mix_r[smoother][ptr + stencil[Param::bottom_right]] = row;
                    A_Zebra_Mix_c[smoother][ptr + stencil[Param::bottom_right]] = col2;
                    A_Zebra_Mix_v[smoother][ptr + stencil[Param::bottom_right]] += val2;
                }
        }
    }

    /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! RADIAL   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         */
    if (smoother > 1) {
        int diff  = ntheta_int % 3;
        int odd_d = delete_circles % 2;

        for (int i = 0; i < ntheta_int; i++) {
            {
                int im        = imoins(i, ntheta_int);
                int ip        = iplus(i, ntheta_int);
                int im2       = imoins(im, ntheta_int);
                int ip2       = iplus(ip, ntheta_int);
                int odd_i     = i % 2;
                int ind_m_odd = i, ind_p_odd = i;
                if (odd_i) {
                    ind_m_odd = im;
                    ind_p_odd = ip;
                }

                int s_cur     = (i % 2) + 2;
                int s_next    = ((i + 1) % 2) + 2;
                kt            = thetaplus_per[i + 1];
                ktmin1        = thetaplus_per[i];
                smoother_cur  = s_cur;
                smoother_prev = s_next;
                smoother_next = s_next;
                row_vect_prev = get_row_i_glob(nr, im, smoother_prev, extrapol);
                row_vect      = get_row_i_glob(nr, i, smoother_cur, extrapol);
                row_vect_next = get_row_i_glob(nr, ip, smoother_next, extrapol);
                gyro::arr_att_art(r, theta[i], arr_vect, att_vect, art_vect, 0);

                /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                         * !!!!!!!!!!! Link Circle-Radial (radial) !!!!!!!!!!
                         * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                         */
                j                        = delete_circles - 1;
                int smoother_prev_circle = get_smoother(i, j - 1);
                int smoother_cur_left    = get_smoother(im, j);
                int smoother_cur_circle  = get_smoother(i, j);
                int smoother_cur_right   = get_smoother(ip, j);
                int smoother_next_circle = get_smoother(i, j + 1);
                ptr_vect_prev            = get_ptr_sc(j - 1, smoother_prev_circle, 1);
                ptr_vect                 = get_ptr_sc(j, smoother_cur_circle, 1);
                ptr_vect_next            = get_ptr_sc(j + 1, smoother_next_circle, 1);
                stencil_prev             = get_stencil_sc(j - 1, smoother_prev_circle, 1);
                stencil_cur_left         = get_stencil_sc(j, smoother_cur_left, 1);
                stencil_cur              = get_stencil_sc(j, smoother_cur_circle, 1);
                stencil_cur_right        = get_stencil_sc(j, smoother_cur_right, 1);
                stencil_next             = get_stencil_sc(j + 1, smoother_next_circle, 1);
                if (!(extrapol && smoother == 2 && j % 2 == 1)) {
                    hs     = hplus[j];
                    hsmin1 = hplus[j - 1];
                    if ((smoother == 3 && i % 2 == 1) || (smoother == 2 && i % 2 == 0)) {
                        stencil                                               = stencil_next;
                        ptr                                                   = ptr_vect_next[i];
                        row                                                   = row_vect[j + 1];
                        col                                                   = j * ntheta_int + i;
                        coeff                                                 = 0.5 * (kt + ktmin1) * arr_vect[j] / hs;
                        val                                                   = -coeff;
                        A_Zebra_Mix_r[smoother][ptr + stencil[Param::bottom]] = row;
                        A_Zebra_Mix_c[smoother][ptr + stencil[Param::bottom]] = col;
                        A_Zebra_Mix_v[smoother][ptr + stencil[Param::bottom]] += -coeff;

                        if (gyro::icntl[Param::mod_pk] > 0) {
                            stencil                                                    = stencil_next;
                            ptr                                                        = ptr_vect_next[i];
                            row                                                        = row_vect[j + 1];
                            col                                                        = j * ntheta_int + im;
                            val                                                        = -art_vect[j];
                            A_Zebra_Mix_r[smoother][ptr + stencil[Param::bottom_left]] = row;
                            A_Zebra_Mix_c[smoother][ptr + stencil[Param::bottom_left]] = col;
                            A_Zebra_Mix_v[smoother][ptr + stencil[Param::bottom_left]] += val;
                            col2                                                        = j * ntheta_int + ip;
                            val2                                                        = art_vect[j];
                            A_Zebra_Mix_r[smoother][ptr + stencil[Param::bottom_right]] = row;
                            A_Zebra_Mix_c[smoother][ptr + stencil[Param::bottom_right]] = col2;
                            A_Zebra_Mix_v[smoother][ptr + stencil[Param::bottom_right]] += val2;
                        }
                    }
                }

                /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             * !!!!!!!!!!! Link Radial-Circle (radial) !!!!!!!!!!
             * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             */

                j                    = delete_circles;
                smoother_prev_circle = get_smoother(i, j - 1);
                smoother_cur_left    = get_smoother(im, j);
                smoother_cur_circle  = get_smoother(i, j);
                smoother_cur_right   = get_smoother(ip, j);
                smoother_next_circle = get_smoother(i, j + 1);
                ptr_vect_prev        = get_ptr_sc(j - 1, smoother_prev_circle, 1);
                ptr_vect             = get_ptr_sc(j, smoother_cur_circle, 1);
                ptr_vect_next        = get_ptr_sc(j + 1, smoother_next_circle, 1);
                stencil_prev         = get_stencil_sc(j - 1, smoother_prev_circle, 1);
                stencil_cur_left     = get_stencil_sc(j, smoother_cur_left, 1);
                stencil_cur          = get_stencil_sc(j, smoother_cur_circle, 1);
                stencil_cur_right    = get_stencil_sc(j, smoother_cur_right, 1);
                stencil_next         = get_stencil_sc(j + 1, smoother_next_circle, 1);
                hs                   = hplus[j];
                hsmin1               = hplus[j - 1];
                if (!(extrapol && smoother == 2 && j % 2 == 0)) {
                    if ((smoother == 3 && i % 2 == 1) || (smoother == 2 && i % 2 == 0)) {
                        // Update (j, i)
                        stencil = stencil_cur;
                        ptr     = ptr_vect[i];
                        row     = row_vect[j];
                        coeff   = 0.5 * (kt + ktmin1) * arr_vect[j];
                        col     = (j - 1) * ntheta_int + i;
                        val     = -coeff / hsmin1;
                        A_Zebra_Mix_v[smoother][ptr + stencil[Param::bottom]] += val;

                        coeff2 = 0.5 * (hs + hsmin1) * att_vect[j];
                        col    = j * ntheta_int + im;
                        val    = -coeff2 / ktmin1;
                        A_Zebra_Mix_v[smoother][ptr + stencil[Param::left]] += val;
                        col2 = j * ntheta_int + ip;
                        val2 = -coeff2 / kt;
                        A_Zebra_Mix_v[smoother][ptr + stencil[Param::right]] += val2;
                    }
                    if ((smoother == 3 && i % 2 == 0) || (smoother == 2 && i % 2 == 1)) {
                        stencil = stencil_cur;
                        // Update (j, i+1)
                        stencil                                             = stencil_cur_right;
                        ptr                                                 = ptr_vect[ip];
                        row                                                 = row_vect_next[j];
                        col                                                 = j * ntheta_int + i;
                        coeff                                               = 0.5 * (hs + hsmin1) * att_vect[j] / kt;
                        val                                                 = -coeff;
                        A_Zebra_Mix_r[smoother][ptr + stencil[Param::left]] = row;
                        A_Zebra_Mix_c[smoother][ptr + stencil[Param::left]] = col;
                        A_Zebra_Mix_v[smoother][ptr + stencil[Param::left]] += val;

                        // Update (j, i-1)
                        stencil = stencil_cur_left;
                        ptr     = ptr_vect[im];
                        row     = row_vect_prev[j];
                        col     = j * ntheta_int + i;
                        coeff   = 0.5 * (hs + hsmin1) * att_vect[j] / ktmin1;
                        val     = -coeff;
                        A_Zebra_Mix_r[smoother][ptr + stencil[Param::right]] = row;
                        A_Zebra_Mix_c[smoother][ptr + stencil[Param::right]] = col;
                        A_Zebra_Mix_v[smoother][ptr + stencil[Param::right]] += -coeff;
                    }
                }
                // Update (j+1, i) (Not in DB_ext)
                if (extrapol && smoother == 2 && i % 2 == 0) {
                    if (j % 2 == 1) {
                        stencil                                            = stencil_cur;
                        ptr                                                = ptr_vect[i];
                        row                                                = row_vect[j];
                        coeff                                              = 0.5 * (kt + ktmin1) * arr_vect[j];
                        col                                                = (j + 1) * ntheta_int + i;
                        val                                                = -coeff / hs;
                        A_Zebra_Mix_r[smoother][ptr + stencil[Param::top]] = row;
                        A_Zebra_Mix_c[smoother][ptr + stencil[Param::top]] = col;
                        A_Zebra_Mix_v[smoother][ptr + stencil[Param::top]] += val;
                    }
                    else {
                        stencil                                               = stencil_next;
                        ptr                                                   = ptr_vect_next[i];
                        row                                                   = row_vect[j + 1];
                        col                                                   = j * ntheta_int + i;
                        coeff                                                 = 0.5 * (kt + ktmin1) * arr_vect[j] / hs;
                        val                                                   = -coeff;
                        A_Zebra_Mix_r[smoother][ptr + stencil[Param::bottom]] = row;
                        A_Zebra_Mix_c[smoother][ptr + stencil[Param::bottom]] = col;
                        A_Zebra_Mix_v[smoother][ptr + stencil[Param::bottom]] += val;
                    }
                }
                if (gyro::icntl[Param::mod_pk] > 0) {
                    if (!(extrapol && smoother == 2 && j % 2 == 0)) {
                        if ((smoother == 3 && i % 2 == 0) || (smoother == 2 && i % 2 == 1)) {
                            stencil                                                    = stencil_cur_right;
                            ptr                                                        = ptr_vect[ip];
                            row                                                        = row_vect_next[j];
                            col                                                        = (j - 1) * ntheta_int + i;
                            val                                                        = -art_vect[j];
                            A_Zebra_Mix_r[smoother][ptr + stencil[Param::bottom_left]] = row;
                            A_Zebra_Mix_c[smoother][ptr + stencil[Param::bottom_left]] = col;
                            A_Zebra_Mix_v[smoother][ptr + stencil[Param::bottom_left]] += val;
                            col2                                                    = (j + 1) * ntheta_int + i;
                            val2                                                    = art_vect[j];
                            A_Zebra_Mix_r[smoother][ptr + stencil[Param::top_left]] = row;
                            A_Zebra_Mix_c[smoother][ptr + stencil[Param::top_left]] = col2;
                            A_Zebra_Mix_v[smoother][ptr + stencil[Param::top_left]] += val2;

                            // Update (j, i-1)
                            stencil                                                     = stencil_cur_left;
                            ptr                                                         = ptr_vect[im];
                            row2                                                        = row_vect_prev[j];
                            A_Zebra_Mix_r[smoother][ptr + stencil[Param::bottom_right]] = row2;
                            A_Zebra_Mix_c[smoother][ptr + stencil[Param::bottom_right]] = col;
                            A_Zebra_Mix_v[smoother][ptr + stencil[Param::bottom_right]] -= val;

                            A_Zebra_Mix_r[smoother][ptr + stencil[Param::top_right]] = row2;
                            A_Zebra_Mix_c[smoother][ptr + stencil[Param::top_right]] = col2;
                            A_Zebra_Mix_v[smoother][ptr + stencil[Param::top_right]] -= val2;
                        }
                    }
                    if (!(extrapol && smoother == 2 && j % 2 == 1)) {
                        if ((smoother == 3 && i % 2 == 1) || (smoother == 2 && i % 2 == 0)) {
                            // Update (j+1, i)
                            stencil                                                    = stencil_next;
                            ptr                                                        = ptr_vect_next[i];
                            row                                                        = row_vect[j + 1];
                            col                                                        = j * ntheta_int + im;
                            val                                                        = -art_vect[j];
                            A_Zebra_Mix_r[smoother][ptr + stencil[Param::bottom_left]] = row;
                            A_Zebra_Mix_c[smoother][ptr + stencil[Param::bottom_left]] = col;
                            A_Zebra_Mix_v[smoother][ptr + stencil[Param::bottom_left]] += val;
                            col2                                                        = j * ntheta_int + ip;
                            val2                                                        = art_vect[j];
                            A_Zebra_Mix_r[smoother][ptr + stencil[Param::bottom_right]] = row;
                            A_Zebra_Mix_c[smoother][ptr + stencil[Param::bottom_right]] = col2;
                            A_Zebra_Mix_v[smoother][ptr + stencil[Param::bottom_right]] += val2;
                        }
                    }
                }

                /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!  Interior nodes   !!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     */
                for (j = delete_circles + 1; j < nr_int - 1; j++) {
                    smoother_prev_circle = get_smoother(i, j - 1);
                    smoother_cur_left    = get_smoother(im, j);
                    smoother_cur_circle  = get_smoother(i, j);
                    smoother_cur_right   = get_smoother(ip, j);
                    smoother_next_circle = get_smoother(i, j + 1);
                    ptr_vect_prev        = get_ptr_sc(j - 1, smoother_prev_circle, 1);
                    ptr_vect             = get_ptr_sc(j, smoother_cur_circle, 1);
                    ptr_vect_next        = get_ptr_sc(j + 1, smoother_next_circle, 1);
                    stencil_prev         = get_stencil_sc(j - 1, smoother_prev_circle, 1);
                    stencil_cur_left     = get_stencil_sc(j, smoother_cur_left, 1);
                    stencil_cur          = get_stencil_sc(j, smoother_cur_circle, 1);
                    stencil_cur_right    = get_stencil_sc(j, smoother_cur_right, 1);
                    stencil_next         = get_stencil_sc(j + 1, smoother_next_circle, 1);
                    hs                   = hplus[j];
                    hsmin1               = hplus[j - 1];
                    if (smoother > 1 && !(extrapol && smoother == 2 && j % 2 == 0)) {
                        if ((smoother == 3 && i % 2 == 1) || (smoother == 2 && i % 2 == 0)) {
                            // Update (j, i)
                            stencil = stencil_cur;
                            ptr     = ptr_vect[i];
                            row     = row_vect[j];
                            // coeff  = 0.5 * (kt + ktmin1) * arr_vect[j];
                            coeff2 = 0.5 * (hs + hsmin1) * att_vect[j];
                            col    = j * ntheta_int + im;
                            val    = -coeff2 / ktmin1;
                            A_Zebra_Mix_v[smoother][ptr + stencil[Param::left]] += val;
                            col2 = j * ntheta_int + ip;
                            val2 = -coeff2 / kt;
                            A_Zebra_Mix_v[smoother][ptr + stencil[Param::right]] += val2;
                        }
                        if ((smoother == 3 && i % 2 == 0) || (smoother == 2 && i % 2 == 1)) {
                            // Update (j, i+1)
                            stencil                                             = stencil_cur_right;
                            ptr                                                 = ptr_vect[ip];
                            row                                                 = row_vect_next[j];
                            col                                                 = j * ntheta_int + i;
                            coeff                                               = 0.5 * (hs + hsmin1) * att_vect[j];
                            val                                                 = -coeff / kt;
                            A_Zebra_Mix_r[smoother][ptr + stencil[Param::left]] = row;
                            A_Zebra_Mix_c[smoother][ptr + stencil[Param::left]] = col;
                            A_Zebra_Mix_v[smoother][ptr + stencil[Param::left]] += val;

                            // Update (j, i-1)
                            stencil                                              = stencil_cur_left;
                            ptr                                                  = ptr_vect[im];
                            row                                                  = row_vect_prev[j];
                            col                                                  = j * ntheta_int + i;
                            coeff                                                = 0.5 * (hs + hsmin1) * att_vect[j];
                            val                                                  = -coeff / ktmin1;
                            A_Zebra_Mix_r[smoother][ptr + stencil[Param::right]] = row;
                            A_Zebra_Mix_c[smoother][ptr + stencil[Param::right]] = col;
                            A_Zebra_Mix_v[smoother][ptr + stencil[Param::right]] += val;
                        }
                    }
                    if (extrapol && smoother == 2 && i % 2 == 0) {
                        if (j % 2 == 1) {
                            // Update (j, i)
                            stencil = stencil_cur;
                            ptr     = ptr_vect[i];
                            row     = row_vect[j];
                            coeff   = 0.5 * (kt + ktmin1) * arr_vect[j];
                            col     = (j + 1) * ntheta_int + i;
                            val     = -coeff / hs;
                            A_Zebra_Mix_v[smoother][ptr + stencil[Param::top]] += val;
                            col2 = (j - 1) * ntheta_int + i;
                            val2 = -coeff / hsmin1;
                            A_Zebra_Mix_v[smoother][ptr + stencil[Param::bottom]] += val2;
                        }
                        else {

                            // Update (j-1, i)
                            stencil = stencil_prev;
                            ptr     = ptr_vect_prev[i];
                            row     = row_vect[j - 1];
                            col     = j * ntheta_int + i;
                            coeff   = 0.5 * (kt + ktmin1) * arr_vect[j] / hsmin1;
                            val     = -coeff;
                            A_Zebra_Mix_r[smoother][ptr + stencil[Param::top]] = row;
                            A_Zebra_Mix_c[smoother][ptr + stencil[Param::top]] = col;
                            A_Zebra_Mix_v[smoother][ptr + stencil[Param::top]] += val;

                            // Update (j+1, i)
                            stencil = stencil_next;
                            ptr     = ptr_vect_next[i];
                            row     = row_vect[j + 1];
                            col     = j * ntheta_int + i;
                            coeff   = 0.5 * (kt + ktmin1) * arr_vect[j] / hs;
                            val     = -coeff;
                            A_Zebra_Mix_r[smoother][ptr + stencil[Param::bottom]] = row;
                            A_Zebra_Mix_c[smoother][ptr + stencil[Param::bottom]] = col;
                            A_Zebra_Mix_v[smoother][ptr + stencil[Param::bottom]] += val;
                        }
                    }
                    if (gyro::icntl[Param::mod_pk] > 0) {
                        if (smoother > 1 && !(extrapol && smoother == 2 && j % 2 == 1)) {
                            if ((smoother == 3 && i % 2 == 1) || (smoother == 2 && i % 2 == 0)) {
                                // Update (j-1, i)
                                stencil                                                 = stencil_prev;
                                ptr                                                     = ptr_vect_prev[i];
                                row                                                     = row_vect[j - 1];
                                col                                                     = j * ntheta_int + im;
                                val                                                     = art_vect[j];
                                A_Zebra_Mix_r[smoother][ptr + stencil[Param::top_left]] = row;
                                A_Zebra_Mix_c[smoother][ptr + stencil[Param::top_left]] = col;
                                A_Zebra_Mix_v[smoother][ptr + stencil[Param::top_left]] += val;
                                col2                                                     = j * ntheta_int + ip;
                                val2                                                     = -art_vect[j];
                                A_Zebra_Mix_r[smoother][ptr + stencil[Param::top_right]] = row;
                                A_Zebra_Mix_c[smoother][ptr + stencil[Param::top_right]] = col2;
                                A_Zebra_Mix_v[smoother][ptr + stencil[Param::top_right]] += val2;

                                // Update (j+1, i)
                                stencil                                                    = stencil_next;
                                ptr                                                        = ptr_vect_next[i];
                                row2                                                       = row_vect[j + 1];
                                A_Zebra_Mix_r[smoother][ptr + stencil[Param::bottom_left]] = row2;
                                A_Zebra_Mix_c[smoother][ptr + stencil[Param::bottom_left]] = col;
                                A_Zebra_Mix_v[smoother][ptr + stencil[Param::bottom_left]] -= val;

                                A_Zebra_Mix_r[smoother][ptr + stencil[Param::bottom_right]] = row2;
                                A_Zebra_Mix_c[smoother][ptr + stencil[Param::bottom_right]] = col2;
                                A_Zebra_Mix_v[smoother][ptr + stencil[Param::bottom_right]] -= val2;
                            }
                        }
                        if (smoother > 1 && !(extrapol && smoother == 2 && j % 2 == 0)) {
                            if ((smoother == 3 && i % 2 == 0) || (smoother == 2 && i % 2 == 1)) {
                                // Update (j, i+1)
                                stencil                                                    = stencil_cur_right;
                                ptr                                                        = ptr_vect[ip];
                                row                                                        = row_vect_next[j];
                                col                                                        = (j - 1) * ntheta_int + i;
                                val                                                        = -art_vect[j];
                                A_Zebra_Mix_r[smoother][ptr + stencil[Param::bottom_left]] = row;
                                A_Zebra_Mix_c[smoother][ptr + stencil[Param::bottom_left]] = col;
                                A_Zebra_Mix_v[smoother][ptr + stencil[Param::bottom_left]] += val;
                                col2                                                    = (j + 1) * ntheta_int + i;
                                val2                                                    = art_vect[j];
                                A_Zebra_Mix_r[smoother][ptr + stencil[Param::top_left]] = row;
                                A_Zebra_Mix_c[smoother][ptr + stencil[Param::top_left]] = col2;
                                A_Zebra_Mix_v[smoother][ptr + stencil[Param::top_left]] += val2;

                                // Update (j, i-1)
                                stencil                                                     = stencil_cur_left;
                                ptr                                                         = ptr_vect[im];
                                row                                                         = row_vect_prev[j];
                                col                                                         = (j - 1) * ntheta_int + i;
                                val                                                         = -art_vect[j];
                                A_Zebra_Mix_r[smoother][ptr + stencil[Param::bottom_right]] = row;
                                A_Zebra_Mix_c[smoother][ptr + stencil[Param::bottom_right]] = col;
                                A_Zebra_Mix_v[smoother][ptr + stencil[Param::bottom_right]] -= val;
                                col2                                                     = (j + 1) * ntheta_int + i;
                                val2                                                     = art_vect[j];
                                A_Zebra_Mix_r[smoother][ptr + stencil[Param::top_right]] = row;
                                A_Zebra_Mix_c[smoother][ptr + stencil[Param::top_right]] = col2;
                                A_Zebra_Mix_v[smoother][ptr + stencil[Param::top_right]] -= val2;
                            }
                        }
                    }
                }
                /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!        Last lines       !!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     */
                // DB_ext updates (~~~ Interior - (j+1, i) + DB)
                j                    = nr_int - 1;
                smoother_prev_circle = get_smoother(i, j - 1);
                smoother_cur_left    = get_smoother(im, j);
                smoother_cur_circle  = get_smoother(i, j);
                smoother_cur_right   = get_smoother(ip, j);
                smoother_next_circle = get_smoother(i, j + 1);
                ptr_vect_prev        = get_ptr_sc(j - 1, smoother_prev_circle, 1);
                ptr_vect             = get_ptr_sc(j, smoother_cur_circle, 1);
                ptr_vect_next        = get_ptr_sc(j + 1, smoother_next_circle, 1);
                stencil_prev         = get_stencil_sc(j - 1, smoother_prev_circle, 1);
                stencil_cur_left     = get_stencil_sc(j, smoother_cur_left, 1);
                stencil_cur          = get_stencil_sc(j, smoother_cur_circle, 1);
                stencil_cur_right    = get_stencil_sc(j, smoother_cur_right, 1);
                stencil_next         = get_stencil_sc(j + 1, smoother_next_circle, 1);
                hs                   = hplus[j];
                hsmin1               = hplus[j - 1];
                if (smoother > 1 && !(extrapol && smoother == 2 && j % 2 == 0)) {
                    // Update (j, i)
                    if ((smoother == 3 && i % 2 == 1) || (smoother == 2 && i % 2 == 0)) {
                        stencil = stencil_cur;
                        ptr     = ptr_vect[i];
                        row     = row_vect[j];
                        // Contribution to middle (top) from DB
                        coeff3 = 0.5 * (hs + hsmin1) * att_vect[j];
                        col    = j * ntheta_int + im;
                        val    = -coeff3 / ktmin1;
                        A_Zebra_Mix_v[smoother][ptr + stencil[Param::left]] += val;
                        col2 = j * ntheta_int + ip;
                        val2 = -coeff3 / kt;
                        A_Zebra_Mix_v[smoother][ptr + stencil[Param::right]] += val2;
                    }
                    if ((smoother == 3 && i % 2 == 0) || (smoother == 2 && i % 2 == 1)) {
                        stencil = stencil_cur;
                        // Update (j, i+1)
                        stencil                                             = stencil_cur_right;
                        ptr                                                 = ptr_vect[ip];
                        row                                                 = row_vect_next[j];
                        col                                                 = j * ntheta_int + i;
                        coeff                                               = 0.5 * (hs + hsmin1) * att_vect[j];
                        val                                                 = -coeff / kt;
                        A_Zebra_Mix_r[smoother][ptr + stencil[Param::left]] = row;
                        A_Zebra_Mix_c[smoother][ptr + stencil[Param::left]] = col;
                        A_Zebra_Mix_v[smoother][ptr + stencil[Param::left]] += val;

                        // Update (j, i-1)
                        stencil                                              = stencil_cur_left;
                        ptr                                                  = ptr_vect[im];
                        row                                                  = row_vect_prev[j];
                        col                                                  = j * ntheta_int + i;
                        coeff                                                = 0.5 * (hs + hsmin1) * att_vect[j];
                        val                                                  = -coeff / ktmin1;
                        A_Zebra_Mix_r[smoother][ptr + stencil[Param::right]] = row;
                        A_Zebra_Mix_c[smoother][ptr + stencil[Param::right]] = col;
                        A_Zebra_Mix_v[smoother][ptr + stencil[Param::right]] += val;
                    }
                }
                if (extrapol && smoother == 2 && i % 2 == 0) {
                    if (j % 2 == 0) {
                        // Update (j-1, i)
                        stencil                                            = stencil_prev;
                        ptr                                                = ptr_vect_prev[i];
                        row                                                = row_vect[j - 1];
                        col                                                = j * ntheta_int + i;
                        coeff                                              = 0.5 * (kt + ktmin1) * arr_vect[j] / hsmin1;
                        val                                                = -coeff;
                        A_Zebra_Mix_r[smoother][ptr + stencil[Param::top]] = row;
                        A_Zebra_Mix_c[smoother][ptr + stencil[Param::top]] = col;
                        A_Zebra_Mix_v[smoother][ptr + stencil[Param::top]] += val;
                    }
                    else {
                        // Update (j, i)
                        stencil = stencil_cur;
                        ptr     = ptr_vect[i];
                        row     = row_vect[j];
                        coeff   = 0.5 * (kt + ktmin1) * arr_vect[j];
                        col     = (j - 1) * ntheta_int + i;
                        val     = -coeff / hsmin1;
                        A_Zebra_Mix_v[smoother][ptr + stencil[Param::bottom]] += val;
                    }
                }
                if (gyro::icntl[Param::mod_pk] > 0) {
                    if (smoother > 1 && !(extrapol && smoother == 2 && j % 2 == 1)) {
                        if ((smoother == 3 && i % 2 == 1) || (smoother == 2 && i % 2 == 0)) {
                            // Update (j-1, i)
                            stencil                                                 = stencil_prev;
                            ptr                                                     = ptr_vect_prev[i];
                            row                                                     = row_vect[j - 1];
                            col                                                     = j * ntheta_int + im;
                            val                                                     = art_vect[j];
                            A_Zebra_Mix_r[smoother][ptr + stencil[Param::top_left]] = row;
                            A_Zebra_Mix_c[smoother][ptr + stencil[Param::top_left]] = col;
                            A_Zebra_Mix_v[smoother][ptr + stencil[Param::top_left]] += val;
                            col2                                                     = j * ntheta_int + ip;
                            val2                                                     = -art_vect[j];
                            A_Zebra_Mix_r[smoother][ptr + stencil[Param::top_right]] = row;
                            A_Zebra_Mix_c[smoother][ptr + stencil[Param::top_right]] = col2;
                            A_Zebra_Mix_v[smoother][ptr + stencil[Param::top_right]] += val2;
                        }
                    }
                    if (smoother > 1 && !(extrapol && smoother == 2 && j % 2 == 0)) {
                        if ((smoother == 3 && i % 2 == 0) || (smoother == 2 && i % 2 == 1)) {
                            stencil = stencil_cur;
                            // Update (j, i+1)
                            stencil                                                    = stencil_cur_right;
                            ptr                                                        = ptr_vect[ip];
                            row                                                        = row_vect_next[j];
                            col                                                        = (j - 1) * ntheta_int + i;
                            val                                                        = -art_vect[j];
                            A_Zebra_Mix_r[smoother][ptr + stencil[Param::bottom_left]] = row;
                            A_Zebra_Mix_c[smoother][ptr + stencil[Param::bottom_left]] = col;
                            A_Zebra_Mix_v[smoother][ptr + stencil[Param::bottom_left]] += val;

                            // Update (j, i-1)
                            stencil                                                     = stencil_cur_left;
                            ptr                                                         = ptr_vect[im];
                            row                                                         = row_vect_prev[j];
                            col                                                         = (j - 1) * ntheta_int + i;
                            val                                                         = -art_vect[j];
                            A_Zebra_Mix_r[smoother][ptr + stencil[Param::bottom_right]] = row;
                            A_Zebra_Mix_c[smoother][ptr + stencil[Param::bottom_right]] = col;
                            A_Zebra_Mix_v[smoother][ptr + stencil[Param::bottom_right]] -= val;
                        }
                    }
                }
            }
        }
    }
} /* ----- end of level::build_Asc_ortho ----- */

/*! \brief Applies the matrix A_sc_ortho for a smoother explicitely on the level l
 *
 * Asc_ortho corresponds to all lines of a smoother s with colour c and the columns not in Asc
 * 
 * \param smoother_todo: the smoother of this Asc_ortho matrix
*/
void level::apply_Asc_ortho(std::vector<double>& Au, std::vector<double>& u, int smoother_todo, int v, int c,
                            int* dep_Asc_cur, int* dep_Asc_prev, int* dep_Asc1, int* dep_Asc_ortho_cur)
{
    int start_j;
    int extrapol = gyro::icntl[Param::extrapolation] == 1 && l == 0;

    int smoother = smoother_todo;
    int base_prec;
    if (smoother == 1)
        base_prec = delete_circles + 1;
    else if (smoother == 0)
        base_prec = -(delete_circles + 1);
    else if (smoother == 3)
        base_prec = ntheta_int;
    else if (smoother == 2)
        base_prec = -ntheta_int;

    /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!      Circle Smoother     !!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     */
    // Take boundary condition into account: Dirichlet-RB
    if (gyro::icntl[Param::DirBC_Interior]) { // (r[0],0) is on Dirichlet boundary
        dep_Asc_ortho[0][0] = 1;
        dep_Asc_ortho[1][0] = 1;
        start_j             = 1;
    }
    else { // (r[0],0) is not on Dirichlet boundary
        start_j = 0;
    }

    if (smoother < 2) {
        /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!! Link Radial-Circle (circle) !!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     */
        int j     = delete_circles;
        int odd_j = j % 2;

#pragma omp task shared(u, Au) firstprivate(smoother, j, odd_j) depend(in                                              \
                                                                       : dep_Asc_prev[j - 1])                          \
    depend(out                                                                                                         \
           : dep_Asc_ortho_cur[j])
        {
            int start_loop, shift_loop;
            int row, col, col2;
            int smoother_prev;
            double coeff, kt, ktmin1, hs, hsmin1, val, val2;
            std::vector<int> row_vect_prev, row_vect, row_vect_next;
            std::vector<int> smoother_vect_prev, smoother_vect, smoother_vect_next;
            std::vector<double> arr_vect, arr_vect2, att_vect, art_vect;

            smoother_prev      = get_smoother(0, j - 1);
            smoother_vect      = get_smoother_radial(j);
            smoother_vect_next = get_smoother_radial(j + 1);
            row_vect_prev      = get_row(j - 1, smoother_prev, extrapol, 0, 1);
            row_vect           = get_row(j, 2, extrapol, 0, 1);
            row_vect_next      = get_row(j + 1, 2, extrapol, 0, 1);
            hs                 = hplus[j];
            hsmin1             = hplus[j - 1];
            gyro::arr_att_art(r[j], theta, sin_theta, cos_theta, ntheta_int, arr_vect, att_vect, art_vect, 0);
            if (smoother_prev == smoother) {
                if (extrapol && smoother == 0) {
                    start_loop = 1;
                    shift_loop = 2;
                }
                else {
                    start_loop = 0;
                    shift_loop = 1;
                }
                for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                    kt     = thetaplus_per[i + 1];
                    ktmin1 = thetaplus_per[i];

                    // Update (j-1, i)
                    row = row_vect_prev[i];
                    col = j * ntheta_int + i;

                    coeff = 0.5 * (kt + ktmin1) * arr_vect[i] / hsmin1;
                    val   = -coeff;

                    Au[row] -= val * u[col];
                }
                if (gyro::icntl[Param::mod_pk] > 0) {
                    if (extrapol && smoother == 0) {
                        start_loop = 1;
                        shift_loop = 2;
                    }
                    else {
                        start_loop = 0;
                        shift_loop = 1;
                    }
                    for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                        kt     = thetaplus_per[i + 1];
                        ktmin1 = thetaplus_per[i];

                        // Update (j-1, i)
                        row  = row_vect_prev[i];
                        col  = j * ntheta_int + imoins(i, ntheta_int);
                        val  = art_vect[i];
                        col2 = j * ntheta_int + iplus(i, ntheta_int);
                        val2 = -art_vect[i];

                        Au[row] -= val * u[col] + val2 * u[col2];
                    }
                }
            }
            dep_Asc_ortho_cur[j] = 1;
        }

        /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!! Link Circle-Radial !!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     */
        j     = delete_circles - 1;
        odd_j = j % 2;
#pragma omp task shared(u, Au) firstprivate(smoother, j, odd_j) depend(in                                              \
                                                                       : dep_Asc_prev[j - odd_j])                      \
    depend(in                                                                                                          \
           : dep_Asc_ortho_cur[j + 1]) depend(out                                                                      \
                                              : dep_Asc_ortho_cur[j])
        {
            int start_loop, shift_loop;
            int row, row2, col, col2;
            int smoother_prev, smoother_cur;
            double coeff, coeff2, kt, ktmin1, hs, hsmin1, val, val2;
            std::vector<int> row_vect_prev, row_vect, row_vect_next;
            std::vector<int> smoother_vect_prev, smoother_vect, smoother_vect_next;
            std::vector<double> arr_vect, arr_vect2, att_vect, art_vect;

            if (delete_circles > 2 || !gyro::icntl[Param::DirBC_Interior]) {
                smoother_prev      = get_smoother(0, j - 1);
                smoother_cur       = get_smoother(0, j);
                smoother_vect_next = get_smoother_radial(j + 1);
                row_vect_prev      = get_row(j - 1, smoother_prev, extrapol, 0, 1);
                row_vect           = get_row(j, smoother_cur, extrapol, 0, 1);
                row_vect_next      = get_row(j + 1, 2, extrapol, 0, 1);
                hs                 = hplus[j];
                hsmin1             = hplus[j - 1];
                gyro::arr_att_art(r[j], theta, sin_theta, cos_theta, ntheta_int, arr_vect, att_vect, art_vect, 0);
                if (smoother_cur == smoother) {
                    if (extrapol && smoother == 0) {
                        start_loop = 1;
                        shift_loop = 2;
                    }
                    else {
                        start_loop = 0;
                        shift_loop = 1;
                    }
                    for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                        kt     = thetaplus_per[i + 1];
                        ktmin1 = thetaplus_per[i];

                        // Update (j, i)
                        row   = row_vect[i];
                        coeff = 0.5 * (kt + ktmin1) * arr_vect[i];
                        col   = (j + 1) * ntheta_int + i;
                        val   = -coeff / hs;
                        col2  = (j - 1) * ntheta_int + i;
                        val2  = -coeff / hsmin1;

                        Au[row] -= val * u[col] + val2 * u[col2];
                    }
                    if (extrapol && smoother == 0) {
                        start_loop = 0;
                        shift_loop = 2;
                        for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                            kt     = thetaplus_per[i + 1];
                            ktmin1 = thetaplus_per[i];

                            // Update (j, i+1)
                            row   = row_vect[iplus(i, ntheta_int)];
                            col   = j * ntheta_int + i;
                            coeff = 0.5 * (hs + hsmin1) * att_vect[i];
                            val   = -coeff / kt;

                            Au[row] -= val * u[col];

                            // Update (j, i-1)
                            row = row_vect[imoins(i, ntheta_int)];
                            val = -coeff / ktmin1;

                            Au[row] -= val * u[col];

                            kt     = thetaplus_per[i + 2];
                            ktmin1 = thetaplus_per[i + 1];
                            row    = row_vect[i + 1];
                            coeff2 = 0.5 * (hs + hsmin1) * att_vect[i + 1];
                            col    = j * ntheta_int + imoins(i + 1, ntheta_int);
                            val    = -coeff2 / ktmin1;
                            col2   = j * ntheta_int + iplus(i + 1, ntheta_int);
                            val2   = -coeff2 / kt;

                            Au[row] -= val * u[col] + val2 * u[col2];
                        }
                    }
                }
                else if (smoother_prev == smoother) {
                    if (extrapol && smoother == 0) {
                        start_loop = 1;
                        shift_loop = 2;
                    }
                    else {
                        start_loop = 0;
                        shift_loop = 1;
                    }
                    for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                        kt     = thetaplus_per[i + 1];
                        ktmin1 = thetaplus_per[i];

                        // Update (j-1, i)
                        row   = row_vect_prev[i];
                        col   = j * ntheta_int + i;
                        coeff = 0.5 * (kt + ktmin1) * arr_vect[i] / hsmin1;
                        val   = -coeff;

                        Au[row] -= val * u[col];
                    }
                }
                if (gyro::icntl[Param::mod_pk] > 0) {
                    if (smoother_cur == smoother) {
                        start_loop = 0;
                        if (extrapol && smoother == 0)
                            shift_loop = 2;
                        else
                            shift_loop = 1;
                        for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                            kt     = thetaplus_per[i + 1];
                            ktmin1 = thetaplus_per[i];

                            // Update (j, i+1)
                            row  = row_vect[iplus(i, ntheta_int)];
                            row2 = row_vect[imoins(i, ntheta_int)];
                            col  = (j - 1) * ntheta_int + i;
                            val  = -art_vect[i];
                            col2 = (j + 1) * ntheta_int + i;
                            val2 = art_vect[i];

                            val = val * u[col] + val2 * u[col2];
                            Au[row] -= val;

                            // Update (j, i-1)
                            Au[row2] -= -val;
                        }
                    }
                    else if (smoother_prev == smoother) {
                        if (extrapol && smoother == 0) {
                            start_loop = 1;
                            shift_loop = 2;
                        }
                        else {
                            start_loop = 0;
                            shift_loop = 1;
                        }
                        for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                            kt     = thetaplus_per[i + 1];
                            ktmin1 = thetaplus_per[i];

                            // Update (j-1, i)
                            row  = row_vect_prev[i];
                            col  = j * ntheta_int + imoins(i, ntheta_int);
                            val  = art_vect[i];
                            col2 = j * ntheta_int + iplus(i, ntheta_int);
                            val2 = -art_vect[i];

                            Au[row] -= val * u[col] + val2 * u[col2];
                        }
                    }
                }
            }
            dep_Asc_ortho_cur[j] = 1;
        }

        /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!  Interior nodes (1)   !!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     */
        for (int j = delete_circles - 4; j > start_j; j -= 3) {
            odd_j = j % 2;
#pragma omp task shared(u, Au) firstprivate(smoother, j, odd_j)                                                        \
    depend(in                                                                                                          \
           : dep_Asc_prev[j - odd_j], dep_Asc_prev[j + odd_j]) depend(out                                              \
                                                                      : dep_Asc_ortho_cur[j])
            {
                int start_loop, shift_loop;
                int row, row2, col, col2;
                int smoother_prev, smoother_cur, smoother_next;
                double coeff, coeff2, kt, ktmin1, hs, hsmin1, val, val2;
                std::vector<int> row_vect_prev, row_vect, row_vect_next;
                std::vector<int> smoother_vect_prev, smoother_vect, smoother_vect_next;
                std::vector<double> arr_vect, arr_vect2, att_vect, art_vect;

                smoother_prev = get_smoother(0, j - 1);
                smoother_cur  = get_smoother(0, j);
                smoother_next = get_smoother(0, j + 1);
                row_vect_prev = get_row(j - 1, smoother_prev, extrapol, 0, 1);
                row_vect      = get_row(j, smoother_cur, extrapol, 0, 1);
                row_vect_next = get_row(j + 1, smoother_next, extrapol, 0, 1);
                hs            = hplus[j];
                hsmin1        = hplus[j - 1];
                gyro::arr_att_art(r[j], theta, sin_theta, cos_theta, ntheta_int, arr_vect, att_vect, art_vect, 0);
                if (smoother_cur == smoother) {
                    if (extrapol && smoother == 0) {
                        start_loop = 1;
                        shift_loop = 2;
                    }
                    else {
                        start_loop = 0;
                        shift_loop = 1;
                    }
                    for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                        kt     = thetaplus_per[i + 1];
                        ktmin1 = thetaplus_per[i];

                        // Update (j, i)
                        row   = row_vect[i];
                        coeff = 0.5 * (kt + ktmin1) * arr_vect[i];
                        col   = (j + 1) * ntheta_int + i;
                        val   = -coeff / hs;
                        col2  = (j - 1) * ntheta_int + i;
                        val2  = -coeff / hsmin1;

                        Au[row] -= val * u[col] + val2 * u[col2];
                    }
                    if (extrapol && smoother == 0) {
                        start_loop = 0;
                        shift_loop = 2;
                        for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                            kt     = thetaplus_per[i + 1];
                            ktmin1 = thetaplus_per[i];

                            // Update (j, i+1)
                            row   = row_vect[iplus(i, ntheta_int)];
                            col   = j * ntheta_int + i;
                            coeff = 0.5 * (hs + hsmin1) * att_vect[i];
                            val   = -coeff / kt;

                            Au[row] -= val * u[col];

                            // Update (j, i-1)
                            row = row_vect[imoins(i, ntheta_int)];
                            val = -coeff / ktmin1;

                            Au[row] -= val * u[col];

                            // Update (j, i)
                            kt     = thetaplus_per[i + 2];
                            ktmin1 = thetaplus_per[i + 1];
                            row    = row_vect[i + 1];
                            coeff2 = 0.5 * (hs + hsmin1) * att_vect[i + 1];
                            col    = j * ntheta_int + imoins(i + 1, ntheta_int);
                            val    = -coeff2 / ktmin1;
                            col2   = j * ntheta_int + iplus(i + 1, ntheta_int);
                            val2   = -coeff2 / kt;

                            Au[row] -= val * u[col] + val2 * u[col2];
                        }
                    }
                }
                else {
                    if (extrapol && smoother == 0) {
                        start_loop = 1;
                        shift_loop = 2;
                    }
                    else {
                        start_loop = 0;
                        shift_loop = 1;
                    }
                    for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                        kt     = thetaplus_per[i + 1];
                        ktmin1 = thetaplus_per[i];

                        // Update (j-1, i)
                        row   = row_vect_prev[i];
                        col   = j * ntheta_int + i;
                        coeff = 0.5 * (kt + ktmin1) * arr_vect[i] / hsmin1;
                        val   = -coeff;

                        Au[row] -= val * u[col];

                        // Update (j+1, i)
                        row = row_vect_next[i];
                        col = j * ntheta_int + i;

                        coeff = 0.5 * (kt + ktmin1) * arr_vect[i] / hs;
                        val   = -coeff;

                        Au[row] -= val * u[col];
                    }
                }
                if (gyro::icntl[Param::mod_pk] > 0) {
                    if (smoother_cur == smoother) {
                        start_loop = 0;
                        if (extrapol && smoother == 0)
                            shift_loop = 2;
                        else
                            shift_loop = 1;
                        for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                            kt     = thetaplus_per[i + 1];
                            ktmin1 = thetaplus_per[i];

                            // Update (j, i+1)
                            row  = row_vect[iplus(i, ntheta_int)];
                            row2 = row_vect[imoins(i, ntheta_int)];
                            col  = (j - 1) * ntheta_int + i;
                            val  = -art_vect[i];
                            col2 = (j + 1) * ntheta_int + i;
                            val2 = art_vect[i];

                            val = val * u[col] + val2 * u[col2];
                            Au[row] -= val;

                            // Update (j, i-1)
                            Au[row2] -= -val;
                        }
                    }
                    else {
                        if (extrapol && smoother == 0) {
                            start_loop = 1;
                            shift_loop = 2;
                        }
                        else {
                            start_loop = 0;
                            shift_loop = 1;
                        }
                        for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                            kt     = thetaplus_per[i + 1];
                            ktmin1 = thetaplus_per[i];

                            // Update (j-1, i)
                            row  = row_vect_prev[i];
                            col  = j * ntheta_int + imoins(i, ntheta_int);
                            val  = art_vect[i];
                            col2 = j * ntheta_int + iplus(i, ntheta_int);
                            val2 = -art_vect[i];

                            Au[row] -= val * u[col] + val2 * u[col2];

                            // Update (j+1, i)
                            row  = row_vect_next[i];
                            val  = -art_vect[i];
                            val2 = art_vect[i];

                            Au[row] -= val * u[col] + val2 * u[col2];
                        }
                    }
                }
            }
            dep_Asc_ortho_cur[j] = 1;
        }

        /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!  Interior nodes (2)   !!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     */
        for (int j = delete_circles - 2; j > start_j; j -= 3) {
            odd_j = j % 2;
#pragma omp task shared(u, Au) firstprivate(smoother, j, odd_j)                                                        \
    depend(in                                                                                                          \
           : dep_Asc_prev[j - odd_j], dep_Asc_prev[j + odd_j], dep_Asc_ortho_cur[j + 1])                               \
        depend(in                                                                                                      \
               : dep_Asc_ortho_cur[j - 2]) depend(out                                                                  \
                                                  : dep_Asc_ortho_cur[j])
            {
                int start_loop, shift_loop;
                int row, row2, col, col2;
                int smoother_prev, smoother_cur, smoother_next;
                double coeff, coeff2, kt, ktmin1, hs, hsmin1, val, val2;
                std::vector<int> row_vect_prev, row_vect, row_vect_next;
                std::vector<int> smoother_vect_prev, smoother_vect, smoother_vect_next;
                std::vector<double> arr_vect, arr_vect2, att_vect, art_vect;

                smoother_prev = get_smoother(0, j - 1);
                smoother_cur  = get_smoother(0, j);
                smoother_next = get_smoother(0, j + 1);
                row_vect_prev = get_row(j - 1, smoother_prev, extrapol, 0, 1);
                row_vect      = get_row(j, smoother_cur, extrapol, 0, 1);
                row_vect_next = get_row(j + 1, smoother_next, extrapol, 0, 1);
                hs            = hplus[j];
                hsmin1        = hplus[j - 1];
                gyro::arr_att_art(r[j], theta, sin_theta, cos_theta, ntheta_int, arr_vect, att_vect, art_vect, 0);
                if (smoother_cur == smoother) {
                    if (extrapol && smoother == 0) {
                        start_loop = 1;
                        shift_loop = 2;
                    }
                    else {
                        start_loop = 0;
                        shift_loop = 1;
                    }
                    for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                        kt     = thetaplus_per[i + 1];
                        ktmin1 = thetaplus_per[i];

                        // Update (j, i)
                        row   = row_vect[i];
                        coeff = 0.5 * (kt + ktmin1) * arr_vect[i];
                        col   = (j + 1) * ntheta_int + i;
                        val   = -coeff / hs;
                        col2  = (j - 1) * ntheta_int + i;
                        val2  = -coeff / hsmin1;

                        Au[row] -= val * u[col] + val2 * u[col2];
                    }
                    if (extrapol && smoother == 0) {
                        start_loop = 0;
                        shift_loop = 2;
                        for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                            kt     = thetaplus_per[i + 1];
                            ktmin1 = thetaplus_per[i];

                            // Update (j, i+1)
                            row   = row_vect[iplus(i, ntheta_int)];
                            col   = j * ntheta_int + i;
                            coeff = 0.5 * (hs + hsmin1) * att_vect[i];
                            val   = -coeff / kt;

                            Au[row] -= val * u[col];

                            // Update (j, i-1)
                            row = row_vect[imoins(i, ntheta_int)];
                            val = -coeff / ktmin1;

                            Au[row] -= val * u[col];

                            // Update (j, i)
                            kt     = thetaplus_per[i + 2];
                            ktmin1 = thetaplus_per[i + 1];
                            row    = row_vect[i + 1];
                            coeff2 = 0.5 * (hs + hsmin1) * att_vect[i + 1];
                            col    = j * ntheta_int + imoins(i + 1, ntheta_int);
                            val    = -coeff2 / ktmin1;
                            col2   = j * ntheta_int + iplus(i + 1, ntheta_int);
                            val2   = -coeff2 / kt;

                            Au[row] -= val * u[col] + val2 * u[col2];
                        }
                    }
                }
                else {
                    if (extrapol && smoother == 0) {
                        start_loop = 1;
                        shift_loop = 2;
                    }
                    else {
                        start_loop = 0;
                        shift_loop = 1;
                    }
                    for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                        kt     = thetaplus_per[i + 1];
                        ktmin1 = thetaplus_per[i];

                        // Update (j-1, i)
                        row   = row_vect_prev[i];
                        col   = j * ntheta_int + i;
                        coeff = 0.5 * (kt + ktmin1) * arr_vect[i] / hsmin1;
                        val   = -coeff;

                        Au[row] -= val * u[col];

                        // Update (j+1, i)
                        row = row_vect_next[i];
                        col = j * ntheta_int + i;

                        coeff = 0.5 * (kt + ktmin1) * arr_vect[i] / hs;
                        val   = -coeff;

                        Au[row] -= val * u[col];
                    }
                }
                if (gyro::icntl[Param::mod_pk] > 0) {
                    if (smoother_cur == smoother) {
                        start_loop = 0;
                        if (extrapol && smoother == 0)
                            shift_loop = 2;
                        else
                            shift_loop = 1;
                        for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                            kt     = thetaplus_per[i + 1];
                            ktmin1 = thetaplus_per[i];

                            // Update (j, i+1)
                            row  = row_vect[iplus(i, ntheta_int)];
                            row2 = row_vect[imoins(i, ntheta_int)];
                            col  = (j - 1) * ntheta_int + i;
                            val  = -art_vect[i];
                            col2 = (j + 1) * ntheta_int + i;
                            val2 = art_vect[i];

                            val = val * u[col] + val2 * u[col2];
                            Au[row] -= val;

                            // Update (j, i-1)
                            Au[row2] -= -val;
                        }
                    }
                    else {
                        if (extrapol && smoother == 0) {
                            start_loop = 1;
                            shift_loop = 2;
                        }
                        else {
                            start_loop = 0;
                            shift_loop = 1;
                        }
                        for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                            kt     = thetaplus_per[i + 1];
                            ktmin1 = thetaplus_per[i];

                            // Update (j-1, i)
                            row  = row_vect_prev[i];
                            col  = j * ntheta_int + imoins(i, ntheta_int);
                            val  = art_vect[i];
                            col2 = j * ntheta_int + iplus(i, ntheta_int);
                            val2 = -art_vect[i];

                            Au[row] -= val * u[col] + val2 * u[col2];

                            // Update (j+1, i)
                            row  = row_vect_next[i];
                            val  = -art_vect[i];
                            val2 = art_vect[i];

                            Au[row] -= val * u[col] + val2 * u[col2];
                        }
                    }
                }
            }
            dep_Asc_ortho_cur[j] = 1;
        }

        /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!  Interior nodes (3)   !!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     */
        for (int j = delete_circles - 3; j > start_j; j -= 3) {
            odd_j = j % 2;
#pragma omp task shared(u, Au) firstprivate(smoother, j, odd_j)                                                        \
    depend(in                                                                                                          \
           : dep_Asc_prev[j - odd_j], dep_Asc_prev[j + odd_j], dep_Asc_ortho_cur[j + 1])                               \
        depend(in                                                                                                      \
               : dep_Asc_ortho_cur[j - 2]) depend(out                                                                  \
                                                  : dep_Asc_ortho_cur[j])
            {
                int start_loop, shift_loop;
                int row, row2, col, col2;
                int smoother_prev, smoother_cur, smoother_next;
                double coeff, coeff2, kt, ktmin1, hs, hsmin1, val, val2;
                std::vector<int> row_vect_prev, row_vect, row_vect_next;
                std::vector<int> smoother_vect_prev, smoother_vect, smoother_vect_next;
                std::vector<double> arr_vect, arr_vect2, att_vect, art_vect;

                smoother_prev = get_smoother(0, j - 1);
                smoother_cur  = get_smoother(0, j);
                smoother_next = get_smoother(0, j + 1);
                row_vect_prev = get_row(j - 1, smoother_prev, extrapol, 0, 1);
                row_vect      = get_row(j, smoother_cur, extrapol, 0, 1);
                row_vect_next = get_row(j + 1, smoother_next, extrapol, 0, 1);
                hs            = hplus[j];
                hsmin1        = hplus[j - 1];
                gyro::arr_att_art(r[j], theta, sin_theta, cos_theta, ntheta_int, arr_vect, att_vect, art_vect, 0);
                if (smoother_cur == smoother) {
                    if (extrapol && smoother == 0) {
                        start_loop = 1;
                        shift_loop = 2;
                    }
                    else {
                        start_loop = 0;
                        shift_loop = 1;
                    }
                    for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                        kt     = thetaplus_per[i + 1];
                        ktmin1 = thetaplus_per[i];

                        // Update (j, i)
                        row   = row_vect[i];
                        coeff = 0.5 * (kt + ktmin1) * arr_vect[i];
                        col   = (j + 1) * ntheta_int + i;
                        val   = -coeff / hs;
                        col2  = (j - 1) * ntheta_int + i;
                        val2  = -coeff / hsmin1;

                        Au[row] -= val * u[col] + val2 * u[col2];
                    }
                    if (extrapol && smoother == 0) {
                        start_loop = 0;
                        shift_loop = 2;
                        for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                            kt     = thetaplus_per[i + 1];
                            ktmin1 = thetaplus_per[i];

                            // Update (j, i+1)
                            row   = row_vect[iplus(i, ntheta_int)];
                            col   = j * ntheta_int + i;
                            coeff = 0.5 * (hs + hsmin1) * att_vect[i];
                            val   = -coeff / kt;

                            Au[row] -= val * u[col];

                            // Update (j, i-1)
                            row = row_vect[imoins(i, ntheta_int)];
                            val = -coeff / ktmin1;

                            Au[row] -= val * u[col];

                            // Update (j, i)
                            kt     = thetaplus_per[i + 2];
                            ktmin1 = thetaplus_per[i + 1];
                            row    = row_vect[i + 1];
                            coeff2 = 0.5 * (hs + hsmin1) * att_vect[i + 1];
                            col    = j * ntheta_int + imoins(i + 1, ntheta_int);
                            val    = -coeff2 / ktmin1;
                            col2   = j * ntheta_int + iplus(i + 1, ntheta_int);
                            val2   = -coeff2 / kt;

                            Au[row] -= val * u[col] + val2 * u[col2];
                        }
                    }
                }
                else {
                    if (extrapol && smoother == 0) {
                        start_loop = 1;
                        shift_loop = 2;
                    }
                    else {
                        start_loop = 0;
                        shift_loop = 1;
                    }
                    for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                        kt     = thetaplus_per[i + 1];
                        ktmin1 = thetaplus_per[i];

                        // Update (j-1, i)
                        row   = row_vect_prev[i];
                        col   = j * ntheta_int + i;
                        coeff = 0.5 * (kt + ktmin1) * arr_vect[i] / hsmin1;
                        val   = -coeff;

                        Au[row] -= val * u[col];

                        // Update (j+1, i)
                        row = row_vect_next[i];
                        col = j * ntheta_int + i;

                        coeff = 0.5 * (kt + ktmin1) * arr_vect[i] / hs;
                        val   = -coeff;

                        Au[row] -= val * u[col];
                    }
                }
                if (gyro::icntl[Param::mod_pk] > 0) {
                    if (smoother_cur == smoother) {
                        start_loop = 0;
                        if (extrapol && smoother == 0)
                            shift_loop = 2;
                        else
                            shift_loop = 1;
                        for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                            kt     = thetaplus_per[i + 1];
                            ktmin1 = thetaplus_per[i];

                            // Update (j, i+1)
                            row  = row_vect[iplus(i, ntheta_int)];
                            row2 = row_vect[imoins(i, ntheta_int)];
                            col  = (j - 1) * ntheta_int + i;
                            val  = -art_vect[i];
                            col2 = (j + 1) * ntheta_int + i;
                            val2 = art_vect[i];

                            val = val * u[col] + val2 * u[col2];
                            Au[row] -= val;

                            // Update (j, i-1)
                            Au[row2] -= -val;
                        }
                    }
                    else {
                        if (extrapol && smoother == 0) {
                            start_loop = 1;
                            shift_loop = 2;
                        }
                        else {
                            start_loop = 0;
                            shift_loop = 1;
                        }
                        for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                            kt     = thetaplus_per[i + 1];
                            ktmin1 = thetaplus_per[i];

                            // Update (j-1, i)
                            row  = row_vect_prev[i];
                            col  = j * ntheta_int + imoins(i, ntheta_int);
                            val  = art_vect[i];
                            col2 = j * ntheta_int + iplus(i, ntheta_int);
                            val2 = -art_vect[i];

                            Au[row] -= val * u[col] + val2 * u[col2];

                            // Update (j+1, i)
                            row  = row_vect_next[i];
                            val  = -art_vect[i];
                            val2 = art_vect[i];

                            Au[row] -= val * u[col] + val2 * u[col2];
                        }
                    }
                }
            }
            dep_Asc_ortho_cur[j] = 1;
        }

        /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 * !!!!!!!!!!!    First lines    !!!!!!!!!!!
                 * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 */
        j     = start_j;
        odd_j = j % 2;
#pragma omp task shared(u, Au) firstprivate(smoother, j, odd_j)                                                        \
    depend(in                                                                                                          \
           : dep_Asc_prev[j - odd_j], dep_Asc_prev[j + odd_j], dep_Asc_ortho_cur[j + 1])                               \
        depend(in                                                                                                      \
               : dep_Asc_ortho_cur[j + 2]) depend(out                                                                  \
                                                  : dep_Asc_ortho_cur[j])
        {
            int start_loop, shift_loop;
            int row, col, col2;
            int smoother_cur, smoother_next;
            double coeff, coeff2, kt, ktmin1, hs, hsmin1, val, val2;
            std::vector<int> row_vect_prev, row_vect, row_vect_next;
            std::vector<int> smoother_vect_prev, smoother_vect, smoother_vect_next;
            std::vector<double> arr_vect, arr_vect2, att_vect, art_vect;

            smoother_cur = j;
            row_vect     = get_row(j, smoother_cur, extrapol, 0, 1);
            if (j < delete_circles - 1) {
                smoother_next = (j + 1) % 2;
                row_vect_next = get_row(j + 1, smoother_next, extrapol, 0, 1);
            }
            else if (j == delete_circles - 1) {
                smoother_vect_next = get_smoother_radial(j + 1);
                smoother_next      = smoother_vect_next[0];
                row_vect_next      = get_row(j + 1, 2, extrapol, 0, 1);
            }

            hs = hplus[j];
            if (j > 0) {
                hsmin1 = hplus[j - 1];
            }
            else {
                hsmin1 = 2 * r[0]; // across the origin
            }
            // Across and DB_int updates
            gyro::arr_att_art(r[j], theta, sin_theta, cos_theta, ntheta_int, arr_vect, att_vect, art_vect, 0);
            if (smoother_cur == smoother) {
                if (extrapol && smoother == 0) {
                    start_loop = 1;
                    shift_loop = 2;
                }
                else {
                    start_loop = 0;
                    shift_loop = 1;
                }
                for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                    kt     = thetaplus_per[i + 1];
                    ktmin1 = thetaplus_per[i];

                    // Update (j, i)
                    row   = row_vect[i];
                    coeff = 0.5 * (kt + ktmin1) * arr_vect[i];
                    col   = (j + 1) * ntheta_int + i;
                    val   = -coeff / hs;

                    Au[row] -= val * u[col];
                }
                if (extrapol && smoother == 0) {
                    start_loop = 0;
                    shift_loop = 2;
                    for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                        kt     = thetaplus_per[i + 1];
                        ktmin1 = thetaplus_per[i];

                        // Update (j, i)
                        row   = row_vect[iplus(i, ntheta_int)];
                        col   = j * ntheta_int + i;
                        coeff = 0.5 * (hs + hsmin1) * att_vect[i];
                        val   = -coeff / kt;

                        Au[row] -= val * u[col];

                        // Update (j, i-1)
                        row = row_vect[imoins(i, ntheta_int)];
                        val = -coeff / ktmin1;

                        Au[row] -= val * u[col];

                        // Update (j, i)
                        kt     = thetaplus_per[i + 2];
                        ktmin1 = thetaplus_per[i + 1];
                        row    = row_vect[i + 1];
                        coeff2 = 0.5 * (hs + hsmin1) * att_vect[i + 1];
                        col    = j * ntheta_int + imoins(i + 1, ntheta_int);
                        val    = -coeff2 / ktmin1;
                        col2   = j * ntheta_int + iplus(i + 1, ntheta_int);
                        val2   = -coeff2 / kt;

                        Au[row] -= val * u[col] + val2 * u[col2];
                    }
                }
                if (gyro::icntl[Param::mod_pk] > 0) {
                    start_loop = 0;
                    if (extrapol && smoother == 0)
                        shift_loop = 2;
                    else
                        shift_loop = 1;
                    for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                        kt     = thetaplus_per[i + 1];
                        ktmin1 = thetaplus_per[i];

                        // Update (j, i+1)
                        row = row_vect[iplus(i, ntheta_int)];
                        col = (j + 1) * ntheta_int + i;
                        val = art_vect[i];

                        Au[row] -= val * u[col];

                        // Update (j, i-1)
                        row = row_vect[imoins(i, ntheta_int)];

                        Au[row] -= -val * u[col];
                    }
                }
            }

            // Update (j+1, i)
            if (j < delete_circles - 1 && smoother_next == smoother) {
                if (extrapol && smoother == 0) {
                    start_loop = 1;
                    shift_loop = 2;
                }
                else {
                    start_loop = 0;
                    shift_loop = 1;
                }
                for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                    kt     = thetaplus_per[i + 1];
                    ktmin1 = thetaplus_per[i];

                    row   = row_vect_next[i];
                    col   = j * ntheta_int + i;
                    coeff = 0.5 * (kt + ktmin1) * arr_vect[i] / hs;

                    val = -coeff;

                    Au[row] -= val * u[col];
                }
                if (gyro::icntl[Param::mod_pk] > 0)
                    for (int i = start_loop; i < ntheta_int; i += shift_loop) {
                        row  = row_vect_next[i];
                        col  = j * ntheta_int + imoins(i, ntheta_int);
                        val  = -art_vect[i];
                        col2 = j * ntheta_int + iplus(i, ntheta_int);
                        val2 = art_vect[i];

                        Au[row] -= val * u[col] + val2 * u[col2];
                    }
            }
            dep_Asc_ortho_cur[j] = 1;
        }
    }

    /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! RADIAL   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         */
    if (smoother > 1) {
        int diff  = ntheta_int % 3;
        int odd_d = delete_circles % 2;

        /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                * !!!!!!!!!!!!!!!!!!!!!!!            (0))           !!!!!!!!!!!!!!!!!!!!!!
                * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                */

        for (int i = 0; i < diff; i++) {
            {
                int im        = imoins(i, ntheta_int);
                int ip        = iplus(i, ntheta_int);
                int im2       = imoins(im, ntheta_int);
                int ip2       = iplus(ip, ntheta_int);
                int odd_i     = i % 2;
                int ind_m_odd = i, ind_p_odd = i;
                if (odd_i) {
                    ind_m_odd = im;
                    ind_p_odd = ip;
                }
#pragma omp task shared(u, Au) firstprivate(smoother, im, i, ip, im2, ip2, odd_i, odd_d)                               \
    depend(in                                                                                                          \
           : dep_Asc_prev[ind_m_odd]) depend(in                                                                        \
                                             : dep_Asc_prev[ind_p_odd], dep_Asc1[delete_circles - 1 - odd_d])          \
        depend(in                                                                                                      \
               : dep_Asc_ortho_cur[i - 1]) depend(out                                                                  \
                                                  : dep_Asc_ortho_cur[i])
                {
                    int row, row2, col, col2;
                    int smoother_prev, smoother_cur, smoother_next;
                    double coeff, coeff2, coeff3, kt, ktmin1, hs, hsmin1, val, val2;
                    std::vector<int> row_vect_prev, row_vect, row_vect_next;
                    std::vector<int> smoother_vect_prev, smoother_vect, smoother_vect_next;
                    std::vector<double> arr_vect, arr_vect2, att_vect, art_vect;

                    int j;
                    int s_cur     = (i % 2) + 2;
                    int s_next    = ((i + 1) % 2) + 2;
                    kt            = thetaplus_per[i + 1];
                    ktmin1        = thetaplus_per[i];
                    smoother_cur  = s_cur;
                    smoother_prev = s_next;
                    smoother_next = s_next;
                    row_vect_prev = get_row_i_glob(nr, im, smoother_prev, extrapol);
                    row_vect      = get_row_i_glob(nr, i, smoother_cur, extrapol);
                    row_vect_next = get_row_i_glob(nr, ip, smoother_next, extrapol);
                    gyro::arr_att_art(r, theta[i], arr_vect, att_vect, art_vect, 0);

                    /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                         * !!!!!!!!!!! Link Circle-Radial (radial) !!!!!!!!!!
                         * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                         */
                    j = delete_circles - 1;
                    if (!(extrapol && smoother == 2 && j % 2 == 1)) {
                        hs     = hplus[j];
                        hsmin1 = hplus[j - 1];
                        if ((smoother == 3 && i % 2 == 1) || (smoother == 2 && i % 2 == 0)) {
                            row   = row_vect[j + 1];
                            col   = j * ntheta_int + i;
                            coeff = 0.5 * (kt + ktmin1) * arr_vect[j] / hs;
                            val   = -coeff;

                            Au[row] -= val * u[col];
                            if (gyro::icntl[Param::mod_pk] > 0) {
                                row  = row_vect[j + 1];
                                col  = j * ntheta_int + im;
                                val  = -art_vect[j];
                                col2 = j * ntheta_int + ip;
                                val2 = art_vect[j];

                                Au[row] -= val * u[col] + val2 * u[col2];
                            }
                        }
                    }

                    /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             * !!!!!!!!!!! Link Radial-Circle (radial) !!!!!!!!!!
             * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             */

                    j      = delete_circles;
                    hs     = hplus[j];
                    hsmin1 = hplus[j - 1];
                    if (!(extrapol && smoother == 2 && j % 2 == 0)) {
                        if ((smoother == 3 && i % 2 == 1) || (smoother == 2 && i % 2 == 0)) {
                            // Update (j, i)
                            row   = row_vect[j];
                            coeff = 0.5 * (kt + ktmin1) * arr_vect[j];
                            col   = (j - 1) * ntheta_int + i;
                            val   = -coeff / hsmin1;

                            Au[row] -= val * u[col];

                            coeff2 = 0.5 * (hs + hsmin1) * att_vect[j];
                            col    = j * ntheta_int + im;
                            val    = -coeff2 / ktmin1;
                            col2   = j * ntheta_int + ip;
                            val2   = -coeff2 / kt;

                            Au[row] -= val * u[col] + val2 * u[col2];
                        }
                        if ((smoother == 3 && i % 2 == 0) || (smoother == 2 && i % 2 == 1)) {
                            // Update (j, i+1)
                            row   = row_vect_next[j];
                            col   = j * ntheta_int + i;
                            coeff = 0.5 * (hs + hsmin1) * att_vect[j] / kt;
                            val   = -coeff;

                            Au[row] -= val * u[col];

                            // Update (j, i-1)
                            row   = row_vect_prev[j];
                            col   = j * ntheta_int + i;
                            coeff = 0.5 * (hs + hsmin1) * att_vect[j] / ktmin1;
                            val   = -coeff;

                            Au[row] -= val * u[col];
                        }
                    }
                    if (extrapol && smoother == 2 && i % 2 == 0) {
                        if (j % 2 == 1) {
                            row   = row_vect[j];
                            coeff = 0.5 * (kt + ktmin1) * arr_vect[j];
                            col   = (j + 1) * ntheta_int + i;
                            val   = -coeff / hs;

                            Au[row] -= val * u[col];
                        }
                        else {
                            row   = row_vect[j + 1];
                            col   = j * ntheta_int + i;
                            coeff = 0.5 * (kt + ktmin1) * arr_vect[j] / hs;
                            val   = -coeff;

                            Au[row] -= val * u[col];
                        }
                    }
                    if (gyro::icntl[Param::mod_pk] > 0) {
                        if (!(extrapol && smoother == 2 && j % 2 == 0)) {
                            if ((smoother == 3 && i % 2 == 0) || (smoother == 2 && i % 2 == 1)) {
                                row  = row_vect_next[j];
                                col  = (j - 1) * ntheta_int + i;
                                val  = -art_vect[j];
                                row2 = row_vect_prev[j];
                                col2 = (j + 1) * ntheta_int + i;
                                val2 = art_vect[j];

                                val = val * u[col] + val2 * u[col2];
                                Au[row] -= val;

                                // Update (j, i-1)
                                Au[row2] -= -val;
                            }
                        }
                        if (!(extrapol && smoother == 2 && j % 2 == 1)) {
                            if ((smoother == 3 && i % 2 == 1) || (smoother == 2 && i % 2 == 0)) {
                                // Update (j+1, i)
                                row  = row_vect[j + 1];
                                col  = j * ntheta_int + im;
                                val  = -art_vect[j];
                                col2 = j * ntheta_int + ip;
                                val2 = art_vect[j];

                                Au[row] -= val * u[col] + val2 * u[col2];
                            }
                        }
                    }

                    /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!  Interior nodes   !!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     */
                    for (j = delete_circles + 1; j < nr_int - 1; j++) {
                        hs     = hplus[j];
                        hsmin1 = hplus[j - 1];
                        if (smoother > 1 && !(extrapol && smoother == 2 && j % 2 == 0)) {
                            if ((smoother == 3 && i % 2 == 1) || (smoother == 2 && i % 2 == 0)) {
                                // Update (j, i)
                                row = row_vect[j];
                                // coeff  = 0.5 * (kt + ktmin1) * arr_vect[j];
                                coeff2 = 0.5 * (hs + hsmin1) * att_vect[j];
                                col    = j * ntheta_int + im;
                                val    = -coeff2 / ktmin1;
                                col2   = j * ntheta_int + ip;
                                val2   = -coeff2 / kt;

                                Au[row] -= val * u[col] + val2 * u[col2];
                            }
                            if ((smoother == 3 && i % 2 == 0) || (smoother == 2 && i % 2 == 1)) {
                                // Update (j, i+1)
                                row   = row_vect_next[j];
                                col   = j * ntheta_int + i;
                                coeff = 0.5 * (hs + hsmin1) * att_vect[j];
                                val   = -coeff / kt;

                                Au[row] -= val * u[col];

                                // Update (j, i-1)
                                row   = row_vect_prev[j];
                                col   = j * ntheta_int + i;
                                coeff = 0.5 * (hs + hsmin1) * att_vect[j];
                                val   = -coeff / ktmin1;

                                Au[row] -= val * u[col];
                            }
                        }
                        if (extrapol && smoother == 2 && i % 2 == 0) {
                            if (j % 2 == 1) {
                                // Update (j, i)
                                row   = row_vect[j];
                                coeff = 0.5 * (kt + ktmin1) * arr_vect[j];
                                col   = (j + 1) * ntheta_int + i;
                                val   = -coeff / hs;
                                col2  = (j - 1) * ntheta_int + i;
                                val2  = -coeff / hsmin1;

                                Au[row] -= val * u[col] + val2 * u[col2];
                            }
                            else {
                                // Update (j-1, i)
                                row   = row_vect[j - 1];
                                col   = j * ntheta_int + i;
                                coeff = 0.5 * (kt + ktmin1) * arr_vect[j] / hsmin1;
                                val   = -coeff;

                                Au[row] -= val * u[col];

                                // Update (j+1, i)
                                row   = row_vect[j + 1];
                                col   = j * ntheta_int + i;
                                coeff = 0.5 * (kt + ktmin1) * arr_vect[j] / hs;
                                val   = -coeff;

                                Au[row] -= val * u[col];
                            }
                        }
                        if (gyro::icntl[Param::mod_pk] > 0) {
                            if (smoother > 1 && !(extrapol && smoother == 2 && j % 2 == 1)) {
                                if ((smoother == 3 && i % 2 == 1) || (smoother == 2 && i % 2 == 0)) {
                                    // Update (j-1, i)
                                    row  = row_vect[j - 1];
                                    row2 = row_vect[j + 1];
                                    col  = j * ntheta_int + im;
                                    val  = art_vect[j];
                                    col2 = j * ntheta_int + ip;
                                    val2 = -art_vect[j];

                                    val = val * u[col] + val2 * u[col2];
                                    Au[row] -= val;

                                    // Update (j+1, i)
                                    Au[row2] -= -val;
                                }
                            }
                            if (smoother > 1 && !(extrapol && smoother == 2 && j % 2 == 0)) {
                                if ((smoother == 3 && i % 2 == 0) || (smoother == 2 && i % 2 == 1)) {
                                    // Update (j, i+1)
                                    row  = row_vect_next[j];
                                    col  = (j - 1) * ntheta_int + i;
                                    val  = -art_vect[j];
                                    col2 = (j + 1) * ntheta_int + i;
                                    val2 = art_vect[j];

                                    val = val * u[col] + val2 * u[col2];
                                    Au[row] -= val;

                                    // Update (j, i-1)
                                    row  = row_vect_prev[j];
                                    col  = (j - 1) * ntheta_int + i;
                                    val  = -art_vect[j];
                                    col2 = (j + 1) * ntheta_int + i;
                                    val2 = art_vect[j];

                                    val = val * u[col] + val2 * u[col2];
                                    Au[row] -= -val;
                                }
                            }
                        }
                    }
                    /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!        Last lines       !!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     */
                    // DB_ext updates (~~~ Interior - (j+1, i) + DB)
                    j      = nr_int - 1;
                    hs     = hplus[j];
                    hsmin1 = hplus[j - 1];
                    if (smoother > 1 && !(extrapol && smoother == 2 && j % 2 == 0)) {
                        if ((smoother == 3 && i % 2 == 1) || (smoother == 2 && i % 2 == 0)) {
                            row = row_vect[j];
                            // Contribution to middle (top) from DB
                            coeff3 = 0.5 * (hs + hsmin1) * att_vect[j];
                            col    = j * ntheta_int + im;
                            val    = -coeff3 / ktmin1;
                            col2   = j * ntheta_int + ip;
                            val2   = -coeff3 / kt;

                            Au[row] -= val * u[col] + val2 * u[col2];
                        }
                        if ((smoother == 3 && i % 2 == 0) || (smoother == 2 && i % 2 == 1)) {
                            // Update (j, i+1)
                            row   = row_vect_next[j];
                            col   = j * ntheta_int + i;
                            coeff = 0.5 * (hs + hsmin1) * att_vect[j];
                            val   = -coeff / kt;

                            Au[row] -= val * u[col];

                            // Update (j, i-1)
                            row   = row_vect_prev[j];
                            col   = j * ntheta_int + i;
                            coeff = 0.5 * (hs + hsmin1) * att_vect[j];
                            val   = -coeff / ktmin1;

                            Au[row] -= val * u[col];
                        }
                    }
                    if (extrapol && smoother == 2 && i % 2 == 0) {
                        if (j % 2 == 0) {
                            // Update (j-1, i)
                            row   = row_vect[j - 1];
                            col   = j * ntheta_int + i;
                            coeff = 0.5 * (kt + ktmin1) * arr_vect[j] / hsmin1;
                            val   = -coeff;

                            Au[row] -= val * u[col];
                        }
                        else {
                            // Update (j, i)
                            row   = row_vect[j];
                            coeff = 0.5 * (kt + ktmin1) * arr_vect[j];
                            col   = (j - 1) * ntheta_int + i;
                            val   = -coeff / hsmin1;

                            Au[row] -= val * u[col];
                        }
                    }
                    if (gyro::icntl[Param::mod_pk] > 0) {
                        if (smoother > 1 && !(extrapol && smoother == 2 && j % 2 == 1)) {
                            if ((smoother == 3 && i % 2 == 1) || (smoother == 2 && i % 2 == 0)) {
                                // Update (j-1, i)
                                row  = row_vect[j - 1];
                                col  = j * ntheta_int + im;
                                val  = art_vect[j];
                                col2 = j * ntheta_int + ip;
                                val2 = -art_vect[j];

                                Au[row] -= val * u[col] + val2 * u[col2];
                            }
                        }
                        if (smoother > 1 && !(extrapol && smoother == 2 && j % 2 == 0)) {
                            if ((smoother == 3 && i % 2 == 0) || (smoother == 2 && i % 2 == 1)) {
                                row = row_vect_next[j];
                                col = (j - 1) * ntheta_int + i;
                                val = -art_vect[j];

                                Au[row] -= val * u[col];

                                // Update (j, i-1)
                                row = row_vect_prev[j];
                                col = (j - 1) * ntheta_int + i;
                                val = -art_vect[j];

                                Au[row] -= -val * u[col];
                            }
                        }
                    }
                }
                dep_Asc_ortho_cur[i] = 1;
            }
        }

        /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                * !!!!!!!!!!!!!!!!!!!!!!!            (1))           !!!!!!!!!!!!!!!!!!!!!!
                * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                */
        for (int i = diff; i < ntheta_int; i += 3) {
            {
                int im        = imoins(i, ntheta_int);
                int ip        = iplus(i, ntheta_int);
                int im2       = imoins(im, ntheta_int);
                int ip2       = iplus(ip, ntheta_int);
                int odd_i     = i % 2;
                int ind_m_odd = i, ind_p_odd = i;
                if (odd_i) {
                    ind_m_odd = im;
                    ind_p_odd = ip;
                }
#pragma omp task shared(u, Au) firstprivate(smoother, im, i, ip, im2, ip2, odd_i) depend(in                            \
                                                                                         : dep_Asc_prev[ind_m_odd])    \
    depend(in                                                                                                          \
           : dep_Asc_prev[ind_p_odd], dep_Asc_ortho_cur[diff - 2], dep_Asc_ortho_cur[diff - 1])                        \
        depend(out                                                                                                     \
               : dep_Asc_ortho_cur[i])
                {
                    int row, row2, col, col2;
                    int smoother_prev, smoother_cur, smoother_next;
                    double coeff, coeff2, coeff3, kt, ktmin1, hs, hsmin1, val, val2;
                    std::vector<int> row_vect_prev, row_vect, row_vect_next;
                    std::vector<int> smoother_vect_prev, smoother_vect, smoother_vect_next;
                    std::vector<double> arr_vect, arr_vect2, att_vect, art_vect;

                    int j;
                    int s_cur     = (i % 2) + 2;
                    int s_next    = ((i + 1) % 2) + 2;
                    kt            = thetaplus_per[i + 1];
                    ktmin1        = thetaplus_per[i];
                    smoother_cur  = s_cur;
                    smoother_prev = s_next;
                    smoother_next = s_next;
                    row_vect_prev = get_row_i_glob(nr, im, smoother_prev, extrapol);
                    row_vect      = get_row_i_glob(nr, i, smoother_cur, extrapol);
                    row_vect_next = get_row_i_glob(nr, ip, smoother_next, extrapol);
                    gyro::arr_att_art(r, theta[i], arr_vect, att_vect, art_vect, 0);

                    /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                         * !!!!!!!!!!! Link Circle-Radial (radial) !!!!!!!!!!
                         * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                         */
                    j = delete_circles - 1;
                    if (!(extrapol && smoother == 2 && j % 2 == 1)) {
                        hs     = hplus[j];
                        hsmin1 = hplus[j - 1];
                        if ((smoother == 3 && i % 2 == 1) || (smoother == 2 && i % 2 == 0)) {
                            row   = row_vect[j + 1];
                            col   = j * ntheta_int + i;
                            coeff = 0.5 * (kt + ktmin1) * arr_vect[j] / hs;
                            val   = -coeff;

                            Au[row] -= val * u[col];
                            if (gyro::icntl[Param::mod_pk] > 0) {
                                row  = row_vect[j + 1];
                                col  = j * ntheta_int + im;
                                val  = -art_vect[j];
                                col2 = j * ntheta_int + ip;
                                val2 = art_vect[j];

                                Au[row] -= val * u[col] + val2 * u[col2];
                            }
                        }
                    }

                    /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             * !!!!!!!!!!! Link Radial-Circle (radial) !!!!!!!!!!
             * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             */

                    j      = delete_circles;
                    hs     = hplus[j];
                    hsmin1 = hplus[j - 1];
                    if (!(extrapol && smoother == 2 && j % 2 == 0)) {
                        if ((smoother == 3 && i % 2 == 1) || (smoother == 2 && i % 2 == 0)) {
                            // Update (j, i)
                            row   = row_vect[j];
                            coeff = 0.5 * (kt + ktmin1) * arr_vect[j];
                            col   = (j - 1) * ntheta_int + i;
                            val   = -coeff / hsmin1;

                            Au[row] -= val * u[col];

                            coeff2 = 0.5 * (hs + hsmin1) * att_vect[j];
                            col    = j * ntheta_int + im;
                            val    = -coeff2 / ktmin1;
                            col2   = j * ntheta_int + ip;
                            val2   = -coeff2 / kt;

                            Au[row] -= val * u[col] + val2 * u[col2];
                        }
                        if ((smoother == 3 && i % 2 == 0) || (smoother == 2 && i % 2 == 1)) {
                            // Update (j, i+1)
                            row   = row_vect_next[j];
                            col   = j * ntheta_int + i;
                            coeff = 0.5 * (hs + hsmin1) * att_vect[j] / kt;
                            val   = -coeff;

                            Au[row] -= val * u[col];

                            // Update (j, i-1)
                            row   = row_vect_prev[j];
                            col   = j * ntheta_int + i;
                            coeff = 0.5 * (hs + hsmin1) * att_vect[j] / ktmin1;
                            val   = -coeff;

                            Au[row] -= val * u[col];
                        }
                    }
                    if (extrapol && smoother == 2 && i % 2 == 0) {
                        if (j % 2 == 1) {
                            row   = row_vect[j];
                            coeff = 0.5 * (kt + ktmin1) * arr_vect[j];
                            col   = (j + 1) * ntheta_int + i;
                            val   = -coeff / hs;

                            Au[row] -= val * u[col];
                        }
                        else {
                            row   = row_vect[j + 1];
                            col   = j * ntheta_int + i;
                            coeff = 0.5 * (kt + ktmin1) * arr_vect[j] / hs;
                            val   = -coeff;

                            Au[row] -= val * u[col];
                        }
                    }
                    if (gyro::icntl[Param::mod_pk] > 0) {
                        if (!(extrapol && smoother == 2 && j % 2 == 0)) {
                            if ((smoother == 3 && i % 2 == 0) || (smoother == 2 && i % 2 == 1)) {
                                row  = row_vect_next[j];
                                col  = (j - 1) * ntheta_int + i;
                                val  = -art_vect[j];
                                row2 = row_vect_prev[j];
                                col2 = (j + 1) * ntheta_int + i;
                                val2 = art_vect[j];

                                val = val * u[col] + val2 * u[col2];
                                Au[row] -= val;

                                // Update (j, i-1)
                                Au[row2] -= -val;
                            }
                        }
                        if (!(extrapol && smoother == 2 && j % 2 == 1)) {
                            if ((smoother == 3 && i % 2 == 1) || (smoother == 2 && i % 2 == 0)) {
                                // Update (j+1, i)
                                row  = row_vect[j + 1];
                                col  = j * ntheta_int + im;
                                val  = -art_vect[j];
                                col2 = j * ntheta_int + ip;
                                val2 = art_vect[j];

                                Au[row] -= val * u[col] + val2 * u[col2];
                            }
                        }
                    }

                    /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!  Interior nodes   !!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     */
                    for (j = delete_circles + 1; j < nr_int - 1; j++) {
                        hs     = hplus[j];
                        hsmin1 = hplus[j - 1];
                        if (smoother > 1 && !(extrapol && smoother == 2 && j % 2 == 0)) {
                            if ((smoother == 3 && i % 2 == 1) || (smoother == 2 && i % 2 == 0)) {
                                // Update (j, i)
                                row = row_vect[j];
                                // coeff  = 0.5 * (kt + ktmin1) * arr_vect[j];
                                coeff2 = 0.5 * (hs + hsmin1) * att_vect[j];
                                col    = j * ntheta_int + im;
                                val    = -coeff2 / ktmin1;
                                col2   = j * ntheta_int + ip;
                                val2   = -coeff2 / kt;

                                Au[row] -= val * u[col] + val2 * u[col2];
                            }
                            if ((smoother == 3 && i % 2 == 0) || (smoother == 2 && i % 2 == 1)) {
                                // Update (j, i+1)
                                row   = row_vect_next[j];
                                col   = j * ntheta_int + i;
                                coeff = 0.5 * (hs + hsmin1) * att_vect[j];
                                val   = -coeff / kt;

                                Au[row] -= val * u[col];

                                // Update (j, i-1)
                                row   = row_vect_prev[j];
                                col   = j * ntheta_int + i;
                                coeff = 0.5 * (hs + hsmin1) * att_vect[j];
                                val   = -coeff / ktmin1;

                                Au[row] -= val * u[col];
                            }
                        }
                        if (extrapol && smoother == 2 && i % 2 == 0) {
                            if (j % 2 == 1) {
                                // Update (j, i)
                                row   = row_vect[j];
                                coeff = 0.5 * (kt + ktmin1) * arr_vect[j];
                                col   = (j + 1) * ntheta_int + i;
                                val   = -coeff / hs;
                                col2  = (j - 1) * ntheta_int + i;
                                val2  = -coeff / hsmin1;

                                Au[row] -= val * u[col] + val2 * u[col2];
                            }
                            else {
                                // Update (j-1, i)
                                row   = row_vect[j - 1];
                                col   = j * ntheta_int + i;
                                coeff = 0.5 * (kt + ktmin1) * arr_vect[j] / hsmin1;
                                val   = -coeff;

                                Au[row] -= val * u[col];

                                // Update (j+1, i)
                                row   = row_vect[j + 1];
                                col   = j * ntheta_int + i;
                                coeff = 0.5 * (kt + ktmin1) * arr_vect[j] / hs;
                                val   = -coeff;

                                Au[row] -= val * u[col];
                            }
                        }
                        if (gyro::icntl[Param::mod_pk] > 0) {
                            if (smoother > 1 && !(extrapol && smoother == 2 && j % 2 == 1)) {
                                if ((smoother == 3 && i % 2 == 1) || (smoother == 2 && i % 2 == 0)) {
                                    // Update (j-1, i)
                                    row  = row_vect[j - 1];
                                    row2 = row_vect[j + 1];
                                    col  = j * ntheta_int + im;
                                    val  = art_vect[j];
                                    col2 = j * ntheta_int + ip;
                                    val2 = -art_vect[j];

                                    val = val * u[col] + val2 * u[col2];
                                    Au[row] -= val;

                                    // Update (j+1, i)
                                    Au[row2] -= -val;
                                }
                            }
                            if (smoother > 1 && !(extrapol && smoother == 2 && j % 2 == 0)) {
                                if ((smoother == 3 && i % 2 == 0) || (smoother == 2 && i % 2 == 1)) {
                                    // Update (j, i+1)
                                    row  = row_vect_next[j];
                                    col  = (j - 1) * ntheta_int + i;
                                    val  = -art_vect[j];
                                    col2 = (j + 1) * ntheta_int + i;
                                    val2 = art_vect[j];

                                    val = val * u[col] + val2 * u[col2];
                                    Au[row] -= val;

                                    // Update (j, i-1)
                                    row  = row_vect_prev[j];
                                    col  = (j - 1) * ntheta_int + i;
                                    val  = -art_vect[j];
                                    col2 = (j + 1) * ntheta_int + i;
                                    val2 = art_vect[j];

                                    val = val * u[col] + val2 * u[col2];
                                    Au[row] -= -val;
                                }
                            }
                        }
                    }
                    /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!        Last lines       !!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     */
                    // DB_ext updates (~~~ Interior - (j+1, i) + DB)
                    j      = nr_int - 1;
                    hs     = hplus[j];
                    hsmin1 = hplus[j - 1];
                    if (smoother > 1 && !(extrapol && smoother == 2 && j % 2 == 0)) {
                        if ((smoother == 3 && i % 2 == 1) || (smoother == 2 && i % 2 == 0)) {
                            row = row_vect[j];
                            // Contribution to middle (top) from DB
                            coeff3 = 0.5 * (hs + hsmin1) * att_vect[j];
                            col    = j * ntheta_int + im;
                            val    = -coeff3 / ktmin1;
                            col2   = j * ntheta_int + ip;
                            val2   = -coeff3 / kt;

                            Au[row] -= val * u[col] + val2 * u[col2];
                        }
                        if ((smoother == 3 && i % 2 == 0) || (smoother == 2 && i % 2 == 1)) {
                            // Update (j, i+1)
                            row   = row_vect_next[j];
                            col   = j * ntheta_int + i;
                            coeff = 0.5 * (hs + hsmin1) * att_vect[j];
                            val   = -coeff / kt;

                            Au[row] -= val * u[col];

                            // Update (j, i-1)
                            row   = row_vect_prev[j];
                            col   = j * ntheta_int + i;
                            coeff = 0.5 * (hs + hsmin1) * att_vect[j];
                            val   = -coeff / ktmin1;

                            Au[row] -= val * u[col];
                        }
                    }
                    if (extrapol && smoother == 2 && i % 2 == 0) {
                        if (j % 2 == 0) {
                            // Update (j-1, i)
                            row   = row_vect[j - 1];
                            col   = j * ntheta_int + i;
                            coeff = 0.5 * (kt + ktmin1) * arr_vect[j] / hsmin1;
                            val   = -coeff;

                            Au[row] -= val * u[col];
                        }
                        else {
                            // Update (j, i)
                            row   = row_vect[j];
                            coeff = 0.5 * (kt + ktmin1) * arr_vect[j];
                            col   = (j - 1) * ntheta_int + i;
                            val   = -coeff / hsmin1;

                            Au[row] -= val * u[col];
                        }
                    }
                    if (gyro::icntl[Param::mod_pk] > 0) {
                        if (smoother > 1 && !(extrapol && smoother == 2 && j % 2 == 1)) {
                            if ((smoother == 3 && i % 2 == 1) || (smoother == 2 && i % 2 == 0)) {
                                // Update (j-1, i)
                                row  = row_vect[j - 1];
                                col  = j * ntheta_int + im;
                                val  = art_vect[j];
                                col2 = j * ntheta_int + ip;
                                val2 = -art_vect[j];

                                Au[row] -= val * u[col] + val2 * u[col2];
                            }
                        }
                        if (smoother > 1 && !(extrapol && smoother == 2 && j % 2 == 0)) {
                            if ((smoother == 3 && i % 2 == 0) || (smoother == 2 && i % 2 == 1)) {
                                row = row_vect_next[j];
                                col = (j - 1) * ntheta_int + i;
                                val = -art_vect[j];

                                Au[row] -= val * u[col];

                                // Update (j, i-1)
                                row = row_vect_prev[j];
                                col = (j - 1) * ntheta_int + i;
                                val = -art_vect[j];

                                Au[row] -= -val * u[col];
                            }
                        }
                    }
                }
                dep_Asc_ortho_cur[i] = 1;
            }
        }

        /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                * !!!!!!!!!!!!!!!!!!!!!!!            (2))           !!!!!!!!!!!!!!!!!!!!!!
                * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                */
        for (int i = diff + 1; i < ntheta_int; i += 3) {
            {
                int im        = imoins(i, ntheta_int);
                int ip        = iplus(i, ntheta_int);
                int im2       = imoins(im, ntheta_int);
                int ip2       = iplus(ip, ntheta_int);
                int odd_i     = i % 2;
                int ind_m_odd = i, ind_p_odd = i;
                if (odd_i) {
                    ind_m_odd = im;
                    ind_p_odd = ip;
                }
#pragma omp task shared(u, Au) firstprivate(smoother, im, i, ip, im2, ip2, odd_i) depend(in                            \
                                                                                         : dep_Asc_prev[ind_m_odd])    \
    depend(in                                                                                                          \
           : dep_Asc_prev[ind_p_odd], dep_Asc_ortho_cur[im], dep_Asc_ortho_cur[ip2]) depend(out                        \
                                                                                            : dep_Asc_ortho_cur[i])
                {
                    int row, row2, col, col2;
                    int smoother_prev, smoother_cur, smoother_next;
                    double coeff, coeff2, coeff3, kt, ktmin1, hs, hsmin1, val, val2;
                    std::vector<int> row_vect_prev, row_vect, row_vect_next;
                    std::vector<int> smoother_vect_prev, smoother_vect, smoother_vect_next;
                    std::vector<double> arr_vect, arr_vect2, att_vect, art_vect;

                    int j;
                    int s_cur     = (i % 2) + 2;
                    int s_next    = ((i + 1) % 2) + 2;
                    kt            = thetaplus_per[i + 1];
                    ktmin1        = thetaplus_per[i];
                    smoother_cur  = s_cur;
                    smoother_prev = s_next;
                    smoother_next = s_next;
                    row_vect_prev = get_row_i_glob(nr, im, smoother_prev, extrapol);
                    row_vect      = get_row_i_glob(nr, i, smoother_cur, extrapol);
                    row_vect_next = get_row_i_glob(nr, ip, smoother_next, extrapol);
                    gyro::arr_att_art(r, theta[i], arr_vect, att_vect, art_vect, 0);

                    /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                         * !!!!!!!!!!! Link Circle-Radial (radial) !!!!!!!!!!
                         * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                         */
                    j = delete_circles - 1;
                    if (!(extrapol && smoother == 2 && j % 2 == 1)) {
                        hs     = hplus[j];
                        hsmin1 = hplus[j - 1];
                        if ((smoother == 3 && i % 2 == 1) || (smoother == 2 && i % 2 == 0)) {
                            row   = row_vect[j + 1];
                            col   = j * ntheta_int + i;
                            coeff = 0.5 * (kt + ktmin1) * arr_vect[j] / hs;
                            val   = -coeff;

                            Au[row] -= val * u[col];
                            if (gyro::icntl[Param::mod_pk] > 0) {
                                row  = row_vect[j + 1];
                                col  = j * ntheta_int + im;
                                val  = -art_vect[j];
                                col2 = j * ntheta_int + ip;
                                val2 = art_vect[j];

                                Au[row] -= val * u[col] + val2 * u[col2];
                            }
                        }
                    }

                    /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             * !!!!!!!!!!! Link Radial-Circle (radial) !!!!!!!!!!
             * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             */

                    j      = delete_circles;
                    hs     = hplus[j];
                    hsmin1 = hplus[j - 1];
                    if (!(extrapol && smoother == 2 && j % 2 == 0)) {
                        if ((smoother == 3 && i % 2 == 1) || (smoother == 2 && i % 2 == 0)) {
                            // Update (j, i)
                            row   = row_vect[j];
                            coeff = 0.5 * (kt + ktmin1) * arr_vect[j];
                            col   = (j - 1) * ntheta_int + i;
                            val   = -coeff / hsmin1;

                            Au[row] -= val * u[col];

                            coeff2 = 0.5 * (hs + hsmin1) * att_vect[j];
                            col    = j * ntheta_int + im;
                            val    = -coeff2 / ktmin1;
                            col2   = j * ntheta_int + ip;
                            val2   = -coeff2 / kt;

                            Au[row] -= val * u[col] + val2 * u[col2];
                        }
                        if ((smoother == 3 && i % 2 == 0) || (smoother == 2 && i % 2 == 1)) {
                            // Update (j, i+1)
                            row   = row_vect_next[j];
                            col   = j * ntheta_int + i;
                            coeff = 0.5 * (hs + hsmin1) * att_vect[j] / kt;
                            val   = -coeff;

                            Au[row] -= val * u[col];

                            // Update (j, i-1)
                            row   = row_vect_prev[j];
                            col   = j * ntheta_int + i;
                            coeff = 0.5 * (hs + hsmin1) * att_vect[j] / ktmin1;
                            val   = -coeff;

                            Au[row] -= val * u[col];
                        }
                    }
                    if (extrapol && smoother == 2 && i % 2 == 0) {
                        if (j % 2 == 1) {
                            row   = row_vect[j];
                            coeff = 0.5 * (kt + ktmin1) * arr_vect[j];
                            col   = (j + 1) * ntheta_int + i;
                            val   = -coeff / hs;

                            Au[row] -= val * u[col];
                        }
                        else {
                            row   = row_vect[j + 1];
                            col   = j * ntheta_int + i;
                            coeff = 0.5 * (kt + ktmin1) * arr_vect[j] / hs;
                            val   = -coeff;

                            Au[row] -= val * u[col];
                        }
                    }
                    if (gyro::icntl[Param::mod_pk] > 0) {
                        if (!(extrapol && smoother == 2 && j % 2 == 0)) {
                            if ((smoother == 3 && i % 2 == 0) || (smoother == 2 && i % 2 == 1)) {
                                row  = row_vect_next[j];
                                col  = (j - 1) * ntheta_int + i;
                                val  = -art_vect[j];
                                row2 = row_vect_prev[j];
                                col2 = (j + 1) * ntheta_int + i;
                                val2 = art_vect[j];

                                val = val * u[col] + val2 * u[col2];
                                Au[row] -= val;

                                // Update (j, i-1)
                                Au[row2] -= -val;
                            }
                        }
                        if (!(extrapol && smoother == 2 && j % 2 == 1)) {
                            if ((smoother == 3 && i % 2 == 1) || (smoother == 2 && i % 2 == 0)) {
                                // Update (j+1, i)
                                row  = row_vect[j + 1];
                                col  = j * ntheta_int + im;
                                val  = -art_vect[j];
                                col2 = j * ntheta_int + ip;
                                val2 = art_vect[j];

                                Au[row] -= val * u[col] + val2 * u[col2];
                            }
                        }
                    }

                    /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!  Interior nodes   !!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     */
                    for (j = delete_circles + 1; j < nr_int - 1; j++) {
                        hs     = hplus[j];
                        hsmin1 = hplus[j - 1];
                        if (smoother > 1 && !(extrapol && smoother == 2 && j % 2 == 0)) {
                            if ((smoother == 3 && i % 2 == 1) || (smoother == 2 && i % 2 == 0)) {
                                // Update (j, i)
                                row = row_vect[j];
                                // coeff  = 0.5 * (kt + ktmin1) * arr_vect[j];
                                coeff2 = 0.5 * (hs + hsmin1) * att_vect[j];
                                col    = j * ntheta_int + im;
                                val    = -coeff2 / ktmin1;
                                col2   = j * ntheta_int + ip;
                                val2   = -coeff2 / kt;

                                Au[row] -= val * u[col] + val2 * u[col2];
                            }
                            if ((smoother == 3 && i % 2 == 0) || (smoother == 2 && i % 2 == 1)) {
                                // Update (j, i+1)
                                row   = row_vect_next[j];
                                col   = j * ntheta_int + i;
                                coeff = 0.5 * (hs + hsmin1) * att_vect[j];
                                val   = -coeff / kt;

                                Au[row] -= val * u[col];

                                // Update (j, i-1)
                                row   = row_vect_prev[j];
                                col   = j * ntheta_int + i;
                                coeff = 0.5 * (hs + hsmin1) * att_vect[j];
                                val   = -coeff / ktmin1;

                                Au[row] -= val * u[col];
                            }
                        }
                        if (extrapol && smoother == 2 && i % 2 == 0) {
                            if (j % 2 == 1) {
                                // Update (j, i)
                                row   = row_vect[j];
                                coeff = 0.5 * (kt + ktmin1) * arr_vect[j];
                                col   = (j + 1) * ntheta_int + i;
                                val   = -coeff / hs;
                                col2  = (j - 1) * ntheta_int + i;
                                val2  = -coeff / hsmin1;

                                Au[row] -= val * u[col] + val2 * u[col2];
                            }
                            else {
                                // Update (j-1, i)
                                row   = row_vect[j - 1];
                                col   = j * ntheta_int + i;
                                coeff = 0.5 * (kt + ktmin1) * arr_vect[j] / hsmin1;
                                val   = -coeff;

                                Au[row] -= val * u[col];

                                // Update (j+1, i)
                                row   = row_vect[j + 1];
                                col   = j * ntheta_int + i;
                                coeff = 0.5 * (kt + ktmin1) * arr_vect[j] / hs;
                                val   = -coeff;

                                Au[row] -= val * u[col];
                            }
                        }
                        if (gyro::icntl[Param::mod_pk] > 0) {
                            if (smoother > 1 && !(extrapol && smoother == 2 && j % 2 == 1)) {
                                if ((smoother == 3 && i % 2 == 1) || (smoother == 2 && i % 2 == 0)) {
                                    // Update (j-1, i)
                                    row  = row_vect[j - 1];
                                    row2 = row_vect[j + 1];
                                    col  = j * ntheta_int + im;
                                    val  = art_vect[j];
                                    col2 = j * ntheta_int + ip;
                                    val2 = -art_vect[j];

                                    val = val * u[col] + val2 * u[col2];
                                    Au[row] -= val;

                                    // Update (j+1, i)
                                    Au[row2] -= -val;
                                }
                            }
                            if (smoother > 1 && !(extrapol && smoother == 2 && j % 2 == 0)) {
                                if ((smoother == 3 && i % 2 == 0) || (smoother == 2 && i % 2 == 1)) {
                                    // Update (j, i+1)
                                    row  = row_vect_next[j];
                                    col  = (j - 1) * ntheta_int + i;
                                    val  = -art_vect[j];
                                    col2 = (j + 1) * ntheta_int + i;
                                    val2 = art_vect[j];

                                    val = val * u[col] + val2 * u[col2];
                                    Au[row] -= val;

                                    // Update (j, i-1)
                                    row  = row_vect_prev[j];
                                    col  = (j - 1) * ntheta_int + i;
                                    val  = -art_vect[j];
                                    col2 = (j + 1) * ntheta_int + i;
                                    val2 = art_vect[j];

                                    val = val * u[col] + val2 * u[col2];
                                    Au[row] -= -val;
                                }
                            }
                        }
                    }
                    /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!        Last lines       !!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     */
                    // DB_ext updates (~~~ Interior - (j+1, i) + DB)
                    j      = nr_int - 1;
                    hs     = hplus[j];
                    hsmin1 = hplus[j - 1];
                    if (smoother > 1 && !(extrapol && smoother == 2 && j % 2 == 0)) {
                        if ((smoother == 3 && i % 2 == 1) || (smoother == 2 && i % 2 == 0)) {
                            row = row_vect[j];
                            // Contribution to middle (top) from DB
                            coeff3 = 0.5 * (hs + hsmin1) * att_vect[j];
                            col    = j * ntheta_int + im;
                            val    = -coeff3 / ktmin1;
                            col2   = j * ntheta_int + ip;
                            val2   = -coeff3 / kt;

                            Au[row] -= val * u[col] + val2 * u[col2];
                        }
                        if ((smoother == 3 && i % 2 == 0) || (smoother == 2 && i % 2 == 1)) {
                            // Update (j, i+1)
                            row   = row_vect_next[j];
                            col   = j * ntheta_int + i;
                            coeff = 0.5 * (hs + hsmin1) * att_vect[j];
                            val   = -coeff / kt;

                            Au[row] -= val * u[col];

                            // Update (j, i-1)
                            row   = row_vect_prev[j];
                            col   = j * ntheta_int + i;
                            coeff = 0.5 * (hs + hsmin1) * att_vect[j];
                            val   = -coeff / ktmin1;

                            Au[row] -= val * u[col];
                        }
                    }
                    if (extrapol && smoother == 2 && i % 2 == 0) {
                        if (j % 2 == 0) {
                            // Update (j-1, i)
                            row   = row_vect[j - 1];
                            col   = j * ntheta_int + i;
                            coeff = 0.5 * (kt + ktmin1) * arr_vect[j] / hsmin1;
                            val   = -coeff;

                            Au[row] -= val * u[col];
                        }
                        else {
                            // Update (j, i)
                            row   = row_vect[j];
                            coeff = 0.5 * (kt + ktmin1) * arr_vect[j];
                            col   = (j - 1) * ntheta_int + i;
                            val   = -coeff / hsmin1;

                            Au[row] -= val * u[col];
                        }
                    }
                    if (gyro::icntl[Param::mod_pk] > 0) {
                        if (smoother > 1 && !(extrapol && smoother == 2 && j % 2 == 1)) {
                            if ((smoother == 3 && i % 2 == 1) || (smoother == 2 && i % 2 == 0)) {
                                // Update (j-1, i)
                                row  = row_vect[j - 1];
                                col  = j * ntheta_int + im;
                                val  = art_vect[j];
                                col2 = j * ntheta_int + ip;
                                val2 = -art_vect[j];

                                Au[row] -= val * u[col] + val2 * u[col2];
                            }
                        }
                        if (smoother > 1 && !(extrapol && smoother == 2 && j % 2 == 0)) {
                            if ((smoother == 3 && i % 2 == 0) || (smoother == 2 && i % 2 == 1)) {
                                row = row_vect_next[j];
                                col = (j - 1) * ntheta_int + i;
                                val = -art_vect[j];

                                Au[row] -= val * u[col];

                                // Update (j, i-1)
                                row = row_vect_prev[j];
                                col = (j - 1) * ntheta_int + i;
                                val = -art_vect[j];

                                Au[row] -= -val * u[col];
                            }
                        }
                    }
                }
                dep_Asc_ortho_cur[i] = 1;
            }
        }

        /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                * !!!!!!!!!!!!!!!!!!!!!!!            (3))           !!!!!!!!!!!!!!!!!!!!!!
                * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                */
        for (int i = diff + 2; i < ntheta_int; i += 3) {
            {
                int im        = i - 1;
                int ip        = iplus(i, ntheta_int);
                int im2       = imoins(im, ntheta_int);
                int ip2       = iplus(ip, ntheta_int);
                int odd_i     = i % 2;
                int ind_m_odd = i, ind_p_odd = i;
                if (odd_i) {
                    ind_m_odd = im;
                    ind_p_odd = ip;
                }
#pragma omp task shared(u, Au) firstprivate(smoother, im, i, ip, im2, ip2, odd_i) depend(in                            \
                                                                                         : dep_Asc_prev[ind_m_odd])    \
    depend(in                                                                                                          \
           : dep_Asc_prev[ind_p_odd], dep_Asc_ortho_cur[im], dep_Asc_ortho_cur[ip2]) depend(out                        \
                                                                                            : dep_Asc_ortho_cur[i])
                {
                    int row, row2, col, col2;
                    int smoother_prev, smoother_cur, smoother_next;
                    double coeff, coeff2, coeff3, kt, ktmin1, hs, hsmin1, val, val2;
                    std::vector<int> row_vect_prev, row_vect, row_vect_next;
                    std::vector<int> smoother_vect_prev, smoother_vect, smoother_vect_next;
                    std::vector<double> arr_vect, arr_vect2, att_vect, art_vect;

                    int j;
                    int s_cur     = (i % 2) + 2;
                    int s_next    = ((i + 1) % 2) + 2;
                    kt            = thetaplus_per[i + 1];
                    ktmin1        = thetaplus_per[i];
                    smoother_cur  = s_cur;
                    smoother_prev = s_next;
                    smoother_next = s_next;
                    row_vect_prev = get_row_i_glob(nr, im, smoother_prev, extrapol);
                    row_vect      = get_row_i_glob(nr, i, smoother_cur, extrapol);
                    row_vect_next = get_row_i_glob(nr, ip, smoother_next, extrapol);
                    gyro::arr_att_art(r, theta[i], arr_vect, att_vect, art_vect, 0);

                    /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                         * !!!!!!!!!!! Link Circle-Radial (radial) !!!!!!!!!!
                         * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                         */
                    j = delete_circles - 1;
                    if (!(extrapol && smoother == 2 && j % 2 == 1)) {
                        hs     = hplus[j];
                        hsmin1 = hplus[j - 1];
                        if ((smoother == 3 && i % 2 == 1) || (smoother == 2 && i % 2 == 0)) {
                            row   = row_vect[j + 1];
                            col   = j * ntheta_int + i;
                            coeff = 0.5 * (kt + ktmin1) * arr_vect[j] / hs;
                            val   = -coeff;

                            Au[row] -= val * u[col];
                            if (gyro::icntl[Param::mod_pk] > 0) {
                                row  = row_vect[j + 1];
                                col  = j * ntheta_int + im;
                                val  = -art_vect[j];
                                col2 = j * ntheta_int + ip;
                                val2 = art_vect[j];

                                Au[row] -= val * u[col] + val2 * u[col2];
                            }
                        }
                    }

                    /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             * !!!!!!!!!!! Link Radial-Circle (radial) !!!!!!!!!!
             * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             */

                    j      = delete_circles;
                    hs     = hplus[j];
                    hsmin1 = hplus[j - 1];
                    if (!(extrapol && smoother == 2 && j % 2 == 0)) {
                        if ((smoother == 3 && i % 2 == 1) || (smoother == 2 && i % 2 == 0)) {
                            // Update (j, i)
                            row   = row_vect[j];
                            coeff = 0.5 * (kt + ktmin1) * arr_vect[j];
                            col   = (j - 1) * ntheta_int + i;
                            val   = -coeff / hsmin1;

                            Au[row] -= val * u[col];

                            coeff2 = 0.5 * (hs + hsmin1) * att_vect[j];
                            col    = j * ntheta_int + im;
                            val    = -coeff2 / ktmin1;
                            col2   = j * ntheta_int + ip;
                            val2   = -coeff2 / kt;

                            Au[row] -= val * u[col] + val2 * u[col2];
                        }
                        if ((smoother == 3 && i % 2 == 0) || (smoother == 2 && i % 2 == 1)) {
                            // Update (j, i+1)
                            row   = row_vect_next[j];
                            col   = j * ntheta_int + i;
                            coeff = 0.5 * (hs + hsmin1) * att_vect[j] / kt;
                            val   = -coeff;

                            Au[row] -= val * u[col];

                            // Update (j, i-1)
                            row   = row_vect_prev[j];
                            col   = j * ntheta_int + i;
                            coeff = 0.5 * (hs + hsmin1) * att_vect[j] / ktmin1;
                            val   = -coeff;

                            Au[row] -= val * u[col];
                        }
                    }
                    if (extrapol && smoother == 2 && i % 2 == 0) {
                        if (j % 2 == 1) {
                            row   = row_vect[j];
                            coeff = 0.5 * (kt + ktmin1) * arr_vect[j];
                            col   = (j + 1) * ntheta_int + i;
                            val   = -coeff / hs;

                            Au[row] -= val * u[col];
                        }
                        else {
                            row   = row_vect[j + 1];
                            col   = j * ntheta_int + i;
                            coeff = 0.5 * (kt + ktmin1) * arr_vect[j] / hs;
                            val   = -coeff;

                            Au[row] -= val * u[col];
                        }
                    }
                    if (gyro::icntl[Param::mod_pk] > 0) {
                        if (!(extrapol && smoother == 2 && j % 2 == 0)) {
                            if ((smoother == 3 && i % 2 == 0) || (smoother == 2 && i % 2 == 1)) {
                                row  = row_vect_next[j];
                                col  = (j - 1) * ntheta_int + i;
                                val  = -art_vect[j];
                                row2 = row_vect_prev[j];
                                col2 = (j + 1) * ntheta_int + i;
                                val2 = art_vect[j];

                                val = val * u[col] + val2 * u[col2];
                                Au[row] -= val;

                                // Update (j, i-1)
                                Au[row2] -= -val;
                            }
                        }
                        if (!(extrapol && smoother == 2 && j % 2 == 1)) {
                            if ((smoother == 3 && i % 2 == 1) || (smoother == 2 && i % 2 == 0)) {
                                // Update (j+1, i)
                                row  = row_vect[j + 1];
                                col  = j * ntheta_int + im;
                                val  = -art_vect[j];
                                col2 = j * ntheta_int + ip;
                                val2 = art_vect[j];

                                Au[row] -= val * u[col] + val2 * u[col2];
                            }
                        }
                    }

                    /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!  Interior nodes   !!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     */
                    for (j = delete_circles + 1; j < nr_int - 1; j++) {
                        hs     = hplus[j];
                        hsmin1 = hplus[j - 1];
                        if (smoother > 1 && !(extrapol && smoother == 2 && j % 2 == 0)) {
                            if ((smoother == 3 && i % 2 == 1) || (smoother == 2 && i % 2 == 0)) {
                                // Update (j, i)
                                row = row_vect[j];
                                // coeff  = 0.5 * (kt + ktmin1) * arr_vect[j];
                                coeff2 = 0.5 * (hs + hsmin1) * att_vect[j];
                                col    = j * ntheta_int + im;
                                val    = -coeff2 / ktmin1;
                                col2   = j * ntheta_int + ip;
                                val2   = -coeff2 / kt;

                                Au[row] -= val * u[col] + val2 * u[col2];
                            }
                            if ((smoother == 3 && i % 2 == 0) || (smoother == 2 && i % 2 == 1)) {
                                // Update (j, i+1)
                                row   = row_vect_next[j];
                                col   = j * ntheta_int + i;
                                coeff = 0.5 * (hs + hsmin1) * att_vect[j];
                                val   = -coeff / kt;

                                Au[row] -= val * u[col];

                                // Update (j, i-1)
                                row   = row_vect_prev[j];
                                col   = j * ntheta_int + i;
                                coeff = 0.5 * (hs + hsmin1) * att_vect[j];
                                val   = -coeff / ktmin1;

                                Au[row] -= val * u[col];
                            }
                        }
                        if (extrapol && smoother == 2 && i % 2 == 0) {
                            if (j % 2 == 1) {
                                // Update (j, i)
                                row   = row_vect[j];
                                coeff = 0.5 * (kt + ktmin1) * arr_vect[j];
                                col   = (j + 1) * ntheta_int + i;
                                val   = -coeff / hs;
                                col2  = (j - 1) * ntheta_int + i;
                                val2  = -coeff / hsmin1;

                                Au[row] -= val * u[col] + val2 * u[col2];
                            }
                            else {
                                // Update (j-1, i)
                                row   = row_vect[j - 1];
                                col   = j * ntheta_int + i;
                                coeff = 0.5 * (kt + ktmin1) * arr_vect[j] / hsmin1;
                                val   = -coeff;

                                Au[row] -= val * u[col];

                                // Update (j+1, i)
                                row   = row_vect[j + 1];
                                col   = j * ntheta_int + i;
                                coeff = 0.5 * (kt + ktmin1) * arr_vect[j] / hs;
                                val   = -coeff;

                                Au[row] -= val * u[col];
                            }
                        }
                        if (gyro::icntl[Param::mod_pk] > 0) {
                            if (smoother > 1 && !(extrapol && smoother == 2 && j % 2 == 1)) {
                                if ((smoother == 3 && i % 2 == 1) || (smoother == 2 && i % 2 == 0)) {
                                    // Update (j-1, i)
                                    row  = row_vect[j - 1];
                                    row2 = row_vect[j + 1];
                                    col  = j * ntheta_int + im;
                                    val  = art_vect[j];
                                    col2 = j * ntheta_int + ip;
                                    val2 = -art_vect[j];

                                    val = val * u[col] + val2 * u[col2];
                                    Au[row] -= val;

                                    // Update (j+1, i)
                                    Au[row2] -= -val;
                                }
                            }
                            if (smoother > 1 && !(extrapol && smoother == 2 && j % 2 == 0)) {
                                if ((smoother == 3 && i % 2 == 0) || (smoother == 2 && i % 2 == 1)) {
                                    // Update (j, i+1)
                                    row  = row_vect_next[j];
                                    col  = (j - 1) * ntheta_int + i;
                                    val  = -art_vect[j];
                                    col2 = (j + 1) * ntheta_int + i;
                                    val2 = art_vect[j];

                                    val = val * u[col] + val2 * u[col2];
                                    Au[row] -= val;

                                    // Update (j, i-1)
                                    row  = row_vect_prev[j];
                                    col  = (j - 1) * ntheta_int + i;
                                    val  = -art_vect[j];
                                    col2 = (j + 1) * ntheta_int + i;
                                    val2 = art_vect[j];

                                    val = val * u[col] + val2 * u[col2];
                                    Au[row] -= -val;
                                }
                            }
                        }
                    }
                    /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!        Last lines       !!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     */
                    // DB_ext updates (~~~ Interior - (j+1, i) + DB)
                    j      = nr_int - 1;
                    hs     = hplus[j];
                    hsmin1 = hplus[j - 1];
                    if (smoother > 1 && !(extrapol && smoother == 2 && j % 2 == 0)) {
                        if ((smoother == 3 && i % 2 == 1) || (smoother == 2 && i % 2 == 0)) {
                            row = row_vect[j];
                            // Contribution to middle (top) from DB
                            coeff3 = 0.5 * (hs + hsmin1) * att_vect[j];
                            col    = j * ntheta_int + im;
                            val    = -coeff3 / ktmin1;
                            col2   = j * ntheta_int + ip;
                            val2   = -coeff3 / kt;

                            Au[row] -= val * u[col] + val2 * u[col2];
                        }
                        if ((smoother == 3 && i % 2 == 0) || (smoother == 2 && i % 2 == 1)) {
                            // Update (j, i+1)
                            row   = row_vect_next[j];
                            col   = j * ntheta_int + i;
                            coeff = 0.5 * (hs + hsmin1) * att_vect[j];
                            val   = -coeff / kt;

                            Au[row] -= val * u[col];

                            // Update (j, i-1)
                            row   = row_vect_prev[j];
                            col   = j * ntheta_int + i;
                            coeff = 0.5 * (hs + hsmin1) * att_vect[j];
                            val   = -coeff / ktmin1;

                            Au[row] -= val * u[col];
                        }
                    }
                    if (extrapol && smoother == 2 && i % 2 == 0) {
                        if (j % 2 == 0) {
                            // Update (j-1, i)
                            row   = row_vect[j - 1];
                            col   = j * ntheta_int + i;
                            coeff = 0.5 * (kt + ktmin1) * arr_vect[j] / hsmin1;
                            val   = -coeff;

                            Au[row] -= val * u[col];
                        }
                        else {
                            // Update (j, i)
                            row   = row_vect[j];
                            coeff = 0.5 * (kt + ktmin1) * arr_vect[j];
                            col   = (j - 1) * ntheta_int + i;
                            val   = -coeff / hsmin1;

                            Au[row] -= val * u[col];
                        }
                    }
                    if (gyro::icntl[Param::mod_pk] > 0) {
                        if (smoother > 1 && !(extrapol && smoother == 2 && j % 2 == 1)) {
                            if ((smoother == 3 && i % 2 == 1) || (smoother == 2 && i % 2 == 0)) {
                                // Update (j-1, i)
                                row  = row_vect[j - 1];
                                col  = j * ntheta_int + im;
                                val  = art_vect[j];
                                col2 = j * ntheta_int + ip;
                                val2 = -art_vect[j];

                                Au[row] -= val * u[col] + val2 * u[col2];
                            }
                        }
                        if (smoother > 1 && !(extrapol && smoother == 2 && j % 2 == 0)) {
                            if ((smoother == 3 && i % 2 == 0) || (smoother == 2 && i % 2 == 1)) {
                                row = row_vect_next[j];
                                col = (j - 1) * ntheta_int + i;
                                val = -art_vect[j];

                                Au[row] -= val * u[col];

                                // Update (j, i-1)
                                row = row_vect_prev[j];
                                col = (j - 1) * ntheta_int + i;
                                val = -art_vect[j];

                                Au[row] -= -val * u[col];
                            }
                        }
                    }
                }
                dep_Asc_ortho_cur[i] = 1;
            }
        }
    }
} /* ----- end of level::apply_Asc_ortho ----- */
