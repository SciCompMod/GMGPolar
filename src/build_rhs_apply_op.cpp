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
 * \file build_rhs_apply_op.cpp
 * \brief Implementation of the 9p FD operators A and RHS
 * \author M. Kuehn, C. Kruse, P. Leleux, C. Schwarz
 * \version 0.0
 */
#include "level.h"

/*!
 *  \brief Build the operator A
 *
 * Builds the matrix A of the system based on 9p FD. Relations with Dirichlet 
 * boundary condition nodes are shifted to the RHS to keep a symmetric operator.
 *
 * everything expressed in polar coordinates
 * allows anisotropy in r
 * - r: array of discretized r values
 * - ntheta: number of discretization intervals in angle phi (also variable m)
 *
 * solves the transformed DGL:  r*u_rr + u_r + 1/r*u_phiphi = rf fuer r>0, r<=1,
 * fuer r=0, approximate: u_0 = 1/m sum_j u(h,jk)-(h/2)^2*f
 * sorting of nodes ringwise. first node origin, last ntheta
 * Node all the way out
 */
void level::build_A()
{
    int start_j;
    int* dep = new int[nr];

    /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!       DB first line      !!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     */
    // Take boundary condition into account: Dirichlet-RB
    if (gyro::icntl[Param::DirBC_Interior]) { // (r[0],0) is on Dirichlet boundary
        start_j = 1;
    }
    else { // (r[0],0) is not on Dirichlet boundary
        start_j = 0;
    }

#pragma omp parallel shared(dep)
    {
#pragma omp single
        {
            if (gyro::icntl[Param::DirBC_Interior]) { // (r[0],0) is on Dirichlet boundary
#pragma omp task shared(dep, start_j) depend(out : dep[0])
                {
                    for (int i = 0; i < ntheta_int; i++) {
                        row_indices[i] = i;
                        col_indices[i] = i;

                        vals[i] += 1.0;
                    }
                } //end of task and parallel
            }

#pragma omp task shared(dep, start_j) depend(out : dep[nr_int])
            {
                // int ptr, row; // To be removed ?
                // Take boundary condition into account: Dirichlet-RB
                for (int i = 0; i < ntheta_int; i++) {
                    // row                              = m - ntheta_int + i;
                    // ptr                              = nz - ntheta_int + i;
                    row_indices[nz - ntheta_int + i] = m - ntheta_int + i;
                    col_indices[nz - ntheta_int + i] = m - ntheta_int + i;

                    vals[nz - ntheta_int + i] += 1.0;
                }
            } //end of task and parallel

#pragma omp task shared(dep, start_j) depend(out : dep[start_j])
            {
                int i, j, ptr, row, col;
                double coeff, coeff2, val, kt, ktmin1, hs, hsmin1;
                std::vector<int> ptr_vect_prev, ptr_vect, ptr_vect_next;
                std::vector<int> stencil_prev, stencil_cur, stencil_next, stencil;
                std::vector<double> arr_vect, arr_vect2, att_vect, art_vect;
                /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!       First lines       !!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     */
                j             = start_j;
                ptr_vect      = get_ptr(j);
                ptr_vect_next = get_ptr(j + 1);
                stencil_cur   = get_stencil(j);
                stencil_next  = get_stencil(j + 1);
                stencil       = stencil_cur;

                for (int i = 0; i < ntheta_int; i++) {
                    ptr = ptr_vect[i + 1];
                    row = j * ntheta_int + i;
                    val = betaVec[row];

                    vals[ptr + stencil[Param::middle]] += val;
                }

                hs       = hplus[j];
                arr_vect = gyro::arr(r[j], theta, sin_theta, cos_theta, ntheta_int, 0);
                // Across: bottom update
                if (!gyro::icntl[Param::DirBC_Interior]) {
                    hsmin1 = 2 * r[0];
                    // Accross the origin theta and arr
                    arr_vect2 = gyro::arr(r[j], theta_PI, sin_theta_PI, cos_theta_PI, ntheta_int, 0);

                    for (i = 0; i < ntheta_int; i++) {
                        // Define the index, position and interval size of current and previous node
                        // - in theta
                        kt     = thetaplus_per[i + 1];
                        ktmin1 = thetaplus_per[i];

                        ptr = ptr_vect[i + 1];
                        row = j * ntheta_int + i;
                        if (i + 1 <= ntheta_int / 2) // first half of circle: half a turn further
                            col = row + ntheta_int / 2;
                        else // second half of circle: half a turn back
                            col = row - ntheta_int / 2;
                        row_indices[ptr + stencil[Param::bottom]] = row;
                        col_indices[ptr + stencil[Param::bottom]] = col;

                        coeff  = 0.5 * (kt + ktmin1) * arr_vect2[i] / hsmin1;
                        coeff2 = 0.5 * (kt + ktmin1) * arr_vect[i] / hsmin1;

                        vals[ptr + stencil[Param::bottom]] += -coeff - coeff2;

                        vals[ptr + stencil[Param::middle]] += coeff;
                    }
                }
                else {
                    hsmin1 = hplus[j - 1];
                    // DB contribution arr (r(0))
                    arr_vect2 = gyro::arr(r[j - 1], theta, sin_theta, cos_theta, ntheta_int, 0);

                    for (i = 0; i < ntheta_int; i++) {
                        // Define the index, position and interval size of current and previous node
                        // - in theta
                        kt     = thetaplus_per[i + 1];
                        ktmin1 = thetaplus_per[i];

                        ptr   = ptr_vect[i + 1];
                        coeff = 0.5 * (kt + ktmin1) * arr_vect2[i] / hsmin1;
                        row   = j * ntheta_int + i;

                        vals[ptr + stencil[Param::middle]] += coeff;
                    }
                }

                // Across and DB_int updates (~~~ Interior - (j-1, i))
                gyro::arr_att_art(r[j], theta, sin_theta, cos_theta, ntheta_int, arr_vect, att_vect, art_vect, 0);
                for (i = 0; i < ntheta_int; i++) {
                    // Define the index, position and interval size of current and previous node
                    // - in theta
                    kt     = thetaplus_per[i + 1];
                    ktmin1 = thetaplus_per[i];

                    // Update (j, i) (=== Interior - bottom)
                    ptr                                       = ptr_vect[i + 1];
                    stencil                                   = stencil_cur;
                    row                                       = j * ntheta_int + i;
                    row_indices[ptr + stencil[Param::middle]] = row;
                    col_indices[ptr + stencil[Param::middle]] = row;

                    coeff = 0.5 * (kt + ktmin1) * arr_vect[i];
                    col   = (j + 1) * ntheta_int + i;

                    vals[ptr + stencil[Param::top]] += -coeff / hs;

                    coeff2 = 0.5 * (hs + hsmin1) * att_vect[i];
                    col    = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);

                    vals[ptr + stencil[Param::left]] += -coeff2 / ktmin1;
                    col = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);

                    vals[ptr + stencil[Param::right]] += -coeff2 / kt;

                    vals[ptr + stencil[Param::middle]] += coeff / hs + coeff / hsmin1 + coeff2 / ktmin1 + coeff2 / kt;

                    // Update (j, i+1) (=== Interior - bottom_left)
                    ptr                                     = ptr_vect[i + 2];
                    row                                     = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);
                    col                                     = j * ntheta_int + i;
                    row_indices[ptr + stencil[Param::left]] = row;
                    col_indices[ptr + stencil[Param::left]] = col;

                    coeff = 0.5 * (hs + hsmin1) * att_vect[i] / kt;

                    vals[ptr + stencil[Param::left]] += -coeff;

                    vals[ptr + stencil[Param::middle]] += coeff;

                    // Update (j, i-1) (=== Interior - bottom_right)
                    ptr                                      = ptr_vect[i];
                    row                                      = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);
                    col                                      = j * ntheta_int + i;
                    row_indices[ptr + stencil[Param::right]] = row;
                    col_indices[ptr + stencil[Param::right]] = col;

                    coeff = 0.5 * (hs + hsmin1) * att_vect[i] / ktmin1;

                    vals[ptr + stencil[Param::right]] += -coeff;

                    vals[ptr + stencil[Param::middle]] += coeff;

                    // Update (j+1, i) (=== Interior)
                    ptr                                       = ptr_vect_next[i + 1];
                    stencil                                   = stencil_next;
                    row                                       = (j + 1) * ntheta_int + i;
                    col                                       = j * ntheta_int + i;
                    row_indices[ptr + stencil[Param::bottom]] = row;
                    col_indices[ptr + stencil[Param::bottom]] = col;

                    coeff = 0.5 * (kt + ktmin1) * arr_vect[i] / hs;

                    vals[ptr + stencil[Param::bottom]] += -coeff;

                    vals[ptr + stencil[Param::middle]] += coeff;
                }
                if (gyro::icntl[Param::mod_pk] > 0)
                    for (i = 0; i < ntheta_int; i++) {
                        // Define the index, position and interval size of current and previous node
                        // - in theta
                        kt     = thetaplus_per[i + 1];
                        ktmin1 = thetaplus_per[i];

                        // Update (j, i+1) (=== Interior - bottom_left)
                        ptr     = ptr_vect[i + 2];
                        stencil = stencil_cur;
                        row     = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);
                        col     = (j + 1) * ntheta_int + i;
                        row_indices[ptr + stencil[Param::top_left]] = row;
                        col_indices[ptr + stencil[Param::top_left]] = col;

                        vals[ptr + stencil[Param::top_left]] += art_vect[i];

                        // Update (j, i-1) (=== Interior - bottom_right)
                        ptr = ptr_vect[i];
                        row = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);
                        col = (j + 1) * ntheta_int + i;
                        row_indices[ptr + stencil[Param::top_right]] = row;
                        col_indices[ptr + stencil[Param::top_right]] = col;

                        vals[ptr + stencil[Param::top_right]] += -art_vect[i];

                        // Update (j+1, i) (=== Interior)
                        ptr     = ptr_vect_next[i + 1];
                        stencil = stencil_next;
                        row     = (j + 1) * ntheta_int + i;
                        col     = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);
                        row_indices[ptr + stencil[Param::bottom_left]] = row;
                        col_indices[ptr + stencil[Param::bottom_left]] = col;

                        vals[ptr + stencil[Param::bottom_left]] += -art_vect[i];
                        col = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);
                        row_indices[ptr + stencil[Param::bottom_right]] = row;
                        col_indices[ptr + stencil[Param::bottom_right]] = col;

                        vals[ptr + stencil[Param::bottom_right]] += art_vect[i];
                    }
            } // end of task

            /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!      Interior nodes (1)      !!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     */
            for (int j = start_j + 3; j < nr_int - 1; j += 3) {
#pragma omp task shared(dep, start_j) firstprivate(j) depend(out : dep[j])
                {
                    int ptr, row, col;
                    double coeff, coeff2, val, kt, ktmin1, hs, hsmin1;
                    std::vector<int> ptr_vect_prev, ptr_vect, ptr_vect_next;
                    std::vector<int> stencil_prev, stencil_cur, stencil_next, stencil;
                    std::vector<double> arr_vect, arr_vect2, att_vect, art_vect;

                    ptr_vect_prev = get_ptr(j - 1);
                    ptr_vect      = get_ptr(j);
                    ptr_vect_next = get_ptr(j + 1);
                    stencil_prev  = get_stencil(j - 1);
                    stencil_cur   = get_stencil(j);
                    stencil_next  = get_stencil(j + 1);

                    stencil = stencil_cur;
                    for (int i = 0; i < ntheta_int; i++) {
                        ptr = ptr_vect[i + 1];
                        row = j * ntheta_int + i;
                        val = betaVec[row];

                        vals[ptr + stencil[Param::middle]] += val;
                    }

                    // - in r
                    hs     = hplus[j];
                    hsmin1 = hplus[j - 1];
                    gyro::arr_att_art(r[j], theta, sin_theta, cos_theta, ntheta_int, arr_vect, att_vect, art_vect, 0);

                    for (int i = 0; i < ntheta_int; i++) {
                        // Define the index, position and interval size of current and previous node
                        // - in theta
                        kt     = thetaplus_per[i + 1];
                        ktmin1 = thetaplus_per[i];

                        // Update (j-1, i) (Not in DB_int)
                        ptr                                    = ptr_vect_prev[i + 1];
                        stencil                                = stencil_prev;
                        row                                    = (j - 1) * ntheta_int + i;
                        col                                    = j * ntheta_int + i;
                        row_indices[ptr + stencil[Param::top]] = row;
                        col_indices[ptr + stencil[Param::top]] = col;

                        coeff = 0.5 * (kt + ktmin1) * arr_vect[i] / hsmin1;

                        vals[ptr + stencil[Param::top]] += -coeff;

                        vals[ptr + stencil[Param::middle]] += coeff;

                        // Update (j, i)
                        ptr                                       = ptr_vect[i + 1];
                        stencil                                   = stencil_cur;
                        row                                       = j * ntheta_int + i;
                        row_indices[ptr + stencil[Param::middle]] = row;
                        col_indices[ptr + stencil[Param::middle]] = row;

                        coeff = 0.5 * (kt + ktmin1) * arr_vect[i];
                        col   = (j + 1) * ntheta_int + i;

                        vals[ptr + stencil[Param::top]] += -coeff / hs;
                        col = (j - 1) * ntheta_int + i;

                        vals[ptr + stencil[Param::bottom]] += -coeff / hsmin1;

                        coeff2 = 0.5 * (hs + hsmin1) * att_vect[i];
                        col    = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);

                        vals[ptr + stencil[Param::left]] += -coeff2 / ktmin1;
                        col = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);

                        vals[ptr + stencil[Param::right]] += -coeff2 / kt;

                        vals[ptr + stencil[Param::middle]] +=
                            coeff / hs + coeff / hsmin1 + coeff2 / ktmin1 + coeff2 / kt;

                        // Update (j, i+1)
                        ptr = ptr_vect[i + 2];
                        row = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);
                        col = j * ntheta_int + i;
                        row_indices[ptr + stencil[Param::left]] = row;
                        col_indices[ptr + stencil[Param::left]] = col;

                        coeff = 0.5 * (hs + hsmin1) * att_vect[i] / kt;

                        vals[ptr + stencil[Param::left]] += -coeff;

                        vals[ptr + stencil[Param::middle]] += coeff;

                        // Update (j, i-1)
                        ptr = ptr_vect[i];
                        row = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);
                        col = j * ntheta_int + i;
                        row_indices[ptr + stencil[Param::right]] = row;
                        col_indices[ptr + stencil[Param::right]] = col;

                        coeff = 0.5 * (hs + hsmin1) * att_vect[i] / ktmin1;

                        vals[ptr + stencil[Param::right]] += -coeff;

                        vals[ptr + stencil[Param::middle]] += coeff;

                        // Update (j+1, i) (Not in DB_ext)
                        ptr                                       = ptr_vect_next[i + 1];
                        stencil                                   = stencil_next;
                        row                                       = (j + 1) * ntheta_int + i;
                        col                                       = j * ntheta_int + i;
                        row_indices[ptr + stencil[Param::bottom]] = row;
                        col_indices[ptr + stencil[Param::bottom]] = col;

                        coeff = 0.5 * (kt + ktmin1) * arr_vect[i] / hs;

                        vals[ptr + stencil[Param::bottom]] += -coeff;

                        vals[ptr + stencil[Param::middle]] += coeff;
                    }
                    if (gyro::icntl[Param::mod_pk] > 0)
                        for (int i = 0; i < ntheta_int; i++) {
                            // Define the index, position and interval size of current and previous node
                            // - in theta
                            kt     = thetaplus_per[i + 1];
                            ktmin1 = thetaplus_per[i];

                            // Update (j-1, i) (Not in DB_int)
                            ptr     = ptr_vect_prev[i + 1];
                            stencil = stencil_prev;
                            row     = (j - 1) * ntheta_int + i;
                            col     = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);
                            row_indices[ptr + stencil[Param::top_left]] = row;
                            col_indices[ptr + stencil[Param::top_left]] = col;

                            vals[ptr + stencil[Param::top_left]] += art_vect[i];
                            col = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);
                            row_indices[ptr + stencil[Param::top_right]] = row;
                            col_indices[ptr + stencil[Param::top_right]] = col;

                            vals[ptr + stencil[Param::top_right]] += -art_vect[i];

                            // Update (j, i+1)
                            ptr     = ptr_vect[i + 2];
                            stencil = stencil_cur;
                            row     = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);
                            col     = (j - 1) * ntheta_int + i;
                            row_indices[ptr + stencil[Param::bottom_left]] = row;
                            col_indices[ptr + stencil[Param::bottom_left]] = col;

                            vals[ptr + stencil[Param::bottom_left]] += -art_vect[i];
                            col                                         = (j + 1) * ntheta_int + i;
                            row_indices[ptr + stencil[Param::top_left]] = row;
                            col_indices[ptr + stencil[Param::top_left]] = col;

                            vals[ptr + stencil[Param::top_left]] += art_vect[i];

                            // Update (j, i-1)
                            ptr = ptr_vect[i];
                            row = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);
                            col = (j - 1) * ntheta_int + i;
                            row_indices[ptr + stencil[Param::bottom_right]] = row;
                            col_indices[ptr + stencil[Param::bottom_right]] = col;

                            vals[ptr + stencil[Param::bottom_right]] += art_vect[i];
                            col                                          = (j + 1) * ntheta_int + i;
                            row_indices[ptr + stencil[Param::top_right]] = row;
                            col_indices[ptr + stencil[Param::top_right]] = col;

                            vals[ptr + stencil[Param::top_right]] += -art_vect[i];

                            // Update (j+1, i) (Not in DB_ext)
                            ptr     = ptr_vect_next[i + 1];
                            stencil = stencil_next;
                            row     = (j + 1) * ntheta_int + i;
                            col     = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);
                            row_indices[ptr + stencil[Param::bottom_left]] = row;
                            col_indices[ptr + stencil[Param::bottom_left]] = col;

                            vals[ptr + stencil[Param::bottom_left]] += -art_vect[i];
                            col = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);
                            row_indices[ptr + stencil[Param::bottom_right]] = row;
                            col_indices[ptr + stencil[Param::bottom_right]] = col;

                            vals[ptr + stencil[Param::bottom_right]] += art_vect[i];
                        }
                } // end of task
            }
            /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!      Interior nodes (2)      !!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     */
            for (int j = start_j + 1; j < nr_int - 1; j += 3) {
#pragma omp task shared(dep, start_j) firstprivate(j) depend(in                                                        \
                                                             : dep[j - 1]) depend(in                                   \
                                                                                  : dep[j + 2]) depend(out             \
                                                                                                       : dep[j])
                {
                    int ptr, row, col;
                    double coeff, coeff2, val, kt, ktmin1, hs, hsmin1;
                    std::vector<int> ptr_vect_prev, ptr_vect, ptr_vect_next;
                    std::vector<int> stencil_prev, stencil_cur, stencil_next, stencil;
                    std::vector<double> arr_vect, arr_vect2, att_vect, art_vect;

                    ptr_vect_prev = get_ptr(j - 1);
                    ptr_vect      = get_ptr(j);
                    ptr_vect_next = get_ptr(j + 1);
                    stencil_prev  = get_stencil(j - 1);
                    stencil_cur   = get_stencil(j);
                    stencil_next  = get_stencil(j + 1);

                    stencil = stencil_cur;
                    for (int i = 0; i < ntheta_int; i++) {
                        ptr = ptr_vect[i + 1];
                        row = j * ntheta_int + i;
                        val = betaVec[row];

                        vals[ptr + stencil[Param::middle]] += val;
                    }

                    // - in r
                    hs     = hplus[j];
                    hsmin1 = hplus[j - 1];
                    gyro::arr_att_art(r[j], theta, sin_theta, cos_theta, ntheta_int, arr_vect, att_vect, art_vect, 0);

                    for (int i = 0; i < ntheta_int; i++) {
                        // Define the index, position and interval size of current and previous node
                        // - in theta
                        kt     = thetaplus_per[i + 1];
                        ktmin1 = thetaplus_per[i];

                        // Update (j-1, i) (Not in DB_int)
                        ptr                                    = ptr_vect_prev[i + 1];
                        stencil                                = stencil_prev;
                        row                                    = (j - 1) * ntheta_int + i;
                        col                                    = j * ntheta_int + i;
                        row_indices[ptr + stencil[Param::top]] = row;
                        col_indices[ptr + stencil[Param::top]] = col;

                        coeff = 0.5 * (kt + ktmin1) * arr_vect[i] / hsmin1;

                        vals[ptr + stencil[Param::top]] += -coeff;

                        vals[ptr + stencil[Param::middle]] += coeff;

                        // Update (j, i)
                        ptr                                       = ptr_vect[i + 1];
                        stencil                                   = stencil_cur;
                        row                                       = j * ntheta_int + i;
                        row_indices[ptr + stencil[Param::middle]] = row;
                        col_indices[ptr + stencil[Param::middle]] = row;

                        coeff = 0.5 * (kt + ktmin1) * arr_vect[i];
                        col   = (j + 1) * ntheta_int + i;

                        vals[ptr + stencil[Param::top]] += -coeff / hs;
                        col = (j - 1) * ntheta_int + i;

                        vals[ptr + stencil[Param::bottom]] += -coeff / hsmin1;

                        coeff2 = 0.5 * (hs + hsmin1) * att_vect[i];
                        col    = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);

                        vals[ptr + stencil[Param::left]] += -coeff2 / ktmin1;
                        col = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);

                        vals[ptr + stencil[Param::right]] += -coeff2 / kt;

                        vals[ptr + stencil[Param::middle]] +=
                            coeff / hs + coeff / hsmin1 + coeff2 / ktmin1 + coeff2 / kt;

                        // Update (j, i+1)
                        ptr = ptr_vect[i + 2];
                        row = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);
                        col = j * ntheta_int + i;
                        row_indices[ptr + stencil[Param::left]] = row;
                        col_indices[ptr + stencil[Param::left]] = col;

                        coeff = 0.5 * (hs + hsmin1) * att_vect[i] / kt;

                        vals[ptr + stencil[Param::left]] += -coeff;

                        vals[ptr + stencil[Param::middle]] += coeff;

                        // Update (j, i-1)
                        ptr = ptr_vect[i];
                        row = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);
                        col = j * ntheta_int + i;
                        row_indices[ptr + stencil[Param::right]] = row;
                        col_indices[ptr + stencil[Param::right]] = col;

                        coeff = 0.5 * (hs + hsmin1) * att_vect[i] / ktmin1;

                        vals[ptr + stencil[Param::right]] += -coeff;

                        vals[ptr + stencil[Param::middle]] += coeff;

                        // Update (j+1, i) (Not in DB_ext)
                        ptr                                       = ptr_vect_next[i + 1];
                        stencil                                   = stencil_next;
                        row                                       = (j + 1) * ntheta_int + i;
                        col                                       = j * ntheta_int + i;
                        row_indices[ptr + stencil[Param::bottom]] = row;
                        col_indices[ptr + stencil[Param::bottom]] = col;

                        coeff = 0.5 * (kt + ktmin1) * arr_vect[i] / hs;

                        vals[ptr + stencil[Param::bottom]] += -coeff;

                        vals[ptr + stencil[Param::middle]] += coeff;
                    }
                    if (gyro::icntl[Param::mod_pk] > 0)
                        for (int i = 0; i < ntheta_int; i++) {
                            // Define the index, position and interval size of current and previous node
                            // - in theta
                            kt     = thetaplus_per[i + 1];
                            ktmin1 = thetaplus_per[i];

                            // Update (j-1, i) (Not in DB_int)
                            ptr     = ptr_vect_prev[i + 1];
                            stencil = stencil_prev;
                            row     = (j - 1) * ntheta_int + i;
                            col     = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);
                            row_indices[ptr + stencil[Param::top_left]] = row;
                            col_indices[ptr + stencil[Param::top_left]] = col;

                            vals[ptr + stencil[Param::top_left]] += art_vect[i];
                            col = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);
                            row_indices[ptr + stencil[Param::top_right]] = row;
                            col_indices[ptr + stencil[Param::top_right]] = col;

                            vals[ptr + stencil[Param::top_right]] += -art_vect[i];

                            // Update (j, i+1)
                            ptr     = ptr_vect[i + 2];
                            stencil = stencil_cur;
                            row     = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);
                            col     = (j - 1) * ntheta_int + i;
                            row_indices[ptr + stencil[Param::bottom_left]] = row;
                            col_indices[ptr + stencil[Param::bottom_left]] = col;

                            vals[ptr + stencil[Param::bottom_left]] += -art_vect[i];
                            col                                         = (j + 1) * ntheta_int + i;
                            row_indices[ptr + stencil[Param::top_left]] = row;
                            col_indices[ptr + stencil[Param::top_left]] = col;

                            vals[ptr + stencil[Param::top_left]] += art_vect[i];

                            // Update (j, i-1)
                            ptr = ptr_vect[i];
                            row = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);
                            col = (j - 1) * ntheta_int + i;
                            row_indices[ptr + stencil[Param::bottom_right]] = row;
                            col_indices[ptr + stencil[Param::bottom_right]] = col;

                            vals[ptr + stencil[Param::bottom_right]] += art_vect[i];
                            col                                          = (j + 1) * ntheta_int + i;
                            row_indices[ptr + stencil[Param::top_right]] = row;
                            col_indices[ptr + stencil[Param::top_right]] = col;

                            vals[ptr + stencil[Param::top_right]] += -art_vect[i];

                            // Update (j+1, i) (Not in DB_ext)
                            ptr     = ptr_vect_next[i + 1];
                            stencil = stencil_next;
                            row     = (j + 1) * ntheta_int + i;
                            col     = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);
                            row_indices[ptr + stencil[Param::bottom_left]] = row;
                            col_indices[ptr + stencil[Param::bottom_left]] = col;

                            vals[ptr + stencil[Param::bottom_left]] += -art_vect[i];
                            col = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);
                            row_indices[ptr + stencil[Param::bottom_right]] = row;
                            col_indices[ptr + stencil[Param::bottom_right]] = col;

                            vals[ptr + stencil[Param::bottom_right]] += art_vect[i];
                        }
                } // end of task
            }
            /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!      Interior nodes (3)      !!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     */
            for (int j = start_j + 2; j < nr_int - 1; j += 3) {
#pragma omp task shared(dep, start_j) firstprivate(j) depend(in                                                        \
                                                             : dep[j - 1]) depend(in                                   \
                                                                                  : dep[j + 2]) depend(out             \
                                                                                                       : dep[j])
                {
                    int ptr, row, col;
                    double coeff, coeff2, val, kt, ktmin1, hs, hsmin1;
                    std::vector<int> ptr_vect_prev, ptr_vect, ptr_vect_next;
                    std::vector<int> stencil_prev, stencil_cur, stencil_next, stencil;
                    std::vector<double> arr_vect, arr_vect2, att_vect, art_vect;

                    ptr_vect_prev = get_ptr(j - 1);
                    ptr_vect      = get_ptr(j);
                    ptr_vect_next = get_ptr(j + 1);
                    stencil_prev  = get_stencil(j - 1);
                    stencil_cur   = get_stencil(j);
                    stencil_next  = get_stencil(j + 1);

                    stencil = stencil_cur;
                    for (int i = 0; i < ntheta_int; i++) {
                        ptr = ptr_vect[i + 1];
                        row = j * ntheta_int + i;
                        val = betaVec[row];

                        vals[ptr + stencil[Param::middle]] += val;
                    }

                    // - in r
                    hs     = hplus[j];
                    hsmin1 = hplus[j - 1];
                    gyro::arr_att_art(r[j], theta, sin_theta, cos_theta, ntheta_int, arr_vect, att_vect, art_vect, 0);

                    for (int i = 0; i < ntheta_int; i++) {
                        // Define the index, position and interval size of current and previous node
                        // - in theta
                        kt     = thetaplus_per[i + 1];
                        ktmin1 = thetaplus_per[i];

                        // Update (j-1, i) (Not in DB_int)
                        ptr                                    = ptr_vect_prev[i + 1];
                        stencil                                = stencil_prev;
                        row                                    = (j - 1) * ntheta_int + i;
                        col                                    = j * ntheta_int + i;
                        row_indices[ptr + stencil[Param::top]] = row;
                        col_indices[ptr + stencil[Param::top]] = col;

                        coeff = 0.5 * (kt + ktmin1) * arr_vect[i] / hsmin1;

                        vals[ptr + stencil[Param::top]] += -coeff;

                        vals[ptr + stencil[Param::middle]] += coeff;

                        // Update (j, i)
                        ptr                                       = ptr_vect[i + 1];
                        stencil                                   = stencil_cur;
                        row                                       = j * ntheta_int + i;
                        row_indices[ptr + stencil[Param::middle]] = row;
                        col_indices[ptr + stencil[Param::middle]] = row;

                        coeff = 0.5 * (kt + ktmin1) * arr_vect[i];
                        col   = (j + 1) * ntheta_int + i;

                        vals[ptr + stencil[Param::top]] += -coeff / hs;
                        col = (j - 1) * ntheta_int + i;

                        vals[ptr + stencil[Param::bottom]] += -coeff / hsmin1;

                        coeff2 = 0.5 * (hs + hsmin1) * att_vect[i];
                        col    = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);

                        vals[ptr + stencil[Param::left]] += -coeff2 / ktmin1;
                        col = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);

                        vals[ptr + stencil[Param::right]] += -coeff2 / kt;

                        vals[ptr + stencil[Param::middle]] +=
                            coeff / hs + coeff / hsmin1 + coeff2 / ktmin1 + coeff2 / kt;

                        // Update (j, i+1)
                        ptr = ptr_vect[i + 2];
                        row = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);
                        col = j * ntheta_int + i;
                        row_indices[ptr + stencil[Param::left]] = row;
                        col_indices[ptr + stencil[Param::left]] = col;

                        coeff = 0.5 * (hs + hsmin1) * att_vect[i] / kt;

                        vals[ptr + stencil[Param::left]] += -coeff;

                        vals[ptr + stencil[Param::middle]] += coeff;

                        // Update (j, i-1)
                        ptr = ptr_vect[i];
                        row = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);
                        col = j * ntheta_int + i;
                        row_indices[ptr + stencil[Param::right]] = row;
                        col_indices[ptr + stencil[Param::right]] = col;

                        coeff = 0.5 * (hs + hsmin1) * att_vect[i] / ktmin1;

                        vals[ptr + stencil[Param::right]] += -coeff;

                        vals[ptr + stencil[Param::middle]] += coeff;

                        // Update (j+1, i) (Not in DB_ext)
                        ptr                                       = ptr_vect_next[i + 1];
                        stencil                                   = stencil_next;
                        row                                       = (j + 1) * ntheta_int + i;
                        col                                       = j * ntheta_int + i;
                        row_indices[ptr + stencil[Param::bottom]] = row;
                        col_indices[ptr + stencil[Param::bottom]] = col;

                        coeff = 0.5 * (kt + ktmin1) * arr_vect[i] / hs;

                        vals[ptr + stencil[Param::bottom]] += -coeff;

                        vals[ptr + stencil[Param::middle]] += coeff;
                    }
                    if (gyro::icntl[Param::mod_pk] > 0)
                        for (int i = 0; i < ntheta_int; i++) {
                            // Define the index, position and interval size of current and previous node
                            // - in theta
                            kt     = thetaplus_per[i + 1];
                            ktmin1 = thetaplus_per[i];

                            // Update (j-1, i) (Not in DB_int)
                            ptr     = ptr_vect_prev[i + 1];
                            stencil = stencil_prev;
                            row     = (j - 1) * ntheta_int + i;
                            col     = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);
                            row_indices[ptr + stencil[Param::top_left]] = row;
                            col_indices[ptr + stencil[Param::top_left]] = col;

                            vals[ptr + stencil[Param::top_left]] += art_vect[i];
                            col = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);
                            row_indices[ptr + stencil[Param::top_right]] = row;
                            col_indices[ptr + stencil[Param::top_right]] = col;

                            vals[ptr + stencil[Param::top_right]] += -art_vect[i];

                            // Update (j, i+1)
                            ptr     = ptr_vect[i + 2];
                            stencil = stencil_cur;
                            row     = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);
                            col     = (j - 1) * ntheta_int + i;
                            row_indices[ptr + stencil[Param::bottom_left]] = row;
                            col_indices[ptr + stencil[Param::bottom_left]] = col;

                            vals[ptr + stencil[Param::bottom_left]] += -art_vect[i];
                            col                                         = (j + 1) * ntheta_int + i;
                            row_indices[ptr + stencil[Param::top_left]] = row;
                            col_indices[ptr + stencil[Param::top_left]] = col;

                            vals[ptr + stencil[Param::top_left]] += art_vect[i];

                            // Update (j, i-1)
                            ptr = ptr_vect[i];
                            row = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);
                            col = (j - 1) * ntheta_int + i;
                            row_indices[ptr + stencil[Param::bottom_right]] = row;
                            col_indices[ptr + stencil[Param::bottom_right]] = col;

                            vals[ptr + stencil[Param::bottom_right]] += art_vect[i];
                            col                                          = (j + 1) * ntheta_int + i;
                            row_indices[ptr + stencil[Param::top_right]] = row;
                            col_indices[ptr + stencil[Param::top_right]] = col;

                            vals[ptr + stencil[Param::top_right]] += -art_vect[i];

                            // Update (j+1, i) (Not in DB_ext)
                            ptr     = ptr_vect_next[i + 1];
                            stencil = stencil_next;
                            row     = (j + 1) * ntheta_int + i;
                            col     = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);
                            row_indices[ptr + stencil[Param::bottom_left]] = row;
                            col_indices[ptr + stencil[Param::bottom_left]] = col;

                            vals[ptr + stencil[Param::bottom_left]] += -art_vect[i];
                            col = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);
                            row_indices[ptr + stencil[Param::bottom_right]] = row;
                            col_indices[ptr + stencil[Param::bottom_right]] = col;

                            vals[ptr + stencil[Param::bottom_right]] += art_vect[i];
                        }
                } // end of task
            }

#pragma omp task shared(dep, start_j) depend(in                                                                        \
                                             : dep[nr_int - 2]) depend(in                                              \
                                                                       : dep[nr_int - 3]) depend(out                   \
                                                                                                 : dep[nr_int - 1])
            {
                /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!        Last lines       !!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     */
                // DB_ext updates (~~~ Interior - (j+1, i) + DB)
                int j, ptr, row, col;
                double coeff, coeff2, val, coeff3, kt, ktmin1, hs, hsmin1;
                std::vector<int> ptr_vect_prev, ptr_vect, ptr_vect_next;
                std::vector<int> stencil_prev, stencil_cur, stencil_next, stencil;
                std::vector<double> arr_vect, arr_vect2, att_vect, art_vect;

                j             = nr_int - 1;
                ptr_vect_prev = get_ptr(j - 1);
                ptr_vect      = get_ptr(j);
                stencil_prev  = get_stencil(j - 1);
                stencil_cur   = get_stencil(j);

                stencil = stencil_cur;
                for (int i = 0; i < ntheta_int; i++) {
                    ptr = ptr_vect[i + 1];
                    row = j * ntheta_int + i;
                    val = betaVec[row];

                    vals[ptr + stencil[Param::middle]] += val;
                }

                hs     = hplus[j];
                hsmin1 = hplus[j - 1];
                gyro::arr_att_art(r[j], theta, sin_theta, cos_theta, ntheta_int, arr_vect, att_vect, art_vect, 0);
                arr_vect2 = gyro::arr(r[j + 1], theta, sin_theta, cos_theta, ntheta_int, 0);
                for (int i = 0; i < ntheta_int; i++) {
                    // Define the index, position and interval size of current and previous node
                    // - in theta
                    kt     = thetaplus_per[i + 1];
                    ktmin1 = thetaplus_per[i];

                    // Update (j-1, i) (=== Interior)
                    ptr                                    = ptr_vect_prev[i + 1];
                    stencil                                = stencil_prev;
                    row                                    = (j - 1) * ntheta_int + i;
                    col                                    = j * ntheta_int + i;
                    row_indices[ptr + stencil[Param::top]] = row;
                    col_indices[ptr + stencil[Param::top]] = col;

                    coeff = 0.5 * (kt + ktmin1) * arr_vect[i] / hsmin1;

                    vals[ptr + stencil[Param::top]] += -coeff;

                    vals[ptr + stencil[Param::middle]] += coeff;

                    // Update (j, i) (=== Interior - top)
                    ptr                                       = ptr_vect[i + 1];
                    stencil                                   = stencil_cur;
                    row                                       = j * ntheta_int + i;
                    col                                       = row;
                    row_indices[ptr + stencil[Param::middle]] = row;
                    col_indices[ptr + stencil[Param::middle]] = col;

                    coeff = 0.5 * (kt + ktmin1) * arr_vect[i];
                    col   = (j - 1) * ntheta_int + i;

                    vals[ptr + stencil[Param::bottom]] += -coeff / hsmin1;

                    // Contribution to middle (top) from DB
                    coeff2 = 0.5 * (kt + ktmin1) * arr_vect2[i] / hs;

                    coeff3 = 0.5 * (hs + hsmin1) * att_vect[i];
                    col    = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);

                    vals[ptr + stencil[Param::left]] += -coeff3 / ktmin1;
                    col = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);

                    vals[ptr + stencil[Param::right]] += -coeff3 / kt;

                    vals[ptr + stencil[Param::middle]] +=
                        coeff / hsmin1 + coeff / hs + coeff2 + coeff3 / ktmin1 + coeff3 / kt;

                    // Update (j, i+1) (=== Interior - top_left)
                    ptr                                     = ptr_vect[i + 2];
                    row                                     = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);
                    col                                     = j * ntheta_int + i;
                    row_indices[ptr + stencil[Param::left]] = row;
                    col_indices[ptr + stencil[Param::left]] = col;

                    coeff = 0.5 * (hs + hsmin1) * att_vect[i] / kt;

                    vals[ptr + stencil[Param::left]] += -coeff;

                    vals[ptr + stencil[Param::middle]] += coeff;

                    // Update (j, i-1) (=== Interior - top_right)
                    ptr                                      = ptr_vect[i];
                    row                                      = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);
                    col                                      = j * ntheta_int + i;
                    row_indices[ptr + stencil[Param::right]] = row;
                    col_indices[ptr + stencil[Param::right]] = col;

                    coeff = 0.5 * (hs + hsmin1) * att_vect[i] / ktmin1;

                    vals[ptr + stencil[Param::right]] += -coeff;

                    vals[ptr + stencil[Param::middle]] += coeff;
                }
                if (gyro::icntl[Param::mod_pk] > 0)
                    for (int i = 0; i < ntheta_int; i++) {
                        // Define the index, position and interval size of current and previous node
                        // - in theta
                        kt     = thetaplus_per[i + 1];
                        ktmin1 = thetaplus_per[i];

                        // Update (j-1, i) (=== Interior)
                        ptr     = ptr_vect_prev[i + 1];
                        stencil = stencil_prev;
                        row     = (j - 1) * ntheta_int + i;
                        col     = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);
                        row_indices[ptr + stencil[Param::top_left]] = row;
                        col_indices[ptr + stencil[Param::top_left]] = col;

                        vals[ptr + stencil[Param::top_left]] += art_vect[i];
                        col = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);
                        row_indices[ptr + stencil[Param::top_right]] = row;
                        col_indices[ptr + stencil[Param::top_right]] = col;

                        vals[ptr + stencil[Param::top_right]] += -art_vect[i];

                        // Update (j, i+1) (=== Interior - top_left)
                        ptr = ptr_vect[i + 2];
                        row = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);
                        col = (j - 1) * ntheta_int + i;
                        row_indices[ptr + stencil[Param::bottom_left]] = row;
                        col_indices[ptr + stencil[Param::bottom_left]] = col;

                        vals[ptr + stencil[Param::bottom_left]] += -art_vect[i];

                        // Update (j, i-1) (=== Interior - top_right)
                        ptr = ptr_vect[i];
                        row = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);
                        col = (j - 1) * ntheta_int + i;
                        row_indices[ptr + stencil[Param::bottom_right]] = row;
                        col_indices[ptr + stencil[Param::bottom_right]] = col;

                        vals[ptr + stencil[Param::bottom_right]] += art_vect[i];
                    }
            } //end of task
        } //end of single
    } //end of parallel

    delete[] dep;
} /* ----- end of level::build_A ----- */

/*!
 *  \brief Build the RHS
 *
 * Builds the RHS of the system based on 9p FD. Relations with Dirichlet boundary 
 * condition nodes are shifted to the RHS to keep a symmetric operator.
 *
 */
void level::build_rhs()
{
    // double tol_bound_check = 1e-8;
    // double Rmax            = gyro::dcntl[Param::R];
    int start_j;
    int* dep = new int[nr];

    /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!       DB first line      !!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     */
    // Take boundary condition into account: Dirichlet-RB
    if (gyro::icntl[Param::DirBC_Interior]) { // (r[0],0) is on Dirichlet boundary
        start_j = 1;
    }
    else // (r[0],0) is not on Dirichlet boundary
        start_j = 0;

#pragma omp parallel shared(dep)
    {
#pragma omp single
        {
            if (gyro::icntl[Param::DirBC_Interior]) { // (r[0],0) is on Dirichlet boundary
#pragma omp task shared(dep, start_j) depend(out : dep[0])
                {
                    std::vector<double> sol = gyro::def_solution_rt(r[0], theta, sin_theta, cos_theta, ntheta_int, 0);
                    for (int i = 0; i < ntheta_int; i++) {

                        fVec[i] = sol[i];
                    }
                } //end of task
            }

            /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!       Fill base RHS      !!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     */
            // across the origin (treat separately because of hsmin1)
            // std::vector<double> coeff_b = gyro::coeff_beta(r, 0);
            for (int j = start_j; j < nr_int; j++) {
#pragma omp task shared(dep, start_j) firstprivate(j) depend(out : dep[j])
                {
                    int i, row;
                    double kt, ktmin1, hs, hsmin1;
                    if (j == 0)
                        hsmin1 = 2 * r[0];
                    else
                        hsmin1 = hplus[j - 1];
                    std::vector<double> det = gyro::detDFinv(r[j], theta, sin_theta, cos_theta, ntheta_int, 0);
                    std::vector<double> rhs = gyro::eval_def_rhs(r[j], theta, sin_theta, cos_theta, ntheta_int, 0);
                    // Define the index, position and interval size of current and previous node
                    // - in r
                    hs = hplus[j];
                    for (i = 0; i < ntheta_int; i++) {
                        row = j * ntheta_int + i;
                        // - in theta
                        kt     = thetaplus_per[i + 1];
                        ktmin1 = thetaplus_per[i];

                        double fac_hk = 0.25 * (kt + ktmin1) * (hs + hsmin1);

                        fVec[row] = fac_hk * fabs(det[i]) * rhs[i];
                    }
                    hsmin1 = hs;
                }
            } //end of task

/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!          DB_int          !!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     */
#pragma omp task shared(dep, start_j) depend(in : dep[1])
            {
                // int hs, ;
                int j = 1;
                if (gyro::icntl[Param::DirBC_Interior]) {
                    int i, row;
                    double kt, ktmin1;
                    std::vector<double> sol =
                        gyro::def_solution_rt(r[j - 1], theta, sin_theta, cos_theta, ntheta_int, 0);
                    if (gyro::icntl[Param::mod_pk] > 0) {
                        std::vector<double> art_vect_prev =
                            gyro::art(r[j - 1], theta_per, sin_theta_per, cos_theta_per, ntheta_int + 1, 0);
                        std::vector<double> art_vect =
                            gyro::art(r[j], theta_per, sin_theta_per, cos_theta_per, ntheta_int + 1, 0);
                        for (i = 1; i < ntheta_int - 1; i++) {
                            row = j * ntheta_int + i;
                            // Define the index, position and interval size of current and previous node
                            // - in theta
                            kt     = thetaplus_per[i + 1];
                            ktmin1 = thetaplus_per[i];

                            // bottom_left
                            fVec[row] += (art_vect[i] + art_vect_prev[i + 1]) * sol[i - 1];
                            // bottom_right
                            fVec[row] -= (art_vect_prev[i + 1] + art_vect[i + 2]) * sol[i + 1];
                        }
                        // i=0
                        i   = 0;
                        row = j * ntheta_int + i;

                        fVec[row] += (art_vect[ntheta_int] + art_vect_prev[1]) * sol[ntheta_int - 1];

                        fVec[row] -= (art_vect_prev[1] + art_vect[2]) * sol[1];

                        // i=ntheta_int-1
                        i   = ntheta_int - 1;
                        row = j * ntheta_int + i;

                        fVec[row] += (art_vect[ntheta_int - 1] + art_vect_prev[ntheta_int]) * sol[ntheta_int - 2];

                        fVec[row] -= (art_vect_prev[ntheta_int] + art_vect[1]) * sol[0];
                    }
                    std::vector<double> arr_vect_prev = gyro::arr(r[j - 1], theta, sin_theta, cos_theta, ntheta_int, 0);
                    std::vector<double> arr_vect      = gyro::arr(r[j], theta, sin_theta, cos_theta, ntheta_int, 0);
                    for (i = 0; i < ntheta_int; i++) {
                        row = j * ntheta_int + i;
                        // Define the index, position and interval size of current and previous node
                        // - in theta
                        kt     = thetaplus_per[i + 1];
                        ktmin1 = thetaplus_per[i];
                        // - in r
                        double hsmin1 = hplus[j - 1];

                        // bottom

                        fVec[row] += 0.5 * (kt + ktmin1) * (arr_vect_prev[i] + arr_vect[i]) * sol[i] / hsmin1;
                    }
                }
            } //end of task

/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!          DB_ext          !!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     */
#pragma omp task shared(dep, start_j) depend(in : dep[nr_int - 1])
            {
                int i, j, row;
                double kt, ktmin1, hs;
                // double hsmin1;
                j                       = nr_int - 1;
                std::vector<double> sol = gyro::def_solution_rt(r[j + 1], theta, sin_theta, cos_theta, ntheta_int, 0);
                if (gyro::icntl[Param::mod_pk] > 0) {
                    std::vector<double> art_vect =
                        gyro::art(r[j], theta_per, sin_theta_per, cos_theta_per, ntheta_int + 1, 0);
                    std::vector<double> art_vect_next =
                        gyro::art(r[j + 1], theta_per, sin_theta_per, cos_theta_per, ntheta_int + 1, 0);
                    for (i = 1; i < ntheta_int - 1; i++) {
                        row = j * ntheta_int + i;
                        // Define the index, position and interval size of current and previous node
                        // - in theta
                        kt     = thetaplus_per[i + 1];
                        ktmin1 = thetaplus_per[i];
                        // - in r
                        hs = hplus[j];

                        // top_left

                        fVec[row] -= (art_vect[i] + art_vect_next[i + 1]) * sol[i - 1];
                        // top_right

                        fVec[row] += (art_vect_next[i + 1] + art_vect[i + 2]) * sol[i + 1];
                    }
                    // i=0
                    i   = 0;
                    row = j * ntheta_int + i;

                    fVec[row] -= (art_vect[ntheta_int] + art_vect_next[1]) * sol[ntheta_int - 1];

                    fVec[row] += (art_vect_next[1] + art_vect[2]) * sol[1];

                    // i=ntheta_int-1
                    i   = ntheta_int - 1;
                    row = j * ntheta_int + i;

                    fVec[row] -= (art_vect[ntheta_int - 1] + art_vect_next[ntheta_int]) * sol[ntheta_int - 2];

                    fVec[row] += (art_vect_next[ntheta_int] + art_vect[1]) * sol[0];
                }
                std::vector<double> arr_vect      = gyro::arr(r[j], theta, sin_theta, cos_theta, ntheta_int, 0);
                std::vector<double> arr_vect_next = gyro::arr(r[j + 1], theta, sin_theta, cos_theta, ntheta_int, 0);
                for (i = 0; i < ntheta_int; i++) {
                    row = j * ntheta_int + i;
                    // Define the index, position and interval size of current and previous node
                    // - in theta
                    kt     = thetaplus_per[i + 1];
                    ktmin1 = thetaplus_per[i];
                    // - in r
                    hs = hplus[j];

                    // top

                    fVec[row] += (0.5 * (kt + ktmin1) * (arr_vect[i] + arr_vect_next[i])) * sol[i] / hs;
                }
            } //end of task

/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!       DB last line      !!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     */
#pragma omp task shared(dep, start_j) depend(out : dep[nr_int])
            {
                std::vector<double> sol = gyro::def_solution_rt(r[nr_int], theta, sin_theta, cos_theta, ntheta_int, 0);
                for (int i = 0; i < ntheta_int; i++) {

                    fVec[m - ntheta_int + i] = sol[i];
                }
            }
        } //end of single
    } //end of parallel

    delete[] dep;
} /* ----- end of level::build_rhs ----- */

/*!
 *  \brief Build the RHS
 *
 * Builds the RHS of the system based on 9p FD. Relations with Dirichlet boundary 
 * condition nodes are shifted to the RHS to keep a symmetric operator.
 *
 */
void level::build_betaVec()
{
    // #pragma omp parallel shared(dep)
    {
        // #pragma omp single
        {
            // across the origin (treat separately because of hsmin1)
            for (int j = 0; j < nr_int; j++) {
                // #pragma omp task shared(dep, start_j) firstprivate(j) depend(out : dep[j])
                {
                    int i, row;
                    double kt, ktmin1, hs, hsmin1;
                    if (j == 0)
                        hsmin1 = 2 * r[0];
                    else
                        hsmin1 = hplus[j - 1];
                    std::vector<double> det = gyro::detDFinv(r[j], theta, sin_theta, cos_theta, ntheta_int, 0);
                    // Define the index, position and interval size of current and previous node
                    // - in r
                    hs = hplus[j];
                    for (i = 0; i < ntheta_int; i++) {
                        row = j * ntheta_int + i;
                        // - in theta
                        kt             = thetaplus_per[i + 1];
                        ktmin1         = thetaplus_per[i];
                        double coeff_b = gyro::coeff_beta(r[j], 0);

                        double fac_hk = 0.25 * (kt + ktmin1) * (hs + hsmin1);
                        betaVec[row]  = fac_hk * coeff_b * fabs(det[i]);
                    }
                    hsmin1 = hs;
                }
            } //end of task
        } //end of single
    } //end of parallel
} /* ----- end of level::build_betaVec ----- */

/*!
 *  \brief Applies the operator A without construction
 *
 * Applies the matrix A of the system based on 9p FD. Relations with Dirichlet 
 * boundary condition nodes are shifted to the RHS to keep a symmetric operator.
 *
 * everything expressed in polar coordinates
 * allows anisotropy in r
 * - r: array of discretized r values
 * - ntheta: number of discretization intervals in angle phi (also variable m)
 *
 * solves the transformed DGL:  r*u_rr + u_r + 1/r*u_phiphi = rf fuer r>0, r<=1,
 * fuer r=0, approximate: u_0 = 1/m sum_j u(h,jk)-(h/2)^2*f
 * sorting of nodes ringwise. first node origin, last ntheta
 * Node all the way out
 *
 * \param u: vector on which to apply A
 * \param Au: result vector = A*u
 *
 */
void level::apply_A(std::vector<double> u, std::vector<double>& Au)
{
    int start_j;
    int* dep = new int[nr];

    /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!       DB first line      !!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     */
    // Take boundary condition into account: Dirichlet-RB
    if (gyro::icntl[Param::DirBC_Interior]) { // (r[0],0) is on Dirichlet boundary
        start_j = 1;
    }
    else { // (r[0],0) is not on Dirichlet boundary
        start_j = 0;
    }

#pragma omp parallel shared(dep, Au)
    {
#pragma omp single
        {
            if (gyro::icntl[Param::DirBC_Interior]) { // (r[0],0) is on Dirichlet boundary
#pragma omp task depend(out : dep[0])
                {
                    for (int i = 0; i < ntheta_int; i++) {

                        Au[i] += u[i];
                    }
                } //end of task
            }

#pragma omp task depend(out : dep[nr_int])
            {
                // Take boundary condition into account: Dirichlet-RB
                for (int i = 0; i < ntheta_int; i++) {
                    int row = m - ntheta_int + i;

                    Au[row] += u[row];
                }
            } //end of task and parallel

#pragma omp task depend(out : dep[start_j])
            {
                double coeff, coeff2, val, kt, ktmin1, hs, hsmin1;
                int row, col;
                std::vector<double> arr_vect, arr_vect2, att_vect, art_vect;

                /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!       First lines       !!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     */
                int j = start_j;

                for (int i = 0; i < ntheta_int; i++) {
                    row = j * ntheta_int + i;
                    val = betaVec[row];

                    Au[row] += val * u[row];
                }

                hs       = hplus[j];
                arr_vect = gyro::arr(r[j], theta, sin_theta, cos_theta, ntheta_int, 0);
                // Across: bottom update
                if (!gyro::icntl[Param::DirBC_Interior]) {
                    hsmin1 = 2 * r[0];
                    // Accross the origin theta and arr
                    arr_vect2 = gyro::arr(r[j], theta_PI, sin_theta_PI, cos_theta_PI, ntheta_int, 0);

                    for (int i = 0; i < ntheta_int; i++) {
                        // Define the index, position and interval size of current and previous node
                        // - in theta
                        kt     = thetaplus_per[i + 1];
                        ktmin1 = thetaplus_per[i];

                        row = j * ntheta_int + i;
                        if (i + 1 <= ntheta_int / 2) // first half of circle: half a turn further
                            col = row + ntheta_int / 2;
                        else // second half of circle: half a turn back
                            col = row - ntheta_int / 2;

                        coeff  = 0.5 * (kt + ktmin1) * arr_vect2[i] / hsmin1;
                        coeff2 = 0.5 * (kt + ktmin1) * arr_vect[i] / hsmin1;
                        // bottom
                        val = -coeff - coeff2;

                        Au[row] += val * u[col];
                        // middle
                        val = coeff;

                        Au[row] += val * u[row];
                    }
                }
                else {
                    hsmin1 = hplus[j - 1];
                    // DB contribution arr (r(0))
                    arr_vect2 = gyro::arr(r[j - 1], theta, sin_theta, cos_theta, ntheta_int, 0);

                    for (int i = 0; i < ntheta_int; i++) {
                        // Define the index, position and interval size of current and previous node
                        // - in theta
                        kt     = thetaplus_per[i + 1];
                        ktmin1 = thetaplus_per[i];

                        coeff = 0.5 * (kt + ktmin1) * arr_vect2[i] / hsmin1;
                        row   = j * ntheta_int + i;
                        // middle
                        val = coeff;

                        Au[row] += val * u[row];
                    }
                }

                // Across and DB_int updates (~~~ Interior - (j-1, i))
                gyro::arr_att_art(r[j], theta, sin_theta, cos_theta, ntheta_int, arr_vect, att_vect, art_vect, 0);
                for (int i = 0; i < ntheta_int; i++) {
                    // Define the index, position and interval size of current and previous node
                    // - in theta
                    kt     = thetaplus_per[i + 1];
                    ktmin1 = thetaplus_per[i];

                    // Update (j, i) (=== Interior - bottom)
                    row = j * ntheta_int + i;

                    coeff = 0.5 * (kt + ktmin1) * arr_vect[i];
                    col   = (j + 1) * ntheta_int + i;
                    // top
                    val = -coeff / hs;

                    Au[row] += val * u[col];

                    coeff2 = 0.5 * (hs + hsmin1) * att_vect[i];
                    col    = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);
                    // left
                    val = -coeff2 / ktmin1;

                    Au[row] += val * u[col];
                    col = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);
                    // right
                    val = -coeff2 / kt;

                    Au[row] += val * u[col];

                    // middle
                    val = coeff / hs + coeff / hsmin1 + coeff2 / ktmin1 + coeff2 / kt;

                    Au[row] += val * u[row];

                    // Update (j, i+1) (=== Interior - bottom_left)
                    row = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);
                    col = j * ntheta_int + i;

                    coeff = 0.5 * (hs + hsmin1) * att_vect[i] / kt;
                    // left
                    val = -coeff;

                    Au[row] += val * u[col];
                    // middle
                    val = coeff;

                    Au[row] += val * u[row];

                    // Update (j, i-1) (=== Interior - bottom_right)
                    row = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);
                    col = j * ntheta_int + i;

                    coeff = 0.5 * (hs + hsmin1) * att_vect[i] / ktmin1;
                    // right
                    val = -coeff;

                    Au[row] += val * u[col];
                    // middle
                    val = coeff;

                    Au[row] += val * u[row];

                    // Update (j+1, i) (=== Interior)
                    row = (j + 1) * ntheta_int + i;
                    col = j * ntheta_int + i;

                    coeff = 0.5 * (kt + ktmin1) * arr_vect[i] / hs;
                    // bottom
                    val = -coeff;

                    Au[row] += val * u[col];
                    // middle
                    val = coeff;

                    Au[row] += val * u[row];
                }

                if (gyro::icntl[Param::mod_pk] > 0) {
                    for (int i = 0; i < ntheta_int; i++) {
                        // Define the index, position and interval size of current and previous node
                        // - in theta
                        kt     = thetaplus_per[i + 1];
                        ktmin1 = thetaplus_per[i];

                        // Update (j, i+1) (=== Interior - bottom_left)
                        row = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);
                        col = (j + 1) * ntheta_int + i;
                        // top_left
                        val = art_vect[i];

                        Au[row] += val * u[col];

                        // Update (j, i-1) (=== Interior - bottom_right)
                        row = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);
                        col = (j + 1) * ntheta_int + i;
                        // top_right
                        val = -art_vect[i];

                        Au[row] += val * u[col];

                        // Update (j+1, i) (=== Interior)
                        row = (j + 1) * ntheta_int + i;
                        col = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);
                        // bottom_left
                        val = -art_vect[i];

                        Au[row] += val * u[col];
                        col = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);
                        // bottom_right
                        val = art_vect[i];

                        Au[row] += val * u[col];
                    }
                }
            } // end of task --------------------------------------------------------------------------------------------------

            /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!      Interior nodes  (1)    !!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     */
            for (int j = start_j + 3; j < nr_int - 1; j += 3) {
#pragma omp task firstprivate(j) depend(out : dep[j])
                {
                    double coeff, coeff2, val, kt, ktmin1, hs, hsmin1;
                    int row, col;
                    std::vector<double> arr_vect, arr_vect2, att_vect, art_vect;

                    for (int i = 0; i < ntheta_int; i++) {
                        row = j * ntheta_int + i;
                        val = betaVec[row];

                        Au[row] += val * u[row];
                    }

                    // - in r
                    hs     = hplus[j];
                    hsmin1 = hplus[j - 1];
                    gyro::arr_att_art(r[j], theta, sin_theta, cos_theta, ntheta_int, arr_vect, att_vect, art_vect, 0);

                    for (int i = 0; i < ntheta_int; i++) {
                        // Define the index, position and interval size of current and previous node
                        // - in theta
                        kt     = thetaplus_per[i + 1];
                        ktmin1 = thetaplus_per[i];

                        // Update (j-1, i) (Not in DB_int)
                        row = (j - 1) * ntheta_int + i;
                        col = j * ntheta_int + i;

                        coeff = 0.5 * (kt + ktmin1) * arr_vect[i] / hsmin1;
                        // top
                        val = -coeff;

                        Au[row] += val * u[col];
                        // middle
                        val = coeff;

                        Au[row] += val * u[row];

                        // Update (j, i)
                        row = j * ntheta_int + i;

                        coeff = 0.5 * (kt + ktmin1) * arr_vect[i];
                        col   = (j + 1) * ntheta_int + i;
                        // top
                        val = -coeff / hs;

                        Au[row] += val * u[col];
                        col = (j - 1) * ntheta_int + i;
                        // bottom
                        val = -coeff / hsmin1;

                        Au[row] += val * u[col];

                        coeff2 = 0.5 * (hs + hsmin1) * att_vect[i];
                        col    = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);
                        // left
                        val = -coeff2 / ktmin1;
                        Au[row] += val * u[col];

                        col = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);
                        // right
                        val = -coeff2 / kt;

                        Au[row] += val * u[col];

                        // middle
                        val = coeff / hs + coeff / hsmin1 + coeff2 / ktmin1 + coeff2 / kt;

                        Au[row] += val * u[row];

                        // Update (j, i+1)
                        row = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);
                        col = j * ntheta_int + i;

                        coeff = 0.5 * (hs + hsmin1) * att_vect[i] / kt;
                        // left
                        val = -coeff;

                        Au[row] += val * u[col];
                        // middle
                        val = coeff;

                        Au[row] += val * u[row];

                        // Update (j, i-1)
                        row = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);
                        col = j * ntheta_int + i;

                        coeff = 0.5 * (hs + hsmin1) * att_vect[i] / ktmin1;
                        // right
                        val = -coeff;

                        Au[row] += val * u[col];
                        // middle
                        val = coeff;

                        Au[row] += val * u[row];

                        // Update (j+1, i) (Not in DB_ext)
                        row = (j + 1) * ntheta_int + i;
                        col = j * ntheta_int + i;

                        coeff = 0.5 * (kt + ktmin1) * arr_vect[i] / hs;
                        // bottom
                        val = -coeff;

                        Au[row] += val * u[col];
                        // middle
                        val = coeff;

                        Au[row] += val * u[row];
                    }
                    if (gyro::icntl[Param::mod_pk] > 0) {
                        for (int i = 0; i < ntheta_int; i++) {
                            // Define the index, position and interval size of current and previous node
                            // - in theta
                            kt     = thetaplus_per[i + 1];
                            ktmin1 = thetaplus_per[i];

                            // Update (j-1, i) (Not in DB_int)
                            row = (j - 1) * ntheta_int + i;
                            col = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);
                            // top_left
                            val = art_vect[i];

                            Au[row] += val * u[col];
                            col = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);
                            // top_right
                            val = -art_vect[i];

                            Au[row] += val * u[col];

                            // Update (j, i+1)
                            row = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);
                            col = (j - 1) * ntheta_int + i;
                            // bottom_left
                            val = -art_vect[i];

                            Au[row] += val * u[col];
                            col = (j + 1) * ntheta_int + i;
                            // top_left
                            val = art_vect[i];

                            Au[row] += val * u[col];

                            // Update (j, i-1)
                            row = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);
                            col = (j - 1) * ntheta_int + i;
                            // bottom_right
                            val = art_vect[i];

                            Au[row] += val * u[col];
                            col = (j + 1) * ntheta_int + i;
                            // top_right
                            val = -art_vect[i];

                            Au[row] += val * u[col];

                            // Update (j+1, i) (Not in DB_ext)
                            row = (j + 1) * ntheta_int + i;
                            col = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);
                            // bottom_left
                            val = -art_vect[i];

                            Au[row] += val * u[col];
                            col = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);
                            // bottom_right
                            val = art_vect[i];

                            Au[row] += val * u[col];
                        }
                    }
                } // end of task -----------------------------------------------------------------------------------------------
            } //end of for

            /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!      Interior nodes (2)     !!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     */
            for (int j = start_j + 1; j < nr_int - 1; j += 3) {
#pragma omp task firstprivate(j) depend(in : dep[j - 1]) depend(in : dep[j + 2]) depend(out : dep[j])
                {
                    double coeff, coeff2, val, kt, ktmin1, hs, hsmin1;
                    int row, col;
                    std::vector<double> arr_vect, arr_vect2, att_vect, art_vect;

                    for (int i = 0; i < ntheta_int; i++) {
                        row = j * ntheta_int + i;
                        val = betaVec[row];

                        Au[row] += val * u[row];
                    }

                    // - in r
                    hs     = hplus[j];
                    hsmin1 = hplus[j - 1];
                    gyro::arr_att_art(r[j], theta, sin_theta, cos_theta, ntheta_int, arr_vect, att_vect, art_vect, 0);

                    for (int i = 0; i < ntheta_int; i++) {
                        // Define the index, position and interval size of current and previous node
                        // - in theta
                        kt     = thetaplus_per[i + 1];
                        ktmin1 = thetaplus_per[i];

                        // Update (j-1, i) (Not in DB_int)
                        row = (j - 1) * ntheta_int + i;
                        col = j * ntheta_int + i;

                        coeff = 0.5 * (kt + ktmin1) * arr_vect[i] / hsmin1;
                        // top
                        val = -coeff;

                        Au[row] += val * u[col];
                        // middle
                        val = coeff;

                        Au[row] += val * u[row];

                        // Update (j, i)
                        row = j * ntheta_int + i;

                        coeff = 0.5 * (kt + ktmin1) * arr_vect[i];
                        col   = (j + 1) * ntheta_int + i;
                        // top
                        val = -coeff / hs;

                        Au[row] += val * u[col];
                        col = (j - 1) * ntheta_int + i;
                        // bottom
                        val = -coeff / hsmin1;

                        Au[row] += val * u[col];

                        coeff2 = 0.5 * (hs + hsmin1) * att_vect[i];
                        col    = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);
                        // left
                        val = -coeff2 / ktmin1;

                        Au[row] += val * u[col];
                        col = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);
                        // right
                        val = -coeff2 / kt;

                        Au[row] += val * u[col];

                        // middle
                        val = coeff / hs + coeff / hsmin1 + coeff2 / ktmin1 + coeff2 / kt;

                        Au[row] += val * u[row];

                        // Update (j, i+1)
                        row = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);
                        col = j * ntheta_int + i;

                        coeff = 0.5 * (hs + hsmin1) * att_vect[i] / kt;
                        // left
                        val = -coeff;

                        Au[row] += val * u[col];
                        // middle
                        val = coeff;

                        Au[row] += val * u[row];

                        // Update (j, i-1)
                        row = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);
                        col = j * ntheta_int + i;

                        coeff = 0.5 * (hs + hsmin1) * att_vect[i] / ktmin1;
                        // right
                        val = -coeff;

                        Au[row] += val * u[col];
                        // middle
                        val = coeff;

                        Au[row] += val * u[row];

                        // Update (j+1, i) (Not in DB_ext)
                        row = (j + 1) * ntheta_int + i;
                        col = j * ntheta_int + i;

                        coeff = 0.5 * (kt + ktmin1) * arr_vect[i] / hs;
                        // bottom
                        val = -coeff;

                        Au[row] += val * u[col];
                        // middle
                        val = coeff;

                        Au[row] += val * u[row];
                    }
                    if (gyro::icntl[Param::mod_pk] > 0) {
                        for (int i = 0; i < ntheta_int; i++) {
                            // Define the index, position and interval size of current and previous node
                            // - in theta
                            kt     = thetaplus_per[i + 1];
                            ktmin1 = thetaplus_per[i];

                            // Update (j-1, i) (Not in DB_int)
                            row = (j - 1) * ntheta_int + i;
                            col = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);
                            // top_left
                            val = art_vect[i];

                            Au[row] += val * u[col];
                            col = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);
                            // top_right
                            val = -art_vect[i];

                            Au[row] += val * u[col];

                            // Update (j, i+1)
                            row = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);
                            col = (j - 1) * ntheta_int + i;
                            // bottom_left
                            val = -art_vect[i];

                            Au[row] += val * u[col];
                            col = (j + 1) * ntheta_int + i;
                            // top_left
                            val = art_vect[i];

                            Au[row] += val * u[col];

                            // Update (j, i-1)
                            row = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);
                            col = (j - 1) * ntheta_int + i;
                            // bottom_right
                            val = art_vect[i];

                            Au[row] += val * u[col];
                            col = (j + 1) * ntheta_int + i;
                            // top_right
                            val = -art_vect[i];

                            Au[row] += val * u[col];

                            // Update (j+1, i) (Not in DB_ext)
                            row = (j + 1) * ntheta_int + i;
                            col = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);
                            // bottom_left
                            val = -art_vect[i];

                            Au[row] += val * u[col];
                            col = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);
                            // bottom_right
                            val = art_vect[i];

                            Au[row] += val * u[col];
                        }
                    }
                } // end of task -------------------------------------------------------------------------------------------
            } //end of for

            /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!      Interior nodes (2)     !!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     */
            for (int j = start_j + 2; j < nr_int - 1; j += 3) {
#pragma omp task depend(in : dep[j - 1]) depend(in : dep[j + 2]) depend(out : dep[j])
                {
                    double coeff, coeff2, val, kt, ktmin1, hs, hsmin1;
                    int row, col;
                    std::vector<double> arr_vect, arr_vect2, att_vect, art_vect;

                    for (int i = 0; i < ntheta_int; i++) {
                        row = j * ntheta_int + i;
                        val = betaVec[row];

                        Au[row] += val * u[row];
                    }

                    // - in r
                    hs     = hplus[j];
                    hsmin1 = hplus[j - 1];
                    gyro::arr_att_art(r[j], theta, sin_theta, cos_theta, ntheta_int, arr_vect, att_vect, art_vect, 0);

                    for (int i = 0; i < ntheta_int; i++) {
                        // Define the index, position and interval size of current and previous node
                        // - in theta
                        kt     = thetaplus_per[i + 1];
                        ktmin1 = thetaplus_per[i];

                        // Update (j-1, i) (Not in DB_int)
                        row = (j - 1) * ntheta_int + i;
                        col = j * ntheta_int + i;

                        coeff = 0.5 * (kt + ktmin1) * arr_vect[i] / hsmin1;
                        // top
                        val = -coeff;

                        Au[row] += val * u[col];
                        // middle
                        val = coeff;

                        Au[row] += val * u[row];

                        // Update (j, i)
                        row = j * ntheta_int + i;

                        coeff = 0.5 * (kt + ktmin1) * arr_vect[i];
                        col   = (j + 1) * ntheta_int + i;
                        // top
                        val = -coeff / hs;

                        Au[row] += val * u[col];
                        col = (j - 1) * ntheta_int + i;
                        // bottom
                        val = -coeff / hsmin1;

                        Au[row] += val * u[col];

                        coeff2 = 0.5 * (hs + hsmin1) * att_vect[i];
                        col    = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);
                        // left
                        val = -coeff2 / ktmin1;

                        Au[row] += val * u[col];
                        col = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);
                        // right
                        val = -coeff2 / kt;

                        Au[row] += val * u[col];

                        // middle
                        val = coeff / hs + coeff / hsmin1 + coeff2 / ktmin1 + coeff2 / kt;

                        Au[row] += val * u[row];

                        // Update (j, i+1)
                        row = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);
                        col = j * ntheta_int + i;

                        coeff = 0.5 * (hs + hsmin1) * att_vect[i] / kt;
                        // left
                        val = -coeff;

                        Au[row] += val * u[col];
                        // middle
                        val = coeff;

                        Au[row] += val * u[row];

                        // Update (j, i-1)
                        row = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);
                        col = j * ntheta_int + i;

                        coeff = 0.5 * (hs + hsmin1) * att_vect[i] / ktmin1;
                        // right
                        val = -coeff;

                        Au[row] += val * u[col];
                        // middle
                        val = coeff;

                        Au[row] += val * u[row];

                        // Update (j+1, i) (Not in DB_ext)
                        row = (j + 1) * ntheta_int + i;
                        col = j * ntheta_int + i;

                        coeff = 0.5 * (kt + ktmin1) * arr_vect[i] / hs;
                        // bottom
                        val = -coeff;

                        Au[row] += val * u[col];
                        // middle
                        val = coeff;

                        Au[row] += val * u[row];
                    }
                    if (gyro::icntl[Param::mod_pk] > 0) {
                        for (int i = 0; i < ntheta_int; i++) {
                            // Define the index, position and interval size of current and previous node
                            // - in theta
                            kt     = thetaplus_per[i + 1];
                            ktmin1 = thetaplus_per[i];

                            // Update (j-1, i) (Not in DB_int)
                            row = (j - 1) * ntheta_int + i;
                            col = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);
                            // top_left
                            val = art_vect[i];

                            Au[row] += val * u[col];
                            col = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);
                            // top_right
                            val = -art_vect[i];

                            Au[row] += val * u[col];

                            // Update (j, i+1)
                            row = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);
                            col = (j - 1) * ntheta_int + i;
                            // bottom_left
                            val = -art_vect[i];

                            Au[row] += val * u[col];
                            col = (j + 1) * ntheta_int + i;
                            // top_left
                            val = art_vect[i];

                            Au[row] += val * u[col];

                            // Update (j, i-1)
                            row = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);
                            col = (j - 1) * ntheta_int + i;
                            // bottom_right
                            val = art_vect[i];

                            Au[row] += val * u[col];
                            col = (j + 1) * ntheta_int + i;
                            // top_right
                            val = -art_vect[i];

                            Au[row] += val * u[col];

                            // Update (j+1, i) (Not in DB_ext)
                            row = (j + 1) * ntheta_int + i;
                            col = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);
                            // bottom_left
                            val = -art_vect[i];

                            Au[row] += val * u[col];
                            col = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);
                            // bottom_right
                            val = art_vect[i];

                            Au[row] += val * u[col];
                        }
                    }
                } // end of task -------------------------------------------------------------------------------------------
            } //end of for

#pragma omp task depend(in : dep[nr_int - 2]) depend(in : dep[nr_int - 3]) depend(out : dep[nr_int - 1])
            {
                double coeff, coeff2, coeff3, val, kt, ktmin1, hs, hsmin1;
                int row, col;
                std::vector<double> arr_vect, arr_vect2, att_vect, art_vect;

                /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!        Last lines (1)      !!!!!!!!!!!!!!!!!!!!!!!!
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     */
                // DB_ext updates (~~~ Interior - (j+1, i) + DB)
                int j = nr_int - 1;
                for (int i = 0; i < ntheta_int; i++) {
                    row = j * ntheta_int + i;
                    val = betaVec[row];

                    Au[row] += val * u[row];
                }

                hs     = hplus[j];
                hsmin1 = hplus[j - 1];
                gyro::arr_att_art(r[j], theta, sin_theta, cos_theta, ntheta_int, arr_vect, att_vect, art_vect, 0);
                arr_vect2 = gyro::arr(r[j + 1], theta, sin_theta, cos_theta, ntheta_int, 0);
                for (int i = 0; i < ntheta_int; i++) {
                    // Define the index, position and interval size of current and previous node
                    // - in theta
                    kt     = thetaplus_per[i + 1];
                    ktmin1 = thetaplus_per[i];

                    // Update (j-1, i) (=== Interior)
                    row = (j - 1) * ntheta_int + i;
                    col = j * ntheta_int + i;

                    coeff = 0.5 * (kt + ktmin1) * arr_vect[i] / hsmin1;
                    // top
                    val = -coeff;

                    Au[row] += val * u[col];
                    // middle
                    val = coeff;

                    Au[row] += val * u[row];

                    // Update (j, i) (=== Interior - top)
                    row = j * ntheta_int + i;
                    col = row;

                    coeff = 0.5 * (kt + ktmin1) * arr_vect[i];
                    col   = (j - 1) * ntheta_int + i;
                    // bottom
                    val = -coeff / hsmin1;

                    Au[row] += val * u[col];

                    // Contribution to middle (top) from DB
                    coeff2 = 0.5 * (kt + ktmin1) * arr_vect2[i] / hs;

                    coeff3 = 0.5 * (hs + hsmin1) * att_vect[i];
                    col    = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);
                    // left
                    val = -coeff3 / ktmin1;

                    Au[row] += val * u[col];
                    col = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);
                    // right
                    val = -coeff3 / kt;

                    Au[row] += val * u[col];

                    // middle
                    val = coeff / hsmin1 + coeff / hs + coeff2 + coeff3 / ktmin1 + coeff3 / kt;

                    Au[row] += val * u[row];

                    // Update (j, i+1) (=== Interior - top_left)
                    row = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);
                    col = j * ntheta_int + i;

                    coeff = 0.5 * (hs + hsmin1) * att_vect[i] / kt;
                    // left
                    val = -coeff;

                    Au[row] += val * u[col];
                    // middle
                    val = coeff;

                    Au[row] += val * u[row];

                    // Update (j, i-1) (=== Interior - top_right)
                    row = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);
                    col = j * ntheta_int + i;

                    coeff = 0.5 * (hs + hsmin1) * att_vect[i] / ktmin1;
                    // right
                    val = -coeff;

                    Au[row] += val * u[col];
                    // middle
                    val = coeff;

                    Au[row] += val * u[row];
                }
                if (gyro::icntl[Param::mod_pk] > 0) {
                    for (int i = 0; i < ntheta_int; i++) {
                        // Define the index, position and interval size of current and previous node
                        // - in theta
                        kt     = thetaplus_per[i + 1];
                        ktmin1 = thetaplus_per[i];

                        // Update (j-1, i) (=== Interior)
                        row = (j - 1) * ntheta_int + i;
                        col = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);
                        // top_left
                        val = art_vect[i];

                        Au[row] += val * u[col];
                        col = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);
                        // top_right
                        val = -art_vect[i];

                        Au[row] += val * u[col];

                        // Update (j, i+1) (=== Interior - top_left)
                        row = j * ntheta_int + ((i + 1 > ntheta_int - 1) ? 0 : i + 1);
                        col = (j - 1) * ntheta_int + i;
                        // bottom_left
                        val = -art_vect[i];

                        Au[row] += val * u[col];

                        // Update (j, i-1) (=== Interior - top_right)
                        row = j * ntheta_int + ((i - 1 < 0) ? ntheta_int - 1 : i - 1);
                        col = (j - 1) * ntheta_int + i;
                        // bottom_right
                        val = art_vect[i];

                        Au[row] += val * u[col];
                    }
                }
            } //end of task
        } //end of single
    } //end of parallel

    delete[] dep;
} /* ----- end of level::apply_A ----- */
