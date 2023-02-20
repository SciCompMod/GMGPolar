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
 * \file build_bi_aniso_rdir_phiper.cpp
 * \brief Implementation of the prolongation operators (deprecated)
 * \author M. Kuehn, C. Kruse, P. Leleux, C. Schwarz
 * \version 0.0
 */
#include "level.h"

/*!
 *  \brief Applies the bilinear interpolation
 *
 * Applies bilinear interpolation for Dirichlet boundary conditions in
 * variable r and periodic boundary conditions in phi. Can be used for a
 * disk as well as for an annulus. (might be used for other geometries with
 * Dirichlet boundary conditions in 'x' and periodic boundary conditions
 * in 'y' as well; not tested!)
 *
 * uses Anisotropic Bilinear Interpolation stencil
 * for isotropic mesh, anisotropic stencil reduces to isotropic definition
 *      ]1  2   1[
 *  1/4 ]2  4   2[
 *      ]1  2   1[
 *
 * \param u: vector to prolongate/restrict
 * \param mc: coarse size of prolongator (level->mc)
 * \param ncoarse: fine size of prolongator (level->m)
 * \param coarse_r: vector(m), contains -1 if fine nodes, else the coordinate r of coarse nodes (level->coarse_nodes_list_r)
 * \param coarse_theta: idem with theta coordinates (level->coarse_nodes_list_theta)
 * \param trans: should we prolongate (P) or restrict (R=Pt)
 *
 */
std::vector<double> level::apply_prolongation_bi0(std::vector<double> u, int mc, int ncoarse, std::vector<int> coarse_r,
                                                  std::vector<int> coarse_theta, int trans)
{
    std::vector<double> Pu;
    if (trans == 0) { //Prolongation (P * u = Pu), (m x mc) * mc = m
        Pu = std::vector<double>(m, 0);
    }
    else if (trans == 1) { //Restriction (P^T * u = Pu), (mc x m) * m = mc
        Pu = std::vector<double>(mc, 0);
    }

    int index               = 0;
    int start_node_loop     = 0;
    int max_row_indices_act = 0;
    int coarse_index        = 0;

    int row_index;
    int col_index_base;
    int ri, ci;
    double vi;

    for (int k = start_node_loop; k < ncoarse; k++) {

        if (coarse_r[k] > -1) { // Coarse Node; no interpolation necessary

            row_index      = max_row_indices_act;
            col_index_base = coarse_index;
            ri             = row_index;
            ci             = col_index_base; // u(r,phi)
            vi             = 1; // 4/4

            if (trans) {
                Pu[ci] += vi * u[ri];
            }
            else
                Pu[ri] += vi * u[ci];
            index++;
            max_row_indices_act++;
            coarse_index++;
        }
        else { // only fine node
            // alt: r0=0/origin_NOT_coarse
            int phi_index_prev, r_index_1_prev, r_index_1_next;
            int coarse_nodes_prevCirc, coarse_nodes_nexCirc, ncoarse_2prevCirc;
            double phi_next, h_prev, h_next, k_prev, k_next;

            if (coarse_r[k - 1] > -1 &&
                ((k + 1 < ncoarse && coarse_r[k - 1] - coarse_r[k + 1] == 0) ||
                 (k - ntheta_int + 1 >= 0 && coarse_r[k - 1] - coarse_r[k - ntheta_int + 1] == 0))) {
                // previous (k-1) and following (either k+1 or k-nphi+1) nodes are coarse and have the same
                // radius r; interpolation between u(phi-k) and u(phi+k)
                // no periodicity workaround for phi-l since first node (r,0)
                // on each circle (where coarse nodes are present) is always coarse!
                // periodicity for phi+k handled below
                // r_index_1_bef = c(k-1,1)+1; % index in r (starting with 1); r_indices are 0 based otherwise.

                phi_index_prev = coarse_theta[k - 1]; // index one before phi_k
                if (k + 1 < ncoarse && coarse_r[k - 1] - coarse_r[k + 1] == 0)
                    phi_next = theta[coarse_theta[k + 1]]; // index one after phi_k
                else
                    phi_next = 2 * PI; // index one after phi_k (periodic bc come in)

                k_prev = theta[phi_index_prev + 1] - theta[phi_index_prev];
                k_next = phi_next - theta[phi_index_prev + 1];

                row_index = max_row_indices_act;

                // u(phi-k)
                col_index_base = coarse_index;

                ri = row_index;
                ci = col_index_base - 1; // u(r,phi-k) (voriger Coarse Node)
                vi = k_next / (k_prev + k_next); // 1/2

                max_row_indices_act = row_index + 1;
                if (trans)
                    Pu[ci] += vi * u[ri];
                else
                    Pu[ri] += vi * u[ci];
                index++;

                // u(phi+k)
                ri = row_index;
                if (k + 1 < ncoarse && coarse_r[k - 1] - coarse_r[k + 1] == 0)
                    ci = col_index_base; // u(r,phi+k) (next or currently 'considered' Coarse Node)
                else if (coarse_r[k - 1] - coarse_r[k - ntheta_int + 1] == 0)
                    ci = col_index_base - ceil(ntheta_int / 2); // periodicity
                else {
                    std::cout << "In build_bi_aniso_rDirphiPer (k=" << k << "): should not be entered.\n";
                    throw std::runtime_error("In build_bi_aniso_rDirphiPer: should not be entered.");
                }
                vi = k_prev / (k_prev + k_next); // 1/2

                if (trans)
                    Pu[ci] += vi * u[ri];
                else
                    Pu[ri] += vi * u[ci];
                index++;
            }
            // HERE
            else if (coarse_r[k - ntheta_int] > -1 &&
                     coarse_theta[k - ntheta_int] - coarse_theta[k + ntheta_int] == 0) {

                // corresp. coarse nodes with same phi are at r-h and r+h
                // interpolation between u(r-h) and u(r+h)
                // previous index in r (starting with 1); r_indices are 0 based otherwise.
                r_index_1_prev = coarse_r[k - ntheta_int];
                // next index in r (starting with 1); r_indices are 0 based otherwise.
                r_index_1_next = coarse_r[k + ntheta_int];

                if (r_index_1_prev < 0)
                    h_prev = r[r_index_1_prev + 1] - gyro::dcntl[Param::r0_DB];
                else
                    h_prev = r[r_index_1_prev + 1] - r[r_index_1_prev];
                if (r_index_1_next + 1 > nr_int)
                    h_next = gyro::dcntl[Param::R] - r[r_index_1_prev + 1];
                else
                    h_next = r[r_index_1_next] - r[r_index_1_prev + 1];

                row_index      = max_row_indices_act;
                col_index_base = coarse_index;

                // u(r-h)
                // do not introduce interactions between interior nodes and boundary nodes.
                // Correction of boundary to interior zero anyway (since residual is zero);
                // introduction of interior values to boundary nodes is not intended, so do not do it !
                if (r_index_1_prev > 0) {
                    coarse_nodes_prevCirc = 0;
                    for (int i = k - ntheta_int; i < k + 1; i++)
                        if (coarse_r[i] > -1)
                            coarse_nodes_prevCirc += 1;

                    ri = row_index;
                    ci = col_index_base - coarse_nodes_prevCirc; // Coarse Node u(r-h,phi)
                    // not touching boundary, standard case %% if/else with changed values should be bullshit; removing later
                    vi                  = h_next / (h_prev + h_next); // 1/2
                    max_row_indices_act = row_index + 1; // update max
                    if (trans)
                        Pu[ci] += vi * u[ri];
                    else
                        Pu[ri] += vi * u[ci];
                    index++;
                }

                // u(r+h)
                // do not introduce interactions between interior nodes and boundary nodes.
                // Correction of boundary to interior zero anyway (since residual is zero);
                // introduction of interior values to boundary nodes is not intended, so do not do it !
                if (r_index_1_next + 1 < nr_int) {
                    coarse_nodes_nexCirc = 0;
                    for (int i = k; i < k + ntheta_int + 1; i++)
                        if (coarse_r[i] > -1)
                            coarse_nodes_nexCirc += 1;
                    ri = row_index;
                    ci = col_index_base + coarse_nodes_nexCirc - 1; // u(r+h,phi)
                    // not touching boundary, standard case %% if/else with changed values should be bullshit; removing later
                    vi = h_prev / (h_prev + h_next); // 1/2

                    max_row_indices_act = row_index + 1; // update max
                    if (trans)
                        Pu[ci] += vi * u[ri];
                    else
                        Pu[ri] += vi * u[ci];
                    index++;
                }
                // weighing with 1/4 from all 'diagonally adjacent' coarse node
                // for standard coarsening this is in fact 0.5x the previous interpolated fine node
                // plus 0.5x the next (...taking periodicity into account!) interpolated fine node
                // ...thus take formula from previous elseif and change
                // correspondingly.
            }
            else {
                // previous index in r (starting with 1); r_indices are 0 based otherwise.

                r_index_1_prev = coarse_r[k - ntheta_int - 1];
                // next index in r (starting with 1); r_indices are 0 based otherwise
                // (take -1 instead of +1 (attention to periodic BC!) to not go one circle to far)
                r_index_1_next = coarse_r[k + ntheta_int - 1];

                // index one before phi_k  (access is possible via 'k-nphi-1' since subsequent if should never be entered)
                phi_index_prev = coarse_theta[k - ntheta_int - 1];

                // means that first node on r-circle with (r,phi)=(r,0) is a fine node and does not have adjacent coarse nodes at (r-h,0) and (r+h,0)
                // [jump across periodic bc had to be respected; not implemented; must and shall not appear with current coarsening]
                if (std::min(coarse_theta[k - ntheta_int - 1], coarse_theta[k - ntheta_int + 1]) > -1 &&
                    coarse_theta[k - ntheta_int - 1] > coarse_theta[k - ntheta_int + 1])
                    std::cout << "WARNING: Coarsening strategy has to be adapted for bilinear interpolation\n";
                if (phi_index_prev + 1 < ntheta_int - 1)
                    phi_next = theta[coarse_theta[k + ntheta_int + 1]]; // index one after phi_k
                else
                    phi_next = 2 * PI; // index one after phi_k (periodic bc come in)

                if (r_index_1_prev < 0)
                    h_prev = r[r_index_1_prev + 1] - gyro::dcntl[Param::r0_DB]; // h_{i-1}
                else
                    h_prev = r[r_index_1_prev + 1] - r[r_index_1_prev]; // h_{i-1}
                if (r_index_1_next + 1 > nr_int)
                    h_next = gyro::dcntl[Param::R] - r[r_index_1_prev + 1]; // h_{i}
                else
                    h_next = r[r_index_1_next] - r[r_index_1_prev + 1]; // h_{i}

                k_prev = theta[phi_index_prev + 1] - theta[phi_index_prev];
                k_next = phi_next - theta[phi_index_prev + 1];

                row_index      = max_row_indices_act;
                col_index_base = coarse_index;

                // do not introduce interactions between interior nodes and boundary nodes.
                // Correction of boundary to interior zero anyway (since residual is zero);
                // introduction of interior values to boundary nodes is not intended, so do not do it !
                if (r_index_1_prev > 0) {
                    // u(r-h,phi-h)
                    coarse_nodes_prevCirc = 0;
                    for (int i = k - ntheta_int - 1; i < k; i++)
                        if (coarse_r[i] > -1)
                            coarse_nodes_prevCirc += 1;

                    ri = row_index;
                    ci = col_index_base - coarse_nodes_prevCirc; // Coarse Node u(r-h,phi-h)
                    vi = h_next * k_next / ((k_prev + k_next) * (h_prev + h_next)); // isotrop: 1/4

                    max_row_indices_act = row_index + 1; // update max
                    if (trans)
                        Pu[ci] += vi * u[ri];
                    else
                        Pu[ri] += vi * u[ci];
                    index++;

                    // u(r-h,phi+h)
                    ri = row_index;
                    // coarse node at (r-h,phi-h) is at phi=nphi-1, so 'next' coarse
                    // node at (r-h,phi+h) is at (r-h,0) since phi=2pi-h
                    // %%%%%%%%%%%%%%%%%%%%%%%%% Hier Fehler k-2*nphi:k=0,....
                    if (coarse_theta[k - 1 - ntheta_int] + 1 == ntheta_int - 1) {
                        ncoarse_2prevCirc = 0;
                        for (int i = k - 2 * ntheta_int + 1; i < k - 1; i++)
                            if (coarse_r[i] > -1)
                                ncoarse_2prevCirc += 1;
                        ci = col_index_base - ncoarse_2prevCirc; // Coarse Node u(r-h,phi+h)
                    }
                    else
                        ci = col_index_base - coarse_nodes_prevCirc + 1; // Coarse Node u(r-h,phi+h)
                    // not touching boundary, standard case
                    vi = h_next * k_prev / ((k_prev + k_next) * (h_prev + h_next)); // isotrop: 1/4

                    if (trans)
                        Pu[ci] += vi * u[ri];
                    else
                        Pu[ri] += vi * u[ci];
                    index++;
                }

                // do not introduce interactions between interior nodes and boundary nodes.
                // Correction of boundary to interior zero anyway (since residual is zero);
                // introduction of interior values to boundary nodes is not intended, so do not do it !
                if (r_index_1_next + 1 < nr_int) {
                    // u(r+h,phi-h)
                    coarse_nodes_nexCirc = 0;
                    for (int i = k + 1; i < k + ntheta_int + 1; i++) {
                        if (coarse_r[i] > -1)
                            coarse_nodes_nexCirc += 1;
                    }

                    ri = row_index;
                    ci = col_index_base + coarse_nodes_nexCirc - 1; // u(r+h,phi-h)
                    vi = h_prev * k_next / ((k_prev + k_next) * (h_prev + h_next)); // isotrop: 1/4

                    max_row_indices_act = row_index + 1; // update max
                    if (trans)
                        Pu[ci] += vi * u[ri];
                    else
                        Pu[ri] += vi * u[ci];
                    index++;

                    // u(r+h,phi+h)
                    ri = row_index;
                    // 'next' coarse node at (..,phi+h) is at (..,0) since phi=2pi-h;
                    // corresp. index equals coarse_index since this is the index of the next coarse node
                    if (coarse_theta[k - 1 - ntheta_int] + 1 == ntheta_int - 1)
                        ci = col_index_base; // u(r+h,phi-h)
                    else
                        ci = col_index_base + coarse_nodes_nexCirc; // u(r+h,phi-h)
                    vi = h_prev * k_prev / ((k_prev + k_next) * (h_prev + h_next)); // isotrop: 1/4

                    if (trans)
                        Pu[ci] += vi * u[ri];
                    else
                        Pu[ri] += vi * u[ci];
                    index++;
                }
            }
        }
    }

    return Pu;
} /* ----- end of level::apply_prolongation_bi0 ----- */

/*!
 *  \brief Applies the injection
 *
 * Applies injection: directly inject values of coarse nodes.
 * Just uses values 1.0 for nodes of type 0.
 *
 * \param u: vector to prolongate/restrict
 * \param mc: coarse size of prolongator (level->mc)
 * \param ncoarse: fine size of prolongator (level->m)
 * \param coarse_r: vector(m), contains -1 if fine nodes, else the coordinate r of coarse nodes (level->coarse_nodes_list_r)
 * \param coarse_theta: idem with theta coordinates (level->coarse_nodes_list_theta)
 * \param trans: should we prolongate (P) or restrict (R=Pt)
 *
 */
std::vector<double> level::apply_prolongation_inj0(std::vector<double> u, int mc, int ncoarse,
                                                   std::vector<int> coarse_r, std::vector<int> coarse_theta, int trans)
{
    std::vector<double> Pu;
    if (trans == 0) { //Prolongation (P * u = Pu), (m x mc) * mc = m
        Pu = std::vector<double>(m, 0);
    }
    else if (trans == 1) { //Restriction (P^T * u = Pu), (mc x m) * m = mc
        Pu = std::vector<double>(mc, 0);
    }

    // Beware ! we only consider the m-th first elements in ncoarse_lists_*:
    // ncoarse_lists{i+1}(nodes_remain{i},:)

    int index = 0;
    // alt: r0=0/origin_NOT_coarse
    int start_node_loop     = 0;
    int max_row_indices_act = 0;
    int coarse_index        = 0;

    int row_index;
    int col_index_base;
    int ri, ci;
    double vi;

    for (int k = start_node_loop; k < ncoarse; k++) {
        //! CASE 1: Coarse Node; no interpolation necessary
        if (coarse_r[k] > -1) { // Coarse Node; no interpolation necessary

            row_index      = max_row_indices_act;
            col_index_base = coarse_index;
            ri             = row_index;
            ci             = col_index_base; // u(r,phi)
            vi             = 1.0; // 4/4

            if (trans)
                Pu[ci] += vi * u[ri];
            else
                Pu[ri] += vi * u[ci];
            index++;
            max_row_indices_act++;
            coarse_index++;
        }
        else { // only fine node
            // alt: r0=0/origin_NOT_coarse
            int phi_index_prev, r_index_1_prev, r_index_1_next;
            int coarse_nodes_prevCirc, coarse_nodes_nexCirc, ncoarse_2prevCirc;
            double phi_next, h_prev, h_next, k_prev, k_next;

            //! CASE 2
            if (coarse_r[k - 1] > -1 &&
                ((k + 1 < ncoarse && coarse_r[k - 1] - coarse_r[k + 1] == 0) ||
                 (k - ntheta_int + 1 >= 0 && coarse_r[k - 1] - coarse_r[k - ntheta_int + 1] == 0))) {
                // previous (k-1) and following (either k+1 or k-nphi+1) nodes are coarse and have the same
                // radius r; interpolation between u(phi-k) and u(phi+k)
                // no periodicity workaround for phi-l since first node (r,0)
                // on each circle (where coarse nodes are present) is always coarse!
                // periodicity for phi+k handled below
                // r_index_1_bef = c(k-1,1)+1; % index in r (starting with 1); r_indices are 0 based otherwise.

                phi_index_prev = coarse_theta[k - 1]; // index one before phi_k
                if (k + 1 < ncoarse && coarse_r[k - 1] - coarse_r[k + 1] == 0)
                    phi_next = theta[coarse_theta[k + 1]]; // index one after phi_k
                else
                    phi_next = 2 * PI; // index one after phi_k (periodic bc come in)

                k_prev = theta[phi_index_prev + 1] - theta[phi_index_prev];
                k_next = phi_next - theta[phi_index_prev + 1];

                row_index = max_row_indices_act;

                //! u(phi-k)
                col_index_base = coarse_index;

                ri = row_index;
                ci = col_index_base - 1; // u(r,phi-k) (voriger Coarse Node)
                vi = k_next / (k_prev + k_next); // 1/2

                max_row_indices_act = row_index + 1;
                index++;

                //! u(phi+k)
                ri = row_index;
                if (k + 1 < ncoarse && coarse_r[k - 1] - coarse_r[k + 1] == 0)
                    ci = col_index_base; // u(r,phi+k) (next or currently 'considered' Coarse Node)
                else if (coarse_r[k - 1] - coarse_r[k - ntheta_int + 1] == 0)
                    ci = col_index_base - ceil(ntheta_int / 2); // periodicity
                else {
                    std::cout << "In build_bi_aniso_rDirphiPer (k=" << k << "): should not be entered.\n";
                    throw std::runtime_error("In build_bi_aniso_rDirphiPer: should not be entered.");
                }
                vi = k_prev / (k_prev + k_next); // 1/2

                index++;
            }
            //! CASE 3
            else if (coarse_r[k - ntheta_int] > -1 &&
                     coarse_theta[k - ntheta_int] - coarse_theta[k + ntheta_int] == 0) {

                // corresp. coarse nodes with same phi are at r-h and r+h
                // interpolation between u(r-h) and u(r+h)
                // previous index in r (starting with 1); r_indices are 0 based otherwise.
                r_index_1_prev = coarse_r[k - ntheta_int];
                // next index in r (starting with 1); r_indices are 0 based otherwise.
                r_index_1_next = coarse_r[k + ntheta_int];

                if (r_index_1_prev < 0)
                    h_prev = r[r_index_1_prev + 1] - gyro::dcntl[Param::r0_DB];
                else
                    h_prev = r[r_index_1_prev + 1] - r[r_index_1_prev];
                if (r_index_1_next + 1 > nr_int)
                    h_next = gyro::dcntl[Param::R] - r[r_index_1_prev + 1];
                else
                    h_next = r[r_index_1_next] - r[r_index_1_prev + 1];

                row_index      = max_row_indices_act;
                col_index_base = coarse_index;

                //! u(r-h)
                // do not introduce interactions between interior nodes and boundary nodes.
                // Correction of boundary to interior zero anyway (since residual is zero);
                // introduction of interior values to boundary nodes is not intended, so do not do it !
                if (r_index_1_prev > 0) {
                    coarse_nodes_prevCirc = 0;
                    for (int i = k - ntheta_int; i < k + 1; i++)
                        if (coarse_r[i] > -1)
                            coarse_nodes_prevCirc += 1;
                    ri = row_index;
                    ci = col_index_base - coarse_nodes_prevCirc; // Coarse Node u(r-h,phi)
                    // not touching boundary, standard case %% if/else with changed values should be bullshit; removing later
                    vi                  = h_next / (h_prev + h_next); // 1/2
                    max_row_indices_act = row_index + 1; // update max
                    index++;
                }

                //! u(r+h)
                // do not introduce interactions between interior nodes and boundary nodes.
                // Correction of boundary to interior zero anyway (since residual is zero);
                // introduction of interior values to boundary nodes is not intended, so do not do it !
                if (r_index_1_next + 1 < nr_int) {
                    coarse_nodes_nexCirc = 0;
                    for (int i = k; i < k + ntheta_int + 1; i++)
                        if (coarse_r[i] > -1)
                            coarse_nodes_nexCirc += 1;
                    ri = row_index;
                    ci = col_index_base + coarse_nodes_nexCirc - 1; // u(r+h,phi)
                    // not touching boundary, standard case %% if/else with changed values should be bullshit; removing later
                    vi = h_prev / (h_prev + h_next); // 1/2

                    max_row_indices_act = row_index + 1; // update max
                    index++;
                }
                // weighing with 1/4 from all 'diagonally adjacent' coarse node
                // for standard coarsening this is in fact 0.5x the previous interpolated fine node
                // plus 0.5x the next (...taking periodicity into account!) interpolated fine node
                // ...thus take formula from previous elseif and change
                // correspondingly.
            }
            //! CASE 4
            else {
                // previous index in r (starting with 1); r_indices are 0 based otherwise.

                r_index_1_prev = coarse_r[k - ntheta_int - 1];
                // next index in r (starting with 1); r_indices are 0 based otherwise
                // (take -1 instead of +1 (attention to periodic BC!) to not go one circle to far)
                r_index_1_next = coarse_r[k + ntheta_int - 1];

                // index one before phi_k  (access is possible via 'k-nphi-1' since subsequent if should never be entered)
                phi_index_prev = coarse_theta[k - ntheta_int - 1];

                // means that first node on r-circle with (r,phi)=(r,0) is a fine node and does not have adjacent coarse nodes at (r-h,0) and (r+h,0)
                // [jump across periodic bc had to be respected; not implemented; must and shall not appear with current coarsening]
                if (std::min(coarse_theta[k - ntheta_int - 1], coarse_theta[k - ntheta_int + 1]) > -1 &&
                    coarse_theta[k - ntheta_int - 1] > coarse_theta[k - ntheta_int + 1])
                    std::cout << "WARNING: Coarsening strategy has to be adapted for bilinear interpolation\n";
                if (phi_index_prev + 1 < ntheta_int - 1)
                    phi_next = theta[coarse_theta[k + ntheta_int + 1]]; // index one after phi_k
                else
                    phi_next = 2 * PI; // index one after phi_k (periodic bc come in)

                if (r_index_1_prev < 0)
                    h_prev = r[r_index_1_prev + 1] - gyro::dcntl[Param::r0_DB]; // h_{i-1}
                else
                    h_prev = r[r_index_1_prev + 1] - r[r_index_1_prev]; // h_{i-1}
                if (r_index_1_next + 1 > nr_int)
                    h_next = gyro::dcntl[Param::R] - r[r_index_1_prev + 1]; // h_{i}
                else
                    h_next = r[r_index_1_next] - r[r_index_1_prev + 1]; // h_{i}

                k_prev = theta[phi_index_prev + 1] - theta[phi_index_prev];
                k_next = phi_next - theta[phi_index_prev + 1];

                row_index      = max_row_indices_act;
                col_index_base = coarse_index;

                // do not introduce interactions between interior nodes and boundary nodes.
                // Correction of boundary to interior zero anyway (since residual is zero);
                // introduction of interior values to boundary nodes is not intended, so do not do it !
                if (r_index_1_prev > 0) {
                    // u(r-h,phi-h)
                    coarse_nodes_prevCirc = 0;
                    for (int i = k - ntheta_int - 1; i < k; i++)
                        if (coarse_r[i] > -1)
                            coarse_nodes_prevCirc += 1;

                    ri = row_index;
                    ci = col_index_base - coarse_nodes_prevCirc; // Coarse Node u(r-h,phi-h)
                    vi = h_next * k_next / ((k_prev + k_next) * (h_prev + h_next)); // isotrop: 1/4

                    max_row_indices_act = row_index + 1; // update max
                    index++;

                    // u(r-h,phi+h)
                    ri = row_index;
                    // coarse node at (r-h,phi-h) is at phi=nphi-1, so 'next' coarse
                    // node at (r-h,phi+h) is at (r-h,0) since phi=2pi-h
                    // %%%%%%%%%%%%%%%%%%%%%%%%% Hier Fehler k-2*nphi:k=0,....
                    if (coarse_theta[k - 1 - ntheta_int] + 1 == ntheta_int - 1) {
                        ncoarse_2prevCirc = 0;
                        for (int i = k - 2 * ntheta_int + 1; i < k - 1; i++)
                            if (coarse_r[i] > -1)
                                ncoarse_2prevCirc += 1;
                        ci = col_index_base - ncoarse_2prevCirc; // Coarse Node u(r-h,phi+h)
                    }
                    else
                        ci = col_index_base - coarse_nodes_prevCirc + 1; // Coarse Node u(r-h,phi+h)
                    // not touching boundary, standard case
                    vi = h_next * k_prev / ((k_prev + k_next) * (h_prev + h_next)); // isotrop: 1/4

                    index++;
                }

                // do not introduce interactions between interior nodes and boundary nodes.
                // Correction of boundary to interior zero anyway (since residual is zero);
                // introduction of interior values to boundary nodes is not intended, so do not do it !
                if (r_index_1_next + 1 < nr_int) {
                    // u(r+h,phi-h)
                    coarse_nodes_nexCirc = 0;
                    for (int i = k + 1; i < k + ntheta_int + 1; i++) {
                        if (coarse_r[i] > -1)
                            coarse_nodes_nexCirc += 1;
                    }

                    ri = row_index;
                    ci = col_index_base + coarse_nodes_nexCirc - 1; // u(r+h,phi-h)
                    vi = h_prev * k_next / ((k_prev + k_next) * (h_prev + h_next)); // isotrop: 1/4

                    max_row_indices_act = row_index + 1; // update max
                    index++;

                    // u(r+h,phi+h)
                    ri = row_index;
                    // 'next' coarse node at (..,phi+h) is at (..,0) since phi=2pi-h;
                    // corresp. index equals coarse_index since this is the index of the next coarse node
                    if (coarse_theta[k - 1 - ntheta_int] + 1 == ntheta_int - 1)
                        ci = col_index_base; // u(r+h,phi-h)
                    else
                        ci = col_index_base + coarse_nodes_nexCirc; // u(r+h,phi-h)
                    vi = h_prev * k_prev / ((k_prev + k_next) * (h_prev + h_next)); // isotrop: 1/4
                    index++;
                }
            }
        }
    }

    //std::cout << "end of applyP_inj \n";

    return Pu;
} /* ----- end of level::apply_prolongation_inj0 ----- */

/*!
 *  \brief Applies the extrapolation operator
 *
 * Applies the prolongation operator for implicit extrapolation.
 *
 * The stencil is the same as the isotropic bilinear interpolation stencil,
 * except for fine points with diagonal coarse neighboors (types 5,6,7) where only
 * top_left and top_right values are used.
 *
 * \param u: vector to prolongate/restrict
 * \param mv: coarse size of prolongator (level->mc)
 * \param ncoarse: fine size of prolongator (level->m)
 * \param coarse_r: vector(m), contains -1 if fine nodes, else the coordinate r of coarse nodes (level->coarse_nodes_list_r)
 * \param coarse_theta: idem with theta coordinates (level->coarse_nodes_list_theta)
 * \param trans: should we prolongate (P) or restrict (R=Pt)
 *
 */
std::vector<double> level::apply_prolongation_ex0(std::vector<double> u, int mc, int ncoarse, std::vector<int> coarse_r,
                                                  std::vector<int> coarse_theta, int trans)
{
    std::vector<double> Pu;
    if (trans == 0) { //Prolongation (P * u = Pu), (m x mc) * mc = m
        Pu = std::vector<double>(m, 0);
    }
    else if (trans == 1) { //Restriction (P^T * u = Pu), (mc x m) * m = mc
        Pu = std::vector<double>(mc, 0);
    }

    // Beware ! we only consider the m-th first elements in ncoarse_lists_*:
    // ncoarse_lists{i+1}(nodes_remain{i},:)

    int index               = 0;
    int start_node_loop     = 0;
    int max_row_indices_act = 0;
    int coarse_index        = 0;

    int row_index;
    int col_index_base;
    int ri, ci;
    double vi;

    for (int k = start_node_loop; k < ncoarse; k++) {
        //! CASE 1: Coarse Node; no interpolation necessary (coarse node = fine node)
        if (coarse_r[k] > -1) { // Coarse Node; no interpolation necessary

            row_index      = max_row_indices_act;
            col_index_base = coarse_index;
            ri             = row_index;
            ci             = col_index_base; // u(r,phi)
            vi             = 1.0; // 4/4

            if (trans)
                Pu[ci] += vi * u[ri];
            else
                Pu[ri] += vi * u[ci];
            index++;
            max_row_indices_act++;
            coarse_index++;
        }
        else { // only fine nodes
            //check if the fine node has a diagonal link to a coarse node
            int diagonal_link = 0; //indicates wheather a fine node has a diagonal link or not
            int i             = k % ntheta_int;
            int j             = floor(k / ntheta_int);
            if (i % 2 == 1 && j % 2 == 1)
                diagonal_link = 1;

            // alt: r0=0/origin_NOT_coarse
            int phi_index_prev, r_index_1_prev, r_index_1_next;
            int coarse_nodes_prevCirc, coarse_nodes_nexCirc, ncoarse_2prevCirc;
            double phi_next, h_prev, h_next, k_prev, k_next;

            //! CASE 2: fine nodes which have coarse points to the right and to the left
            if (coarse_r[k - 1] > -1 &&
                ((k + 1 < ncoarse && coarse_r[k - 1] - coarse_r[k + 1] == 0) ||
                 (k - ntheta_int + 1 >= 0 && coarse_r[k - 1] - coarse_r[k - ntheta_int + 1] == 0))) {
                // previous (k-1) and following (either k+1 or k-nphi+1) nodes are coarse and have the same
                // radius r; interpolation between u(phi-k) and u(phi+k)
                // no periodicity workaround for phi-l since first node (r,0)
                // on each circle (where coarse nodes are present) is always coarse!
                // periodicity for phi+k handled below
                // r_index_1_bef = c(k-1,1)+1; % index in r (starting with 1); r_indices are 0 based otherwise.

                phi_index_prev = coarse_theta[k - 1]; // index one before phi_k
                if (k + 1 < ncoarse && coarse_r[k - 1] - coarse_r[k + 1] == 0)
                    phi_next = theta[coarse_theta[k + 1]]; // index one after phi_k
                else
                    phi_next = 2 * PI; // index one after phi_k (periodic bc come in)

                k_prev = theta[phi_index_prev + 1] - theta[phi_index_prev];
                k_next = phi_next - theta[phi_index_prev + 1];

                row_index = max_row_indices_act;

                //! u(phi-k) (left)
                col_index_base = coarse_index;

                ri = row_index;
                ci = col_index_base - 1; // u(r,phi-k) (voriger Coarse Node)
                vi = 0.5;

                max_row_indices_act = row_index + 1;
                if (diagonal_link == 0) {
                    if (trans)
                        Pu[ci] += vi * u[ri];
                    else
                        Pu[ri] += vi * u[ci];
                }
                index++;

                // u(phi+k) (right)
                ri = row_index;
                if (k + 1 < ncoarse && coarse_r[k - 1] - coarse_r[k + 1] == 0)
                    ci = col_index_base; // u(r,phi+k) (next or currently 'considered' Coarse Node)
                else if (coarse_r[k - 1] - coarse_r[k - ntheta_int + 1] == 0)
                    ci = col_index_base - ceil(ntheta_int / 2); // periodicity
                else {
                    std::cout << "In build_bi_aniso_rDirphiPer (k=" << k << "): should not be entered.\n";
                    throw std::runtime_error("In build_bi_aniso_rDirphiPer: should not be entered.");
                }
                vi = 0.5;

                if (diagonal_link == 0) {
                    if (trans)
                        Pu[ci] += vi * u[ri];
                    else
                        Pu[ri] += vi * u[ci];
                }
                index++;
            }
            //! CASE 3: fine nodes which have coarse points at the top and bottom side
            else if (coarse_r[k - ntheta_int] > -1 &&
                     coarse_theta[k - ntheta_int] - coarse_theta[k + ntheta_int] == 0) {

                // corresp. coarse nodes with same phi are at r-h and r+h
                // interpolation between u(r-h) and u(r+h)
                // previous index in r (starting with 1); r_indices are 0 based otherwise.
                r_index_1_prev = coarse_r[k - ntheta_int];
                // next index in r (starting with 1); r_indices are 0 based otherwise.
                r_index_1_next = coarse_r[k + ntheta_int];

                if (r_index_1_prev < 0)
                    h_prev = r[r_index_1_prev + 1] - gyro::dcntl[Param::r0_DB];
                else
                    h_prev = r[r_index_1_prev + 1] - r[r_index_1_prev];
                if (r_index_1_next + 1 > nr_int)
                    h_next = gyro::dcntl[Param::R] - r[r_index_1_prev + 1];
                else
                    h_next = r[r_index_1_next] - r[r_index_1_prev + 1];

                row_index      = max_row_indices_act;
                col_index_base = coarse_index;

                //! u(r-h) (bottom)
                // do not introduce interactions between interior nodes and boundary nodes.
                // Correction of boundary to interior zero anyway (since residual is zero);
                // introduction of interior values to boundary nodes is not intended, so do not do it !
                if (r_index_1_prev > 0) {
                    coarse_nodes_prevCirc = 0;
                    for (int i = k - ntheta_int; i < k + 1; i++)
                        if (coarse_r[i] > -1)
                            coarse_nodes_prevCirc += 1;

                    ri = row_index;
                    ci = col_index_base - coarse_nodes_prevCirc; // Coarse Node u(r-h,phi)
                    vi = 0.5;

                    max_row_indices_act = row_index + 1; // update max
                    if (diagonal_link == 0) {
                        if (trans)
                            Pu[ci] += vi * u[ri];
                        else
                            Pu[ri] += vi * u[ci];
                    }
                    index++;
                }

                //! u(r+h) (top)
                // do not introduce interactions between interior nodes and boundary nodes.
                // Correction of boundary to interior zero anyway (since residual is zero);
                // introduction of interior values to boundary nodes is not intended, so do not do it !
                if (r_index_1_next + 1 < nr_int) {
                    coarse_nodes_nexCirc = 0;
                    for (int i = k; i < k + ntheta_int + 1; i++)
                        if (coarse_r[i] > -1)
                            coarse_nodes_nexCirc += 1;
                    ri = row_index;
                    ci = col_index_base + coarse_nodes_nexCirc - 1; // u(r+h,phi)
                    vi = h_prev / (h_prev + h_next); // 1/2
                    vi = 0.5;

                    max_row_indices_act = row_index + 1; // update max
                    if (diagonal_link == 0) {
                        if (trans)
                            Pu[ci] += vi * u[ri];
                        else
                            Pu[ri] += vi * u[ci];
                    }
                    index++;
                }
                // weighing with 1/4 from all 'diagonally adjacent' coarse node
                // for standard coarsening this is in fact 0.5x the previous interpolated fine node
                // plus 0.5x the next (...taking periodicity into account!) interpolated fine node
                // ...thus take formula from previous elseif and change
                // correspondingly.
            }
            //! CASE 4: fine nodes which have coarse points only diagonally
            else {
                // previous index in r (starting with 1); r_indices are 0 based otherwise.

                r_index_1_prev = coarse_r[k - ntheta_int - 1];
                // next index in r (starting with 1); r_indices are 0 based otherwise
                // (take -1 instead of +1 (attention to periodic BC!) to not go one circle to far)
                r_index_1_next = coarse_r[k + ntheta_int - 1];

                // index one before phi_k  (access is possible via 'k-nphi-1' since subsequent if should never be entered)
                phi_index_prev = coarse_theta[k - ntheta_int - 1];

                // means that first node on r-circle with (r,phi)=(r,0) is a fine node and does not have adjacent coarse nodes at (r-h,0) and (r+h,0)
                // [jump across periodic bc had to be respected; not implemented; must and shall not appear with current coarsening]
                if (std::min(coarse_theta[k - ntheta_int - 1], coarse_theta[k - ntheta_int + 1]) > -1 &&
                    coarse_theta[k - ntheta_int - 1] > coarse_theta[k - ntheta_int + 1])
                    std::cout << "WARNING: Coarsening strategy has to be adapted for bilinear interpolation\n";
                if (phi_index_prev + 1 < ntheta_int - 1)
                    phi_next = theta[coarse_theta[k + ntheta_int + 1]]; // index one after phi_k
                else
                    phi_next = 2 * PI; // index one after phi_k (periodic bc come in)

                if (r_index_1_prev < 0)
                    h_prev = r[r_index_1_prev + 1] - gyro::dcntl[Param::r0_DB]; // h_{i-1}
                else
                    h_prev = r[r_index_1_prev + 1] - r[r_index_1_prev]; // h_{i-1}
                if (r_index_1_next + 1 > nr_int)
                    h_next = gyro::dcntl[Param::R] - r[r_index_1_prev + 1]; // h_{i}
                else
                    h_next = r[r_index_1_next] - r[r_index_1_prev + 1]; // h_{i}

                k_prev = theta[phi_index_prev + 1] - theta[phi_index_prev];
                k_next = phi_next - theta[phi_index_prev + 1];

                row_index      = max_row_indices_act;
                col_index_base = coarse_index;

                // do not introduce interactions between interior nodes and boundary nodes.
                // Correction of boundary to interior zero anyway (since residual is zero);
                // introduction of interior values to boundary nodes is not intended, so do not do it !
                if (r_index_1_prev > 0) {
                    //! u(r-h,phi-h) (bottom left)
                    coarse_nodes_prevCirc = 0;
                    for (int i = k - ntheta_int - 1; i < k; i++)
                        if (coarse_r[i] > -1)
                            coarse_nodes_prevCirc += 1;

                    ri = row_index;
                    ci = col_index_base - coarse_nodes_prevCirc; // Coarse Node u(r-h,phi-h)
                    vi = 0.5;

                    max_row_indices_act = row_index + 1; // update max
                    if (diagonal_link == -1) { //should not be treated
                        if (trans)
                            Pu[ci] += vi * u[ri];
                        else
                            Pu[ri] += vi * u[ci];
                    }
                    index++;

                    //! u(r-h,phi+h) (bottom right)
                    ri = row_index;
                    // coarse node at (r-h,phi-h) is at phi=nphi-1, so 'next' coarse
                    // node at (r-h,phi+h) is at (r-h,0) since phi=2pi-h
                    // %%%%%%%%%%%%%%%%%%%%%%%%% Hier Fehler k-2*nphi:k=0,....
                    if (coarse_theta[k - 1 - ntheta_int] + 1 == ntheta_int - 1) {
                        ncoarse_2prevCirc = 0;
                        for (int i = k - 2 * ntheta_int + 1; i < k - 1; i++)
                            if (coarse_r[i] > -1)
                                ncoarse_2prevCirc += 1;
                        ci = col_index_base - ncoarse_2prevCirc; // Coarse Node u(r-h,phi+h)
                    }
                    else
                        ci = col_index_base - coarse_nodes_prevCirc + 1; // Coarse Node u(r-h,phi+h)
                    // not touching boundary, standard case
                    vi = 0.5;

                    if (diagonal_link == 1) {
                        if (trans) {
                            Pu[ci] += vi * u[ri];
                        }
                        else {
                            Pu[ri] += vi * u[ci];
                        }
                    }
                    index++;
                }

                // do not introduce interactions between interior nodes and boundary nodes.
                // Correction of boundary to interior zero anyway (since residual is zero);
                // introduction of interior values to boundary nodes is not intended, so do not do it !
                if (r_index_1_next + 1 < nr_int) {
                    //! u(r+h,phi-h) (top left)
                    coarse_nodes_nexCirc = 0;
                    for (int i = k + 1; i < k + ntheta_int + 1; i++) {
                        if (coarse_r[i] > -1)
                            coarse_nodes_nexCirc += 1;
                    }

                    ri = row_index;
                    ci = col_index_base + coarse_nodes_nexCirc - 1; // u(r+h,phi-h)
                    vi = 0.5;

                    max_row_indices_act = row_index + 1; // update max
                    if (diagonal_link == 1) {
                        if (trans) {
                            Pu[ci] += vi * u[ri];
                        }
                        else {
                            Pu[ri] += vi * u[ci];
                        }
                    }
                    index++;

                    //! u(r+h,phi+h) (top right)
                    ri = row_index;
                    // 'next' coarse node at (..,phi+h) is at (..,0) since phi=2pi-h;
                    // corresp. index equals coarse_index since this is the index of the next coarse node
                    if (coarse_theta[k - 1 - ntheta_int] + 1 == ntheta_int - 1)
                        ci = col_index_base; // u(r+h,phi-h)
                    else
                        ci = col_index_base + coarse_nodes_nexCirc; // u(r+h,phi-h)
                    vi = 0.5;

                    if (diagonal_link == -1) { //should not be treated
                        if (trans)
                            Pu[ci] += vi * u[ri];
                        else
                            Pu[ri] += vi * u[ci];
                    }
                    index++;
                }
            }
        }
    }

    return Pu;
} /* ----- end of level::apply_prolongation_ex0 ----- */