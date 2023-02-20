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
 * \file build_fd9star_anisotr_scaled.cpp
 * \brief Implementation of the 9p FD operators A and RHS (deprecated)
 * \author M. Kuehn, C. Kruse, P. Leleux, C. Schwarz
 * \version 0.0
 */
#include "level.h"

/*!
 *  \brief Build the operator and RHS (deprecated)
 *
 * Builds the matrix A and the RHS of the system based on 9p FD. Relations with 
 * Dirichlet boundary condition nodes are shifted to the RHS to keep a symmetric 
 * operator.
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
 */
void level::build_A0()
{
    // if boundary is only defined on radius, checks for wrongly detected boundary
    // nodes due to rounding errors (double checking...)
    double tol_bound_check = gyro::dcntl[Param::tol_bound_check];
    // double x, y;
    int index;

    // Take boundary condition into account
    // In the operator
    // gyro::trafo(r[0], theta[0], x, y, 0);
    // alt: r0=0/neumBound
    // dirBound([r(1),0])==0 means that (r(1),0) is on Dirichlet boundary
    if (fabs(gyro::distBoundary(r[0], theta[0], 0)) < tol_bound_check) {
        // Dirichlet-RB (intérieur ; premier noeud dans la numérotation des disques)
        for (int i = 0; i < ntheta_int; i++) {
            row_indices.push_back(i);
            col_indices.push_back(i);
            vals.push_back(1.0);
        }
        index = ntheta_int;
    }
    else {
        index = 0;
    }

    // gyro::trafo(r[0], theta[0], x, y, 0);
    // alt: r0=0/neumBound
    int start_j = 0;
    // dirBound([r(1),0])==0 means that (r(1),0) is on Dirichlet boundary
    if (fabs(gyro::distBoundary(r[0], theta[0], 0)) < tol_bound_check) {
        start_j = 1;
    }
    // dirBound([r(1),0])>0 means that (r(1),0) is not on Dirichlet boundary
    else if (r[0] > 0 && fabs(gyro::distBoundary(r[0], theta[0], 0)) > tol_bound_check) {
        start_j = 0;
    }

    // int i1 = 0, i2 = 0;
    for (int j = start_j; j < nr_int; j++) { // nodes of Disk K(eps,1-eps)
        for (int i = 0; i < ntheta_int; i++) {
            // alt: r0=0
            int row_index = j * ntheta_int + i;

            // % ==================== %
            double kt = thetaplus[i];
            double thetamin1, ktmin1, hs, hsmin1;
            if (i > 0) {
                thetamin1 = theta[i - 1];
                ktmin1    = thetaplus[i - 1];
            }
            else {
                thetamin1 = theta[ntheta_int - 1]; // equal to 2*PI-ktmin1
                ktmin1    = thetaplus[ntheta_int - 1];
            }
            // % -------------------- %
            // r(j+1) is the current r
            hs = hplus[j];
            if (j > 0)
                hsmin1 = hplus[j - 1];
            else
                hsmin1 = 2 * r[0]; // across the origin
            // % ==================== %

            // 9-Point Stencil (attention order; first bottom, then right, then top right, top left, bottom left, middle)
            // r-h- (h- = x_i-x_i-1)
            // if (j == 1)
            //     gyro::trafo(r[j - 1], theta[i], x, y, 0);

            // j=1 means r-h is on the boundary, dirBound([r(j),0])==0 means that that Dirichlet BC are set at (r(j),0))
            // --> to symmetrize, put it on the right hand side
            if (j == 1 && r[j - 1] > 0 && fabs(gyro::distBoundary(r[j - 1], theta[i], 0)) < tol_bound_check) {
            }
            else {
                // row_indices[index] = row_index;
                row_indices.push_back(row_index);
                col_indices.push_back(-42);
                vals.push_back(-42);
                // second circle (at least one 'real' circle with ntheta_int nodes which is closer to the origin)
                if (j > 1 || (j == 1 && r[0] > 0))
                    // to u at r-h geometrically recorded (i.e. -ntheta_int nodes at slice sorting)
                    col_indices[index] = row_index - ntheta_int;
                else { // j=0 or j=1
                    if (r[0] == 0) { // all refering to the single origin node! (j=0 not possible in this case)
                        // Reference to u at r-h=0 geometrically recorded (for all i the origin, because here r=h)
                        col_indices[index] = 1;
                    }
                    else if (j == 0) { // across the origin
                        if (i + 1 > ntheta_int / 2) {
                            col_indices[index] = row_index - ntheta_int / 2; // half a turn back
                        }
                        else {
                            col_indices[index] = row_index + ntheta_int / 2; // half a turn further
                        }
                    }
                }

                if (j > 0) {
                    // alt: neumBC
                    // a l interieur sans contact aux conditions aux limites
                    vals[index] =
                        -0.5 * kt / hsmin1 * (gyro::arr(r[j - 1], theta[i], 0) + gyro::arr(r[j], theta[i], 0)) -
                        0.5 * ktmin1 / hsmin1 * (gyro::arr(r[j - 1], theta[i], 0) + gyro::arr(r[j], theta[i], 0));
                }
                else {
                    // j==0; across the origin
                    // just use r_{s-1}=0
                    // Use r(0):=r(j) to go to the other side of the origin
                    vals[index] =
                        -0.5 * kt / hsmin1 * (gyro::arr(r[j], theta[i] + PI, 0) + gyro::arr(r[j], theta[i], 0)) -
                        0.5 * ktmin1 / hsmin1 * (gyro::arr(r[j], theta[i] + PI, 0) + gyro::arr(r[j], theta[i], 0));
                }
                index++;
            }

            // (r-h-,phi-k) (en bas a gauche)
            if (gyro::icntl[Param::mod_pk] > 0) {
                // if (j == 1) // coords necessary for potential boundary conditions
                //     gyro::trafo(r[j - 1], thetamin1, x, y, 0); // r(j) is PREVIOUS r, actual is r(J+1) !
                // j=1 means r-h is on the boundary, dirBound([r(j),0])==0 means that that Dirichlet BC are set at (r(j),0))
                // --> to symmetrize, put it on the right hand side
                if (j == 1 && r[j - 1] > 0 && fabs(gyro::distBoundary(r[j - 1], thetamin1, 0)) < tol_bound_check) {
                }
                else {
                    // second circle (at least one 'real' circle with ntheta_int nodes which is closer to the origin)
                    // but no dirichlet bc on innermost circle
                    if (j > 1 || (j == 1 && r[0] > 0)) {
                        row_indices.push_back(row_index);
                        if (i > 0) // next node in theta direction but one circle before
                            col_indices.push_back(row_index - ntheta_int - 1);
                        else // periodicity condition
                            col_indices.push_back(row_index - 1);
                    }

                    //    if j > 0
                    // alt: neumBound
                    if (j > 1 || (j == 1 && r[0] > 0)) {
                        // a l interieur sans contact aux conditions aux limites
                        //         if (j > 1) // pas de reference vers l origine ici
                        vals.push_back(-(gyro::art(r[j], thetamin1, 0) + gyro::art(r[j - 1], theta[i], 0)));
                        index++;
                    }
                }
            }

            // phi-k
            // row_indices[index] = row_index;
            row_indices.push_back(row_index);
            if (i > 0) // previous node in phi direction
                col_indices.push_back(row_index - 1);
            else // periodicity condition
                col_indices.push_back(row_index + ntheta_int - 1);
            vals.push_back(-0.5 * hsmin1 / ktmin1 * (gyro::att(r[j], thetamin1, 0) + gyro::att(r[j], theta[i], 0)) -
                           0.5 * hs / ktmin1 * (gyro::att(r[j], thetamin1, 0) + gyro::att(r[j], theta[i], 0)));
            index++;

            // (r+h+,phi-k) (en haut a droite)
            if (gyro::icntl[Param::mod_pk] > 0) {
                double theta_eval_prelast = 0;
                double r_tmp, theta_tmp;
                if (i > 0) { // not the first node on the corresponding circle; we can take theta(i-1)
                    // gyro::trafo(r[j + 1], theta[i - 1], x, y, 0);
                    r_tmp     = r[j + 1];
                    theta_tmp = theta[i - 1];
                }
                else {
                    theta_eval_prelast = theta[ntheta_int - 1];
                    // gyro::trafo(r[j + 1], theta_eval_prelast, x, y, 0);
                    r_tmp     = r[j + 1];
                    theta_tmp = theta_eval_prelast;
                }

                // means that r[j+1] is not on the Dirichlet boundary
                if (j < nr_int - 1 ||
                    (j == nr_int - 1 && fabs(gyro::distBoundary(r_tmp, theta_tmp, 0)) > tol_bound_check)) {
                    // row_indices[index] = row_index;
                    row_indices.push_back(row_index);
                    if (i > 0) // previous node in phi direction but at r+h
                        col_indices.push_back(row_index + ntheta_int - 1);
                    else // first node on corresponding circle; pay gyro::attention to periodicity condition!
                        col_indices.push_back(row_index + 2 * ntheta_int - 1);
                    vals.push_back(gyro::art(r[j], thetamin1, 0) + gyro::art(r[j + 1], theta[i], 0));
                    index++;
                }
            }

            // // r + h+ (h+ = x_i+1-x_i)
            // gyro::trafo(r[j + 1], theta[i], x, y, 0);
            // means that r[j+1] is not on the Dirichlet boundary
            if (j < nr_int - 1 ||
                (j == nr_int - 1 && fabs(gyro::distBoundary(r[j + 1], theta[i], 0)) > tol_bound_check)) {
                row_indices.push_back(row_index);
                col_indices.push_back(row_index + ntheta_int); // zu u bei r+h
                vals.push_back(-0.5 * kt / hs * (gyro::arr(r[j], theta[i], 0) + gyro::arr(r[j + 1], theta[i], 0)) -
                               0.5 * ktmin1 / hs * (gyro::arr(r[j], theta[i], 0) + gyro::arr(r[j + 1], theta[i], 0)));
                index++;
            }

            // (r+h+,phi+k) (en haut a gauche)
            if (gyro::icntl[Param::mod_pk] > 0) {
                double thetap1 = 0;
                if (i + 1 < ntheta_int) // Define theta(i+1), this might be on the periodicity condition.
                    thetap1 = theta[i + 1];
                else
                    thetap1 = 2 * PI;

                // gyro::trafo(r[j + 1], thetap1, x, y, 0); ???
                // means that r[j+1] is not on the Dirichlet boundary
                if (j < nr_int - 1 ||
                    (j == nr_int - 1 && fabs(gyro::distBoundary(r[j + 1], theta[i], 0)) > tol_bound_check)) {
                    row_indices.push_back(row_index);
                    if (i + 1 < ntheta_int) // previous node in phi direction but at r+h
                        col_indices.push_back(row_index + ntheta_int + 1);
                    else // first node on corresponding circle; pay gyro::attention to periodicity condition!
                        col_indices.push_back(row_index + 1);
                    vals.push_back(-gyro::art(r[j + 1], theta[i], 0) - gyro::art(r[j], thetap1, 0));
                    index++;
                }
            }

            // phi+k
            double thetap1, r_tmp, theta_tmp;
            // row_indices[index] = row_index;
            row_indices.push_back(row_index);
            if (i + 1 < ntheta_int) // next node in theta direction
                col_indices.push_back(row_index + 1);
            else // periodicity condition
                col_indices.push_back(row_index - ntheta_int + 1);
            if (i + 1 < ntheta_int) // Define theta(i+1), this might be on the periodicity condition.
                thetap1 = theta[i + 1];
            else
                thetap1 = 2 * PI;
            vals.push_back(-0.5 * hs / kt * (gyro::att(r[j], theta[i], 0) + gyro::att(r[j], thetap1, 0)) -
                           0.5 * hsmin1 / kt * (gyro::att(r[j], theta[i], 0) + gyro::att(r[j], thetap1, 0)));
            index++;

            // (r-h-,phi+k) (en bas a gauche)
            if (gyro::icntl[Param::mod_pk] > 0) {
                if (j == 1) { // coords necessary for potential boundary conditions
                    if (i + 1 < ntheta_int) {
                        // gyro::trafo(r[j - 1], theta[i + 1], x, y, 0);
                        r_tmp     = r[j - 1];
                        theta_tmp = theta[i + 1];
                    }
                    else {
                        double pi2 = 2 * PI;
                        // gyro::trafo(r[j - 1], pi2, x, y, 0);
                        r_tmp     = r[j - 1];
                        theta_tmp = pi2;
                    }
                }
                if (i + 1 < ntheta_int) // Define theta(i+1), this might be on the periodicity condition.
                    thetap1 = theta[i + 1];
                else
                    thetap1 = 2 * PI;
                // j=1 means r-h is on the boundary, gyro::distBoundary([r(j),0])==0 means that that Dirichlet BC are set at (r(j),0)) --> to symmetrize, put it on the right hand side
                if (j == 1 && r[j - 1] > 0 && fabs(gyro::distBoundary(r_tmp, theta_tmp, 0)) < tol_bound_check) {
                }
                else {
                    // second circle (at least one 'real' circle with ntheta_int nodes which is closer to the origin) but no dirichlet bc on innermost circle
                    if (j > 1 || (j == 1 && r[0] > 0)) {
                        row_indices.push_back(row_index);
                        if (i + 1 < ntheta_int) // next node in theta direction but one circle before
                            col_indices.push_back(row_index - ntheta_int + 1);
                        else { // periodicity condition
                            col_indices.push_back(row_index - 2 * ntheta_int + 1);
                        }
                    }

                    // alt: neumBC
                    if (j > 1 || (j == 1 && r[0] > 0)) {
                        // a l'interieur sans contact aux conditions aux limites
                        //        if j > 1 // pas de reference vers l'origine ici
                        vals.push_back(gyro::art(r[j - 1], theta[i], 0) + gyro::art(r[j], thetap1, 0));
                        index++;
                    }
                }
            }

            // (r,phi)
            row_indices.push_back(row_index);
            col_indices.push_back(row_index);

            if (i + 1 < ntheta_int) // Define theta(i+1), this might be on the periodicity condition.
                thetap1 = theta[i + 1];
            else
                thetap1 = 2 * PI;
            if (j > 0) {
                // alt: neumBound

                vals.push_back(0.5 * kt / hs * (gyro::arr(r[j], theta[i], 0) + gyro::arr(r[j + 1], theta[i], 0)) +
                               0.5 * hs / kt * (gyro::att(r[j], theta[i], 0) + gyro::att(r[j], thetap1, 0)) +
                               0.5 * kt / hsmin1 * (gyro::arr(r[j - 1], theta[i], 0) + gyro::arr(r[j], theta[i], 0)) +
                               0.5 * hsmin1 / kt * (gyro::att(r[j], theta[i], 0) + gyro::att(r[j], thetap1, 0)) +
                               0.5 * ktmin1 / hsmin1 *
                                   (gyro::arr(r[j - 1], theta[i], 0) + gyro::arr(r[j], theta[i], 0)) +
                               0.5 * hsmin1 / ktmin1 * (gyro::att(r[j], thetamin1, 0) + gyro::att(r[j], theta[i], 0)) +
                               0.5 * ktmin1 / hs * (gyro::arr(r[j], theta[i], 0) + gyro::arr(r[j + 1], theta[i], 0)) +
                               0.5 * hs / ktmin1 * (gyro::att(r[j], thetamin1, 0) + gyro::att(r[j], theta[i], 0)));
            }
            else { //across the origin; j=0

                vals.push_back(0.5 * kt / hs * (gyro::arr(r[j], theta[i], 0) + gyro::arr(r[j + 1], theta[i], 0)) +
                               0.5 * hs / kt * (gyro::att(r[j], theta[i], 0) + gyro::att(r[j], thetap1, 0)) +
                               0.5 * hsmin1 / kt * (gyro::att(r[j], theta[i], 0) + gyro::att(r[j], thetap1, 0)) +
                               0.5 * hsmin1 / ktmin1 * (gyro::att(r[j], thetamin1, 0) + gyro::att(r[j], theta[i], 0)) +
                               0.5 * ktmin1 / hs * (gyro::arr(r[j], theta[i], 0) + gyro::arr(r[j + 1], theta[i], 0)) +
                               0.5 * hs / ktmin1 * (gyro::att(r[j], thetamin1, 0) + gyro::att(r[j], theta[i], 0)));
                if (i + 1 > ntheta_int / 2)
                    vals[index] +=
                        0.5 * kt / hsmin1 * (gyro::arr(r[j], theta[i] + PI, 0) + gyro::arr(r[j], theta[i], 0)) +
                        0.5 * ktmin1 / hsmin1 * (gyro::arr(r[j], theta[i] + PI, 0) + gyro::arr(r[j], theta[i], 0));
                else // to make it symmetric, take values from second half circle
                    vals[index] +=
                        0.5 * kt / hsmin1 * (gyro::arr(r[j], theta[i], 0) + gyro::arr(r[j], theta[i] + PI, 0)) +
                        0.5 * ktmin1 / hsmin1 * (gyro::arr(r[j], theta[i], 0) + gyro::arr(r[j], theta[i] + PI, 0));
            }

            vals[index] += betaVec[row_index];
            index++;
        }
    }

    int size_ri          = row_indices.size();
    int realsize_indices = 0;
    for (int i = 0; i < size_ri; i++)
        if (row_indices[i] > realsize_indices)
            realsize_indices = row_indices[i];

    // Dirichlet-RB (outside; last node in slice numbering)
    double t_tmp = 0;
    // gyro::trafo(r[nr_int], t_tmp, x, y, 0);
    if (fabs(gyro::distBoundary(r[nr_int], t_tmp, 0)) <
        tol_bound_check) { // dirBound([r(n+1),0])==0 means that (r(n+1),0) is on Dirichlet boundary
        // if (bound_def_rad && r[nr_int] != Rmax)
        //     throw std::runtime_error("Node on the boundary erroneously detected as node with boundary conditions "
        //                              "(although there are none...)");
        // alt: r0=0 ?
        double ri_end = row_indices[size_ri - 1];
        for (int i = 0; i < ntheta_int; i++) {
            row_indices.push_back(ri_end + i + 1);
            col_indices.push_back(ri_end + i + 1); // Matlab format....
            vals.push_back(1.0);
        }
    }
} /* ----- end of level::build_A0 ----- */

/*!
 *  \brief Build the RHS (deprecated)
 *
 * Builds the RHS of the system based on 9p FD. Relations with Dirichlet boundary 
 * condition nodes are shifted to the RHS to keep a symmetric operator.
 *
 */
void level::build_rhs0()
{
    // Take boundary condition into account

    // created transformed vector from original right side
    // zu -Laplace(u) = f to K(0,1) after transformed to polar coordinates, then
    // right side for r in [eps, 1-eps]: r*f(r,phi) (values in loop below)
    // D_scal_rhs = ones(global_size,1); % save h and k dependent scaling of the right side
    // alt: r0=0/neumBound
    int start_j;
    // dirBound([r(1),0])==0 means that (r(1),0) is on Dirichlet boundary
    if (gyro::icntl[Param::DirBC_Interior]) {
        start_j = 1;
    }
    // dirBound([r(1),0])>0 means that (r(1),0) is not on Dirichlet boundary
    else if (r[0] > 0 && !gyro::icntl[Param::DirBC_Interior]) {
        start_j = 0;
    }
    // dirBound([r(1),0])==0 means that (r(1),0) is on Dirichlet boundary
    if (gyro::icntl[Param::DirBC_Interior]) {
        for (int i = 0; i < ntheta_int; i++) {
            // Dirichlet-RB inside % one value for each point inside (Dirichlet-RB), i.e. ntheta_int-many
            fVec.push_back(gyro::def_solution_rt(r[0], theta[i], 0));
        }
    }

    // int i1 = 0, i2 = 0;
    for (int j = start_j; j < nr_int; j++) { // nodes of Disk K(eps,1-eps)
        for (int i = 0; i < ntheta_int; i++) {
            // alt: r0=0
            int row_index = j * ntheta_int + i;

            // % ==================== %
            double kt = thetaplus[i];
            double thetamin1, ktmin1, hs, hsmin1;
            if (i > 0) {
                thetamin1 = theta[i - 1];
                ktmin1    = thetaplus[i - 1];
            }
            else {
                thetamin1 = theta[ntheta_int - 1]; // equal to 2*PI-ktmin1
                ktmin1    = thetaplus[ntheta_int - 1];
            }
            // % -------------------- %
            // r(j+1) is the current r
            hs = hplus[j];
            if (j > 0)
                hsmin1 = hplus[j - 1];
            else
                hsmin1 = 2 * r[0]; // across the origin
            // % ==================== %

            // right side (r=r(j+1) because r(1)=0 or start of j at 0 if r(1)>0)
            //             fac_hk=1; %% Factor now in diagonal matrix
            double fac_hk = 0.25 * (kt + ktmin1) * (hs + hsmin1);
            fVec.push_back(fac_hk * fabs(gyro::detDFinv(r[j], theta[i], 0)) *
                           gyro::eval_def_rhs(r[j], theta[i], 0)); // right side for r\in(eps,1-eps): r*f(r,phi)

            // 9-Point Stencil (attention order; first bottom, then right, then top right, top left, bottom left, middle)
            // r-h- (h- = x_i-x_i-1)

            // j=1 means r-h is on the boundary, dirBound([r(j),0])==0 means that that Dirichlet BC are set at (r(j),0))
            // --> to symmetrize, put it on the right hand side
            if (j == 1 && r[j - 1] > 0 && gyro::icntl[Param::DirBC_Interior]) {
                // conditions aux limites Dirichlet
                fVec[row_index] +=
                    (0.5 * kt / hsmin1 * (gyro::arr(r[j - 1], theta[i], 0) + gyro::arr(r[j], theta[i], 0)) +
                     0.5 * ktmin1 / hsmin1 * (gyro::arr(r[j - 1], theta[i], 0) + gyro::arr(r[j], theta[i], 0))) *
                    gyro::def_solution_rt(r[j - 1], theta[i], 0);
            }

            // (r-h-,phi-k) (en bas a gauche)
            if (gyro::icntl[Param::mod_pk] > 0) {
                // j=1 means r-h is on the boundary, dirBound([r(j),0])==0 means that that Dirichlet BC are set at (r(j),0))
                // --> to symmetrize, put it on the right hand side
                if (j == 1 && r[j - 1] > 0 && gyro::icntl[Param::DirBC_Interior]) {
                    // conditions aux limites Dirichlet
                    fVec[row_index] +=
                        (gyro::art(r[j], thetamin1, 0) + gyro::art(r[j - 1], theta[i], 0)) *
                        gyro::def_solution_rt(r[j - 1], thetamin1,
                                              0); // factor is positive on left hand side, so it is negative here...
                }
            }

            // phi-k

            // (r+h+,phi-k) (en haut a droite)
            double theta_eval_prelast = 0;
            if (gyro::icntl[Param::mod_pk] > 0) {
                if (i <= 0) {
                    theta_eval_prelast = theta[ntheta_int - 1];
                }

                // j=nr_int-1 means r+h is on the boundary --> to symmetrize, put it on the right hand side
                if (j == nr_int - 1) {
                    if (i > 0) { // not the first node on the corresponding circle; we can take theta(i-1)
                        fVec[row_index] -= (gyro::art(r[j], thetamin1, 0) + gyro::art(r[j + 1], theta[i], 0)) *
                                           gyro::def_solution_rt(r[j + 1], theta[i - 1], 0);
                    }
                    else {
                        fVec[row_index] -= (gyro::art(r[j], thetamin1, 0) + gyro::art(r[j + 1], theta[i], 0)) *
                                           gyro::def_solution_rt(r[j + 1], theta_eval_prelast, 0);
                    }
                }
            }

            // r + h+ (h+ = x_i+1-x_i)
            // j=nr_int-1 means r+h is on the boundary --> to symmetrize, put it on the right hand side
            if (j == nr_int - 1) {
                fVec[row_index] +=
                    (0.5 * kt / hs * (gyro::arr(r[j], theta[i], 0) + gyro::arr(r[j + 1], theta[i], 0)) +
                     0.5 * ktmin1 / hs * (gyro::arr(r[j], theta[i], 0) + gyro::arr(r[j + 1], theta[i], 0))) *
                    gyro::def_solution_rt(r[j + 1], theta[i], 0);
            }

            ////// INSERT HERE IN THE UPPER LEFT CORNER
            // (r+h+,phi+k) (en haut a gauche)

            double thetap1 = 0;
            if (gyro::icntl[Param::mod_pk] > 0) {
                if (i + 1 < ntheta_int) // Define theta(i+1), this might be on the periodicity condition.
                    thetap1 = theta[i + 1];
                else
                    thetap1 = 2 * PI;

                // j=nr_int-1 means r+h is on the boundary --> to symmetrize, put it on the right hand side
                if (j == nr_int - 1) {
                    fVec[row_index] += (gyro::art(r[j + 1], theta[i], 0) + gyro::art(r[j], thetap1, 0)) *
                                       gyro::def_solution_rt(r[j + 1], thetap1, 0);
                }
            }

            // phi+k
            double r_tmp, theta_tmp;

            // (r-h-,phi+k) (en bas a gauche)
            if (gyro::icntl[Param::mod_pk] > 0) {
                if (j == 1) { // coords necessary for potential boundary conditions
                    if (i + 1 < ntheta_int) {
                        // gyro::trafo(r[j - 1], theta[i + 1], x, y, 0);
                        r_tmp     = r[j - 1];
                        theta_tmp = theta[i + 1];
                    }
                    else {
                        double pi2 = 2 * PI;
                        // gyro::trafo(r[j - 1], pi2, x, y, 0);
                        r_tmp     = r[j - 1];
                        theta_tmp = pi2;
                    }
                }
                if (i + 1 < ntheta_int) // Define theta(i+1), this might be on the periodicity condition.
                    thetap1 = theta[i + 1];
                else
                    thetap1 = 2 * PI;
                // j=1 means r-h is on the boundary, gyro::distBoundary([r(j),0])==0 means that that Dirichlet BC are set at (r(j),0)) --> to symmetrize, put it on the right hand side
                if (j == 1 && r[j - 1] > 0 && gyro::icntl[Param::DirBC_Interior]) {
                    // conditions aux limites Dirichlet
                    // factor is positive on left hand side, so it is negative here...
                    fVec[row_index] -= (gyro::art(r[j - 1], theta[i], 0) + gyro::art(r[j], thetap1, 0)) *
                                       gyro::def_solution_rt(r_tmp, theta_tmp, 0);
                }
            }

            // (r,phi)
        }
    }

    for (int i = 0; i < ntheta_int; i++) {
        // Dirichlet-RB outside % one value for each point outside (Dirichlet-RB), i.e. ntheta-many,
        fVec.push_back(gyro::def_solution_rt(r[nr_int], theta[i], 0));
    }
} /* ----- end of level::build_rhs0 ----- */

/*!
 *  \brief Applies the operator A without construction (deprecated)
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
 */
void level::apply_A0(std::vector<double> u, std::vector<double>& Au)
{
    // if boundary is only defined on radius, checks for wrongly detected boundary
    // nodes due to rounding errors (double checking...)
    double tol_bound_check = gyro::dcntl[Param::tol_bound_check];

    // double x, y;
    int index;

    // Take boundary condition into account
    // In the operator
    // gyro::trafo(r[0], theta[0], x, y, 0);
    // alt: r0=0/neumBound
    // dirBound([r(1),0])==0 means that (r(1),0) is on Dirichlet boundary
    if (fabs(gyro::distBoundary(r[0], theta[0], 0)) < tol_bound_check) {
        // Dirichlet-RB (intérieur ; premier noeud dans la numérotation des disques)
        for (int i = 0; i < ntheta_int; i++) {
            Au[i] += u[i];
        }
        index = ntheta_int;
    }
    else {
        index = 0;
    }

    int row, col;
    double v;

    // gyro::trafo(r[0], theta[0], x, y, 0);
    // alt: r0=0/neumBound
    int start_j = 0;
    // dirBound([r(1),0])==0 means that (r(1),0) is on Dirichlet boundary
    if (fabs(gyro::distBoundary(r[0], theta[0], 0)) < tol_bound_check) {
        // alt: neumBC
        start_j = 1;
    }
    // dirBound([r(1),0])>0 means that (r(1),0) is not on Dirichlet boundary
    else if (r[0] > 0 && fabs(gyro::distBoundary(r[0], theta[0], 0)) > tol_bound_check) {
        // alt: neumBC
        start_j = 0;
    }
    // int i1 = 0, i2 = 0;
    for (int j = start_j; j < nr_int; j++) { // nodes of Disk K(eps,1-eps)
        for (int i = 0; i < ntheta_int; i++) {
            // alt: r0=0
            int row_index = j * ntheta_int + i;

            // % ==================== %
            double kt = thetaplus[i];
            double thetamin1, ktmin1, hs, hsmin1;
            if (i > 0) {
                thetamin1 = theta[i - 1];
                ktmin1    = thetaplus[i - 1];
            }
            else {
                thetamin1 = theta[ntheta_int - 1]; // equal to 2*PI-ktmin1
                ktmin1    = thetaplus[ntheta_int - 1];
            }
            // % -------------------- %
            // r(j+1) is the current r
            hs = hplus[j];
            if (j > 0)
                hsmin1 = hplus[j - 1];
            else
                hsmin1 = 2 * r[0]; // across the origin
            // % ==================== %

            // 9-Point Stencil (attention order; first bottom, then right, then top right, top left, bottom left, middle)
            // r-h- (h- = x_i-x_i-1)
            // if (j == 1)
            //     gyro::trafo(r[j - 1], theta[i], x, y, 0);
            // j=1 means r-h is on the boundary, dirBound([r(j),0])==0 means that that Dirichlet BC are set at (r(j),0))
            // --> to symmetrize, put it on the right hand side
            if (j == 1 && r[j - 1] > 0 && fabs(gyro::distBoundary(r[j - 1], theta[i], 0)) < tol_bound_check) {
            }
            else {
                // row_indices[index] = row_index;
                row = row_index;
                col = -42;
                v   = -42;
                // second circle (at least one 'real' circle with ntheta_int nodes which is closer to the origin)
                if (j > 1 || (j == 1 && r[0] > 0)) {
                    // to u at r-h geometrically recorded (i.e. -ntheta_int nodes at slice sorting)
                    col = row_index - ntheta_int;
                }
                else { // j=0 or j=1
                    if (r[0] == 0) { // all refering to the single origin node! (j=0 not possible in this case)
                        // Reference to u at r-h=0 geometrically recorded (for all i the origin, because here r=h)
                        col = 1;
                    }
                    else if (j == 0) { // across the origin
                        if (i + 1 > ntheta_int / 2) {
                            col = row_index - ntheta_int / 2; // half a turn back
                        }
                        else {
                            col = row_index + ntheta_int / 2; // half a turn further
                        }
                    }
                }

                if (j > 0) {
                    // alt: neumBC
                    // a l interieur sans contact aux conditions aux limites
                    v = -0.5 * kt / hsmin1 * (gyro::arr(r[j - 1], theta[i], 0) + gyro::arr(r[j], theta[i], 0)) -
                        0.5 * ktmin1 / hsmin1 * (gyro::arr(r[j - 1], theta[i], 0) + gyro::arr(r[j], theta[i], 0));
                }
                else {
                    // j==0; across the origin
                    // just use r_{s-1}=0
                    if (i + 1 > ntheta_int / 2)
                        // Use r(0):=r(j) to go to the other side of the origin
                        v = -0.5 * kt / hsmin1 * (gyro::arr(r[j], theta[i] + PI, 0) + gyro::arr(r[j], theta[i], 0)) -
                            0.5 * ktmin1 / hsmin1 * (gyro::arr(r[j], theta[i] + PI, 0) + gyro::arr(r[j], theta[i], 0));
                    else // to make it symmetric, take values from second half circle
                        v = -0.5 * kt / hsmin1 * (gyro::arr(r[j], theta[i], 0) + gyro::arr(r[j], theta[i] + PI, 0)) -
                            0.5 * ktmin1 / hsmin1 * (gyro::arr(r[j], theta[i], 0) + gyro::arr(r[j], theta[i] + PI, 0));
                }
                Au[row] += v * u[col];
                index++;
            }
            // (r-h-,phi-k) (en bas a gauche)
            if (gyro::icntl[Param::mod_pk] > 0) {
                // if (j == 1) // coords necessary for potential boundary conditions
                //     gyro::trafo(r[j - 1], thetamin1, x, y, 0); // r(j) is PREVIOUS r, actual is r(J+1) !
                // j=1 means r-h is on the boundary, dirBound([r(j),0])==0 means that that Dirichlet BC are set at (r(j),0))
                // --> to symmetrize, put it on the right hand side
                if (j == 1 && r[j - 1] > 0 && fabs(gyro::distBoundary(r[j - 1], thetamin1, 0)) < tol_bound_check) {
                }
                else {
                    // second circle (at least one 'real' circle with ntheta_int nodes which is closer to the origin)
                    // but no dirichlet bc on innermost circle
                    if (j > 1 || (j == 1 && r[0] > 0)) {
                        row = row_index;
                        if (i > 0) // next node in theta direction but one circle before
                            col = row_index - ntheta_int - 1;
                        else // periodicity condition
                            col = row_index - 1;
                    }

                    //    if j > 0
                    // alt: neumBound
                    if (j > 1 || (j == 1 && r[0] > 0)) {
                        // a l interieur sans contact aux conditions aux limites
                        //         if (j > 1) // pas de reference vers l origine ici
                        v = -(gyro::art(r[j], thetamin1, 0) + gyro::art(r[j - 1], theta[i], 0));

                        Au[row] += v * u[col];
                        index++;
                    }
                }
            }

            // phi-k
            // row_indices[index] = row_index;
            row = row_index;
            if (i > 0) // previous node in phi direction
                col = row_index - 1;
            else // periodicity condition
                col = row_index + ntheta_int - 1;
            if (j > 0)
                v = -0.5 * hsmin1 / ktmin1 * (gyro::att(r[j], thetamin1, 0) + gyro::att(r[j], theta[i], 0)) -
                    0.5 * hs / ktmin1 * (gyro::att(r[j], thetamin1, 0) + gyro::att(r[j], theta[i], 0));

            else // j=0; take r(0)=0 //
                v = -0.5 * hsmin1 / ktmin1 * (gyro::att(r[j], thetamin1, 0) + gyro::att(r[j], theta[i], 0)) -
                    0.5 * hs / ktmin1 * (gyro::att(r[j], thetamin1, 0) + gyro::att(r[j], theta[i], 0));

            Au[row] += v * u[col];
            index++;

            // (r+h+,phi-k) (en haut a droite)
            double theta_eval_prelast = 0;
            double r_tmp, theta_tmp;
            if (gyro::icntl[Param::mod_pk] > 0) {
                if (i > 0) { // not the first node on the corresponding circle; we can take theta(i-1)
                    // gyro::trafo(r[j + 1], theta[i - 1], x, y, 0);
                    r_tmp     = r[j + 1];
                    theta_tmp = theta[i - 1];
                }
                else {
                    theta_eval_prelast = theta[ntheta_int - 1];
                    // gyro::trafo(r[j + 1], theta_eval_prelast, x, y, 0);
                    r_tmp     = r[j + 1];
                    theta_tmp = theta_eval_prelast;
                }

                // means that r[j+1] is not on the Dirichlet boundary
                if (j < nr_int - 1 ||
                    (j == nr_int - 1 && fabs(gyro::distBoundary(r_tmp, theta_tmp, 0)) > tol_bound_check)) {
                    // row_indices[index] = row_index;
                    row = row_index;
                    if (i > 0) // previous node in phi direction but at r+h
                        col = row_index + ntheta_int - 1;
                    else // first node on corresponding circle; pay gyro::attention to periodicity condition!
                        col = row_index + 2 * ntheta_int - 1;
                    v = gyro::art(r[j], thetamin1, 0) + gyro::art(r[j + 1], theta[i], 0);
                    Au[row] += v * u[col];
                    index++;
                }
            }

            // r + h+ (h+ = x_i+1-x_i)
            // gyro::trafo(r[j + 1], theta[i], x, y, 0);
            // means that r[j+1] is not on the Dirichlet boundary
            if (j < nr_int - 1 ||
                (j == nr_int - 1 && fabs(gyro::distBoundary(r[j + 1], theta[i], 0)) > tol_bound_check)) {
                row = row_index;
                col = row_index + ntheta_int; // zu u bei r+h
                v   = -0.5 * kt / hs * (gyro::arr(r[j], theta[i], 0) + gyro::arr(r[j + 1], theta[i], 0)) -
                    0.5 * ktmin1 / hs * (gyro::arr(r[j], theta[i], 0) + gyro::arr(r[j + 1], theta[i], 0));
                Au[row] += v * u[col];
                index++;
            }
            // j=nr_int-1 means r+h is on the boundary --> to symmetrize, put it on the right hand side
            else if (j == nr_int - 1 && fabs(gyro::distBoundary(r[j + 1], theta[i], 0)) < tol_bound_check) {
            }

            // (r+h+,phi+k) (en haut a gauche)
            double thetap1 = 0;
            if (gyro::icntl[Param::mod_pk] > 0) {
                if (i + 1 < ntheta_int) // Define theta(i+1), this might be on the periodicity condition.
                    thetap1 = theta[i + 1];
                else
                    thetap1 = 2 * PI;

                // means that r[j+1] is not on the Dirichlet boundary
                if (j < nr_int - 1 ||
                    (j == nr_int - 1 && fabs(gyro::distBoundary(r[j + 1], theta[i], 0)) > tol_bound_check)) {
                    row = row_index;
                    if (i + 1 < ntheta_int) // previous node in phi direction but at r+h
                        col = row_index + ntheta_int + 1;
                    else // first node on corresponding circle; pay gyro::attention to periodicity condition!
                        col = row_index + 1;
                    v = -gyro::art(r[j + 1], theta[i], 0) - gyro::art(r[j], thetap1, 0);

                    Au[row] += v * u[col];
                    index++;
                }
                // j=nr_int-1 means r+h is on the boundary --> to symmetrize, put it on the right hand side
                else if (j == nr_int - 1 && fabs(gyro::distBoundary(r[j + 1], theta[i], 0)) < tol_bound_check) {
                    // gyro::trafo(r[j + 1], thetap1, x, y, 0);
                }
            }

            // phi+k
            // row_indices[index] = row_index;
            row = row_index;
            if (i + 1 < ntheta_int) // next node in theta direction
                col = row_index + 1;
            else // periodicity condition
                col = row_index - ntheta_int + 1;
            if (i + 1 < ntheta_int) // Define theta(i+1), this might be on the periodicity condition.
                thetap1 = theta[i + 1];
            else
                thetap1 = 2 * PI;
            if (j > 0)
                v = -0.5 * hs / kt * (gyro::att(r[j], theta[i], 0) + gyro::att(r[j], thetap1, 0)) -
                    0.5 * hsmin1 / kt * (gyro::att(r[j], theta[i], 0) + gyro::att(r[j], thetap1, 0));
            else // j=0;
                v = -0.5 * hs / kt * (gyro::att(r[j], theta[i], 0) + gyro::att(r[j], thetap1, 0)) -
                    0.5 * hsmin1 / kt * (gyro::att(r[j], theta[i], 0) + gyro::att(r[j], thetap1, 0));

            Au[row] += v * u[col];
            index++;

            // (r-h-,phi+k) (en bas a gauche)
            if (gyro::icntl[Param::mod_pk] > 0) {
                double theta_tmp;
                if (j == 1) { // coords necessary for potential boundary conditions
                    if (i + 1 < ntheta_int) {
                        // gyro::trafo(r[j - 1], theta[i + 1], x, y, 0);
                        r_tmp     = r[j - 1];
                        theta_tmp = theta[i + 1];
                    }
                    else {
                        double pi2 = 2 * PI;
                        // gyro::trafo(r[j - 1], pi2, x, y, 0);
                        r_tmp     = r[j - 1];
                        theta_tmp = pi2;
                    }
                }
                if (i + 1 < ntheta_int) // Define theta(i+1), this might be on the periodicity condition.
                    thetap1 = theta[i + 1];
                else
                    thetap1 = 2 * PI;
                // j=1 means r-h is on the boundary, gyro::distBoundary([r(j),0])==0 means that that Dirichlet BC are set at (r(j),0)) --> to symmetrize, put it on the right hand side
                if (j == 1 && r[j - 1] > 0 && fabs(gyro::distBoundary(r_tmp, theta_tmp, 0)) < tol_bound_check) {
                }
                else {
                    // second circle (at least one 'real' circle with ntheta_int nodes which is closer to the origin) but no dirichlet bc on innermost circle
                    if (j > 1 || (j == 1 && r[0] > 0)) {
                        row = row_index;
                        if (i + 1 < ntheta_int) // next node in theta direction but one circle before
                            col = row_index - ntheta_int + 1;
                        else { // periodicity condition
                            col = row_index - 2 * ntheta_int + 1;
                            // theta_eval_prelast = theta[ntheta_int];
                        }
                    }

                    // alt: neumBC
                    if (j > 1 || (j == 1 && r[0] > 0)) {
                        // a l'interieur sans contact aux conditions aux limites
                        //        if j > 1 // pas de reference vers l'origine ici
                        v = gyro::art(r[j - 1], theta[i], 0) + gyro::art(r[j], thetap1, 0);
                        Au[row] += v * u[col];
                        index++;
                    }
                }
            }

            // (r,phi)
            row = row_index;
            col = row_index;

            if (i + 1 < ntheta_int) // Define theta(i+1), this might be on the periodicity condition.
                thetap1 = theta[i + 1];
            else
                thetap1 = 2 * PI;
            if (j > 0)
                // alt: neumBound
                v = 0.5 * kt / hs * (gyro::arr(r[j], theta[i], 0) + gyro::arr(r[j + 1], theta[i], 0)) +
                    0.5 * hs / kt * (gyro::att(r[j], theta[i], 0) + gyro::att(r[j], thetap1, 0)) +
                    0.5 * kt / hsmin1 * (gyro::arr(r[j - 1], theta[i], 0) + gyro::arr(r[j], theta[i], 0)) +
                    0.5 * hsmin1 / kt * (gyro::att(r[j], theta[i], 0) + gyro::att(r[j], thetap1, 0)) +
                    0.5 * ktmin1 / hsmin1 * (gyro::arr(r[j - 1], theta[i], 0) + gyro::arr(r[j], theta[i], 0)) +
                    0.5 * hsmin1 / ktmin1 * (gyro::att(r[j], thetamin1, 0) + gyro::att(r[j], theta[i], 0)) +
                    0.5 * ktmin1 / hs * (gyro::arr(r[j], theta[i], 0) + gyro::arr(r[j + 1], theta[i], 0)) +
                    0.5 * hs / ktmin1 * (gyro::att(r[j], thetamin1, 0) + gyro::att(r[j], theta[i], 0));
            else { //across the origin; j=0
                v = 0.5 * kt / hs * (gyro::arr(r[j], theta[i], 0) + gyro::arr(r[j + 1], theta[i], 0)) +
                    0.5 * hs / kt * (gyro::att(r[j], theta[i], 0) + gyro::att(r[j], thetap1, 0)) +
                    0.5 * hsmin1 / kt * (gyro::att(r[j], theta[i], 0) + gyro::att(r[j], thetap1, 0)) +
                    0.5 * hsmin1 / ktmin1 * (gyro::att(r[j], thetamin1, 0) + gyro::att(r[j], theta[i], 0)) +
                    0.5 * ktmin1 / hs * (gyro::arr(r[j], theta[i], 0) + gyro::arr(r[j + 1], theta[i], 0)) +
                    0.5 * hs / ktmin1 * (gyro::att(r[j], thetamin1, 0) + gyro::att(r[j], theta[i], 0));
                if (i + 1 > ntheta_int / 2)
                    v += 0.5 * kt / hsmin1 * (gyro::arr(r[j], theta[i] + PI, 0) + gyro::arr(r[j], theta[i], 0)) +
                         0.5 * ktmin1 / hsmin1 * (gyro::arr(r[j], theta[i] + PI, 0) + gyro::arr(r[j], theta[i], 0));
                else // to make it symmetric, take values from second half circle
                    v += 0.5 * kt / hsmin1 * (gyro::arr(r[j], theta[i], 0) + gyro::arr(r[j], theta[i] + PI, 0)) +
                         0.5 * ktmin1 / hsmin1 * (gyro::arr(r[j], theta[i], 0) + gyro::arr(r[j], theta[i] + PI, 0));
            }
            Au[row] += v * u[col];

            double val = betaVec[row];
            Au[row] += val * u[col];
            index++;
        }
    }

    // Dirichlet-RB (outside; last node in slice numbering)
    double t_tmp = 0;
    // gyro::trafo(r[nr_int], t_tmp, x, y, 0);
    if (fabs(gyro::distBoundary(r[nr_int], t_tmp, 0)) <
        tol_bound_check) { // dirBound([r(n+1),0])==0 means that (r(n+1),0) is on Dirichlet boundary
        // alt: r0=0 ?
        for (int i = 0; i < ntheta_int; i++) {
            Au[row + i + 1] += u[row + i + 1];
        }
    }
} /* ----- end of level::apply_A0 ----- */
