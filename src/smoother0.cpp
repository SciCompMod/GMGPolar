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
 * \file multigrid_iter.cpp
 * \brief Implementation of the smoother (deprecated)
 * \author M. Kuehn, C. Kruse, P. Leleux, C. Schwarz
 * \version 0.0
 */

#include "level.h"

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
void level::multigrid_smoothing0(int smoother)
{
    //std::cout << "\n---------------- smoother: " << smoother << std::endl;

    double t, t_smoothing_tmp;
    TIC;
    t_smoothing_tmp = t;

    //! solving for u_sc (solution vector)
    std::vector<double> u_sc;

    //create f_sc
    std::vector<double> f_sc;
    build_fsc0(f_sc, smoother);

    t_f_sc += TOC;
    TIC;

    //compute f_total = f_sc - A_sc_ortho * u
    //size of u_sc for the definition of size of A_sc_ortho = size of f_sc
    std::vector<double> f_total(m_sc[smoother], 0);
    apply_Asc_ortho0(f_total, smoother);

    for (long unsigned int i = 0; i < f_total.size(); ++i) {
        f_total[i] = f_sc[i] - f_total[i];
    }

    t_Asc_ortho += TOC;
    TIC;

#ifdef GMGPOLAR_USE_MUMPS
    if (gyro::icntl[Param::optimized] == 0) {
#endif
        u_sc =
            solve_gaussian_elimination(A_Zebra_r_LU[smoother], A_Zebra_c_LU[smoother], A_Zebra_v_LU[smoother], f_total);
#ifdef GMGPOLAR_USE_MUMPS
    }
    else
        u_sc = solve_mumps(mumps_A_Zebra[smoother], f_total);
#endif

    if (gyro::icntl[Param::verbose] > 5)
        gyro::disp(u_sc, "u_sc");

    t_Asc += TOC;
    TIC;

    //-------------------------------------------------------------------------------------------
    //! insert the solution vector u_sc into the whole vector u, depending on the smoother/colour
    int extrapol = gyro::icntl[Param::extrapolation];
    //only for the circle black smoother, don't change the variable of the class itself
    int ntheta_int_local = ntheta_int;
    if (extrapol == 1 && l == 0 && smoother == 0) { //circle, black
        ntheta_int_local = ntheta_int / 2;
    }
    int n_lines_radial_b = ceil((double)ntheta_int / 2);
    int n_lines_radial_w = floor((double)ntheta_int / 2);

    //computation of indices in the total vector u corresponding to the indices in u_sc
    for (long unsigned int i = 0; i < u_sc.size(); ++i) {
        int row;
        int col;

        if (smoother < 2) { //circle
            row = i / ntheta_int_local; //row within the smoother
            row = row * 2; //row within the whole matrix
            col = i % ntheta_int_local;
            if (smoother == 1) { //white
                row++;
            }
            if (extrapol == 1 && l == 0 && smoother == 0) { //black
                col = col * 2 + 1; //augment col in case of extrapolation
            }
        }
        else { //radial
            if (smoother == 2) { //black
                row = i / n_lines_radial_b; //row within the smoother
                col = i % n_lines_radial_b; //col within the smoother
                col = col * 2; //col within the whole matrix
                if (extrapol == 1 && l == 0) {
                    row = row * 2; //augment row in case of extrapolation
                    if (delete_circles % 2 == 0) { //delete_circles = even
                        row = row + 1;
                    }
                }
            }
            else { //white
                row = i / n_lines_radial_w; //row within the smoother
                col = i % n_lines_radial_w; //col within the smoother
                col = col * 2 + 1; //col within the whole matrix
            }
            row += delete_circles; //row within the whole matrix
        }

        int index = row * ntheta_int + col;

        if (gyro::icntl[Param::smoother] == 13) {
            if (smoother < 2) { //circle
                u_previous_c[index] = u_sc[i];
            }
            else { //radial
                u_previous_r[index] = u_sc[i];
            }
        }
        else {
            u[index] = u_sc[i];
        }
    }

    t_f_sc += TOC;
    TIC;

    t = t_smoothing_tmp;
    t_smoothing += TOC;
    //std::cout << "smoothing end \n";
}

/*! \brief Create the RHS part corresponding to Asc_ortho for a smoother (on level l)
 *
 * Create f_sc, i.e. the RHS part corresponding to Asc_ortho for a smoother (on level l)
 * 
 * \param f_sc: the RHS part (out)
 * \param smoother: the smoother
*/
void level::build_fsc0(std::vector<double>& f_sc, int smoother)
{
    int n_indices_circle = delete_circles * ntheta_int; //delete_circles = index of radius of border between smoothers
    int extrapol         = gyro::icntl[Param::extrapolation];

    if (smoother == 0) { //circle black smoother
        if (extrapol == 1 && l == 0) { //extrapolation
            //skip every second element
            for (int ind = 1; ind < n_indices_circle; ind = ind + 2) { //iteration over all elements in the smoother
                int r_index = ind / ntheta_int; //index in r-direction
                //check if even r_index
                if (!(r_index % 2)) {
                    f_sc.push_back(fVec[ind]); //insert elements of f corresponding to the colour
                }
            }
        }
        else { //no extrapolation
            for (int ind = 0; ind < n_indices_circle; ++ind) { //iteration over all elements in the smoother
                int r_index = ind / ntheta_int; //index in r-direction
                //check if even r_index
                if (!(r_index % 2)) {
                    f_sc.push_back(fVec[ind]); //insert elements of f corresponding to the colour
                }
            }
        }
    }
    else if (smoother == 1) { //circle white smoother
        for (int ind = 0; ind < n_indices_circle; ++ind) { //iteration over all elements in the smoother
            int ind_local = ind;
            int r_index   = ind_local / ntheta_int; //index in r-direction
            //check odd r_index
            if (r_index % 2) {
                f_sc.push_back(fVec[ind]); //insert elements of f corresponding to the colour
            }
        }
    }
    else if (smoother == 2) { //radial black smoother
        //skip every second row
        for (int ind = n_indices_circle; ind < m; ++ind) { //iteration over all elements in the smoother
            int r_index     = ind / ntheta_int; //index in r-direction
            int theta_index = ind - r_index * ntheta_int; //index in theta-direction

            //check if black (even theta_index) or white (odd theta_index)
            if (extrapol == 1 && l == 0) { //extrapolation
                if (!(theta_index % 2) && r_index % 2) { //radial black, even theta_index, uneven row_index
                    f_sc.push_back(fVec[ind]); //insert elements of f corresponding to the colour
                }
            }
            else { //no extrapolation
                if (!(theta_index % 2)) { //radial black, even
                    f_sc.push_back(fVec[ind]); //insert elements of f corresponding to the colour
                }
            }
        }
    }
    else { //radial white smoother
        for (int ind = n_indices_circle; ind < m; ++ind) { //iteration over all elements in the smoother
            int r_index     = ind / ntheta_int; //index in r-direction
            int theta_index = ind - r_index * ntheta_int; //index in theta-direction
            if (theta_index % 2) { //radial white, odd
                f_sc.push_back(fVec[ind]); //insert elements of f corresponding to the colour
            }
        }
    }
} /* ----- end of level::build_subvectors ----- */

/*! \brief Build the matrix A_sc explicitely on the level l (deprecated)
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
void level::build_Asc0()
{
    //smoother=0 (circle), smoother=1 (radial), (0: white, 1: black)

    // if boundary is only defined on radius, checks for wrongly detected boundary
    // nodes due to rounding errors (double checking...)
    double tol_bound_check = 1e-8;
    int smoother           = 0;
    int extrapol           = gyro::icntl[Param::extrapolation];

    A_Zebra_r.assign(4, std::vector<int>());
    A_Zebra_c.assign(4, std::vector<int>());
    A_Zebra_v.assign(4, std::vector<double>());

    int n_cols_radial = ntheta_int / 2;
    int col;
    // double x, y;
    double val;
    int extrapol_fine_node; //indicates if we have extrapolation on level 0 AND a fine node
    // to count the fine nodes, to know the index without the coarse nodes
    std::vector<double> count_nodes = {0, 0, 0, 0}; //cb, cw, rb, rw

    //! Take boundary condition into account in the operator; inner circle is part of the circle smoother (black)
    // gyro::trafo(r[0], theta[0], x, y, 0);
    //check if r0 is on the boundary, if so, we have Diriclet BC
    if (fabs(gyro::distBoundary(r[0], theta[0], 0)) < tol_bound_check) {
        // if (r[0] != gyro::dcntl[Param::r0_DB]) {
        //     throw std::runtime_error("Node on the boundary erroneously detected as node with boundary conditions "
        //                              "(although there are none...)");
        // }
        // Dirichlet-RB (the inner circle, that is on the Diriclet boundary, just insert "ones" on the diagonal)
        for (int i = 0; i < ntheta_int; i++) {
            extrapol_fine_node = 0;
            //only in case of extrapolation on level 0
            if (extrapol == 1 && l == 0) {
                if (coarse_nodes_list_theta[i] == -1 || coarse_nodes_list_r[0] == -1) {
                    //if the index is -1 in both lists, we have a fine node, otherwise the node is coarse
                    extrapol_fine_node = 1;
                }
                else {
                    continue; //the node is coarse, so we can skip this index, continue with next i in for-loop
                }
            }

            //if we have no extrapolation
            //or if we have extrapolation but are not on level 0
            //or if we have extrapolation on level 0, and the node is fine

            count_nodes[smoother]++;
            int index = count_nodes[smoother] - 1;

            A_Zebra_r[smoother].push_back(index); //just take indices of first line (r0), (j=0)
            A_Zebra_c[smoother].push_back(index);

            A_Zebra_v[smoother].push_back(1.0);
        }
    }

    //! check if we start at j=0 (across the origin) or j=1 (Diriclet boundary)
    // gyro::trafo(r[0], theta[0], x, y, 0);
    int start_j;
    if (fabs(gyro::distBoundary(r[0], theta[0], 0)) < tol_bound_check) { //check if r0 is on the boundary
        // if (fabs(gyro::distBoundary(x, y, 0)) < tol_bound_check && r[0] != gyro::dcntl[Param::r0_DB])
        //     throw std::runtime_error("Node on the boundary erroneously detected as node with boundary conditions "
        //                              "(although there are none...)");
        start_j = 1; //Diriclet
    }
    // dirBound([r(1),0])>0 means that (r(1),0) is not on Dirichlet boundary
    else {
        // if (fabs(gyro::distBoundary(x, y, 0)) > tol_bound_check && r[0] == gyro::dcntl[Param::r0_DB]) {
        //     throw std::runtime_error("Node on the boundary erroneously detected as node with boundary conditions "
        //                              "(although there are none...)");
        // }
        start_j = 0; //across the origin
    }

    // % -------------------------------------------------------------------------------------------------- %
    //!Iteration over all points!
    for (int j = start_j; j < nr_int; j++) { // start with second circle (j=1), as we already treated the origin(j=0)
        for (int i = 0; i < ntheta_int; i++) {
            //!extrapolation: check if we have a fine or a coarse node
            extrapol_fine_node = 0;
            if (extrapol == 1 && l == 0) { //only do this in case of extrapolation and level 0
                if (coarse_nodes_list_theta[i] == -1 || coarse_nodes_list_r[j] == -1) {
                    extrapol_fine_node = 1; //the node is fine
                }
                else {
                    continue; //the node is coarse, so we can skip this index, go on with next i in for-loop
                }
            }

            //!check, which smoother the point belongs to: 0 circle black, 1 circle white, 2 radial black, 3 radial white
            if (j < delete_circles) {
                if (j % 2 == 0) {
                    smoother = 0; //circle black (even)
                }
                else {
                    smoother = 1; //circle white (odd)
                }
            }
            else {
                if (i % 2 == 0) {
                    smoother = 2; //radial black (even)
                        //odd number of points in theta-direction, we have one additional black colom
                    if (ntheta_int % 2) {
                        n_cols_radial = ntheta_int / 2 + 1;
                    }
                }
                else {
                    smoother = 3; //radial white (odd)
                }
            }

            //!compute the index of the point in the grid depending on the smoother (the index in the matrix)
            count_nodes[smoother]++;
            int index = count_nodes[smoother] - 1; //index of the point (i,j)

            //!variables for the calculation of vals
            double kt = thetaplus[i]; //k_j
            double thetamin1; //the coordinate of theta to the left
            double thetap1 = 0;
            double ktmin1; //k_j-1
            double hs; //h_i
            double hsmin1; //h_i-1

            if (i > 0) { //normal case
                thetamin1 = theta[i - 1];
                ktmin1    = thetaplus[i - 1];
            }
            else { //i==0, on the left of the domain, take periodic BCs into acount
                thetamin1 = theta[ntheta_int - 1]; // equal to 2*PI-ktmin1
                ktmin1    = thetaplus[ntheta_int - 1];
            }
            hs = hplus[j];
            if (j > 0) {
                hsmin1 = hplus[j - 1];
            }
            else {
                hsmin1 = 2 * r[0]; // across the origin
            }

            // % -------------------------------------------------------------------------------------------------- %
            // 9-Point Stencil (attention order: bottom, bottom left, left, top left, top, top right, right, bottom right, middle)

            //check if we have a black point --> fine black point has two coarse neighbours and thus we only treat the middle
            //skip bottom/left/top/right if we have a fine black point
            int fine_black_point = 0;
            if (extrapol_fine_node == 1 && (smoother == 0 || smoother == 2)) {
                fine_black_point = 1;
            }
            //! bottom (r - h-) (only radial smoother) (and circle smoother with across-the-origin)
            // (as every fine point has two coarse nodes as neighbours, this entry is shifted to the rhs --> Asc_ortho)
            if ((fine_black_point == 0 && smoother > 1 && j != delete_circles) || (smoother == 0 && j == 0)) {
                //std::cout << ", bottom";
                //first row of radial smoother: link to the circle smoother --> A_sc_ortho
                //across the origin (for j=0): link of points of inner circle --> A_sc

                // if (j == 1) { // j=1 means r-h is on the boundary, for Diriclet
                //     gyro::trafo(r[j - 1], theta[i], x, y, 0);
                // }
                if (j == 1 && fabs(gyro::distBoundary(r[j - 1], theta[i], 0)) < tol_bound_check) {
                    //empty
                    //for r1: r0 is on the boundary (for Diriclet BC)
                    //we don't treat the bottom point for r1 as it is connected to the boundary (r0) and we thus bring it to the rhs
                }
                else {
                    // second circle (at least one 'real' circle with ntheta_int nodes which is closer to the origin)
                    if (j > 0) {
                        // to u at r-h geometrically recorded (i.e. -ntheta_int nodes at slice sorting)
                        col = index - n_cols_radial;
                    }
                    else { //j=0: across the origin
                        int shift;
                        if (extrapol == 1 && l == 0)
                            shift = n_cols_radial / 2;
                        else
                            shift = ntheta_int / 2;
                        if (i + 1 > ntheta_int / 2) {
                            col = index - shift; // half a turn back
                        }
                        else {
                            col = index + shift; // half a turn further
                        }
                    }
                    if (j > 0) {

                        val = -0.5 * kt / hsmin1 * (gyro::arr(r[j - 1], theta[i], 0) + gyro::arr(r[j], theta[i], 0)) -
                              0.5 * ktmin1 / hsmin1 * (gyro::arr(r[j - 1], theta[i], 0) + gyro::arr(r[j], theta[i], 0));
                    }
                    else { // j==0; across the origin, just use r_{s-1}=0

                        val = -0.5 * (kt + ktmin1) / hsmin1 *
                              (gyro::arr(r[j], theta[i], 0) + gyro::arr(r[j], theta[i] + PI, 0));
                    }

                    if (val != 0) {
                        A_Zebra_v[smoother].push_back(val);
                        A_Zebra_r[smoother].push_back(index);
                        A_Zebra_c[smoother].push_back(col);
                    }
                }
            }
            if (fine_black_point == 0) { //Extrapolation: skip, if we have a fine point (only for the black smoother)
                //! left (phi-k) (only circle smoother)
                if (smoother < 2) {
                    //std::cout << ", left";
                    val = -0.5 * hsmin1 / ktmin1 * (gyro::att(r[j], thetamin1, 0) + gyro::att(r[j], theta[i], 0)) -
                          0.5 * hs / ktmin1 * (gyro::att(r[j], thetamin1, 0) + gyro::att(r[j], theta[i], 0));

                    if (val != 0) {
                        A_Zebra_v[smoother].push_back(val);
                        A_Zebra_r[smoother].push_back(index);

                        if (i > 0) { // previous node in phi direction
                            col = index - 1;
                        }
                        else { // periodicity condition
                            col = index + ntheta_int - 1;
                        }
                        A_Zebra_c[smoother].push_back(col);
                    }
                }

                //! top (r + h+) (only radial smoother)
                if (smoother > 1) {
                    //std::cout << ", top";
                    // gyro::trafo(r[j + 1], theta[i], x, y, 0);
                    // means that r[j+1] is not on the Dirichlet boundary
                    if (j < nr_int - 1) {
                        // means that r[j+1] is not on the upper Dirichlet boundary
                        val = -0.5 * kt / hs * (gyro::arr(r[j], theta[i], 0) + gyro::arr(r[j + 1], theta[i], 0)) -
                              0.5 * ktmin1 / hs * (gyro::arr(r[j], theta[i], 0) + gyro::arr(r[j + 1], theta[i], 0));

                        if (val != 0) {
                            A_Zebra_v[smoother].push_back(val);
                            A_Zebra_r[smoother].push_back(index);
                            col = index + n_cols_radial;
                            A_Zebra_c[smoother].push_back(col);
                        }
                    }
                }

                //! right (phi+k) (only circle smoother)
                if (smoother < 2) {
                    //std::cout << ", right";
                    if (i + 1 < ntheta_int) { // next node in theta direction
                        thetap1 = theta[i + 1];
                    }
                    else { // periodicity condition
                        thetap1 = 2 * PI;
                    }
                    val = -0.5 * hs / kt * (gyro::att(r[j], theta[i], 0) + gyro::att(r[j], thetap1, 0)) -
                          0.5 * hsmin1 / kt * (gyro::att(r[j], theta[i], 0) + gyro::att(r[j], thetap1, 0));

                    if (val != 0) {
                        A_Zebra_v[smoother].push_back(val);
                        A_Zebra_r[smoother].push_back(index);

                        if (i + 1 < ntheta_int) { // next node in theta direction
                            col = index + 1;
                        }
                        else { // periodicity condition
                            col = index - ntheta_int + 1;
                        }
                        A_Zebra_c[smoother].push_back(col);
                    }
                }
            } //treat middle also in the case of a fine black point

            //! middle (r,phi) (both smoothers)
            //std::cout << ", middle" << std::endl;
            if (i + 1 < ntheta_int) { // Define theta(i+1), this might be on the periodicity condition.
                thetap1 = theta[i + 1];
            }
            else {
                thetap1 = 2 * PI;
            }

            col = index;
            if (j > 0) {

                val = 0.5 * kt / hs * (gyro::arr(r[j], theta[i], 0) + gyro::arr(r[j + 1], theta[i], 0)) +
                      0.5 * hs / kt * (gyro::att(r[j], theta[i], 0) + gyro::att(r[j], thetap1, 0)) +
                      0.5 * kt / hsmin1 * (gyro::arr(r[j - 1], theta[i], 0) + gyro::arr(r[j], theta[i], 0)) +
                      0.5 * hsmin1 / kt * (gyro::att(r[j], theta[i], 0) + gyro::att(r[j], thetap1, 0)) +
                      0.5 * ktmin1 / hsmin1 * (gyro::arr(r[j - 1], theta[i], 0) + gyro::arr(r[j], theta[i], 0)) +
                      0.5 * hsmin1 / ktmin1 * (gyro::att(r[j], thetamin1, 0) + gyro::att(r[j], theta[i], 0)) +
                      0.5 * ktmin1 / hs * (gyro::arr(r[j], theta[i], 0) + gyro::arr(r[j + 1], theta[i], 0)) +
                      0.5 * hs / ktmin1 * (gyro::att(r[j], thetamin1, 0) + gyro::att(r[j], theta[i], 0));
            }
            else { //across the origin; j=0

                val = 0.5 * kt / hs * (gyro::arr(r[j], theta[i], 0) + gyro::arr(r[j + 1], theta[i], 0)) +
                      0.5 * hs / kt * (gyro::att(r[j], theta[i], 0) + gyro::att(r[j], thetap1, 0)) +
                      0.5 * hsmin1 / kt * (gyro::att(r[j], theta[i], 0) + gyro::att(r[j], thetap1, 0)) +
                      0.5 * hsmin1 / ktmin1 * (gyro::att(r[j], thetamin1, 0) + gyro::att(r[j], theta[i], 0)) +
                      0.5 * ktmin1 / hs * (gyro::arr(r[j], theta[i], 0) + gyro::arr(r[j + 1], theta[i], 0)) +
                      0.5 * hs / ktmin1 * (gyro::att(r[j], thetamin1, 0) + gyro::att(r[j], theta[i], 0)) +
                      0.5 * kt / hsmin1 * (gyro::arr(r[j], theta[i] + PI, 0) + gyro::arr(r[j], theta[i], 0)) +
                      0.5 * ktmin1 / hsmin1 * (gyro::arr(r[j], theta[i] + PI, 0) + gyro::arr(r[j], theta[i], 0));
            }

            val += betaVec[j * ntheta_int + i];

            if (val != 0) {
                A_Zebra_v[smoother].push_back(val);
                A_Zebra_r[smoother].push_back(index);
                A_Zebra_c[smoother].push_back(index);
            }
        }
    }

    // % -------------------------------------------------------------------------------------------------- %

    //!treat the outer circle, which is on the Diriclet boundary (radial smoother)
    double tmp = 0;
    // gyro::trafo(r[nr_int], tmp, x, y, 0);
    if (fabs(gyro::distBoundary(r[nr_int], tmp, 0)) < tol_bound_check) { //check if R_max is on boundary
        // if (r[nr_int] != Rmax) {
        //     throw std::runtime_error("Node on the boundary erroneously detected as node with boundary conditions "
        //                              "(although there are none...)");
        // }

        for (int i = 0; i < ntheta_int; i++) { // iterate over all points of that circle
            //check if the point is black or white
            if (i % 2 == 0) {
                smoother = 2;
            }
            else {
                smoother = 3;
            }

            //check if the node is coarse or fine
            extrapol_fine_node = 0;
            //only in case of extrapolation and level 0
            if (extrapol == 1 && l == 0) {
                if (coarse_nodes_list_theta[i] == -1 || coarse_nodes_list_r[nr_int] == -1) {
                    extrapol_fine_node = 1; //the node is fine
                }
                else {
                    continue; //skip the index if the node is coarse
                }
            }
            /*else {
                int row = nr - delete_circles - 1; //last row = r_max-delete circles (-1 because indexing starts at zero)
                int col = i / 2;
                index = row * n_cols_radial + col;
            }*/

            count_nodes[smoother]++;
            int index = count_nodes[smoother] - 1;

            A_Zebra_r[smoother].push_back(index);
            A_Zebra_c[smoother].push_back(index);
            A_Zebra_v[smoother].push_back(1.0); //just insert "ones" on the diagonal
        }
    }
} /* ----- end of level::build_Asc0 ----- */

/*! \brief Applies the matrix A_sc_ortho for a smoother explicitely on the level l (deprecated)
 *
 * Asc_ortho corresponds to all lines of a smoother s with colour c and the columns not in Asc
 * 
 * \param smoother_todo: the smoother of this Asc_ortho matrix
*/
void level::apply_Asc_ortho0(std::vector<double>& Au, int smoother)
{
    //smoother: 0 (circle black), 1 (circle white), 2 (radial black), 3(radial white)

    // if boundary is only defined on radius, checks for wrongly detected boundary
    // nodes due to rounding errors (double checking...)
    double tol_bound_check = 1e-8;
    int extrapol           = gyro::icntl[Param::extrapolation];

    // double x, y;
    int col;
    double val;
    int n_cols_radial = ntheta_int / 2;

    int extrapol_fine_node; //indicates if we have extrapolation on level 0 AND a fine node
    // to count the fine nodes, to know the index without the coarse nodes
    std::vector<double> count_nodes = {0, 0, 0, 0}; //cb, cw, rb, rw

    //! check if we start at j=0 (across the origin) or j=1 (Diriclet boundary)
    // gyro::trafo(r[0], theta[0], x, y, 0);
    // alt: r0=0/neumBound
    int start_j;
    // dirBound([r(1),0])==0 means that (r(1),0) is on Dirichlet boundary
    if (fabs(gyro::distBoundary(r[0], theta[0], 0)) < tol_bound_check) { //check if r0 is on boundary
        // if (fabs(gyro::distBoundary(x, y, 0)) < tol_bound_check && r[0] != gyro::dcntl[Param::r0_DB]) {
        //     throw std::runtime_error("Node on the boundary erroneously detected as node with boundary conditions "
        //                              "(although there are none...)");
        // }
        start_j = 1;
        if (extrapol == 1 && l == 0)
            count_nodes[0] = n_cols_radial; //in case of the circle/black smoother, we start at the next line
        else
            count_nodes[0] = ntheta_int;
    }
    // dirBound([r(1),0])>0 means that (r(1),0) is not on Dirichlet boundary
    else {
        // if (fabs(gyro::distBoundary(x, y, 0)) > tol_bound_check && r[0] == gyro::dcntl[Param::r0_DB]) {
        //     throw std::runtime_error("Node on the boundary erroneously detected as node with boundary conditions "
        //                              "(although there are none...)");
        // }
        start_j = 0;
    }

    // % -------------------------------------------------------------------------------------------------- %
    //!Iteration over all points!
    for (int j = start_j; j < nr_int; j++) { // start with second circle (j=1), as we already treated the origin(j=0)
        for (int i = 0; i < ntheta_int; i++) {
            //! if the index and the smoother don't fit together, just skip this point!
            if ((j < delete_circles && ((j % 2 == 0 && smoother != 0) || (j % 2 != 0 && smoother != 1))) ||
                (j >= delete_circles && ((i % 2 == 0 && smoother != 2) || (i % 2 != 0 && smoother != 3)))) {
                continue;
            }

            //!index of the point in the overall total mxm matrix
            //to indicate the col-index of the matrix A_sc_ortho of size (m_sc x m)
            int base_row_index = j * ntheta_int + i;

            //!extrapolation: check if we have a fine or a coarse node
            extrapol_fine_node = 0;
            if (extrapol == 1 && l == 0) { //only do this in case of extrapolation and level 0
                if (coarse_nodes_list_theta[i] == -1 || coarse_nodes_list_r[j] == -1) {
                    extrapol_fine_node = 1; //the node is fine
                }
                else {
                    continue; //skip the index, if the node is coarse
                }
            }

            //!compute the index of the point in the grid depending on the smoother, and thus the row_index
            count_nodes[smoother]++;
            int row_index = count_nodes[smoother] - 1; //index of the point (i,j)

            double kt = thetaplus[i]; //k_j
            double thetamin1; //the coordinate of theta to the left
            double thetap1 = 0; //the coordinate of theta to the right
            double ktmin1; //k_j-1
            double hs; //h_i
            double hsmin1; //h_i-1

            if (i > 0) { //normal case
                thetamin1 = theta[i - 1];
                ktmin1    = thetaplus[i - 1];
            }
            else { //i==0, on the left of the domain, take periodic BCs into acount
                thetamin1 = theta[ntheta_int - 1]; // equal to 2*PI-ktmin1
                ktmin1    = thetaplus[ntheta_int - 1];
            }
            hs = hplus[j];
            if (j > 0) {
                hsmin1 = hplus[j - 1];
            }
            else {
                hsmin1 = 2 * r[0]; // across the origin
            }

            // % -------------------------------------------------------------------------------------------------- %
            // 9-Point Stencil (attention order: bottom, bottom left, left, top left, top, top right, right, bottom right, middle)

            //check if we have a black point
            // --> fine black point has two coarse neighbours, thus we additionally need to treat the bottom/top or left/right
            //
            int fine_black_circle = 0;
            int fine_black_radial = 0;
            if (extrapol_fine_node == 1) {
                if (smoother == 0) {
                    //the point belongs to the black/circle smoother --> additionally treat left/right
                    fine_black_circle = 1;
                }
                else if (smoother == 2) {
                    //the point belongs to the black/radial smoother --> additionally treat bottom/top
                    fine_black_radial = 1;
                }
            }

            //! bottom (r - h-) (only circle smoother)
            //Extrapolation: also for radial smoother in case of a fine black node
            if ((smoother < 2 && j < delete_circles && j > 0) || (smoother > 1 && j == delete_circles) ||
                fine_black_radial == 1 || (fine_black_circle == 1 && j > 0)) {
                //for the circle smoother (but not the inner circle)
                //and for the radial smoother for the first line, which has a bottom link to the circle smoother

                // if (j == 1) {
                //     gyro::trafo(r[j - 1], theta[i], x, y, 0);
                // }
                // j=1 means r-h is on the boundary
                if (j == 1 && fabs(gyro::distBoundary(r[j - 1], theta[i], 0)) < tol_bound_check) {
                    //empty
                    //for r1: r0 is on the boundary for Diriclet BC
                    //we don't treat the bottom point for r1 as it is connected to the boundary (r0) and we thus bring it to the rhs
                }
                else {
                    // second circle (at least one 'real' circle with ntheta_int nodes which is closer to the origin)
                    col = base_row_index - ntheta_int;
                    val = -0.5 * kt / hsmin1 * (gyro::arr(r[j - 1], theta[i], 0) + gyro::arr(r[j], theta[i], 0)) -
                          0.5 * ktmin1 / hsmin1 * (gyro::arr(r[j - 1], theta[i], 0) + gyro::arr(r[j], theta[i], 0));
                    if (val != 0) {
                        if (gyro::icntl[Param::smoother] == 13) {
                            if (smoother < 2) { //circle

                                Au[row_index] += val * u_previous_c[col];
                            }
                            else { //radial

                                Au[row_index] += val * u_previous_r[col];
                            }
                        }
                        else {

                            Au[row_index] += val * u[col];
                        }
                    }
                }
            }

            //! bottom left (r-h-,phi-k) (for both smoothers)
            if (j > 0) { //not for j=0
                // if (j == 1) { // coords necessary for potential boundary conditions
                //     gyro::trafo(r[j - 1], thetamin1, x, y, 0); // r(j) is PREVIOUS r, actual is r(J+1) !
                // }
                // j=1 means r-h is on the boundary
                if (j == 1 && fabs(gyro::distBoundary(r[j - 1], thetamin1, 0)) < tol_bound_check) {
                    //empty
                    //no treatment of bottom point for r1 in the case of Diriclet
                }
                else {
                    if (i > 0) { // next node in theta direction but one circle before
                        col = base_row_index - ntheta_int - 1;
                    }
                    else { // periodicity condition
                        col = base_row_index - 1;
                    }
                    val = -(gyro::art(r[j], thetamin1, 0) + gyro::art(r[j - 1], theta[i], 0));

                    if (val != 0) {
                        if (gyro::icntl[Param::smoother] == 13) {
                            if (smoother < 2) { //circle

                                Au[row_index] += val * u_previous_c[col];
                            }
                            else { //radial

                                Au[row_index] += val * u_previous_r[col];
                            }
                        }
                        else {

                            Au[row_index] += val * u[col];
                        }
                    }
                }
            }

            //! left (phi-k) (radial smoother only)
            //Extrapolation: also for circle smoother in case of a fine node
            if ((smoother > 1 && j >= delete_circles) || fine_black_circle == 1) {
                if (i > 0) { // previous node in phi direction
                    col = base_row_index - 1;
                }
                else { // periodicity condition
                    col = base_row_index + ntheta_int - 1;
                }
                val = -0.5 * hsmin1 / ktmin1 * (gyro::att(r[j], thetamin1, 0) + gyro::att(r[j], theta[i], 0)) -
                      0.5 * hs / ktmin1 * (gyro::att(r[j], thetamin1, 0) + gyro::att(r[j], theta[i], 0));

                if (val != 0) {
                    if (gyro::icntl[Param::smoother] == 13) {
                        if (smoother < 2) { //circle

                            Au[row_index] += val * u_previous_c[col];
                        }
                        else { //radial

                            Au[row_index] += val * u_previous_r[col];
                        }
                    }
                    else {

                        Au[row_index] += val * u[col];
                    }
                }
            }

            //! top left (r+h+,phi-k) (for both smoothers)
            double theta_eval_prelast = 0;
            if (i > 0) { // not the first node on the corresponding circle; we can take theta(i-1)
                // gyro::trafo(r[j + 1], theta[i - 1], x, y, 0);
            }
            else {
                theta_eval_prelast = theta[ntheta_int - 1];
                // gyro::trafo(r[j + 1], theta_eval_prelast, x, y, 0);
            }

            if (j < nr_int - 1) {
                if (i > 0) { // previous node in phi direction but at r+h
                    col = base_row_index + ntheta_int - 1;
                }
                else { // first node on corresponding circle; pay gyro::attention to periodicity condition!
                    col = base_row_index + 2 * ntheta_int - 1;
                }
                val = gyro::art(r[j], thetamin1, 0) + gyro::art(r[j + 1], theta[i], 0);

                if (val != 0) {
                    if (gyro::icntl[Param::smoother] == 13) {
                        if (smoother < 2) { //circle
                            Au[row_index] += val * u_previous_c[col];
                        }
                        else { //radial
                            Au[row_index] += val * u_previous_r[col];
                        }
                    }
                    else {

                        Au[row_index] += val * u[col];
                    }
                }
            }

            //! top (r + h+) (circle smoother only)
            //Extrapolation: also for radial smoother in case of a fine black node
            if ((smoother < 2 && j < delete_circles) || fine_black_circle == 1 || fine_black_radial == 1) {
                // gyro::trafo(r[j + 1], theta[i], x, y, 0);
                // means that r[j+1] is not on the Dirichlet boundary
                if (j < nr_int - 1) {
                    col = base_row_index + ntheta_int;
                    val = -0.5 * kt / hs * (gyro::arr(r[j], theta[i], 0) + gyro::arr(r[j + 1], theta[i], 0)) -
                          0.5 * ktmin1 / hs * (gyro::arr(r[j], theta[i], 0) + gyro::arr(r[j + 1], theta[i], 0));

                    if (val != 0) {
                        if (gyro::icntl[Param::smoother] == 13) {
                            if (smoother < 2) { //circle

                                Au[row_index] += val * u_previous_c[col];
                            }
                            else { //radial

                                Au[row_index] += val * u_previous_r[col];
                            }
                        }
                        else {

                            Au[row_index] += val * u[col];
                        }
                    }
                }
            }

            //! top right (r+h+,phi+k) (for both smoothers)
            if (i + 1 < ntheta_int) { // previous node in phi direction but at r+h
                thetap1 = theta[i + 1];
            }
            else { // first node on corresponding circle; pay gyro::attention to periodicity condition!
                thetap1 = 2 * PI;
            }

            // gyro::trafo(r[j + 1], thetap1, x, y, 0);
            // means that r[j+1] is not on the Dirichlet boundary
            if (j < nr_int - 1) {
                if (i + 1 < ntheta_int) { // previous node in phi direction but at r+h
                    col = base_row_index + ntheta_int + 1;
                }
                else { // first node on corresponding circle; pay gyro::attention to periodicity condition!
                    col = base_row_index + 1;
                }
                val = -gyro::art(r[j + 1], theta[i], 0) - gyro::art(r[j], thetap1, 0);

                if (val != 0) {
                    if (gyro::icntl[Param::smoother] == 13) {
                        if (smoother < 2) { //circle
                            Au[row_index] += val * u_previous_c[col];
                        }
                        else { //radial
                            Au[row_index] += val * u_previous_r[col];
                        }
                    }
                    else {

                        Au[row_index] += val * u[col];
                    }
                }
            }

            //! right (phi+k) (radial smoother only)
            //Extrapolation: also for circle smoother in case of a fine black node
            if ((smoother > 1 && j >= delete_circles) || fine_black_circle == 1) {
                if (i + 1 < ntheta_int) { // next node in theta direction
                    col     = base_row_index + 1;
                    thetap1 = theta[i + 1];
                }
                else { // periodicity condition
                    col     = base_row_index - ntheta_int + 1;
                    thetap1 = 2 * PI;
                }
                val = -0.5 * hs / kt * (gyro::att(r[j], theta[i], 0) + gyro::att(r[j], thetap1, 0)) -
                      0.5 * hsmin1 / kt * (gyro::att(r[j], theta[i], 0) + gyro::att(r[j], thetap1, 0));

                if (val != 0) {
                    if (gyro::icntl[Param::smoother] == 13) {
                        if (smoother < 2) { //circle

                            Au[row_index] += val * u_previous_c[col];
                        }
                        else { //radial

                            Au[row_index] += val * u_previous_r[col];
                        }
                    }
                    else {

                        Au[row_index] += val * u[col];
                    }
                }
            }

            //! bottom right (r-h-,phi+k) (for both smoothers)
            if (j > 0) { //not for j=0
                double r_tmp, theta_tmp;
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
                if (i + 1 < ntheta_int) {
                    thetap1 = theta[i + 1];
                }
                else {
                    thetap1 = 2 * PI;
                }

                if (j == 1 && fabs(gyro::distBoundary(r_tmp, theta_tmp, 0)) < tol_bound_check) {
                    //empty
                    //no treatment of bottom point for r1 in the case of Diriclet
                }
                else {
                    if (i + 1 < ntheta_int) { // next node in theta direction but one circle before
                        col = base_row_index - ntheta_int + 1;
                    }
                    else { // periodicity condition
                        col = base_row_index - 2 * ntheta_int + 1;
                    }
                    val = gyro::art(r[j - 1], theta[i], 0) + gyro::art(r[j], thetap1, 0);

                    if (val != 0) {
                        if (gyro::icntl[Param::smoother] == 13) {
                            if (smoother < 2) { //circle

                                Au[row_index] += val * u_previous_c[col];
                            }
                            else { //radial

                                Au[row_index] += val * u_previous_r[col];
                            }
                        }
                        else {

                            Au[row_index] += val * u[col];
                        }
                    }
                }
            }
        }
    }
} /* ----- end of level::apply_Asc_ortho0 ----- */
