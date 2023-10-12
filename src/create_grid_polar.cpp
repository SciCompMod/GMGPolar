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
 * \file create_grid_polar.cpp
 * \brief Implementation of the creation of the polar grid (r and theta)
 * \author M. Kuehn, C. Kruse, P. Leleux, C. Schwarz
 * \version 0.0
 */
#include "gmgpolar.h"
#include "level.h"

/*!
 *  \brief Create the polar grids
 *
 *  Create the polar finest grids
 * 
 */
void gmgpolar::create_grid_polar()
{
    if (gyro::icntl[Param::verbose] > 3)
        std::cout << "Creating polar grid...\n";

    level* new_level = new level(0);
    v_level.push_back(new_level);
    if (gyro::f_grid_r.empty() || gyro::icntl[Param::write_radii_angles] == 1) {
        v_level[0]->build_r();
    }
    else {
        v_level[0]->read_grid_r();
    }
    if (gyro::f_grid_theta.empty() || gyro::icntl[Param::write_radii_angles] == 1) {
        v_level[0]->build_theta();
    }
    else {
        v_level[0]->read_grid_theta();
    }

    if (v_level[0]->ntheta < 2)
        throw std::runtime_error("Choose size theta > 3!");

    gyro::dcntl[Param::R0]     = v_level[0]->r[0];
    gyro::dcntl[Param::R]      = v_level[0]->r[v_level[0]->nr - 1];
    gyro::dcntl[Param::THETA0] = v_level[0]->theta[0];
    gyro::dcntl[Param::THETA]  = v_level[0]->theta[v_level[0]->ntheta - 1];

    //split grid intervals in r and theta direction recursively in the middle to get a finer grid
    if (gyro::icntl[Param::divideBy2] != 0) {
        create_grid_polar_divide();
    }

    if (gyro::icntl[Param::verbose] > 1) {
        std::cout << "System of size (nr x ntheta) = (" << v_level[0]->nr << " x " << v_level[0]->ntheta << ")\n";
        std::cout << "on the coordinates (r x theta): (" << gyro::dcntl[Param::R0] << ", " << gyro::dcntl[Param::R]
                  << ") x (" << gyro::dcntl[Param::THETA0] << ", " << gyro::dcntl[Param::THETA] << ")\n";
    }
    if (gyro::icntl[Param::verbose] > 5) {
        v_level[0]->display_r();
        v_level[0]->display_theta();
        std::cout << "\n";
    }

    if (!gyro::f_grid_r.empty() && gyro::icntl[Param::write_radii_angles] == 1) {
        v_level[0]->write_grid_r();
    }
    if (!gyro::f_grid_theta.empty() && gyro::icntl[Param::write_radii_angles] == 1) {
        v_level[0]->write_grid_theta();
    }
} /* ----- end of gmgpolar::create_grid_polar ----- */

/*!
 *  \brief Build the array r
 *
 *  Build the array r, containing the coordinates in the direction r,
 *  following an iso- or aniso-tropic division.
 */
void level::build_r()
{
    int nr_exp = gyro::icntl[Param::nr_exp];
    int aniso  = gyro::icntl[Param::fac_ani];
    double R   = gyro::dcntl[Param::R];
    double R0  = gyro::dcntl[Param::R0];

    if (l > 0)
        throw std::runtime_error("Function build_r should only be called on the finest grid");
    std::vector<double> r_tmp;
    if (!aniso) {
        nr         = pow(2, nr_exp - 1);
        double fac = (R - R0) / nr;
        r_tmp      = std::vector<double>(nr + 1);
        for (int i = 0; i < nr; i++) {
            r_tmp[i] = R0 + i * fac;
        }
        r_tmp[nr] = R;
        nr++;
    }
    else {
        // 1) uniform division with nr=2^dummy_lognr - 2^aniso
        // 2) remaining nodes are added by refining the part centered around 2/3 of r

        std::set<double, std::greater<double>>::iterator itr, itr_p1;
        // very ugly anisotropy hack.... dividing recursively smaller and smaller number of cells

        /* uniform division of r in 2^nr_exp - 2^aniso */
        int dummy_lognr  = nr_exp;
        int n_elems_equi = pow(2, dummy_lognr) - pow(2, aniso);
        if (aniso < 0 || n_elems_equi <= 0) {
            throw std::runtime_error("Please choose anisotropy factor a such that 2^fac_ani < 2^nr_exp.\n");
        }

        if ((aniso % 2) == 1) // odd number of elements on an open circular disk is desired because of coarsening
            n_elems_equi++;
        double fac                 = (R - R0) / n_elems_equi;
        nr                         = n_elems_equi + 1;
        std::vector<double> r_tmp2 = std::vector<double>(nr);
        for (int i = 0; i < nr - 1; i++)
            r_tmp2[i] = R0 + i * fac;
        r_tmp2[nr - 1] = R;

        /* refine around 2/3 of r */
        int n_elems_refined = pow(2, aniso);

        // edge
        int se;

        // allow refining of the grid at r_jump, the center point of the 
        // drop of the diffusion coefficient alpha.
        double r_jump;
        if (gyro::icntl[Param::alpha_coeff] == SONNENDRUCKER) {
            // The center of the coefficient jump lies at 0.6888697651782026
            // for backward stability with previous runs and the Matlab code, 
            // we use 0.66 though.
            r_jump = 0.66;
        } else if (gyro::icntl[Param::alpha_coeff] == ZONI) {
            r_jump = 0.4837;
        } else if (gyro::icntl[Param::alpha_coeff] == ZONI_SHIFTED) {
            // Choose center point of descent.
            // a) - ln(0.5 * (alpha(0) - alpha(Rmax))):
            //    - ln(0.5 * (np.exp(-np.tanh(-14)) - np.exp(-np.tanh(6)))) = 0.16143743821247852
            // b) r_center = Rmax * (np.arctanh(0.16143743821247852) + 14) / 20 = 0.7081431124450334 Rmax
            r_jump = 0.7081;
        } else if (gyro::icntl[Param::alpha_coeff] == POISSON) {
            r_jump = 0.5; // There is no jump for Poisson so this is an arbitrary choice
        } else {
            throw std::runtime_error("Unknown alpha coeff");
        }
        se = floor(nr * r_jump) - n_elems_refined / 2;
        int ee = se + n_elems_refined;
        // takeout
        int st = ceil((double)n_elems_refined / 4.0 + 1) - 1;
        int et = floor(3 * ((double)n_elems_refined / 4.0));

        std::set<double> r_set;
        std::set<double> r_set_p1;
        int count = 0;
        for (int i = 0; i < n_elems_refined; i++) {
            r_set_p1.insert(r_tmp2[se + i]);
            count++;
        }
        double half = fac / 2.0;
        for (int k = 0; k < aniso; k++) {
            std::set<double> r_set_p1_tmp;
            itr_p1     = r_set_p1.begin();
            int r_size = count;
            count      = 0;
            for (int i = 0; i < r_size - 1; i++) {
                r_set.insert((*itr_p1) + half);
                if (k < aniso - 1 && i >= st && i < et) {
                    r_set_p1_tmp.insert(*(itr_p1));
                    r_set_p1_tmp.insert(*(itr_p1) + half);
                    count += 2;
                }
                itr_p1++;
            }
            r_set_p1 = r_set_p1_tmp;
            half *= 0.5;
        }

        // such that the total size is 8*x+1 (or we do not refine)
        nr        = nr + r_set.size();
        int shift = 0;
        shift     = std::min(nr % 8 - 1, (int)r_set.size());
        itr       = r_set.begin();
        std::advance(itr, shift);
        r_set.erase(r_set.begin(), itr);
        for (int i = 0; i < n_elems_refined; i++)
            r_set.insert(r_tmp2[se + i]);

        // group all in r_tmp
        nr = n_elems_equi - n_elems_refined + r_set.size() + 1;

        r_tmp = std::vector<double>(nr);
        for (int i = 0; i < se; i++)
            r_tmp[i] = r_tmp2[i];
        itr = r_set.begin();
        for (int i = 0; i < (int)r_set.size(); i++) {
            r_tmp[se + i] = *itr;
            itr++;
        }
        for (int i = 0; i < n_elems_equi - ee + 1; i++)
            r_tmp[se + r_set.size() + i] = r_tmp2[ee + i];
    }

    // refine again by the middle
    nr = 2 * nr - 1;
    r  = std::vector<double>(nr);
    for (int i = 0; i < nr; i++) {
        if (!(i % 2))
            r[i] = r_tmp[i / 2];
        else
            r[i] = 0.5 * (r_tmp[(i - 1) / 2] + r_tmp[(i + 1) / 2]);
    }
} /* ----- end of level::build_r ----- */

/*!
 *  \brief Build the array theta
 *
 *  Build the array theta, containing the coordinates in the direction theta,
 *  following an iso- or aniso-tropic division.
 */
void level::build_theta()
{
    // int ntheta_exp = gyro::icntl[Param::ntheta_exp];
    int periodic = gyro::icntl[Param::periodic];
    int aniso    = gyro::icntl[Param::theta_aniso];

    if (l > 0)
        throw std::runtime_error("Function build_theta should only be called on the finest grid");
    if (aniso < 0 || aniso > 1)
        throw std::runtime_error("Please choose anisotropy on theta as 0 or 1.\n");
    if (!aniso) {
        ntheta     = pow(2, ceil(log2(nr)));
        double fac = 2 * PI / ntheta;

        // If not periodic BC on theta, last point is present (not implicit)
        if (!periodic)
            ntheta++;

        theta = std::vector<double>(ntheta);
        for (int i = 0; i < ntheta; i++)
            theta[i] = fac * i;
    }
    else {
        if (gyro::icntl[Param::verbose] > 2) {
            std::cout << "Anisotropy chosen in theta.\n";
        }
        ntheta     = pow(2, ceil(log2(nr)) - 1);
        double fac = 2 * PI / ntheta;

        ntheta++;
        std::vector<double> theta_tmp = std::vector<double>(ntheta);
        for (int i = 0; i < ntheta; i++)
            theta_tmp[i] = fac * i;
        for (int i = 1; i < ntheta - 3; i += 3)
            theta_tmp[i] = theta_tmp[i] + 0.3 * (theta_tmp[i + 1] - theta_tmp[i]);
        for (int i = 1; i < ntheta - 3; i += 5)
            theta_tmp[i] = theta_tmp[i] + 0.5 * (theta_tmp[i + 1] - theta_tmp[i]);

        ntheta = 2 * ntheta - 1;
        // If periodic BC on theta, last point is implicit
        if (periodic)
            ntheta--;

        theta = std::vector<double>(ntheta);
        for (int i = 0; i < ntheta; i++)
            if (!(i % 2))
                theta[i] = theta_tmp[i / 2];
            else
                theta[i] = 0.5 * (theta_tmp[(i - 1) / 2] + theta_tmp[(i + 1) / 2]);
    }
} /* ----- end of level::build_theta ----- */

/*!
 *  \brief Divides the intervals of the finest grid in half
 *
 * For a previsouly created finest polar grid, divides the intervals in half
 * in both r and theta directions. Mainly used in order to compute order
 * of approximations in tests (e.G. "implicitly extrapolated ..." paper).
 *
 */
void gmgpolar::create_grid_polar_divide()
{
    if (gyro::icntl[Param::verbose] > 3)
        std::cout << "Dividing a coarser grid...\n";

    int divide = gyro::icntl[Param::divideBy2];
    std::vector<double> r_tmp;
    std::vector<double> theta_tmp;

    r_tmp       = v_level[0]->r;
    theta_tmp   = v_level[0]->theta;
    double last = v_level[0]->theta[v_level[0]->ntheta - 1];

    //! divide the initial coarse grid at the center of the intervals
    for (int k = 0; k < divide; k++) { //do this 'divide' times
        std::vector<double> r_tmp2;
        std::vector<double> theta_tmp2;

        //r
        v_level[0]->nr = 2 * v_level[0]->nr - 1; //double the size of nr (fill gaps with new points)
        r_tmp2         = std::vector<double>(v_level[0]->nr); //new list for the coordinates (of larger size)
        for (int i = 0; i < v_level[0]->nr; i++) {
            if (!(i % 2))
                r_tmp2[i] = r_tmp[i / 2]; //use the available point
            else
                r_tmp2[i] = 0.5 * (r_tmp[(i - 1) / 2] + r_tmp[(i + 1) / 2]); //take the middle value
        }
        r_tmp = r_tmp2;

        //theta
        if ((!fabs(last - 2 * PI)) < 1e-10 && gyro::icntl[Param::periodic]) {
            v_level[0]->ntheta++;
            theta_tmp2 = std::vector<double>(v_level[0]->ntheta);
            for (int i = 0; i < v_level[0]->ntheta - 1; ++i) {
                theta_tmp2[i] = theta_tmp[i];
            }
            theta_tmp2[v_level[0]->ntheta - 1] = 2 * PI; //add 2*PI to the theta array
            theta_tmp                          = theta_tmp2;
        }

        v_level[0]->ntheta = 2 * v_level[0]->ntheta - 1; //double the size of ntheta (fill gaps with new points)
        theta_tmp2         = std::vector<double>(v_level[0]->ntheta);
        for (int i = 0; i < v_level[0]->ntheta; i++) {
            if (!(i % 2))
                theta_tmp2[i] = theta_tmp[i / 2]; //use the available point
            else
                theta_tmp2[i] = 0.5 * (theta_tmp[(i - 1) / 2] + theta_tmp[(i + 1) / 2]); //take the middle value
        }

        if (gyro::icntl[Param::periodic]) {
            v_level[0]->ntheta--;
            theta_tmp = std::vector<double>(v_level[0]->ntheta);
            for (int i = 0; i < v_level[0]->ntheta; ++i) {
                theta_tmp[i] = theta_tmp2[i]; //leave out last value (2*PI) again
            }
        }
        else {
            theta_tmp = theta_tmp2;
        }
    }

    v_level[0]->r     = r_tmp;
    v_level[0]->theta = theta_tmp;

    gyro::dcntl[Param::R0]     = v_level[0]->r[0];
    gyro::dcntl[Param::R]      = v_level[0]->r[v_level[0]->nr - 1];
    gyro::dcntl[Param::THETA0] = v_level[0]->theta[0];
    gyro::dcntl[Param::THETA]  = v_level[0]->theta[v_level[0]->ntheta - 1];
} /* ----- end of level::create_grid_polar_divide ----- */
