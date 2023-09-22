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
 * \file define_coarse_nodes.cpp
 * \brief Implementation of the coarsening
 * \author M. Kuehn, C. Kruse, P. Leleux, C. Schwarz
 * \version 0.0
 */
#include "gmgpolar.h"
#include "level.h"

/*!
 *  \brief Define the coarse nodes on all levels
 *
 * Defines coarse nodes on all level (two columns per level) by recursively
 * calling define_coarse_nodes_onelevel(); using a vector of arrays r/theta
 * to store the different levels.
 * - nodes are ordered (per level) circlewise starting with r_min (here r_0, in 
 * article denoted r_1)
 * - [-1, -1]: means that this node is not a coarse node
 *   [m, n]: means that this node (with r_m+1 and theta_n) is coarse node
 *
 */
void gmgpolar::define_coarse_nodes()
{
    // int periodic = gyro::icntl[Param::periodic];

    for (int i = 1; i < levels; i++) {
        level* new_level = new level(i);
        v_level.push_back(new_level);
    }
    levels_orig = levels;
    for (int i = 0; i < levels; i++) {
        if (gyro::icntl[Param::verbose] > 3)
            std::cout << "=> Defining coarse nodes on level " << i << "\n";
        v_level[i]->store_theta_n_co();
        if (i < levels - 1) {
            v_level[i]->define_coarse_nodes_onelevel(v_level[i + 1]);
        }
        if (gyro::icntl[Param::verbose] > 3) {
            std::cout << "level: " << i;
            if (i != levels - 1)
                std::cout << ", coarse_nodes: " << v_level[i]->coarse_nodes;
            std::cout << std::endl;
            std::cout << "nr: " << v_level[i]->nr << ", ntheta: " << v_level[i]->ntheta << "\n";
        }
        if (i != levels - 1 &&
            (v_level[i]->nr % 2 != 1 || v_level[i]->nr < 5 || v_level[i]->ntheta % 2 != 0 || v_level[i]->ntheta < 8)) {
            std::cout << "WARNING: Cannot reach the desired number of levels (" << levels << "). Replacing it with "
                      << i + 1 << "...\n";
            ;
            levels = i + 1;
            break;
        }
    }
} /* ----- end of gmgpolar::define_coarse_nodes ----- */

/*!
 *  \brief Define the coarse nodes on one level
 *
 * defines coarse nodes on one level (two columns)
 * - nodes are ordered circlewise starting with r_min (here r_0, in article 
 * denoted r_1)
 * - [-1, -1]: means that this node is not a coarse node
 *   [m, n]: means that this node (with r_m+1 and theta_n) is coarse node
 * 
 * \param coarser: the next coarser level
 *
 */
void level::define_coarse_nodes_onelevel(level* coarser)
{
    // int periodic = gyro::icntl[Param::periodic];

    // int count_coarse = 0, type;
    std::set<int> indices_r;
    std::set<int> indices_theta;
    coarse_nodes            = (nr_int + 1) * ntheta_int;
    coarse_nodes_list_r     = std::vector<int>(coarse_nodes);
    coarse_nodes_list_theta = std::vector<int>(coarse_nodes);
    for (int j = 0; j < nr_int + 1; j++) {
        for (int i = 0; i < ntheta_int; i++) {

            if (i % 2 == 0 && j % 2 == 0) {
                coarse_nodes_list_r[j * ntheta_int + i]     = j;
                coarse_nodes_list_theta[j * ntheta_int + i] = i;

                indices_r.insert(j);
                indices_theta.insert(i);
            }
            else {
                coarse_nodes_list_r[j * ntheta_int + i]     = -1;
                coarse_nodes_list_theta[j * ntheta_int + i] = -1;
            }
        }
    }

    std::set<int, std::greater<int>>::iterator itr_r, itr_theta;

    coarser->nr     = (int)indices_r.size();
    coarser->r      = std::vector<double>(coarser->nr);
    coarser->ntheta = (int)indices_theta.size();
    coarser->theta  = std::vector<double>(coarser->ntheta);

    itr_r = indices_r.begin();
    for (int i = 0; i < coarser->nr; i++) {
        coarser->r[i] = r[*itr_r];
        itr_r++;
    }

    itr_theta = indices_theta.begin();
    for (int i = 0; i < coarser->ntheta; i++) {
        coarser->theta[i] = theta[*itr_theta];
        itr_theta++;
    }

    if (gyro::icntl[Param::verbose] > 5) {
        display_r();
        display_theta();
        std::cout << "Coarse grid: \n";
        gyro::disp(coarse_nodes_list_r, "coarse_nodes_list_r");
        gyro::disp(coarse_nodes_list_theta, "coarse_nodes_list_theta");
    }
} /* ----- end of level::define_coarse_nodes_onelevel ----- */

/*!
 *  \brief Define the coarse nodes on one level
 *
 * defines coarse nodes on one level (two columns)
 * - nodes are ordered circlewise starting with r_min (here r_0, in article 
 * denoted r_1)
 * - [-1, -1]: means that this node is not a coarse node
 *   [m, n]: means that this node (with r_m+1 and theta_n) is coarse node
 *
 */
void level::store_theta_n_co()
{
    // Define number of intervals
    nr_int     = nr - 1;
    ntheta_int = ntheta;
    if (fabs(theta[ntheta - 1] - 2 * PI) < 1e-10) // check if end is 2*Pi (for periodic boundary conditions)
        ntheta_int--;
    if (ntheta_int % 2)
        std::cout << "WARNING: The number of thetas is odd. Please use even numbers only.\n";

    // To take first/last theta into account without conditions:
    // - theta_per = [theta_last, theta, 2PI]
    theta_per    = std::vector<double>(ntheta_int + 2);
    theta_per[0] = theta[ntheta_int - 1];
    for (int i = 0; i < ntheta_int; i++)
        theta_per[i + 1] = theta[i];
    theta_per[ntheta_int + 1] = 2 * PI;

    cos_theta_per = std::vector<double>(ntheta_int + 2);
    sin_theta_per = std::vector<double>(ntheta_int + 2);
    for (int i = 0; i < ntheta_int + 2; i++) {
        cos_theta_per[i] = cos(theta_per[i]);
        sin_theta_per[i] = sin(theta_per[i]);
    }

    // Store cosines/sines (expensive)
    cos_theta = std::vector<double>(ntheta);
    for (int i = 0; i < ntheta; i++) {
        cos_theta[i] = cos(theta[i]);
    }
    sin_theta = std::vector<double>(ntheta);
    for (int i = 0; i < ntheta; i++) {
        sin_theta[i] = sin(theta[i]);
    }

    // Store theta+PI (only O(ntheta)) and associated cosines/sines (expensive)
    theta_PI = std::vector<double>(ntheta);
    for (int i = 0; i < ntheta; i++) {
        theta_PI[i] = theta[i] + PI;
    }
    cos_theta_PI = std::vector<double>(ntheta);
    for (int i = 0; i < ntheta; i++) {
        cos_theta_PI[i] = cos(theta_PI[i]);
    }
    sin_theta_PI = std::vector<double>(ntheta);
    for (int i = 0; i < ntheta; i++) {
        sin_theta_PI[i] = sin(theta_PI[i]);
    }

    // Size of intervals in theta
    // thetaplus = [k_last, thetaplus, k_last]
    thetaplus = std::vector<double>(ntheta_int);
    for (int i = 0; i < ntheta_int - 1; i++)
        thetaplus[i] = theta[i + 1] - theta[i];
    thetaplus[ntheta_int - 1] = 2 * PI - theta[ntheta_int - 1];

    // Size of intervals in theta
    // thetaplus_per = [k_last, thetaplus, k_last]
    thetaplus_per = std::vector<double>(ntheta_int + 1);
    for (int i = 1; i < ntheta_int; i++)
        thetaplus_per[i] = theta_per[i + 1] - theta_per[i];
    thetaplus_per[ntheta_int] = 2 * PI - theta_per[ntheta_int];
    thetaplus_per[0]          = thetaplus_per[ntheta_int];

    // Size of intervals in r
    hplus = std::vector<double>(nr_int);
    for (int i = 0; i < nr_int; i++)
        hplus[i] = r[i + 1] - r[i];

    if (gyro::icntl[Param::verbose] > 5) {
        gyro::disp(theta_per, "theta_per");
        gyro::disp(cos_theta, "cos_theta");
        gyro::disp(sin_theta, "sin_theta");
        gyro::disp(theta_PI, "theta_PI");
        gyro::disp(cos_theta_PI, "cos_theta_PI");
        gyro::disp(sin_theta_PI, "sin_theta_PI");
        gyro::disp(thetaplus, "thetaplus");
        gyro::disp(thetaplus_per, "thetaplus_per");
        gyro::disp(hplus, "hplus");
    }
} /* ----- end of level::store_theta_n_co ----- */
