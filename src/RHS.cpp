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
 * \file RHS.cpp
 * \brief Implementation of the right hand side
 * \author M. Kuehn, C. Kruse, P. Leleux, C. Schwarz
 * \version 0.0
 */
#include "gyro.h"

/*!
 *  \brief Defines the RHS (single)
 *
 * Defines the RHS on (r, theta)
 *
 */
double gyro::eval_def_rhs(double r, double theta, int verbose)
{
    double kappa_eps = dcntl[Param::kappa_eps];
    double delta_e   = dcntl[Param::delta_e];
    double Rmax      = dcntl[Param::R];

    double rhs_val = 0;
    if (r > 0) {
        rhs_val = functions->rho_glob(r, theta, kappa_eps, delta_e, Rmax);
    }
    else {
        rhs_val = functions->rho_pole(r, theta, kappa_eps, delta_e, Rmax);
    }
    if (verbose) {
        std::cout << "RHS(" << r << ", " << theta << "): " << rhs_val << "\n";
    }
    return rhs_val;
} /* ----- end of gyro::eval_def_rhs ----- */

/*!
 *  \brief Defines the RHS (vector)
 *
 * Defines the RHS for a whole radius
 *
 */
std::vector<double> gyro::eval_def_rhs(double r, std::vector<double> theta, std::vector<double> sin_theta,
                                       std::vector<double> cos_theta, int ntheta, int verbose)
{
    double kappa_eps = dcntl[Param::kappa_eps];
    double delta_e   = dcntl[Param::delta_e];
    double Rmax      = dcntl[Param::R];

    std::vector<double> rhs_val(ntheta);
    if (r > 0) {
        functions->rho_glob(r, theta, kappa_eps, delta_e, Rmax, rhs_val, sin_theta, cos_theta);
    }
    else {
        functions->rho_pole(r, theta, kappa_eps, delta_e, Rmax, rhs_val, sin_theta, cos_theta);
    }
    if (verbose) {
        for (int i = 0; i < ntheta; i++)
            std::cout << "RHS(" << r << ", " << theta[i] << "): " << rhs_val[i] << "\n";
    }
    return rhs_val;
} /* ----- end of gyro::eval_def_rhs ----- */
