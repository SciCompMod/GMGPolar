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
 * \file gyro.cpp
 * \brief Header for the class gyro
 * \author M. Kuehn, C. Kruse, P. Leleux, C. Schwarz
 * \version 0.0
 */
#include "gyro.h"
#include "CartesianR2GyroSonnendruckerCircular.h"
#include "CartesianR2GyroSonnendruckerShafranov.h"
#include "CartesianR2GyroSonnendruckerTriangular.h"
#include "PolarR6GyroSonnendruckerCircular.h"
#include "PolarR6GyroSonnendruckerShafranov.h"
#include "PolarR6GyroSonnendruckerTriangular.h"
#include "CartesianR6GyroSonnendruckerCircular.h"
#include "CartesianR6GyroSonnendruckerShafranov.h"
#include "CartesianR6GyroSonnendruckerTriangular.h"
#include "CartesianR2GyroZoniCircular.h"
#include "CartesianR2GyroZoniShafranov.h"
#include "CartesianR2GyroZoniTriangular.h"
#include "PolarR6GyroZoniCircular.h"
#include "PolarR6GyroZoniShafranov.h"
#include "PolarR6GyroZoniTriangular.h"
#include "CartesianR6GyroZoniCircular.h"
#include "CartesianR6GyroZoniShafranov.h"
#include "CartesianR6GyroZoniTriangular.h"
#include "CartesianR2GyroZoniShiftedCircular.h"
#include "CartesianR2GyroZoniShiftedShafranov.h"
#include "CartesianR2GyroZoniShiftedTriangular.h"
#include "PolarR6GyroZoniShiftedCircular.h"
#include "PolarR6GyroZoniShiftedShafranov.h"
#include "PolarR6GyroZoniShiftedTriangular.h"
#include "PolarR6GyroZoniShiftedCulham.h"
#include "CartesianR6GyroZoniShiftedCircular.h"
#include "CartesianR6GyroZoniShiftedShafranov.h"
#include "CartesianR6GyroZoniShiftedTriangular.h"
#include "CartesianR2SonnendruckerCircular.h"
#include "CartesianR2SonnendruckerShafranov.h"
#include "CartesianR2SonnendruckerTriangular.h"
#include "PolarR6SonnendruckerCircular.h"
#include "PolarR6SonnendruckerShafranov.h"
#include "PolarR6SonnendruckerTriangular.h"
#include "CartesianR6SonnendruckerCircular.h"
#include "CartesianR6SonnendruckerShafranov.h"
#include "CartesianR6SonnendruckerTriangular.h"
#include "CartesianR2ZoniCircular.h"
#include "CartesianR2ZoniShafranov.h"
#include "CartesianR2ZoniTriangular.h"
#include "PolarR6ZoniCircular.h"
#include "PolarR6ZoniShafranov.h"
#include "PolarR6ZoniTriangular.h"
#include "CartesianR6ZoniCircular.h"
#include "CartesianR6ZoniShafranov.h"
#include "CartesianR6ZoniTriangular.h"
#include "CartesianR2ZoniShiftedCircular.h"
#include "CartesianR2ZoniShiftedShafranov.h"
#include "CartesianR2ZoniShiftedTriangular.h"
#include "PolarR6ZoniShiftedCircular.h"
#include "PolarR6ZoniShiftedShafranov.h"
#include "PolarR6ZoniShiftedTriangular.h"
#include "CartesianR6ZoniShiftedCircular.h"
#include "CartesianR6ZoniShiftedShafranov.h"
#include "CartesianR6ZoniShiftedTriangular.h"
#include "CartesianR2PoissonCircular.h"
#include "CartesianR2PoissonShafranov.h"
#include "CartesianR2PoissonTriangular.h"
#include "PolarR6PoissonCircular.h"
#include "PolarR6PoissonShafranov.h"
#include "PolarR6PoissonTriangular.h"
#include "CartesianR6PoissonCircular.h"
#include "CartesianR6PoissonShafranov.h"
#include "CartesianR6PoissonTriangular.h"
#include "RefinedGyroZoniShiftedCircular.h"
#include "RefinedGyroZoniShiftedShafranov.h"
#include "RefinedGyroZoniShiftedTriangular.h"
#include "RefinedGyroZoniShiftedCulham.h"

namespace gyro
{
/*******************************************************************************
 * Attributes
 ******************************************************************************/
/***************************************************************************
     * Controls and Informations
     **************************************************************************/
/* controls initialization */
std::vector<int> icntl(40, 0);
std::vector<double> dcntl(30, 0);

/* infos initialization */
std::vector<int> info(10, 0);
std::vector<double> dinfo(10, 0);
std::unique_ptr<ExactFuncs> functions;

std::string f_grid_r     = "";
std::string f_grid_theta = "";
std::string f_sol_in     = "";
std::string f_sol_out    = "";
} // namespace gyro

/*******************************************************************************
 * Methods
 ******************************************************************************/
/***************************************************************************
     * Parameters
     **************************************************************************/
/*!
 *  \brief Initialize the chosen parameters
 *
 *  Initialize the chosen parameters
 *
 */
void gyro::init_params()
{
    icntl[Param::optimized]          = 1;
    icntl[Param::matrix_free]        = 1;
    icntl[Param::debug]              = 0;
    icntl[Param::nr_exp]             = 3;
    icntl[Param::ntheta_exp]         = 3;
    icntl[Param::fac_ani]            = 0; // 0 or 1
    icntl[Param::v1]                 = 1;
    icntl[Param::v2]                 = 1;
    icntl[Param::cycle]              = 1; // 1 or 2
    icntl[Param::mod_pk]             = 0; // 0 or 1
    icntl[Param::compute_rho]        = 0; // 0 or 1
    icntl[Param::level]              = -1;
    icntl[Param::plotit]             = 1;
    icntl[Param::solveit]            = 1;
    icntl[Param::maxiter]            = 150;
    icntl[Param::periodic]           = 1;
    icntl[Param::origin_NOT_coarse]  = 0;
    icntl[Param::theta_aniso]        = 0;
    icntl[Param::smoother]           = 3; // originally 4 intended (optimized)
    icntl[Param::discr]              = 3;
    icntl[Param::extrapolation]      = 2;
    icntl[Param::DirBC_Interior]     = 0; // 0 or 1
    icntl[Param::paraview]           = 0;
    icntl[Param::divideBy2]          = 0;
    icntl[Param::prob]               = 5;
    icntl[Param::alpha_coeff]        = 0;
    icntl[Param::beta_coeff]         = 0;
    icntl[Param::verbose]            = 2;
    icntl[Param::openmp]             = 1;
    icntl[Param::res_norm]           = 3;
    icntl[Param::write_radii_angles] = 0;
    icntl[Param::check_error]        = 0;

    dcntl[Param::r0_DB]           = -1e6;
    dcntl[Param::R0]              = 0.1;
    dcntl[Param::R]               = 1.3;
    dcntl[Param::THETA0]          = 0.1;
    dcntl[Param::THETA]           = 1.3;
    dcntl[Param::kappa_eps]       = 0;
    dcntl[Param::delta_e]         = 0;
    dcntl[Param::tol_bound_check] = 1e-8;
    dcntl[Param::rel_red_conv]    = 1e-8;
    dcntl[Param::t_coeff]         = 0;
    dcntl[Param::t_arr_art_att]   = 0;
    dcntl[Param::t_sol]           = 0;
    dcntl[Param::t_detDFinv]      = 0;
    dcntl[Param::t_trafo]         = 0;
} /* ----- end of gyro::init_params ----- */

/*!
 *  \brief Shows the chosen parameters
 *
 *  Shows the chosen parameters
 *
 */
void gyro::show_params()
{
    std::cout << "optimized: " << icntl[Param::optimized] << "\n";
    std::cout << "matrix_free: " << icntl[Param::matrix_free] << "\n";
    std::cout << "debug: " << icntl[Param::debug] << "\n";
    std::cout << "nr_exp: " << icntl[Param::nr_exp] << "\n";
    std::cout << "ntheta_exp: " << icntl[Param::ntheta_exp] << "\n";
    std::cout << "fac_ani: " << icntl[Param::fac_ani] << "\n";
    std::cout << "theta_aniso: " << icntl[Param::theta_aniso] << "\n";
    std::cout << "v1: " << icntl[Param::v1] << "\n";
    std::cout << "v2: " << icntl[Param::v2] << "\n";
    std::cout << "cycle: " << icntl[Param::cycle] << "\n";
    std::cout << "mod_pk: " << icntl[Param::mod_pk] << "\n";
    std::cout << "compute_rho: " << icntl[Param::compute_rho] << "\n";
    std::cout << "level: " << icntl[Param::level] << "\n";
    std::cout << "plotit: " << icntl[Param::plotit] << "\n";
    std::cout << "solveit: " << icntl[Param::solveit] << "\n";
    std::cout << "maxiter: " << icntl[Param::maxiter] << "\n";
    std::cout << "periodic: " << icntl[Param::periodic] << "\n";
    std::cout << "origin_NOT_coarse: " << icntl[Param::origin_NOT_coarse] << "\n";
    std::cout << "smoother: " << icntl[Param::smoother] << "\n";
    std::cout << "discr: " << icntl[Param::discr] << "\n";
    std::cout << "extrapolation: " << icntl[Param::extrapolation] << "\n";
    std::cout << "DirBC_Interior: " << icntl[Param::DirBC_Interior] << "\n";
    std::cout << "paraview: " << icntl[Param::paraview] << "\n";
    std::cout << "divideBy2: " << icntl[Param::divideBy2] << "\n";
    std::cout << "prob: " << icntl[Param::prob] << "\n";
    std::cout << "alpha_coeff: " << icntl[Param::alpha_coeff] << "\n";
    std::cout << "beta_coeff: " << icntl[Param::beta_coeff] << "\n";
    std::cout << "verbose: " << icntl[Param::verbose] << "\n";
    std::cout << "openmp: " << icntl[Param::openmp] << "\n";
    std::cout << "res_norm: " << icntl[Param::res_norm] << "\n";
    std::cout << "write_radii_angles: " << icntl[Param::write_radii_angles] << "\n";
    std::cout << "check_error: " << icntl[Param::check_error] << "\n";

    std::cout << "R0: " << dcntl[Param::R0] << "\n";
    std::cout << "R: " << dcntl[Param::R] << "\n";
    std::cout << "kappa_eps: " << dcntl[Param::kappa_eps] << "\n";
    std::cout << "delta_e: " << dcntl[Param::delta_e] << "\n";
    std::cout << "tol_bound_check: " << dcntl[Param::tol_bound_check] << "\n";
    std::cout << "rel_red_conv: " << dcntl[Param::rel_red_conv] << "\n";
    std::cout << "\n";
} /* ----- end of gyro::show_params ----- */

/***************************************************************************
     * Boundary and solution
     **************************************************************************/
/************************
 * Single
 ************************/
/*!
 *  \brief Computes the distance of a node (r, theta) to the Dirichlet boundary.
 * The function assumes a boundary given by r_0=Param::r0_DB and R=Param::R.
 *
 * Attention: As we compute (r-r_0)*(r-R), the output of this function is highly dependent on r_0.
 * If r_0 is set "very close" to zero, roundoff errors can appear for function values on the innermost circles.
 * However, as dcntl[Param::r0_DB] is advised to be chosen small if icntl[Param::DirBC_Interior]=0 is set,
 * this function should be handled with care. Modestly small values like 1e-5 or 1e-8 should not create a problem.
 * 
 * For more details on the Across-The-Origin heuristic, which is implemented with icntl[Param::DirBC_Interior]=0,
 * see Kühn, Kruse, Rüde, Journal of Scientific Computing, 91 (28), (2022).
 * 
 * \param r_i: the r coordinate of the node
 * \param theta_j: the theta coordinate of the node
 * \param verbose: verbose level for debug
 *
 * \return the distance
 *
 */
double gyro::distBoundary(double r_i, double theta_j, int verbose)
{
    // // alt: rectangle/periodic
    // double r_i, theta_j;
    // gyro::trafo_back(r_i, theta_j, x, y, 0);
    // x = r_i;
    // y = theta_j;

    double boundDefDim1 = fabs(r_i - dcntl[Param::r0_DB]) * fabs(r_i - dcntl[Param::R]);
    double boundDefDim2 = 1;

    if (verbose > 5) {
        std::cout << "DISTBOUNDARY (" << r_i << ", " << theta_j << "): " << boundDefDim1 * boundDefDim2
                  << " (boundDefDim1: " << boundDefDim1 << ", boundDefDim2: " << boundDefDim2 << ")\n";
    }
    return boundDefDim1 * boundDefDim2;
} /* ----- end of gyro::distBoundary ----- */

/*!
 *  \brief Sign function
 *
 *  Sign function returns -1 (< 0) or 1 (>= 0)
 * 
 * \param n: number which sign is to compute
 *
 * \return the sign
 *
 */
int gyro::sign(double n)
{
    return (n < 0) ? -1 : 1;
}

/*!
 *  \brief Evaluates the solution
 *
 *  Evaluates the solution from the polar coordinates
 * 
 * \param r_i: the r coordinate of the node
 * \param theta_j: the theta coordinate of the node
 * \param verbose: verbose level for debug
 *
 * \return the solution value
 *
 */
double gyro::def_solution_rt(double r_i, double theta_j, int verbose)
{
    double t;
    TIC;

    double kappa_eps = dcntl[Param::kappa_eps];
    double delta_e   = dcntl[Param::delta_e];
    double Rmax      = dcntl[Param::R];

    double sol = functions->phi_exact(r_i, theta_j, kappa_eps, delta_e, Rmax);

    dcntl[Param::t_sol] += TOC;

    if (verbose > 5) {
        std::cout << "SOL (" << r_i << ", " << theta_j << "): " << sol << "\n";
    }
    return sol;
} /* ----- end of gyro::eval_def_solution_vec ----- */

/************************
 * Vector
 ************************/
/*!
 *  \brief Evaluates the solution
 *
 *  Evaluates the solution from the r coordinate on all theta positions
 * 
 * \param r_i: the r coordinate of the node
 * \param theta: vector theta (0, ntheta_int)
 * \param sin_theta: sines of theta
 * \param cos_theta: cosines of theta
 * \param ntheta: number of values in theta
 * \param verbose: verbose level for debug
 *
 * \return the solution vector
 */
std::vector<double> gyro::def_solution_rt(double r_i, std::vector<double> theta, std::vector<double> sin_theta,
                                          std::vector<double> cos_theta, int ntheta, int verbose)
{
    double t;
    TIC;

    double kappa_eps = dcntl[Param::kappa_eps];
    double delta_e   = dcntl[Param::delta_e];
    double Rmax      = dcntl[Param::R];
    // double R0        = dcntl[Param::R0];

    std::vector<double> sol(ntheta);
    functions->phi_exact(r_i, theta, kappa_eps, delta_e, Rmax, sol, sin_theta, cos_theta);

    dcntl[Param::t_sol] += TOC;

    if (verbose > 5) {
        for (int i = 0; i < ntheta; i++)
            std::cout << "SOL (" << r_i << ", " << theta[i] << "): " << sol[i] << "\n";
    }
    return sol;
} /* ----- end of gyro::eval_def_solution_vec ----- */

/***************************************************************************
     * Diffusivity and operator
     **************************************************************************/
/************************
 * Single
 ************************/
/*!
 *  \brief Diffusivity coefficient
 *
 * Computes the diffusivity coefficient depending on the radius r
 * - for prob 3 or 5: alpha(r)=(2.0/(2.6+3.14)) * (1.3 + atan((1-1.3*r/Rmax) / 0.09))
 * - else: 1
 * 
 * \param r: the r coordinate
 * \param verbose: verbose level for debug
 * 
 * \return the coefficient
 *
 */
double gyro::coeff(double r, int verbose)
{
    double t;
    TIC;

    double Rmax = dcntl[Param::R];

    double coeff_a = functions->coeffs1(r, Rmax);

    dcntl[Param::t_coeff] += TOC;

    if (verbose > 5) {
        std::cout << "COEFF (" << r << "): " << coeff << "\n";
    }
    return coeff_a;
} /* ----- end of level::coeff ----- */

/*!
 *  \brief Diffusivity coefficient
 *
 * Computes the diffusivity coefficient depending on the radius r
 * - for prob 3 or 5: alpha(r)=(2.0/(2.6+3.14)) * (1.3 + atan((1-1.3*r/Rmax) / 0.09))
 * - else: 1
 * 
 * \param r: the r coordinate
 * \param verbose: verbose level for debug
 * 
 * \return the coefficient
 *
 */
double gyro::coeff_beta(double r, int verbose)
{
    double t;
    TIC;

    double Rmax = dcntl[Param::R];

    double coeff_b = functions->coeffs2(r, Rmax);

    dcntl[Param::t_coeff] += TOC;

    return coeff_b;
} /* ----- end of level::coeff ----- */

/*!
 *  \brief Diffusivity coefficient
 *
 * Computes the diffusivity coefficient depending on the radius r
 * - for prob 3 or 5: alpha(r)=(2.0/(2.6+3.14)) * (1.3 + atan((1-1.3*r/Rmax) / 0.09))
 * - else: 1
 * 
 * \param r: the r coordinate
 * \param verbose: verbose level for debug
 * 
 * \return the coefficient
 *
 */
std::vector<double> gyro::coeff(std::vector<double> r, int verbose)
{
    double t;
    TIC;

    double Rmax = dcntl[Param::R];

    std::vector<double> coeff_a(r.size());
    functions->coeffs1(r, Rmax, coeff_a);

    dcntl[Param::t_coeff] += TOC;

    if (gyro::icntl[Param::verbose] > 5) {
        disp(r, "r");
        disp(coeff_a, "coeff_a");
    }
    return coeff_a;
} /* ----- end of level::coeff ----- */

/*!
 *  \brief Diffusivity coefficient
 *
 * Computes the diffusivity coefficient depending on the radius r
 * - for prob 3 or 5: alpha(r)=(2.0/(2.6+3.14)) * (1.3 + atan((1-1.3*r/Rmax) / 0.09))
 * - else: 1
 * 
 * \param r: the r coordinate
 * \param verbose: verbose level for debug
 * 
 * \return the coefficient
 *
 */
std::vector<double> gyro::coeff_beta(std::vector<double> r, int verbose)
{
    double t;
    TIC;

    double Rmax = dcntl[Param::R];

    std::vector<double> coeff_b(r.size());
    functions->coeffs2(r, Rmax, coeff_b);

    dcntl[Param::t_coeff] += TOC;
    if (gyro::icntl[Param::verbose] > 5) {
        disp(r, "r");
        disp(coeff_b, "coeff_b");
    }
    return coeff_b;
} /* ----- end of level::coeff_beta ----- */

/*!
 *  \brief detDFinv/arr/art/att/arr_att_art (single)
 *
 * Computes values from the polar coordinates
 * - computing the determinant of the derivative of F_inverse (det(DFinv))
 * - computing the coefficient a_r,r
 * - computing the coefficient a_r,theta
 * - computing the coefficient a_theta,theta
 * - computing all 3 coefficients at once to reduce flops
 * 
 * \param r: the r coordinate of the node
 * \param theta: the theta coordinate of the node
 * \param verbose: verbose level for debug
 * 
 * \return the computed value
 *
 */
double gyro::detDFinv(double r, double theta, int verbose)
{
    double t;
    TIC;

    double kappa_eps = dcntl[Param::kappa_eps];
    double delta_e   = dcntl[Param::delta_e];
    double Rmax      = dcntl[Param::R];

    double Jrr        = functions->J_rr(r, theta, kappa_eps, delta_e, Rmax);
    double Jrt        = functions->J_rt(r, theta, kappa_eps, delta_e, Rmax);
    double Jtr        = functions->J_tr(r, theta, kappa_eps, delta_e, Rmax);
    double Jtt        = functions->J_tt(r, theta, kappa_eps, delta_e, Rmax);
    double detDFinv_r = Jrr * Jtt - Jrt * Jtr;

    dcntl[Param::t_detDFinv] += TOC;

    if (gyro::icntl[Param::verbose] > 5) {
        std::cout << "Value of detDFinv (" << r << ", " << theta << "): " << detDFinv_r << "\n";
    }
    return detDFinv_r;
} /* ----- end of level::detDFinv ----- */
double gyro::arr(double r, double theta, int verbose)
{
    double t;
    TIC;

    double kappa_eps = dcntl[Param::kappa_eps];
    double delta_e   = dcntl[Param::delta_e];
    double Rmax      = dcntl[Param::R];

    double detDFinv_r = detDFinv(r, theta, verbose);
    double coeff_r    = coeff(r, verbose);
    double Jrt        = functions->J_rt(r, theta, kappa_eps, delta_e, Rmax);
    double Jtt        = functions->J_tt(r, theta, kappa_eps, delta_e, Rmax);
    double arr_r      = 0.5 * (Jtt * Jtt + Jrt * Jrt) * coeff_r / fabs(detDFinv_r);

    dcntl[Param::t_arr_art_att] += TOC;

    if (gyro::icntl[Param::verbose] > 5) {
        std::cout << "Value of arr (" << r << ", " << theta << "): " << arr_r << "\n";
    }
    return arr_r;
} /* ----- end of level::arr ----- */
double gyro::art(double r, double theta, int verbose)
{
    double t;
    TIC;

    double kappa_eps = dcntl[Param::kappa_eps];
    double delta_e   = dcntl[Param::delta_e];
    double Rmax      = dcntl[Param::R];

    double detDFinv_r = detDFinv(r, theta, verbose);
    double coeff_r    = coeff(r, verbose);
    double Jrr        = functions->J_rr(r, theta, kappa_eps, delta_e, Rmax);
    double Jrt        = functions->J_rt(r, theta, kappa_eps, delta_e, Rmax);
    double Jtr        = functions->J_tr(r, theta, kappa_eps, delta_e, Rmax);
    double Jtt        = functions->J_tt(r, theta, kappa_eps, delta_e, Rmax);
    double art_r      = -0.25 * (Jtt * Jtr + Jrt * Jrr) * coeff_r / fabs(detDFinv_r);

    dcntl[Param::t_arr_art_att] += TOC;

    if (gyro::icntl[Param::verbose] > 5) {
        std::cout << "Value of art (" << r << ", " << theta << "): " << art_r << "\n";
    }
    return art_r;
} /* ----- end of level::art ----- */
double gyro::att(double r, double theta, int verbose)
{
    double t;
    TIC;

    double kappa_eps = dcntl[Param::kappa_eps];
    double delta_e   = dcntl[Param::delta_e];
    double Rmax      = dcntl[Param::R];

    double detDFinv_r = detDFinv(r, theta, verbose);
    double coeff_r    = coeff(r, verbose);
    double Jrr        = functions->J_rr(r, theta, kappa_eps, delta_e, Rmax);
    double Jtr        = functions->J_tr(r, theta, kappa_eps, delta_e, Rmax);
    double att_r      = 0.5 * (Jtr * Jtr + Jrr * Jrr) * coeff_r / fabs(detDFinv_r);

    dcntl[Param::t_arr_art_att] += TOC;

    if (gyro::icntl[Param::verbose] > 5) {
        std::cout << "Value of att (" << r << ", " << theta << "): " << att_r << "\n";
    }
    return att_r;
} /* ----- end of level::att ----- */
void gyro::arr_att_art(double r, double theta, double& arr, double& att, double& art, int verbose)
{
    double t;
    TIC;

    double kappa_eps = dcntl[Param::kappa_eps];
    double delta_e   = dcntl[Param::delta_e];
    double Rmax      = dcntl[Param::R];

    double detDFinv_r = detDFinv(r, theta, verbose);
    double coeff_r    = coeff(r, verbose);
    double Jrr        = functions->J_rr(r, theta, kappa_eps, delta_e, Rmax);
    double Jrt        = functions->J_rt(r, theta, kappa_eps, delta_e, Rmax);
    double Jtr        = functions->J_tr(r, theta, kappa_eps, delta_e, Rmax);
    double Jtt        = functions->J_tt(r, theta, kappa_eps, delta_e, Rmax);
    double coeff      = coeff_r / fabs(detDFinv_r);
    arr               = 0.5 * (Jtt * Jtt + Jrt * Jrt) * coeff;
    art               = -0.25 * (Jtt * Jtr + Jrt * Jrr) * coeff;
    att               = 0.5 * (Jtr * Jtr + Jrr * Jrr) * coeff;

    dcntl[Param::t_arr_art_att] += TOC;

    if (gyro::icntl[Param::verbose] > 5) {
        std::cout << "Value of arr (" << r << ", " << theta << "): " << arr << "\n";
        std::cout << "Value of att (" << r << ", " << theta << "): " << att << "\n";
        std::cout << "Value of art (" << r << ", " << theta << "): " << art << "\n";
    }
} /* ----- end of level::arr_att_art ----- */

/*!
 *  \brief detDFinv/arr/art/att/arr_att_art (vector)
 *
 * Computes values from the r coordinate on all theta positions
 * - computing the determinant of the derivative of F_inverse (det(DFinv))
 * - computing the coefficient a_r,r
 * - computing the coefficient a_r,theta
 * - computing the coefficient a_theta,theta
 * - computing all 3 coefficients at once to reduce flops
 * 
 * \param r: the r coordinate of the node
 * \param theta: the theta coordinate of the node
 * \param verbose: verbose level for debug
 * 
 * \return the computed vector
 *
 */
std::vector<double> gyro::detDFinv(double r, std::vector<double> theta, std::vector<double> sin_theta,
                                   std::vector<double> cos_theta, int ntheta, int verbose)
{
    double t;
    TIC;

    double kappa_eps = dcntl[Param::kappa_eps];
    double delta_e   = dcntl[Param::delta_e];
    double Rmax      = dcntl[Param::R];

    std::vector<double> detDFinv_r(ntheta);
    std::vector<double> Jrr(ntheta);
    std::vector<double> Jrt(ntheta);
    std::vector<double> Jtr(ntheta);
    std::vector<double> Jtt(ntheta);
    functions->J_rr(r, theta, kappa_eps, delta_e, Rmax, Jrr, sin_theta, cos_theta);
    functions->J_rt(r, theta, kappa_eps, delta_e, Rmax, Jrt, sin_theta, cos_theta);
    functions->J_tr(r, theta, kappa_eps, delta_e, Rmax, Jtr, sin_theta, cos_theta);
    functions->J_tt(r, theta, kappa_eps, delta_e, Rmax, Jtt, sin_theta, cos_theta);
    for (int i = 0; i < ntheta; i++) {
        detDFinv_r[i] = Jrr[i] * Jtt[i] - Jrt[i] * Jtr[i];
    }

    dcntl[Param::t_detDFinv] += TOC;

    if (gyro::icntl[Param::verbose] > 5) {
        disp(theta, "theta");
        disp(cos_theta, "cos_theta");
        disp(detDFinv_r, "detDFinv_r");
    }
    return detDFinv_r;
} /* ----- end of level::detDFinv ----- */

/*!
 *  \brief detDFinv/arr/art/att/arr_att_art (vector)
 *
 * Computes values from the r coordinate on all theta positions
 * - computing the determinant of the derivative of F_inverse (det(DFinv))
 * - computing the coefficient a_r,r
 * - computing the coefficient a_r,theta
 * - computing the coefficient a_theta,theta
 * - computing all 3 coefficients at once to reduce flops
 * 
 * \param r: the r coordinate of the node
 * \param theta: the theta coordinate of the node
 * \param verbose: verbose level for debug
 * 
 * \return the computed vector
 *
 */
std::vector<double> gyro::detDFinv(std::vector<double> r, double theta, int verbose)
{
    double t;
    TIC;

    double kappa_eps = dcntl[Param::kappa_eps];
    double delta_e   = dcntl[Param::delta_e];
    double Rmax      = dcntl[Param::R];

    std::size_t nr = r.size();
    std::vector<double> detDFinv_r(nr);
    std::vector<double> Jrr(nr);
    std::vector<double> Jrt(nr);
    std::vector<double> Jtr(nr);
    std::vector<double> Jtt(nr);
    functions->J_rr(r, theta, kappa_eps, delta_e, Rmax, Jrr);
    functions->J_rt(r, theta, kappa_eps, delta_e, Rmax, Jrt);
    functions->J_tr(r, theta, kappa_eps, delta_e, Rmax, Jtr);
    functions->J_tt(r, theta, kappa_eps, delta_e, Rmax, Jtt);
    for (std::size_t i = 0; i < r.size(); i++) {
        detDFinv_r[i] = Jrr[i] * Jtt[i] - Jrt[i] * Jtr[i];
    }

    dcntl[Param::t_detDFinv] += TOC;

    if (gyro::icntl[Param::verbose] > 5) {
        disp(r, "r");
        disp(detDFinv_r, "detDFinv_r");
    }
    return detDFinv_r;
} /* ----- end of level::detDFinv ----- */

std::vector<double> gyro::arr(double r, std::vector<double> theta, std::vector<double> sin_theta,
                              std::vector<double> cos_theta, int ntheta, int verbose)
{
    double t;
    TIC;

    double kappa_eps = dcntl[Param::kappa_eps];
    double delta_e   = dcntl[Param::delta_e];
    double Rmax      = dcntl[Param::R];

    std::vector<double> arr_r(ntheta);
    std::vector<double> detDFinv_r = detDFinv(r, theta, sin_theta, cos_theta, ntheta, verbose);
    std::vector<double> Jrt(ntheta);
    std::vector<double> Jtt(ntheta);
    functions->J_rt(r, theta, kappa_eps, delta_e, Rmax, Jrt, sin_theta, cos_theta);
    functions->J_tt(r, theta, kappa_eps, delta_e, Rmax, Jtt, sin_theta, cos_theta);
    double coeff_r = coeff(r, verbose);
    ;
    for (int i = 0; i < ntheta; i++) {
        arr_r[i] = 0.5 * (Jtt[i] * Jtt[i] + Jrt[i] * Jrt[i]) * coeff_r / fabs(detDFinv_r[i]);
    }

    dcntl[Param::t_arr_art_att] += TOC;

    if (gyro::icntl[Param::verbose] > 5) {
        disp(theta, "theta");
        disp(sin_theta, "sin_theta");
        disp(cos_theta, "cos_theta");
        disp(arr_r, "arr_r");
    }
    return arr_r;
} /* ----- end of level::arr ----- */
std::vector<double> gyro::art(double r, std::vector<double> theta, std::vector<double> sin_theta,
                              std::vector<double> cos_theta, int ntheta, int verbose)
{
    double t;
    TIC;

    double kappa_eps = dcntl[Param::kappa_eps];
    double delta_e   = dcntl[Param::delta_e];
    double Rmax      = dcntl[Param::R];

    std::vector<double> art_r(ntheta);
    std::vector<double> detDFinv_r = detDFinv(r, theta, sin_theta, cos_theta, ntheta, verbose);
    std::vector<double> Jrr(ntheta);
    std::vector<double> Jrt(ntheta);
    std::vector<double> Jtr(ntheta);
    std::vector<double> Jtt(ntheta);
    functions->J_rr(r, theta, kappa_eps, delta_e, Rmax, Jrr, sin_theta, cos_theta);
    functions->J_rt(r, theta, kappa_eps, delta_e, Rmax, Jrt, sin_theta, cos_theta);
    functions->J_tr(r, theta, kappa_eps, delta_e, Rmax, Jtr, sin_theta, cos_theta);
    functions->J_tt(r, theta, kappa_eps, delta_e, Rmax, Jtt, sin_theta, cos_theta);
    double coeff_r = coeff(r, verbose);
    ;
    for (int i = 0; i < ntheta; i++) {
        art_r[i] = -0.25 * (Jtt[i] * Jtr[i] + Jrt[i] * Jrr[i]) * coeff_r / fabs(detDFinv_r[i]);
    }

    dcntl[Param::t_arr_art_att] += TOC;

    if (gyro::icntl[Param::verbose] > 5) {
        disp(theta, "theta");
        disp(sin_theta, "sin_theta");
        disp(cos_theta, "cos_theta");
        disp(art_r, "art_r");
    }
    return art_r;
} /* ----- end of level::art ----- */
std::vector<double> gyro::att(double r, std::vector<double> theta, std::vector<double> sin_theta,
                              std::vector<double> cos_theta, int ntheta, int verbose)
{
    double t;
    TIC;

    double kappa_eps = dcntl[Param::kappa_eps];
    double delta_e   = dcntl[Param::delta_e];
    double Rmax      = dcntl[Param::R];

    std::vector<double> att_r(ntheta);
    std::vector<double> detDFinv_r = detDFinv(r, theta, sin_theta, cos_theta, ntheta, verbose);
    std::vector<double> Jrr(ntheta);
    std::vector<double> Jrt(ntheta);
    std::vector<double> Jtr(ntheta);
    std::vector<double> Jtt(ntheta);
    functions->J_rr(r, theta, kappa_eps, delta_e, Rmax, Jrr, sin_theta, cos_theta);
    functions->J_tr(r, theta, kappa_eps, delta_e, Rmax, Jtr, sin_theta, cos_theta);
    double coeff_r = coeff(r, verbose);
    ;
    for (int i = 0; i < ntheta; i++) {
        att_r[i] = 0.5 * (Jtr[i] * Jtr[i] + Jrr[i] * Jrr[i]) * coeff_r / fabs(detDFinv_r[i]);
    }

    dcntl[Param::t_arr_art_att] += TOC;

    if (gyro::icntl[Param::verbose] > 5) {
        disp(theta, "theta");
        disp(sin_theta, "sin_theta");
        disp(cos_theta, "cos_theta");
        disp(att_r, "att_r");
    }
    return att_r;
} /* ----- end of level::att ----- */
void gyro::arr_att_art(double r, std::vector<double> theta, std::vector<double> sin_theta,
                       std::vector<double> cos_theta, int ntheta, std::vector<double>& arr, std::vector<double>& att,
                       std::vector<double>& art, int verbose)
{
    double t;
    TIC;

    double kappa_eps = dcntl[Param::kappa_eps];
    double delta_e   = dcntl[Param::delta_e];
    double Rmax      = dcntl[Param::R];

    arr                            = std::vector<double>(ntheta);
    att                            = std::vector<double>(ntheta);
    art                            = std::vector<double>(ntheta);
    std::vector<double> detDFinv_r = detDFinv(r, theta, sin_theta, cos_theta, ntheta, verbose);
    std::vector<double> Jrr(ntheta);
    std::vector<double> Jrt(ntheta);
    std::vector<double> Jtr(ntheta);
    std::vector<double> Jtt(ntheta);
    functions->J_rr(r, theta, kappa_eps, delta_e, Rmax, Jrr, sin_theta, cos_theta);
    functions->J_rt(r, theta, kappa_eps, delta_e, Rmax, Jrt, sin_theta, cos_theta);
    functions->J_tr(r, theta, kappa_eps, delta_e, Rmax, Jtr, sin_theta, cos_theta);
    functions->J_tt(r, theta, kappa_eps, delta_e, Rmax, Jtt, sin_theta, cos_theta);
    double coeff_r = coeff(r, verbose);
    for (int i = 0; i < ntheta; i++) {
        double coeff = coeff_r / fabs(detDFinv_r[i]);
        arr[i]       = 0.5 * (Jtt[i] * Jtt[i] + Jrt[i] * Jrt[i]) * coeff;
        art[i]       = -0.25 * (Jtt[i] * Jtr[i] + Jrt[i] * Jrr[i]) * coeff;
        att[i]       = 0.5 * (Jtr[i] * Jtr[i] + Jrr[i] * Jrr[i]) * coeff;
    }

    dcntl[Param::t_arr_art_att] += TOC;

    if (gyro::icntl[Param::verbose] > 5) {
        disp(theta, "theta");
        disp(sin_theta, "sin_theta");
        disp(cos_theta, "cos_theta");
        disp(arr, "arr");
        disp(att, "att");
        disp(art, "art");
    }
} /* ----- end of level::arr_att_art ----- */
void gyro::arr_att_art(std::vector<double> r, double theta, std::vector<double>& arr_r, std::vector<double>& att_r,
                       std::vector<double>& art_r, int verbose)
{
    double t;
    TIC;

    double kappa_eps = dcntl[Param::kappa_eps];
    double delta_e   = dcntl[Param::delta_e];
    double Rmax      = dcntl[Param::R];

    int size                       = r.size();
    arr_r                          = std::vector<double>(size);
    att_r                          = std::vector<double>(size);
    art_r                          = std::vector<double>(size);
    std::vector<double> detDFinv_r = detDFinv(r, theta, verbose);
    std::vector<double> Jrr(size);
    std::vector<double> Jrt(size);
    std::vector<double> Jtr(size);
    std::vector<double> Jtt(size);
    std::vector<double> coeff_r = coeff(r, verbose);
    functions->J_rr(r, theta, kappa_eps, delta_e, Rmax, Jrr);
    functions->J_rt(r, theta, kappa_eps, delta_e, Rmax, Jrt);
    functions->J_tr(r, theta, kappa_eps, delta_e, Rmax, Jtr);
    functions->J_tt(r, theta, kappa_eps, delta_e, Rmax, Jtt);
    for (std::size_t j = 0; j < r.size(); j++) {
        double coeff = coeff_r[j] / fabs(detDFinv_r[j]);
        arr_r[j]     = 0.5 * (Jtt[j] * Jtt[j] + Jrt[j] * Jrt[j]) * coeff;
        art_r[j]     = -0.25 * (Jtt[j] * Jtr[j] + Jrt[j] * Jrr[j]) * coeff;
        att_r[j]     = 0.5 * (Jtr[j] * Jtr[j] + Jrr[j] * Jrr[j]) * coeff;
    }

    dcntl[Param::t_arr_art_att] += TOC;

    if (gyro::icntl[Param::verbose] > 5) {
        disp(r, "r");
        disp(detDFinv_r, "detDFinv_r");
        disp(coeff_r, "coeff_r");
        disp(arr_r, "arr_r");
        disp(att_r, "att_r");
        disp(art_r, "art_r");
    }
} /* ----- end of level::arr_att_art ----- */

/***************************************************************************
     * Polar to cartesian and back
     **************************************************************************/
/************************
 * Single
 ************************/
/*!
 *  \brief Transform from polar to cartesian
 *
 *  Transform one couple r_i,theta_j from polar to cartesian
 * 
 * \param r_i: the r coordinate of the node
 * \param theta_j: the theta coordinate of the node
 * \param x: the x coordinate of the node (out)
 * \param y: the y coordinate of the node (out)
 * \param verbose: verbose level for debug
 *
 */
void gyro::trafo(double& r_i, double& theta_j, double& x, double& y, int verbose)
{
    double t;
    TIC;

    double kappa_eps = gyro::dcntl[Param::kappa_eps];
    double delta_e   = gyro::dcntl[Param::delta_e];
    double Rmax      = dcntl[Param::R];

    x = functions->x(r_i, theta_j, kappa_eps, delta_e, Rmax);
    y = functions->y(r_i, theta_j, kappa_eps, delta_e, Rmax);

    dcntl[Param::t_trafo] += TOC;

    // if (verbose)
    //     std::cout << "TRAFO (" << r_i << ", " << theta_j << "): (" << x << ", " << y << ")\n";
} /* ----- end of gyro::trafo ----- */

/*!
 *  \brief Transform from polar to cartesian
 *
 *  TODO: This function could be reimplemented and used for simple geometries
 *        For more general geometries, the inverse mapping is not available 
 *        and a check on the geometry would be needed.
 *         
 *  Transform one couple r_i,theta_j from cartesian to polar
 * 
 * \param r_i: the r coordinate of the node (out)
 * \param theta_j: the theta coordinate of the node (out)
 * \param x: the x coordinate of the node
 * \param y: the y coordinate of the node
 * \param verbose: verbose level for debug
 *
 */
void gyro::trafo_back(double& r_i, double& theta_j, double& x, double& y, int verbose)
{
    double t;
    TIC;

    //double kappa_eps = gyro::dcntl[Param::kappa_eps];
    //double delta_e   = gyro::dcntl[Param::delta_e];

    throw std::runtime_error("No general inverse mapping available.");
    /*
    if (mod_pk == geometry::CIRCULAR) {
        r_i     = sqrt(x * x + y * y);
        theta_j = atan2(y, x);
    }
    else if (mod_pk == geometry::SHAFRANOV) {
        // Back transformer on (r,phi) of (x,y), see paper https://arxiv.org/pdf/1712.02201.pdf p.4
        r_i     = sqrt(2 * ((x * x / pow(1 - kappa_eps, 2) + y * y / pow(1 + kappa_eps, 2))) /
                       (1 - 2 * delta_e * x / pow(1 - kappa_eps, 2) +
                    sqrt(pow(1 - 2 * delta_e * x / pow(1 - kappa_eps, 2), 2) -
                             4 * pow(delta_e, 2) / pow(1 - kappa_eps, 2) *
                                 (x * x / pow(1 - kappa_eps, 2) + y * y / pow(1 + kappa_eps, 2)))));
        theta_j = atan2(y / (1 + kappa_eps), (x + delta_e * r_i * r_i) / (1 - kappa_eps));
    }
    else if (mod_pk == geometry::TRIANGULAR) {
        throw std::runtime_error("No inverse mapping available for the Czarny geometry.");
        r_i     = 0;
        theta_j = 0;
    }*/

    dcntl[Param::t_trafo] += TOC;

    // if (verbose)
    //     std::cout << "TRAFO_BACK (" << x << ", " << y << "): (" << r_i << ", " << theta_j << ")\n";
} /* ----- end of gyro::trafo_back ----- */

/************************
 * Vector
 ************************/
/*!
 *  \brief Transform from polar to cartesian
 *
 *  Transform one couple r_i,theta_j from polar to cartesian
 * 
 * \param r_i: the r coordinate of the node
 * \param theta: vector theta (0, ntheta_int)
 * \param sin_theta: sines of theta
 * \param cos_theta: cosines of theta
 * \param ntheta: number of values in theta
 * \param x: vector x
 * \param y: vector y
 * \param verbose: verbose level for debug
 *
 */
void gyro::trafo(double r_i, std::vector<double> theta, std::vector<double> sin_theta, std::vector<double> cos_theta,
                 int ntheta, std::vector<double>& x, std::vector<double>& y, int verbose)
{
    double t;
    TIC;

    double kappa_eps = gyro::dcntl[Param::kappa_eps];
    double delta_e   = gyro::dcntl[Param::delta_e];
    double Rmax      = dcntl[Param::R];
    x                = std::vector<double>(ntheta);
    y                = std::vector<double>(ntheta);

    functions->x(r_i, theta, kappa_eps, delta_e, Rmax, x, sin_theta, cos_theta);
    functions->y(r_i, theta, kappa_eps, delta_e, Rmax, y, sin_theta, cos_theta);

    dcntl[Param::t_trafo] += TOC;

    // if (verbose)
    //     for (int i = 0; i < ntheta; i++) {
    //         std::cout << "TRAFO (" << r_i << ", " << theta[i] << "): (" << x[i] << ", " << y[i] << ")\n";
    //     }
} /* ----- end of gyro::trafo ----- */

/*!
 *  \brief Sparse matrix-vector product
 *
 *  Sparse matrix-vector product:
 *  - (trans=0) y := alpha*A*x + beta*y
 *  - (trans=1) y := alpha*A^T*x + beta*y
 *
 */
void gyro::sp_dgemv(int trans, int m, int n, double alpha, std::vector<int> row_indices, std::vector<int> col_indices,
                    std::vector<double> vals, int lda, std::vector<double> x, int incx, double beta,
                    std::vector<double>& y, int incy)
{
    for (std::size_t i = 0; i < row_indices.size(); i++) {
        y[i] = beta * (y[i] + (double)incy);
    }
    for (std::size_t i = 0; i < row_indices.size(); i++) {
        y[row_indices[i]] += vals[i] * (x[col_indices[i]] + (double)incx);
    }
}

/*!
 *  \brief 2-norm of a vector
 *
 *  Computes the 2-norm of a vector
 *
 */
double gyro::norm(std::vector<double> x)
{
    double nrmres = 0;
    for (std::size_t i = 0; i < x.size(); i++) {
        nrmres += pow(x[i], 2);
    }
    return sqrt(nrmres);
}

/*!
 *  \brief 2-norm of a vector
 *
 *  Computes the 2-norm of a matrix
 *
 */
double gyro::A_norm(std::vector<double> x, int m, std::vector<int> row_indices, std::vector<int> col_indices,
                    std::vector<double> vals)
{
    double nrmres = 0;
    std::vector<double> Ax(m, 0.0);
    sp_dgemv(0, m, m, 1.0, row_indices, col_indices, vals, m, x, 0, 1.0, Ax, 0);
    for (std::size_t i = 0; i < x.size(); i++) {
        nrmres += pow(x[i] * Ax[i], 2);
    }
    return sqrt(nrmres);
}

void gyro::get_geometry_coeffs(geometry_type geom)
{
    switch (geom) {
    case CIRCULAR:
    case CULHAM:
        gyro::dcntl[Param::kappa_eps] = 0;
        gyro::dcntl[Param::delta_e]   = 0;
        break;
    case SHAFRANOV:
        gyro::dcntl[Param::kappa_eps] = 0.3;
        gyro::dcntl[Param::delta_e]   = 0.2;
        break;
    case TRIANGULAR:
        gyro::dcntl[Param::kappa_eps] = 0.3;
        gyro::dcntl[Param::delta_e]   = 1.4;
        break;
    }
}

void gyro::select_functions_class(int alpha_coeff, int beta_coeff, int geometry, int problem)
{
    if (alpha_coeff < 0 or alpha_coeff > 3)
        throw std::runtime_error("Unknown alpha coeff");
    if (beta_coeff < 0 or beta_coeff > 1)
        throw std::runtime_error("Unknown beta coeff");
    if (geometry < 0 or geometry > 3)
        throw std::runtime_error("Unknown geometry");
    if ((problem != 1 and problem < 3) or problem > 7)
        throw std::runtime_error("Unknown problem");
    if (problem == 1)
        throw std::runtime_error("Flat not implemented");

    alpha_val alpha    = (alpha_val)alpha_coeff;
    geometry_type geom = (geometry_type)geometry;
    problem_type prob  = (problem_type)problem;
    std::cout << "GEOMETRY: " << geom << "\n";
    if (beta_coeff) {
        switch (alpha) {
        case SONNENDRUCKER:
            switch (prob) {
            case CARTESIAN_R2:
                switch (geom) {
                case CIRCULAR:
                    gyro::functions = std::make_unique<CartesianR2GyroSonnendruckerCircular>();
                    break;
                case SHAFRANOV:
                    gyro::functions = std::make_unique<CartesianR2GyroSonnendruckerShafranov>();
                    break;
                case TRIANGULAR:
                    gyro::functions = std::make_unique<CartesianR2GyroSonnendruckerTriangular>();
                    break;
                default:
                    throw std::runtime_error("Wrong choice for the geometry\n");
                    break;
                }
                break;
            case POLAR_R6:
                switch (geom) {
                case CIRCULAR:
                    gyro::functions = std::make_unique<PolarR6GyroSonnendruckerCircular>();
                    break;
                case SHAFRANOV:
                    gyro::functions = std::make_unique<PolarR6GyroSonnendruckerShafranov>();
                    break;
                case TRIANGULAR:
                    gyro::functions = std::make_unique<PolarR6GyroSonnendruckerTriangular>();
                    break;
                default:
                    throw std::runtime_error("Wrong choice for the geometry\n");
                    break;
                }
                break;
            case CARTESIAN_R6:
                switch (geom) {
                case CIRCULAR:
                    gyro::functions = std::make_unique<CartesianR6GyroSonnendruckerCircular>();
                    break;
                case SHAFRANOV:
                    gyro::functions = std::make_unique<CartesianR6GyroSonnendruckerShafranov>();
                    break;
                case TRIANGULAR:
                    gyro::functions = std::make_unique<CartesianR6GyroSonnendruckerTriangular>();
                    break;
                default:
                    throw std::runtime_error("Wrong choice for the geometry\n");
                    break;
                }
                break;
                //case FLAT:
                //    throw std::runtime_error("Beta coeff cannot be 1/0");
                //    break;
            default:
                throw std::runtime_error("Wrong choice for the problem\n");
                break;
            }
            break;
        case ZONI:
            switch (prob) {
            case CARTESIAN_R2:
                switch (geom) {
                case CIRCULAR:
                    gyro::functions = std::make_unique<CartesianR2GyroZoniCircular>();
                    break;
                case SHAFRANOV:
                    gyro::functions = std::make_unique<CartesianR2GyroZoniShafranov>();
                    break;
                case TRIANGULAR:
                    gyro::functions = std::make_unique<CartesianR2GyroZoniTriangular>();
                    break;
                default:
                    throw std::runtime_error("Wrong choice for the geometry\n");
                    break;
                }
                break;
            case POLAR_R6:
                switch (geom) {
                case CIRCULAR:
                    gyro::functions = std::make_unique<PolarR6GyroZoniCircular>();
                    break;
                case SHAFRANOV:
                    gyro::functions = std::make_unique<PolarR6GyroZoniShafranov>();
                    break;
                case TRIANGULAR:
                    gyro::functions = std::make_unique<PolarR6GyroZoniTriangular>();
                    break;
                default:
                    throw std::runtime_error("Wrong choice for the geometry\n");
                    break;
                }
                break;
            case CARTESIAN_R6:
                switch (geom) {
                case CIRCULAR:
                    gyro::functions = std::make_unique<CartesianR6GyroZoniCircular>();
                    break;
                case SHAFRANOV:
                    gyro::functions = std::make_unique<CartesianR6GyroZoniShafranov>();
                    break;
                case TRIANGULAR:
                    gyro::functions = std::make_unique<CartesianR6GyroZoniTriangular>();
                    break;
                default:
                    throw std::runtime_error("Wrong choice for the geometry\n");
                    break;
                }
                break;
                //case FLAT:
                //    throw std::runtime_error("Beta coeff cannot be 1/0");
                //    break;
            default:
                throw std::runtime_error("Wrong choice for the problem\n");
                break;
            }
            break;
        case ZONI_SHIFTED:
            switch (prob) {
            case CARTESIAN_R2:
                switch (geom) {
                case CIRCULAR:
                    gyro::functions = std::make_unique<CartesianR2GyroZoniShiftedCircular>();
                    break;
                case SHAFRANOV:
                    gyro::functions = std::make_unique<CartesianR2GyroZoniShiftedShafranov>();
                    break;
                case TRIANGULAR:
                    gyro::functions = std::make_unique<CartesianR2GyroZoniShiftedTriangular>();
                    break;
                default:
                    throw std::runtime_error("Wrong choice for the geometry\n");
                    break;
                }
                break;
            case POLAR_R6:
                switch (geom) {
                case CIRCULAR:
                    gyro::functions = std::make_unique<PolarR6GyroZoniShiftedCircular>();
                    break;
                case SHAFRANOV:
                    gyro::functions = std::make_unique<PolarR6GyroZoniShiftedShafranov>();
                    break;
                case TRIANGULAR:
                    gyro::functions = std::make_unique<PolarR6GyroZoniShiftedTriangular>();
                    break;
                case CULHAM:
                    gyro::functions = std::make_unique<PolarR6GyroZoniShiftedCulham>();
                    break;
                default:
                    throw std::runtime_error("Wrong choice for the geometry\n");
                    break;
                }
                break;
            case CARTESIAN_R6:
                switch (geom) {
                case CIRCULAR:
                    gyro::functions = std::make_unique<CartesianR6GyroZoniShiftedCircular>();
                    break;
                case SHAFRANOV:
                    gyro::functions = std::make_unique<CartesianR6GyroZoniShiftedShafranov>();
                    break;
                case TRIANGULAR:
                    gyro::functions = std::make_unique<CartesianR6GyroZoniShiftedTriangular>();
                    break;
                default:
                    throw std::runtime_error("Wrong choice for the geometry\n");
                    break;
                }
                break;
            case REFINED_RADIUS:
                switch (geom) {
                case CIRCULAR:
                    gyro::functions = std::make_unique<RefinedGyroZoniShiftedCircular>();
                    break;
                case SHAFRANOV:
                    gyro::functions = std::make_unique<RefinedGyroZoniShiftedShafranov>();
                    break;
                case TRIANGULAR:
                    gyro::functions = std::make_unique<RefinedGyroZoniShiftedTriangular>();
                    break;
                case CULHAM:
                    gyro::functions = std::make_unique<RefinedGyroZoniShiftedCulham>();
                    break;
                default:
                    throw std::runtime_error("Wrong choice for the geometry\n");
                    break;
                }
                break;
                //case FLAT:
                //    throw std::runtime_error("Beta coeff cannot be 1/0");
                //    break;
            default:
                throw std::runtime_error("Wrong choice for the problem\n");
                break;
            }
            break;
        case POISSON:
            throw std::runtime_error("Beta coeff cannot be 1/0");
            break;
        default:
            throw std::runtime_error("Wrong choice for the alpha coefficient\n");
            break;
        }
    }
    else {
        switch (alpha) {
        case SONNENDRUCKER:
            switch (prob) {
            case CARTESIAN_R2:
                switch (geom) {
                case CIRCULAR:
                    gyro::functions = std::make_unique<CartesianR2SonnendruckerCircular>();
                    break;
                case SHAFRANOV:
                    gyro::functions = std::make_unique<CartesianR2SonnendruckerShafranov>();
                    break;
                case TRIANGULAR:
                    gyro::functions = std::make_unique<CartesianR2SonnendruckerTriangular>();
                    break;
                default:
                    throw std::runtime_error("Wrong choice for the geometry\n");
                    break;
                }
                break;
            case POLAR_R6:
                switch (geom) {
                case CIRCULAR:
                    gyro::functions = std::make_unique<PolarR6SonnendruckerCircular>();
                    break;
                case SHAFRANOV:
                    gyro::functions = std::make_unique<PolarR6SonnendruckerShafranov>();
                    break;
                case TRIANGULAR:
                    gyro::functions = std::make_unique<PolarR6SonnendruckerTriangular>();
                    break;
                default:
                    throw std::runtime_error("Wrong choice for the geometry\n");
                    break;
                }
                break;
            case CARTESIAN_R6:
                switch (geom) {
                case CIRCULAR:
                    gyro::functions = std::make_unique<CartesianR6SonnendruckerCircular>();
                    break;
                case SHAFRANOV:
                    gyro::functions = std::make_unique<CartesianR6SonnendruckerShafranov>();
                    break;
                case TRIANGULAR:
                    gyro::functions = std::make_unique<CartesianR6SonnendruckerTriangular>();
                    break;
                default:
                    throw std::runtime_error("Wrong choice for the geometry\n");
                    break;
                }
                break;
                //case FLAT:
                //    switch (geom) {
                //        case CIRCULAR: gyro::functions   = std::make_unique<FlatSonnendruckerCircular>(); break;
                //        case SHAFRANOV: gyro::functions  = std::make_unique<FlatSonnendruckerShafranov>(); break;
                //        case TRIANGULAR: gyro::functions = std::make_unique<FlatSonnendruckerTriangular>(); break;
                //    }
                //    break;
            default:
                throw std::runtime_error("Wrong choice for the problem\n");
                break;
            }
            break;
        case ZONI:
            switch (prob) {
            case CARTESIAN_R2:
                switch (geom) {
                case CIRCULAR:
                    gyro::functions = std::make_unique<CartesianR2ZoniCircular>();
                    break;
                case SHAFRANOV:
                    gyro::functions = std::make_unique<CartesianR2ZoniShafranov>();
                    break;
                case TRIANGULAR:
                    gyro::functions = std::make_unique<CartesianR2ZoniTriangular>();
                    break;
                default:
                    throw std::runtime_error("Wrong choice for the geometry\n");
                    break;
                }
                break;
            case POLAR_R6:
                switch (geom) {
                case CIRCULAR:
                    gyro::functions = std::make_unique<PolarR6ZoniCircular>();
                    break;
                case SHAFRANOV:
                    gyro::functions = std::make_unique<PolarR6ZoniShafranov>();
                    break;
                case TRIANGULAR:
                    gyro::functions = std::make_unique<PolarR6ZoniTriangular>();
                    break;
                default:
                    throw std::runtime_error("Wrong choice for the geometry\n");
                    break;
                }
                break;
            case CARTESIAN_R6:
                switch (geom) {
                case CIRCULAR:
                    gyro::functions = std::make_unique<CartesianR6ZoniCircular>();
                    break;
                case SHAFRANOV:
                    gyro::functions = std::make_unique<CartesianR6ZoniShafranov>();
                    break;
                case TRIANGULAR:
                    gyro::functions = std::make_unique<CartesianR6ZoniTriangular>();
                    break;
                default:
                    throw std::runtime_error("Wrong choice for the geometry\n");
                    break;
                }
                break;
                //case FLAT:
                //    switch (geom) {
                //        case CIRCULAR: gyro::functions   = std::make_unique<FlatZoniCircular>(); break;
                //        case SHAFRANOV: gyro::functions  = std::make_unique<FlatZoniShafranov>(); break;
                //        case TRIANGULAR: gyro::functions = std::make_unique<FlatZoniTriangular>(); break;
                //    }
                //    break;
            default:
                throw std::runtime_error("Wrong choice for the problem\n");
                break;
            }
            break;
        case ZONI_SHIFTED:
            switch (prob) {
            case CARTESIAN_R2:
                switch (geom) {
                case CIRCULAR:
                    gyro::functions = std::make_unique<CartesianR2ZoniShiftedCircular>();
                    break;
                case SHAFRANOV:
                    gyro::functions = std::make_unique<CartesianR2ZoniShiftedShafranov>();
                    break;
                case TRIANGULAR:
                    gyro::functions = std::make_unique<CartesianR2ZoniShiftedTriangular>();
                    break;
                default:
                    throw std::runtime_error("Wrong choice for the geometry\n");
                    break;
                }
                break;
            case POLAR_R6:
                switch (geom) {
                case CIRCULAR:
                    gyro::functions = std::make_unique<PolarR6ZoniShiftedCircular>();
                    break;
                case SHAFRANOV:
                    gyro::functions = std::make_unique<PolarR6ZoniShiftedShafranov>();
                    break;
                case TRIANGULAR:
                    gyro::functions = std::make_unique<PolarR6ZoniShiftedTriangular>();
                    break;
                default:
                    throw std::runtime_error("Wrong choice for the geometry\n");
                    break;
                }
                break;
            case CARTESIAN_R6:
                switch (geom) {
                case CIRCULAR:
                    gyro::functions = std::make_unique<CartesianR6ZoniShiftedCircular>();
                    break;
                case SHAFRANOV:
                    gyro::functions = std::make_unique<CartesianR6ZoniShiftedShafranov>();
                    break;
                case TRIANGULAR:
                    gyro::functions = std::make_unique<CartesianR6ZoniShiftedTriangular>();
                    break;
                default:
                    throw std::runtime_error("Wrong choice for the geometry\n");
                    break;
                }
                break;
                //case FLAT:
                //    switch (geom) {
                //        case CIRCULAR: gyro::functions   = std::make_unique<FlatZoniShiftedCircular>(); break;
                //        case SHAFRANOV: gyro::functions  = std::make_unique<FlatZoniShiftedShafranov>(); break;
                //        case TRIANGULAR: gyro::functions = std::make_unique<FlatZoniShiftedTriangular>(); break;
                //    }
                //    break;
            default:
                throw std::runtime_error("Wrong choice for the problem\n");
                break;
            }
            break;
        case POISSON:
            switch (prob) {
            case CARTESIAN_R2:
                switch (geom) {
                case CIRCULAR:
                    gyro::functions = std::make_unique<CartesianR2PoissonCircular>();
                    break;
                case SHAFRANOV:
                    gyro::functions = std::make_unique<CartesianR2PoissonShafranov>();
                    break;
                case TRIANGULAR:
                    gyro::functions = std::make_unique<CartesianR2PoissonTriangular>();
                    break;
                default:
                    throw std::runtime_error("Wrong choice for the geometry\n");
                    break;
                }
                break;
            case POLAR_R6:
                switch (geom) {
                case CIRCULAR:
                    gyro::functions = std::make_unique<PolarR6PoissonCircular>();
                    break;
                case SHAFRANOV:
                    gyro::functions = std::make_unique<PolarR6PoissonShafranov>();
                    break;
                case TRIANGULAR:
                    gyro::functions = std::make_unique<PolarR6PoissonTriangular>();
                    break;
                default:
                    throw std::runtime_error("Wrong choice for the geometry\n");
                    break;
                }
                break;
            case CARTESIAN_R6:
                switch (geom) {
                case CIRCULAR:
                    gyro::functions = std::make_unique<CartesianR6PoissonCircular>();
                    break;
                case SHAFRANOV:
                    gyro::functions = std::make_unique<CartesianR6PoissonShafranov>();
                    break;
                case TRIANGULAR:
                    gyro::functions = std::make_unique<CartesianR6PoissonTriangular>();
                    break;
                default:
                    throw std::runtime_error("Wrong choice for the geometry\n");
                    break;
                }
                break;
                //case FLAT:
                //    switch (geom) {
                //        case CIRCULAR: gyro::functions   = std::make_unique<FlatPoissonCircular>(); break;
                //        case SHAFRANOV: gyro::functions  = std::make_unique<FlatPoissonShafranov>(); break;
                //        case TRIANGULAR: gyro::functions = std::make_unique<FlatPoissonTriangular>(); break;
                //    }
                //    break;
            default:
                throw std::runtime_error("Wrong choice for the problem\n");
                break;
            }
            break;
        default:
            throw std::runtime_error("Wrong choice for the alpha coefficient\n");
            break;
        }
    }
}
