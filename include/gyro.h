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

// The class for all util fonctions and parameters including
// - parameters
// - the dirichlet boundary conditions
// - trafo functions

/*!
 * \file gyro.h
 * \brief Header for the class gyro
 * \author M. Kuehn, C. Kruse, P. Leleux
 * \version 0.0
 */
#ifndef GYRO_HXX_
#define GYRO_HXX_

#include <omp.h>
#include <stdlib.h>
#include <cstdio>
#include <iostream>
#include <memory>
#include <vector>
#include <stdexcept>
#include <iterator>

#include "constants.h"
#include "exact_funcs.h"

// #define TIC t = clock()
// #define TOC ((double)(clock() - t)) / CLOCKS_PER_SEC
#define TIC t = omp_get_wtime()
#define TOC (omp_get_wtime() - t)

namespace gyro
{
// /*******************************************************************************
//  * Timings
//  ******************************************************************************/
// extern double t_coeff; // compute coeff alpha
// extern double t_arr_art_att; // compute arr/art/Att
// extern double t_sol; // compute the exact solution (for Dirichlet, RHS, and error)
// extern double t_detDFinv; // compute the inverse of the determinant from the mapping of curvilinear coordinates
// extern double t_trafo; // compute transformation cartesian-curvilinear

/*******************************************************************************
 * Attributes
 ******************************************************************************/
/***************************************************************************
     * Controls and Informations
     **************************************************************************/
/*! The integer control array, see Controls::icontrols for the
     *  possible values
     */
extern std::vector<int> icntl;

/*! The real control array, see Controls::dcontrols for the
     *  possible values
     */
extern std::vector<double> dcntl;

/*! The integer info output array, see Controls::info */
extern std::vector<int> info;
/*! The real info output array, see Controls::dinfo */
extern std::vector<double> dinfo;
/*! The exact functions for the configuration*/
extern std::unique_ptr<ExactFuncs> functions;

/* File names containing the grid
 * Comma separated list of radii
 */
extern std::string f_grid_r, f_grid_theta;

/* File name containing the solution as input
 * 1 value per line
 */
extern std::string f_sol_in, f_sol_out;

/*******************************************************************************
 * Methods
 ******************************************************************************/
/***************************************************************************
     * Parameters
     **************************************************************************/
void init_params();
void show_params();
void get_geometry_coeffs(geometry_type geom);
void select_functions_class(int alpha_coeff, int beta_coeff, int geom, int prob);

/***************************************************************************
     * Boundary and solution
     **************************************************************************/
/************************
 * Single
 ************************/
double distBoundary(double x, double y, int verbose);
double eval_def_rhs(double r, double theta, int verbose);
int sign(double n);
double def_solution_rt(double r, double t, int verbose);
/************************
 * Vector
 ************************/
std::vector<double> eval_def_rhs(double r, std::vector<double> theta, std::vector<double> sin_theta,
                                 std::vector<double> cos_theta, int ntheta, int verbose);
std::vector<double> def_solution_rt(double r_i, std::vector<double> theta, std::vector<double> sin_theta,
                                    std::vector<double> cos_theta, int ntheta, int verbose);

/***************************************************************************
     * Diffusivity and operator
     **************************************************************************/
/************************
 * Single
 ************************/
double coeff(double r, int verbose);
double coeff_beta(double r, int verbose);
double detDFinv(double r, double theta, int verbose);
double arr(double r, double theta, int verbose);
double art(double r, double theta, int verbose);
double att(double r, double theta, int verbose);
void arr_att_art(double r, double theta, double& arr, double& att, double& art, int verbose);
/************************
 * Vector
 ************************/
std::vector<double> coeff(std::vector<double> r, int verbose);
std::vector<double> coeff_beta(std::vector<double> r, int verbose);
std::vector<double> detDFinv(double r, std::vector<double> theta, std::vector<double> sin_theta,
                             std::vector<double> cos_theta, int ntheta, int verbose);
std::vector<double> detDFinv(std::vector<double> r, double theta, int verbose);
std::vector<double> arr(double r, std::vector<double> theta, std::vector<double> sin_theta,
                        std::vector<double> cos_theta, int ntheta, int verbose);
std::vector<double> art(double r, std::vector<double> theta, std::vector<double> sin_theta,
                        std::vector<double> cos_theta, int ntheta, int verbose);
std::vector<double> att(double r, std::vector<double> theta, std::vector<double> cos_theta,
                        std::vector<double> sin_theta, int ntheta, int verbose);
void arr_att_art(double r, std::vector<double> theta, std::vector<double> cos_theta, std::vector<double> sin_theta,
                 int ntheta, std::vector<double>& arr, std::vector<double>& att, std::vector<double>& art, int verbose);
void arr_att_art(std::vector<double> r, double theta, std::vector<double>& arr, std::vector<double>& att,
                 std::vector<double>& art, int verbose);

/***************************************************************************
     * Polar to cartesian and back
     **************************************************************************/
/************************
 * Single
 ************************/
void trafo(double& r_i, double& theta_j, double& x, double& y, int verbose);
void trafo_back(double& r_i, double& theta_j, double& x, double& y, int verbose);
/************************
 * Vector
 ************************/
void trafo(double r_i, std::vector<double> theta, std::vector<double> sin_theta, std::vector<double> cos_theta,
           int ntheta, std::vector<double>& x, std::vector<double>& y, int verbose);

/***************************************************************************
* Display arrays or vectors
**************************************************************************/
/*!
 *  \brief Forwards all values inside an array to std::cout.
 * 
 * \param na: Size of the array.
 * \param a: The array.
 * \param s_a: Name of the array (to be printed).
 *
 */
 template <typename T>
void disp(int na, T* a, const std::string& s_a)
{
    std::cout << s_a << "(" << na << "): ";
    for (int i = 0; i < na; i++)
        std::cout << a[i] << " ";
    std::cout << "\n";
} 

/*!
 *  \brief Forwards all values inside a vector to std::cout.
 * 
 * \param a: The vector.
 * \param s_a: Name of the vector (to be printed).
 *
 */
template <typename T>
void disp(std::vector<T> a, const std::string& s_a)
{
    std::cout << s_a << "(" << a.size() << "): ";
    for (std::size_t i = 0; i < a.size(); i++)
        std::cout << a[i] << " ";
    std::cout << "\n";
} 

/***************************************************************************
     * Matrix operations
     **************************************************************************/
void sp_dgemv(int trans, int m, int n, double alpha, std::vector<int> row_indices, std::vector<int> col_indices,
              std::vector<double> vals, int lda, std::vector<double> x, int incx, double beta, std::vector<double>& y,
              int incy);
double norm(std::vector<double> x);
double A_norm(std::vector<double> x, int m, std::vector<int> row_indices, std::vector<int> col_indices,
              std::vector<double> vals);
} // namespace gyro

#endif // GYRO_HXX_
