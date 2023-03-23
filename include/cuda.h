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
 * \file cuda.h
 * \brief Header for the CUDA functions
 * \author M. Kuehn, C. Kruse, P. Leleux, C. Schwarz
 * \version 0.0
 */
#ifndef CUDA_HXX_
#define CUDA_HXX_

#include <stdlib.h>
#include <cstdio>
#include <iostream>
#include <vector>
#include <cusolverSp.h>
#include <cusparse.h>

// void compute_residual_gpu(int m, double* Au, double* fVec, double* res);
__device__ void apply_A_gpu(int m, double* u, double* Au, int nr_int, double* r, double* hplus, int ntheta_int,
                            double* theta, double* thetaplus, int DirBC_Interior, int mod_pk, int prob,
                            double kappa_eps, double delta_e, double Rmax);
__global__ void apply_A_SA(int m, double* u, double* Au, int nr_int, double* r, double* hplus, int ntheta_int,
                           double* theta, double* thetaplus, int DirBC_Interior, int mod_pk, int prob, double kappa_eps,
                           double delta_e, double Rmax);
__global__ void compute_residual_gpu(int m, double* u, double* Au, double* fVec, double* res, int nr_int, double* r,
                                     double* hplus, int ntheta_int, double* theta, double* thetaplus,
                                     int DirBC_Interior, int mod_pk, int prob, double kappa_eps, double delta_e,
                                     double Rmax);

__device__ double coeff_gpu(double r, double theta, int verbose, int mod_pk, int prob, double kappa_eps, double delta_e,
                            double Rmax);
__device__ double detDFinv_gpu(double r, double theta, int verbose, int mod_pk, int prob, double kappa_eps,
                               double delta_e, double Rmax);
__device__ double arr_gpu(double r, double theta, int verbose, int mod_pk, int prob, double kappa_eps, double delta_e,
                          double Rmax);
__device__ double art_gpu(double r, double theta, int verbose, int mod_pk, int prob, double kappa_eps, double delta_e,
                          double Rmax);
__device__ double att_gpu(double r, double theta, int verbose, int mod_pk, int prob, double kappa_eps, double delta_e,
                          double Rmax);

__global__ void compute_error_gpu(double* error, int m, int nr, double* r, double* theta, int ntheta_int, double* u,
                                  double Rmax, int mod_pk, double kappa_eps, double delta_e);
__device__ void def_solution_rt_gpu(double* sol, int m, double r_j, double theta_i, double Rmax, int mod_pk,
                                    double kappa_eps, double delta_e);
__device__ void trafo_gpu(double r_j, double theta_i, double* x, double* y, int mod_pk, double kappa_eps,
                          double delta_e);
__global__ void apply_restriction_bi_gpu(
    double* Pu, double* u, int m, int mc, double* thetaplus, int ntheta_int, double* hplus,
    int* coarse_nodes_list_type_shift, int* coarse_nodes_list_type_start, int k0, int k1, int k2, int k3, int k4,
    int k5, int k6, int k7, int* coarse_nodes_list_type0, int* coarse_nodes_list_type1, int* coarse_nodes_list_type2,
    int* coarse_nodes_list_type3, int* coarse_nodes_list_type4, int* coarse_nodes_list_type5,
    int* coarse_nodes_list_type6, int* coarse_nodes_list_type7, int* coarse_nodes_list_i0, int* coarse_nodes_list_i1,
    int* coarse_nodes_list_i2, int* coarse_nodes_list_i3, int* coarse_nodes_list_i4, int* coarse_nodes_list_i5,
    int* coarse_nodes_list_i6, int* coarse_nodes_list_i7, int* coarse_nodes_list_j0, int* coarse_nodes_list_j1,
    int* coarse_nodes_list_j2, int* coarse_nodes_list_j3, int* coarse_nodes_list_j4, int* coarse_nodes_list_j5,
    int* coarse_nodes_list_j6, int* coarse_nodes_list_j7);
__global__ void apply_prolongation_bi_gpu(
    double* Pu, double* u, int m, int mc, double* thetaplus, int ntheta_int, double* hplus,
    int* coarse_nodes_list_type_shift, int* coarse_nodes_list_type_start, int k0, int k1, int k2, int k3, int k4,
    int k5, int k6, int k7, int* coarse_nodes_list_type0, int* coarse_nodes_list_type1, int* coarse_nodes_list_type2,
    int* coarse_nodes_list_type3, int* coarse_nodes_list_type4, int* coarse_nodes_list_type5,
    int* coarse_nodes_list_type6, int* coarse_nodes_list_type7, int* coarse_nodes_list_i0, int* coarse_nodes_list_i1,
    int* coarse_nodes_list_i2, int* coarse_nodes_list_i3, int* coarse_nodes_list_i4, int* coarse_nodes_list_i5,
    int* coarse_nodes_list_i6, int* coarse_nodes_list_i7, int* coarse_nodes_list_j0, int* coarse_nodes_list_j1,
    int* coarse_nodes_list_j2, int* coarse_nodes_list_j3, int* coarse_nodes_list_j4, int* coarse_nodes_list_j5,
    int* coarse_nodes_list_j6, int* coarse_nodes_list_j7);
__global__ void build_fsc_gpu(double* f_sc, int smoother, int m, int nr, int ntheta_int, int delete_circles,
                              double* fVec);
__global__ void apply_Asc_ortho_gpu(int m, double* Au, double* u, int smoother, int nr_int, double* r, double* hplus,
                                    int* coarse_nodes_list_r, int ntheta_int, double* theta, double* thetaplus,
                                    int* coarse_nodes_list_theta, int delete_circles, int DirBC_Interior, int mod_pk,
                                    double r0_DB, int prob, double kappa_eps, double delta_e, double Rmax);
__device__ int get_ptr_sc(int i, int j, int smoother, int ortho, int delete_circles, int nr_int, int ntheta_int);
__global__ void empty_vect_gpu(int m, double* vect);
__global__ void usc_to_u(int msc, double* usc, double* u, int ntheta_int, int delete_circles, int smoother);
void facto_gaussian_elimination_gpu(int* row_ind, int* col_ind, double* val, int nnz, int m, double* rhs, double* sol);

#endif // CUDA_HXX_
