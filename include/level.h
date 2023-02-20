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

// The class for a level in the gmgpolar, containing the disk-like shape
// geometry of the problem, the operator, etc.

/*!
 * \file level.h
 * \brief Header for the class level
 * \author M. Kuehn, C. Kruse, P. Leleux
 * \version 0.0
 */
#ifndef LEVEL_HXX_
#define LEVEL_HXX_

#include <stdlib.h>
#include <cstdio>
#include <iostream>
#include <vector>
#include <iterator>
#include <set>
#include <cmath>
#include <algorithm>
#include <map>
#include <fstream>
#include <iomanip>
#include <sys/stat.h>
#include "gyro.h"

#ifdef USE_MUMPS
    #include "mpi.h"
    #include "dmumps_c.h"

    #define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */
    #define CNTL(I) cntl[(I)-1] /* macro s.t. indices match documentation */
#endif

class level
{
public:
    /*******************************************************************************
 * Attributes
 ******************************************************************************/
    /*! the level (0=finest) */
    int l;

    /*! size of r and theta */
    int nr, ntheta; // Number of nodes
    int ntheta_int, nr_int; // Number of intervals

    /* execution times */
    double t_smoothing, t_f_sc, t_Asc_ortho, t_Asc;
    double t_get_ptr, t_get_stencil, t_get_smoother, t_get_row;

    /*! Coordinates in r and theta directions */
    std::vector<double> r;
    std::vector<int> is_bound;
    std::vector<double> theta, thetaplus, hplus;
    std::vector<double> theta_per, thetaplus_per, cos_theta, sin_theta, cos_theta_per, sin_theta_per;
    std::vector<double> theta_PI, cos_theta_PI, sin_theta_PI;

    /*! Coarse nodes */
    int coarse_nodes;
    std::vector<int> coarse_nodes_list_r;
    std::vector<int> coarse_nodes_list_theta;

    /*! Operator */
    int m, nz;
    std::vector<int> row_indices, col_indices;
    std::vector<double> vals;
    // Factorization of the coarse operator
    // - using in-house solver
    std::vector<int> row_Ac_LU, col_Ac_LU;
    std::vector<double> vals_Ac_LU;
#ifdef USE_MUMPS
    // - using MUMPS
    DMUMPS_STRUC_C mumps_Ac;
    DMUMPS_STRUC_C mumps_across;
#endif

    /*! Beta coefficient update */
    std::vector<double> betaVec;

    /*! RHS */
    std::vector<double> fVec;
    std::vector<double> fVec_initial;

    /*! Solution */
    std::vector<double> sol_in;

    /*! Prolongation */
    int mc, nz_P, nz_P_inj, nz_P_ex;
    std::vector<int> ri_prol, ci_prol;
    std::vector<double> v_prol;
    std::vector<int> ri_prol_inj, ci_prol_inj;
    std::vector<double> v_prol_inj;
    std::vector<int> ri_prol_ex, ci_prol_ex;
    std::vector<double> v_prol_ex;

    /*! Smoother */
    int delete_circles;
    std::vector<int> zebra1, zebra2;
    std::vector<std::vector<int>> coloring;
    std::vector<int> zebra_BnW;
    // size and number of entries in Asc/sc_ortho
    std::vector<int> m_sc;
    std::vector<int> nz_sc;
    std::vector<int> nz_sc_ortho;
    // Asc and Asc_ortho matrices stored in 6 2D vectors:
    // - smoother(Circle-Radial) + color(Black/White) (array)
    // - ij (vect)
    // - row/col/val (vect)
    std::vector<std::vector<int>> A_Zebra_r;
    std::vector<std::vector<int>> A_Zebra_c;
    std::vector<std::vector<double>> A_Zebra_v;
    std::vector<std::vector<int>> A_Zebra_r_row[4];
    std::vector<std::vector<int>> A_Zebra_c_row[4];
    std::vector<std::vector<double>> A_Zebra_v_row[4];

    // Number of blocks per smoother (and maximum unmber of block as 5th element)
    std::vector<int> nblocks;
    std::vector<std::vector<int>> A_Zebra_Mix_r;
    std::vector<std::vector<int>> A_Zebra_Mix_c;
    std::vector<std::vector<double>> A_Zebra_Mix_v;
    std::vector<std::vector<int>> A_Zebra_Mix_r_row[4];
    std::vector<std::vector<int>> A_Zebra_Mix_c_row[4];
    std::vector<std::vector<double>> A_Zebra_Mix_v_row[4];
    // Vectors necessary for the parallel application of Asc_ortho in multigrid_smoothing
    std::vector<int> shift_vect_s2, shift_vect_s3, ptr_vect_s2, ptr_vect_s3;
    // Factorization of Asc
    // - using in-house solver
    std::vector<std::vector<int>> A_Zebra_r_LU;
    std::vector<std::vector<int>> A_Zebra_c_LU;
    std::vector<std::vector<double>> A_Zebra_v_LU;
    std::vector<std::vector<int>> A_Zebra_r_LU_row[4];
    std::vector<std::vector<int>> A_Zebra_c_LU_row[4];
    std::vector<std::vector<double>> A_Zebra_v_LU_row[4];
#ifdef USE_MUMPS
    // - using MUMPS
    DMUMPS_STRUC_C mumps_A_Zebra[4];

    std::vector<DMUMPS_STRUC_C> mumps_A_Zebra_row[4];
#endif

    /* Dependencies */
    // int *dep_Asc_ortho, *dep_Asc, *dep_u;
    std::vector<int*> dep_Asc_ortho, dep_Asc;
    std::vector<int> size_Asc_ortho, size_Asc;

    /*! Solution */
    std::vector<double> u;
    // u from the previous smoothing procedure
    // - for the circular smoother
    std::vector<double> u_previous_c;
    // - for the radial smoother
    std::vector<double> u_previous_r;
    std::vector<double> res; // residual

    /*******************************************************************************
 * Methods
 ******************************************************************************/
    level(int l_);
    ~level();
    void reset_timers();

    /***************************************************************************
     * Geometry
     **************************************************************************/
    void build_r();
    void read_grid_r();
    void write_grid_r();
    void build_theta();
    void read_grid_theta();
    void write_grid_theta();
    void display_r();
    void display_theta();

    /***************************************************************************
     * gmgpolar
     **************************************************************************/
    void build_bound();
    void define_coarse_nodes_onelevel(level* finer);
    void store_theta_n_co();
    // void define_colors();
    void define_line_splitting();

    /* System */
    // original versions (deprecated)
    void build_A0();
    void apply_A0(std::vector<double> u, std::vector<double>& Au);
    void build_rhs0();
    // optimized
    void define_nz();
    int get_ptr(int i, int j);
    std::vector<int> get_ptr(int j);
    std::vector<int> get_stencil(int j);
    void build_A();
    void apply_A(std::vector<double> u, std::vector<double>& Au);
    void build_rhs();
    void build_betaVec();
    void read_sol();
    void write_sol();

    /* Prolongator */
    // original versions (deprecated)
    std::vector<double> apply_prolongation_bi0(std::vector<double> u, int mc, int ncoarse, std::vector<int> coarse_r,
                                               std::vector<int> coarse_theta, int trans);
    std::vector<double> apply_prolongation_inj0(std::vector<double> u, int mc, int ncoarse, std::vector<int> coarse_r,
                                                std::vector<int> coarse_theta, int trans);
    std::vector<double> apply_prolongation_ex0(std::vector<double> u, int mc, int ncoarse, std::vector<int> coarse_r,
                                               std::vector<int> coarse_theta, int trans);
    // optimized
    void define_nz_P();
    void build_prolongation_bi();
    std::vector<double> apply_prolongation_bi(std::vector<double> u);
    std::vector<double> apply_restriction_bi(std::vector<double> u);
    void build_prolongation_inj();
    std::vector<double> apply_prolongation_inj(std::vector<double> u);
    std::vector<double> apply_restriction_inj(std::vector<double> u);
    void build_prolongation_ex();
    std::vector<double> apply_prolongation_ex(std::vector<double> u);
    std::vector<double> apply_restriction_ex(std::vector<double> u);

    /* Smoothing */
    void multigrid_smoothing(int smoother, int v, std::vector<double>& f_Asc_u, int nblocks, int c, int* dep_Asc_cur,
                             int* dep_Asc_prev, int* dep_Asc1, int* dep_Asc_ortho_cur);
    // void build_fsc(std::vector<double>& f_sc, int smoother);
    // void build_fsc(std::vector<double>& f_sc, std::vector<double>& f, int smoother, int loc_to_glob);
    void build_fsc(std::vector<double>& f_sc, std::vector<double>& f, int smoother, int loc_to_glob, int start,
                   int end);
    // original versions (deprecated)
    void multigrid_smoothing0(int smoother);
    void build_fsc0(std::vector<double>& f_sc, int smoother);
    void build_Asc0();
    void apply_Asc_ortho0(std::vector<double>& Au, int smoother);
    // optimized
    void define_m_nz_Asc();
    int define_nz_Asc_ij(int smoother, int ij, int ortho);
    std::vector<int> get_ptr_sc(int j, int smoother, int ortho);
    int get_smoother(int i, int j);
    std::vector<int> get_smoother_circle(int i);
    std::vector<int> get_smoother_radial(int j);
    std::vector<int> get_stencil_sc(int j, int smoother, int ortho);
    int get_row(int i, int j, int smoother, int extrapol);
    int mapping_usc_to_u(int ind_sc, int smoother);
    std::vector<int> mapping_usc_to_u(int ind_sc_start, int ind_sc_end, int smoother);
    std::vector<int> get_row(int j, int smoother, int extrapol, int local, int col_wise);
    std::vector<int> get_row_i(int i, int size, int smoother, int extrapol);
    std::vector<int> get_row_i_glob(int i, int size, int smoother, int extrapol);
    void build_Asc();
    void build_Asc_ortho(int smoother);
    void apply_Asc_ortho(std::vector<double>& Au, std::vector<double>& u, int smoother, int v, int c, int* dep_Asc_cur,
                         int* dep_Asc_prev, int* dep_Asc1, int* dep_Asc_ortho_cur);
    void apply_Asc_ortho2(std::vector<double>& Au, std::vector<double> u, int smoother);
    void apply_Asc_ortho_ij(int _ij, std::vector<double>& Au, std::vector<double> u, int smoother_todo);

    /* Direct solver */
    // original versions (deprecated)
    std::vector<double> solve_gaussian_elimination_fb_subst(std::vector<int> A_row_indices,
                                                            std::vector<int> A_col_indices, std::vector<double> A_vals,
                                                            std::vector<double> f);
    double get_element(std::vector<int> A_row_indices, std::vector<int> A_col_indices, std::vector<double> A_vals,
                       int row_index, int col_index);
    void set_element(std::vector<int>& A_row_indices, std::vector<int>& A_col_indices, std::vector<double>& A_vals,
                     int row_index, int col_index, double value);
    // optimized
    void facto_gaussian_elimination(std::vector<int>& A_row_indices, std::vector<int>& A_col_indices,
                                    std::vector<double>& A_vals, int m_solution);
    std::vector<double> solve_gaussian_elimination(std::vector<int> A_row_indices, std::vector<int> A_col_indices,
                                                   std::vector<double> A_vals, std::vector<double> f);
#ifdef USE_MUMPS
    void init_mumps(DMUMPS_STRUC_C& mumps);
    void facto_mumps(DMUMPS_STRUC_C& mumps, std::vector<int> A_row_indices, std::vector<int> A_col_indices,
                     std::vector<double> A_vals, int m_solution);
    std::vector<double> solve_mumps(DMUMPS_STRUC_C& mumps, std::vector<double> f);
    void finalize_mumps(DMUMPS_STRUC_C& mumps);
#endif
    // specialized
    void fill_in_circle(int ij, int smoother);
    std::vector<double> solve_diag(std::vector<double> A_vals, std::vector<double> f);
    void facto_circle(std::vector<int>& A_row_indices, std::vector<int>& A_col_indices, std::vector<double>& A_vals,
                      int m_solution);
    std::vector<double> solve_circle(std::vector<int> A_row_indices, std::vector<int> A_col_indices,
                                     std::vector<double> A_vals, std::vector<double> f);
    void facto_radial(std::vector<int>& A_row_indices, std::vector<int>& A_col_indices, std::vector<double>& A_vals,
                      int m_solution);
    std::vector<double> solve_radial(std::vector<int> A_row_indices, std::vector<int> A_col_indices,
                                     std::vector<double> A_vals, std::vector<double> f);

private:
    /*******************************************************************************
 * Attributes
 ******************************************************************************/

    /*******************************************************************************
 * Methods
 ******************************************************************************/
};

#endif // LEVEL_HXX
