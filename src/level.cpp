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
 * \file level.cpp
 * \brief Header for the class level
 * \author M. Kuehn, C. Kruse, P. Leleux, C. Schwarz
 * \version 0.0
 */
#include "level.h"

/*!
 *  \brief Default Constructor of level class
 *
 *  Default constructor of the level class, sets default values for
 *  attributes.
 *
 */
level::level(int l_)
{
    l      = l_;
    nr     = 0;
    ntheta = 0;
    reset_timers();
    delete_circles = 0;

#ifdef USE_MUMPS
    init_mumps(mumps_Ac);
    if (gyro::icntl[Param::optimized] == 0) {
        for (int i = 0; i < 4; i++) {
            init_mumps(mumps_A_Zebra[i]);
        }
    }
    else {
        init_mumps(mumps_across);
    }
#endif
} /* ----- end of constructor level:level ----- */

/*!
 *  \brief Default Destructor of level class
 *
 *  Default destructor of the level class.
 *
 */
level::~level()
{
#ifdef USE_MUMPS
    finalize_mumps(mumps_Ac);
    if (gyro::icntl[Param::optimized] == 0) {
        for (int i = 0; i < 4; i++) {
            finalize_mumps(mumps_A_Zebra[i]);
        }
    }
    else {
        finalize_mumps(mumps_across);
    }
#endif
    if (delete_circles > 0) {
        for (int smoother = 0; smoother < 4; smoother++) {
            delete[] dep_Asc_ortho[smoother];
            delete[] dep_Asc[smoother];
        }
        // delete[] dep_u;
    }
} /* ----- end of destructor level::~level ----- */

/*!
 *  \brief Sets execution times to 0
 *
 *  Sets execution times to 0.
 *
 */
void level::reset_timers()
{
    t_smoothing    = 0;
    t_f_sc         = 0;
    t_Asc_ortho    = 0;
    t_Asc          = 0;
    t_get_ptr      = 0;
    t_get_stencil  = 0;
    t_get_smoother = 0;
    t_get_row      = 0;
} /* ----- end of level::reset_timers ----- */

/*!
 *  \brief Display the array theta
 *
 *  Display the array theta
 *
 */
void level::display_r()
{
    gyro::disp(r, "Coordinates r");
} /* ----- end of level::display_r ----- */

/*!
 *  \brief Display the array theta
 *
 *  Display the array theta
 *
 */
void level::display_theta()
{
    gyro::disp(theta, "Coordinates theta");
} /* ----- end of level::display_theta ----- */

/*!
 *  \brief Defines which nodes are on the boundary
 *
 *  returns ndistance to Dirichlet boundary in binary (0: o boundary, >0: not on boundary)
 *  uses tolerance of 1e-10 to define boundary
 *
 */
void level::build_bound()
{
    double tol_bound_check = gyro::dcntl[Param::tol_bound_check];
    double r0_DB           = gyro::dcntl[Param::r0_DB];
    double R               = gyro::dcntl[Param::R];

    is_bound = std::vector<int>(m);

    for (int j = 0; j < nr; j++) {
        int is_DB = fabs((r[j] - r0_DB) * (r[j] - R)) < tol_bound_check;
        for (int i = 0; i < ntheta_int; i++) {
            is_bound[j * ntheta_int + i] = is_DB;
        }
    }

    if (gyro::icntl[Param::verbose] > 5)
        for (int j = 0; j < nr; j++)
            for (int i = 0; i < ntheta_int; i++)
                std::cout << "DISTBOUNDARY (" << r[j] << ", " << theta[j] << "): " << is_bound[j * ntheta_int + i]
                          << "\n";
} /* ----- end of gyro::build_bound ----- */

/*!
 *  \brief Defines the number of entries in A
 *
 *  returns the numbr of entries in the matrix A
 *
 */
void level::define_nz()
{
    nz          = 0;
    int nr_left = nr;
    int nb_DB   = (gyro::icntl[Param::DirBC_Interior]) ? 2 : 1;
    // - Dirichlet BC nodes
    nz += nb_DB * ntheta_int;
    nr_left -= nb_DB;
    nz += nb_DB * 4 * ntheta_int; // Nodes linked to Dirichlet BC nodes
    if (gyro::icntl[Param::mod_pk] > 0)
        nz += nb_DB * 2 * ntheta_int; // Only take diagonal values if deformed circle
    nr_left -= nb_DB;
    // - Across the origin nodes
    if (!gyro::icntl[Param::DirBC_Interior]) {
        nz += 5 * ntheta_int; // accross the origin
        if (gyro::icntl[Param::mod_pk] > 0)
            nz += 2 * ntheta_int; // Only take diagonal values if deformed circle
        nr_left--;
    }
    // - Interior nodes
    nz += 5 * ntheta_int * nr_left;
    if (gyro::icntl[Param::mod_pk] > 0)
        nz += 4 * ntheta_int * nr_left; // internal circle
} /* ----- end of gyro::define_nz ----- */

/*!
 *  \brief Defines the index of an entry in A
 *
 *  returns the index of the first entry in the row corresponding to the node (r_j, theta_i)
 *
 *  \param i: the index of the theta coordinate
 *  \param j: the index of the r coordinate
 * 
 *  \return the index
 */
int level::get_ptr(int i, int j)
{
    int ptr     = 0;
    int nr_left = (j > nr_int - 2) ? nr_int - 1 : j;
    i           = (i + ntheta_int) % ntheta_int;

    // First circle
    int index_theta = (j == 0) ? i : ntheta_int;
    // - Dirichlet
    if (gyro::icntl[Param::DirBC_Interior])
        ptr += index_theta;
    // - Across the origin
    else {
        ptr += 5 * index_theta; // accross the origin
        if (gyro::icntl[Param::mod_pk] > 0)
            ptr += 2 * index_theta; // Only take diagonal values if deformed circle
    }
    nr_left--;
    if (j == 0)
        return ptr;

    // Second circle if Dirichlet
    index_theta = (j == 1) ? i : ntheta_int;
    if (gyro::icntl[Param::DirBC_Interior]) {
        ptr += 4 * index_theta; // Nodes linked to Dirichlet BC nodes
        if (gyro::icntl[Param::mod_pk] > 0)
            ptr += 2 * index_theta; // Only take diagonal values if deformed circle
        nr_left--;
    }
    if (gyro::icntl[Param::DirBC_Interior] && j == 1)
        return ptr;

    // Interior nodes
    index_theta = (j < nr_int - 1) ? i : 0;
    ptr += 5 * (ntheta_int * nr_left + index_theta);
    if (gyro::icntl[Param::mod_pk] > 0)
        ptr += 4 * (ntheta_int * nr_left + index_theta); // internal circle
    if (j < nr_int - 1)
        return ptr;

    // Penultian circle
    index_theta = (j == nr_int - 1) ? i : ntheta_int;
    // - Dirichlet
    ptr += 4 * index_theta; // Nodes linked to Dirichlet BC nodes
    if (gyro::icntl[Param::mod_pk] > 0)
        ptr += 2 * index_theta; // Only take diagonal values if deformed circle
    if (j == nr_int - 1)
        return ptr;

    // Last circle
    index_theta = (j == nr_int) ? i : ntheta_int;
    // - Dirichlet
    ptr += index_theta;

    return ptr;
} /* ----- end of gyro::get_ptr ----- */

/*!
 *  \brief Defines the index of entry for a whole radius in A
 *
 *  returns the index of all entries in the row corresponding to the nodes in r_j
 *
 *  \param j: the index of the r coordinate
 * 
 *  \return the vector of index
 *
 */
std::vector<int> level::get_ptr(int j)
{
    int ptr = 0;
    int shift;
    int nr_left = (j > nr_int - 2) ? nr_int - 1 : j;
    std::vector<int> ptr_vect(ntheta_int + 2);

    // First circle
    // - Dirichlet
    if (gyro::icntl[Param::DirBC_Interior])
        shift = 1;
    // - Across the origin
    else {
        shift = 5; // accross the origin
        if (gyro::icntl[Param::mod_pk] > 0)
            shift = 7; // Only take diagonal values if deformed circle
    }
    if (j > 0) {
        ptr += shift * ntheta_int;
        nr_left--;

        // DB_int
        if (gyro::icntl[Param::DirBC_Interior]) {
            shift = 4; // Nodes linked to Dirichlet BC nodes
            if (gyro::icntl[Param::mod_pk] > 0)
                shift = 6; // Only take diagonal values if deformed circle
            if (j > 1) {
                ptr += shift * ntheta_int;
                nr_left--;
            }
        }
    }
    // Interior nodes
    if (j > 1 || (j > 0 && !gyro::icntl[Param::DirBC_Interior])) {
        shift = 5;
        if (gyro::icntl[Param::mod_pk] > 0)
            shift = 9;
        ptr += shift * ntheta_int * nr_left;
    }
    // Penultian circle
    if (j >= nr_int - 1) {
        // - DB_ext
        shift = 4;
        if (gyro::icntl[Param::mod_pk] > 0)
            shift = 6;
    }
    // Last circle
    if (j > nr_int - 1) {
        ptr += shift * ntheta_int;
        // - Dirichlet
        shift = 1;
    }

    for (int i = 0; i < ntheta_int; i++)
        ptr_vect[i + 1] = ptr + i * shift;
    ptr_vect[0]              = ptr_vect[ntheta_int];
    ptr_vect[ntheta_int + 1] = ptr_vect[1];
    return ptr_vect;
} /* ----- end of gyro::get_ptr ----- */

/*!
 *  \brief Defines shifting of each element in the stencil for A
 *
 *  For each node, there are a certain number of entries in the matrix.
 * Here we define, for the 9 point stencil, what is the position of each entry.
 * -1 means that this entry is not part of the stencil.
 * We consider all nodes on a same row have an identical stencil.
 *
 *  \param j: the index of the r coordinate
 * 
 *  \return the corresponding stencil
 *
 */
std::vector<int> level::get_stencil(int j)
{
    std::vector<int> stencil;
    std::vector<int> stencil_DB{-1, -1, -1, -1, 0, -1, -1, -1, -1};
    std::vector<int> stencil_DB_int_next, stencil_DB_ext_next, stencil_accross, stencil_interior;
    if (gyro::icntl[Param::mod_pk] == 0) {
        stencil_DB_int_next = std::vector<int>{-1, -1, -1, 0, 1, 2, -1, 3, -1};
        stencil_DB_ext_next = std::vector<int>{-1, 0, -1, 1, 2, 3, -1, -1, -1};
        stencil_accross     = std::vector<int>{-1, 0, -1, 1, 2, 3, -1, 4, -1};
        stencil_interior    = std::vector<int>{-1, 0, -1, 1, 2, 3, -1, 4, -1};
    }
    else {
        stencil_DB_int_next = std::vector<int>{-1, -1, -1, 0, 1, 2, 3, 4, 5};
        stencil_DB_ext_next = std::vector<int>{0, 1, 2, 3, 4, 5, -1, -1, -1};
        stencil_accross     = std::vector<int>{-1, 0, -1, 1, 2, 3, 4, 5, 6};
        stencil_interior    = std::vector<int>{0, 1, 2, 3, 4, 5, 6, 7, 8};
    }

    if ((j == 0 && gyro::icntl[Param::DirBC_Interior]) || j == nr_int)
        stencil = stencil_DB;
    else if (j == 0 && !gyro::icntl[Param::DirBC_Interior])
        stencil = stencil_accross;
    else if (j == 1 && gyro::icntl[Param::DirBC_Interior])
        stencil = stencil_DB_int_next;
    else if (j == nr_int - 1)
        stencil = stencil_DB_ext_next;
    else
        stencil = stencil_interior;

    return stencil;
} /* ----- end of gyro::get_stencil ----- */

/*!
 *  \brief Defines the number of entries in P
 *
 * Each coarse node prolongation on a 9-pt stencil
 * except the first and last lines which only propagates on the same radius (3-pt)
 *  returns the numbr of entries in the matrix P
 *
 */
void level::define_nz_P()
{
    // nz_P     = 9 * ntheta_int * (nr_int / 2 - 1) / 2 + 3 * ntheta_int;
    nz_P = (nr_int / 2 + 1) * (ntheta_int / 2) // coarse nodes
           + 2 * (nr_int / 2 + 1) * ntheta_int / 2 // same theta
           + 2 * (nr_int / 2) * ntheta_int / 2 // same r
           + 4 * (nr_int / 2) * ntheta_int / 2; // diagonal relations

    nz_P_inj = mc;
    // nz_P_ex  = 7 * ntheta_int * (nr_int / 2 - 1) / 2 + 3 * ntheta_int;
    nz_P_ex = (nr_int / 2 + 1) * (ntheta_int / 2) // coarse nodes
              + 2 * (nr_int / 2 + 1) * ntheta_int / 2 // same theta
              + 2 * (nr_int / 2) * ntheta_int / 2 // same r
              + 2 * (nr_int / 2) * ntheta_int / 2; // diagonal relations

    // 9 * ntheta_int * nr_int / 4
    // - 9 * ntheta_int / 2
    // + 3 * ntheta_int;
} /* ----- end of gyro::define_nz_P ----- */

/*!
 *  \brief Defines the size and number of entries in Asc/Asc_ortho
 *
 *  returns the size and number of entries in the matrix Asc/Asc_ortho
 *
 */
void level::define_m_nz_Asc()
{
    int extrapol = gyro::icntl[Param::extrapolation] == 1 && l == 0;
    int nz_DB, nz_accross, nz_circle, nz_radial_circle, nz_radial, nz_DB_ext;
    int nz_DB_noextr, nz_accross_noextr, nz_circle_noextr, nz_radial_circle_noextr, nz_radial_noextr, nz_DB_ext_noextr;
    int nz_1st, nz_2nd;
    double half = (extrapol) ? 0.5 : 1.0;

    // Define the size of Asc
    // - not optimized: complete size
    // - optimized: just the size for 1 row
    m_sc = std::vector<int>(4);
    m_sc = std::vector<int>(4);
    if (gyro::icntl[Param::optimized] == 0) {
        m_sc[0] = ceil(delete_circles * 0.5) * ntheta_int * half;
        m_sc[1] = floor(delete_circles * 0.5) * ntheta_int;
        m_sc[2] = ntheta_int * floor((nr_int - delete_circles + 1) * half) * 0.5;
        m_sc[3] = ntheta_int * (nr_int - delete_circles + 1) * 0.5;
    }
    else {
        m_sc[0] = ntheta_int * half;
        m_sc[1] = ntheta_int;
        m_sc[2] = (nr - delete_circles + extrapol * ((nr % 2 == 0) - (delete_circles % 2 == 0))) * half;
        m_sc[3] = (nr_int - delete_circles + 1);
    }

    // Define the number of entries in Asc_ortho
    if (gyro::icntl[Param::mod_pk] == 0) {
        nz_DB_noextr            = 0;
        nz_accross_noextr       = 1;
        nz_circle_noextr        = 2;
        nz_radial_circle_noextr = 3;
        nz_radial_noextr        = 2;
        nz_DB_ext_noextr        = 2;
        nz_accross              = (extrapol) ? 3 : nz_accross_noextr;
        nz_circle               = (extrapol) ? 4 : nz_circle_noextr;
        nz_radial_circle        = (extrapol) ? 4 : nz_radial_circle_noextr;
        nz_radial               = (extrapol) ? 4 : nz_radial_noextr;
        nz_DB_ext               = (extrapol) ? 3 : nz_DB_ext_noextr;
    }
    else {
        nz_DB_noextr            = 0;
        nz_accross_noextr       = 3;
        nz_circle_noextr        = 6;
        nz_radial_circle_noextr = 7;
        nz_radial_noextr        = 6;
        nz_DB_ext_noextr        = 4;
        nz_accross              = (extrapol) ? 5 : nz_accross_noextr;
        nz_circle               = (extrapol) ? 8 : nz_circle_noextr;
        nz_radial_circle        = (extrapol) ? 8 : nz_radial_circle_noextr;
        nz_radial               = (extrapol) ? 8 : nz_radial_noextr;
        nz_DB_ext               = (extrapol) ? 5 : nz_DB_ext_noextr;
    }
    nz_radial_circle = (extrapol && delete_circles % 2 == 0) ? 0 : nz_radial_circle;
    nz_DB            = (extrapol) ? 0 : nz_DB_noextr;
    nz_1st           = (gyro::icntl[Param::DirBC_Interior]) ? nz_DB_noextr : nz_accross;
    // DB_int and Across give the same number of entries here
    nz_2nd         = (gyro::icntl[Param::DirBC_Interior]) ? nz_accross_noextr : nz_circle_noextr;
    nz_sc_ortho    = std::vector<int>(4);
    nz_sc_ortho[0] = (nz_1st + nz_circle * (ceil(delete_circles * 0.5) - 1)) * ntheta_int * half;
    nz_sc_ortho[1] = (nz_2nd + nz_circle_noextr * (floor(delete_circles * 0.5) - 1)) * ntheta_int;
    nz_sc_ortho[2] = (nz_radial_circle + floor(nz_radial * (nr_int - delete_circles - 2) * half) + nz_DB_ext + nz_DB) *
                     ntheta_int * 0.5;
    nz_sc_ortho[3] =
        (nz_radial_circle_noextr + nz_radial_noextr * (nr_int - delete_circles - 2) + nz_DB_ext_noextr + nz_DB_noextr) *
        ntheta_int * 0.5;
} /* ----- end of gyro::define_m_nz_Asc ----- */

/*!
 *  \brief Defines the size and number of entries in Asc/Asc_ortho
 *
 *  returns the size and number of entries in the matrix Asc/Asc_ortho
 *
 */
int level::define_nz_Asc_ij(int smoother, int ij, int ortho)
{
    int extrapol = gyro::icntl[Param::extrapolation] == 1 && l == 0;
    int nz_DB, nz_accross, nz_circle, nz_radial_circle, nz_radial, nz_DB_ext;
    int nz_DB_noextr, nz_accross_noextr, nz_circle_noextr, nz_radial_circle_noextr, nz_radial_noextr, nz_DB_ext_noextr;
    int nz_1st, nz_2nd;
    int nz      = 0;
    double half = (extrapol) ? 0.5 : 1.0;

    if (!ortho) {
        nz_DB_noextr            = 1;
        nz_accross_noextr       = 4;
        nz_circle_noextr        = 3;
        nz_radial_circle_noextr = 2;
        nz_radial_noextr        = 3;
        nz_DB_ext_noextr        = 2;
        nz_DB                   = (extrapol) ? 0 : nz_DB_noextr;
        nz_accross              = (extrapol) ? 2 : nz_accross_noextr;
        nz_circle               = (extrapol) ? 1 : nz_circle_noextr;
        nz_radial_circle        = (extrapol) ? 1 : nz_radial_circle_noextr;
        nz_radial_circle        = (extrapol && delete_circles % 2 == 0) ? 0 : nz_radial_circle;
        nz_radial               = (extrapol) ? 1 : nz_radial_noextr;
        nz_DB_ext               = (extrapol) ? 1 : nz_DB_ext_noextr;
        nz_1st                  = (gyro::icntl[Param::DirBC_Interior]) ? nz_DB_noextr : nz_accross;
        nz_2nd                  = nz_circle_noextr;
    }
    else {
        if (gyro::icntl[Param::mod_pk] == 0) {
            nz_DB_noextr            = 0;
            nz_accross_noextr       = 1;
            nz_circle_noextr        = 2;
            nz_radial_circle_noextr = 3;
            nz_radial_noextr        = 2;
            nz_DB_ext_noextr        = 2;
            nz_accross              = (extrapol) ? 3 : nz_accross_noextr;
            nz_circle               = (extrapol) ? 4 : nz_circle_noextr;
            nz_radial_circle        = (extrapol) ? 4 : nz_radial_circle_noextr;
            nz_radial               = (extrapol) ? 4 : nz_radial_noextr;
            nz_DB_ext               = (extrapol) ? 3 : nz_DB_ext_noextr;
        }
        else {
            nz_DB_noextr            = 0;
            nz_accross_noextr       = 3;
            nz_circle_noextr        = 6;
            nz_radial_circle_noextr = 7;
            nz_radial_noextr        = 6;
            nz_DB_ext_noextr        = 4;
            nz_accross              = (extrapol) ? 5 : nz_accross_noextr;
            nz_circle               = (extrapol) ? 8 : nz_circle_noextr;
            nz_radial_circle        = (extrapol) ? 8 : nz_radial_circle_noextr;
            nz_radial               = (extrapol) ? 8 : nz_radial_noextr;
            nz_DB_ext               = (extrapol) ? 5 : nz_DB_ext_noextr;
        }
        nz_radial_circle = (extrapol && delete_circles % 2 == 0) ? 0 : nz_radial_circle;
        nz_DB            = (extrapol) ? 0 : nz_DB_noextr;
        nz_1st           = (gyro::icntl[Param::DirBC_Interior]) ? nz_DB_noextr : nz_accross;
        // DB_int and Across give the same number of entries here
        nz_2nd = (gyro::icntl[Param::DirBC_Interior]) ? nz_accross_noextr : nz_circle_noextr;
    }
    if (smoother == 0) {
        if (ij == 0)
            nz = nz_1st * ntheta_int * half;
        else
            nz = nz_circle * ntheta_int * half;
    }
    else if (smoother == 1) {
        if (ij == 1)
            nz = nz_2nd * ntheta_int;
        else
            nz = nz_circle_noextr * ntheta_int;
    }
    else if (smoother == 2) {
        nz = nz_radial_circle + floor(nz_radial * (nr_int - delete_circles - 2) * half) + nz_DB_ext + nz_DB;
    }
    else if (smoother == 3) {
        nz = nz_radial_circle_noextr + nz_radial_noextr * (nr_int - delete_circles - 2) + nz_DB_ext_noextr +
             nz_DB_noextr;
    }
    return nz;
} /* ----- end of gyro::define_m_nz_Asc_ij ----- */

/*!
 *  \brief Defines the index of entry for a whole radius in Asc/Asc_ortho
 *
 *  returns the index of all entries in the row corresponding to the nodes in r_j
 *
 *  \param j: the index of the r coordinate
 *  \param smoother: the current smoother
 *  \param ortho: 0: Asc, 1: Asc_ortho
 * 
 *  \return the vector of index
 *
 */
std::vector<int> level::get_ptr_sc(int j, int smoother, int ortho)
{
    double t;
    TIC;

    if ((j >= delete_circles && smoother < 2) || j < delete_circles && smoother > 1)
        throw std::runtime_error("(get_ptr_sc) Incompatible radius and smoother.");

    int extrapol = gyro::icntl[Param::extrapolation] == 1 && l == 0;
    // For the radial smoother, we need to separate smoother 2 and 3 (ptr/ptr3 and shift/shift3)
    // since we explore radius per raidus
    int ptr = 0, ptr3 = 0;
    int shift, shift3, nr_left;
    int nz_DB, nz_accross, nz_circle, nz_radial_circle, nz_radial, nz_DB_ext;
    int nz_DB_noextr, nz_accross_noextr, nz_circle_noextr, nz_radial_circle_noextr, nz_radial_noextr, nz_DB_ext_noextr;
    int nz_1st, nz_2nd;
    double half = (extrapol && smoother != 1) ? 0.5 : 1.0;
    std::vector<int> ptr_vect(ntheta_int);

    if (!ortho) {
        nz_DB_noextr            = 1;
        nz_accross_noextr       = 4;
        nz_circle_noextr        = 3;
        nz_radial_circle_noextr = 2;
        nz_radial_noextr        = 3;
        nz_DB_ext_noextr        = 2;
        nz_accross              = (extrapol) ? 2 : nz_accross_noextr;
        nz_circle               = (extrapol) ? 1 : nz_circle_noextr;
        nz_radial_circle        = (extrapol) ? 1 : nz_radial_circle_noextr;
        nz_radial               = (extrapol) ? 1 : nz_radial_noextr;
        nz_DB_ext               = (extrapol) ? 1 : nz_DB_ext_noextr;
        nz_2nd                  = nz_circle_noextr;
    }
    else {
        if (gyro::icntl[Param::mod_pk] == 0) {
            nz_DB_noextr            = 0;
            nz_accross_noextr       = 1;
            nz_circle_noextr        = 2;
            nz_radial_circle_noextr = 3;
            nz_radial_noextr        = 2;
            nz_DB_ext_noextr        = 2;
            nz_accross              = (extrapol) ? 3 : nz_accross_noextr;
            nz_circle               = (extrapol) ? 4 : nz_circle_noextr;
            nz_radial_circle        = (extrapol) ? 4 : nz_radial_circle_noextr;
            nz_radial               = (extrapol) ? 4 : nz_radial_noextr;
            nz_DB_ext               = (extrapol) ? 3 : nz_DB_ext_noextr;
        }
        else {
            nz_DB_noextr            = 0;
            nz_accross_noextr       = 3;
            nz_circle_noextr        = 6;
            nz_radial_circle_noextr = 7;
            nz_radial_noextr        = 6;
            nz_DB_ext_noextr        = 4;
            nz_accross              = (extrapol) ? 5 : nz_accross_noextr;
            nz_circle               = (extrapol) ? 8 : nz_circle_noextr;
            nz_radial_circle        = (extrapol) ? 8 : nz_radial_circle_noextr;
            nz_radial               = (extrapol) ? 8 : nz_radial_noextr;
            nz_DB_ext               = (extrapol) ? 5 : nz_DB_ext_noextr;
        }
        // DB_int and Across give the same number of entries here
        nz_2nd = (gyro::icntl[Param::DirBC_Interior]) ? nz_accross_noextr : nz_circle_noextr;
    }
    nz_1st           = (gyro::icntl[Param::DirBC_Interior]) ? nz_DB_noextr : nz_accross;
    nz_radial_circle = (extrapol && delete_circles % 2 == 0) ? 0 : nz_radial_circle;
    nz_DB            = (extrapol) ? 0 : nz_DB_noextr;

    if (!ortho) {
        if (j < delete_circles) {
            if (smoother == 0) {
                shift = nz_1st;
                if (j > 0) {
                    shift = nz_circle;
                }
            }
            else if (smoother == 1) {
                shift = nz_2nd;
                if (j > 1) {
                    shift = nz_circle_noextr;
                }
            }
            // In case of extrapol and smoother==0, there are only half of points in the grid: take twice
            // the same values for each in order to respect the theta index i
            for (int i = 0; i < ntheta_int; i++)
                ptr_vect[i] = ptr + floor(i * half) * shift;
        }
        else {
            shift  = nz_radial_circle;
            shift3 = nz_radial_circle_noextr;
            if (j > delete_circles) {
                ptr += shift;
                ptr3 += shift3;

                // Radial nodes
                shift  = nz_radial;
                shift3 = nz_radial_noextr;
                ptr += shift * floor(std::min(j - delete_circles - 1, nr_int - delete_circles - 2) * half);
                ptr3 += shift3 * std::min(j - delete_circles - 1, nr_int - delete_circles - 2);
            }
            // - DB_ext
            if (j >= nr_int - 1) {
                shift  = nz_DB_ext;
                shift3 = nz_DB_ext_noextr;
            }
            // Last circle
            if (j > nr_int - 1) {
                ptr += shift;
                ptr3 += shift3;
                shift  = nz_DB;
                shift3 = nz_DB_noextr;
            }

            for (int i = 0; i < ntheta_int; i += 2) {
                ptr_vect[i]     = ptr;
                ptr_vect[i + 1] = ptr3;
            }
        }
    }
    else {
        if (j < delete_circles) {
            if (smoother == 0) {
                shift = nz_1st;
                if (j > 0) {
                    ptr += shift * ntheta_int * half;

                    shift   = nz_circle;
                    nr_left = ceil((double)(j - 1) * 0.5) - 1;
                    ptr += shift * nr_left * ntheta_int * half;
                }
            }
            else if (smoother == 1) {
                shift = nz_2nd;
                if (j > 1) {
                    ptr += shift * ntheta_int * half;

                    shift   = nz_circle_noextr;
                    nr_left = floor((double)(j - 1) * 0.5) - 1;
                    ptr += shift * nr_left * ntheta_int;
                }
            }
            // In case of extrapol and smoother==0, there are only half of points in the grid: take twice
            // the same values for each in order to respect the theta index i
            for (int i = 0; i < ntheta_int; i++)
                ptr_vect[i] = ptr + floor(i * half) * shift;
        }
        else {
            // Radial-Circle nodes
            shift  = nz_radial_circle;
            shift3 = nz_radial_circle_noextr;
            if (j > delete_circles) {
                ptr += shift * ntheta_int * 0.5;
                ptr3 += shift3 * ntheta_int * 0.5;

                // Radial nodes
                shift  = nz_radial;
                shift3 = nz_radial_noextr;
                ptr += shift * floor(std::min(j - delete_circles - 1, nr_int - delete_circles - 2) * half) *
                       ntheta_int / 2;
                ptr3 += shift3 * std::min(j - delete_circles - 1, nr_int - delete_circles - 2) * ntheta_int / 2;
            }
            // - DB_ext
            if (j >= nr_int - 1) {
                shift  = nz_DB_ext;
                shift3 = nz_DB_ext_noextr;
            }
            // Last circle
            if (j > nr_int - 1) {
                ptr += shift * ntheta_int / 2;
                ptr3 += shift3 * ntheta_int / 2;
                shift  = nz_DB;
                shift3 = nz_DB_noextr;
            }

            for (int i = 0; i < ntheta_int; i += 2) {
                ptr_vect[i]     = ptr + i / 2 * shift;
                ptr_vect[i + 1] = ptr3 + i / 2 * shift3;
            }
        }
    }

#pragma omp atomic
    t_get_ptr += TOC;

    return ptr_vect;
} /* ----- end of gyro::get_ptr_sc ----- */

/*!
 *  \brief Defines the smoother on node (i, j)
 *
 *  Defines the smoother on node (i, j)
 *
 *  \param i: the index of the theta coordinate
 *  \param j: the index of the r coordinate
 * 
 *  \return the corresponding smoother
 *
 */
int level::get_smoother(int i, int j)
{
    double t;
    TIC;

    int smoother = 0;
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
        }
        else {
            smoother = 3; //radial white (odd)
        }
    }

#pragma omp atomic
    t_get_smoother += TOC;

    return smoother;
} /* ----- end of gyro::get_smoother ----- */

/*!
 *  \brief Defines the smoother for a whole theta line
 *
 *  Defines the smoother on theta i
 *
 *  \param i: the index of the theta coordinate
 * 
 *  \return the corresponding smoothers
 *
 */
std::vector<int> level::get_smoother_circle(int i)
{
    double t;
    TIC;

    std::vector<int> smoother(delete_circles, 0);
    for (int j = 1; j < delete_circles; j += 2)
        smoother[j] = 1;

#pragma omp atomic
    t_get_smoother += TOC;
    return smoother;
} /* ----- end of gyro::get_smoother_circle ----- */

/*!
 *  \brief Defines the smoother for a whole radius
 *
 *  Defines the smoother on r j
 *
 *  \param j: the index of the r coordinate
 * 
 *  \return the corresponding smoothers
 *
 */
std::vector<int> level::get_smoother_radial(int j)
{
    double t;
    TIC;

    std::vector<int> smoother(ntheta_int, 2);
    for (int i = 1; i < ntheta_int; i += 2)
        smoother[i] = 3;

#pragma omp atomic
    t_get_smoother += TOC;
    return smoother;
} /* ----- end of gyro::get_smoother_radial ----- */

/*!
 *  \brief Defines shifting of each element in the stencil for Asc/Asc_ortho
 *
 *  For each node, there are a certain number of entries in the matrix Asc/Asc_ortho.
 * Here we define, for the 9 point stencil, what is the position of each entry.
 * -1 means that this entry is not part of the stencil.
 *
 *  \param j: the index of the r coordinate
 *  \param smoother: the current smoother
 *  \param ortho: 0: Asc, 1: Asc_ortho
 * 
 *  \return the corresponding stencil
 *
 */
std::vector<int> level::get_stencil_sc(int j, int smoother, int ortho)
{
    double t;
    TIC;

    int extrapol = gyro::icntl[Param::extrapolation] == 1 && l == 0 && (smoother == 0 || smoother == 2);
    std::vector<int> stencil, stencil_DB, stencil_DB_int, stencil_accross, stencil_circle, stencil_radial_circle,
        stencil_radial, stencil_DB_ext;

    if (!ortho) {
        if (!extrapol) {
            stencil_DB            = std::vector<int>{-1, -1, -1, -1, 0, -1, -1, -1, -1};
            stencil_circle        = std::vector<int>{-1, -1, -1, 0, 1, 2, -1, -1, -1};
            stencil_accross       = std::vector<int>{-1, 0, -1, 1, 2, 3, -1, -1, -1};
            stencil_radial_circle = std::vector<int>{-1, -1, -1, -1, 0, -1, -1, 1, -1};
            stencil_radial        = std::vector<int>{-1, 0, -1, -1, 1, -1, -1, 2, -1};
            stencil_DB_ext        = std::vector<int>{-1, 0, -1, -1, 1, -1, -1, -1, -1};
        }
        else {
            stencil_DB            = std::vector<int>{-1, -1, -1, -1, 0, -1, -1, -1, -1};
            stencil_accross       = std::vector<int>{-1, 0, -1, -1, 1, -1, -1, -1, -1};
            stencil_circle        = stencil_DB;
            stencil_radial_circle = stencil_DB;
            stencil_radial        = stencil_DB;
            stencil_DB_ext        = stencil_DB;
        }
        stencil_DB_int = stencil_circle;
    }
    else {
        if (gyro::icntl[Param::mod_pk] == 0) {
            if (!extrapol) {
                stencil_DB            = std::vector<int>{-1, -1, -1, -1, -1, -1, -1, -1, -1};
                stencil_accross       = std::vector<int>{-1, -1, -1, -1, -1, -1, -1, 0, -1};
                stencil_circle        = std::vector<int>{-1, 0, -1, -1, -1, -1, -1, 1, -1};
                stencil_radial_circle = std::vector<int>{-1, 0, -1, 1, -1, 2, -1, -1, -1};
                stencil_radial        = std::vector<int>{-1, -1, -1, 0, -1, 1, -1, -1, -1};
                stencil_DB_ext        = std::vector<int>{-1, -1, -1, 0, -1, 1, -1, -1, -1};
            }
            else {
                stencil_DB            = std::vector<int>{-1, -1, -1, -1, -1, -1, -1, -1, -1};
                stencil_accross       = std::vector<int>{-1, -1, -1, 0, -1, 1, -1, 2, -1};
                stencil_circle        = std::vector<int>{-1, 0, -1, 1, -1, 2, -1, 3, -1};
                stencil_radial_circle = stencil_circle;
                stencil_radial        = stencil_circle;
                stencil_DB_ext        = std::vector<int>{-1, 0, -1, 1, -1, 2, -1, -1, -1};
            }
        }
        else {
            if (!extrapol) {
                stencil_DB            = std::vector<int>{-1, -1, -1, -1, -1, -1, -1, -1, -1};
                stencil_accross       = std::vector<int>{-1, -1, -1, -1, -1, -1, 0, 1, 2};
                stencil_circle        = std::vector<int>{0, 1, 2, -1, -1, -1, 3, 4, 5};
                stencil_radial_circle = std::vector<int>{0, 1, 2, 3, -1, 4, 5, -1, 6};
                stencil_radial        = std::vector<int>{0, -1, 1, 2, -1, 3, 4, -1, 5};
                stencil_DB_ext        = std::vector<int>{0, -1, 1, 2, -1, 3, -1, -1, -1};
            }
            else {
                stencil_DB            = std::vector<int>{-1, -1, -1, -1, -1, -1, -1, -1, -1};
                stencil_accross       = std::vector<int>{-1, -1, -1, 0, -1, 1, 2, 3, 4};
                stencil_circle        = std::vector<int>{0, 1, 2, 3, -1, 4, 5, 6, 7};
                stencil_radial_circle = stencil_circle;
                stencil_radial        = stencil_circle;
                stencil_DB_ext        = std::vector<int>{0, 1, 2, 3, -1, 4, -1, -1, -1};
            }
        }
        stencil_DB_int = stencil_accross;
    }

    if (j == 0)
        stencil = (gyro::icntl[Param::DirBC_Interior]) ? stencil_DB : stencil_accross;
    else if (j == 1 && gyro::icntl[Param::DirBC_Interior])
        stencil = stencil_DB_int;
    else if (j < delete_circles)
        stencil = stencil_circle;
    else if (j == delete_circles)
        stencil = stencil_radial_circle;
    else if (j < nr_int - 1)
        stencil = stencil_radial;
    else if (j == nr_int - 1)
        stencil = stencil_DB_ext;
    else if (j == nr_int)
        stencil = stencil_DB;
    else
        throw std::runtime_error("(get_stencil_sc) Stencil not recognized.");

#pragma omp atomic
    t_get_stencil += TOC;
    return stencil;
} /* ----- end of level::get_stencil_sc ----- */

/*!
 *  \brief Row/Column for the point (i, j) for Asc
 *
 *  Returns the row/column for the point (i, j) with smoother and extrapol parameters for Asc
 *
 *  \param i: the index of the theta coordinate
 *  \param j: the index of the r coordinate
 *  \param smoother: the current smoother
 *  \param extrapol: level=0 and we use implicit extrapolation
 * 
 *  \return the corresponding row
 *
 */
int level::get_row(int i, int j, int smoother, int extrapol)
{
    double t;
    TIC;

    int row = 0;
    if (smoother < 2) {
        j   = smoother; // row per row
        row = ((j - smoother) * ntheta_int * 0.5 + i + extrapol * (smoother - 1)) / (1 + extrapol * (1 - smoother));
    }
    else if (smoother > 1) {
        // // Row-wise (mixed B n W)
        // row = (j - delete_circles - extrapol * (3 - smoother) * (delete_circles % 2 == 0)) * ntheta_int * 0.5 /
        //           (1 + extrapol * (3 - smoother)) +
        //       (i + 2 - smoother) / 2;
        // Col-wise (line by line)
        // row = ((j - delete_circles - extrapol * (3 - smoother) * (delete_circles % 2 == 0)) +
        //        (i + 2 - smoother) * (nr - delete_circles + 1 - extrapol * (3 - smoother) * (nr + 1) % 2) * 0.5) /
        //       (1 + extrapol * (3 - smoother));
        i   = 0; // column per column
        row = (j - delete_circles - extrapol * (3 - smoother) * (delete_circles % 2 == 0)) /
                  (1 + extrapol * (3 - smoother)) +
              (i + 2 - smoother) * 0.5 *
                  (nr - delete_circles + extrapol * (3 - smoother) * ((nr % 2 == 0) - (delete_circles % 2 == 0))) /
                  (1 + extrapol * (3 - smoother));
    }

#pragma omp atomic
    t_get_row += TOC;
    return row;
} /* ----- end of level::get_row ----- */

/*!
 *  \brief Row/Column for a whole radius for Asc
 *
 *  Returns the row/column for the whole radius j with smoother and extrapol parameters for Asc
 *
 *  \param j: the index of the r coordinate
 *  \param smoother: the current smoother
 *  \param extrapol: level=0 and we use implicit extrapolation
 * 
 *  \return the vector of row indices
 *
 */
std::vector<int> level::get_row(int j, int smoother, int extrapol, int local, int col_wise)
{
    double t;
    TIC;

    double coeff, coeff2, coeff3, shift, shift2, shift3;

    std::vector<int> row(ntheta_int);
    if (smoother < 2) {
        if (local)
            j = smoother; // row per row
        coeff = 1.0 / (double)(1 + extrapol * (1 - smoother));
        shift = ((j - smoother) * ntheta_int * 0.5 + extrapol * (smoother - 1));
        for (int i = 0; i < ntheta_int; i++)
            row[i] = coeff * (i + shift);
    }
    else if (smoother > 1) {
        // Row-wise (mixed B n W)
        if (!col_wise) {
            coeff = 0.5;
            shift2 =
                (j - delete_circles - extrapol * (delete_circles % 2 == 0)) * ntheta_int * 0.5 / (double)(1 + extrapol);
            shift3 = ((j - delete_circles) * ntheta_int - 1) * 0.5;
            for (int i = 0; i < ntheta_int; i += 2)
                row[i] = coeff * i + shift2;
            for (int i = 1; i < ntheta_int; i += 2)
                row[i] = coeff * i + shift3;
        }
        else {
            if (local) {
                // Col-wise (line by line)
                shift2 = (j - delete_circles - extrapol * (delete_circles % 2 == 0)) / (1 + extrapol);
                shift3 = (j - delete_circles); // - 0.5 * (nr - delete_circles);
                coeff2 = 0.5 * (nr - delete_circles + extrapol * ((nr % 2 == 0) - (delete_circles % 2 == 0))) /
                         (1 + extrapol);
                coeff3 = 0.5 * (nr - delete_circles);
                for (int i = 0; i < ntheta_int; i += 2)
                    row[i] = shift2;
                for (int i = 1; i < ntheta_int; i += 2)
                    row[i] = shift3;
            }
            else {
                shift2 = (j - delete_circles - extrapol * (delete_circles % 2 == 0)) / (1 + extrapol);
                shift3 = (j - delete_circles) - 0.5 * (nr - delete_circles);
                coeff2 = 0.5 * (nr - delete_circles + extrapol * ((nr % 2 == 0) - (delete_circles % 2 == 0))) /
                         (1 + extrapol);
                coeff3 = 0.5 * (nr - delete_circles);
                for (int i = 0; i < ntheta_int; i += 2)
                    row[i] = coeff2 * i + shift2;
                for (int i = 1; i < ntheta_int; i += 2)
                    row[i] = coeff3 * i + shift3;
            }
        }
    }

#pragma omp atomic
    t_get_row += TOC;
    return row;
} /* ----- end of level::get_row ----- */

/*!
 *  \brief Row/Column for a whole radius for Asc
 *
 *  Returns the row/column for the whole radius j with smoother and extrapol parameters for Asc
 *
 *  \param j: the index of the r coordinate
 *  \param smoother: the current smoother
 *  \param extrapol: level=0 and we use implicit extrapolation
 * 
 *  \return the vector of row indices
 *
 */
std::vector<int> level::get_row_i(int size, int i, int smoother, int extrapol)
{
    double t;
    TIC;

    double coeff, shift;

    std::vector<int> row(size);
    if (smoother < 2) {
        std::cout << "SHOULD NOT BE CALLED WITH THIS SMOOTHER\n";
    }
    else if (smoother > 1) {
        coeff = 1.0 / (1 + (3 - smoother) * extrapol);
        shift = -delete_circles - (3 - smoother) * extrapol * (delete_circles % 2 == 0);
        for (int j = 0; j < size; j++)
            row[j] = coeff * (j + shift);
    }

#pragma omp atomic
    t_get_row += TOC;
    return row;
} /* ----- end of level::get_row ----- */

/*!
 *  \brief Row/Column for a whole radius for Asc
 *
 *  Returns the row/column for the whole radius j with smoother and extrapol parameters for Asc
 *
 *  \param j: the index of the r coordinate
 *  \param smoother: the current smoother
 *  \param extrapol: level=0 and we use implicit extrapolation
 * 
 *  \return the vector of row indices
 *
 */
std::vector<int> level::get_row_i_glob(int size, int i, int smoother, int extrapol)
{
    double t;
    TIC;

    double coeff, shift;

    std::vector<int> row(size);
    if (smoother < 2) {
        std::cout << "SHOULD NOT BE CALLED WITH THIS SMOOTHER\n";
    }
    else if (smoother > 1) {
        coeff = 1.0 / (1 + (3 - smoother) * extrapol);
        shift = (-delete_circles - extrapol * (3 - smoother) * (delete_circles % 2 == 0)) +
                (i + 2 - smoother) * 0.5 *
                    (nr - delete_circles + extrapol * (3 - smoother) * ((nr % 2 == 0) - (delete_circles % 2 == 0)));
        for (int j = 0; j < size; j++)
            row[j] = coeff * (j + shift);
    }

#pragma omp atomic
    t_get_row += TOC;
    return row;
} /* ----- end of level::get_row ----- */

/*!
 *  \brief Row/Column for a whole radius for Asc
 *
 *  Returns the row/column for the whole radius j with smoother and extrapol parameters for Asc
 *
 *  \param j: the index of the r coordinate
 *  \param smoother: the current smoother
 *  \param extrapol: level=0 and we use implicit extrapolation
 * 
 *  \return the vector of row indices
 *
 */
int level::mapping_usc_to_u(int ind_sc, int smoother)
{
    double t;
    TIC;

    int extrapol = gyro::icntl[Param::extrapolation];
    extrapol     = extrapol == 1 && l == 0;
    //only for the circle black smoother, don't change the variable of the class itself

    //computation of indices in the total vector u corresponding to the indices in u_sc
    // for (long unsigned int ind_sc = 0; ind_sc < u_sc.size(); ++i) {
    int row;
    int col;

    if (smoother < 2) { //circle
        int ntheta_int_local = ntheta_int;
        if (extrapol && smoother == 0) { //circle, black
            ntheta_int_local = ntheta_int / 2;
        }
        row = 2 * (ind_sc / ntheta_int_local); //row within the smoother
        col = ind_sc % ntheta_int_local;
        if (smoother == 1) { //white
            row++;
        }
        if (extrapol && smoother == 0) { //black
            col = col * 2 + 1; //augment col in case of extrapolation
        }
    }
    else { //radial
        // // Row-wise (mixed B n W)
        // if (smoother == 2) { //black
        //     int n_lines_radial_b = ceil((double)ntheta_int / 2);
        //     row                  = ind_sc / n_lines_radial_b; //row within the smoother
        //     col                  = 2 * (ind_sc % n_lines_radial_b); //col within the smoother
        //     if (extrapol) {
        //         row = row * 2; //augment row in case of extrapolation
        //         if (delete_circles % 2 == 0) { //delete_circles = even
        //             row = row + 1;
        //         }
        //     }
        // }
        // else { //white
        //     int n_lines_radial_w = floor((double)ntheta_int / 2);
        //     row                  = ind_sc / n_lines_radial_w; //row within the smoother
        //     col                  = 2 * (ind_sc % n_lines_radial_w) + 1; //col within the smoother
        // }
        // Col-wise (line by line)
        int n_rows = (nr - delete_circles + extrapol * (3 - smoother) * ((nr % 2 == 0) - (delete_circles % 2 == 0))) /
                     (1 + extrapol * (3 - smoother));
        col = 2 * (ind_sc / n_rows);
        if (smoother == 3)
            col++;
        row = ind_sc % n_rows;
        if (extrapol && smoother == 2)
            row = 2 * row + (delete_circles + 1) % 2;

        row += delete_circles; //row within the whole matrix
    }

    int ind = row * ntheta_int + col;

#pragma omp atomic
    t_get_row += TOC;

    return ind;
} /* ----- end of level::mapping_usc_to_u ----- */

/*!
 *  \brief Row/Column for a whole radius for Asc
 *
 *  Returns the row/column for the whole radius j with smoother and extrapol parameters for Asc
 *
 *  \param j: the index of the r coordinate
 *  \param smoother: the current smoother
 *  \param extrapol: level=0 and we use implicit extrapolation
 * 
 *  \return the vector of row indices
 *
 */
std::vector<int> level::mapping_usc_to_u(int ind_sc_start, int ind_sc_end, int smoother)
{
    double t;
    TIC;

    int extrapol = gyro::icntl[Param::extrapolation];
    extrapol     = extrapol == 1 && l == 0;
    //only for the circle black smoother, don't change the variable of the class itself

    //computation of indices in the total vector u corresponding to the indices in u_sc
    // for (long unsigned int ind_sc = 0; ind_sc < u_sc.size(); ++i) {
    int row;
    int col;
    std::vector<int> ind(ind_sc_end - ind_sc_start);

    if (smoother < 2) { //circle
        int ntheta_int_local = ntheta_int;
        int is_extrapol      = 0;
        if (extrapol && smoother == 0) { //circle, black
            ntheta_int_local = ntheta_int / 2;
            is_extrapol      = 1;
        }
        for (int i = ind_sc_start; i < ind_sc_end; i++) {
            row                   = 2 * (i / ntheta_int_local) + (smoother == 1); //row within the smoother
            col                   = i % ntheta_int_local;
            col                   = col * (1 + is_extrapol) + is_extrapol;
            ind[i - ind_sc_start] = row * ntheta_int + col;
        }
    }
    else { //radial
        // Col-wise (line by line)
        int n_rows = (nr - delete_circles + extrapol * (3 - smoother) * ((nr % 2 == 0) - (delete_circles % 2 == 0))) /
                     (1 + extrapol * (3 - smoother));
        int is_extrapol = 0;
        if (extrapol && smoother == 2)
            is_extrapol = 1;
        for (int i = ind_sc_start; i < ind_sc_end; i++) {
            col = 2 * (i / n_rows) + (smoother == 3);
            row = i % n_rows;
            row = row * (1 + is_extrapol) + is_extrapol * ((delete_circles + 1) % 2);
            row += delete_circles;
            ind[i - ind_sc_start] = row * ntheta_int + col;
        }
    }

#pragma omp atomic
    t_get_row += TOC;

    return ind;
} /* ----- end of level::mapping_usc_to_u ----- */

/*!
 *  \brief Creates grid division in r from a file
 *
 *  Creates grid division in r from a file containing comma-separated radii
 *
 */
void level::read_grid_r()
{
    std::fstream f;

    f.open(gyro::f_grid_r.c_str(), std::ios::in);
    if (!f) {
        std::cout << "No such file: " << gyro::f_grid_r.c_str() << "\n";
    }
    else {
        while (1) {
            double r_val;
            f >> r_val;
            if (f.eof())
                break;
            r.push_back(r_val);
        }
    }
    f.close();
    nr = r.size();
} /* ----- end of level::read_grid ----- */

/*!
 *  \brief Creates grid division in theta from a file
 *
 *  Creates grid division in theta from a file containing comma-separated radii
 *
 */
void level::read_grid_theta()
{
    std::fstream f;

    f.open(gyro::f_grid_theta.c_str(), std::ios::in);
    if (!f) {
        std::cout << "No such file: " << gyro::f_grid_theta.c_str() << "\n";
    }
    else {
        while (1) {
            double theta_val;
            f >> theta_val;
            if (f.eof())
                break;
            theta.push_back(theta_val);
        }
    }
    f.close();
    ntheta = theta.size();
} /* ----- end of level::read_grid ----- */

/*!
 *  \brief Creates grid division in r from a file
 *
 *  Creates grid division in r from a file containing comma-separated radii
 *
 */
void level::write_grid_r()
{
    std::fstream f;
    std::stringstream folder;

    size_t pos = 0;
    std::string token;
    std::string original = gyro::f_grid_r;
    while ((pos = original.find("/")) != std::string::npos) {
        token = original.substr(0, pos);
        folder << token << "/";
        original.erase(0, pos + 1);
        const int dir_err0 = mkdir(folder.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    }

    f.open(gyro::f_grid_r.c_str(), std::ios::out);
    if (!f) {
        std::cout << "No such file: " << gyro::f_grid_r.c_str() << "\n";
    }
    else {
        for (int i = 0; i < r.size(); i++) {
            f << std::setprecision(17) << r[i] << "\n";
        }
    }
    f.close();
} /* ----- end of level::read_grid ----- */

/*!
 *  \brief Creates grid division in theta from a file
 *
 *  Creates grid division in theta from a file containing comma-separated radii
 *
 */
void level::write_grid_theta()
{
    std::fstream f;
    std::stringstream folder;

    size_t pos = 0;
    std::string token;
    std::string original = gyro::f_grid_theta;
    while ((pos = original.find("/")) != std::string::npos) {
        token = original.substr(0, pos);
        folder << token << "/";
        original.erase(0, pos + 1);
        const int dir_err0 = mkdir(folder.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    }

    f.open(gyro::f_grid_theta.c_str(), std::ios::out);
    if (!f) {
        std::cout << "No such file: " << gyro::f_grid_theta.c_str() << "\n";
    }
    else {
        for (int i = 0; i < theta.size(); i++) {
            f << std::setprecision(17) << theta[i] << "\n";
        }
    }
    f.close();
} /* ----- end of level::read_grid ----- */

/*!
 *  \brief Reads the theoretical solution
 *
 *  Reads the theoretical solution from an input file with 1 entry per line
 *
 */
void level::read_sol()
{
    std::fstream f;

    f.open(gyro::f_sol_in.c_str(), std::ios::in);
    if (!f) {
        std::cout << "No such file: " << gyro::f_sol_in.c_str() << "\n";
    }
    else {
        while (1) {
            double r_val;
            f >> r_val;
            if (f.eof())
                break;
            sol_in.push_back(r_val);
        }
    }
    f.close();
    if (sol_in.size() != (size_t)m)
        throw std::runtime_error("The provided theoretical solution does not have the same size m as the system.");
} /* ----- end of level::read_grid ----- */

/*!
 *  \brief Creates grid division in theta from a file
 *
 *  Creates grid division in theta from a file containing comma-separated radii
 *
 */
void level::write_sol()
{
    std::fstream f;
    std::stringstream folder;

    size_t pos = 0;
    std::string token;
    std::string original = gyro::f_sol_out;
    while ((pos = original.find("/")) != std::string::npos) {
        token = original.substr(0, pos);
        folder << token << "/";
        original.erase(0, pos + 1);
        const int dir_err0 = mkdir(folder.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    }

    f.open(gyro::f_sol_out.c_str(), std::ios::out);
    if (!f) {
        std::cout << "No such file: " << gyro::f_sol_out.c_str() << "\n";
    }
    else {
        f << "# nr : " << r.size() << "  ntheta : " << theta.size() << std::endl;
        double kappa_eps = gyro::dcntl[Param::kappa_eps];
        double delta_e   = gyro::dcntl[Param::delta_e];
        double Rmax      = gyro::dcntl[Param::R];
        for (int i = 0; i < u.size(); i++) {
            double r_i     = r[int(i / theta.size())];
            double theta_i = theta[i % theta.size()];
            double x       = gyro::functions->x(r_i, theta_i, kappa_eps, delta_e, Rmax);
            double y       = gyro::functions->y(r_i, theta_i, kappa_eps, delta_e, Rmax);
            f << std::setprecision(17) << r_i << " " << theta_i << " " << x << " " << y << " " << u[i] << "\n";
        }
    }
    f.close();
} /* ----- end of level::read_grid ----- */
