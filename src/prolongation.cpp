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
 * \file prolongation.cpp
 * \brief Implementation of the prolongation operators
 * \author M. Kuehn, C. Kruse, P. Leleux, C. Schwarz
 * \version 0.0
 */
#include "level.h"

/*!
 *  \brief Builds the bilinear interpolation
 *
 * Builds bilinear interpolation for Dirichlet boundary conditions in
 * variable r and periodic boundary conditions in phi. Can be used for a
 * disk as well as for an annulus. (might be used for other geometries with
 * Dirichlet boundary conditions in 'x' and periodic boundary conditions
 * in 'y' as well; not tested!)
 *
 * uses Anisotropic Bilinear Interpolation stencil
 * for isotropic mesh, anisotropic stencil reduces to isotropic definition
 *      ]1  2   1[
 *  1/4 ]2  4   2[
 *      ]1  2   1[
 *
 */
void level::build_prolongation_bi()
{
    int index = 0;
    int row, col;
    double h_prev, h_next, k_prev, k_next, denum, val;

    std::vector<int> iplus_vect(ntheta_int);
    for (int i = 0; i < ntheta_int - 1; i++)
        iplus_vect[i] = i + 1;
    iplus_vect[ntheta_int - 1] = 0;

    // Coarse nodes
    // => nb_nodes = ntheta_int * (nr_int / 2 + 1) / 2;
    for (int j = 0; j < nr; j += 2) {
        for (int i = 0; i < ntheta_int; i += 2) {
            row            = j * ntheta_int + i;
            col            = j * ntheta_int / 4 + i / 2;
            val            = 1.0;
            ri_prol[index] = row;
            ci_prol[index] = col;
            v_prol[index]  = val;
            index++;
        }
    }
    // Coarse nodes in same r (2i, 2j+1)
    // => nb_nodes = ntheta_int * (nr_int / 2 + 1) / 2;
    for (int j = 0; j < nr; j += 2) {
        for (int i = 1; i < ntheta_int; i += 2) {
            row    = j * ntheta_int + i;
            k_prev = thetaplus[(i + ntheta_int - 1) % ntheta_int];
            k_next = thetaplus[i];
            denum  = 1 / (k_prev + k_next);

            // Previous coarse node (left)
            col            = j * ntheta_int / 4 + (i - 1) / 2;
            val            = k_next * denum; // 1/2
            ri_prol[index] = row;
            ci_prol[index] = col;
            v_prol[index]  = val;
            index++;

            // Next coarse node (right)
            col            = j * ntheta_int / 4 + iplus_vect[i] / 2;
            val            = k_prev * denum; // 1/2
            ri_prol[index] = row;
            ci_prol[index] = col;
            v_prol[index]  = val;
            index++;
        }
    }
    // Coarse nodes in same theta (2i+1, 2j)
    // => nb_nodes = ntheta_int * (nr_int / 2 - 2) / 2; (TODO correct comment)
    for (int j = 1; j <= nr_int - 1; j += 2) {
        for (int i = 0; i < ntheta_int; i += 2) {
            row    = j * ntheta_int + i;
            h_prev = hplus[j - 1];
            h_next = hplus[j];
            denum  = 1 / (h_prev + h_next);

            // Previous coarse node (bottom)
            col            = (j - 1) * ntheta_int / 4 + i / 2;
            val            = h_next * denum; // 1/2
            ri_prol[index] = row;
            ci_prol[index] = col;
            v_prol[index]  = val;
            index++;

            // Next coarse node (top)
            col            = (j + 1) * ntheta_int / 4 + i / 2;
            val            = h_prev * denum; // 1/2
            ri_prol[index] = row;
            ci_prol[index] = col;
            v_prol[index]  = val;
            index++;
        }
    }
    // Coarse nodes in diagonals (2i+1, 2j+1)
    // => nb_nodes = ntheta_int * (nr_int / 2 - 2) / 2; (TODO correct comment)
    for (int j = 1; j <= nr_int - 1; j += 2) {
        for (int i = 1; i < ntheta_int; i += 2) {
            row    = j * ntheta_int + i;
            k_prev = thetaplus[i - 1];
            k_next = thetaplus[i];
            h_prev = hplus[j - 1];
            h_next = hplus[j];
            denum  = 1 / ((k_prev + k_next) * (h_prev + h_next));

            // top_right
            col            = (j + 1) * ntheta_int / 4 + iplus_vect[i] / 2;
            val            = h_prev * k_prev * denum; // isotrop: 1/4
            ri_prol[index] = row;
            ci_prol[index] = col;
            v_prol[index]  = val;
            index++;
            // top_left
            col            = (j + 1) * ntheta_int / 4 + (i - 1) / 2;
            val            = h_prev * k_next * denum; // isotrop: 1/4
            ri_prol[index] = row;
            ci_prol[index] = col;
            v_prol[index]  = val;
            index++;

            // bottom_right
            col            = (j - 1) * ntheta_int / 4 + iplus_vect[i] / 2;
            val            = h_next * k_prev * denum; // isotrop: 1/4
            ri_prol[index] = row;
            ci_prol[index] = col;
            v_prol[index]  = val;
            index++;
            // bottom_left
            col            = (j - 1) * ntheta_int / 4 + (i - 1) / 2;
            val            = h_next * k_next * denum; // isotrop: 1/4
            ri_prol[index] = row;
            ci_prol[index] = col;
            v_prol[index]  = val;
            index++;
        }
    }
} /* ----- end of level::build_prolongation_bi ----- */

/*!
 *  \brief Applies the bilinear interpolation
 *
 * Applies bilinear interpolation for Dirichlet boundary conditions in
 * variable r and periodic boundary conditions in phi. Can be used for a
 * disk as well as for an annulus. (might be used for other geometries with
 * Dirichlet boundary conditions in 'x' and periodic boundary conditions
 * in 'y' as well; not tested!)
 *
 * uses Anisotropic Bilinear Interpolation stencil
 * for isotropic mesh, anisotropic stencil reduces to isotropic definition
 *      ]1  2   1[
 *  1/4 ]2  4   2[
 *      ]1  2   1[
 *
 */
std::vector<double> level::apply_prolongation_bi(std::vector<double> u)
{
    //Prolongation (P * u = Pu), (m x mc) * mc = m
    std::vector<double> Pu = std::vector<double>(m, 0);

    std::vector<int> iplus_vect(ntheta_int);
    for (int i = 0; i < ntheta_int - 1; i++)
        iplus_vect[i] = i + 1;
    iplus_vect[ntheta_int - 1] = 0;

#pragma omp parallel shared(iplus_vect)
    {
#pragma omp single
        {
            // Coarse nodes
            // => nb_nodes = ntheta_int * (nr_int / 2 + 1) / 2;
            for (int j = 0; j < nr; j += 2) {
#pragma omp task firstprivate(j) shared(iplus_vect, Pu)
                {
                    int row, col;
                    double val;
                    for (int i = 0; i < ntheta_int; i += 2) {
                        row = j * ntheta_int + i;
                        col = j * ntheta_int / 4 + i / 2;
                        val = 1.0;
                        Pu[row] += val * u[col];
                    }
                } //end of task
            }
            // Coarse nodes in same r (2i, 2j+1)
            // => nb_nodes = ntheta_int * (nr_int / 2 + 1) / 2;
            for (int j = 0; j < nr; j += 2) {
#pragma omp task firstprivate(j) shared(iplus_vect, Pu)
                {
                    int row, col;
                    double k_prev, k_next, denum, val;
                    for (int i = 1; i < ntheta_int; i += 2) {
                        row    = j * ntheta_int + i;
                        k_prev = thetaplus[(i + ntheta_int - 1) % ntheta_int];
                        k_next = thetaplus[i];
                        denum  = 1 / (k_prev + k_next);

                        // Previous coarse node (left)
                        col = j * ntheta_int / 4 + (i - 1) / 2;
                        val = k_next * denum; // 1/2
                        Pu[row] += val * u[col];

                        // Next coarse node (right)
                        col = j * ntheta_int / 4 + iplus_vect[i] / 2;
                        val = k_prev * denum; // 1/2
                        Pu[row] += val * u[col];
                    }
                } //end of task
            }
            // Coarse nodes in same theta (2i+1, 2j)
            // => nb_nodes = ntheta_int * (nr_int / 2 - 2) / 2; (TODO: correct comment)
            for (int j = 1; j <= nr_int - 1; j += 2) {
#pragma omp task firstprivate(j) shared(iplus_vect, Pu)
                {
                    int row, col;
                    double h_prev, h_next, denum, val;
                    for (int i = 0; i < ntheta_int; i += 2) {
                        row    = j * ntheta_int + i;
                        h_prev = hplus[j - 1];
                        h_next = hplus[j];
                        denum  = 1 / (h_prev + h_next);

                        // Previous coarse node (bottom)
                        col = (j - 1) * ntheta_int / 4 + i / 2;
                        val = h_next * denum; // 1/2
                        Pu[row] += val * u[col];

                        // Next coarse node (top)
                        col = (j + 1) * ntheta_int / 4 + i / 2;
                        val = h_prev * denum; // 1/2
                        Pu[row] += val * u[col];
                    }
                } //end of task
            }
            // Coarse nodes in diagonals (2i+1, 2j+1)
            // => nb_nodes = ntheta_int * (nr_int / 2 - 2) / 2; (TODO: correct comment)
            for (int j = 1; j <= nr_int - 1; j += 2) {
#pragma omp task firstprivate(j) shared(iplus_vect, Pu)
                {
                    int row, col;
                    double k_prev, k_next, h_prev, h_next, denum, val;
                    for (int i = 1; i < ntheta_int; i += 2) {
                        row    = j * ntheta_int + i;
                        k_prev = thetaplus[i - 1];
                        k_next = thetaplus[i];
                        h_prev = hplus[j - 1];
                        h_next = hplus[j];
                        denum  = 1 / ((k_prev + k_next) * (h_prev + h_next));

                        // top_right
                        col = (j + 1) * ntheta_int / 4 + iplus_vect[i] / 2;
                        val = h_prev * k_prev * denum; // isotrop: 1/4
                        Pu[row] += val * u[col];
                        // top_left
                        col = (j + 1) * ntheta_int / 4 + (i - 1) / 2;
                        val = h_prev * k_next * denum; // isotrop: 1/4
                        Pu[row] += val * u[col];

                        // bottom_right
                        col = (j - 1) * ntheta_int / 4 + iplus_vect[i] / 2;
                        val = h_next * k_prev * denum; // isotrop: 1/4
                        Pu[row] += val * u[col];
                        // bottom_left
                        col = (j - 1) * ntheta_int / 4 + (i - 1) / 2;
                        val = h_next * k_next * denum; // isotrop: 1/4
                        Pu[row] += val * u[col];
                    }
                } //end of task
            }
        } //end of single
    } //end of parallel

    return Pu;
} /* ----- end of level::apply_prolongation_bi ----- */

/*!
 *  \brief Applies the bilinear interpolation
 *
 * Applies bilinear interpolation for Dirichlet boundary conditions in
 * variable r and periodic boundary conditions in phi. Can be used for a
 * disk as well as for an annulus. (might be used for other geometries with
 * Dirichlet boundary conditions in 'x' and periodic boundary conditions
 * in 'y' as well; not tested!)
 *
 * uses Anisotropic Bilinear Interpolation stencil
 * for isotropic mesh, anisotropic stencil reduces to isotropic definition
 *      ]1  2   1[
 *  1/4 ]2  4   2[
 *      ]1  2   1[
 *
 */
std::vector<double> level::apply_restriction_bi(std::vector<double> u)
{
    //Restriction (P^T * u = Pu), (mc x m) * m = mc
    std::vector<double> Pu = std::vector<double>(mc, 0);

    std::vector<int> iplus_vect(ntheta_int);
    for (int i = 0; i < ntheta_int - 1; i++)
        iplus_vect[i] = i + 1;
    iplus_vect[ntheta_int - 1] = 0;
    std::vector<int> imoins_vect(ntheta_int);
    for (int i = 1; i < ntheta_int; i++)
        imoins_vect[i] = i - 1;
    imoins_vect[0] = ntheta_int - 1;

#pragma omp parallel shared(iplus_vect)
    {
#pragma omp single
        {
            // // Loop through all coarse nodes
            // First and last radius
            // => nb_nodes = ntheta_int;
            int jc = 0;
#pragma omp task firstprivate(jc) shared(iplus_vect, imoins_vect, Pu)
            {
                int row, col, i, j, ic, i_f, j_f;
                double k_prev, k_next, h_prev, h_next, denum, val;
                for (ic = 0; ic < ntheta_int / 2; ic++) {
                    col = jc * ntheta_int / 2 + ic;
                    i_f = ic * 2;
                    j_f = jc * 2;

                    // // Coarse nodes
                    i   = i_f;
                    j   = j_f;
                    row = j * ntheta_int + i;
                    val = 1.0;
                    Pu[col] += val * u[row];

                    // // Fine nodes in same r (2i, 2j+1)
                    // Next fine node (right)
                    i      = iplus_vect[i_f];
                    j      = j_f;
                    row    = j * ntheta_int + i;
                    k_prev = thetaplus[imoins_vect[i]];
                    k_next = thetaplus[i];
                    denum  = 1 / (k_prev + k_next);
                    val    = k_next * denum; // 1/2
                    Pu[col] += val * u[row];
                    // Previous fine node (left)
                    i      = imoins_vect[i_f];
                    j      = j_f;
                    row    = j * ntheta_int + i;
                    k_prev = thetaplus[imoins_vect[i]];
                    k_next = thetaplus[i];
                    denum  = 1 / (k_prev + k_next);
                    val    = k_prev * denum; // 1/2
                    Pu[col] += val * u[row];

                    // // Fine nodes in same theta (2i+1, 2j)
                    // Next fine node (top)
                    i      = i_f;
                    j      = j_f + 1;
                    row    = j * ntheta_int + i;
                    h_prev = hplus[j - 1];
                    h_next = hplus[j];
                    denum  = 1 / (h_prev + h_next);
                    val    = h_next * denum; // 1/2
                    Pu[col] += val * u[row];

                    // // Fine nodes in diagonals (2i+1, 2j+1)
                    // top_left
                    i      = imoins_vect[i_f];
                    j      = j_f + 1;
                    row    = j * ntheta_int + i;
                    k_prev = thetaplus[imoins_vect[i]];
                    k_next = thetaplus[i];
                    h_prev = hplus[j - 1];
                    h_next = hplus[j];
                    denum  = 1 / ((k_prev + k_next) * (h_prev + h_next));
                    val    = h_next * k_prev * denum; // isotrop: 1/4
                    Pu[col] += val * u[row];
                    // top_right
                    i      = iplus_vect[i_f];
                    j      = j_f + 1;
                    row    = j * ntheta_int + i;
                    k_prev = thetaplus[imoins_vect[i]];
                    k_next = thetaplus[i];
                    h_prev = hplus[j - 1];
                    h_next = hplus[j];
                    denum  = 1 / ((k_prev + k_next) * (h_prev + h_next));
                    val    = h_next * k_next * denum; // isotrop: 1/4
                    Pu[col] += val * u[row];
                }
            } //end of task
            jc = nr_int / 2;
#pragma omp task firstprivate(jc) shared(iplus_vect, imoins_vect, Pu)
            {
                int row, col, i, j, ic, i_f, j_f;
                double k_prev, k_next, h_prev, h_next, denum, val;
                for (ic = 0; ic < ntheta_int / 2; ic++) {
                    col = jc * ntheta_int / 2 + ic;
                    i_f = ic * 2;
                    j_f = jc * 2;

                    // // Coarse nodes
                    i   = i_f;
                    j   = j_f;
                    row = j * ntheta_int + i;
                    val = 1.0;
                    Pu[col] += val * u[row];

                    // // Fine nodes in same r (2i, 2j+1)
                    // Next fine node (right)
                    i      = iplus_vect[i_f];
                    j      = j_f;
                    row    = j * ntheta_int + i;
                    k_prev = thetaplus[imoins_vect[i]];
                    k_next = thetaplus[i];
                    denum  = 1 / (k_prev + k_next);
                    val    = k_next * denum; // 1/2
                    Pu[col] += val * u[row];
                    // Previous fine node (left)
                    i      = imoins_vect[i_f];
                    j      = j_f;
                    row    = j * ntheta_int + i;
                    k_prev = thetaplus[imoins_vect[i]];
                    k_next = thetaplus[i];
                    denum  = 1 / (k_prev + k_next);
                    val    = k_prev * denum; // 1/2
                    Pu[col] += val * u[row];

                    // // Fine nodes in same theta (2i+1, 2j)
                    // Previous fine node (bottom)
                    i      = i_f;
                    j      = j_f - 1;
                    row    = j * ntheta_int + i;
                    h_prev = hplus[j - 1];
                    h_next = hplus[j];
                    denum  = 1 / (h_prev + h_next);
                    val    = h_prev * denum; // 1/2
                    Pu[col] += val * u[row];

                    // // Fine nodes in diagonals (2i+1, 2j+1)
                    // bottom_left
                    i      = imoins_vect[i_f];
                    j      = j_f - 1;
                    row    = j * ntheta_int + i;
                    k_prev = thetaplus[imoins_vect[i]];
                    k_next = thetaplus[i];
                    h_prev = hplus[j - 1];
                    h_next = hplus[j];
                    denum  = 1 / ((k_prev + k_next) * (h_prev + h_next));
                    val    = h_prev * k_prev * denum; // isotrop: 1/4
                    Pu[col] += val * u[row];
                    // bottom_right
                    i      = iplus_vect[i_f];
                    j      = j_f - 1;
                    row    = j * ntheta_int + i;
                    k_prev = thetaplus[imoins_vect[i]];
                    k_next = thetaplus[i];
                    h_prev = hplus[j - 1];
                    h_next = hplus[j];
                    denum  = 1 / ((k_prev + k_next) * (h_prev + h_next));
                    val    = h_prev * k_next * denum; // isotrop: 1/4
                    Pu[col] += val * u[row];
                }
            } //end of task
            // Interior
            // => nb_nodes = ntheta_int * (nr_int / 2 + 1) / 2; TODO: Correct comment
            for (int jc = 1; jc < nr_int / 2; jc++) {
#pragma omp task firstprivate(jc) shared(iplus_vect, imoins_vect, Pu)
                {
                    int row, col, i, j, ic, i_f, j_f;
                    double k_prev, k_next, h_prev, h_next, denum, val;
                    for (ic = 0; ic < ntheta_int / 2; ic++) {
                        col = jc * ntheta_int / 2 + ic;
                        i_f = ic * 2;
                        j_f = jc * 2;

                        // // Coarse nodes
                        i   = i_f;
                        j   = j_f;
                        row = j * ntheta_int + i;
                        val = 1.0;
                        Pu[col] += val * u[row];

                        // // Fine nodes in same r (2i, 2j+1)
                        // Next fine node (right)
                        i      = iplus_vect[i_f];
                        j      = j_f;
                        row    = j * ntheta_int + i;
                        k_prev = thetaplus[imoins_vect[i]];
                        k_next = thetaplus[i];
                        denum  = 1 / (k_prev + k_next);
                        val    = k_next * denum; // 1/2
                        Pu[col] += val * u[row];
                        // Previous fine node (left)
                        i      = imoins_vect[i_f];
                        j      = j_f;
                        row    = j * ntheta_int + i;
                        k_prev = thetaplus[imoins_vect[i]];
                        k_next = thetaplus[i];
                        denum  = 1 / (k_prev + k_next);
                        val    = k_prev * denum; // 1/2
                        Pu[col] += val * u[row];

                        // // Coarse nodes in same theta (2i+1, 2j)
                        // Next fine node (top)
                        i      = i_f;
                        j      = j_f + 1;
                        row    = j * ntheta_int + i;
                        h_prev = hplus[j - 1];
                        h_next = hplus[j];
                        denum  = 1 / (h_prev + h_next);
                        val    = h_next * denum; // 1/2
                        Pu[col] += val * u[row];
                        // Previous fine node (bottom)
                        i      = i_f;
                        j      = j_f - 1;
                        row    = j * ntheta_int + i;
                        h_prev = hplus[j - 1];
                        h_next = hplus[j];
                        denum  = 1 / (h_prev + h_next);
                        val    = h_prev * denum; // 1/2
                        Pu[col] += val * u[row];

                        // // Coarse nodes in diagonals (2i+1, 2j+1)
                        // bottom_left
                        i      = imoins_vect[i_f];
                        j      = j_f - 1;
                        row    = j * ntheta_int + i;
                        k_prev = thetaplus[imoins_vect[i]];
                        k_next = thetaplus[i];
                        h_prev = hplus[j - 1];
                        h_next = hplus[j];
                        denum  = 1 / ((k_prev + k_next) * (h_prev + h_next));
                        val    = h_prev * k_prev * denum; // isotrop: 1/4
                        Pu[col] += val * u[row];
                        // bottom_right
                        i      = iplus_vect[i_f];
                        j      = j_f - 1;
                        row    = j * ntheta_int + i;
                        k_prev = thetaplus[imoins_vect[i]];
                        k_next = thetaplus[i];
                        h_prev = hplus[j - 1];
                        h_next = hplus[j];
                        denum  = 1 / ((k_prev + k_next) * (h_prev + h_next));
                        val    = h_prev * k_next * denum; // isotrop: 1/4
                        Pu[col] += val * u[row];
                        // top_left
                        i      = imoins_vect[i_f];
                        j      = j_f + 1;
                        row    = j * ntheta_int + i;
                        k_prev = thetaplus[imoins_vect[i]];
                        k_next = thetaplus[i];
                        h_prev = hplus[j - 1];
                        h_next = hplus[j];
                        denum  = 1 / ((k_prev + k_next) * (h_prev + h_next));
                        val    = h_next * k_prev * denum; // isotrop: 1/4
                        Pu[col] += val * u[row];
                        // top_right
                        i      = iplus_vect[i_f];
                        j      = j_f + 1;
                        row    = j * ntheta_int + i;
                        k_prev = thetaplus[imoins_vect[i]];
                        k_next = thetaplus[i];
                        h_prev = hplus[j - 1];
                        h_next = hplus[j];
                        denum  = 1 / ((k_prev + k_next) * (h_prev + h_next));
                        val    = h_next * k_next * denum; // isotrop: 1/4
                        Pu[col] += val * u[row];
                    }
                } //end of task
            }
        } //end of single
    } //end of parallel

    return Pu;
} /* ----- end of level::apply_restriction_bi ----- */

/*!
 *  \brief Builds the injection
 *
 * Builds injection: directly inject values of coarse nodes.
 * Just uses values 1.0 for nodes of type 0.
 *
 */
void level::build_prolongation_inj()
{
    int index = 0;

    std::vector<int> iplus_vect(ntheta_int);
    for (int i = 0; i < ntheta_int - 1; i++)
        iplus_vect[i] = i + 1;
    iplus_vect[ntheta_int - 1] = 0;

    // Coarse nodes
    // => nb_nodes = ntheta_int * (nr_int / 2 + 1) / 2;
    for (int j = 0; j < nr; j += 2) {
        int row, col;
        double val;
        for (int i = 0; i < ntheta_int; i += 2) {
            row                = j * ntheta_int + i;
            col                = j * ntheta_int / 4 + i / 2;
            val                = 1.0;
            ri_prol_inj[index] = row;
            ci_prol_inj[index] = col;
            v_prol_inj[index]  = val;
            index++;
        }
    }
} /* ----- end of level::build_prolongation_inj ----- */

/*!
 *  \brief Applies the injection
 *
 * Applies injection: directly inject values of coarse nodes.
 * Just uses values 1.0 for nodes of type 0.
 *
 */
std::vector<double> level::apply_prolongation_inj(std::vector<double> u)
{
    std::vector<double> Pu;
    //Prolongation (P * u = Pu), (m x mc) * mc = m
    Pu = std::vector<double>(m, 0);

    std::vector<int> iplus_vect(ntheta_int);
    for (int i = 0; i < ntheta_int - 1; i++)
        iplus_vect[i] = i + 1;
    iplus_vect[ntheta_int - 1] = 0;

#pragma omp parallel shared(iplus_vect)
    {
#pragma omp single
        {
            // Coarse nodes
            // => nb_nodes = ntheta_int * (nr_int / 2 + 1) / 2;
            for (int j = 0; j < nr; j += 2) {
#pragma omp task firstprivate(j) shared(iplus_vect, Pu)
                {
                    int row, col;
                    double val;
                    for (int i = 0; i < ntheta_int; i += 2) {
                        row = j * ntheta_int + i;
                        col = j * ntheta_int / 4 + i / 2;
                        val = 1.0;
                        Pu[row] += val * u[col];
                    }
                } //end of task
            }
        } //end of single
    } //end of parallel

    return Pu;
} /* ----- end of level::apply_prolongation_inj ----- */

/*!
 *  \brief Applies the injection
 *
 * Applies injection: directly inject values of coarse nodes.
 * Just uses values 1.0 for nodes of type 0.
 *
 */
std::vector<double> level::apply_restriction_inj(std::vector<double> u)
{
    //Restriction (P^T * u = Pu), (mc x m) * m = mc
    std::vector<double> Pu = std::vector<double>(mc, 0);

#pragma omp parallel
    {
#pragma omp single
        {
            // Loop through all coarse nodes
            // => nb_nodes = ntheta_int * (nr_int / 2 + 1) / 2;
            for (int jc = 0; jc < nr_int / 2 + 1; jc++) {
#pragma omp task firstprivate(jc) shared(Pu)
                {
                    int row, col, i, j, i_f, j_f;
                    double val;
                    for (int ic = 0; ic < ntheta_int / 2; ic++) {
                        col = jc * ntheta_int / 2 + ic;
                        i_f = ic * 2;
                        j_f = jc * 2;

                        // // Coarse nodes
                        i   = i_f;
                        j   = j_f;
                        row = j * ntheta_int + i;
                        val = 1.0;
                        Pu[col] += val * u[row];
                    }
                } //end of task
            }
        } //end of single
    } //end of parallel

    return Pu;
} /* ----- end of level::apply_restriction_inj ----- */

/*!
 *  \brief Builds the extrapolation operator
 *
 * Builds the prolongation operator for implicit extrapolation.
 * 
 * The stencil is the same as the isotropic bilinear interpolation stencil,
 * except for fine points with diagonal coarse neighboors (types 5,6,7) where only
 * top_left and top_right values are used.
 *
 */
void level::build_prolongation_ex()
{
    int index = 0;

    std::vector<int> iplus_vect(ntheta_int);
    for (int i = 0; i < ntheta_int - 1; i++)
        iplus_vect[i] = i + 1;
    iplus_vect[ntheta_int - 1] = 0;

    // Coarse nodes
    // => nb_nodes = ntheta_int * (nr_int / 2 + 1) / 2;
    for (int j = 0; j < nr; j += 2) {
        int row, col;
        double val;
        for (int i = 0; i < ntheta_int; i += 2) {
            row               = j * ntheta_int + i;
            col               = j * ntheta_int / 4 + i / 2;
            val               = 1.0;
            ri_prol_ex[index] = row;
            ci_prol_ex[index] = col;
            v_prol_ex[index]  = val;
            index++;
        }
    }
    // Coarse nodes in same r (2i, 2j+1)
    // => nb_nodes = ntheta_int * (nr_int / 2 + 1) / 2;
    for (int j = 0; j < nr; j += 2) {
        int row, col;
        double val;
        for (int i = 1; i < ntheta_int; i += 2) {
            row = j * ntheta_int + i;

            // Previous coarse node (left)
            col               = j * ntheta_int / 4 + (i - 1) / 2;
            val               = 0.5; // 1/2
            ri_prol_ex[index] = row;
            ci_prol_ex[index] = col;
            v_prol_ex[index]  = val;
            index++;

            // Next coarse node (right)
            col               = j * ntheta_int / 4 + iplus_vect[i] / 2;
            val               = 0.5; // 1/2
            ri_prol_ex[index] = row;
            ci_prol_ex[index] = col;
            v_prol_ex[index]  = val;
            index++;
        }
    }
    // Coarse nodes in same theta (2i+1, 2j)
    // => nb_nodes = ntheta_int * (nr_int / 2 - 2) / 2; TODO: Correct comment
    int row, col;
    double val;
    for (int j = 1; j <= nr_int - 1; j += 2) {
        for (int i = 0; i < ntheta_int; i += 2) {
            row = j * ntheta_int + i;

            // Previous coarse node (bottom)
            col               = (j - 1) * ntheta_int / 4 + i / 2;
            val               = 0.5; // 1/2
            ri_prol_ex[index] = row;
            ci_prol_ex[index] = col;
            v_prol_ex[index]  = val;
            index++;

            // Next coarse node (top)
            col               = (j + 1) * ntheta_int / 4 + i / 2;
            val               = 0.5; // 1/2
            ri_prol_ex[index] = row;
            ci_prol_ex[index] = col;
            v_prol_ex[index]  = val;
            index++;
        }
    }
    // Coarse nodes in diagonals (2i+1, 2j+1)
    // => nb_nodes = ntheta_int * (nr_int / 2 - 2) / 2; TODO: Correct comment
    for (int j = 1; j <= nr_int - 1; j += 2) {
        for (int i = 1; i < ntheta_int; i += 2) {
            row = j * ntheta_int + i;

            // no top_right
            // top_left
            col               = (j + 1) * ntheta_int / 4 + (i - 1) / 2;
            val               = 0.5; // 1/2
            ri_prol_ex[index] = row;
            ci_prol_ex[index] = col;
            v_prol_ex[index]  = val;
            index++;
            // bottom_right
            col               = (j - 1) * ntheta_int / 4 + iplus_vect[i] / 2;
            val               = 0.5; // 1/2
            ri_prol_ex[index] = row;
            ci_prol_ex[index] = col;
            v_prol_ex[index]  = val;
            index++;
            // no bottom_left
        }
    }
} /* ----- end of level::build_prolongation_ex ----- */

/*!
 *  \brief Applies the extrapolation operator
 *
 * Applies the prolongation operator for implicit extrapolation.
 * 
 * The stencil is the same as the isotropic bilinear interpolation stencil,
 * except for fine points with diagonal coarse neighboors (types 5,6,7) where only
 * top_left and top_right values are used.
 *
 */
std::vector<double> level::apply_prolongation_ex(std::vector<double> u)
{
    std::vector<double> Pu;
    //Prolongation (P * u = Pu), (m x mc) * mc = m
    Pu = std::vector<double>(m, 0);

    std::vector<int> iplus_vect(ntheta_int);
    for (int i = 0; i < ntheta_int - 1; i++)
        iplus_vect[i] = i + 1;
    iplus_vect[ntheta_int - 1] = 0;

#pragma omp parallel shared(iplus_vect)
    {
#pragma omp single
        {
            // Coarse nodes
            // => nb_nodes = ntheta_int * (nr_int / 2 + 1) / 2;
            for (int j = 0; j < nr; j += 2) {
#pragma omp task firstprivate(j) shared(iplus_vect, Pu)
                {
                    int row, col;
                    double val;
                    for (int i = 0; i < ntheta_int; i += 2) {
                        row = j * ntheta_int + i;
                        col = j * ntheta_int / 4 + i / 2;
                        val = 1.0;
                        Pu[row] += val * u[col];
                    }
                } //end of task
            }
            // Coarse nodes in same r (2i, 2j+1)
            // => nb_nodes = ntheta_int * (nr_int / 2 + 1) / 2;
            for (int j = 0; j < nr; j += 2) {
#pragma omp task firstprivate(j) shared(iplus_vect, Pu)
                {
                    int row, col;
                    double val;
                    for (int i = 1; i < ntheta_int; i += 2) {
                        row = j * ntheta_int + i;

                        // Previous coarse node (left)
                        col = j * ntheta_int / 4 + (i - 1) / 2;
                        val = 0.5; // 1/2
                        Pu[row] += val * u[col];

                        // Next coarse node (right)
                        col = j * ntheta_int / 4 + iplus_vect[i] / 2;
                        val = 0.5; // 1/2
                        Pu[row] += val * u[col];
                    }
                } //end of task
            }
            // Coarse nodes in same theta (2i+1, 2j)
            // => nb_nodes = ntheta_int * (nr_int / 2 - 2) / 2; TODO: Correct comment
            for (int j = 1; j <= nr_int - 1; j += 2) {
#pragma omp task firstprivate(j) shared(iplus_vect, Pu)
                {
                    int row, col;
                    double val;
                    for (int i = 0; i < ntheta_int; i += 2) {
                        row = j * ntheta_int + i;

                        // Previous coarse node (bottom)
                        col = (j - 1) * ntheta_int / 4 + i / 2;
                        val = 0.5; // 1/2
                        Pu[row] += val * u[col];

                        // Next coarse node (top)
                        col = (j + 1) * ntheta_int / 4 + i / 2;
                        val = 0.5; // 1/2
                        Pu[row] += val * u[col];
                    }
                } //end of task
            }
            // Coarse nodes in diagonals (2i+1, 2j+1)
            // => nb_nodes = ntheta_int * (nr_int / 2 - 2) / 2; TODO: Correct comment
            for (int j = 1; j <= nr_int - 1; j += 2) {
#pragma omp task firstprivate(j) shared(iplus_vect, Pu)
                {
                    int row, col;
                    double val;
                    for (int i = 1; i < ntheta_int; i += 2) {
                        row = j * ntheta_int + i;

                        // no top_right
                        // top_left
                        col = (j + 1) * ntheta_int / 4 + (i - 1) / 2;
                        val = 0.5; // isotrop: 1/4
                        Pu[row] += val * u[col];

                        // bottom_right
                        col = (j - 1) * ntheta_int / 4 + iplus_vect[i] / 2;
                        val = 0.5; // isotrop: 1/4
                        Pu[row] += val * u[col];
                        // no bottom_left
                    }
                } //end of task
            }
        } //end of single
    } //end of parallel

    return Pu;
} /* ----- end of level::apply_prolongation_ex ----- */

/*!
 *  \brief Applies the extrapolation operator
 *
 * Applies the prolongation operator for implicit extrapolation.
 * 
 * The stencil is the same as the isotropic bilinear interpolation stencil,
 * except for fine points with diagonal coarse neighboors (types 5,6,7) where only
 * top_left and top_right values are used.
 *
 */
std::vector<double> level::apply_restriction_ex(std::vector<double> u)
{
    //Restriction (P^T * u = Pu), (mc x m) * m = mc
    std::vector<double> Pu = std::vector<double>(mc, 0);

    std::vector<int> iplus_vect(ntheta_int);
    for (int i = 0; i < ntheta_int - 1; i++)
        iplus_vect[i] = i + 1;
    iplus_vect[ntheta_int - 1] = 0;
    std::vector<int> imoins_vect(ntheta_int);
    for (int i = 1; i < ntheta_int; i++)
        imoins_vect[i] = i - 1;
    imoins_vect[0] = ntheta_int - 1;

#pragma omp parallel shared(iplus_vect)
    {
#pragma omp single
        {
            // // Loop through all coarse nodes
            int jc = 0;
#pragma omp task firstprivate(jc) shared(iplus_vect, imoins_vect, Pu)
            {
                int row, col, i, j, i_f, j_f;
                double val;
                for (int ic = 0; ic < ntheta_int / 2; ic++) {
                    col = jc * ntheta_int / 2 + ic;
                    i_f = ic * 2;
                    j_f = jc * 2;

                    // // Coarse nodes
                    i   = i_f;
                    j   = j_f;
                    row = j * ntheta_int + i;
                    val = 1.0;
                    Pu[col] += val * u[row];

                    // // Fine nodes in same r (2i, 2j+1)
                    // Next fine node (right)
                    i   = iplus_vect[i_f];
                    j   = j_f;
                    row = j * ntheta_int + i;
                    val = 0.5; // 1/2
                    Pu[col] += val * u[row];
                    // Previous fine node (left)
                    i   = imoins_vect[i_f];
                    j   = j_f;
                    row = j * ntheta_int + i;
                    val = 0.5; // 1/2
                    Pu[col] += val * u[row];

                    // // Fine nodes in same theta (2i+1, 2j)
                    // Next fine node (top)
                    i   = i_f;
                    j   = j_f + 1;
                    row = j * ntheta_int + i;
                    val = 0.5; // 1/2
                    Pu[col] += val * u[row];

                    // // Fine nodes in diagonals (2i+1, 2j+1)
                    // no bottom_left
                    // no bottom_right
                    // top_left
                    i   = imoins_vect[i_f];
                    j   = j_f + 1;
                    row = j * ntheta_int + i;
                    val = 0.5; // isotrop: 1/4
                    Pu[col] += val * u[row];
                    // no top_right
                }
            } //end of task
            jc = nr_int / 2;
#pragma omp task firstprivate(jc) shared(iplus_vect, imoins_vect, Pu)
            {
                int row, col, i, j, i_f, j_f;
                double val;
                for (int ic = 0; ic < ntheta_int / 2; ic++) {
                    col = jc * ntheta_int / 2 + ic;
                    i_f = ic * 2;
                    j_f = jc * 2;

                    // // Coarse nodes
                    i   = i_f;
                    j   = j_f;
                    row = j * ntheta_int + i;
                    val = 1.0;
                    Pu[col] += val * u[row];

                    // // Fine nodes in same r (2i, 2j+1)
                    // Next fine node (right)
                    i   = iplus_vect[i_f];
                    j   = j_f;
                    row = j * ntheta_int + i;
                    val = 0.5; // 1/2
                    Pu[col] += val * u[row];
                    // Previous fine node (left)
                    i   = imoins_vect[i_f];
                    j   = j_f;
                    row = j * ntheta_int + i;
                    val = 0.5; // 1/2
                    Pu[col] += val * u[row];

                    // // Fine nodes in same theta (2i+1, 2j)
                    // Previous fine node (bottom)
                    i   = i_f;
                    j   = j_f - 1;
                    row = j * ntheta_int + i;
                    val = 0.5; // 1/2
                    Pu[col] += val * u[row];

                    // // Fine nodes in diagonals (2i+1, 2j+1)
                    // no bottom_left
                    // bottom_right
                    i   = iplus_vect[i_f];
                    j   = j_f - 1;
                    row = j * ntheta_int + i;
                    val = 0.5; // isotrop: 1/4
                    Pu[col] += val * u[row];
                    // no top_left
                    // no top_right
                }
            } //end of task
            // => nb_nodes = ntheta_int * (nr_int / 2 + 1) / 2; TODO: Correct comment
            for (int jc = 1; jc < nr_int / 2; jc++) {
#pragma omp task firstprivate(jc) shared(iplus_vect, imoins_vect, Pu)
                {
                    int row, col, i, j, i_f, j_f;
                    double val;
                    for (int ic = 0; ic < ntheta_int / 2; ic++) {
                        col = jc * ntheta_int / 2 + ic;
                        i_f = ic * 2;
                        j_f = jc * 2;

                        // // Coarse nodes
                        i   = i_f;
                        j   = j_f;
                        row = j * ntheta_int + i;
                        val = 1.0;
                        Pu[col] += val * u[row];

                        // // Fine nodes in same r (2i, 2j+1)
                        // Next fine node (right)
                        i   = iplus_vect[i_f];
                        j   = j_f;
                        row = j * ntheta_int + i;
                        val = 0.5; // 1/2
                        Pu[col] += val * u[row];
                        // Previous fine node (left)
                        i   = imoins_vect[i_f];
                        j   = j_f;
                        row = j * ntheta_int + i;
                        val = 0.5; // 1/2
                        Pu[col] += val * u[row];

                        // // Fine nodes in same theta (2i+1, 2j)
                        // Next fine node (top)
                        i   = i_f;
                        j   = j_f + 1;
                        row = j * ntheta_int + i;
                        val = 0.5; // 1/2
                        Pu[col] += val * u[row];
                        // Previous fine node (bottom)
                        i   = i_f;
                        j   = j_f - 1;
                        row = j * ntheta_int + i;
                        val = 0.5; // 1/2
                        Pu[col] += val * u[row];

                        // // Fine nodes in diagonals (2i+1, 2j+1)
                        // no bottom_left
                        // bottom_right
                        i   = iplus_vect[i_f];
                        j   = j_f - 1;
                        row = j * ntheta_int + i;
                        val = 0.5; // isotrop: 1/4
                        Pu[col] += val * u[row];
                        // top_left
                        i   = imoins_vect[i_f];
                        j   = j_f + 1;
                        row = j * ntheta_int + i;
                        val = 0.5; // isotrop: 1/4
                        Pu[col] += val * u[row];
                        // no top_right
                    }
                } //end of task
            }
        } //end of single
    } //end of parallel

    return Pu;
} /* ----- end of level::apply_restriction_ex ----- */
