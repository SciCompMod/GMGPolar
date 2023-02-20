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
 * \file direct_solve.cpp
 * \brief Implementation of the direct solver for the coarse grid operator and the smoother
 * \author M. Kuehn, C. Kruse, P. Leleux, C. Schwarz
 * \version 0.0
 */
#include "level.h"

#ifdef USE_MUMPS
/*!
 *  \brief Initialization of a MUMPS structure
 *
 * Initialize the MUMPS structure: calls dmumps_c with job=-1
 *
 * \param mumps: the structure
 */
void level::init_mumps(DMUMPS_STRUC_C& mumps)
{
    mumps.sym          = 0; // symmetric
    mumps.par          = 1; // working host
    mumps.comm_fortran = -987654; // MPI_COMM_WORLD
    mumps.job          = -1; // Initialization

    dmumps_c(&mumps);
} /* ----- end of level::init_mumps ----- */

/*!
 *  \brief Analysis and factorization with MUMPS (no RHS)
 *
 * Analysis and factorization with MUMPS (no RHS): calls dmumps_c with job=4
 *
 * \param mumps: the structure
 * \param A_row_indices: vector of row indices for the matrix
 * \param A_col_indices: vector of column indices for the matrix
 * \param A_vals: vector of valus for the matrix
 * \param m_solution: size of the solutions in MUMPS
 */
void level::facto_mumps(DMUMPS_STRUC_C& mumps, std::vector<int> A_row_indices, std::vector<int> A_col_indices,
                        std::vector<double> A_vals, int m_solution)
{
    int nz    = A_row_indices.size();
    mumps.job = 4;
    mumps.n   = m_solution;
    mumps.nz  = nz;
    mumps.irn = new int[nz];
    mumps.jcn = new int[nz];
    mumps.a   = new double[nz];
    for (int j = 0; j < nz; j++) {
        mumps.irn[j] = A_row_indices[j] + 1;
        mumps.jcn[j] = A_col_indices[j] + 1;
        mumps.a[j]   = A_vals[j];
    }
    mumps.lrhs     = mumps.n;
    mumps.nrhs     = 1;
    mumps.rhs      = new double[mumps.lrhs];
    mumps.ICNTL(1) = 0;
    mumps.ICNTL(2) = 0;
    mumps.ICNTL(3) = 0;
    mumps.ICNTL(4) = 0;
    mumps.ICNTL(7) = 5;
    dmumps_c(&(mumps));
    mumps.job = 3;
} /* ----- end of level::facto_mumps ----- */

/*!
 *  \brief Solve with MUMPS on a specific RHS
 *
 * Solve with MUMPS on a specific RHS: calls dmumps_c with job=3
 *
 * \param mumps: the structure
 * \param f: the RHS
 */
std::vector<double> level::solve_mumps(DMUMPS_STRUC_C& mumps, std::vector<double> f)
{
    std::vector<double> sol(mumps.lrhs);

    for (int i = 0; i < mumps.lrhs; i++) {
        mumps.rhs[i] = f[i];
    }
    dmumps_c(&(mumps));
    for (int i = 0; i < mumps.lrhs; i++) {
        sol[i] = mumps.rhs[i];
    }

    return sol;
} /* ----- end of level::solve_mumps ----- */

/*!
 *  \brief Finalization of a MUMPS structure
 *
 * Finalization the MUMPS structure: calls dmumps_c with job=-2
 *
 * \param mumps: the structure
 */
void level::finalize_mumps(DMUMPS_STRUC_C& mumps)
{
    int job_cur    = mumps.job;
    mumps.job      = -2;
    mumps.ICNTL(1) = 0;
    mumps.ICNTL(2) = 0;
    mumps.ICNTL(3) = 0;
    mumps.ICNTL(4) = 0;
    dmumps_c(&mumps);
    if (job_cur == 3) {
        delete[] mumps.a;
        delete[] mumps.irn;
        delete[] mumps.jcn;
        delete[] mumps.rhs;
    }
} /* ----- end of level::finalize_mumps ----- */
#endif

/*!
 *  \brief Factorization of A using an in-house sparse solver
 *
 * Factorization of A using an in-house sparse solver
 *
 * \param A_row_indices: vector of row indices for the matrix
 * \param A_col_indices: vector of column indices for the matrix
 * \param A_vals: vector of valus for the matrix
 * \param m_solution: size of the solutions in MUMPS
 */
void level::facto_gaussian_elimination(std::vector<int>& A_row_indices, std::vector<int>& A_col_indices,
                                       std::vector<double>& A_vals, int m_solution)

{
    // Factorization with Gaussian elimination
    int row = A_row_indices[0];
    double pivot;
    std::vector<int> cols;
    std::vector<double> vals;
    std::vector<int> rows_L;
    std::vector<double> vals_L;

    std::vector<int> A_row_ptr(m_solution + 1, 0);
    int count = 1;
    for (std::size_t k = 0; k < A_row_indices.size() + 1; ++k) {
        if (k == A_row_indices.size() || A_row_indices[k] != row) {
            A_row_ptr[count] = k;
            count++;
            if (k != A_row_indices.size())
                row = A_row_indices[k];
        }
    }

    row = A_row_indices[0];
    // Go through the list and accumulate the elements of a same row (assume increasing rows)
    for (std::size_t k = 0; k < A_row_indices.size() + 1; ++k) {

        // If new row, construct col of L and row of U
        if (k == A_row_indices.size() || row != A_row_indices[k]) {
            // Update the column of L with L(j,k) = U(j,k) / U(k,k)
            // i.e. to save storage: A(j,k) = A(j,k) / A(k,k), j > k
            for (std::size_t j = k; j < A_row_indices.size(); ++j) {
                // col == k and row > col
                if (A_col_indices[j] == row && A_row_indices[j] > A_col_indices[j]) {
                    A_vals[j] /= pivot;
                    rows_L.push_back(A_row_indices[j]);
                    vals_L.push_back(A_vals[j]);
                }
            }

            // Update already existing entries
            std::vector<int> add_rows, add_cols, add_index;
            std::vector<double> add_vals;
            for (std::size_t j = k; j < A_row_indices.size(); ++j) {
                if (A_row_indices[j] < row)
                    throw std::runtime_error("Loop through previous row in Gaussian elimination, should not happen if "
                                             "rows sorted in increasing order in the matrix.");

                // If row == k or col <= k, continue searching
                if (A_row_indices[j] == row || A_col_indices[j] <= row)
                    continue;

                //  col not in pivot row (U_ki = 0) or row not in L_k (L_jk = 0), continue searching
                auto col_ind = std::find(cols.begin(), cols.end(), A_col_indices[j]);
                auto row_ind = std::find(rows_L.begin(), rows_L.end(), A_row_indices[j]);
                if (col_ind == cols.end() || row_ind == rows_L.end())
                    continue;
                // Else, update the row of U with U(j,i) = U(j,i) - L(j,k) * U(k,i)
                // i.e. to save storage: A(j,i) = A(j,i) - A(j,k) * A(k,i), j > k
                A_vals[j] -= vals_L[row_ind - rows_L.begin()] * vals[col_ind - cols.begin()];
                add_rows.push_back(A_row_indices[j]);
                add_cols.push_back(A_col_indices[j]);
                add_vals.push_back(A_vals[j]);
            }

            // Take fill-in into account: values at the intersection of entries L_jk and U_ki
            std::vector<int> add_rows_FI_tot, add_cols_FI_tot;
            std::vector<double> add_vals_FI_tot;
            for (std::size_t j = 0; j < rows_L.size(); ++j) {
                std::vector<int> add_rows_FI, add_cols_FI;
                std::vector<double> add_vals_FI;
                for (std::size_t i = 0; i < cols.size(); ++i) {
                    // Only update the upper part of U
                    if (cols[i] <= row || rows_L[j] <= row)
                        continue;

                    // Check if entry exists: was already updated
                    int exists = 0;
                    for (std::size_t k = 0; k < add_rows.size(); k++)
                        if (add_rows[k] == rows_L[j] && add_cols[k] == cols[i]) {
                            exists = 1;
                            break;
                        }
                    if (exists)
                        continue;

                    // Else, update the row of U with U(j,i) = U(j,i) - L(j,k) * U(k,i)
                    // i.e. to save storage: A(j,i) = A(j,i) - A(j,k) * A(k,i), j > k
                    add_rows_FI.push_back(rows_L[j]);
                    add_cols_FI.push_back(cols[i]);
                    add_vals_FI.push_back(-vals_L[j] * vals[i]);
                }
                A_row_indices.insert(A_row_indices.begin() + A_row_ptr[rows_L[j] + 1], add_rows_FI.begin(),
                                     add_rows_FI.end());
                A_col_indices.insert(A_col_indices.begin() + A_row_ptr[rows_L[j] + 1], add_cols_FI.begin(),
                                     add_cols_FI.end());
                A_vals.insert(A_vals.begin() + A_row_ptr[rows_L[j] + 1], add_vals_FI.begin(), add_vals_FI.end());

                for (std::size_t r = rows_L[j] + 1; r < A_row_ptr.size(); r++)
                    A_row_ptr[r] += add_rows_FI.size();
            }

            // New row:
            if (k < A_row_indices.size()) {
                row = A_row_indices[k];
                // - empty columns/rows and values for pivot row and L_k
                cols.clear();
                vals.clear();
                rows_L.clear();
                vals_L.clear();
                // - loop again on this entry
                if (k != A_row_indices.size())
                    k--;
            }
        }
        // If same row, store the values and col indices
        else {
            // Do not consider values from L
            if (A_col_indices[k] < row)
                continue;
            if (A_col_indices[k] == row) {
                pivot = A_vals[k];
                if (pivot == 0)
                    throw std::runtime_error("Null pivot found in Gaussian elimination.");
            }
            cols.push_back(A_col_indices[k]);
            vals.push_back(A_vals[k]);
        }
    }
} /* ----- end of level::facto_gaussian_elimination ----- */

/*!
 *  \brief Solve on a specific RHS using an in-house sparse solver
 *
 * Solve on a specific RHS using an in-house sparse solver
 *
 * \param A_row_indices: vector of row indices for the matrix
 * \param A_col_indices: vector of column indices for the matrix
 * \param A_vals: vector of valus for the matrix
 * \param f: the RHS
 */
std::vector<double> level::solve_gaussian_elimination(std::vector<int> A_row_indices, std::vector<int> A_col_indices,
                                                      std::vector<double> A_vals, std::vector<double> f)
{
    // int m_solution = f.size();

    // std::cout << "Forward substitution \n";
    // Solve Ly = b for y (Forward Substitution)
    // (assumes rows in increasing order)
    int row = A_row_indices[0];
    std::vector<double> y(f);
    for (std::size_t j = 0; j < A_row_indices.size(); ++j) {
        if (A_col_indices[j] >= A_row_indices[j])
            continue;
        y[A_row_indices[j]] -= y[A_col_indices[j]] * A_vals[j];
    }

    // std::cout << "Backward substitution \n";
    //Solve Ux = y for x (Backward Substitution)
    row = A_row_indices[A_row_indices.size() - 1];
    int jind;
    double U_jj;
    std::vector<double> x(y);
    for (int j = A_row_indices.size() - 1; j >= -1; --j) {
        if (j != -1 && A_col_indices[j] < A_row_indices[j])
            continue;

        // If new row, we can apply the diagonal element of U
        if (j == -1 || A_row_indices[j] != row) {
            x[row] /= U_jj;

            if (j == -1)
                break;

            row = A_row_indices[j];
        }

        if (A_row_indices[j] == A_col_indices[j]) {
            U_jj = A_vals[j];
            jind = j;
            continue;
        }

        x[row] -= x[A_col_indices[j]] * A_vals[j];
    }

    return x;

} /* ----- end of level::solve_gaussian_elimination ----- */

/*!
 *  \brief Factorization and solve on a specific RHS using an in-house sparse solver (deprecated)
 *
 * Factorization and solve on a specific RHS using an in-house sparse solver
 *
 * \param A_row_indices: vector of row indices for the matrix
 * \param A_col_indices: vector of column indices for the matrix
 * \param A_vals: vector of valus for the matrix
 * \param f: the RHS
 */
std::vector<double> level::solve_gaussian_elimination_fb_subst(std::vector<int> A_row_indices,
                                                               std::vector<int> A_col_indices,
                                                               std::vector<double> A_vals, std::vector<double> f)
{
    //std::cout << "Gaussian Elimination!\n";

    int m_solution = f.size();
    std::vector<double> x(m_solution, 0.0); //x is the solution vector, it has the same size than f, fill with zeros

    //Gaussian Elimination
    //create matrices U and L
    //U = A
    std::vector<int> U_row_indices(A_row_indices);
    std::vector<int> U_col_indices(A_col_indices);
    std::vector<double> U_vals(A_vals);

    //L = I
    std::vector<int> L_row_indices; //=[0,1,2,3,4...m-1]
    std::vector<int> L_col_indices; //=[0,1,2,3,4...m-1]
    std::vector<double> L_vals; //=[1,1,1,1...] m-times
    for (int i = 0; i < m_solution; ++i) {
        L_col_indices.push_back(i);
        L_row_indices.push_back(i);
        L_vals.push_back(1);
    }

    //Gaussian Elimination
    for (int k = 0; k < m_solution - 1; ++k) { //k=col index
        for (int j = k + 1; j < m_solution; ++j) { //j=row index

            double U_jk = get_element(U_row_indices, U_col_indices, U_vals, j, k); //get the element (j,k) of U

            //if the U_jk is zero, we don't need to treat it and can just go on with the next element
            if (U_jk != 0) {
                double U_kk = get_element(U_row_indices, U_col_indices, U_vals, k, k); //get the element (k,k) of U
                double L_jk = U_jk / U_kk; //step of GE

                //insert the new element L_jk into the matrix L
                //as it has to be a new value, which is not yet there, just append new value  to the end
                L_vals.push_back(L_jk);
                L_row_indices.push_back(j);
                L_col_indices.push_back(k);

                for (int i = k; i < m_solution; ++i) {

                    double U_ki = get_element(U_row_indices, U_col_indices, U_vals, k, i); //get the element (k,i) of U

                    if (U_ki != 0) { //we can skip this part, if U_ki=0
                        //U(j,i) = U(j,i) - L(j,k) * U(k,i)    --> step of GE
                        double value = get_element(U_row_indices, U_col_indices, U_vals, j, i) - L_jk * U_ki;
                        set_element(U_row_indices, U_col_indices, U_vals, j, i, value);
                    }
                }
            }
        }
    }

    //std::cout << "backward substitution \n";
    //Forward-Backward-Substitution

    //Solve Ly=b for y (Forward Substitution)
    std::vector<double> y(m_solution, 0.0);
    for (int j = 0; j < m_solution; ++j) {

        double L_jj = get_element(L_row_indices, L_col_indices, L_vals, j, j); //get the element (j,j) of U
        double sum  = 0;

        for (int k = 0; k < j; ++k) {
            double L_jk = get_element(L_row_indices, L_col_indices, L_vals, j, k); //get the element (j,k) of U
            sum += y[k] * L_jk;
        }
        y[j] = (f[j] - sum) / L_jj;
    }

    //Solve Ux=y for x (Backward Substitution)
    for (int j = m_solution - 1; j >= 0; --j) {

        double U_jj = get_element(U_row_indices, U_col_indices, U_vals, j, j); //get the element (j,j) of U
        double sum  = 0;

        for (int k = j; k < m_solution; ++k) {
            double U_jk = get_element(U_row_indices, U_col_indices, U_vals, j, k); //get the element (j,k) of U
            sum += x[k] * U_jk;
        }
        x[j] = (y[j] - sum) / U_jj;
    }

    return x;
} /* ----- end of level::solve_gaussian_elimination_fb_subst ----- */

/*!
 *  \brief Finds value of entry (row, col) in the sparse matrix (deprecated)
 *
 * Finds value of entry (row, col) in the sparse matrix
 *
 * \param A_row_indices: vector of row indices for the matrix
 * \param A_col_indices: vector of column indices for the matrix
 * \param A_vals: vector of valus for the matrix
 * \param row_index: row to find
 * \param col_index: col to find
 */
double level::get_element(std::vector<int> A_row_indices, std::vector<int> A_col_indices, std::vector<double> A_vals,
                          int row_index, int col_index)
{
    //get the element (i,j) of A
    double A_ij = 0;

    //second version of code
    int min   = 0;
    int max   = A_vals.size();
    int count = log2(A_vals.size() / 10);

    //try to make the search region smaller, values in matrix are sorted by ascending rows
    for (int z = 0; z < count; ++z) {
        int middle = min + (max - min) / 2;
        if (row_index > middle) {
            min = middle;
        }
        else {
            max = middle;
        }
    }

    for (int v = min; v < (int)A_vals.size(); ++v) {
        if (A_row_indices[v] == row_index) {
            if (A_col_indices[v] == col_index) {
                A_ij = A_vals[v];
                break; //break the for loop when we found the value
            }
        }
    }

    return A_ij;
} /* ----- end of level::get_element ----- */

/*!
 *  \brief Sets value of entry (row, col) in the sparse matrix (deprecated)
 *
 * Sets value of entry (row, col) in the sparse matrix
 *
 * \param A_row_indices: vector of row indices for the matrix
 * \param A_col_indices: vector of column indices for the matrix
 * \param A_vals: vector of valus for the matrix
 * \param row_index: row to find
 * \param col_index: col to find
 */
void level::set_element(std::vector<int>& A_row_indices, std::vector<int>& A_col_indices, std::vector<double>& A_vals,
                        int row_index, int col_index, double value)
{

    //set the element (i,j) of A to 'value'
    //if there exists already an element, it is overwritten with 'value'
    //if the element does not exist yet, just append it to the end

    int min   = 0;
    int max   = A_vals.size();
    int count = log2(A_vals.size() / 10);

    //try to make the search region smaller, values in matrix are sorted by ascending rows
    for (int z = 0; z < count; ++z) {
        int middle = min + (max - min) / 2;
        if (row_index > middle) {
            min = middle;
        }
        else {
            max = middle;
        }
    }

    int existing_element = 0;
    for (int v = min; v < (int)A_vals.size(); ++v) { //check if there exists already an element
        if (A_row_indices[v] == row_index) {
            if (A_col_indices[v] == col_index) {
                existing_element = 1;
                if (value != 0) { //if value is not zero: set it
                    A_vals[v] = value;
                }
                else { //if value is zero: delete element from sparse matrix
                    A_vals.erase(A_vals.begin() + v);
                    A_row_indices.erase(A_row_indices.begin() + v);
                    A_col_indices.erase(A_col_indices.begin() + v);
                }
                break; //if we found the element, we can break the for-loop
            }
        }
    }
    if (existing_element == 0) { //if there exists no element yet, just append it
        A_vals.push_back(value);
        A_row_indices.push_back(row_index);
        A_col_indices.push_back(col_index);
    }
} /* ----- end of level::set_element ----- */

/*!
 *  \brief Solve diagonal system on a specific RHS
 *
 * Solve diagonal system on a specific RHS
 *
 * \param A_vals: vector of values for the diagonal matrix
 * \param f: the RHS
 */
std::vector<double> level::solve_diag(std::vector<double> A_vals, std::vector<double> f)
{
    int m = f.size();
    std::vector<double> x(f);
    // #pragma omp parallel for shared(x, A_vals)
    for (int i = 0; i < m; i++) {
        x[i] /= A_vals[i];
    }
    return x;
} /* ----- end of level::solve_gaussian_elimination_fb_subst ----- */

/*!
 *  \brief Solve system for Circle smoother on a specific RHS
 *
 * Solve system for Circle smoother on a specific RHS
 * The size of the system should be > 2
 *
 * \param A_row_indices: vector of row indices for the matrix
 * \param A_col_indices: vector of column indices for the matrix
 * \param A_vals: vector of valus for the matrix
 * \param f: the RHS
 */
std::vector<double> level::solve_circle(std::vector<int> A_row_indices, std::vector<int> A_col_indices,
                                        std::vector<double> A_vals, std::vector<double> f)
{
    int m = f.size();
    std::vector<double> x(f);

    /**************************************
     * Forward substitution (Ly=f)
     *************************************/
    // interior
    for (int i = 0; i < m - 2; i++) {
        int ind = 4 * i + 3;
        x[i + 1] -= A_vals[ind] * x[i];
        int ind2 = 4 * (m - 3) + 3 + 3;
        x[m - 1] -= A_vals[ind2 + i] * x[i];
    }
    // Penultian line
    int i    = m - 2;
    int ind2 = 4 * (m - 3) + 3 + 3;
    x[m - 1] -= A_vals[ind2 + i] * x[i];

    /**************************************
     * Backward substitution (Ux=y)
     *************************************/
    // Last line
    int ind = 4 * (m - 3) + 3 + 3;
    x[m - 1] /= A_vals[ind + m - 1];
    ind = 4 * (m - 3) + 3;
    x[m - 2] -= A_vals[ind + 2] * x[m - 1];
    // Penultian line
    i   = m - 3;
    ind = 4 * i + 3;
    x[i + 1] /= A_vals[ind + 1];
    x[i] -= A_vals[ind - 2] * x[i + 1];

    // Interior
    for (int i = m - 4; i >= 0; i--) {
        ind = 4 * i + 3;
        x[i + 1] -= A_vals[ind + 3] * x[m - 1];
        x[i + 1] /= A_vals[ind + 1];
        x[i] -= A_vals[ind - 2] * x[i + 1];
    }
    // first
    x[0] -= A_vals[2] * x[m - 1];
    x[0] /= A_vals[0];

    return x;
} /* ----- end of level::solve_gaussian_elimination_fb_subst ----- */

/*!
 *  \brief Solve system for Circle smoother on a specific RHS
 *
 * Solve system for Circle smoother on a specific RHS
 * The size of the system should be > 2
 *
 * \param A_row_indices: vector of row indices for the matrix
 * \param A_col_indices: vector of column indices for the matrix
 * \param A_vals: vector of valus for the matrix
 * \param f: the RHS
 */
std::vector<double> level::solve_radial(std::vector<int> A_row_indices, std::vector<int> A_col_indices,
                                        std::vector<double> A_vals, std::vector<double> f)
{
    int i, ind, m = f.size();
    std::vector<double> x(f);

    /**************************************
     * Forward substitution (Ly=f)
     *************************************/
    // interior
    for (i = 0; i < m - 2; i++) {
        ind = 3 * i + 2;
        x[i + 1] -= A_vals[ind] * x[i];
    }

    /**************************************
     * Backward substitution (Ux=y)
     *************************************/
    ind = 3 * (m - 3) + 4;
    x[m - 1] /= A_vals[ind];
    for (i = m - 3; i >= 0; i--) {
        ind = 3 * i + 2;
        x[i + 1] /= A_vals[ind + 1];
        x[i] -= A_vals[ind - 1] * x[i + 1];
    }
    // first lines
    x[0] /= A_vals[0];

    return x;
} /* ----- end of level::solve_gaussian_elimination_fb_subst ----- */

/*!
 *  \brief Factorization of A using an in-house sparse solver
 *
 * Factorization of A using an in-house sparse solver
 *
 * \param A_row_indices: vector of row indices for the matrix
 * \param A_col_indices: vector of column indices for the matrix
 * \param A_vals: vector of valus for the matrix
 * \param m_solution: size of the solutions in MUMPS
 */
void level::facto_radial(std::vector<int>& A_row_indices, std::vector<int>& A_col_indices, std::vector<double>& A_vals,
                         int m)

{
    int i, ind;

    // Interior
    for (i = 0; i < m - 2; i++) {
        ind = 3 * i + 2;
        A_vals[ind] /= A_vals[ind - 2]; // L
        A_vals[ind + 1] -= A_vals[ind - 1] * A_vals[ind]; // U
    }
} /* ----- end of level::facto_gaussian_elimination ----- */
/*!
 *  \brief Factorization of A using an in-house sparse solver
 *
 * Factorization of A using an in-house sparse solver
 *
 * \param A_row_indices: vector of row indices for the matrix
 * \param A_col_indices: vector of column indices for the matrix
 * \param A_vals: vector of valus for the matrix
 * \param m_solution: size of the solutions in MUMPS
 */
void level::facto_circle(std::vector<int>& A_row_indices, std::vector<int>& A_col_indices, std::vector<double>& A_vals,
                         int m)

{
    int i, ind, ind2;

    // First line
    ind = 3;
    A_vals[ind] /= A_vals[ind - 3]; // L
    A_vals[ind + 1] -= A_vals[ind - 2] * A_vals[ind]; // U (diag)
    A_vals[ind + 3] -= A_vals[ind - 1] * A_vals[ind]; // U (last col=fill-in)
    // Update last line
    ind2 = 4 * (m - 3) + 3 + 3;
    A_vals[ind2] /= A_vals[ind - 3]; // L
    A_vals[ind2 + 1] -= A_vals[ind - 2] * A_vals[ind2]; // L (last line=fill-in)
    A_vals[ind2 + m - 1] -= A_vals[ind - 1] * A_vals[ind2]; // U (diag)
    // Interior
    for (i = 1; i < m - 3; i++) {
        ind = 4 * i + 3;
        A_vals[ind] /= A_vals[ind - 3]; // L
        A_vals[ind + 1] -= A_vals[ind - 2] * A_vals[ind]; // U (diag)
        A_vals[ind + 3] -= A_vals[ind - 1] * A_vals[ind]; // U (last col=fill-in)
        // Update last line
        ind2 = 4 * (m - 3) + 3 + 3;
        A_vals[ind2 + i] /= A_vals[ind - 3]; // L
        A_vals[ind2 + i + 1] -= A_vals[ind - 2] * A_vals[ind2 + i]; // L (last line=fill-in)
        A_vals[ind2 + m - 1] -= A_vals[ind - 1] * A_vals[ind2 + i]; // U (diag)
    }
    // Penultian line
    i   = (m - 3);
    ind = 4 * i + 3;
    A_vals[ind] /= A_vals[ind - 3]; // L
    A_vals[ind + 1] -= A_vals[ind - 2] * A_vals[ind]; // U (diag)
    A_vals[ind + 2] -= A_vals[ind - 1] * A_vals[ind]; // U (last col)
    // Update last line
    ind2 = 4 * (m - 3) + 3 + 3;
    A_vals[ind2 + i] /= A_vals[ind - 3]; // L
    A_vals[ind2 + i + 1] -= A_vals[ind - 2] * A_vals[ind2 + i]; // U (diag)
    A_vals[ind2 + m - 1] -= A_vals[ind - 1] * A_vals[ind2 + i]; // U (last col)
    // Last line
    A_vals[ind2 + m - 2] /= A_vals[ind + 1]; // L
    A_vals[ind2 + m - 1] -= A_vals[ind + 2] * A_vals[ind2 + m - 2]; // L
} /* ----- end of level::facto_gaussian_elimination ----- */

/*!
 *  \brief Initialize the LU factors for circle
 *
 *  Initialize the LU factors for the circle smoother based on the matrix
 *  while taking into account the fill-in.
 *
 *  \param ij: block number for the current smoother
 *  \param smoother: the current smoother
 *
 */
void level::fill_in_circle(int ij, int smoother)
{
    int ind, ind2;
    int msc     = m_sc[smoother];
    int size_LU = 3 * msc + 2 * (msc - 3);
    // std::cout << "ij: " << ij << ", smoother: " << smoother << ", m_sc: " << msc << ", size_LU: " << size_LU << "\n";
    // std::cout << "A_Zebra_r_row[smoother][ij].size(): " << A_Zebra_r_row[smoother][ij].size() << "\n";
    // first line
    A_Zebra_r_LU_row[smoother][ij]    = std::vector<int>(size_LU);
    A_Zebra_c_LU_row[smoother][ij]    = std::vector<int>(size_LU);
    A_Zebra_v_LU_row[smoother][ij]    = std::vector<double>(size_LU);
    A_Zebra_r_LU_row[smoother][ij][0] = A_Zebra_r_row[smoother][ij][1];
    A_Zebra_c_LU_row[smoother][ij][0] = A_Zebra_c_row[smoother][ij][1];
    A_Zebra_v_LU_row[smoother][ij][0] = A_Zebra_v_row[smoother][ij][1];
    A_Zebra_r_LU_row[smoother][ij][1] = A_Zebra_r_row[smoother][ij][2];
    A_Zebra_c_LU_row[smoother][ij][1] = A_Zebra_c_row[smoother][ij][2];
    A_Zebra_v_LU_row[smoother][ij][1] = A_Zebra_v_row[smoother][ij][2];
    A_Zebra_r_LU_row[smoother][ij][2] = A_Zebra_r_row[smoother][ij][0];
    A_Zebra_c_LU_row[smoother][ij][2] = A_Zebra_c_row[smoother][ij][0];
    A_Zebra_v_LU_row[smoother][ij][2] = A_Zebra_v_row[smoother][ij][0];
    // interior
    for (int i = 0; i < msc - 2; i++) {
        ind  = 3 * i + 3;
        ind2 = 4 * i + 3;
        // std::cout << "ind: " << ind << ", ind2: " << ind2 << ", ij: " << ij << "\n";
        // std::cout << "A_Zebra_r_row[smoother][ij][ind]: " << A_Zebra_r_row[smoother][ij][ind] << "\n";
        // std::cout << "A_Zebra_r_LU_row[smoother][ij][ind2]: " << A_Zebra_r_LU_row[smoother][ij][ind2] << "\n";
        A_Zebra_r_LU_row[smoother][ij][ind2]     = A_Zebra_r_row[smoother][ij][ind];
        A_Zebra_r_LU_row[smoother][ij][ind2 + 1] = A_Zebra_r_row[smoother][ij][ind + 1];
        A_Zebra_r_LU_row[smoother][ij][ind2 + 2] = A_Zebra_r_row[smoother][ij][ind + 2];
        A_Zebra_r_LU_row[smoother][ij][ind2 + 3] = A_Zebra_r_row[smoother][ij][ind + 2];
        A_Zebra_c_LU_row[smoother][ij][ind2]     = A_Zebra_c_row[smoother][ij][ind];
        A_Zebra_c_LU_row[smoother][ij][ind2 + 1] = A_Zebra_c_row[smoother][ij][ind + 1];
        A_Zebra_c_LU_row[smoother][ij][ind2 + 2] = A_Zebra_c_row[smoother][ij][ind + 2];
        A_Zebra_c_LU_row[smoother][ij][ind2 + 3] = msc - 1;
        A_Zebra_v_LU_row[smoother][ij][ind2]     = A_Zebra_v_row[smoother][ij][ind];
        A_Zebra_v_LU_row[smoother][ij][ind2 + 1] = A_Zebra_v_row[smoother][ij][ind + 1];
        A_Zebra_v_LU_row[smoother][ij][ind2 + 2] = A_Zebra_v_row[smoother][ij][ind + 2];
        A_Zebra_v_LU_row[smoother][ij][ind2 + 3] = 0;
    }
    // Penultian line
    ind                                      = 3 * (msc - 3) + 3;
    ind2                                     = 4 * (msc - 3) + 3;
    A_Zebra_r_LU_row[smoother][ij][ind2]     = A_Zebra_r_row[smoother][ij][ind];
    A_Zebra_r_LU_row[smoother][ij][ind2 + 1] = A_Zebra_r_row[smoother][ij][ind + 1];
    A_Zebra_r_LU_row[smoother][ij][ind2 + 2] = A_Zebra_r_row[smoother][ij][ind + 2];
    A_Zebra_c_LU_row[smoother][ij][ind2]     = A_Zebra_c_row[smoother][ij][ind];
    A_Zebra_c_LU_row[smoother][ij][ind2 + 1] = A_Zebra_c_row[smoother][ij][ind + 1];
    A_Zebra_c_LU_row[smoother][ij][ind2 + 2] = A_Zebra_c_row[smoother][ij][ind + 2];
    A_Zebra_v_LU_row[smoother][ij][ind2]     = A_Zebra_v_row[smoother][ij][ind];
    A_Zebra_v_LU_row[smoother][ij][ind2 + 1] = A_Zebra_v_row[smoother][ij][ind + 1];
    A_Zebra_v_LU_row[smoother][ij][ind2 + 2] = A_Zebra_v_row[smoother][ij][ind + 2];
    // last line
    ind                                  = 3 * (msc - 3) + 3 + 3;
    ind2                                 = 4 * (msc - 3) + 3 + 3;
    A_Zebra_r_LU_row[smoother][ij][ind2] = A_Zebra_r_row[smoother][ij][ind + 2];
    A_Zebra_c_LU_row[smoother][ij][ind2] = A_Zebra_c_row[smoother][ij][ind + 2];
    A_Zebra_v_LU_row[smoother][ij][ind2] = A_Zebra_v_row[smoother][ij][ind + 2];
    for (int i = 1; i < msc - 2; i++) {
        A_Zebra_r_LU_row[smoother][ij][ind2 + i] = msc - 1;
        A_Zebra_c_LU_row[smoother][ij][ind2 + i] = i;
        A_Zebra_v_LU_row[smoother][ij][ind2 + i] = 0;
    }
    A_Zebra_r_LU_row[smoother][ij][ind2 + msc - 2] = A_Zebra_r_row[smoother][ij][ind];
    A_Zebra_c_LU_row[smoother][ij][ind2 + msc - 2] = A_Zebra_c_row[smoother][ij][ind];
    A_Zebra_v_LU_row[smoother][ij][ind2 + msc - 2] = A_Zebra_v_row[smoother][ij][ind];
    A_Zebra_r_LU_row[smoother][ij][ind2 + msc - 1] = A_Zebra_r_row[smoother][ij][ind + 1];
    A_Zebra_c_LU_row[smoother][ij][ind2 + msc - 1] = A_Zebra_c_row[smoother][ij][ind + 1];
    A_Zebra_v_LU_row[smoother][ij][ind2 + msc - 1] = A_Zebra_v_row[smoother][ij][ind + 1];
} /* ----- end of level::fill_in_circle ----- */
