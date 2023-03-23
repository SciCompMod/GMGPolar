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
 * \file cuda.cu
 * \brief Implementation of the CUDA functions
 * \author M. Kuehn, C. Kruse, P. Leleux, C. Schwarz
 * \version 0.0
 */

#ifdef CUDA

    #include "gmgpolar.h"
    #include "cuda.h"

    #define NTHREADS_CUDA 256

/*!
 *  \brief Computes the residual
 *
 *  Computes the residual based on the approximation res = b-Au
 *
 * \param l: the level (0=finest)
 *
 */
void gmgpolar::compute_residual(int l)
{
    //printf("COMPUTE RESIDUAL *************************************** \n");

    v_level[l]->res = std::vector<double>(v_level[l]->m);

    int n_threads = 32;
    int n_blocks  = ceil((double)v_level[l]->m / (double)n_threads);

    compute_residual_gpu<<<n_blocks, n_threads>>>(
        v_level[l]->m, v_level[l]->dev_u, v_level[l]->dev_Au, v_level[l]->dev_fVec, v_level[l]->dev_res,
        v_level[l]->nr_int, v_level[l]->dev_r, v_level[l]->dev_hplus, v_level[l]->ntheta_int, v_level[l]->dev_theta,
        v_level[l]->dev_thetaplus, gyro::icntl[Param::DirBC_Interior], gyro::icntl[Param::mod_pk],
        gyro::icntl[Param::prob], gyro::dcntl[Param::kappa_eps], gyro::dcntl[Param::delta_e], gyro::dcntl[Param::R]);

    cudaMemcpy((void*)&(v_level[l]->res[0]), (void*)v_level[l]->dev_res, v_level[l]->m * sizeof(double),
               cudaMemcpyDeviceToHost);

    // TODO: return 2-norm of residual

} /* ----- end of gmgpolar::compute_residual ----- */
/*!
 *  \brief Initializes to 0
 *
 *  Initializes to 0
 * 
 * \param m: the size of the vector
 * \param vect: the vector to initialize
 *
 */
__global__ void empty_vect_gpu(int m, double* vect)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if (tid < m)
        vect[tid] = 0;
} /* ----- end of empty_vect_gpu ----- */

/*!
 *  \brief Initializes to 0
 *
 *  Initializes to 0
 * 
 * \param m: the size of the vector
 * \param vect: the vector to initialize
 *
 */
__global__ void usc_to_u(int msc, double* u_sc, double* u, int ntheta_int, int delete_circles, int smoother)
{
    int i                = blockIdx.x * blockDim.x + threadIdx.x;
    int n_lines_radial_b = ceil((double)ntheta_int / 2);
    int n_lines_radial_w = floor((double)ntheta_int / 2);

    //computation of indices in the total vector u corresponding to the indices in u_sc
    if (i < msc) {
        // for (long unsigned int i = 0; i < msc; ++i) {
        int row;
        int col;

        if (smoother < 2) { //circle
            row = i / ntheta_int; //row within the smoother
            row = row * 2; //row within the whole matrix
            col = i % ntheta_int;
            if (smoother == 1) { //white
                row++;
            }
        }
        else { //radial
            if (smoother == 2) { //black
                row = i / n_lines_radial_b; //row within the smoother
                col = i % n_lines_radial_b; //col within the smoother
                col = col * 2; //col within the whole matrix
            }
            else { //white
                row = i / n_lines_radial_w; //row within the smoother
                col = i % n_lines_radial_w; //col within the smoother
                col = col * 2 + 1; //col within the whole matrix
            }
            row += delete_circles; //row within the whole matrix
        }

        int index = row * ntheta_int + col;

        u[index] = u_sc[i];
    }
} /* ----- end of usc_to_u ----- */

/*!
 *  \brief Computes the error
 *
 *  Computes the error compared to the theoretical solution and sclaed by its norm
 *
 * \return the error vector
 *
 */
std::vector<double> gmgpolar::compute_error()
{
    //compute the error vector (on level 0) between the computed and the analytical solution
    std::vector<double> error(v_level[0]->m);
    int mod_pk       = gyro::icntl[Param::mod_pk];
    double Rmax      = gyro::dcntl[Param::R];
    double kappa_eps = gyro::dcntl[Param::kappa_eps];
    double delta_e   = gyro::dcntl[Param::delta_e];
    int nr           = v_level[0]->nr;
    int ntheta       = v_level[0]->ntheta_int;
    int m            = v_level[0]->m;

    compute_error_gpu<<<ceil((double)m / 128), 128>>>(v_level[0]->dev_error, m, nr, v_level[0]->dev_r,
                                                      v_level[0]->dev_theta, ntheta, v_level[0]->dev_u, Rmax, mod_pk,
                                                      kappa_eps, delta_e);

    cudaMemcpy(&error[0], v_level[0]->dev_error, m * sizeof(double), cudaMemcpyDeviceToHost);

    // TODO: return 2-norm of error

    return error;
}
/* ----- end of gmgpolar::compute_error ----- */

/*!
 *  \brief Multigrid iterations on level l
 *
 *  Multigrid iterations on level l
 * 
 * \param l: the level (0=finest)
 *
 */
void gmgpolar::multigrid_cycle_extrapol(int l)
{
    //l = current level we are on (from class level)
    if (gyro::icntl[Param::verbose] > 2) {
        std::cout << "********************************************" << std::endl;
        std::cout << "MULTIGRID ON LEVEL " << l << std::endl;
        std::cout << "nr = " << v_level[l]->nr << ", ntheta = " << v_level[l]->ntheta << std::endl;
        std::cout << "********************************************" << std::endl;
    }

    //time measuring
    double t, t_total_tmp, t_smoothing_tmp;
    TIC;
    t_total_tmp     = t;
    t_smoothing_tmp = t;

    //! v1 presmoothing steps (on level l)
    for (int v = 0; v < gyro::icntl[Param::v1]; ++v) {
        //std::cout << "smoothing ... \n";

        // iterate over the 4 smoothers (0:c/b, 1:c/w, 2:r/b, 3:r/w)
        for (int smoother = 0; smoother < 4; smoother++) {
            TIC;

            //call the smooothing function, the result is u_sc which is directly inserted into u
            v_level[l]->multigrid_smoothing(smoother);

            if (gyro::icntl[Param::verbose] > 2)
                std::cout << "SMOOTHER: " << smoother << " = " << TOC << "\n";
        }
    }

    t = t_smoothing_tmp;
    t_smoothing += TOC;
    TIC;

    //! compute residual (of level l)
    //even if we have extrapolation, compute just normal residual (extrapolation-restriction follows in the next step)
    gmgpolar::compute_residual(l);

    t_residual += TOC;
    TIC;

    v_level[l + 1]->fVec = std::vector<double>(v_level[l + 1]->m, 0);

    //! Restriction of residual (coarsening)
    //the restricted residual of level l becomes the right hand side on level l+1, f(l+1)=res_coarse(l)
    //res(l+1)=P^T*res(l)      //trans=1 (Restriction, P^T)
    //P=(ncoarse*mc), ncoarse=m(size of fine level), mc(size of coarse level), coarse_r=coarse_nodes_list_r
    int nblocks = ceil((double)v_level[l]->m / NTHREADS_CUDA);
    apply_restriction_bi_gpu<<<nblocks, NTHREADS_CUDA>>>(
        v_level[l + 1]->dev_fVec, v_level[l]->dev_res, v_level[l]->m, v_level[l]->mc, v_level[l]->dev_thetaplus,
        v_level[l]->ntheta_int, v_level[l]->dev_hplus, v_level[l]->dev_coarse_nodes_list_type_shift,
        v_level[l]->dev_coarse_nodes_list_type_start, v_level[l]->coarse_nodes_list_type[0].size(),
        v_level[l]->coarse_nodes_list_type[1].size(), v_level[l]->coarse_nodes_list_type[2].size(),
        v_level[l]->coarse_nodes_list_type[3].size(), v_level[l]->coarse_nodes_list_type[4].size(),
        v_level[l]->coarse_nodes_list_type[5].size(), v_level[l]->coarse_nodes_list_type[6].size(),
        v_level[l]->coarse_nodes_list_type[7].size(), v_level[l]->dev_coarse_nodes_list_type0,
        v_level[l]->dev_coarse_nodes_list_type1, v_level[l]->dev_coarse_nodes_list_type2,
        v_level[l]->dev_coarse_nodes_list_type3, v_level[l]->dev_coarse_nodes_list_type4,
        v_level[l]->dev_coarse_nodes_list_type5, v_level[l]->dev_coarse_nodes_list_type6,
        v_level[l]->dev_coarse_nodes_list_type7, v_level[l]->dev_coarse_nodes_list_i0,
        v_level[l]->dev_coarse_nodes_list_i1, v_level[l]->dev_coarse_nodes_list_i2,
        v_level[l]->dev_coarse_nodes_list_i3, v_level[l]->dev_coarse_nodes_list_i4,
        v_level[l]->dev_coarse_nodes_list_i5, v_level[l]->dev_coarse_nodes_list_i6,
        v_level[l]->dev_coarse_nodes_list_i7, v_level[l]->dev_coarse_nodes_list_j0,
        v_level[l]->dev_coarse_nodes_list_j1, v_level[l]->dev_coarse_nodes_list_j2,
        v_level[l]->dev_coarse_nodes_list_j3, v_level[l]->dev_coarse_nodes_list_j4,
        v_level[l]->dev_coarse_nodes_list_j5, v_level[l]->dev_coarse_nodes_list_j6,
        v_level[l]->dev_coarse_nodes_list_j7);

    t_restriction += TOC;
    TIC;

    //! iterative call of multigrid_cycle_extrapol
    v_level[l + 1]->u.assign(v_level[l + 1]->m, 0); //zero u in every iteration
    nblocks = ceil((double)v_level[l + 1]->m / NTHREADS_CUDA);
    empty_vect_gpu<<<nblocks, NTHREADS_CUDA>>>(v_level[l + 1]->m, v_level[l + 1]->dev_u);
    if (l == levels - 2) {
        cudaMemcpy(&(v_level[l + 1]->fVec[0]), v_level[l + 1]->dev_fVec, sizeof(double) * v_level[l + 1]->fVec.size(),
                   cudaMemcpyDeviceToHost);

        //solve exactly on the coarsest level for the error (A * error = res) (use whole A from coarsest level)
        // check for the second coarsest level (levels-2), because we have no smoothing on the coarsest level (and thus no Asc and no Asc_ortho)
        //v_level[l + 1]->u = v_level[l + 1]->solve_gaussian_elimination(
        //    v_level[l + 1]->row_Ac_LU, v_level[l + 1]->col_Ac_LU, v_level[l + 1]->vals_Ac_LU, v_level[l + 1]->fVec);

        // Now a call on the GPU

        facto_gaussian_elimination_gpu(&(v_level[l + 1]->row_indices[0]), &(v_level[l + 1]->col_indices[0]),
                                       &(v_level[l + 1]->vals[0]), v_level[l + 1]->nz, v_level[l + 1]->m,
                                       &(v_level[l + 1]->fVec[0]), &(v_level[l + 1]->u[0]));

        // TODO: only for coarsest
        cudaMemcpy(v_level[l + 1]->dev_u, &(v_level[l + 1]->u[0]), sizeof(double) * v_level[l]->mc,
                   cudaMemcpyHostToDevice);
        t_Ac += TOC;
        TIC;
    }
    else {
        multigrid_cycle_extrapol(l + 1);
        TIC;
    }

    if (gyro::icntl[Param::verbose] > 2) {
        std::cout << "********************************************" << std::endl;
        std::cout << "BACK ON LEVEL " << l << std::endl;
        std::cout << "********************************************" << std::endl;
    }

    std::vector<double> error_fine = std::vector<double>(v_level[l]->m, 0);
    //! Prolongation of error_coarse to error_fine
    //error(l) = P * error(l+1)
    nblocks = ceil((double)v_level[l]->m / NTHREADS_CUDA);
    // gyro::disp(v_level[l]->u, "u");
    // empty_vect_gpu<<<nblocks, NTHREADS_CUDA>>>(v_level[l]->m, v_level[l]->dev_u);
    apply_prolongation_bi_gpu<<<nblocks, NTHREADS_CUDA>>>(
        v_level[l]->dev_u, v_level[l + 1]->dev_u, v_level[l]->m, v_level[l]->mc, v_level[l]->dev_thetaplus,
        v_level[l]->ntheta_int, v_level[l]->dev_hplus, v_level[l]->dev_coarse_nodes_list_type_shift,
        v_level[l]->dev_coarse_nodes_list_type_start, v_level[l]->coarse_nodes_list_i[0].size(),
        v_level[l]->coarse_nodes_list_i[1].size(), v_level[l]->coarse_nodes_list_i[2].size(),
        v_level[l]->coarse_nodes_list_i[3].size(), v_level[l]->coarse_nodes_list_i[4].size(),
        v_level[l]->coarse_nodes_list_i[5].size(), v_level[l]->coarse_nodes_list_i[6].size(),
        v_level[l]->coarse_nodes_list_i[7].size(), v_level[l]->dev_coarse_nodes_list_type0,
        v_level[l]->dev_coarse_nodes_list_type1, v_level[l]->dev_coarse_nodes_list_type2,
        v_level[l]->dev_coarse_nodes_list_type3, v_level[l]->dev_coarse_nodes_list_type4,
        v_level[l]->dev_coarse_nodes_list_type5, v_level[l]->dev_coarse_nodes_list_type6,
        v_level[l]->dev_coarse_nodes_list_type7, v_level[l]->dev_coarse_nodes_list_i0,
        v_level[l]->dev_coarse_nodes_list_i1, v_level[l]->dev_coarse_nodes_list_i2,
        v_level[l]->dev_coarse_nodes_list_i3, v_level[l]->dev_coarse_nodes_list_i4,
        v_level[l]->dev_coarse_nodes_list_i5, v_level[l]->dev_coarse_nodes_list_i6,
        v_level[l]->dev_coarse_nodes_list_i7, v_level[l]->dev_coarse_nodes_list_j0,
        v_level[l]->dev_coarse_nodes_list_j1, v_level[l]->dev_coarse_nodes_list_j2,
        v_level[l]->dev_coarse_nodes_list_j3, v_level[l]->dev_coarse_nodes_list_j4,
        v_level[l]->dev_coarse_nodes_list_j5, v_level[l]->dev_coarse_nodes_list_j6,
        v_level[l]->dev_coarse_nodes_list_j7);

    // std::cout << "solution corrected \n";
    t_prolongation += TOC;
    TIC;

    t_smoothing_tmp = t;

    //! v2 postsmoothing steps (on level l)
    for (int v = 0; v < gyro::icntl[Param::v2]; ++v) {
        //std::cout << "smoothing ... \n";

        // iterate over the 4 smoothers (0:c/b, 1:c/w, 2:r/b, 3:r/w)
        for (int smoother = 0; smoother < 4; smoother++) {
            TIC;

            //call the smooothing function, the result is u_sc which is directly inserted into u
            v_level[l]->multigrid_smoothing(smoother);

            if (gyro::icntl[Param::verbose] > 2)
                std::cout << "SMOOTHER: " << smoother << " = " << TOC << "\n";
        }
    }
    //std::cout << "post-smoothing done \n";
    t = t_smoothing_tmp;
    t_smoothing += TOC;
    TIC;

    t = t_total_tmp;
    t_total += TOC;
} /* ----- end of gmgpolar::multigrid_cycle_extrapol ----- */

/*! \brief Applies smoothing
 *
 * For all lines of the smoother s and the colour c.
 * The matrix A_sc is for all lines of points of the smoother s/c.
 * The result is in u.
 * 
 * if gyro::icntl[Param::smoother]=
 * - 3: block Gauss Seidel for all 4 smoothers.
 * 
 * \param smoother: the smoother
*/
void level::multigrid_smoothing(int smoother)
{
    double t, t_smoothing_tmp;
    TIC;
    t_smoothing_tmp = t;

    int n_threads = 32;
    int n_blocks  = ceil((double)m / (double)n_threads);

    //! solving for u_sc (solution vector)
    std::vector<double> u_sc(m_sc[smoother], 0);

    //create f_sc
    int msc = m_sc[smoother];

    build_fsc_gpu<<<n_blocks, n_threads>>>(dev_fsc[smoother], smoother, m, nr, ntheta_int, delete_circles, dev_fVec);

    t_f_sc += TOC;
    TIC;

    std::vector<double> f_total(m_sc[smoother], 0);

    apply_Asc_ortho_gpu<<<n_blocks, n_threads>>>(
        m, dev_fsc[smoother], dev_u, smoother, nr_int, dev_r, dev_hplus, dev_coarse_nodes_list_r, ntheta_int, dev_theta,
        dev_thetaplus, dev_coarse_nodes_list_theta, delete_circles, gyro::icntl[Param::DirBC_Interior],
        gyro::icntl[Param::mod_pk], gyro::dcntl[Param::r0_DB], gyro::icntl[Param::prob], gyro::dcntl[Param::kappa_eps],
        gyro::dcntl[Param::delta_e], gyro::dcntl[Param::R]);

    cudaMemcpy((void*)&f_total[0], (void*)dev_fsc[smoother], m_sc[smoother] * sizeof(double), cudaMemcpyDeviceToHost);

    t_Asc_ortho += TOC;
    TIC;
    // u_sc = solve_gaussian_elimination(A_Zebra_r_LU[smoother], A_Zebra_c_LU[smoother], A_Zebra_v_LU[smoother], f_total);
    // Call on the GPU
    facto_gaussian_elimination_gpu(&(A_Zebra_r[smoother][0]), &(A_Zebra_c[smoother][0]), &(A_Zebra_v[smoother][0]),
                                   nz_sc[smoother], m_sc[smoother], &f_total[0], &u_sc[0]);

    cudaMemcpy(dev_fsc[smoother], &u_sc[0], m_sc[smoother] * sizeof(double), cudaMemcpyHostToDevice);

    t_Asc += TOC;
    TIC;

    usc_to_u<<<n_blocks, n_threads>>>(msc, dev_fsc[smoother], dev_u, ntheta_int, delete_circles, smoother);
    // usc_to_u<<<1, 1>>>(msc, dev_fsc[smoother], dev_u, ntheta_int, delete_circles, smoother);

    t_indices += TOC;
    TIC;

    t = t_smoothing_tmp;
    t_smoothing += TOC;
    //std::cout << "smoothing end \n";
} /* ----- end of level::multigrid_smoothing ----- */

/*****************************************************************************************************************************
 * ******************   COMPUTE_RESIDUAL
 * ***************************************************************************************************************************/
/*!
 *  \brief Computes the residual
 *
 *  Computes the residual based on the approximation res = b-Au
 * 
 * \param l: the level (0=finest)
 *
 */
__global__ void compute_residual_gpu(int m, double* u, double* Au, double* fVec, double* res, int nr_int, double* r,
                                     double* hplus, int ntheta_int, double* theta, double* thetaplus,
                                     int DirBC_Interior, int mod_pk, int prob, double kappa_eps, double delta_e,
                                     double Rmax)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if (tid < m) {
        Au[tid] = 0;
    }

    apply_A_gpu(m, u, Au, nr_int, r, hplus, ntheta_int, theta, thetaplus, DirBC_Interior, mod_pk, prob, kappa_eps,
                delta_e, Rmax);

    //compute res = f - A * u
    if (tid < m) {
        res[tid] = fVec[tid] - Au[tid];
    }

} /* ----- end of gmgpolar::compute_residual ----- */

/*!
 *  \brief Apply A on the GPU
 *
 *  Apply A on the GPU without residual
 * 
 * \param l: the level (0=finest)
 *
 */
__global__ void apply_A_SA(int m, double* u, double* Au, int nr_int, double* r, double* hplus, int ntheta_int,
                           double* theta, double* thetaplus, int DirBC_Interior, int mod_pk, int prob, double kappa_eps,
                           double delta_e, double Rmax)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if (tid < m)
        Au[tid] = 0;

    apply_A_gpu(m, u, Au, nr_int, r, hplus, ntheta_int, theta, thetaplus, DirBC_Interior, mod_pk, prob, kappa_eps,
                delta_e, Rmax);
} /* ----- end of gmgpolar::compute_residual ----- */

/*!
 *  \brief Applies the operator A without construction (deprecated)
 *
 * Applies the matrix A of the system based on 9p FD. Relations with Dirichlet 
 * boundary condition nodes are shifted to the RHS to keep a symmetric operator.
 *
 * everything expressed in polar coordinates
 * allows anisotropy in r
 * - r: array of discretized r values
 * - ntheta: number of discretization intervals in angle phi (also variable m)
 *
 * solves the transformed DGL:  r*u_rr + u_r + 1/r*u_phiphi = rf fuer r>0, r<=1,
 * fuer r=0, approximate: u_0 = 1/m sum_j u(h,jk)-(h/2)^2*f
 * sorting of nodes ringwise. first node origin, last ntheta
 * Node all the way out
 *
 */
__device__ void apply_A_gpu(int m, double* u, double* Au, int nr_int, double* r, double* hplus, int ntheta_int,
                            double* theta, double* thetaplus, int DirBC_Interior, int mod_pk, int prob,
                            double kappa_eps, double delta_e, double Rmax)
{
    int tid, row, col, i, ip, im, j;
    double v;
    double kt, ktmin1, hs, hsmin1;
    double arr, att, arr_PI, arr_top, art_top, arr_bottom, art_bottom, att_right, art_right, att_left, art_left;

    tid = blockIdx.x * blockDim.x + threadIdx.x;
    row = tid;
    i   = tid % ntheta_int;
    j   = tid / ntheta_int;
    // Dirichlet
    if (j == 0 && DirBC_Interior) {
        Au[tid] += u[tid];
        return;
    }
    else if (j >= nr_int && tid < m) {
        Au[tid] += u[tid];
        return;
    }
    else if (tid < m) {
        ip = (i + 1 > ntheta_int - 1) ? 0 : i + 1;
        im = (i - 1 < 0) ? ntheta_int - 1 : i - 1;

        arr     = arr_gpu(r[j], theta[i], 0, mod_pk, prob, kappa_eps, delta_e, Rmax);
        att     = att_gpu(r[j], theta[i], 0, mod_pk, prob, kappa_eps, delta_e, Rmax);
        arr_top = arr_gpu(r[j + 1], theta[i], 0, mod_pk, prob, kappa_eps, delta_e, Rmax);
        if (j > 0)
            arr_bottom = arr_gpu(r[j - 1], theta[i], 0, mod_pk, prob, kappa_eps, delta_e, Rmax);
        att_right = att_gpu(r[j], theta[ip], 0, mod_pk, prob, kappa_eps, delta_e, Rmax);
        att_left  = att_gpu(r[j], theta[im], 0, mod_pk, prob, kappa_eps, delta_e, Rmax);
        if (j == 0 && !DirBC_Interior)
            arr_PI = arr_gpu(r[j], theta[i] + PI, 0, mod_pk, prob, kappa_eps, delta_e, Rmax);
        if (mod_pk > 0) {
            if (j < nr_int - 1)
                art_top = art_gpu(r[j + 1], theta[i], 0, mod_pk, prob, kappa_eps, delta_e, Rmax);
            if (j > 1 || j > 0 && !DirBC_Interior)
                art_bottom = art_gpu(r[j - 1], theta[i], 0, mod_pk, prob, kappa_eps, delta_e, Rmax);
            art_right = art_gpu(r[j], theta[ip], 0, mod_pk, prob, kappa_eps, delta_e, Rmax);
            art_left  = art_gpu(r[j], theta[im], 0, mod_pk, prob, kappa_eps, delta_e, Rmax);
        }

        kt     = thetaplus[i];
        ktmin1 = thetaplus[im];

        hs = hplus[j];
        if (j > 0)
            hsmin1 = hplus[j - 1];
        else
            hsmin1 = 2 * r[0];

        // 9-Point Stencil (attention order; first bottom, then right, then top right, top left, bottom left, middle)
        // Bottom
        if (j > 1 || !DirBC_Interior) {
            // second circle (at least one 'real' circle with ntheta_int nodes which is closer to the origin)
            if (j > 0) {
                // a l interieur sans contact aux conditions aux limites
                col = (j - 1) * ntheta_int + i;
                v   = -0.5 * (kt + ktmin1) / hsmin1 * (arr_bottom + arr);
            }
            else { // across the origin
                if (i + 1 > ntheta_int / 2)
                    col = row - ntheta_int / 2; // half a turn back
                else
                    col = row + ntheta_int / 2; // half a turn further
                v = -0.5 * (kt + ktmin1) / hsmin1 * (arr_PI + arr);
            }
            Au[row] += v * u[col];
        }
        // Bottom-left
        if (mod_pk > 0 && (j > 1 || j > 0 && !DirBC_Interior)) {
            col = (j - 1) * ntheta_int + im;
            v   = -0.25 * (art_left + art_bottom);

            Au[row] += v * u[col];
        }

        // Left
        col = j * ntheta_int + im;
        v   = -0.5 * (hs + hsmin1) / ktmin1 * (att_left + att);
        Au[row] += v * u[col];

        // Top-Left
        if (mod_pk > 0 && j < nr_int - 1) {
            col = (j + 1) * ntheta_int + im;
            v   = 0.25 * (art_left + art_top);
            Au[row] += v * u[col];
        }

        // Top
        if (j < nr_int - 1) {
            col = (j + 1) * ntheta_int + i;
            v   = -0.5 * (kt + ktmin1) / hs * (arr + arr_top);
            Au[row] += v * u[col];
        }

        // Top-Right
        if (mod_pk > 0 && j < nr_int - 1) {
            col = (j + 1) * ntheta_int + ip;
            v   = -0.25 * (art_top + art_right);
            Au[row] += v * u[col];
        }

        // Right
        col = j * ntheta_int + ip;
        v   = -0.5 * (hs + hsmin1) / kt * (att + att_right);
        Au[row] += v * u[col];

        // Top-Right
        if (mod_pk > 0 && (j > 1 || j > 0 && !DirBC_Interior)) {
            col = (j - 1) * ntheta_int + ip;
            v   = 0.25 * (art_bottom + art_right);
            Au[row] += v * u[col];
        }

        // Middle
        col = row;
        if (j > 0)
            // alt: neumBound
            v = 0.5 * (kt + ktmin1) * (arr + arr_top) / hs + 0.5 * (hs + hsmin1) * (att + att_right) / kt +
                0.5 * (kt + ktmin1) * (arr_bottom + arr) / hsmin1 + 0.5 * (hs + hsmin1) * (att_left + att) / ktmin1;
        else { //across the origin; j=0
            v = 0.5 * (kt + ktmin1) * (arr + arr_top) / hs + 0.5 * (hs + hsmin1) * (att + att_right) / kt +
                0.5 * (hs + hsmin1) * (att_left + att) / ktmin1 + 0.5 * (kt + ktmin1) * (arr_PI + arr) / hsmin1;
        }
        Au[row] += v * u[col];
    }
} /* ----- end of level::apply_A0 ----- */

__device__ double coeff_gpu(double r, double theta, int verbose, int mod_pk, int prob, double kappa_eps, double delta_e,
                            double Rmax)
{
    double c1, c2, c3, c4;
    double coeff;
    if (prob == 3 || prob == 5) {
        c1    = 2.0 / (2.6 + 3.14);
        c2    = 1.3;
        c3    = 1.0 / 0.09;
        c4    = 1.3 / Rmax;
        coeff = c1 * (c2 + atan(c3 * (1 - c4 * r)));
    }
    else {
        c1    = 0;
        c2    = 0;
        c3    = 0;
        c4    = 0;
        coeff = 1.0;
    }
    return coeff;
} /* ----- end of level::coeff ----- */

__device__ double detDFinv_gpu(double r, double theta, int verbose, int mod_pk, int prob, double kappa_eps,
                               double delta_e, double Rmax)
{
    double detDFinv_r = -1;
    if (mod_pk == 1)
        detDFinv_r = r * (1 + kappa_eps) * (1 - kappa_eps - 2 * delta_e * r * cos(theta));
    else if (mod_pk == 0)
        detDFinv_r = r;
    return detDFinv_r;
} /* ----- end of level::detDFinv ----- */

__device__ double arr_gpu(double r, double theta, int verbose, int mod_pk, int prob, double kappa_eps, double delta_e,
                          double Rmax)
{
    double arr_r      = 1.0;
    double detDFinv_r = detDFinv_gpu(r, theta, verbose, mod_pk, prob, kappa_eps, delta_e, Rmax);
    double coeff_r    = coeff_gpu(r, theta, verbose, mod_pk, prob, kappa_eps, delta_e, Rmax);
    if (mod_pk == 1)
        arr_r = 0.5 * r * r * ((pow((1 + kappa_eps) * cos(theta), 2) + pow((1 - kappa_eps) * sin(theta), 2))) *
                coeff_r / abs(detDFinv_r);
    else if (mod_pk == 0)
        arr_r = 0.5 * detDFinv_r * coeff_r;
    return arr_r;
} /* ----- end of level::arr ----- */

__device__ double art_gpu(double r, double theta, int verbose, int mod_pk, int prob, double kappa_eps, double delta_e,
                          double Rmax)
{
    double art_r      = 1.0;
    double detDFinv_r = detDFinv_gpu(r, theta, verbose, mod_pk, prob, kappa_eps, delta_e, Rmax);
    double coeff_r    = coeff_gpu(r, theta, verbose, mod_pk, prob, kappa_eps, delta_e, Rmax);
    if (mod_pk == 1)
        art_r = (-4 * r * kappa_eps * cos(theta) - 2 * delta_e * pow(r, 2) * (1 - kappa_eps)) * coeff_r * sin(theta) /
                abs(detDFinv_r);
    else if (mod_pk == 0)
        art_r = 0;
    return art_r;
} /* ----- end of level::art ----- */

__device__ double att_gpu(double r, double theta, int verbose, int mod_pk, int prob, double kappa_eps, double delta_e,
                          double Rmax)
{
    double att_r      = 1.0;
    double detDFinv_r = detDFinv_gpu(r, theta, verbose, mod_pk, prob, kappa_eps, delta_e, Rmax);
    double coeff_r    = coeff_gpu(r, theta, verbose, mod_pk, prob, kappa_eps, delta_e, Rmax);
    if (mod_pk == 1) {
        att_r = 0.5 * coeff_r *
                (pow((1 + kappa_eps) * sin(theta), 2) + pow((1 - kappa_eps) * cos(theta) - 2 * delta_e * r, 2)) /
                abs(detDFinv_r);
    }
    else if (mod_pk == 0)
        att_r = 0.5 * coeff_r / r;
    return att_r;
} /* ----- end of level::att ----- */

/*****************************************************************************************************************************
 * ******************   COMPUTE_ERROR
 * ***************************************************************************************************************************/
/*!
 *  \brief Computes the residual
 *
 *  Computes the residual based on the approximation res = b-Au
 * 
 * \param l: the level (0=finest)
 *
 */
__global__ void compute_error_gpu(double* error, int m, int nr, double* r, double* theta, int ntheta_int, double* u,
                                  double Rmax, int mod_pk, double kappa_eps, double delta_e)
{
    double computed_solution;
    double theoretical_sol;
    int index, tid, tidr, tidtheta;
    tid      = blockIdx.x * blockDim.x + threadIdx.x;
    tidtheta = tid % ntheta_int;
    tidr     = floor((double)tid / ((double)ntheta_int));

    if (tidr < nr) {
        if (tidtheta < ntheta_int) {
            def_solution_rt_gpu(&theoretical_sol, m, r[tidr], theta[tidtheta], Rmax, mod_pk, kappa_eps, delta_e);
            index             = tidr * ntheta_int + tidtheta; //corresponding index of the solution vector
            computed_solution = u[index];
            error[index]      = abs(theoretical_sol - computed_solution);
        }
    }

    // }

} /* ----- end of gmgpolar::compute_error_gpu ----- */

/*!
 *  \brief Evaluates the solution
 *
 *  Evaluates the solution from the r coordinate on all theta positions
 * 
 * \param r_i: the r coordinate of the node
 * \param theta: vector theta (0, ntheta_int)
 * \param ntheta: number of values in theta
 * \param verbose: verbose level for debug
 *
 * \return the solution vector
 */
__device__ void def_solution_rt_gpu(double* sol, int m, double r_j, double theta_i, double Rmax, int mod_pk,
                                    double kappa_eps, double delta_e)
{
    double x, y;
    if ((Rmax - 1.3) > 1e-7)
        printf("Solution only valid for Rmax=1.3.\n");

    //   throw std::runtime_error("Solution only valid for Rmax=1.3.");
    trafo_gpu(r_j, theta_i, &x, &y, mod_pk, kappa_eps, delta_e);
    *sol = (Rmax * Rmax - r_j * r_j) * cos(2 * PI * x) * sin(2 * PI * y);
    return;
} /* ----- end of eval_def_solution_vec ----- */

/*!
 *  \brief Transform from polar to cartesian
 *
 *  Transform one couple r_i,theta_j from polar to cartesian
 * 
 * \param r_i: the r coordinate of the node
 * \param theta: vector theta (0, ntheta_int)
 * \param ntheta: number of values in theta
 * \param x: vector x
 * \param y: vector y
 * \param verbose: verbose level for debug
 *
 */
__device__ void trafo_gpu(double r_j, double theta_i, double* x, double* y, int mod_pk, double kappa_eps,
                          double delta_e)
{
    if (mod_pk == 1) {
        *x = (1 - kappa_eps) * r_j * cos(theta_i) - delta_e * pow(r_j, 2);
        *y = (1 + kappa_eps) * r_j * sin(theta_i);
    }
    else if (mod_pk == 2) {
        *x = 1 / kappa_eps * (1 - sqrt(1 + kappa_eps * (kappa_eps + 2 * r_j * cos(theta_i))));
        *y = (delta_e / sqrt(1 - pow(kappa_eps, 2) / 4) * r_j * cos(theta_i)) /
             (2 - sqrt(1 + kappa_eps * (kappa_eps + 2 * r_j * cos(theta_i))));
    }
    else {
        *x = r_j * cos(theta_i);
        *y = r_j * sin(theta_i);
    }
} /* ----- end of trafo ----- */

/*****************************************************************************************************************************
 * ******************   MULTIGIRD_CYCLE_EXTRAPOL
 * ***************************************************************************************************************************/
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
__global__ void apply_restriction_bi_gpu(
    double* Pu, double* u, int m, int mc, double* thetaplus, int ntheta_int, double* hplus,
    int* coarse_nodes_list_type_shift, int* coarse_nodes_list_type_start, int k0, int k1, int k2, int k3, int k4,
    int k5, int k6, int k7, int* coarse_nodes_list_type0, int* coarse_nodes_list_type1, int* coarse_nodes_list_type2,
    int* coarse_nodes_list_type3, int* coarse_nodes_list_type4, int* coarse_nodes_list_type5,
    int* coarse_nodes_list_type6, int* coarse_nodes_list_type7, int* coarse_nodes_list_i0, int* coarse_nodes_list_i1,
    int* coarse_nodes_list_i2, int* coarse_nodes_list_i3, int* coarse_nodes_list_i4, int* coarse_nodes_list_i5,
    int* coarse_nodes_list_i6, int* coarse_nodes_list_i7, int* coarse_nodes_list_j0, int* coarse_nodes_list_j1,
    int* coarse_nodes_list_j2, int* coarse_nodes_list_j3, int* coarse_nodes_list_j4, int* coarse_nodes_list_j5,
    int* coarse_nodes_list_j6, int* coarse_nodes_list_j7)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid < mc) {
        Pu[tid] = 0;
        int row, col, i, j, count, shift, start;
        double k_prev, k_next, h_prev, h_next, denum, val;
        // Coarse nodes
        shift = coarse_nodes_list_type_shift[0];
        start = coarse_nodes_list_type_start[0];
        count = 0;
        for (int k = start; k < k0; k += shift) {
            row = coarse_nodes_list_type0[k];
            i   = coarse_nodes_list_i0[count];
            j   = coarse_nodes_list_j0[count];
            count++;
            col = k;
            if (col == tid) {
                val = 1.0;
                Pu[tid] += val * u[row];
            }
        }
        // Coarse nodes in same r (2i, 2j+1)
        shift = coarse_nodes_list_type_shift[1];
        start = coarse_nodes_list_type_start[1];
        count = 0;
        for (int k = start; k < k1; k += shift) {
            row = coarse_nodes_list_type1[k];
            i   = coarse_nodes_list_i1[count];
            j   = coarse_nodes_list_j1[count];
            count++;
            k_prev = thetaplus[(i + ntheta_int - 1) % ntheta_int];
            k_next = thetaplus[i];
            denum  = 1 / (k_prev + k_next);

            // Previous coarse node (left)
            col = coarse_nodes_list_type1[k - 1];
            if (col == tid) {
                val = k_next * denum; // 1/2
                Pu[tid] += val * u[row];
            }

            // Next coarse node (right)
            col = coarse_nodes_list_type1[k + 1];
            if (col == tid) {
                // val = k_prev * denum; // 1/2
                val = k_prev * denum; // 1/2
                Pu[tid] += val * u[row];
            }
        }
        // Coarse nodes in same theta (2i+1, 2j)
        // - j = 1
        shift = coarse_nodes_list_type_shift[3];
        start = coarse_nodes_list_type_start[3];
        count = 0;
        for (int k = start; k < k3; k += shift) {
            row = coarse_nodes_list_type3[k];
            i   = coarse_nodes_list_i3[count];
            j   = coarse_nodes_list_j3[count];
            count++;
            h_prev = hplus[j - 1];
            h_next = hplus[j];
            denum  = 1 / (h_prev + h_next);

            // Next coarse node (top)
            col = coarse_nodes_list_type3[k + 1];
            if (col == tid) {
                // val = h_prev * denum; // 1/2
                val = h_prev * denum; // 1/2
                Pu[tid] += val * u[row];
            }
        }
        // - j = nr_int - 1
        shift = coarse_nodes_list_type_shift[4];
        start = coarse_nodes_list_type_start[4];
        count = 0;
        for (int k = start; k < k4; k += shift) {
            row = coarse_nodes_list_type4[k];
            i   = coarse_nodes_list_i4[count];
            j   = coarse_nodes_list_j4[count];
            count++;
            h_prev = hplus[j - 1];
            h_next = hplus[j];
            denum  = 1 / (h_prev + h_next);

            // Previous coarse node (bottom)
            col = coarse_nodes_list_type4[k - 1];
            if (col == tid) {
                // val = h_next * denum; // 1/2
                val = h_next * denum; // 1/2
                Pu[tid] += val * u[row];
            }
        }
        // - interior
        shift = coarse_nodes_list_type_shift[2];
        start = coarse_nodes_list_type_start[2];
        count = 0;
        for (int k = start; k < k2; k += shift) {
            row = coarse_nodes_list_type2[k];
            i   = coarse_nodes_list_i2[count];
            j   = coarse_nodes_list_j2[count];
            count++;
            h_prev = hplus[j - 1];
            h_next = hplus[j];
            denum  = 1 / (h_prev + h_next);

            // Previous coarse node (bottom)
            col = coarse_nodes_list_type2[k - 1];
            if (col == tid) {
                // val = h_next * denum; // 1/2
                val = h_next * denum; // 1/2
                Pu[tid] += val * u[row];
            }

            // Next coarse node (top)
            col = coarse_nodes_list_type2[k + 1];
            if (col == tid) {
                // val = h_prev * denum; // 1/2
                val = h_prev * denum; // 1/2
                Pu[tid] += val * u[row];
            }
        }
        // Coarse nodes in diagonals (2i+1, 2j+1)
        // - j = 1
        shift = coarse_nodes_list_type_shift[6];
        start = coarse_nodes_list_type_start[6];
        count = 0;
        for (int k = start; k < k6; k += shift) {
            row = coarse_nodes_list_type6[k];
            i   = coarse_nodes_list_i6[count];
            j   = coarse_nodes_list_j6[count];
            count++;
            k_prev = thetaplus[(i + ntheta_int - 1) % ntheta_int];
            k_next = thetaplus[i];
            h_prev = hplus[j - 1];
            h_next = hplus[j];
            denum  = 1 / ((k_prev + k_next) * (h_prev + h_next));

            // bottom_left
            col = coarse_nodes_list_type6[k + 2];
            if (col == tid) {
                val = h_prev * k_prev * denum; // isotrop: 1/4
                Pu[tid] += val * u[row];
            }
            // bottom_right
            col = coarse_nodes_list_type6[k + 1];
            if (col == tid) {
                val = h_prev * k_next * denum; // isotrop: 1/4
                Pu[tid] += val * u[row];
            }
        }
        // - j == nr_int - 1
        shift = coarse_nodes_list_type_shift[7];
        start = coarse_nodes_list_type_start[7];
        count = 0;
        for (int k = start; k < k7; k += shift) {
            row = coarse_nodes_list_type7[k];
            i   = coarse_nodes_list_i7[count];
            j   = coarse_nodes_list_j7[count];
            count++;
            k_prev = thetaplus[(i + ntheta_int - 1) % ntheta_int];
            k_next = thetaplus[i];
            h_prev = hplus[j - 1];
            h_next = hplus[j];
            denum  = 1 / ((k_prev + k_next) * (h_prev + h_next));

            // top_left
            col = coarse_nodes_list_type7[k - 1];
            if (col == tid) {
                val = h_next * k_prev * denum; // isotrop: 1/4
                Pu[tid] += val * u[row];
            }
            // top_right
            col = coarse_nodes_list_type7[k - 2];
            if (col == tid) {
                val = h_next * k_next * denum; // isotrop: 1/4
                Pu[tid] += val * u[row];
            }
        }
        // - interior
        shift = coarse_nodes_list_type_shift[5];
        start = coarse_nodes_list_type_start[5];
        count = 0;
        for (int k = start; k < k5; k += shift) {
            row = coarse_nodes_list_type5[k];
            i   = coarse_nodes_list_i5[count];
            j   = coarse_nodes_list_j5[count];
            count++;
            k_prev = thetaplus[(i + ntheta_int - 1) % ntheta_int];
            k_next = thetaplus[i];
            h_prev = hplus[j - 1];
            h_next = hplus[j];
            denum  = 1 / ((k_prev + k_next) * (h_prev + h_next));

            // bottom_left
            col = coarse_nodes_list_type5[k + 2];
            if (col == tid) {
                val = h_prev * k_prev * denum; // isotrop: 1/4
                Pu[tid] += val * u[row];
            }
            // bottom_right
            col = coarse_nodes_list_type5[k + 1];
            if (col == tid) {
                val = h_prev * k_next * denum; // isotrop: 1/4
                Pu[tid] += val * u[row];
            }

            // top_left
            col = coarse_nodes_list_type5[k - 1];
            if (col == tid) {
                val = h_next * k_prev * denum; // isotrop: 1/4
                Pu[tid] += val * u[row];
            }
            // top_right
            col = coarse_nodes_list_type5[k - 2];
            if (col == tid) {
                val = h_next * k_next * denum; // isotrop: 1/4
                Pu[col] += val * u[row];
            }
        }
    }
} /* ----- end of level::apply_restriction_bi ----- */

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
__global__ void apply_prolongation_bi_gpu(
    double* Pu, double* u, int m, int mc, double* thetaplus, int ntheta_int, double* hplus,
    int* coarse_nodes_list_type_shift, int* coarse_nodes_list_type_start, int k0, int k1, int k2, int k3, int k4,
    int k5, int k6, int k7, int* coarse_nodes_list_type0, int* coarse_nodes_list_type1, int* coarse_nodes_list_type2,
    int* coarse_nodes_list_type3, int* coarse_nodes_list_type4, int* coarse_nodes_list_type5,
    int* coarse_nodes_list_type6, int* coarse_nodes_list_type7, int* coarse_nodes_list_i0, int* coarse_nodes_list_i1,
    int* coarse_nodes_list_i2, int* coarse_nodes_list_i3, int* coarse_nodes_list_i4, int* coarse_nodes_list_i5,
    int* coarse_nodes_list_i6, int* coarse_nodes_list_i7, int* coarse_nodes_list_j0, int* coarse_nodes_list_j1,
    int* coarse_nodes_list_j2, int* coarse_nodes_list_j3, int* coarse_nodes_list_j4, int* coarse_nodes_list_j5,
    int* coarse_nodes_list_j6, int* coarse_nodes_list_j7)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid < m) {
        int k, row, col, i, j, count, shift, start;
        double k_prev, k_next, h_prev, h_next, denum, val;
        // Coarse nodes
        shift = coarse_nodes_list_type_shift[0];
        if (tid < k0) {
            count = tid;
            start = coarse_nodes_list_type_start[0];
            k     = tid * shift + start;
            row   = coarse_nodes_list_type0[k];
            i     = coarse_nodes_list_i0[count];
            j     = coarse_nodes_list_j0[count];
            col   = k;
            val   = 1.0;
            // printf("row: %d, k: %d, Pu: %f, 0\n", row, k, Pu[row]);
            Pu[row] += val * u[col];

            return;
        }
        tid -= k0;
        shift = coarse_nodes_list_type_shift[1];
        if (tid < k1) {
            // Coarse nodes in same r (2i, 2j+1)
            count  = tid;
            start  = coarse_nodes_list_type_start[1];
            k      = tid * shift + start;
            row    = coarse_nodes_list_type1[k];
            i      = coarse_nodes_list_i1[count];
            j      = coarse_nodes_list_j1[count];
            k_prev = thetaplus[(i + ntheta_int - 1) % ntheta_int];
            k_next = thetaplus[i];
            denum  = 1 / (k_prev + k_next);

            // Previous coarse node (left)
            col = coarse_nodes_list_type1[k - 1];
            val = k_next * denum; // 1/2
            // printf("row: %d, k: %d, Pu: %f, 0\n", row, k, Pu[row]);
            Pu[row] += val * u[col];

            // Next coarse node (right)
            col = coarse_nodes_list_type1[k + 1];
            val = k_prev * denum; // 1/2
            Pu[row] += val * u[col];

            return;
        }
        tid -= k1;
        shift = coarse_nodes_list_type_shift[3];
        if (tid < k3) {
            // Coarse nodes in same theta (2i+1, 2j)
            // - j = 1
            count  = tid;
            start  = coarse_nodes_list_type_start[3];
            k      = tid * shift + start;
            row    = coarse_nodes_list_type3[k];
            i      = coarse_nodes_list_i3[count];
            j      = coarse_nodes_list_j3[count];
            h_prev = hplus[j - 1];
            h_next = hplus[j];
            denum  = 1 / (h_prev + h_next);

            // Next coarse node (top)
            col = coarse_nodes_list_type3[k + 1];
            val = h_prev * denum; // 1/2
            // printf("row: %d, k: %d, Pu: %f, 0\n", row, k, Pu[row]);
            Pu[row] += val * u[col];

            return;
        }
        tid -= k3;
        shift = coarse_nodes_list_type_shift[4];
        if (tid < k4) {
            // - j = nr_int - 1
            count  = tid;
            start  = coarse_nodes_list_type_start[4];
            k      = tid * shift + start;
            row    = coarse_nodes_list_type4[k];
            i      = coarse_nodes_list_i4[count];
            j      = coarse_nodes_list_j4[count];
            h_prev = hplus[j - 1];
            h_next = hplus[j];
            denum  = 1 / (h_prev + h_next);

            // Previous coarse node (bottom)
            col = coarse_nodes_list_type4[k - 1];
            val = h_next * denum; // 1/2
            // printf("row: %d, k: %d, Pu: %f, 0\n", row, k, Pu[row]);
            Pu[row] += val * u[col];

            return;
        }
        tid -= k4;
        shift = coarse_nodes_list_type_shift[2];
        if (tid < k2) {
            // - interior
            count  = tid;
            start  = coarse_nodes_list_type_start[2];
            k      = tid * shift + start;
            row    = coarse_nodes_list_type2[k];
            i      = coarse_nodes_list_i2[count];
            j      = coarse_nodes_list_j2[count];
            h_prev = hplus[j - 1];
            h_next = hplus[j];
            denum  = 1 / (h_prev + h_next);

            // Previous coarse node (bottom)
            col = coarse_nodes_list_type2[k - 1];
            val = h_next * denum; // 1/2
            // printf("row: %d, k: %d, Pu: %f, 0\n", row, k, Pu[row]);
            Pu[row] += val * u[col];

            // Next coarse node (top)
            col = coarse_nodes_list_type2[k + 1];
            val = h_prev * denum; // 1/2
            Pu[row] += val * u[col];

            return;
        }
        tid -= k2;
        shift = coarse_nodes_list_type_shift[6];
        if (tid < k6) {
            // Coarse nodes in diagonals (2i+1, 2j+1)
            // - j = 1
            count  = tid;
            start  = coarse_nodes_list_type_start[6];
            k      = tid * shift + start;
            row    = coarse_nodes_list_type6[k];
            i      = coarse_nodes_list_i6[count];
            j      = coarse_nodes_list_j6[count];
            k_prev = thetaplus[(i + ntheta_int - 1) % ntheta_int];
            k_next = thetaplus[i];
            h_prev = hplus[j - 1];
            h_next = hplus[j];
            denum  = 1 / ((k_prev + k_next) * (h_prev + h_next));

            // bottom_left
            col = coarse_nodes_list_type6[k + 2];
            val = h_prev * k_prev * denum; // isotrop: 1/4
            // printf("row: %d, k: %d, Pu: %f, 0\n", row, k, Pu[row]);
            Pu[row] += val * u[col];
            // bottom_right
            col = coarse_nodes_list_type6[k + 1];
            val = h_prev * k_next * denum; // isotrop: 1/4
            Pu[row] += val * u[col];

            return;
        }
        tid -= k6;
        shift = coarse_nodes_list_type_shift[7];
        if (tid < k7) {
            // - j == nr_int - 1
            count  = tid;
            start  = coarse_nodes_list_type_start[7];
            k      = tid * shift + start;
            row    = coarse_nodes_list_type7[k];
            i      = coarse_nodes_list_i7[count];
            j      = coarse_nodes_list_j7[count];
            k_prev = thetaplus[(i + ntheta_int - 1) % ntheta_int];
            k_next = thetaplus[i];
            h_prev = hplus[j - 1];
            h_next = hplus[j];
            denum  = 1 / ((k_prev + k_next) * (h_prev + h_next));

            // top_left
            col = coarse_nodes_list_type7[k - 1];
            val = h_next * k_prev * denum; // isotrop: 1/4
            // printf("row: %d, k: %d, Pu: %f, 0\n", row, k, Pu[row]);
            Pu[row] += val * u[col];
            // top_right
            col = coarse_nodes_list_type7[k - 2];
            val = h_next * k_next * denum; // isotrop: 1/4
            Pu[row] += val * u[col];

            return;
        }
        tid -= k7;
        shift = coarse_nodes_list_type_shift[5];
        if (tid < k5) {
            // - interior
            count  = tid;
            start  = coarse_nodes_list_type_start[5];
            k      = tid * shift + start;
            row    = coarse_nodes_list_type5[k];
            i      = coarse_nodes_list_i5[count];
            j      = coarse_nodes_list_j5[count];
            k_prev = thetaplus[(i + ntheta_int - 1) % ntheta_int];
            k_next = thetaplus[i];
            h_prev = hplus[j - 1];
            h_next = hplus[j];
            denum  = 1 / ((k_prev + k_next) * (h_prev + h_next));

            // bottom_left
            col = coarse_nodes_list_type5[k + 2];
            val = h_prev * k_prev * denum; // isotrop: 1/4
            // printf("row: %d, k: %d, Pu: %f, 0\n", row, k, Pu[row]);
            Pu[row] += val * u[col];
            // bottom_right
            col = coarse_nodes_list_type5[k + 1];
            val = h_prev * k_next * denum; // isotrop: 1/4
            Pu[row] += val * u[col];

            // top_left
            col = coarse_nodes_list_type5[k - 1];
            val = h_next * k_prev * denum; // isotrop: 1/4
            Pu[row] += val * u[col];
            // top_right
            col = coarse_nodes_list_type5[k - 2];
            val = h_next * k_next * denum; // isotrop: 1/4
            Pu[row] += val * u[col];
        }
    }
} /* ----- end of level::apply_prolongation_bi ----- */

/*! \brief Create the RHS part corresponding to Asc_ortho for a smoother (on level l)
 *
 * Create f_sc, i.e. the RHS part corresponding to Asc_ortho for a smoother (on level l)
 * 
 * \param f_sc: the RHS part (out)
 * \param smoother: the smoother
*/
__global__ void build_fsc_gpu(double* f_sc, int smoother, int m, int nr, int ntheta_int, int delete_circles,
                              double* fVec)
{
    int n_indices_circle = delete_circles * ntheta_int; //delete_circles = index of radius of border between smoothers
    int index            = 0;
    int tid              = blockIdx.x * blockDim.x + threadIdx.x;

    if (smoother == 0) { //circle black smoother
        if (tid < n_indices_circle) {
            int r_index = tid / ntheta_int; //index in r-direction
            index       = tid / 2;
            //check if even r_index

            if (!(r_index % 2)) {
                index       = floorf(r_index / 2) * ntheta_int + tid % ntheta_int;
                f_sc[index] = fVec[tid]; //insert elements of f corresponding to the colour
            }
        }
    }
    else if (smoother == 1) { //circle white smoother
        //    for (int ind = 0; ind < n_indices_circle; ++ind) { //iteration over all elements in the smoother
        if (tid < n_indices_circle) {
            int ind_local = tid;
            int r_index   = ind_local / ntheta_int; //index in r-direction
            //check odd r_index

            if (r_index % 2) {
                index       = floorf(r_index / 2) * ntheta_int + tid % ntheta_int;
                f_sc[index] = fVec[tid]; //insert elements of f corresponding to the colour
            }
        }
    }

    else if (smoother == 2) { //radial black smoother
        //skip every second row
        //   for (int ind = n_indices_circle; ind < m; ++ind) { //iteration over all elements in the smoother
        if (tid >= n_indices_circle && tid < m) {
            int r_index     = tid / ntheta_int; //index in r-direction
            int theta_index = tid - r_index * ntheta_int; //index in theta-direction

            //check if black (even theta_index) or white (odd theta_index)
            if (!(theta_index % 2)) { //radial black, even
                index       = (tid - n_indices_circle) / 2;
                f_sc[index] = fVec[tid]; //insert elements of f corresponding to the colour
            }
        }
    }

    else { //radial white smoother
        //   for (int ind = n_indices_circle; ind < m; ++ind) { //iteration over all elements in the smoother
        if (tid >= n_indices_circle && tid < m) {
            int r_index     = tid / ntheta_int; //index in r-direction
            int theta_index = tid - r_index * ntheta_int; //index in theta-direction

            if (theta_index % 2) { //radial white, odd
                index       = (tid - n_indices_circle - 1) / 2;
                f_sc[index] = fVec[tid]; //insert elements of f corresponding to the colour
            }
        }
    }

} /* ----- end of level::build_subvectors ----- */

/*! \brief Applies the matrix A_sc_ortho for a smoother explicitely on the level l (deprecated)
 *
 * Asc_ortho corresponds to all lines of a smoother s with colour c and the columns not in Asc
 * 
 * \param smoother_todo: the smoother of this Asc_ortho matrix
*/
__global__ void apply_Asc_ortho_gpu(int m, double* Au, double* u, int smoother, int nr_int, double* r, double* hplus,
                                    int* coarse_nodes_list_r, int ntheta_int, double* theta, double* thetaplus,
                                    int* coarse_nodes_list_theta, int delete_circles, int DirBC_Interior, int mod_pk,
                                    double r0_DB, int prob, double kappa_eps, double delta_e, double Rmax)
{
    int tid, row, col, i, im, ip, j;
    double kt, ktmin1, hs, hsmin1, val;
    double arr, att, arr_bottom, att_left, att_right, arr_top, art_bottom, art_left, art_right, art_top;

    tid = blockIdx.x * blockDim.x + threadIdx.x;

    if (tid < m) {
        i = tid % ntheta_int;
        j = tid / ntheta_int;
        //! if the index and the smoother don't fit together, just skip this point!
        if (DirBC_Interior && j == 0 || j >= nr_int ||
            (j < delete_circles && ((j % 2 == 0 && smoother != 0) || (j % 2 != 0 && smoother != 1))) ||
            (j >= delete_circles && ((i % 2 == 0 && smoother != 2) || (i % 2 != 0 && smoother != 3)))) {
            return;
        }
        int row_index = get_ptr_sc(i, j, smoother, 1, delete_circles, nr_int, ntheta_int);
        row           = j * ntheta_int + i;

        ip = (i + 1 > ntheta_int - 1) ? 0 : i + 1;
        im = (i - 1 < 0) ? ntheta_int - 1 : i - 1;

        kt     = thetaplus[i];
        ktmin1 = thetaplus[im];

        hs = hplus[j];
        if (j > 0)
            hsmin1 = hplus[j - 1];
        else
            hsmin1 = 2 * r[0];

        if (smoother < 2 && j < delete_circles || smoother > 1 && j == delete_circles)
            arr = arr_gpu(r[j], theta[i], 0, mod_pk, prob, kappa_eps, delta_e, Rmax);
        if (smoother > 1 && j >= delete_circles)
            att = att_gpu(r[j], theta[i], 0, mod_pk, prob, kappa_eps, delta_e, Rmax);
        if ((smoother < 2 && j < delete_circles && j > 0 || smoother > 1 && j == delete_circles) &&
            (j > 1 || !DirBC_Interior))
            arr_bottom = arr_gpu(r[j - 1], theta[i], 0, mod_pk, prob, kappa_eps, delta_e, Rmax);
        if (smoother > 1 && j >= delete_circles)
            att_left = att_gpu(r[j], theta[im], 0, mod_pk, prob, kappa_eps, delta_e, Rmax);
        if (smoother > 1 && j >= delete_circles)
            att_right = att_gpu(r[j], theta[ip], 0, mod_pk, prob, kappa_eps, delta_e, Rmax);
        if (smoother < 2 && j < delete_circles && j < nr_int - 1)
            arr_top = arr_gpu(r[j + 1], theta[i], 0, mod_pk, prob, kappa_eps, delta_e, Rmax);
        if (mod_pk) {
            if (j > 1 || j > 0 && !DirBC_Interior)
                art_bottom = art_gpu(r[j - 1], theta[i], 0, mod_pk, prob, kappa_eps, delta_e, Rmax);
            art_left  = art_gpu(r[j], theta[im], 0, mod_pk, prob, kappa_eps, delta_e, Rmax);
            art_right = art_gpu(r[j], theta[ip], 0, mod_pk, prob, kappa_eps, delta_e, Rmax);
            if (j < nr_int - 1)
                art_top = art_gpu(r[j + 1], theta[i], 0, mod_pk, prob, kappa_eps, delta_e, Rmax);
        }

        // 9-Point Stencil (attention order: bottom, bottom left, left, top left, top, top right, right, bottom right, middle)
        // Bottom
        if (smoother < 2 && j < delete_circles && j > 0 || smoother > 1 && j == delete_circles) {
            //for the circle smoother (but not the inner circle)
            //and for the radial smoother for the first line, which has a bottom link to the circle smoother

            // j=1 means r-h is on the boundary
            if (j > 1 || !DirBC_Interior) {
                // second circle (at least one 'real' circle with ntheta_int nodes which is closer to the origin)
                col = row - ntheta_int;
                val = -0.5 * (kt + ktmin1) * (arr_bottom + arr) / hsmin1;
                Au[row_index] -= val * u[col];
            }
        }

        // Bottom-Left
        if (mod_pk && (j > 1 || j > 0 && !DirBC_Interior)) {
            col = row - (i > 0) * ntheta_int - 1;
            val = -0.25 * (art_left + art_bottom);
            Au[row_index] -= val * u[col];
        }

        // Left
        if (smoother > 1 && j >= delete_circles) {
            col = row + (i == 0) * ntheta_int - 1;
            val = -0.5 * (hsmin1 + hs) * (att_left + att) / ktmin1;
            Au[row_index] -= val * u[col];
        }

        // Top-Left
        if (mod_pk && j < nr_int - 1) {
            col = row + (1 + (i == 0)) * ntheta_int - 1;
            val = 0.25 * (art_left + art_top);
            Au[row_index] -= val * u[col];
        }

        // Top
        if (smoother < 2 && j < delete_circles) {
            // means that r[j+1] is not on the Dirichlet boundary
            if (j < nr_int - 1) {
                col = row + ntheta_int;
                val = -0.5 * (kt + ktmin1) * (arr + arr_top) / hs;
                Au[row_index] -= val * u[col];
            }
        }

        // Top-Right
        // means that r[j+1] is not on the Dirichlet boundary
        if (mod_pk && j < nr_int - 1) {
            col = row + (i + 1 < ntheta_int) * ntheta_int + 1;
            val = -0.25 * (art_top + art_right);
            Au[row_index] -= val * u[col];
        }

        // Right
        if (smoother > 1 && j >= delete_circles) {
            col = row - (i == ntheta_int - 1) * ntheta_int + 1;
            val = -0.5 * (hs + hsmin1) * (att + att_right) / kt;
            Au[row_index] -= val * u[col];
        }

        // Bottom-Right
        if (mod_pk && (j > 1 || j > 0 && !DirBC_Interior)) {
            col = row - (1 + (i == ntheta_int - 1)) * ntheta_int + 1;
            val = 0.25 * (art_bottom + art_right);
            Au[row_index] -= val * u[col];
        }
    }
} /* ----- end of level::apply_Asc_ortho0 ----- */

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
__device__ int get_ptr_sc(int i, int j, int smoother, int ortho, int delete_circles, int nr_int, int ntheta_int)
{
    if (j >= delete_circles && smoother < 2 || j < delete_circles && smoother > 1)
        printf("(get_ptr_sc) Incompatible radius and smoother.");

    int ptr = 0;
    int nr_left;
    if (j < delete_circles) {
        if (smoother == 0) {
            if (j > 0) {
                ptr += ntheta_int;

                nr_left = ceil((double)(j - 1) * 0.5) - 1;
                ptr += nr_left * ntheta_int;
            }
        }
        else if (smoother == 1) {
            if (j > 1) {
                ptr += ntheta_int;

                nr_left = floor((double)(j - 1) * 0.5) - 1;
                ptr += nr_left * ntheta_int;
            }
        }
        // In case of extrapol and smoother==0, there are only half of points in the grid: take twice
        // the same values for each in order to respect the theta index i
        ptr += i;
    }
    else {
        // Radial-Circle nodes
        if (j > delete_circles) {
            ptr += ntheta_int * 0.5;

            // Radial nodes
            int j1   = j - delete_circles - 1;
            int j2   = nr_int - delete_circles - 2;
            int j_nb = (j1 < j2) ? j1 : j2;
            ptr += j_nb * ntheta_int / 2;
        }
        // Last circle
        if (j > nr_int - 1) {
            ptr += ntheta_int / 2;
        }

        ptr += i / 2;
    }

    return ptr;
} /* ----- end of get_ptr_sc ----- */

/*! allocate memory for the GPU */
void level::alloc_cuda()
{
    cudaMalloc((void**)&dev_u, m * sizeof(double));
    cudaMalloc((void**)&dev_Au, m * sizeof(double));
    cudaMalloc((void**)&dev_fVec, m * sizeof(double));
    cudaMalloc((void**)&dev_res, m * sizeof(double));

    cudaMalloc((void**)&dev_r, r.size() * sizeof(double));
    cudaMalloc((void**)&dev_hplus, hplus.size() * sizeof(double));
    cudaMalloc((void**)&dev_theta, theta.size() * sizeof(double));
    cudaMalloc((void**)&dev_thetaplus, thetaplus.size() * sizeof(double));

    cudaMalloc((void**)&dev_error, m * sizeof(double));
    cudaMalloc((void**)&dev_error_fine, m * sizeof(double));

    if (row_Ac_LU.size() == 0) {
        //printf("ALLOC: level %d \n", l);
        cudaMalloc((void**)&dev_fsc[0], m_sc[0] * sizeof(double));
        cudaMalloc((void**)&dev_fsc[1], m_sc[1] * sizeof(double));
        cudaMalloc((void**)&dev_fsc[2], m_sc[2] * sizeof(double));
        cudaMalloc((void**)&dev_fsc[3], m_sc[3] * sizeof(double));

        // // TODO: SIZE=M_SC
        // cudaMalloc((void**)&dev_f_total, m * sizeof(double));

        cudaMalloc((void**)&dev_coarse_nodes_list_r, coarse_nodes_list_r.size() * sizeof(int));
        cudaMalloc((void**)&dev_coarse_nodes_list_theta, coarse_nodes_list_theta.size() * sizeof(int));
        cudaMalloc((void**)&dev_coarse_nodes_list_type_shift, sizeof(int) * coarse_nodes_list_type_shift.size());
        cudaMalloc((void**)&dev_coarse_nodes_list_type_start, sizeof(int) * coarse_nodes_list_type_start.size());

        cudaMalloc((void**)&dev_coarse_nodes_list_type0, sizeof(int) * coarse_nodes_list_type[0].size());
        cudaMalloc((void**)&dev_coarse_nodes_list_type1, sizeof(int) * coarse_nodes_list_type[1].size());
        cudaMalloc((void**)&dev_coarse_nodes_list_type2, sizeof(int) * coarse_nodes_list_type[2].size());
        cudaMalloc((void**)&dev_coarse_nodes_list_type3, sizeof(int) * coarse_nodes_list_type[3].size());
        cudaMalloc((void**)&dev_coarse_nodes_list_type4, sizeof(int) * coarse_nodes_list_type[4].size());
        cudaMalloc((void**)&dev_coarse_nodes_list_type5, sizeof(int) * coarse_nodes_list_type[5].size());
        cudaMalloc((void**)&dev_coarse_nodes_list_type6, sizeof(int) * coarse_nodes_list_type[6].size());
        cudaMalloc((void**)&dev_coarse_nodes_list_type7, sizeof(int) * coarse_nodes_list_type[7].size());

        cudaMalloc((void**)&dev_coarse_nodes_list_i0, sizeof(int) * coarse_nodes_list_i[0].size());
        cudaMalloc((void**)&dev_coarse_nodes_list_i1, sizeof(int) * coarse_nodes_list_i[1].size());
        cudaMalloc((void**)&dev_coarse_nodes_list_i2, sizeof(int) * coarse_nodes_list_i[2].size());
        cudaMalloc((void**)&dev_coarse_nodes_list_i3, sizeof(int) * coarse_nodes_list_i[3].size());
        cudaMalloc((void**)&dev_coarse_nodes_list_i4, sizeof(int) * coarse_nodes_list_i[4].size());
        cudaMalloc((void**)&dev_coarse_nodes_list_i5, sizeof(int) * coarse_nodes_list_i[5].size());
        cudaMalloc((void**)&dev_coarse_nodes_list_i6, sizeof(int) * coarse_nodes_list_i[6].size());
        cudaMalloc((void**)&dev_coarse_nodes_list_i7, sizeof(int) * coarse_nodes_list_i[7].size());
        cudaMalloc((void**)&dev_coarse_nodes_list_j0, sizeof(int) * coarse_nodes_list_j[0].size());
        cudaMalloc((void**)&dev_coarse_nodes_list_j1, sizeof(int) * coarse_nodes_list_j[1].size());
        cudaMalloc((void**)&dev_coarse_nodes_list_j2, sizeof(int) * coarse_nodes_list_j[2].size());
        cudaMalloc((void**)&dev_coarse_nodes_list_j3, sizeof(int) * coarse_nodes_list_j[3].size());
        cudaMalloc((void**)&dev_coarse_nodes_list_j4, sizeof(int) * coarse_nodes_list_j[4].size());
        cudaMalloc((void**)&dev_coarse_nodes_list_j5, sizeof(int) * coarse_nodes_list_j[5].size());
        cudaMalloc((void**)&dev_coarse_nodes_list_j6, sizeof(int) * coarse_nodes_list_j[6].size());
        cudaMalloc((void**)&dev_coarse_nodes_list_j7, sizeof(int) * coarse_nodes_list_j[7].size());
    }
}

/*! Copy vectors to the GPU */
void level::copy_cuda()
{
    cudaMemcpy(dev_r, &(r[0]), r.size() * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_hplus, &(hplus[0]), hplus.size() * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_theta, &(theta[0]), theta.size() * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_thetaplus, &(thetaplus[0]), thetaplus.size() * sizeof(double), cudaMemcpyHostToDevice);
    if (l == 0) {
        int nblocks = ceil((double)m / NTHREADS_CUDA);
        empty_vect_gpu<<<nblocks, NTHREADS_CUDA>>>(m, dev_u);
        cudaMemcpy(dev_fVec, &(fVec[0]), m * sizeof(double), cudaMemcpyHostToDevice);
    }
    if (row_Ac_LU.size() == 0) {
        cudaMemcpy(dev_coarse_nodes_list_r, &coarse_nodes_list_r[0], coarse_nodes_list_r.size() * sizeof(int),
                   cudaMemcpyHostToDevice);
        cudaMemcpy(dev_coarse_nodes_list_theta, &coarse_nodes_list_theta[0],
                   coarse_nodes_list_theta.size() * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(dev_coarse_nodes_list_type_shift, &(coarse_nodes_list_type_shift[0]),
                   sizeof(int) * coarse_nodes_list_type_shift.size(), cudaMemcpyHostToDevice);
        cudaMemcpy(dev_coarse_nodes_list_type_start, &(coarse_nodes_list_type_start[0]),
                   sizeof(int) * coarse_nodes_list_type_start.size(), cudaMemcpyHostToDevice);
        cudaMemcpy(dev_coarse_nodes_list_type0, &(coarse_nodes_list_type[0][0]),
                   sizeof(int) * coarse_nodes_list_type[0].size(), cudaMemcpyHostToDevice);
        cudaMemcpy(dev_coarse_nodes_list_type1, &(coarse_nodes_list_type[1][0]),
                   sizeof(int) * coarse_nodes_list_type[1].size(), cudaMemcpyHostToDevice);
        cudaMemcpy(dev_coarse_nodes_list_type2, &(coarse_nodes_list_type[2][0]),
                   sizeof(int) * coarse_nodes_list_type[2].size(), cudaMemcpyHostToDevice);
        cudaMemcpy(dev_coarse_nodes_list_type3, &(coarse_nodes_list_type[3][0]),
                   sizeof(int) * coarse_nodes_list_type[3].size(), cudaMemcpyHostToDevice);
        cudaMemcpy(dev_coarse_nodes_list_type4, &(coarse_nodes_list_type[4][0]),
                   sizeof(int) * coarse_nodes_list_type[4].size(), cudaMemcpyHostToDevice);
        cudaMemcpy(dev_coarse_nodes_list_type5, &(coarse_nodes_list_type[5][0]),
                   sizeof(int) * coarse_nodes_list_type[5].size(), cudaMemcpyHostToDevice);
        cudaMemcpy(dev_coarse_nodes_list_type6, &(coarse_nodes_list_type[6][0]),
                   sizeof(int) * coarse_nodes_list_type[6].size(), cudaMemcpyHostToDevice);
        cudaMemcpy(dev_coarse_nodes_list_type7, &(coarse_nodes_list_type[7][0]),
                   sizeof(int) * coarse_nodes_list_type[7].size(), cudaMemcpyHostToDevice);
        cudaMemcpy(dev_coarse_nodes_list_i0, &(coarse_nodes_list_i[0][0]), sizeof(int) * coarse_nodes_list_i[0].size(),
                   cudaMemcpyHostToDevice);
        cudaMemcpy(dev_coarse_nodes_list_i1, &(coarse_nodes_list_i[1][0]), sizeof(int) * coarse_nodes_list_i[1].size(),
                   cudaMemcpyHostToDevice);
        cudaMemcpy(dev_coarse_nodes_list_i2, &(coarse_nodes_list_i[2][0]), sizeof(int) * coarse_nodes_list_i[2].size(),
                   cudaMemcpyHostToDevice);
        cudaMemcpy(dev_coarse_nodes_list_i3, &(coarse_nodes_list_i[3][0]), sizeof(int) * coarse_nodes_list_i[3].size(),
                   cudaMemcpyHostToDevice);
        cudaMemcpy(dev_coarse_nodes_list_i4, &(coarse_nodes_list_i[4][0]), sizeof(int) * coarse_nodes_list_i[4].size(),
                   cudaMemcpyHostToDevice);
        cudaMemcpy(dev_coarse_nodes_list_i5, &(coarse_nodes_list_i[5][0]), sizeof(int) * coarse_nodes_list_i[5].size(),
                   cudaMemcpyHostToDevice);
        cudaMemcpy(dev_coarse_nodes_list_i6, &(coarse_nodes_list_i[6][0]), sizeof(int) * coarse_nodes_list_i[6].size(),
                   cudaMemcpyHostToDevice);
        cudaMemcpy(dev_coarse_nodes_list_i7, &(coarse_nodes_list_i[7][0]), sizeof(int) * coarse_nodes_list_i[7].size(),
                   cudaMemcpyHostToDevice);
        cudaMemcpy(dev_coarse_nodes_list_j0, &(coarse_nodes_list_j[0][0]), sizeof(int) * coarse_nodes_list_j[0].size(),
                   cudaMemcpyHostToDevice);
        cudaMemcpy(dev_coarse_nodes_list_j1, &(coarse_nodes_list_j[1][0]), sizeof(int) * coarse_nodes_list_j[1].size(),
                   cudaMemcpyHostToDevice);
        cudaMemcpy(dev_coarse_nodes_list_j2, &(coarse_nodes_list_j[2][0]), sizeof(int) * coarse_nodes_list_j[2].size(),
                   cudaMemcpyHostToDevice);
        cudaMemcpy(dev_coarse_nodes_list_j3, &(coarse_nodes_list_j[3][0]), sizeof(int) * coarse_nodes_list_j[3].size(),
                   cudaMemcpyHostToDevice);
        cudaMemcpy(dev_coarse_nodes_list_j4, &(coarse_nodes_list_j[4][0]), sizeof(int) * coarse_nodes_list_j[4].size(),
                   cudaMemcpyHostToDevice);
        cudaMemcpy(dev_coarse_nodes_list_j5, &(coarse_nodes_list_j[5][0]), sizeof(int) * coarse_nodes_list_j[5].size(),
                   cudaMemcpyHostToDevice);
        cudaMemcpy(dev_coarse_nodes_list_j6, &(coarse_nodes_list_j[6][0]), sizeof(int) * coarse_nodes_list_j[6].size(),
                   cudaMemcpyHostToDevice);
        cudaMemcpy(dev_coarse_nodes_list_j7, &(coarse_nodes_list_j[7][0]), sizeof(int) * coarse_nodes_list_j[7].size(),
                   cudaMemcpyHostToDevice);
    }
} /* ----- end of copy_cuda ----- */

/*!
 *  \brief Default Destructor of level class
 *
 *  Default destructor of the level class.
 *
 */
level::~level()
{

    cudaFree(dev_fVec);
    // cudaFree(dev_f_total);
    cudaFree(dev_u);
    cudaFree(dev_r);
    cudaFree(dev_hplus);
    cudaFree(dev_theta);
    cudaFree(dev_thetaplus);
    cudaFree(dev_Au);
    cudaFree(dev_res);
    cudaFree(dev_error);
    cudaFree(dev_error_fine);

    if (row_Ac_LU.size() == 0) {
        //printf("DESTRUCTOR: level %d \n", l);
        cudaFree(dev_fsc[0]);
        cudaFree(dev_fsc[1]);
        cudaFree(dev_fsc[2]);
        cudaFree(dev_fsc[3]);
        cudaFree(dev_coarse_nodes_list_r);
        cudaFree(dev_coarse_nodes_list_theta);
        cudaFree(dev_coarse_nodes_list_type_shift);
        cudaFree(dev_coarse_nodes_list_type_start);
        cudaFree(dev_coarse_nodes_list_type0);
        cudaFree(dev_coarse_nodes_list_type1);
        cudaFree(dev_coarse_nodes_list_type2);
        cudaFree(dev_coarse_nodes_list_type3);
        cudaFree(dev_coarse_nodes_list_type4);
        cudaFree(dev_coarse_nodes_list_type5);
        cudaFree(dev_coarse_nodes_list_type6);
        cudaFree(dev_coarse_nodes_list_type7);
        cudaFree(dev_coarse_nodes_list_i0);
        cudaFree(dev_coarse_nodes_list_i1);
        cudaFree(dev_coarse_nodes_list_i2);
        cudaFree(dev_coarse_nodes_list_i3);
        cudaFree(dev_coarse_nodes_list_i4);
        cudaFree(dev_coarse_nodes_list_i5);
        cudaFree(dev_coarse_nodes_list_i6);
        cudaFree(dev_coarse_nodes_list_i7);
        cudaFree(dev_coarse_nodes_list_j0);
        cudaFree(dev_coarse_nodes_list_j1);
        cudaFree(dev_coarse_nodes_list_j2);
        cudaFree(dev_coarse_nodes_list_j3);
        cudaFree(dev_coarse_nodes_list_j4);
        cudaFree(dev_coarse_nodes_list_j5);
        cudaFree(dev_coarse_nodes_list_j6);
        cudaFree(dev_coarse_nodes_list_j7);
    }

} /* ----- end of destructor level::~level ----- */

/*!
 *  \brief Solves the problem gyro with multigrid on levels
 *
 *  Solves the problem gyro with multigrid on levels
 *
 */
void gmgpolar::polar_multigrid()
{
    if (gyro::icntl[Param::verbose] > 1)
        std::cout << "Define the coarse grids...\n";
    check_geom();
    define_coarse_nodes();

    if (gyro::icntl[Param::verbose] > 2) {
        std::cout << "\n";
        if (gyro::icntl[Param::prob] == 4 || gyro::icntl[Param::prob] == 1)
            std::cout << "Solving the POISSON EQUATION, PROBLEM " << gyro::icntl[Param::prob] << "\n";
        else if (gyro::icntl[Param::prob] == 5)
            std::cout << "Solving the ARCY EQUATION, PROBLEM " << gyro::icntl[Param::prob]
                      << " with arctan coefficient profile\n";
        if (!gyro::icntl[Param::mod_pk])
            std::cout << "Considering POLAR coordinates with ";
        else
            std::cout << "Considering MODIFIED POLAR coordinates with kappa=" << gyro::dcntl[Param::kappa_eps]
                      << ", delta=" << gyro::dcntl[Param::delta_e] << ", and ";
        std::cout << gyro::dcntl[Param::R0] << " <= r <= " << gyro::dcntl[Param::R] << "\n";
        std::cout << "Using 9 point star FINITE DIFFERENCES.\n";
        std::cout << "Using FULL EXTRAPOLATION\n\n";
    }

    gyro::dcntl[Param::r0_DB] = -1e6;
    if (gyro::icntl[Param::DirBC_Interior])
        gyro::dcntl[Param::r0_DB] = gyro::dcntl[Param::R0];

    if (gyro::icntl[Param::verbose] > 1)
        std::cout
            << "Building discretized system, restriction and interpolation operators, and defining splittings...\n";
    prepare_op_levels();
    //std::cout << "...done.\n\n";

    int m = v_level[0]->m;

    std::string cycle_str = "?";
    if (gyro::icntl[Param::cycle] == 1)
        cycle_str = "V";
    else if (gyro::icntl[Param::cycle] == 2)
        cycle_str = "W";
    //if (gyro::icntl[Param::verbose] > 1)
    std::cout << "\n***** Problem size " << m << " (" << v_level[0]->nr << ", " << v_level[0]->ntheta << "), " << levels
              << " grids, nr_exp=" << gyro::icntl[Param::nr_exp] << ", aniso=" << gyro::icntl[Param::fac_ani]
              << " ***** smoother=" << gyro::icntl[Param::smoother] << ", r0=" << gyro::dcntl[Param::R0]
              << ", mod_pk=" << gyro::icntl[Param::mod_pk] << ", DirBC=" << gyro::icntl[Param::DirBC_Interior]
              << ", divide=" << gyro::icntl[Param::divideBy2] << " *****\n";

    double scaling = 1.0;
    if (gyro::icntl[Param::compute_rho])
        scaling = scaling / (sqrt(m));

    v_level[0]->u.assign(m, scaling); //create an empty vector u

    //allocate memory on the GPU for every level
    for (int i = 0; i < levels; i++) {
        v_level[i]->alloc_cuda();
        // copy_cuda
        // u, fVec, f_sc
        // 0: Au, fsc, error_fine
        v_level[i]->copy_cuda();
    }

    if (gyro::icntl[Param::debug])
        debug();
    // else {
    if (gyro::icntl[Param::verbose] > 1)
        std::cout << "Multigrid iteration....!\n\n";
    multigrid_iter();

    if (gyro::icntl[Param::verbose] > 1)
        for (int i = 0; i < levels; i++) {
            std::cout << "LEVEL " << i << "\n";
            std::cout << "\nt_smoothing: " << v_level[i]->t_smoothing << ", t_f_sc: " << v_level[i]->t_f_sc
                      << ", t_Asc_ortho: " << v_level[i]->t_Asc_ortho << ", t_Asc: " << v_level[i]->t_Asc
                      << ", t_indices: " << v_level[i]->t_indices << "\n";
            std::cout << "\nt_get_ptr: " << v_level[i]->t_get_ptr << ", t_get_stencil: " << v_level[i]->t_get_stencil
                      << ", t_get_smoother: " << v_level[i]->t_get_smoother << ", t_get_row: " << v_level[i]->t_get_row
                      << "\n";
            std::cout << "\n";
        }

    if (gyro::icntl[Param::verbose] > 0) {
        std::cout << "\nt_setup: " << t_setup << ", t_build: " << t_build << ", t_facto_Ac: " << t_facto_Ac
                  << ", t_build_P: " << t_build_P << ", t_build_Asc: " << t_build_Asc
                  << ", t_facto_Asc: " << t_facto_Asc << "\n";
        std::cout << "t_total (fine): " << t_total << ", t_smoothing: " << t_smoothing << ", t_residual: " << t_residual
                  << ", t_restriction: " << t_restriction << ", t_Ac: " << t_Ac
                  << ", t_prolongation: " << t_prolongation << ", t_fine_residual: " << t_fine_residual
                  << ", t_error: " << t_error << "\n";
    }
    // }
} /* ----- end of gmgpolar::polar_multigrid ----- */

void facto_gaussian_elimination_gpu(int* row_ind, int* col_ind, double* val, int nnz, int m, double* rhs, double* sol)
{

    cusolverSpHandle_t cusolver_handle;
    cusolverStatus_t cusolver_status;
    cusolver_status = cusolverSpCreate(&cusolver_handle);
    //  std::cout << "status create cusolver handle: " << cusolver_status << std::endl;

    cusparseHandle_t cusparse_handle;
    cusparseStatus_t cusparse_status;
    cusparse_status = cusparseCreate(&cusparse_handle);
    //  std::cout << "status create cusparse handle: " << cusparse_status << std::endl;

    double tol      = 1e-9;
    int reorder     = 1;
    int singularity = 0;

    double *db, *dcsrVal, *dx;
    int *dcsrColInd, *dcsrRowPtr, *dcooRow;
    cudaMalloc((void**)&db, m * sizeof(double));
    cudaMalloc((void**)&dx, m * sizeof(double));
    cudaMalloc((void**)&dcsrVal, nnz * sizeof(double));
    cudaMalloc((void**)&dcsrColInd, nnz * sizeof(int));
    cudaMalloc((void**)&dcsrRowPtr, (m + 1) * sizeof(int));
    cudaMalloc((void**)&dcooRow, nnz * sizeof(int));

    cudaMemcpy(db, rhs, m * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dcsrVal, val, nnz * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dcsrColInd, col_ind, nnz * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dcooRow, row_ind, nnz * sizeof(int), cudaMemcpyHostToDevice);

    /* Convert matrix from coo format to csr */
    // size of matrix mxm
    // nnz: length of row
    cusparse_status = cusparseXcoo2csr(cusparse_handle, dcooRow, nnz, m, dcsrRowPtr, CUSPARSE_INDEX_BASE_ZERO);
    if (cusparse_status != 0) {
        std::cout << "facto_gaussian_elimination_gpu: cusparseXcoo2csr returned with code:" << cusparse_status
                  << std::endl;
    }
    /* Solve system with QR factorization */

    // First: create a descriptor
    cusparseMatDescr_t descrA = NULL;
    cusparseCreateMatDescr(&descrA);
    cusparseSetMatIndexBase(descrA, CUSPARSE_INDEX_BASE_ZERO); // A is base-1
    cusparseSetMatType(descrA, CUSPARSE_MATRIX_TYPE_GENERAL); // A is a general matrix

    cusolver_status = cusolverSpDcsrlsvqr(cusolver_handle, m, nnz, descrA, dcsrVal, dcsrRowPtr, dcsrColInd, db, tol,
                                          reorder, dx, &singularity);

    if (cusolver_status != 0) {
        std::cout << "facto_gaussian_elimination_gpu: cusolverSpDcsrlsvqr returned with code:" << cusolver_status
                  << std::endl;
    }

    if (singularity != -1) {
        printf("singularity detected at i = %d\n", singularity);
    }
    cudaMemcpy(sol, dx, m * sizeof(double), cudaMemcpyDeviceToHost);

    cusparse_status = cusparseDestroy(cusparse_handle);
    // std::cout << "status destroy cusparse handle: " << cusparse_status << std::endl;

    cusolver_status = cusolverSpDestroy(cusolver_handle);
    //  std::cout << "status destroy cusolver handle: " << cusolver_status << std::endl;

    cudaFree(db);
    cudaFree(dx);
    cudaFree(dcsrVal);
    cudaFree(dcsrColInd);
    cudaFree(dcsrRowPtr);
    cudaFree(dcooRow);

    return;
}

#endif
