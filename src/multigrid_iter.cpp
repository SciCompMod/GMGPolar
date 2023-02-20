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
 * \file multigrid_iter.cpp
 * \brief Implementation of the multigrid scheme
 * \author M. Kuehn, C. Kruse, P. Leleux, C. Schwarz
 * \version 0.0
 */

#include "gmgpolar.h"
#include <unistd.h>

/*!
 *  \brief The multigrid cycle iterations
 *
 *  The multigrid cycle iterations: calls gmgpolar::multigrid_iter_extrapol()
 *
 */
void gmgpolar::multigrid_iter()
{
    double t;
    TIC;

    int extrapol = gyro::icntl[Param::extrapolation];
    int nrm_res  = gyro::icntl[Param::res_norm];

    if (extrapol > 0)
        v_level[1]->fVec_initial = v_level[1]->fVec; //Extrapolation: store the initial fVec on level 1
    if (gyro::icntl[Param::verbose] > 0) {
        if (extrapol == 1)
            std::cout << "WITH IMPLICIT EXTRAPOLATION!!!\n";
        else if (extrapol == 2)
            std::cout << "WITH ALTERNATIVE EXTRAPOLATION!!!\n";
    }

    int it = 0;
    v_level[0]->u.assign(v_level[0]->m, 0); //zero u on level 0 (only once at the beginning)

    TIC;

    if (extrapol < 2) {

        //! compute the initial residual: res0 = f - A*u
        compute_residual(0, extrapol); //compute residual on level 0
        //compute the 2 norm and inf-norm of the residual
        nrm_2_res.push_back(0); //2norm of residual, store the residual norm on level 0 of every iteration
        nrm_inf_res.push_back(0);
        double nrm_inf_err_temp = 0;

        for (long unsigned int i = 0; i < v_level[0]->res.size(); ++i) {
            nrm_2_res[0] += v_level[0]->res[i] * v_level[0]->res[i];

            if (fabs(v_level[0]->res[i]) > nrm_inf_err_temp) {
                nrm_inf_err_temp = fabs(v_level[0]->res[i]);
            }
        }
        nrm_2_res[0]   = sqrt(nrm_2_res[0]);
        nrm_inf_res[0] = nrm_inf_err_temp;

        if (gyro::icntl[Param::verbose] > 0) {
            std::cout << "initial residual: 2-norm = " << nrm_2_res[0] << "\n";
            std::cout << "initial residual: inf-norm = " << nrm_inf_res[0] << "\n";
        }
    }
    //ALternative extrapolation
    else {
        if (gyro::icntl[Param::check_error] == 0)
            throw std::runtime_error("The alternative extrapolation technique requires to check the error, w.r.t. to "
                                     "the theoretical solution, as stopping criterion.");
        //! compute the initial error
        std::vector<double> error = compute_error();
        nrm_2_err.push_back(0);
        for (std::size_t i = 0; i < error.size(); ++i) { //iterate over the grid in polar coordinates
            nrm_2_err[0] += error[i] * error[i];
        }
        nrm_2_err[0] = sqrt(nrm_2_err[0]) / sqrt(v_level[0]->m); //scaling by 1/sqrt(m)
        if (gyro::icntl[Param::verbose] > 0)
            std::cout << "initial error: 2-norm = " << nrm_2_err[0] << "\n";
    }

    double rel_red_conv          = gyro::dcntl[Param::rel_red_conv]; //threshold on relative residual
    double convergence_criterium = 1.0;

    t_fine_residual += TOC;
    TIC;

    //! Start the Multigrid-Iteration
    while (it < gyro::icntl[Param::maxiter] && convergence_criterium > rel_red_conv) {
        it++;

        //call the multigrid_cycle on level 0 (the finest level)
        multigrid_cycle_extrapol(0);

        if (gyro::icntl[Param::verbose] > 3)
            gyro::disp(v_level[0]->u, "u");

        TIC;
        //----------------------------------------------------------------------------------------------------------
        //! compute the convergence criterium
        //use the residual as convergence criterium
        if (extrapol < 2) {
            gmgpolar::compute_residual(0, extrapol); //compute residual on level 0

            if (gyro::icntl[Param::verbose] > 3)
                gyro::disp(v_level[0]->res, "res");

            nrm_2_res.push_back(0);
            nrm_inf_res.push_back(0);
            double nrm_inf_res_temp = 0;

            for (unsigned long int i = 0; i < v_level[0]->res.size(); ++i) {
                nrm_2_res[it] += v_level[0]->res[i] * v_level[0]->res[i];

                if (fabs(v_level[0]->res[i]) > nrm_inf_res_temp) {
                    nrm_inf_res_temp = fabs(v_level[0]->res[i]);
                }
            }
            nrm_2_res[it]   = sqrt(nrm_2_res[it]);
            nrm_inf_res[it] = nrm_inf_res_temp;

            //  * Defines the norm used for the residual in the stopping criterion
            //  * 0:
            //  * 1:
            //  * 2:
            //  * 3:
            if (nrm_res == 0) { // L2 norm scaled by initial res.
                convergence_criterium = nrm_2_res[it] / nrm_2_res[0];
            }
            else if (nrm_res == 1) { // Infinitiy norm scaled by initial res.
                convergence_criterium = nrm_inf_res[it] / nrm_inf_res[0];
            }
            else if (nrm_res == 2) { // L2 norm
                convergence_criterium = nrm_2_res[it];
            }
            else if (nrm_res == 3) { // Infinitiy norm
                convergence_criterium = nrm_inf_res[it];
            }

            if (gyro::icntl[Param::verbose] > 0)
                std::cout << "--> Iteration " << it << ": residual norm = " << nrm_2_res[it]
                          << ", relative residual = " << convergence_criterium << std::endl;
        }
        //Alternative extrapolation: use error instead of residual as convergence criterium
        else {
            if (gyro::icntl[Param::check_error] == 0)
                throw std::runtime_error(
                    "The alternative extrapolation technique requires to check the error, w.r.t. to "
                    "the theoretical solution, as stopping criterion.");

            std::vector<double> error = compute_error();

            if (gyro::icntl[Param::verbose] > 3)
                gyro::disp(error, "error");

            nrm_2_err.push_back(0);
            for (std::size_t i = 0; i < error.size(); ++i) { //iterate over the grid in polar coordinates
                nrm_2_err[it] += error[i] * error[i];
            }
            nrm_2_err[it] = sqrt(nrm_2_err[it]) / sqrt(v_level[0]->m); //scaling by 1/sqrt(m)

            double error_difference = fabs(nrm_2_err[it] - nrm_2_err[it - 1]);
            convergence_criterium   = error_difference / nrm_2_err[0];

            if (gyro::icntl[Param::verbose] > 0)
                std::cout << "--> Iteration " << it << ": error norm = " << nrm_2_err[it]
                          << ", relative error = " << convergence_criterium << std::endl;
        }
        t_fine_residual += TOC;
        TIC;
    }
    // if (gyro::icntl[Param::verbose] > 0) {
    if (it == gyro::icntl[Param::maxiter])
        std::cout << "Multigrid reached maxiter=" << gyro::icntl[Param::maxiter] << "\n";
    else
        std::cout << "Convergence after iteration " << it << std::endl;
    // }
    //----------------------------------------------------------------------------------------------------------
    //!compute mean residual reduction factor rho
    if (extrapol < 2) {
        double rho_mean = std::pow(convergence_criterium, 1.0 / it); //take the it-th root (it=number of iterations)
        std::cout << "mean residual reduction factor: rho = " << rho_mean << std::endl;
    }

    if (!gyro::f_sol_out.empty()) {
        v_level[0]->write_sol();
    }

    TIC;
    if (gyro::icntl[Param::check_error] == 1) {
        //compute the error in the 2-norm and inf-norm
        std::vector<double> error = compute_error();
        double nrm_inf_err        = 0;
        double nrm_2              = 0;
        for (std::size_t i = 0; i < error.size(); ++i) { //iterate over the grid in polar coordinates
            nrm_2 += error[i] * error[i];
            if (fabs(error[i]) > nrm_inf_err) {
                nrm_inf_err = fabs(error[i]);
            }
        }
        nrm_2 = sqrt(nrm_2) / sqrt(v_level[0]->m); //scaling by 1/sqrt(m)

        std::cout << "2-norm of error = " << nrm_2 << std::endl;
        std::cout << "inf-norm of error = " << nrm_inf_err << std::endl;

        if (gyro::icntl[Param::verbose] > 3)
            gyro::disp(error, "error");
    }

    t_error += TOC;
    TIC;

    if (gyro::icntl[Param::verbose] > 3)
        gyro::disp(v_level[0]->u, "u");
} /* ----- end of gmgpolar::multigrid_iter ----- */

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

    int extrapol                             = gyro::icntl[Param::extrapolation];
    int s                                    = 0;
    int c                                    = 0;
    int smoother                             = 0;
    std::vector<std::vector<double>> f_Asc_u = std::vector<std::vector<double>>(4, std::vector<double>());

    TIC;
    //! v1 presmoothing steps (on level l)
    for (int v = 0; v < gyro::icntl[Param::v1]; ++v) {
        for (smoother = 0; smoother < 4; smoother++) {
            int size = 0, size_ortho = 0;
            if (smoother < 2) {
                size_ortho = v_level[l]->delete_circles + 1;
                size       = v_level[l]->delete_circles;
            }
            else if (smoother > 1) {
                size_ortho = v_level[l]->ntheta_int;
                size       = v_level[l]->ntheta_int;
            }
            for (int i = 0; i < size_ortho; i++)
                v_level[l]->dep_Asc_ortho[smoother][i] = 0;
            for (int i = 0; i < size; i++)
                v_level[l]->dep_Asc[smoother][i] = 0;
        }

        if (gyro::icntl[Param::optimized] == 0) {
            // iterate over the 4 smoothers (0:c/b, 1:c/w, 2:r/b, 3:r/w)
            for (int smoother = 0; smoother < 4; smoother++) {
                TIC;

                //call the smooothing function, the result is u_sc which is directly inserted into u
                v_level[l]->multigrid_smoothing0(smoother);

                if (gyro::icntl[Param::verbose] > 2)
                    std::cout << "SMOOTHER: " << smoother << " = " << TOC << "\n";
            }
        }
        else {
            // Initialize vectors to contain Asc_ortho.u (required for OpenMP, else adresses change between threads)
            for (smoother = 0; smoother < 4; smoother++)
                f_Asc_u[smoother] = std::vector<double>(v_level[l]->nblocks[smoother] * v_level[l]->m_sc[smoother], 0);

                //call the smooothing function, the result is u_sc which is directly inserted into u
#pragma omp parallel firstprivate(v, s, c, smoother) shared(f_Asc_u)
            {
#pragma omp single
                {
                    for (s = 0; s < 2; s++) {
                        // iterate over the 4 smoothers (0:c/b, 1:c/w, 2:r/b, 3:r/w)
                        for (c = 0; c < 2; c++) {
                            int smoother = s * 2 + c;
                            TIC;
                            int sm = (smoother < 2) ? 1 - smoother : 5 - smoother;
                            v_level[l]->multigrid_smoothing(
                                smoother, v, f_Asc_u[smoother], v_level[l]->nblocks[smoother], c,
                                v_level[l]->dep_Asc[smoother], v_level[l]->dep_Asc[sm], v_level[l]->dep_Asc[1],
                                v_level[l]->dep_Asc_ortho[smoother]);

                            if (gyro::icntl[Param::verbose] > 2)
                                std::cout << "SMOOTHER: " << smoother << " = " << TOC << "\n";
                        }
                    }
                } // omp single
            } // omp parallel
        }

        if (gyro::icntl[Param::smoother] == 13) {
            //create the current solution u from the two vectors u_previous
            //int number_circle_points = v_level[l]->m_sc[0] + v_level[l]->m_sc[1]; //does not work for extrapolation
            int number_circle_points = v_level[l]->delete_circles * v_level[l]->ntheta;
            for (int i = 0; i < number_circle_points; ++i) {
                v_level[l]->u[i] = v_level[l]->u_previous_c[i]; //insert values from the circle smoother
            }
            for (int i = number_circle_points; i < v_level[l]->m; ++i) {
                v_level[l]->u[i] = v_level[l]->u_previous_r[i]; //insert values from the radial smoother
            }
        }
    }
    //std::cout << "pre-smoothing done \n";

    if (gyro::icntl[Param::verbose] > 3)
        gyro::disp(v_level[l]->u, "u");

    t = t_smoothing_tmp;
    t_smoothing += TOC;
    TIC;

    //! compute residual (of level l)
    //even if we have extrapolation, compute just normal residual (extrapolation-restriction follows in the next step)
    gmgpolar::compute_residual(l, 0);

    if (gyro::icntl[Param::verbose] > 3)
        gyro::disp(v_level[l]->res, "res");

    t_residual += TOC;
    TIC;

    //! Restriction of residual (coarsening)
    //the restricted residual of level l becomes the right hand side on level l+1, f(l+1)=res_coarse(l)
    if (extrapol > 0 && l == 0) { //for extrapol = 1 and extrapol = 2
        std::vector<double> res_1;
        //res_1 = P_ex^T * res
        if (gyro::icntl[Param::matrix_free] == 1) {
            if (gyro::icntl[Param::optimized] == 0) {
                res_1 = v_level[l]->apply_prolongation_ex0(v_level[l]->res, v_level[l]->mc, v_level[l]->m,
                                                           v_level[l]->coarse_nodes_list_r,
                                                           v_level[l]->coarse_nodes_list_theta, 1);
            }
            else {
                res_1 = v_level[l]->apply_restriction_ex(v_level[l]->res);
            }
        }
        else {
            res_1 = std::vector<double>(v_level[l]->mc, 0);
            for (size_t i = 0; i < v_level[l]->ri_prol_ex.size(); i++) {
                res_1[v_level[l]->ci_prol_ex[i]] +=
                    v_level[l]->v_prol_ex[i] * v_level[l]->res[v_level[l]->ri_prol_ex[i]];
            }
        }

        if (gyro::icntl[Param::verbose] > 3)
            gyro::disp(res_1, "res_1");

        //Au_coarse = A(l+1) * u_coarse = A(l+1) * P_inj^T * u(l)
        std::vector<double> Au_coarse(v_level[l + 1]->m, 0); //result of apply_A(u_coarse), empty vector

        std::vector<double> u_coarse;
        if (gyro::icntl[Param::matrix_free] == 1) {
            if (gyro::icntl[Param::optimized] == 0) {
                u_coarse = v_level[l]->apply_prolongation_inj0(v_level[l]->u, v_level[l]->mc, v_level[l]->m,
                                                               v_level[l]->coarse_nodes_list_r,
                                                               v_level[l]->coarse_nodes_list_theta, 1);
            }
            else {
                u_coarse = v_level[l]->apply_restriction_inj(v_level[l]->u);
            }
        }
        else {
            u_coarse = std::vector<double>(v_level[l]->mc, 0);
            for (size_t i = 0; i < v_level[l]->ri_prol_inj.size(); i++) {
                u_coarse[v_level[l]->ci_prol_inj[i]] +=
                    v_level[l]->v_prol_inj[i] * v_level[l]->u[v_level[l]->ri_prol_inj[i]];
            }
        }

        if (gyro::icntl[Param::verbose] > 3)
            gyro::disp(u_coarse, "u_coarse");

        if (gyro::icntl[Param::matrix_free] == 1) {
            if (gyro::icntl[Param::optimized] == 0)
                v_level[l + 1]->apply_A0(u_coarse, Au_coarse);
            else
                v_level[l + 1]->apply_A(u_coarse, Au_coarse);
        }
        else {
            if (gyro::icntl[Param::openmp] == 1) {
                for (size_t i = 0; i < v_level[l + 1]->row_indices.size(); i++) {
                    Au_coarse[v_level[l + 1]->row_indices[i]] +=
                        v_level[l + 1]->vals[i] * u_coarse[v_level[l + 1]->col_indices[i]];
                }
            }
            else {
                int j;
#pragma omp parallel for private(j) firstprivate(u_coarse) shared(Au_coarse)
                for (j = 0; j < v_level[l + 1]->nr; j++) {
                    std::vector<int> ptr_vect = v_level[l + 1]->get_ptr(j);
                    int ptr_start             = ptr_vect[1];
                    int ptr_end               = ptr_start + (ptr_vect[2] - ptr_vect[1]) * v_level[l + 1]->ntheta_int;
                    for (int k = ptr_start; k < ptr_end; k++) {
                        Au_coarse[v_level[l + 1]->row_indices[k]] +=
                            v_level[l + 1]->vals[k] * u_coarse[v_level[l + 1]->col_indices[k]];
                    }
                }
            }
        }

        if (gyro::icntl[Param::verbose] > 3)
            gyro::disp(Au_coarse, "Au_coarse");

        // f(l+1) = 4/3 * res_1 - 1/3 * (f(l+1) - Au_coarse)
        for (unsigned long int i = 0; i < res_1.size(); ++i) {
            v_level[l + 1]->fVec[i] = 4. / 3. * res_1[i] - 1. / 3. * (v_level[l + 1]->fVec_initial[i] - Au_coarse[i]);
        }

        if (gyro::icntl[Param::verbose] > 3)
            gyro::disp(v_level[l + 1]->fVec, "fVec");
    }
    else { //no extrapolation
        //res(l+1)=P^T*res(l)      //trans=1 (Restriction, P^T)
        //P=(ncoarse*mc), ncoarse=m(size of fine level), mc(size of coarse level), coarse_r=coarse_nodes_list_r
        if (gyro::icntl[Param::matrix_free] == 1) {
            if (gyro::icntl[Param::optimized] == 0)
                v_level[l + 1]->fVec = v_level[l]->apply_prolongation_bi0(
                    v_level[l]->res, v_level[l]->mc, v_level[l]->m, v_level[l]->coarse_nodes_list_r,
                    v_level[l]->coarse_nodes_list_theta, 1);
            else
                v_level[l + 1]->fVec = v_level[l]->apply_restriction_bi(v_level[l]->res);
        }
        else {
            v_level[l + 1]->fVec = std::vector<double>(v_level[l]->mc, 0);
            for (size_t i = 0; i < v_level[l]->ri_prol.size(); i++) {
                v_level[l + 1]->fVec[v_level[l]->ci_prol[i]] +=
                    v_level[l]->v_prol[i] * v_level[l]->res[v_level[l]->ri_prol[i]];
            }
        }
    }
    //std::cout << "residual restricted \n";
    t_restriction += TOC;
    TIC;

    //! iterative call of multigrid_cycle_extrapol
    v_level[l + 1]->u.assign(v_level[l + 1]->m, 0); //zero u in every iteration
    std::vector<double> error_coarse;
    if (l == levels - 2) {
//solve exactly on the coarsest level for the error (A * error = res) (use whole A from coarsest level)
// check for the second coarsest level (levels-2), because we have no smoothing on the coarsest level (and thus no Asc and no Asc_ortho)
#ifdef USE_MUMPS
        if (gyro::icntl[Param::optimized] == 0) {
#endif
            error_coarse = v_level[l + 1]->solve_gaussian_elimination(
                v_level[l + 1]->row_Ac_LU, v_level[l + 1]->col_Ac_LU, v_level[l + 1]->vals_Ac_LU, v_level[l + 1]->fVec);
#ifdef USE_MUMPS
        }
        else
            error_coarse = v_level[l + 1]->solve_mumps(v_level[l + 1]->mumps_Ac, v_level[l + 1]->fVec);
#endif
        t_Ac += TOC;
        TIC;
    }
    else {
        multigrid_cycle_extrapol(l + 1);
        error_coarse = v_level[l + 1]->u; //the coarse_error on level l is u on level l+1
        TIC;
    }

    if (gyro::icntl[Param::verbose] > 2) {
        std::cout << "********************************************" << std::endl;
        std::cout << "BACK ON LEVEL " << l << std::endl;
        std::cout << "********************************************" << std::endl;
    }

    if (gyro::icntl[Param::verbose] > 3)
        gyro::disp(error_coarse, "error_coarse");

    TIC;
    //! Prolongation of error_coarse to error_fine
    //error(l) = P * error(l+1)
    std::vector<double> error_fine;
    if (extrapol > 0 && l == 0) { //for extrapol = 1 and extrapol = 2
        if (gyro::icntl[Param::matrix_free] == 1) {
            if (gyro::icntl[Param::optimized] == 0) {
                error_fine = v_level[l]->apply_prolongation_ex0(error_coarse, v_level[l]->mc, v_level[l]->m,
                                                                v_level[l]->coarse_nodes_list_r,
                                                                v_level[l]->coarse_nodes_list_theta, 0);
            }
            else {
                error_fine = v_level[l]->apply_prolongation_ex(error_coarse);
            }
        }
        else {
            error_fine = std::vector<double>(v_level[l]->m, 0);
            for (size_t i = 0; i < v_level[l]->ri_prol_ex.size(); i++) {
                error_fine[v_level[l]->ri_prol_ex[i]] +=
                    v_level[l]->v_prol_ex[i] * error_coarse[v_level[l]->ci_prol_ex[i]];
            }
        }
    }
    else { //no extrapolation
        if (gyro::icntl[Param::matrix_free] == 1) {
            if (gyro::icntl[Param::optimized] == 0)
                error_fine = v_level[l]->apply_prolongation_bi0(error_coarse, v_level[l]->mc, v_level[l]->m,
                                                                v_level[l]->coarse_nodes_list_r,
                                                                v_level[l]->coarse_nodes_list_theta, 0);
            else
                error_fine = v_level[l]->apply_prolongation_bi(error_coarse);
        }
        else {
            error_fine = std::vector<double>(v_level[l]->m, 0);
            for (size_t i = 0; i < v_level[l]->ri_prol.size(); i++) {
                error_fine[v_level[l]->ri_prol[i]] += v_level[l]->v_prol[i] * error_coarse[v_level[l]->ci_prol[i]];
            }
        }
    }
    //std::cout << "error prolongated \n";

    if (gyro::icntl[Param::verbose] > 3)
        gyro::disp(error_fine, "error_fine");

    //! correction of solution (on level l)
    //u += error
    for (long unsigned int i = 0; i < error_fine.size(); ++i) {
        v_level[l]->u[i] += error_fine[i];
    }
    //std::cout << "solution corrected \n";
    t_prolongation += TOC;
    TIC;

    t_smoothing_tmp = t;

    //! v2 postsmoothing steps (on level l)
    for (int v = 0; v < gyro::icntl[Param::v2]; ++v) {
        for (int smoother = 0; smoother < 4; smoother++) {
            int size = 0, size_ortho = 0;
            if (smoother < 2) {
                size_ortho = v_level[l]->delete_circles + 1;
                size       = v_level[l]->delete_circles;
            }
            else if (smoother > 1) {
                size_ortho = v_level[l]->ntheta_int;
                size       = v_level[l]->ntheta_int;
            }
            for (int i = 0; i < size_ortho; i++)
                v_level[l]->dep_Asc_ortho[smoother][i] = 0;
            for (int i = 0; i < size; i++)
                v_level[l]->dep_Asc[smoother][i] = 0;
        }

        //std::cout << "smoothing ... \n";
        if (gyro::icntl[Param::optimized] == 0) {
            // iterate over the 4 smoothers (0:c/b, 1:c/w, 2:r/b, 3:r/w)
            for (int smoother = 0; smoother < 4; smoother++) {
                TIC;

                //call the smooothing function, the result is u_sc which is directly inserted into u
                v_level[l]->multigrid_smoothing0(smoother);

                if (gyro::icntl[Param::verbose] > 2)
                    std::cout << "SMOOTHER: " << smoother << " = " << TOC << "\n";
            }
        }
        else {
            // Initialize vectors to contain Asc_ortho.u (required for OpenMP, else adresses change between threads)
            for (smoother = 0; smoother < 4; smoother++)
                f_Asc_u[smoother] = std::vector<double>(v_level[l]->nblocks[smoother] * v_level[l]->m_sc[smoother], 0);

                //call the smooothing function, the result is u_sc which is directly inserted into u
#pragma omp parallel firstprivate(v, s, c, smoother) shared(f_Asc_u)
            {
#pragma omp single
                {
                    for (s = 0; s < 2; s++) {
                        // iterate over the 4 smoothers (0:c/b, 1:c/w, 2:r/b, 3:r/w)
                        for (c = 0; c < 2; c++) {
                            int smoother = s * 2 + c;
                            TIC;
                            int sm = (smoother < 2) ? 1 - smoother : 5 - smoother;
                            v_level[l]->multigrid_smoothing(
                                smoother, v, f_Asc_u[smoother], v_level[l]->nblocks[smoother], c,
                                v_level[l]->dep_Asc[smoother], v_level[l]->dep_Asc[sm], v_level[l]->dep_Asc[1],
                                v_level[l]->dep_Asc_ortho[smoother]);

                            if (gyro::icntl[Param::verbose] > 2)
                                std::cout << "SMOOTHER: " << smoother << " = " << TOC << "\n";
                        }
                    }
                } // omp single
            } // omp parallel
        }
    }
    //std::cout << "post-smoothing done \n";
    t = t_smoothing_tmp;
    t_smoothing += TOC;
    TIC;

    t = t_total_tmp;
    t_total += TOC;

    if (gyro::icntl[Param::verbose] > 3)
        gyro::disp(v_level[l]->u, "u");
} /* ----- end of gmgpolar::multigrid_cycle_extrapol ----- */

/*!
 *  \brief Computes the residual
 *
 *  Computes the residual based on the approximation res = b-Au
 * 
 * \param l: the level (0=finest)
 * \param extrapol: l==0 and gyro::icntl[Param::extrapolation]
 *
 */
void gmgpolar::compute_residual(int l, int extrapol)
{
    //std::cout << "compute residual function" << std::endl;

    //compute Au = A*u
    std::vector<double> Au(v_level[l]->m, 0);
    if (gyro::icntl[Param::matrix_free] == 1) {
        if (gyro::icntl[Param::optimized] == 0)
            v_level[l]->apply_A0(v_level[l]->u, Au);
        else {
            double start = omp_get_wtime();
            v_level[l]->apply_A(v_level[l]->u, Au);

            double end = omp_get_wtime();
            t_applyA += (end - start);
        }
    }
    else {
        if (gyro::icntl[Param::openmp] == 1) {
            for (size_t i = 0; i < v_level[l]->row_indices.size(); i++) {
                Au[v_level[l]->row_indices[i]] += v_level[l]->vals[i] * v_level[l]->u[v_level[l]->col_indices[i]];
            }
        }
        else {
            int j;
#pragma omp parallel for private(j) shared(Au)
            for (j = 0; j < v_level[l]->nr; j++) {
                std::vector<int> ptr_vect = v_level[l]->get_ptr(j);
                int ptr_start             = ptr_vect[1];
                int ptr_end               = ptr_start + (ptr_vect[2] - ptr_vect[1]) * v_level[l]->ntheta_int;
                for (int k = ptr_start; k < ptr_end; k++) {
                    Au[v_level[l]->row_indices[k]] += v_level[l]->vals[k] * v_level[l]->u[v_level[l]->col_indices[k]];
                }
            }
        }
    }

    //set res(l) = f(l)
    v_level[l]->res = v_level[l]->fVec;

    //Extrapolation (extrapol = 1)
    if (extrapol == 1) {
        //apply P to f(l+1): Pf = prolong_inj * f(l+1)
        std::vector<double> Pf;
        if (gyro::icntl[Param::matrix_free] == 1) {
            if (gyro::icntl[Param::optimized] == 0) {
                Pf = v_level[l]->apply_prolongation_inj0(v_level[l + 1]->fVec_initial, v_level[l]->mc, v_level[l]->m,
                                                         v_level[l]->coarse_nodes_list_r,
                                                         v_level[l]->coarse_nodes_list_theta, 0);
            }
            else {
                Pf = v_level[l]->apply_prolongation_inj(v_level[l + 1]->fVec_initial);
            }
        }
        else {
            Pf = std::vector<double>(v_level[l]->m, 0);
            for (size_t i = 0; i < v_level[l]->ri_prol_inj.size(); i++) {
                Pf[v_level[l]->ri_prol_inj[i]] +=
                    v_level[l]->v_prol_inj[i] * v_level[l + 1]->fVec_initial[v_level[l]->ci_prol_inj[i]];
            }
        }

        //apply P^T to u (restriction)
        std::vector<double> Pu;
        if (gyro::icntl[Param::matrix_free] == 1) {
            if (gyro::icntl[Param::optimized] == 0) {
                Pu = v_level[l]->apply_prolongation_inj0(v_level[l]->u, v_level[l]->mc, v_level[l]->m,
                                                         v_level[l]->coarse_nodes_list_r,
                                                         v_level[l]->coarse_nodes_list_theta, 1);
            }
            else {
                Pu = v_level[l]->apply_restriction_inj(v_level[l]->u);
            }
        }
        else {
            Pu = std::vector<double>(v_level[l]->mc, 0);
            for (size_t i = 0; i < v_level[l]->ri_prol_inj.size(); i++) {
                Pu[v_level[l]->ci_prol_inj[i]] += v_level[l]->v_prol_inj[i] * v_level[l]->u[v_level[l]->ri_prol_inj[i]];
            }
        }

        //apply A(l+1) to Pu
        std::vector<double> APu(v_level[l + 1]->m, 0);
        if (gyro::icntl[Param::matrix_free] == 1) {
            if (gyro::icntl[Param::optimized] == 0)
                v_level[l + 1]->apply_A0(Pu, APu); //APu = A(l+1) * Pu
            else
                v_level[l + 1]->apply_A(Pu, APu); //APu = A(l+1) * Pu
        }
        else {
            if (gyro::icntl[Param::openmp] == 1) {
                for (size_t i = 0; i < v_level[l + 1]->row_indices.size(); i++) {
                    APu[v_level[l + 1]->row_indices[i]] += v_level[l + 1]->vals[i] * Pu[v_level[l + 1]->col_indices[i]];
                }
            }
            else {
                int j;
#pragma omp parallel for private(j) shared(APu)
                for (j = 0; j < v_level[l + 1]->nr; j++) {
                    std::vector<int> ptr_vect = v_level[l + 1]->get_ptr(j);
                    int ptr_start             = ptr_vect[1];
                    int ptr_end               = ptr_start + (ptr_vect[2] - ptr_vect[1]) * v_level[l + 1]->ntheta_int;
                    for (int k = ptr_start; k < ptr_end; k++) {
                        APu[v_level[l + 1]->row_indices[k]] +=
                            v_level[l + 1]->vals[k] * Pu[v_level[l + 1]->col_indices[k]];
                    }
                }
            }
        }

        //apply P to APu (prolongation)
        std::vector<double> PAPu;
        if (gyro::icntl[Param::matrix_free] == 1) {
            if (gyro::icntl[Param::optimized] == 0) {
                PAPu = v_level[l]->apply_prolongation_inj0(APu, v_level[l]->mc, v_level[l]->m,
                                                           v_level[l]->coarse_nodes_list_r,
                                                           v_level[l]->coarse_nodes_list_theta, 0);
            }
            else {
                PAPu = v_level[l]->apply_prolongation_inj(APu);
            }
        }
        else {
            PAPu = std::vector<double>(v_level[l]->m, 0);
            for (size_t i = 0; i < v_level[l]->ri_prol_inj.size(); i++) {
                PAPu[v_level[l]->ri_prol_inj[i]] += v_level[l]->v_prol_inj[i] * APu[v_level[l]->ci_prol_inj[i]];
            }
        }

        //res_ex = 4/3 * f(l) - 1/3 * f_1   -   4/3 * Au - 1/3 *PAPu
        for (int i = 0; i < v_level[l]->m; ++i) {
            v_level[l]->res[i] =
                (4. / 3. * v_level[l]->fVec[i] - 1. / 3. * Pf[i]) - (4. / 3. * Au[i] - 1. / 3. * PAPu[i]);
        }
    }
    //no extrapolation
    else {
        //compute res = f - A * u
        for (unsigned long int i = 0; i < v_level[l]->res.size(); ++i) {
            v_level[l]->res[i] -= Au[i]; //res(l) = res(l) - Au = f(l) - A*u
        }
    }
} /* ----- end of gmgpolar::compute_residual ----- */

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
    std::vector<double> error;

    int ntheta_int = v_level[0]->ntheta;
    if (fabs(v_level[0]->theta[v_level[0]->ntheta - 1] - 2 * PI) < 1e-10) {
        ntheta_int--;
    }

    if (gyro::f_sol_in.empty()) {
        for (int j = 0; j < v_level[0]->nr; ++j) { //iterate over the grid in polar coordinates
            for (int i = 0; i < ntheta_int; ++i) {
                int index = j * ntheta_int + i; //corresponding index of the solution vector
                //double x;
                //double y;
                double theoretical_solution =
                    gyro::def_solution_rt(v_level[0]->r[j], v_level[0]->theta[i], 0); //compute the theoretical solution
                double computed_solution = v_level[0]->u[index];
                error.push_back(fabs(theoretical_solution - computed_solution));
            }
        }
    }
    else {
        for (int j = 0; j < v_level[0]->nr; ++j) { //iterate over the grid in polar coordinates
            for (int i = 0; i < ntheta_int; ++i) {
                int index                   = j * ntheta_int + i; //corresponding index of the solution vector
                double theoretical_solution = v_level[0]->sol_in[index];
                double computed_solution    = v_level[0]->u[index];
                error.push_back(fabs(theoretical_solution - computed_solution));
            }
        }
    }

    return error;
} /* ----- end of gmgpolar::compute_error ----- */

/*!
 *  \brief Compute the backward error on level 0
 *
 *  Compute the backward error on level 0:
 *  ||b-Au||_inf / (||A||_inf ||x||_1 + ||b||_inf)
 * 
 * \return the backward error
 *
 */
double gmgpolar::compute_backwarderror()
{
    gmgpolar::compute_residual(0, 0); //compute residual on level 0, without extrapolation, save in v_level[0]->res

    double nrm_inf_res = 0; //inf norm of res=A*u -fVec
    for (std::size_t i = 0; i < v_level[0]->res.size(); ++i) { //iterate over the grid in polar coordinates
        if (fabs(v_level[0]->res[i]) > nrm_inf_res) {
            nrm_inf_res = fabs(v_level[0]->res[i]);
        }
    }

    double nrm_inf_rhs = 0; //inf norm of the rhs-vector
    for (std::size_t i = 0; i < v_level[0]->fVec.size(); ++i) { //iterate over the grid in polar coordinates
        if (fabs(v_level[0]->fVec[i]) > nrm_inf_rhs) {
            nrm_inf_rhs = fabs(v_level[0]->fVec[i]);
        }
    }

    double nrm_1_u = 0; //1-norm of the solution vector
    for (std::size_t i = 0; i < v_level[0]->u.size(); ++i) { //iterate over the grid in polar coordinates
        nrm_inf_rhs += fabs(v_level[0]->u[i]);
    }

    double nrm_inf_A = 0; //inf norm of the operator A
    double row_sum   = 0;
    int row_index    = 0;
    for (std::size_t i = 0; i < v_level[0]->vals.size(); ++i) { //iterate over all elements in A
        if (v_level[0]->row_indices[i] == row_index) { //we are still in the same row of the matrix
            row_sum += fabs(v_level[0]->vals[i]);
        }
        else { //we just started a new row of the matrix
            if (nrm_inf_A < row_sum) {
                nrm_inf_A = row_sum;
            }
            row_sum   = 0; //set row_sum to zero again
            row_index = v_level[0]->row_indices[i]; //save the last row_index to compare
        }
    }

    double backward_error = nrm_inf_res / (nrm_inf_A * nrm_1_u + nrm_inf_rhs);

    return backward_error;
} /* ----- end of gmgpolar::compute_backwarderror ----- */
