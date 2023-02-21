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

/**
 * \file main.cpp
 * \brief Test program to launch the application from:
 * "Implicitly extrapolated geometric multigrid on disk-like domains for
 the gyrokinetic Poisson equation from fusion plasma applications", Martin
 Joachim Kühn, Carola Kruse, Ulrich Rüde
 * \author M. Kuehn, C. Kruse, P. Leleux, C. Schwarz
 * \version 0.0
 * \date March 31st 2021
 *
 * This program launch the application with the problem and parameters specified in a
 * configuration file (first argument or "config_file.info" by default)
 *
 */
#include <omp.h>
#include <sched.h>
#include <time.h>
#include "cmdline.h"
#include "gmgpolar.h"

#define ICPU sched_getcpu()
#define NOMP omp_get_num_threads()
#define IOMP omp_get_thread_num()

/**
 * \fn int main ()
 * \brief Main function
 *
 * \param 1: configuration file (Optional)
 * \return EXIT_SUCCESS - Arrêt normal du programme.
 */
int main(int argc, char* argv[])
{
    int error = 0;

    ////////////////////////////////////////////////////////////////////////////////
    //    READING PARAMETERS
    ////////////////////////////////////////////////////////////////////////////////
    cmdline::parser a;
    a.add<int>("optimized", '\0', "", false, 1);
    a.add<int>("matrix_free", '\0', "", false, 1);
    a.add<int>("debug", 'D', "", false, 0);
    a.add<int>("nr_exp", 'n', "", false, 4);
    a.add<int>("ntheta_exp", '\0', "", false, 4);
    a.add<int>("fac_ani", 'a', "", false, 3);
    a.add<int>("v1", '\0', "", false, 1);
    a.add<int>("v2", '\0', "", false, 1);
    a.add<int>("cycle", 'c', "", false, 1);
    a.add<int>("mod_pk", '\0', "", false, 0);
    a.add<int>("compute_rho", '\0', "", false, 0);
    a.add<int>("level", 'l', "", false, -1);
    a.add<int>("maxiter", '\0', "", false, 150);
    a.add<int>("theta_aniso", '\0', "", false, 0);
    a.add<int>("smoother", '\0', "", false, 3);
    a.add<int>("extrapolation", 'E', "", false, 0);
    a.add<int>("DirBC_Interior", '\0', "", false, 1);
    a.add<int>("divideBy2", '\0', "", false, 0);
    a.add<int>("prob", '\0', "", false, 5);
    a.add<int>("alpha_coeff", '\0', "", false, 0);
    a.add<int>("beta_coeff", '\0', "", false, 0);
    a.add<int>("verbose", '\0', "", false, 1);
    a.add<int>("openmp", '\0', "", false, 1);
    a.add<int>("res_norm", '\0', "", false, 3);
    a.add<int>("write_radii_angles", '\0', "", false, 0);
    a.add<int>("check_error", '\0', "", false, 1);

    a.add<float>("R0", 'r', "", false, 1e-5);
    a.add<float>("R", 'R', "", false, 1.3);
    a.add<float>("kappa_eps", 'k', "", false, 42);
    a.add<float>("delta_e", 'd', "", false, 42);
    a.add<float>("tol_bound_check", 'e', "", false, 1e-8);
    a.add<float>("rel_red_conv", '\0', "", false, 1e-8);

    a.add<std::string>("f_grid_r", '\0', "", false, "");
    a.add<std::string>("f_grid_theta", '\0', "", false, "");
    a.add<std::string>("f_sol_in", '\0', "", false, "");
    a.add<std::string>("f_sol_out", '\0', "", false, "");

    a.parse_check(argc, argv);

    //std::cout << "Initializing parameters...\n";
    gyro::init_params();
    gyro::icntl[Param::verbose]     = a.get<int>("verbose");
    gyro::icntl[Param::openmp]      = a.get<int>("openmp");
    gyro::icntl[Param::optimized]   = a.get<int>("optimized");
    gyro::icntl[Param::matrix_free] = a.get<int>("matrix_free");
    gyro::icntl[Param::debug]       = a.get<int>("debug");
    // gyro::icntl[Param::divideBy2]    = a.get<int>("divideBy2");
    gyro::icntl[Param::smoother]     = a.get<int>("smoother");
    gyro::dcntl[Param::rel_red_conv] = a.get<float>("rel_red_conv");

    gyro::f_grid_r     = a.get<std::string>("f_grid_r");
    gyro::f_grid_theta = a.get<std::string>("f_grid_theta");
    if (!gyro::f_grid_r.empty() || !gyro::f_grid_theta.empty()) {
        std::cout << "File for grid in r: " << gyro::f_grid_r << ", and theta: " << gyro::f_grid_theta << "\n";
    }
    gyro::f_sol_in  = a.get<std::string>("f_sol_in");
    gyro::f_sol_out = a.get<std::string>("f_sol_out");
    if (!gyro::f_sol_in.empty() && gyro::icntl[Param::check_error]) {
        std::cout << "Warning: an input solution has been provided but the error will not be checked\n";
    }
    else if (!gyro::f_sol_in.empty()) {
        std::cout << "File to read the solution vector: " << gyro::f_sol_in << "\n";
    }
    if (!gyro::f_sol_out.empty()) {
        std::cout << "File to write the solution vector: " << gyro::f_sol_in << "\n";
    }

    ////////////////////////////////////////////////////////////////////////////////
    //    DISPLAY OPENMP INFO
    ////////////////////////////////////////////////////////////////////////////////
    //Set number of threads for the openmp parallelization
    omp_set_num_threads(gyro::icntl[Param::openmp]);

    if (gyro::icntl[Param::verbose] > 0) {
#pragma omp parallel
        {
#pragma omp master
            {
                std::cout << "Number of OpenMP threads: " << NOMP << "\n";
            }
        }
    }
    if (gyro::icntl[Param::verbose] > 0) {

#pragma omp parallel
        {
#pragma omp master
            {
                std::cout << "OMP_thread\tCPU\n";
            }
#pragma omp critical
            {
                std::cout << IOMP << "\t" << ICPU << "\n";
            }
        }
        std::cout << "\n";
    }

    // normal run, NO DEBUGGING
    if (gyro::icntl[Param::debug] == 0) {
        gyro::icntl[Param::nr_exp]     = a.get<int>("nr_exp");
        gyro::icntl[Param::ntheta_exp] = a.get<int>("ntheta_exp");
        gyro::icntl[Param::fac_ani]    = a.get<int>("fac_ani");
        gyro::icntl[Param::v1]         = a.get<int>("v1");
        gyro::icntl[Param::v2]         = a.get<int>("v2");
        gyro::icntl[Param::cycle]      = a.get<int>("cycle");
        gyro::icntl[Param::mod_pk]     = a.get<int>("mod_pk");
        if (gyro::icntl[Param::mod_pk] == 42)
            gyro::icntl[Param::mod_pk] = 0;
        gyro::icntl[Param::compute_rho]        = a.get<int>("compute_rho");
        gyro::icntl[Param::level]              = a.get<int>("level");
        gyro::icntl[Param::maxiter]            = a.get<int>("maxiter");
        gyro::icntl[Param::theta_aniso]        = a.get<int>("theta_aniso");
        gyro::icntl[Param::extrapolation]      = a.get<int>("extrapolation");
        gyro::icntl[Param::DirBC_Interior]     = a.get<int>("DirBC_Interior");
        gyro::icntl[Param::divideBy2]          = a.get<int>("divideBy2");
        gyro::icntl[Param::prob]               = a.get<int>("prob");
        gyro::icntl[Param::alpha_coeff]        = a.get<int>("alpha_coeff");
        gyro::icntl[Param::beta_coeff]         = a.get<int>("beta_coeff");
        gyro::icntl[Param::res_norm]           = a.get<int>("res_norm");
        gyro::icntl[Param::write_radii_angles] = a.get<int>("write_radii_angles");
        gyro::icntl[Param::check_error]        = a.get<int>("check_error");
        if (gyro::icntl[Param::extrapolation] >= 2 && gyro::icntl[Param::check_error] == 0) {
            throw std::runtime_error("The alternative extrapolation technique requires to check the error, w.r.t. to "
                                     "the theoretical solution, as stopping criterion.");
        }

        gyro::dcntl[Param::R0]        = a.get<float>("R0");
        gyro::dcntl[Param::R]         = a.get<float>("R");
        gyro::dcntl[Param::kappa_eps] = a.get<float>("kappa_eps");
        gyro::dcntl[Param::delta_e]   = a.get<float>("delta_e");

        geometry_type geom = (geometry_type)gyro::icntl[Param::mod_pk];
        if (gyro::dcntl[Param::kappa_eps] == 42 && gyro::dcntl[Param::delta_e] == 42) {
            gyro::get_geometry_coeffs(geom);
        }
        if (gyro::icntl[Param::verbose] > 1)
            gyro::show_params();
        gyro::select_functions_class(gyro::icntl[Param::alpha_coeff], gyro::icntl[Param::beta_coeff],
                                     gyro::icntl[Param::mod_pk], gyro::icntl[Param::prob]);

        ////////////////////////////////////////////////////////////////////////////////
        //    LAUNCH
        ////////////////////////////////////////////////////////////////////////////////
        double t;
        TIC;
        gmgpolar gmg;
        try {
            // Create polar grids
            gmg.create_grid_polar();

            if (gyro::icntl[Param::write_radii_angles] == 0)
                // Solve using multigrid
                gmg.polar_multigrid();
        }
        catch (std::runtime_error const& e) {
            std::cout << "Error code : " << e.what() << "\n";
        }
        catch (...) {
            std::cout
                << "I felt a great disturbance in the Force, as if millions of voices suddenly cried out in terror "
                   "and were suddenly silenced. I fear something terrible has happened.\n";
            error = 1;
        }
        // if (gyro::icntl[Param::verbose] > 0)
        std::cout << "Total_execution_time: " << TOC
                  << "\n---------------------------------------------------------------------------\n";
    }

    // DEBUGGING
    else {
        geometry_type geom;
        clock_t t;
        ////////////////////////////////////////////////////////////////////////////////
        //    LOOP PARAMETERS
        ////////////////////////////////////////////////////////////////////////////////
        // debug=1
        // nr_exp=4
        // ntheta_exp=4
        // fac_ani=3
        // divideBy2=0
        gyro::icntl[Param::nr_exp]     = 4;
        gyro::icntl[Param::ntheta_exp] = 4;
        gyro::icntl[Param::fac_ani]    = 3;
        gyro::icntl[Param::divideBy2]  = 0;
        gyro::icntl[Param::verbose]    = 0;
        std::cout << "####################################################\n";
        std::cout << "GLOBAL PARAMETERS // nr_exp" << gyro::icntl[Param::nr_exp]
                  << ", ntheta_exp: " << gyro::icntl[Param::ntheta_exp] << ", fac_ani: " << gyro::icntl[Param::fac_ani]
                  << ", divideBy2: " << gyro::icntl[Param::divideBy2] << ", verbose: " << gyro::icntl[Param::verbose]
                  << "\n";
        std::cout << "####################################################\n";

        // prob5+R1+mod_pk2+alpha1+beta1+optimized1+DirBC_Interior1+smoother3+extrapolation1: 1
        gyro::icntl[Param::prob]           = 5;
        gyro::dcntl[Param::R]              = 1.0;
        gyro::icntl[Param::mod_pk]         = 2;
        gyro::icntl[Param::alpha_coeff]    = 1;
        gyro::icntl[Param::beta_coeff]     = 1;
        gyro::icntl[Param::optimized]      = 1;
        gyro::icntl[Param::DirBC_Interior] = 1;
        gyro::icntl[Param::smoother]       = 3;
        gyro::icntl[Param::extrapolation]  = 1;
        std::cout << "#################################################################################################"
                     "#######\n";
        std::cout << "TESTING DEFAULT CASE // prob" << gyro::icntl[Param::prob] << ", R: " << gyro::dcntl[Param::R]
                  << ", mod_pk: " << gyro::icntl[Param::mod_pk] << ", alpha_coeff: " << gyro::icntl[Param::alpha_coeff]
                  << ", beta_coeff: " << gyro::icntl[Param::beta_coeff] << ", optimized"
                  << gyro::icntl[Param::optimized] << ", DirBC_Interior: " << gyro::icntl[Param::DirBC_Interior]
                  << ", smoother: " << gyro::icntl[Param::smoother]
                  << ", extrapolation: " << gyro::icntl[Param::extrapolation] << "\n";
        std::cout << "#################################################################################################"
                     "#######\n";
        geom = (geometry_type)gyro::icntl[Param::mod_pk];
        gyro::get_geometry_coeffs(geom);
        if (gyro::icntl[Param::verbose] > 1)
            gyro::show_params();
        gyro::select_functions_class(gyro::icntl[Param::alpha_coeff], gyro::icntl[Param::beta_coeff],
                                     gyro::icntl[Param::mod_pk], gyro::icntl[Param::prob]);
        ////////////////////////////////////////////////////////////////////////////////
        //    LAUNCH
        ////////////////////////////////////////////////////////////////////////////////
        TIC;
        gmgpolar gmg;
        try {
            // Create polar grids
            gmg.create_grid_polar();

            // Solve using multigrid
            gmg.polar_multigrid();
        }
        catch (std::runtime_error const& e) {
            std::cout << "Error code : " << e.what() << "\n";
            exit(1);
        }
        catch (...) {
            std::cout << "I felt a great disturbance in the Force, as if millions of voices suddenly "
                         "cried out "
                         "in terror "
                         "and were suddenly silenced. I fear something terrible has happened.\n";
            error = 1;
        }
        if (gyro::icntl[Param::verbose] > 0)
            std::cout << "Total_execution_time: " << TOC << "\n";

        std::cout << "\n\n\n\n\n";

        // prob5+R1+mod_pk2+alpha1+beta1: 6
        // optimized=0-1
        // DirBC_Interior=0-1
        // smoother=3,13
        // extrapolation=0-1
        gyro::icntl[Param::prob]        = 5;
        gyro::dcntl[Param::R]           = 1.0;
        gyro::icntl[Param::mod_pk]      = 2;
        gyro::icntl[Param::alpha_coeff] = 1;
        gyro::icntl[Param::beta_coeff]  = 1;
        std::cout << "#################################################################################################"
                     "#######\n";
        std::cout << "TESTING MG PARAMETERS // prob" << gyro::icntl[Param::prob] << ", R: " << gyro::dcntl[Param::R]
                  << ", mod_pk: " << gyro::icntl[Param::mod_pk] << ", alpha_coeff: " << gyro::icntl[Param::alpha_coeff]
                  << ", beta_coeff: " << gyro::icntl[Param::beta_coeff] << "\n";
        std::cout << "#################################################################################################"
                     "#######\n";
        // for (int optimized = 0; optimized < 2; optimized++) {
        for (int DirBC_Interior = 0; DirBC_Interior < 2; DirBC_Interior++) {
            // for (int smoother = 3; smoother < 23; smoother += 10) {
            for (int extrapolation = 0; extrapolation < 2; extrapolation++) {
                // // int optimized                      = 1;
                // int DirBC_Interior = 0;
                // // int smoother                       = 3;
                // int extrapolation = 1;
                // gyro::icntl[Param::optimized]      = optimized;
                gyro::icntl[Param::DirBC_Interior] = DirBC_Interior;
                // gyro::icntl[Param::smoother]       = smoother;
                gyro::icntl[Param::extrapolation] = extrapolation;

                geom = (geometry_type)gyro::icntl[Param::mod_pk];
                gyro::get_geometry_coeffs(geom);
                if (gyro::icntl[Param::verbose] > 1)
                    gyro::show_params();
                gyro::select_functions_class(gyro::icntl[Param::alpha_coeff], gyro::icntl[Param::beta_coeff],
                                             gyro::icntl[Param::mod_pk], gyro::icntl[Param::prob]);
                std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
                // std::cout << "opti: " << optimized << ", BC: " << DirBC_Interior << ", smoother: " << smoother
                //           << ", extrap: " << extrapolation << "\n";
                std::cout << "BC: " << DirBC_Interior << ", extrap: " << extrapolation << "\n";
                std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";

                ////////////////////////////////////////////////////////////////////////////////
                //    LAUNCH
                ////////////////////////////////////////////////////////////////////////////////
                TIC;
                gmgpolar gmg2;
                try {
                    // Create polar grids
                    gmg2.create_grid_polar();

                    // Solve using multigrid
                    gmg2.polar_multigrid();
                }
                catch (std::runtime_error const& e) {
                    std::cout << "Error code : " << e.what() << "\n";
                    exit(1);
                }
                catch (...) {
                    std::cout << "I felt a great disturbance in the Force, as if millions of voices suddenly "
                                 "cried out "
                                 "in terror "
                                 "and were suddenly silenced. I fear something terrible has happened.\n";
                    error = 1;
                }
                if (gyro::icntl[Param::verbose] > 0)
                    std::cout << "Total_execution_time: " << TOC << "\n";
            }
            // }
        }
        // }

        std::cout << "\n\n\n\n\n";

        // optimized1+DirBC_Interior1+smoother3+extrapolation1: 32
        // R=1.0,1.3
        // prob=5-6
        // mod_pk=0-2(3)
        // alpha_coeff=0-1
        // beta_coeff=0-1
        gyro::icntl[Param::optimized]      = 1;
        gyro::icntl[Param::DirBC_Interior] = 1;
        gyro::icntl[Param::smoother]       = 3;
        gyro::icntl[Param::extrapolation]  = 1;
        std::cout << "#################################################################################################"
                     "#######\n";
        std::cout << "TESTING TEST CASES // optimized" << gyro::icntl[Param::optimized]
                  << ", DirBC_Interior: " << gyro::icntl[Param::DirBC_Interior]
                  << ", smoother: " << gyro::icntl[Param::smoother]
                  << ", extrapolation: " << gyro::icntl[Param::extrapolation] << "\n";
        std::cout << "#################################################################################################"
                     "#######\n";
        for (double R = 1.0; R < 1.6; R += 0.3) {
            for (int prob = 5; prob < 7; prob++) {
                for (int mod_pk = 0; mod_pk < 3; mod_pk++) {
                    for (int alpha_coeff = 0; alpha_coeff < 2; alpha_coeff++) {
                        for (int beta_coeff = 0; beta_coeff < 2; beta_coeff++) {
                            // double R                        = 1.0;
                            // int prob                        = 5;
                            // int mod_pk                      = 1;
                            // int alpha_coeff                 = 1;
                            // int beta_coeff                  = 1;
                            gyro::dcntl[Param::R]           = R;
                            gyro::icntl[Param::prob]        = prob;
                            gyro::icntl[Param::mod_pk]      = mod_pk;
                            gyro::icntl[Param::alpha_coeff] = alpha_coeff;
                            gyro::icntl[Param::beta_coeff]  = beta_coeff;

                            geom = (geometry_type)gyro::icntl[Param::mod_pk];
                            gyro::get_geometry_coeffs(geom);
                            if (gyro::icntl[Param::verbose] > 1)
                                gyro::show_params();
                            gyro::select_functions_class(gyro::icntl[Param::alpha_coeff],
                                                         gyro::icntl[Param::beta_coeff], gyro::icntl[Param::mod_pk],
                                                         gyro::icntl[Param::prob]);
                            std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
                            std::cout << "R: " << R << ", prob: " << prob << ", mod_pk: " << mod_pk
                                      << ", alpha_coeff: " << alpha_coeff << ", beta_coeff: " << beta_coeff << "\n";
                            std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";

                            ////////////////////////////////////////////////////////////////////////////////
                            //    LAUNCH
                            ////////////////////////////////////////////////////////////////////////////////
                            TIC;
                            gmgpolar gmg3;
                            try {
                                // Create polar grids
                                gmg3.create_grid_polar();

                                // Solve using multigrid
                                gmg3.polar_multigrid();
                            }
                            catch (std::runtime_error const& e) {
                                std::cout << "Error code : " << e.what() << "\n";
                                exit(1);
                            }
                            catch (...) {
                                std::cout
                                    << "I felt a great disturbance in the Force, as if millions of voices suddenly "
                                       "cried out "
                                       "in terror "
                                       "and were suddenly silenced. I fear something terrible has happened.\n";
                                error = 1;
                            }
                            if (gyro::icntl[Param::verbose] > 0)
                                std::cout << "Total_execution_time: " << TOC << "\n";
                        }
                    }
                }
            }
        }
    }
    return error;
} /* ----- end of main ----- */
