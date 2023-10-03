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

// The class for the multigrid gmgpolar where each level is from class level

/*!
 * \file gmgpolar.h
 * \brief Header for the class gmgpolar
 * \author M. Kuehn, C. Kruse, P. Leleux
 * \version 0.0
 */
#ifndef GMGPOLAR_HXX_
#define GMGPOLAR_HXX_

#include <stdlib.h>
#include <cstdio>
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <assert.h>
#include <sys/stat.h>
#include <iomanip>
#include <ctime>
#include "config.h"
#include "constants.h"
#include "level.h"
#include "gyro.h"

class gmgpolar
{
public:
/*******************************************************************************
 * Attributes
 ******************************************************************************/
    /***************************************************************************
     * Grid levels
     **************************************************************************/
    /*! Vector of levels for the MG scheme (geometry, operators, etc.) */
    int levels, levels_orig; // Number of levels
    std::vector<level*> v_level;

    /***************************************************************************
     * Multigrid scheme
     **************************************************************************/
    std::vector<double> nrm_2_res;
    std::vector<double> nrm_2_err;
    std::vector<double> nrm_inf_res;

    /* execution times */
    double t_setup; // prepare_op_levels
    double t_build; // build A and RHS
    double t_facto_Ac; // factorization of coarse operator
    double t_build_P; // build coarse nodes, line splitting and P
    double t_build_Asc; // build Asc and Asc_ortho
    double t_facto_Asc; // factorization of Asc
    double t_total; // multigrid_cycle_extrapol
    double t_smoothing; // pre- and post-smoothing (including solve Asc)
    double t_residual; // factorization of Asc
    double t_restriction; // restriction
    double t_Ac; // solve coarse operator
    double t_prolongation; // prolongation and coarse grid correction
    double t_fine_residual; // computation of residual on level 0 (init and after iteration)
    double t_error; // computation of final error
    double t_applyA; // apply the operator A,  method apply_A

    /*******************************************************************************
 * Methods
 ******************************************************************************/
    gmgpolar();
    ~gmgpolar();

    void reset_timers();
    void create_grid_polar();
    void create_grid_polar_divide();
    void polar_multigrid();
    void check_geom();
    void define_coarse_nodes();
    void prepare_op_levels();
    void debug();
    void multigrid_iter();
    void multigrid_cycle_extrapol(int l);
    void compute_residual(int l, int extrapol);
    std::vector<double> compute_error();
    double compute_backwarderror();

private:
    /*******************************************************************************
 * Attributes
 ******************************************************************************/

    /*******************************************************************************
 * Methods
 ******************************************************************************/
};

#endif // GMGPOLAR_HXX
