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
 * \file gmgpolar.cpp
 * \brief Implementation of the class gmgpolar
 * \author M. Kuehn, C. Kruse, P. Leleux, C. Schwarz
 * \version 0.0
 */
#include "gmgpolar.h"

/*!
 *  \brief Default Constructor of gmgpolar class
 *
 *  Default constructor of the gmgpolar class, sets default values for
 *  attributes.
 *
 */
gmgpolar::gmgpolar()
{
    reset_timers();
} /* ----- end of constructor gmgpolar:gmgpolar ----- */

/*!
 *  \brief Default Destructor of gmgpolar class
 *
 *  Default destructor of the gmgpolar class.
 *
 */
gmgpolar::~gmgpolar()
{
    for (int i = 0; i < levels_orig; i++) {
        delete v_level[i];
    }
} /* ----- end of destructor gmgpolar::~gmgpolar ----- */

/*!
 *  \brief Sets execution times to 0
 *
 *  Sets execution times to 0.
 *
 */
void gmgpolar::reset_timers()
{
    t_setup         = 0;
    t_build         = 0;
    t_facto_Ac      = 0;
    t_build_P       = 0;
    t_build_Asc     = 0;
    t_facto_Asc     = 0;
    t_total         = 0;
    t_smoothing     = 0;
    t_residual      = 0;
    t_restriction   = 0;
    t_Ac            = 0;
    t_prolongation  = 0;
    t_applyA        = 0;
    t_fine_residual = 0;
    t_error         = 0;
} /* ----- end of gmgpolar::reset_timers ----- */
