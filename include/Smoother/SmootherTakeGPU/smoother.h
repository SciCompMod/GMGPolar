#pragma once

class LevelCache;
class Level;

#include <chrono>
#include <iostream>
#include <vector>
#include <cuda_runtime.h>
#include <cusparse_v2.h>

#include "../../PolarGrid/polargrid.h"

#include "../../InputFunctions/boundaryConditions.h"
#include "../../InputFunctions/densityProfileCoefficients.h"
#include "../../InputFunctions/domainGeometry.h"
#include "../../InputFunctions/sourceTerm.h"
#include "../../Level/level.h"
#include "../../LinearAlgebra/Vector/gpu_vector.h"
#include "../../LinearAlgebra/Solvers/gpu_symmetric_tridiagonal_solver.h"
#include "../../common/constants.h"



class SmootherTakeGPU
{
public:
    explicit SmootherTakeGPU(const Level& level, const DomainGeometry& domain_geometry,
                      const DensityProfileCoefficients& density_profile_coefficients, const bool DirBC_Interior);
    ~SmootherTakeGPU();

    void smoothingInPlace(GPU_Vector<double>& x, const GPU_Vector<double>& rhs, GPU_Vector<double>& temp);

private:
    /* ------------------- */
    /* Constructor members */
    const Level& level_;
    const DomainGeometry& domain_geometry_;
    const DensityProfileCoefficients& density_profile_coefficients_;
    const bool DirBC_Interior_;

    cusparseHandle_t handle_;

    double* circle_main_diagonals_;
    double* circle_lower_diagonals_;
    double* circle_upper_diagonals_;
    double* sherman_morrison_gammas_;

    double* radial_main_diagonals_;
    double* radial_lower_diagonals_;
    double* radial_upper_diagonals_;
    
    void* pBuffer_;

    void buildAscMatrices();
    void shermannMorrisonAdjustment();
};