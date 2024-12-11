#pragma once

class LevelCache;
class Level;

#include <chrono>
#include <iostream>
#include <vector>
#include <cuda_runtime.h>
#include <cusparse.h>
#include <cusparse_v2.h>
#include <cusolverSp.h>


#include "../../PolarGrid/polargrid.h"
#include "../../InputFunctions/boundaryConditions.h"
#include "../../InputFunctions/densityProfileCoefficients.h"
#include "../../InputFunctions/domainGeometry.h"
#include "../../InputFunctions/sourceTerm.h"
#include "../../Level/level.h"
#include "../../LinearAlgebra/Vector/gpu_vector.h"
#include "../../common/constants.h"

#include "../../LinearAlgebra/Matrix/matrix.h"
#include "dmumps_c.h"
#include "mpi.h"

class ExtrapolatedSmootherTakeGPU
{
public:
    explicit ExtrapolatedSmootherTakeGPU(const Level& level, const DomainGeometry& domain_geometry,
                      const DensityProfileCoefficients& density_profile_coefficients, const bool DirBC_Interior);
    ~ExtrapolatedSmootherTakeGPU();

    void extrapolatedSmoothingInPlace(GPU_Vector<double>& x, const GPU_Vector<double>& rhs, GPU_Vector<double>& temp);

private:
    /* ------------------- */
    /* Constructor members */
    const Level& level_;
    const DomainGeometry& domain_geometry_;
    const DensityProfileCoefficients& density_profile_coefficients_;
    const bool DirBC_Interior_;

    /* Cuda Solvers */
    cusparseHandle_t sparse_handle_;

    /* Circle Tridiagonal Matrices */
    double* circle_main_diagonals_;
    double* circle_lower_diagonals_;
    double* circle_upper_diagonals_;
    double* sherman_morrison_gammas_;
    double* factor_;
    /* Radial Tridiagonal Matrices */
    double* radial_main_diagonals_;
    double* radial_lower_diagonals_;
    double* radial_upper_diagonals_;
    /* Tridiagonal Solver Buffer */
    void* pBuffer_;
    /* Inner Boundary Mumps COO Matrix */
    std::unique_ptr<int[]> inner_boundary_matrix_row_indices_;
    std::unique_ptr<int[]> inner_boundary_matrix_column_indices_;
    std::unique_ptr<double[]> inner_boundary_matrix_values_;
    DMUMPS_STRUC_C inner_boundary_mumps_solver_;
    int inner_boundary_matrix_nnz_;
    int* d_inner_boundary_matrix_row_indices_;
    int* d_inner_boundary_matrix_column_indices_;
    double* d_inner_boundary_matrix_values_;

    /* Build Smoother Matrices we have allocated. */
    void buildAscMatrices();

    void initializeMumps();
    void finalizeMumpsSolver();

    /* The cyclic tridiagonal Matrices need to be adjusted to a system of a non-cyclic tridiagonal matrices. */
    void adjustAscCircle_ShermanMorrison();

    /* temp is needed to prepare for the Sherman-Morrison formula. */
    void applyAscOrtho_BlackCircle(
        GPU_Vector<double>& x, const GPU_Vector<double>& rhs, GPU_Vector<double>& temp, 
        DomainGeometry* device_domain_geometry);
    void applyAscOrtho_WhiteCircle(
        GPU_Vector<double>& x, const GPU_Vector<double>& rhs, GPU_Vector<double>& temp, 
        DomainGeometry* device_domain_geometry);

    void applyAscOrtho_BlackRadial(
        GPU_Vector<double>& x, const GPU_Vector<double>& rhs,
        DomainGeometry* device_domain_geometry);
    void applyAscOrtho_WhiteRadial(
        GPU_Vector<double>& x, const GPU_Vector<double>& rhs,
        DomainGeometry* device_domain_geometry);


    void solveAsc_BlackCircle(GPU_Vector<double>& x, const GPU_Vector<double>& rhs, GPU_Vector<double>& temp);
    void solveAsc_WhiteCircle(GPU_Vector<double>& x, const GPU_Vector<double>& rhs, GPU_Vector<double>& temp);
    void solveCircleDiagonals(GPU_Vector<double>& x, const GPU_Vector<double>& rhs, GPU_Vector<double>& temp);
    void solveCircleTridiagonals(GPU_Vector<double>& x, const GPU_Vector<double>& rhs, GPU_Vector<double>& temp);

    void solveAsc_BlackRadial(GPU_Vector<double>& x, const GPU_Vector<double>& rhs);
    void solveAsc_WhiteRadial(GPU_Vector<double>& x, const GPU_Vector<double>& rhs);
};