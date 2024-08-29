#pragma once

class LevelCache;
class Level;

#include <chrono>
#include <vector>
#include <iostream>

#include "mpi.h" 
#include "dmumps_c.h"   

#include "../PolarGrid/polargrid.h"
#include "../InputFunctions/domainGeometry.h"
#include "../InputFunctions/densityProfileCoefficients.h"
#include "../InputFunctions/boundaryConditions.h"
#include "../InputFunctions/sourceTerm.h"
#include "../LinearAlgebra/vector.h"
#include "../LinearAlgebra/matrix.h"
#include "../LinearAlgebra/vector_operations.h"
#include "../LinearAlgebra/symmetricTridiagonalSolver.h"
#include "../common/constants.h"
#include "../Level/level.h"
#include "../Stencil/stencil.h"

class Smoother {
public:
    explicit Smoother(const PolarGrid& grid, const LevelCache& level_cache, 
                      const DomainGeometry& domain_geometry,
                      bool DirBC_Interior, int num_omp_threads);
    ~Smoother();

    void smoothingInPlace(Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp);

    void smoothingInPlaceSequential(Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp);
    void smoothingInPlaceForLoop(Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp); /* This is the fastest option */
    void smoothingInPlaceTaskLoop(Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp);
    void smoothingInPlaceTaskDependencies(Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp);

private:
    /* ------------------- */
    /* Constructor members */
    const PolarGrid& grid_;
    const std::vector<double>& sin_theta_cache_;
    const std::vector<double>& cos_theta_cache_;
    const std::vector<double>& coeff_alpha_cache_;
    const std::vector<double>& coeff_beta_cache_;
    const DomainGeometry& domain_geometry_;
    const bool DirBC_Interior_;
    const int num_omp_threads_;

    /* ---------------- */
    /* Smoother members */
    SparseMatrix<double> inner_boundary_circle_matrix_;
    DMUMPS_STRUC_C inner_boundary_mumps_solver_;
    std::vector<SymmetricTridiagonalSolver<double>> circle_tridiagonal_solver_;
    std::vector<SymmetricTridiagonalSolver<double>> radial_tridiagonal_solver_;

    const Stencil& getStencil(int i_r) const;
    int getNonZeroCountCircleAsc(const int i_r) const;
    int getNonZeroCountRadialAsc(const int i_theta) const;

    int getCircleAscIndex(const int i_r, const int i_theta) const;
    int getRadialAscIndex(const int i_r, const int i_theta) const;

    void buildAscMatrices();
    void buildAscCircleSection(const int i_r);
    void buildAscRadialSection(const int i_theta);

    void applyAscOrthoCircleSection(const int i_r, const SmootherColor smoother_color, const Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp);
    void applyAscOrthoRadialSection(const int i_theta, const SmootherColor smoother_color, const Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp);
    
    void solveCircleSection(const int i_r, Vector<double>& x, Vector<double>& temp, Vector<double>& solver_storage_1, Vector<double>& solver_storage_2);
    void solveRadialSection(const int i_theta, Vector<double>& x, Vector<double>& temp, Vector<double>& solver_storage);

    void initializeMumpsSolver(DMUMPS_STRUC_C& mumps_solver, const SparseMatrix<double>& solver_matrix);
    void finalizeMumpsSolver(DMUMPS_STRUC_C& mumps_solver);

    /* Paralelization */
    std::vector<int> circle_black_Asc_;
    std::vector<int> circle_white_Asc_;
    std::vector<int> circle_smoother_;
    std::vector<int> radial_black_Asc_;
    std::vector<int> radial_white_Asc_;
    std::vector<int> radial_smoother_;

    bool decrement(std::vector<int>& dependency_counter, const int index);

    void start_radials(Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp);
    void radial_black_Asc0(const int radial_index, Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp);
    void radial_black_Asc1(const int radial_index, Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp);
    void radial_black_Asc2(const int radial_index, Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp);


};

