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
#include "../LinearAlgebra/vector.h"
#include "../LinearAlgebra/matrix.h"
#include "../LinearAlgebra/vector_operations.h"
#include "../LinearAlgebra/symmetricTridiagonalSolver.h"
#include "../LinearAlgebra/diagonalSolver.h"

#include "../common/constants.h"
#include "../Level/level.h"
#include "../Stencil/stencil.h"
#include "../TaskDistribution/taskDistribution.h"

class ExtrapolatedSmoother {
public:
    explicit ExtrapolatedSmoother(const PolarGrid& grid, const LevelCache& level_cache, 
                                  const DomainGeometry& domain_geometry,
                                  bool DirBC_Interior, int num_omp_threads);
    ~ExtrapolatedSmoother();

    void extrapolatedSmoothingInPlace(Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp);

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

    /* ----------------------------- */
    /* Extrapolated Smoother members */
    SparseMatrix<double> inner_boundary_circle_matrix_;
    DMUMPS_STRUC_C inner_boundary_mumps_solver_;
    std::vector<DiagonalSolver<double>> circle_diagonal_solver_;
    std::vector<DiagonalSolver<double>> radial_diagonal_solver_;
    std::vector<SymmetricTridiagonalSolver<double>> circle_tridiagonal_solver_;
    std::vector<SymmetricTridiagonalSolver<double>> radial_tridiagonal_solver_;

    const Stencil& getStencil(int i_r, int i_theta) const;
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
};

