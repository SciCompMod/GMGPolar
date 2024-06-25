#pragma once

class LevelCache;
class Level;

#include "../PolarGrid/polargrid.h"

#include "../InputFunctions/domainGeometry.h"
#include "../InputFunctions/systemParameters.h"

#include "../LinearAlgebra/vector.h"
#include "../LinearAlgebra/matrix.h"
#include "../LinearAlgebra/operations.h"
#include "../LinearAlgebra/symmetricTridiagonalSolver.h"

#include "../common/constants.h"

#include "../Level/level.h"

#include "../Stencil/stencil.h"

#include "../TaskDistribution/taskDistribution.h"

#include "mpi.h" 
#include "dmumps_c.h"   
// #include "mumps_c_types.h"

#include <chrono>
#include <vector>
#include <iostream>

class Smoother {
public:
    explicit Smoother(const PolarGrid& grid, const LevelCache& level_data, 
        const DomainGeometry& domain_geometry, const SystemParameters& system_parameters, const bool DirBC_Interior, 
        const int maxOpenMPThreads, const int openMPTaskThreads);
    ~Smoother();

    void smoothingInPlace(Vector<double>& x, Vector<double>& temp_rhs);

private:
    /* ------------------- */
    /* Constructor members */
    const PolarGrid& grid_;
    /* Level Cache Data */
    const std::vector<double>& sin_theta_;
    const std::vector<double>& cos_theta_;
    const Vector<double>& rhs_f_;
    const std::vector<double>& u_D_;
    const std::vector<double>& u_D_Interior_;

    const DomainGeometry& domain_geometry_;
    const SystemParameters& system_parameters_;
    const bool DirBC_Interior_;

    const int maxOpenMPThreads_;
    const int openMPTaskThreads_;

    /* ---------------- */
    /* Smoother members */
    SparseMatrix<double> inner_boundary_circle_Asc_matrix_;
    DMUMPS_STRUC_C inner_boundary_circle_Asc_mumps_;

    std::vector<SymmetricTridiagonalSolver<double>> circle_symmetric_cyclic_tridiagonal_solver_;
    std::vector<SymmetricTridiagonalSolver<double>> radial_symmetric_tridiagonal_solver_;

    const Stencil& get_stencil(int i_r) const;
    int nnz_circle_Asc(const int i_r) const;
    int nnz_radial_Asc(const int i_theta) const;

    int ptr_nz_index_circle_Asc(const int i_r, const int i_theta) const;
    int ptr_nz_index_radial_Asc(const int i_r, const int i_theta) const;

    void build_Asc_matrices();
    void build_Asc_circle_section(const int i_r);
    void build_Asc_radial_section(const int i_theta);

    void initializeMumps(DMUMPS_STRUC_C& Asc_mumps, const SparseMatrix<double>& Asc_matrix);
    void deleteMumps(DMUMPS_STRUC_C& Asc_mumps);
};

