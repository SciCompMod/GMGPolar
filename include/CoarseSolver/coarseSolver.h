#pragma once

class LevelCache;
class Level;

#include "../PolarGrid/polargrid.h"

#include "../InputFunctions/domainGeometry.h"
#include "../InputFunctions/systemParameters.h"

#include "../LinearAlgebra/vector.h"
#include "../LinearAlgebra/matrix.h"
#include "../LinearAlgebra/operations.h"

#include "../common/constants.h"

#include "../Level/level.h"

#include "../Stencil/stencil.h"

#include "../TaskDistribution/taskDistribution.h"

#include "mpi.h" 
#include "dmumps_c.h"   

#include <chrono>
#include <vector>
#include <iostream>

class CoarseSolver {
public:
    explicit CoarseSolver(const PolarGrid& grid, const LevelCache& level_data, 
        const DomainGeometry& domain_geometry, const SystemParameters& system_parameters, const bool DirBC_Interior, 
        const int maxOpenMPThreads, const int openMPTaskThreads
    );
    ~CoarseSolver();

    void solveInPlace(Vector<double>& x);

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

    /* --------------------- */
    /* Coarse Solver members */
    DMUMPS_STRUC_C mumps_;
    SparseMatrix<double> matrixA_;

    const Stencil& get_stencil(int i_r) const;
    int nnz_matrixA() const;
    int ptr_nz_index_matrixA(const int i_r, const int i_theta) const;

    void buildMatrixA(SparseMatrix<double>& matrixA);
    void subtractSymmetryShift(Vector<double>& rhs);

    void initializeMumps(DMUMPS_STRUC_C& mumps, const SparseMatrix<double>& matrixA);
    void solveMumps(Vector<double>& result_rhs);
    void deleteMumps(DMUMPS_STRUC_C& mumps);
};