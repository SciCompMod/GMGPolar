#pragma once

#include "../directSolver.h"

namespace gmgpolar
{

template <class LevelCacheType>
class DirectSolverTake : public DirectSolver<LevelCacheType>
{
public:
    explicit DirectSolverTake(const PolarGrid& grid, const LevelCacheType& level_cache, bool DirBC_Interior);

    // Note: The rhs (right-hand side) vector gets overwritten during the solution process.
    void solveInPlace(Vector<double> solution) override;
    void updateValues() override;

private:
#ifdef GMGPOLAR_USE_MUMPS
    using SystemMatrix = SparseMatrixCOO<double>;
    using SystemSolver = CooMumpsSolver;
#else
    using SystemMatrix = SparseMatrixCSR<double>;
    using SystemSolver = SparseLUSolver<double>;
#endif
    SystemMatrix system_matrix_;
    SystemSolver system_solver_;

public:
    // Constructs a symmetric solver matrix.
    SystemMatrix buildSolverMatrix();

    // Adjusts the right-hand side vector for symmetry corrections.
    // This modifies the system from
    //    A * solution = rhs
    // to the equivalent system
    //    symmetric_DBc(A) * solution = rhs - applySymmetryShift(rhs).
    // The correction modifies the rhs to account for the influence of the Dirichlet boundary conditions,
    // ensuring that the solution at the boundary is correctly adjusted and maintains the required symmetry.
    void applySymmetryShift(Vector<double> rhs) const;
    void applySymmetryShiftInnerBoundary(Vector<double> rhs) const;
    void applySymmetryShiftOuterBoundary(Vector<double> rhs) const;

private:
    // Shared kernel: fills solver_matrix's entries from the current grid +
    // LevelCache data. Used by buildSolverMatrix() (on freshly allocated
    // storage) and updateValues() (on existing storage, values only).
    void fillSolverMatrix(SystemMatrix& solver_matrix);
};

#include "applySymmetryShift.inl"
#include "buildSolverMatrix.inl"
#include "directSolverTake.inl"

} // namespace gmgpolar
