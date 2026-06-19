#pragma once

#include "../directSolver.h"

namespace gmgpolar
{

template <class LevelCacheType>
class DirectSolverGive : public DirectSolver<LevelCacheType>
{
public:
    explicit DirectSolverGive(const PolarGrid& grid, const LevelCacheType& level_cache, bool DirBC_Interior);

    // Note: The rhs (right-hand side) vector gets overwritten during the solution process.
    void solveInPlace(Vector<double> solution) override;

private:
#ifdef GMGPOLAR_USE_MUMPS
    using SystemMatrix = SparseMatrixCOO<double>;
    using SystemSolver = CooMumpsSolver;
#else
    using SystemMatrix = SparseMatrixCSR<double>;
    using SystemSolver = SparseLUSolver<double>;
    // Stored only for the in-house solver (CSR).
    SystemMatrix system_matrix_;
#endif

    // Solver object (owns matrix if MUMPS, references if in-house solver).
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
    void applySymmetryShiftInnerBoundary(Vector<double> x) const;
    void applySymmetryShiftOuterBoundary(Vector<double> x) const;
};

#include "applySymmetryShift.inl"
#include "buildSolverMatrix.inl"
#include "directSolverGive.inl"

} // namespace gmgpolar
