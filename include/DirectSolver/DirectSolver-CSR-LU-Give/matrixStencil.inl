#pragma once

template <class LevelCacheType>
int DirectSolver_CSR_LU_Give<LevelCacheType>::getNonZeroCountSolverMatrix() const
{
    const PolarGrid& grid     = DirectSolver<LevelCacheType>::grid_;
    const bool DirBC_Interior = DirectSolver<LevelCacheType>::DirBC_Interior_;

    const int size_stencil_inner_boundary      = DirBC_Interior ? 1 : 7;
    const int size_stencil_next_inner_boundary = DirBC_Interior ? 9 : 9;
    const int size_stencil_interior            = 9;
    const int size_stencil_next_outer_boundary = 9;
    const int size_stencil_outer_boundary      = 1;

    assert(grid.nr() >= 4);

    return grid.ntheta() *
           (size_stencil_inner_boundary + size_stencil_next_inner_boundary + (grid.nr() - 4) * size_stencil_interior +
            size_stencil_next_outer_boundary + size_stencil_outer_boundary);
}
