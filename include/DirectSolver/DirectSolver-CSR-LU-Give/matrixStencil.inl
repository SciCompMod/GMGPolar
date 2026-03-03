#pragma once

template <concepts::DomainGeometry DomainGeometry>
int DirectSolver_CSR_LU_Give<DomainGeometry>::getStencilSize(int global_index) const
{
    const PolarGrid& grid           = DirectSolver<DomainGeometry>::grid_;
    const bool       DirBC_Interior = DirectSolver<DomainGeometry>::DirBC_Interior_;

    int i_r, i_theta;
    grid.multiIndex(global_index, i_r, i_theta);

    const int size_stencil_inner_boundary      = DirBC_Interior ? 1 : 7;
    const int size_stencil_next_inner_boundary = DirBC_Interior ? 9 : 9;
    const int size_stencil_interior            = 9;
    const int size_stencil_next_outer_boundary = 9;
    const int size_stencil_outer_boundary      = 1;

    if ((i_r > 1 && i_r < grid.nr() - 2) || (i_r == 1 && !DirBC_Interior)) {
        return size_stencil_interior;
    }
    else if (i_r == 0 && !DirBC_Interior) {
        return size_stencil_inner_boundary;
    }
    else if ((i_r == 0 && DirBC_Interior) || i_r == grid.nr() - 1) {
        return size_stencil_outer_boundary;
    }
    else if (i_r == 1 && DirBC_Interior) {
        return size_stencil_next_inner_boundary;
    }
    else if (i_r == grid.nr() - 2) {
        return size_stencil_next_outer_boundary;
    }
    throw std::out_of_range("Invalid index for stencil");
}

template <concepts::DomainGeometry DomainGeometry>
const Stencil& DirectSolver_CSR_LU_Give<DomainGeometry>::getStencil(int i_r) const
{
    const PolarGrid& grid           = DirectSolver<DomainGeometry>::grid_;
    const bool       DirBC_Interior = DirectSolver<DomainGeometry>::DirBC_Interior_;

    assert(0 <= i_r && i_r < grid.nr());
    assert(grid.nr() >= 4);

    if ((i_r > 1 && i_r < grid.nr() - 2) || (i_r == 1 && !DirBC_Interior)) {
        return stencil_interior_;
    }
    else if (i_r == 0 && !DirBC_Interior) {
        return stencil_across_origin_;
    }
    else if ((i_r == 0 && DirBC_Interior) || i_r == grid.nr() - 1) {
        return stencil_DB_;
    }
    else if (i_r == 1 && DirBC_Interior) {
        return stencil_next_inner_DB_;
    }
    else if (i_r == grid.nr() - 2) {
        return stencil_next_outer_DB_;
    }
    throw std::out_of_range("Invalid index for stencil");
}

template <concepts::DomainGeometry DomainGeometry>
int DirectSolver_CSR_LU_Give<DomainGeometry>::getNonZeroCountSolverMatrix() const
{
    const PolarGrid& grid           = DirectSolver<DomainGeometry>::grid_;
    const bool       DirBC_Interior = DirectSolver<DomainGeometry>::DirBC_Interior_;

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
