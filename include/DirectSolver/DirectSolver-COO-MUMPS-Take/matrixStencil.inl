#pragma once

#ifdef GMGPOLAR_USE_MUMPS

template <concepts::DomainGeometry DomainGeometry>
const Stencil& DirectSolver_COO_MUMPS_Take<DomainGeometry>::getStencil(int i_r) const
{
    const PolarGrid& grid     = DirectSolver<DomainGeometry>::grid_;
    const bool DirBC_Interior = DirectSolver<DomainGeometry>::DirBC_Interior_;

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
int DirectSolver_COO_MUMPS_Take<DomainGeometry>::getNonZeroCountSolverMatrix() const
{
    const PolarGrid& grid     = DirectSolver<DomainGeometry>::grid_;
    const bool DirBC_Interior = DirectSolver<DomainGeometry>::DirBC_Interior_;

    const int size_stencil_inner_boundary      = DirBC_Interior ? 1 : 7;
    const int size_stencil_next_inner_boundary = DirBC_Interior ? 6 : 9;
    const int size_stencil_interior            = 9;
    const int size_stencil_next_outer_boundary = 6;
    const int size_stencil_outer_boundary      = 1;

    assert(grid.nr() >= 4);

    return grid.ntheta() *
           (size_stencil_inner_boundary + size_stencil_next_inner_boundary + (grid.nr() - 4) * size_stencil_interior +
            size_stencil_next_outer_boundary + size_stencil_outer_boundary);
}

/* ----------------------------------------------------------------- */
/* If the indexing is not smoother-based, please adjust the indexing */
template <concepts::DomainGeometry DomainGeometry>
int DirectSolver_COO_MUMPS_Take<DomainGeometry>::getSolverMatrixIndex(const int i_r, const int i_theta) const
{
    const PolarGrid& grid     = DirectSolver<DomainGeometry>::grid_;
    const bool DirBC_Interior = DirectSolver<DomainGeometry>::DirBC_Interior_;

    const int size_stencil_inner_boundary      = DirBC_Interior ? 1 : 7;
    const int size_stencil_next_inner_boundary = DirBC_Interior ? 6 : 9;
    const int size_stencil_interior            = 9;
    const int size_stencil_next_outer_boundary = 6;
    const int size_stencil_outer_boundary      = 1;

    assert(grid.nr() >= 4);
    assert(grid.numberSmootherCircles() >= 2);
    assert(grid.lengthSmootherRadial() >= 2);

    if (1 < i_r && i_r < grid.numberSmootherCircles()) {
        // Interior: Circle index section
        const int prior_inner_boundary_nodes      = grid.ntheta();
        const int prior_next_inner_boundary_nodes = grid.ntheta();
        const int prior_interior_nodes            = (i_r - 2) * grid.ntheta() + i_theta;
        return size_stencil_inner_boundary * prior_inner_boundary_nodes +
               size_stencil_next_inner_boundary * prior_next_inner_boundary_nodes +
               size_stencil_interior * prior_interior_nodes;
    }
    if (i_r >= grid.numberSmootherCircles() && i_r < grid.nr() - 2) {
        // Interior: Radial index section
        const int prior_inner_boundary_nodes      = grid.ntheta();
        const int prior_next_inner_boundary_nodes = grid.ntheta();
        const int prior_interior_nodes            = grid.ntheta() * (grid.numberSmootherCircles() - 2) +
                                         i_theta * (grid.lengthSmootherRadial() - 2) + i_r -
                                         grid.numberSmootherCircles();
        const int prior_next_outer_boundary_nodes = i_theta;
        const int prior_outer_boundary_nodes      = i_theta;
        return size_stencil_inner_boundary * prior_inner_boundary_nodes +
               size_stencil_next_inner_boundary * prior_next_inner_boundary_nodes +
               size_stencil_interior * prior_interior_nodes +
               size_stencil_next_outer_boundary * prior_next_outer_boundary_nodes +
               size_stencil_outer_boundary * prior_outer_boundary_nodes;
    }
    else if (i_r == 0) {
        // Inner boundary
        const int prior_inner_boundary_nodes = i_theta;
        return size_stencil_inner_boundary * prior_inner_boundary_nodes;
    }
    else if (i_r == 1) {
        // Next to inner boundary
        const int prior_inner_boundary_nodes      = grid.ntheta();
        const int prior_next_inner_boundary_nodes = i_theta;
        return size_stencil_inner_boundary * prior_inner_boundary_nodes +
               size_stencil_next_inner_boundary * prior_next_inner_boundary_nodes;
    }
    else if (i_r == grid.nr() - 2) {
        // Next to outer boundary
        const int prior_inner_boundary_nodes      = grid.ntheta();
        const int prior_next_inner_boundary_nodes = grid.ntheta();
        const int prior_interior_nodes =
            grid.ntheta() * (grid.numberSmootherCircles() - 2) + (i_theta + 1) * (grid.lengthSmootherRadial() - 2);
        const int prior_next_outer_boundary_nodes = i_theta;
        const int prior_outer_boundary_nodes      = i_theta;
        return size_stencil_inner_boundary * prior_inner_boundary_nodes +
               size_stencil_next_inner_boundary * prior_next_inner_boundary_nodes +
               size_stencil_interior * prior_interior_nodes +
               size_stencil_next_outer_boundary * prior_next_outer_boundary_nodes +
               size_stencil_outer_boundary * prior_outer_boundary_nodes;
    }
    else if (i_r == grid.nr() - 1) {
        // Outer boundary
        const int prior_inner_boundary_nodes      = grid.ntheta();
        const int prior_next_inner_boundary_nodes = grid.ntheta();
        const int prior_interior_nodes =
            grid.ntheta() * (grid.numberSmootherCircles() - 2) + (i_theta + 1) * (grid.lengthSmootherRadial() - 2);
        const int prior_next_outer_boundary_nodes = i_theta + 1;
        const int prior_outer_boundary_nodes      = i_theta;
        return size_stencil_inner_boundary * prior_inner_boundary_nodes +
               size_stencil_next_inner_boundary * prior_next_inner_boundary_nodes +
               size_stencil_interior * prior_interior_nodes +
               size_stencil_next_outer_boundary * prior_next_outer_boundary_nodes +
               size_stencil_outer_boundary * prior_outer_boundary_nodes;
    }
    throw std::out_of_range("Invalid index for stencil");
}

#endif
