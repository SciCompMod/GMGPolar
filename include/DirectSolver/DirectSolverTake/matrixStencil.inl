#pragma once

namespace direct_solver_take
{

static KOKKOS_INLINE_FUNCTION int getStencilSize(int global_index, const PolarGrid<DefaultMemorySpace>& grid,
                                                 const bool DirBC_Interior)
{
    int i_r, i_theta;
    grid.multiIndex(global_index, i_r, i_theta);

    const int size_stencil_inner_boundary      = DirBC_Interior ? 1 : 7;
    const int size_stencil_next_inner_boundary = DirBC_Interior ? 6 : 9;
    const int size_stencil_interior            = 9;
    const int size_stencil_next_outer_boundary = 6;
    const int size_stencil_outer_boundary      = 1;

    if ((i_r > 1 && i_r < grid.nr() - 2) || (i_r == 1 && !DirBC_Interior)) {
        return size_stencil_interior;
    }
    else if (i_r == 0 && !DirBC_Interior) {
        return size_stencil_inner_boundary;
    }
    else if ((i_r == 0 && DirBC_Interior) || (i_r == grid.nr() - 1)) {
        return size_stencil_outer_boundary;
    }
    else if (i_r == 1 && DirBC_Interior) {
        return size_stencil_next_inner_boundary;
    }
    else if (i_r == grid.nr() - 2) {
        return size_stencil_next_outer_boundary;
    }
    Kokkos::abort("Invalid stencil index");
    return -1;
}

static KOKKOS_INLINE_FUNCTION const Stencil& getStencil(const int i_r, const PolarGrid<DefaultMemorySpace>& grid,
                                                        const bool DirBC_Interior)
{
    KOKKOS_ASSERT(0 <= i_r && i_r < grid.nr());
    KOKKOS_ASSERT(grid.nr() >= 4);

    // clang-format off
    static constexpr Stencil stencil_interior_      = {
        7, 4, 8,
        1, 0, 2,
        5, 3, 6
    };
    static constexpr Stencil stencil_across_origin_ = {
        -1, 4, 6,
        1, 0, 2,
        -1, 3, 5
    };
    static constexpr Stencil stencil_DB_            = {
        -1, -1, -1,
        -1,  0, -1,
        -1, -1, -1
    };
    static constexpr Stencil stencil_next_inner_DB_ = {
        -1,  3,  5,
        -1,  0,  1,
        -1,  2,  4
    };
    static constexpr Stencil stencil_next_outer_DB_ = {
        5,  3, -1,
        1,  0, -1,
        4,  2, -1
    };
    // clang-format on

    if ((i_r > 1 && i_r < grid.nr() - 2) || (i_r == 1 && !DirBC_Interior)) {
        return stencil_interior_;
    }
    else if (i_r == 0 && !DirBC_Interior) {
        return stencil_across_origin_;
    }
    else if ((i_r == 0 && DirBC_Interior) || (i_r == grid.nr() - 1)) {
        return stencil_DB_;
    }
    else if (i_r == 1 && DirBC_Interior) {
        return stencil_next_inner_DB_;
    }
    else {
        //if (i_r == grid.nr() - 2)
        return stencil_next_outer_DB_;
    }
}

static KOKKOS_INLINE_FUNCTION int getNonZeroCountSolverMatrix(const PolarGrid<DefaultMemorySpace>& grid,
                                                              const bool DirBC_Interior)
{
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
static KOKKOS_INLINE_FUNCTION int getSolverMatrixIndex(const int i_r, const int i_theta,
                                                       const PolarGrid<DefaultMemorySpace>& grid,
                                                       const bool DirBC_Interior)
{
    const int size_stencil_inner_boundary      = DirBC_Interior ? 1 : 7;
    const int size_stencil_next_inner_boundary = DirBC_Interior ? 6 : 9;
    const int size_stencil_interior            = 9;
    const int size_stencil_next_outer_boundary = 6;
    const int size_stencil_outer_boundary      = 1;

    assert(grid.nr() >= 4);
    assert(grid.numberSmootherCircles() >= 2);
    assert(grid.lengthRadialSmoother() >= 2);

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
                                         i_theta * (grid.lengthRadialSmoother() - 2) + i_r -
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
            grid.ntheta() * (grid.numberSmootherCircles() - 2) + (i_theta + 1) * (grid.lengthRadialSmoother() - 2);
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
            grid.ntheta() * (grid.numberSmootherCircles() - 2) + (i_theta + 1) * (grid.lengthRadialSmoother() - 2);
        const int prior_next_outer_boundary_nodes = i_theta + 1;
        const int prior_outer_boundary_nodes      = i_theta;
        return size_stencil_inner_boundary * prior_inner_boundary_nodes +
               size_stencil_next_inner_boundary * prior_next_inner_boundary_nodes +
               size_stencil_interior * prior_interior_nodes +
               size_stencil_next_outer_boundary * prior_next_outer_boundary_nodes +
               size_stencil_outer_boundary * prior_outer_boundary_nodes;
    }
    Kokkos::abort("Invalid stencil index");
    return -1;
}

static KOKKOS_INLINE_FUNCTION bool validateSolverMatrixIndexing(const PolarGrid<DefaultMemorySpace>& grid,
                                                                const bool DirBC_Interior)
{
    // 1. Check each node: getSolverMatrixIndex == cumulative sum of prior stencil sizes
    for (int global_index = 0; global_index < grid.numberOfNodes(); ++global_index) {
        int i_r, i_theta;
        grid.multiIndex(global_index, i_r, i_theta);

        int expected = 0;
        for (int prior = 0; prior < global_index; ++prior) {
            expected += getStencilSize(prior, grid, DirBC_Interior);
        }

        if (getSolverMatrixIndex(i_r, i_theta, grid, DirBC_Interior) != expected)
            return false;
        if (getStencilSize(global_index, grid, DirBC_Interior) != getStencil(i_r, grid, DirBC_Interior).size())
            return false;
    }

    // 2. Check total non-zero count
    int total = 0;
    for (int global_index = 0; global_index < grid.numberOfNodes(); ++global_index) {
        total += getStencilSize(global_index, grid, DirBC_Interior);
    }
    if (total != getNonZeroCountSolverMatrix(grid, DirBC_Interior))
        return false;

    return true;
}

} // namespace direct_solver_take
