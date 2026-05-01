#pragma once

/* -------------- */
/* Stencil access */
/* -------------- */

// Select correct stencil depending on the grid position.
static inline const Stencil& getStencil(int i_r, bool DirBC_Interior)
{
    // Only i_r = 0 is implemented.
    // Stencils are only used to obtain offsets to index into COO/CSR matrices.
    // The inner boundary requires a COO/CSR matrix (rather than a tridiagonal one)
    // because it has an additional across-origin coupling.
    assert(i_r == 0);

    /* ------------------- */
    /* Stencil definitions */
    /* ------------------- */

    // The stencil definitions must be defined before the declaration of the inner_boundary_mumps_solver_,
    // since the mumps solver will be build in the member initializer of the Smoother class.

    // Stencils encode neighborhood connectivity for A_sc matrix assembly.
    // It is only used in the construction of COO/CSR matrices.
    // Thus it is only used for the interior boundary matrix and not needed for the tridiagonal matrices.
    // The Stencil class stores the offset for each position.
    // - Non-zero matrix indicesare obtained via `ptr + offset`
    // - A offset value of `-1` means the position is not included in the stencil pattern.
    // - Other values (0, 1, 2, ..., stencil_size - 1) correspond to valid stencil indices.

    // clang-format off
    static const Stencil stencil_DB = {
        -1, -1, -1,
        -1,  0, -1,
        -1, -1, -1
    };
    static const Stencil circle_stencil_across_origin = {
        -1,  3, -1,
         1,  0, -1,
        -1,  2, -1
    };
    // clang-format on

    return DirBC_Interior ? stencil_DB : circle_stencil_across_origin;
}

// Number of nonzero A_sc entries.
static inline int getNonZeroCountCircleAsc(int i_r, const PolarGrid& grid, bool DirBC_Interior)
{
    // Only i_r = 0 is implemented.
    // The number of nonzero elements is only needed to construct COO matrices.
    // The inner boundary requires a COO/CSR matrix (rather than a tridiagonal one)
    // because it has an additional across-origin coupling.
    assert(i_r == 0);

    const int size_stencil_inner_boundary = DirBC_Interior ? 1 : 4;
    return size_stencil_inner_boundary * grid.ntheta();
}

// Obtain a ptr to index into COO matrices.
// It accumulates all stencil sizes within a line up to, but excluding the current node.
static inline int getCircleAscIndex(int i_r, int i_theta, bool DirBC_Interior)
{
    // Only i_r = 0 is implemented.
    // getCircleAscIndex accumulates all stencil sizes within a line up to, but excluding the current node.
    // It is only used to obtain a ptr to index into COO matrices.
    // The inner boundary requires a COO/CSR matrix (rather than a tridiagonal one)
    // because it has an additional across-origin coupling.
    assert(i_r == 0);

    const int size_stencil_inner_boundary = DirBC_Interior ? 1 : 4;
    return size_stencil_inner_boundary * i_theta;
}
