#pragma once

namespace extrapolated_smoother_take
{

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
static KOKKOS_INLINE_FUNCTION const Stencil& getStencil(const int i_r, const int i_theta, const PolarGrid<Kokkos::HostSpace>& grid,
                                                        const bool DirBC_Interior)
{
    // clang-format off
    static constexpr Stencil    stencil_center = {
        -1, -1, -1,
        -1,  0, -1,
        -1, -1, -1
    };
    static constexpr Stencil  stencil_center_left = {
        -1, -1, -1,
        1,  0, -1,
        -1, -1, -1
    };
    // clang-format on

    // Only i_r = 0 is implemented.
    // Stencils are only used to obtain offsets to index into COO/CSR matrices.
    // The inner boundary requires a COO/CSR matrix (rather than a tridiagonal one)
    // because it has an additional across-origin coupling.
    KOKKOS_ASSERT(i_r == 0);

    // The across-origin discretization used by the extrapolated smoother requires
    // ntheta to be divisible by 4. This ensures that nodes mapped across the origin
    // preserve their grid classification, i.e., coarse nodes couple only to coarse
    // nodes and fine nodes couple only to fine nodes.
    // Without this assumption, it is a bit more complex to implement the stencil.
    KOKKOS_ASSERT((grid.ntheta() / 2) % 2 == 0);

    if (i_theta % 2 == 0) {
        return stencil_center;
    }
    else {
        if (!DirBC_Interior) {
            return stencil_center_left;
        }
        else {
            return stencil_center;
        }
    }
}

static KOKKOS_INLINE_FUNCTION int getNonZeroCountCircleAsc(const int i_r, const PolarGrid<Kokkos::HostSpace>& grid,
                                                           const bool DirBC_Interior)
{
    // Only i_r = 0 is implemented.
    // The number of nonzero elements is only needed to construct COO matrices.
    // The inner boundary requires a COO/CSR matrix (rather than a tridiagonal one)
    // because it has an additional across-origin coupling.
    KOKKOS_ASSERT(i_r == 0);

    // The across-origin discretization used by the extrapolated smoother requires
    // ntheta to be divisible by 4. This ensures that nodes mapped across the origin
    // preserve their grid classification, i.e., coarse nodes couple only to coarse
    // nodes and fine nodes couple only to fine nodes.
    // Without this assumption, it is a bit more complex to implement the stencil.
    KOKKOS_ASSERT((grid.ntheta() / 2) % 2 == 0);

    if (!DirBC_Interior) {
        return grid.ntheta() / 2 + 2 * (grid.ntheta() / 2);
    }
    else {
        return grid.ntheta();
    }
}

static KOKKOS_INLINE_FUNCTION int getCircleAscIndex(const int i_r, const int i_theta, const PolarGrid<Kokkos::HostSpace>& grid,
                                                    const bool DirBC_Interior)
{
    // Only i_r = 0 is implemented.
    // getCircleAscIndex accumulates all stencil sizes within a line up to, but excluding the current node.
    // It is only used to obtain a ptr to index into COO matrices.
    // The inner boundary requires a COO/CSR matrix (rather than a tridiagonal one)
    // because it has an additional across-origin coupling.
    KOKKOS_ASSERT(i_r == 0);

    // The across-origin discretization used by the extrapolated smoother requires
    // ntheta to be divisible by 4. This ensures that nodes mapped across the origin
    // preserve their grid classification, i.e., coarse nodes couple only to coarse
    // nodes and fine nodes couple only to fine nodes.
    // Without this assumption, it is a bit more complex to implement the stencil.
    KOKKOS_ASSERT((grid.ntheta() / 2) % 2 == 0);

    if (!DirBC_Interior) {
        if (i_theta % 2 == 0) {
            return 3 * (i_theta / 2);
        }
        else {
            return 3 * (i_theta / 2) + 1;
        }
    }
    else {
        return i_theta;
    }
}

} // namespace extrapolated_smoother_take