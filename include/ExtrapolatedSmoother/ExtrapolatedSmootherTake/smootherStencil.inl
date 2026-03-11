#pragma once

template <concepts::DomainGeometry DomainGeometry>
const Stencil& ExtrapolatedSmootherTake<DomainGeometry>::getStencil(int i_r, int i_theta) const
{
    const PolarGrid& grid     = ExtrapolatedSmoother<DomainGeometry>::grid_;
    const bool DirBC_Interior = ExtrapolatedSmoother<DomainGeometry>::DirBC_Interior_;

    // Only i_r = 0 is implemented.
    // Stencils are only used to obtain offsets to index into COO/CSR matrices.
    // The inner boundary requires a COO/CSR matrix (rather than a tridiagonal one)
    // because it has an additional across-origin coupling.
    assert(i_r == 0);

    // The across-origin discretization used by the extrapolated smoother requires
    // ntheta to be divisible by 4. This ensures that nodes mapped across the origin
    // preserve their grid classification, i.e., coarse nodes couple only to coarse
    // nodes and fine nodes couple only to fine nodes.
    // Without this assumption, it is a bit more complex to implement the stencil.
    assert((grid.ntheta() / 2) % 2 == 0);

    if (i_theta % 2 == 0) {
        return stencil_center_;
    }
    else {
        if (!DirBC_Interior) {
            return stencil_center_left_;
        }
        else {
            return stencil_center_;
        }
    }
}

template <concepts::DomainGeometry DomainGeometry>
int ExtrapolatedSmootherTake<DomainGeometry>::getNonZeroCountCircleAsc(int i_r) const
{
    const PolarGrid& grid     = ExtrapolatedSmoother<DomainGeometry>::grid_;
    const bool DirBC_Interior = ExtrapolatedSmoother<DomainGeometry>::DirBC_Interior_;

    // Only i_r = 0 is implemented.
    // The number of nonzero elements is only needed to construct COO matrices.
    // The inner boundary requires a COO/CSR matrix (rather than a tridiagonal one)
    // because it has an additional across-origin coupling.
    assert(i_r == 0);

    // The across-origin discretization used by the extrapolated smoother requires
    // ntheta to be divisible by 4. This ensures that nodes mapped across the origin
    // preserve their grid classification, i.e., coarse nodes couple only to coarse
    // nodes and fine nodes couple only to fine nodes.
    // Without this assumption, it is a bit more complex to implement the stencil.
    assert((grid.ntheta() / 2) % 2 == 0);

    if (!DirBC_Interior) {
        return grid.ntheta() / 2 + 2 * (grid.ntheta() / 2);
    }
    else {
        return grid.ntheta();
    }
}

template <concepts::DomainGeometry DomainGeometry>
int ExtrapolatedSmootherTake<DomainGeometry>::getCircleAscIndex(int i_r, int i_theta) const
{
    const PolarGrid& grid     = ExtrapolatedSmoother<DomainGeometry>::grid_;
    const bool DirBC_Interior = ExtrapolatedSmoother<DomainGeometry>::DirBC_Interior_;

    // Only i_r = 0 is implemented.
    // getCircleAscIndex accumulates all stencil sizes within a line up to, but excluding the current node.
    // It is only used to obtain a ptr to index into COO matrices.
    // The inner boundary requires a COO/CSR matrix (rather than a tridiagonal one)
    // because it has an additional across-origin coupling.
    assert(i_r == 0);

    // The across-origin discretization used by the extrapolated smoother requires
    // ntheta to be divisible by 4. This ensures that nodes mapped across the origin
    // preserve their grid classification, i.e., coarse nodes couple only to coarse
    // nodes and fine nodes couple only to fine nodes.
    // Without this assumption, it is a bit more complex to implement the stencil.
    assert((grid.ntheta() / 2) % 2 == 0);

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
