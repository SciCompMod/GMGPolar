#pragma once

template <class LevelCacheType>
const Stencil& SmootherGive<LevelCacheType>::getStencil(int i_r) const
{
    // Only i_r = 0 is implemented.
    // Stencils are only used to obtain offsets to index into COO/CSR matrices.
    // The inner boundary requires a COO/CSR matrix (rather than a tridiagonal one)
    // because it has an additional across-origin coupling.
    assert(i_r == 0);

    const bool DirBC_Interior = Smoother<LevelCacheType>::DirBC_Interior_;
    return DirBC_Interior ? stencil_DB_ : circle_stencil_across_origin_;
}

template <class LevelCacheType>
int SmootherGive<LevelCacheType>::getNonZeroCountCircleAsc(int i_r) const
{
    // Only i_r = 0 is implemented.
    // The number of nonzero elements is only needed to construct COO matrices.
    // The inner boundary requires a COO/CSR matrix (rather than a tridiagonal one)
    // because it has an additional across-origin coupling.
    assert(i_r == 0);

    const PolarGrid& grid                 = Smoother<LevelCacheType>::grid_;
    const bool DirBC_Interior             = Smoother<LevelCacheType>::DirBC_Interior_;
    const int size_stencil_inner_boundary = DirBC_Interior ? 1 : 4;
    return size_stencil_inner_boundary * grid.ntheta();
}

template <class LevelCacheType>
int SmootherGive<LevelCacheType>::getCircleAscIndex(int i_r, int i_theta) const
{
    // Only i_r = 0 is implemented.
    // getCircleAscIndex accumulates all stencil sizes within a line up to, but excluding the current node.
    // It is only used to obtain a ptr to index into COO matrices.
    // The inner boundary requires a COO/CSR matrix (rather than a tridiagonal one)
    // because it has an additional across-origin coupling.
    assert(i_r == 0);

    const bool DirBC_Interior             = Smoother<LevelCacheType>::DirBC_Interior_;
    const int size_stencil_inner_boundary = DirBC_Interior ? 1 : 4;
    return size_stencil_inner_boundary * i_theta;
}
