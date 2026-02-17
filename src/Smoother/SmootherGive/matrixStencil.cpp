#include "../../../include/Smoother/SmootherGive/smootherGive.h"

#include "../../../include/Smoother/SmootherTake/smootherTake.h"

const Stencil& SmootherGive::getStencil(int i_r) const
{
    // Only i_r = 0 is implemented.
    // Stencils are only used to obtain offsets to index into COO/CSR matrices.
    // The inner boundary requires a COO/CSR matrix (rather than a tridiagonal one)
    // because it has an additional across-origin coupling.
    assert(i_r == 0);

    return DirBC_Interior_ ? stencil_DB_ : circle_stencil_across_origin_;
}

int SmootherGive::getNonZeroCountCircleAsc(int i_r) const
{
    // Only i_r = 0 is implemented.
    // The number of nonzero elements is only needed to construct COO matrices.
    // The inner boundary requires a COO/CSR matrix (rather than a tridiagonal one)
    // because it has an additional across-origin coupling.
    assert(i_r == 0);

    const int size_stencil_inner_boundary = DirBC_Interior_ ? 1 : 4;
    return size_stencil_inner_boundary * grid_.ntheta();
}

int SmootherGive::getCircleAscIndex(int i_r, int i_theta) const
{
    // Only i_r = 0 is implemented.
    // getCircleAscIndex accumulates all stencil sizes within a line up to, but excluding the current node.
    // It is only used to obtain a ptr to index into COO matrices.
    // The inner boundary requires a COO/CSR matrix (rather than a tridiagonal one)
    // because it has an additional across-origin coupling.
    assert(i_r == 0);

    const int size_stencil_inner_boundary = DirBC_Interior_ ? 1 : 4;
    return size_stencil_inner_boundary * i_theta;
}
