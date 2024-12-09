#include "../../include/Interpolation/interpolation.h"

#include "../../include/LinearAlgebra/Vector/gpu_vector_operations.h"

__global__ void applyInjection_kernel(PolarGrid* coarse_grid, double* result, PolarGrid* fine_grid, double* x) {
    int i_r_coarse = blockIdx.x * blockDim.x + threadIdx.x;
    int i_theta_coarse = blockIdx.y * blockDim.y + threadIdx.y;

    if (i_r_coarse < coarse_grid->nr() && i_theta_coarse < coarse_grid->ntheta()) {
        int i_r = i_r_coarse << 1;
        int i_theta = i_theta_coarse << 1;

        int fine_idx = fine_grid->index(i_r, i_theta);
        int coarse_idx = coarse_grid->index(i_r_coarse, i_theta_coarse);

        result[coarse_idx] = x[fine_idx];
    }
}

/* Remark: This injection is not scaled. */
void Interpolation::applyInjection(const Level& fromLevel, const Level& toLevel, GPU_Vector<double>& result, const GPU_Vector<double>& x) const
{
    assert(toLevel.level() == fromLevel.level() + 1);

    const PolarGrid& fine_grid = fromLevel.grid();
    const PolarGrid& coarse_grid = toLevel.grid();

    assert(x.size() == fine_grid.numberOfNodes());
    assert(result.size() == coarse_grid.numberOfNodes());

    dim3 threadsPerBlock(16, 16);
    dim3 numBlocks((coarse_grid.nr() + threadsPerBlock.x - 1) / threadsPerBlock.x,
                   (coarse_grid.ntheta() + threadsPerBlock.y - 1) / threadsPerBlock.y);

    applyInjection_kernel<<<numBlocks, threadsPerBlock>>>(
        toLevel.device_grid(), result.data(), fromLevel.device_grid(), x.data());

    cudaDeviceSynchronize();
}