#include "../../include/Interpolation/interpolation.h"

#include "../../include/LinearAlgebra/Vector/gpu_vector_operations.h"

__global__ void applyExtrapolatedProlongation_kernel(PolarGrid* fine_grid, double* result, PolarGrid* coarse_grid, double* x) {
    int i_r = blockIdx.x * blockDim.x + threadIdx.x;
    int i_theta = blockIdx.y * blockDim.y + threadIdx.y;

    if (i_r >= fine_grid->nr() || i_theta >= fine_grid->ntheta()) return;

    int i_r_coarse = i_r >> 1;
    int i_theta_coarse = i_theta >> 1;

    double value;
    double multiplier;

    if (i_r & 1) {  // i_r is odd
        if (i_theta & 1) {  // i_theta is odd
            value = x[coarse_grid->index(i_r_coarse+1, i_theta_coarse)] + x[coarse_grid->index(i_r_coarse, i_theta_coarse+1)];
            multiplier = 0.5;
        } else {  // i_theta is even
            value = x[coarse_grid->index(i_r_coarse, i_theta_coarse)] + x[coarse_grid->index(i_r_coarse+1, i_theta_coarse)];
            multiplier = 0.5;
        }
    } else {  // i_r is even
        if (i_theta & 1) {  // i_theta is odd
            value = x[coarse_grid->index(i_r_coarse, i_theta_coarse)] + x[coarse_grid->index(i_r_coarse, i_theta_coarse+1)];
            multiplier = 0.5;
        } else {  // i_theta is even
            value = x[coarse_grid->index(i_r_coarse, i_theta_coarse)];
            multiplier = 1.0;
        }
    }

    result[fine_grid->index(i_r, i_theta)] = multiplier * value;
}

/* Remark: This injection is not scaled. */
void Interpolation::applyExtrapolatedProlongation(const Level& fromLevel, const Level& toLevel, GPU_Vector<double>& result, const GPU_Vector<double>& x) const
{
    assert(toLevel.level() == fromLevel.level() - 1);

    const PolarGrid& coarse_grid = fromLevel.grid();
    const PolarGrid& fine_grid = toLevel.grid();

    assert(x.size() == coarse_grid.numberOfNodes());
    assert(result.size() == fine_grid.numberOfNodes());
    
    dim3 threadsPerBlock(16, 16);
    dim3 numBlocks((fine_grid.nr() + threadsPerBlock.x - 1) / threadsPerBlock.x,
                   (fine_grid.ntheta() + threadsPerBlock.y - 1) / threadsPerBlock.y);

    applyExtrapolatedProlongation_kernel<<<numBlocks, threadsPerBlock>>>(
        toLevel.device_grid(), result.data(), fromLevel.device_grid(), x.data());

    cudaDeviceSynchronize();
}