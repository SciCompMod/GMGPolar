#include "../../include/Interpolation/interpolation.h"

#include "../../include/LinearAlgebra/Vector/gpu_vector_operations.h"

__global__ void applyExtrapolatedRestriction_kernel(PolarGrid* coarse_grid, double* result, PolarGrid* fine_grid, double* x) {
    int i_r_coarse = blockIdx.x * blockDim.x + threadIdx.x;
    int i_theta_coarse = blockIdx.y * blockDim.y + threadIdx.y;

    if (i_r_coarse < coarse_grid->nr() && i_theta_coarse < coarse_grid->ntheta()) {
        int i_r = i_r_coarse << 1;
        int i_theta = i_theta_coarse << 1;

        int i_theta_M1 = fine_grid->wrapThetaIndex(i_theta-1);
        int i_theta_P1 = fine_grid->wrapThetaIndex(i_theta+1);

        // Center, Bottom, Top
        double value = x[fine_grid->index(i_r,i_theta)] +
            0.5 * x[fine_grid->index(i_r, i_theta_M1)] +
            0.5 * x[fine_grid->index(i_r, i_theta_P1)];

        if(i_r_coarse > 0){
            // Left Part
            // Left, Top Left
            value += 
                0.5 * x[fine_grid->index(i_r-1, i_theta)] +
                0.5 * x[fine_grid->index(i_r-1, i_theta_P1)];   
        } 
        if(i_r_coarse < coarse_grid->nr() - 1){
            // Right Part
            // Right, Bottom Right
            value += 
                0.5 * x[fine_grid->index(i_r+1, i_theta)] +
                0.5 * x[fine_grid->index(i_r+1, i_theta_M1)];
        }
        result[coarse_grid->index(i_r_coarse,i_theta_coarse)] = value;
    }
}

void Interpolation::applyExtrapolatedRestriction(const Level& fromLevel, const Level& toLevel, GPU_Vector<double>& result, const GPU_Vector<double>& x) const{
    assert(toLevel.level() == fromLevel.level() + 1);

    const PolarGrid& fine_grid = fromLevel.grid();
    const PolarGrid& coarse_grid = toLevel.grid();

    assert(x.size() == fine_grid.numberOfNodes());
    assert(result.size() == coarse_grid.numberOfNodes());
    
    dim3 threadsPerBlock(16, 16);
    dim3 numBlocks((coarse_grid.nr() + threadsPerBlock.x - 1) / threadsPerBlock.x,
                   (coarse_grid.ntheta() + threadsPerBlock.y - 1) / threadsPerBlock.y);

    applyExtrapolatedRestriction_kernel<<<numBlocks, threadsPerBlock>>>(
        toLevel.device_grid(), result.data(), fromLevel.device_grid(), x.data());

    cudaDeviceSynchronize();
}
