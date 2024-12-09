#include "../../include/Interpolation/interpolation.h"

#include "../../include/LinearAlgebra/Vector/gpu_vector_operations.h"

__global__ void applyRestriction_kernel(PolarGrid* coarse_grid, double* result, PolarGrid* fine_grid, double* x) {
    int i_r_coarse = blockIdx.x * blockDim.x + threadIdx.x;
    int i_theta_coarse = blockIdx.y * blockDim.y + threadIdx.y;

    if (i_r_coarse < coarse_grid->nr() && i_theta_coarse < coarse_grid->ntheta()) {
        int i_r = i_r_coarse << 1;
        int i_theta = i_theta_coarse << 1;

        int i_theta_M2 = fine_grid->wrapThetaIndex(i_theta-2);
        int i_theta_M1 = fine_grid->wrapThetaIndex(i_theta-1);
        int i_theta_P1 = fine_grid->wrapThetaIndex(i_theta+1);
        double k1 = fine_grid->angularSpacing(i_theta_M2);
        double k2 = fine_grid->angularSpacing(i_theta_M1);
        double k3 = fine_grid->angularSpacing(i_theta);
        double k4 = fine_grid->angularSpacing(i_theta_P1);

        // Center, Bottom, Top
        double value = x[fine_grid->index(i_r,i_theta)] +
            k2 * x[fine_grid->index(i_r, i_theta_M1)] / (k1+k2) +
            k3 * x[fine_grid->index(i_r, i_theta_P1)] / (k3+k4);

        if(i_r_coarse > 0){
            // Left Part
            double h1 = fine_grid->radialSpacing(i_r-2);
            double h2 = fine_grid->radialSpacing(i_r-1);
            // Left, Bottom Left, Top Left
            value += h2 * x[fine_grid->index(i_r-1, i_theta)] / (h1+h2) +
                h2*k2 * x[fine_grid->index(i_r-1, i_theta_M1)] / ((h1+h2)*(k1+k2)) +
                h2*k3 * x[fine_grid->index(i_r-1, i_theta_P1)] / ((h1+h2)*(k3+k4));                       
        } 
        if(i_r_coarse < coarse_grid->nr() - 1){
            // Right Part
            double h3 = fine_grid->radialSpacing(i_r);
            double h4 = fine_grid->radialSpacing(i_r+1);
            // Right, Bottom Right, Top Right
            value += h3 * x[fine_grid->index(i_r+1, i_theta)] / (h3+h4) +
                h3*k2 * x[fine_grid->index(i_r+1, i_theta_M1)] / ((h3+h4)*(k1+k2)) +
                h3*k3 * x[fine_grid->index(i_r+1, i_theta_P1)] / ((h3+h4)*(k3+k4));
        }
        result[coarse_grid->index(i_r_coarse,i_theta_coarse)] = value;
    }
}

void Interpolation::applyRestriction(const Level& fromLevel, const Level& toLevel, GPU_Vector<double>& result, const GPU_Vector<double>& x) const{
    assert(toLevel.level() == fromLevel.level() + 1);

    const PolarGrid& fine_grid = fromLevel.grid();
    const PolarGrid& coarse_grid = toLevel.grid();

    assert(x.size() == fine_grid.numberOfNodes());
    assert(result.size() == coarse_grid.numberOfNodes());
    
    dim3 threadsPerBlock(16, 16);
    dim3 numBlocks((coarse_grid.nr() + threadsPerBlock.x - 1) / threadsPerBlock.x,
                   (coarse_grid.ntheta() + threadsPerBlock.y - 1) / threadsPerBlock.y);

    applyRestriction_kernel<<<numBlocks, threadsPerBlock>>>(
        toLevel.device_grid(), result.data(), fromLevel.device_grid(), x.data());

    cudaDeviceSynchronize();
}
