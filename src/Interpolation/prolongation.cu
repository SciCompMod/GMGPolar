#include "../../include/Interpolation/interpolation.h"

#include "../../include/LinearAlgebra/Vector/gpu_vector_operations.h"

__global__ void applyProlongation_kernel(PolarGrid* fine_grid, double* result, PolarGrid* coarse_grid, double* x) {
    int i_r = blockIdx.x * blockDim.x + threadIdx.x;
    int i_theta = blockIdx.y * blockDim.y + threadIdx.y;

    if (i_r >= fine_grid->nr() || i_theta >= fine_grid->ntheta()) return;

    int i_r_coarse = i_r >> 1;
    int i_theta_coarse = i_theta >> 1;

    double value = 0.0;
    double divisor = 1.0;

    if (i_r & 1) {  // i_r is odd
        double h1 = fine_grid->radialSpacing(i_r - 1);
        double h2 = fine_grid->radialSpacing(i_r);
        if (i_theta & 1) {  // i_theta is odd
            double k1 = fine_grid->angularSpacing(i_theta - 1);
            double k2 = fine_grid->angularSpacing(i_theta);
            divisor = (h1 + h2) * (k1 + k2);

            value = h1 * k1 * x[coarse_grid->index(i_r_coarse, i_theta_coarse)] +
                    h2 * k1 * x[coarse_grid->index(i_r_coarse + 1, i_theta_coarse)] +
                    h1 * k2 * x[coarse_grid->index(i_r_coarse, i_theta_coarse + 1)] +
                    h2 * k2 * x[coarse_grid->index(i_r_coarse + 1, i_theta_coarse + 1)];
        } else {  // i_theta is even
            divisor = (h1 + h2);
            value = h1 * x[coarse_grid->index(i_r_coarse, i_theta_coarse)] +
                    h2 * x[coarse_grid->index(i_r_coarse + 1, i_theta_coarse)];
        }
    } else {  // i_r is even
        if (i_theta & 1) {  // i_theta is odd
            double k1 = fine_grid->angularSpacing(i_theta - 1);
            double k2 = fine_grid->angularSpacing(i_theta);
            divisor = (k1 + k2);

            value = k1 * x[coarse_grid->index(i_r_coarse, i_theta_coarse)] +
                    k2 * x[coarse_grid->index(i_r_coarse, i_theta_coarse + 1)];
        } else {  // i_theta is even
            value = x[coarse_grid->index(i_r_coarse, i_theta_coarse)];
        }
    }

    result[fine_grid->index(i_r, i_theta)] = value / divisor;
}

/* Remark: This injection is not scaled. */
void Interpolation::applyProlongation(const Level& fromLevel, const Level& toLevel, GPU_Vector<double>& result, const GPU_Vector<double>& x) const
{
    assert(toLevel.level() == fromLevel.level() - 1);

    const PolarGrid& coarse_grid = fromLevel.grid();
    const PolarGrid& fine_grid = toLevel.grid();

    assert(x.size() == coarse_grid.numberOfNodes());
    assert(result.size() == fine_grid.numberOfNodes());
    
    dim3 threadsPerBlock(16, 16);
    dim3 numBlocks((fine_grid.nr() + threadsPerBlock.x - 1) / threadsPerBlock.x,
                   (fine_grid.ntheta() + threadsPerBlock.y - 1) / threadsPerBlock.y);

    applyProlongation_kernel<<<numBlocks, threadsPerBlock>>>(
        toLevel.device_grid(), result.data(), fromLevel.device_grid(), x.data());

    cudaDeviceSynchronize();
}