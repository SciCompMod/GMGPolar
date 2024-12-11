#include "../../include/GMGPolar/gmgpolar.h"

void GMGPolar::extrapolatedResidual(const int current_level, Vector<double>& residual,
                                    const Vector<double>& residual_next_level)
{
    const PolarGrid& fineGrid   = levels_[current_level].grid();
    const PolarGrid& coarseGrid = levels_[current_level + 1].grid();

    assert(residual.size() == fineGrid.numberOfNodes());
    assert(residual_next_level.size() == coarseGrid.numberOfNodes());

#pragma omp parallel
    {
/* Circluar Indexing Section */
/* For loop matches circular access pattern */
#pragma omp for nowait
        for (int i_r = 0; i_r < fineGrid.numberSmootherCircles(); i_r++) {
            int i_r_coarse = i_r >> 1;
            for (int i_theta = 0; i_theta < fineGrid.ntheta(); i_theta++) {
                int i_theta_coarse = i_theta >> 1;

                if (i_r & 1 || i_theta & 1) {
                    residual[fineGrid.index(i_r, i_theta)] *= 4.0 / 3.0;
                }
                else {
                    int fine_idx       = fineGrid.index(i_r, i_theta);
                    int coarse_idx     = coarseGrid.index(i_r_coarse, i_theta_coarse);
                    residual[fine_idx] = (4.0 * residual[fine_idx] - residual_next_level[coarse_idx]) / 3.0;
                }
            }
        }

/* Radial Indexing Section */
/* For loop matches radial access pattern */
#pragma omp for nowait
        for (int i_theta = 0; i_theta < fineGrid.ntheta(); i_theta++) {
            int i_theta_coarse = i_theta >> 1;
            for (int i_r = fineGrid.numberSmootherCircles(); i_r < fineGrid.nr(); i_r++) {
                int i_r_coarse = i_r >> 1;

                if (i_r & 1 || i_theta & 1) {
                    residual[fineGrid.index(i_r, i_theta)] *= 4.0 / 3.0;
                }
                else {
                    int fine_idx       = fineGrid.index(i_r, i_theta);
                    int coarse_idx     = coarseGrid.index(i_r_coarse, i_theta_coarse);
                    residual[fine_idx] = (4.0 * residual[fine_idx] - residual_next_level[coarse_idx]) / 3.0;
                }
            }
        }
    }
}


__global__ void applyExtrapolatedResidual_kernel(
    double* residual, double* residual_next_level,
    PolarGrid* fineGrid, PolarGrid* coarseGrid) 
{

    int i_r = blockIdx.x * 14 + threadIdx.x - 1;
    int i_theta = blockIdx.y * 14 + threadIdx.y - 1;

    if (i_r < 0 || i_r >= fineGrid->nr() || i_theta < 0 || i_theta >= fineGrid->ntheta()) return;

    int i_r_coarse = i_r >> 1;
    int i_theta_coarse = i_theta >> 1;

    int fine_idx = fineGrid->index(i_r, i_theta);

    if (i_r & 1 || i_theta & 1) {
        residual[fine_idx] *= 4.0 / 3.0;
    }
    else {
        int coarse_idx     = coarseGrid->index(i_r_coarse, i_theta_coarse);
        residual[fine_idx] = (4.0 * residual[fine_idx] - residual_next_level[coarse_idx]) / 3.0;
    }

}

void GMGPolar::extrapolatedResidual(const int current_level, GPU_Vector<double>& residual,
                                    const GPU_Vector<double>& residual_next_level)
{
    const PolarGrid& fineGrid   = levels_[current_level].grid();
    const PolarGrid& coarseGrid = levels_[current_level + 1].grid();

    assert(residual.size() == fineGrid.numberOfNodes());
    assert(residual_next_level.size() == coarseGrid.numberOfNodes());

    dim3 threadsPerBlock(16, 16);
    dim3 numBlocks((fineGrid.nr() + 14 - 1) / 14,
                   (fineGrid.ntheta() + 14 - 1) / 14);

    applyExtrapolatedResidual_kernel<<<numBlocks, threadsPerBlock>>>(
        residual.data(), residual_next_level.data(), 
        levels_[current_level].device_grid(), levels_[current_level+1].device_grid()
    );

    cudaDeviceSynchronize();
}