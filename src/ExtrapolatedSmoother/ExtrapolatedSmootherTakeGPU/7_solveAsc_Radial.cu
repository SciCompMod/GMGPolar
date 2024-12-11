#include "../../../include/ExtrapolatedSmoother/ExtrapolatedSmootherTakeGPU/extrapolatedSmoother.h"


__global__ void solve_black_radial_diagonals_kernel(
    double* x, PolarGrid* grid, double* radial_main_diagonals)
{
    int i_r = grid->numberSmootherCircles() + (blockIdx.x * blockDim.x + threadIdx.x);
    int i_theta = 2 * (blockIdx.y * blockDim.y + threadIdx.y);

    if (grid->numberSmootherCircles() <= i_r && i_r < grid->nr() && i_theta < grid->ntheta()) {
        x[grid->index(i_r, i_theta)] /= radial_main_diagonals[i_theta * grid->lengthSmootherRadial() + i_r - grid->numberSmootherCircles()];
    }
}



void ExtrapolatedSmootherTakeGPU::solveAsc_BlackRadial(GPU_Vector<double>& x, const GPU_Vector<double>& rhs)
{
    const PolarGrid& grid = level_.grid();

    dim3 blockDim(16, 16);
    dim3 gridDim((grid.lengthSmootherRadial() + blockDim.x - 1) / blockDim.x, (grid.ntheta()/2 + blockDim.y - 1) / blockDim.y);
    solve_black_radial_diagonals_kernel<<<gridDim, blockDim>>>(
        x.data(),  level_.device_grid(),
        radial_main_diagonals_
    );
    cudaDeviceSynchronize();
}


void ExtrapolatedSmootherTakeGPU::solveAsc_WhiteRadial(GPU_Vector<double>& x, const GPU_Vector<double>& rhs)
{
    const PolarGrid& grid = level_.grid();

    const int start_white_radial_solver = 1;
    int start = start_white_radial_solver * grid.lengthSmootherRadial();
    int batch_count = grid.ntheta() / 2;
    int batch_stride = 2 * grid.lengthSmootherRadial();
    cusparseDgtsv2StridedBatch(
        sparse_handle_, 
        grid.lengthSmootherRadial(), 
        radial_lower_diagonals_ + start,
        radial_main_diagonals_ + start,
        radial_upper_diagonals_ + start, 
        x.data() + grid.numberCircularSmootherNodes() + start, 
        batch_count, batch_stride, pBuffer_
    );
}