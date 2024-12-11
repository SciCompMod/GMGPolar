#include "../../../include/ExtrapolatedSmoother/ExtrapolatedSmootherTakeGPU/extrapolatedSmoother.h"


__global__ void solve_circle_diagonals_kernel(
    double* x, PolarGrid* grid, double* circle_main_diagonals)
{
    int i_r = 2 + 2 * (blockIdx.x * blockDim.x + threadIdx.x);
    int i_theta = blockIdx.y * blockDim.y + threadIdx.y;

    if (i_r >= 2 && i_r < grid->numberSmootherCircles() && i_theta < grid->ntheta()) {
        x[grid->index(i_r, i_theta)] /= circle_main_diagonals[i_r * grid->ntheta() + i_theta];
    }
}

__global__ void extrapolated_computeCircleFactor_kernel(
    double* x, double* temp, double* factor,
    double* circle_lower_diagonals, double* circle_main_diagonals, double* circle_upper_diagonals,
    double* sherman_morrison_gammas, PolarGrid* grid)
{
    int i_r = 1 + 2 * (blockIdx.x * blockDim.x + threadIdx.x);
    if(i_r < grid->numberSmootherCircles()){
        int matrix_index = i_r * grid->ntheta();
        double alpha = circle_lower_diagonals[matrix_index + 0];
        double gamma = sherman_morrison_gammas[i_r];
        factor[i_r] = (x[matrix_index + 0] + x[matrix_index + grid->ntheta()-1] * alpha / gamma) /
            (1.0 + temp[matrix_index + 0] + temp[matrix_index + grid->ntheta()-1] * alpha / gamma);
    }
}


__global__ void extrapolated_combineCircleSolutions_kernel(
    double* x, double* temp, double* factor, PolarGrid* grid)
{
    int i_r = 1 + 2 * (blockIdx.x * blockDim.x + threadIdx.x);
    int i_theta = blockIdx.y * blockDim.y + threadIdx.y;

    if (i_r < grid->numberSmootherCircles()  && i_theta < grid->ntheta()) {
        int matrix_index = i_r * grid->ntheta();
        x[matrix_index + i_theta] -= factor[i_r] * temp[matrix_index + i_theta];
    }
}



void ExtrapolatedSmootherTakeGPU::solveCircleDiagonals(GPU_Vector<double>& x, const GPU_Vector<double>& rhs, GPU_Vector<double>& temp)
{
    const PolarGrid& grid = level_.grid();

    std::vector<double> host_solution(grid.ntheta());
    cudaMemcpy(host_solution.data(), x.data(), grid.ntheta() * sizeof(double), cudaMemcpyDeviceToHost);

    inner_boundary_mumps_solver_.job    = JOB_COMPUTE_SOLUTION;
    inner_boundary_mumps_solver_.nrhs   = 1; // single rhs vector
    inner_boundary_mumps_solver_.nz_rhs = grid.ntheta(); // non-zeros in rhs
    inner_boundary_mumps_solver_.rhs    = host_solution.data();
    inner_boundary_mumps_solver_.lrhs   = grid.ntheta(); // leading dimension of rhs
    dmumps_c(&inner_boundary_mumps_solver_);
    if (inner_boundary_mumps_solver_.info[0] != 0) {
        std::cerr << "Error solving the system: " << inner_boundary_mumps_solver_.info[0] << std::endl;
    }
    cudaMemcpy(x.data(),  host_solution.data(), grid.ntheta() * sizeof(double), cudaMemcpyHostToDevice);

    dim3 blockDim(16, 16);
    dim3 gridDim(((grid.numberSmootherCircles()-1)/2 + blockDim.x - 1) / blockDim.x, (grid.ntheta() + blockDim.y - 1) / blockDim.y);
    solve_circle_diagonals_kernel<<<gridDim, blockDim>>>(
        x.data(), level_.device_grid(), circle_main_diagonals_
    );
    cudaDeviceSynchronize();

}

void ExtrapolatedSmootherTakeGPU::solveCircleTridiagonals(GPU_Vector<double>& x, const GPU_Vector<double>& rhs, GPU_Vector<double>& temp)
{
    const PolarGrid& grid = level_.grid();

    int start = grid.ntheta();
    int batch_count = grid.numberSmootherCircles() / 2;
    int batch_stride = 2 * grid.ntheta();

    /* cusparseDgtsv2StridedBatch could run in parallel if we use different pBuffer_. */
    cusparseDgtsv2StridedBatch(
        sparse_handle_, 
        grid.ntheta(), 
        circle_lower_diagonals_ + start, 
        circle_main_diagonals_ + start, 
        circle_upper_diagonals_ + start, 
        x.data() + start, 
        batch_count, batch_stride, 
        pBuffer_
    );

    cusparseDgtsv2StridedBatch(
        sparse_handle_, 
        grid.ntheta(), 
        circle_lower_diagonals_ + start, 
        circle_main_diagonals_ + start, 
        circle_upper_diagonals_ + start, 
        temp.data() + start, 
        batch_count, batch_stride, 
        pBuffer_
    );

    int factor_blockSize = 256;
    int factor_numBlocks = (batch_count + factor_blockSize - 1) / factor_blockSize;
    extrapolated_computeCircleFactor_kernel<<<factor_numBlocks, factor_blockSize>>>(
        x.data(), temp.data(), factor_,
        circle_lower_diagonals_, circle_main_diagonals_, circle_upper_diagonals_,
        sherman_morrison_gammas_,
        level_.device_grid()
    );
    cudaDeviceSynchronize();

    dim3 blockDim(16, 16);
    dim3 gridDim((batch_count + blockDim.x - 1) / blockDim.x, (grid.ntheta() + blockDim.y - 1) / blockDim.y);
    extrapolated_combineCircleSolutions_kernel<<<gridDim, blockDim>>>(
        x.data(), temp.data(), factor_, level_.device_grid()
    );
    cudaDeviceSynchronize();
}





void ExtrapolatedSmootherTakeGPU::solveAsc_BlackCircle(GPU_Vector<double>& x, const GPU_Vector<double>& rhs, GPU_Vector<double>& temp)
{
    const PolarGrid& grid = level_.grid();

    const int start_black_circles = (grid.numberSmootherCircles() % 2 == 0) ? 1 : 0;
    if(start_black_circles == 0){
        solveCircleDiagonals(x, rhs, temp);
    }
    else{
        solveCircleTridiagonals(x, rhs, temp);
    }
}

void ExtrapolatedSmootherTakeGPU::solveAsc_WhiteCircle(GPU_Vector<double>& x, const GPU_Vector<double>& rhs, GPU_Vector<double>& temp)
{
    const PolarGrid& grid = level_.grid();

    const int start_white_circles = (grid.numberSmootherCircles() % 2 == 0) ? 0 : 1;
    if(start_white_circles == 0){
        solveCircleDiagonals(x, rhs, temp);
    }
    else{
        solveCircleTridiagonals(x, rhs, temp);
    }
}






// __global__ void extrapolated_computeCicleFactor_kernel(
//     double* x, double* temp, double* factor,
//     double* circle_lower_diagonals, double* circle_main_diagonals, double* circle_upper_diagonals,
//     double* sherman_morrison_gammas,
//     int start_circle_solver, PolarGrid* grid)
// {
//     int i_r = 1 + 2 * (blockIdx.x * blockDim.x + threadIdx.x);
//     if(i_r < grid->numberSmootherCircles()){
//         int matrix_index = i_r/2 * grid->ntheta();
//         double alpha = circle_lower_diagonals[matrix_index + 0];
//         double gamma = sherman_morrison_gammas[i_r/2];
//         factor[i_r/2] = (x[matrix_index + 0] + x[matrix_index + grid->ntheta()-1] * alpha / gamma) /
//             (1.0 + temp[matrix_index + 0] + temp[matrix_index + grid->ntheta()-1] * alpha / gamma);
//     }
// }


// __global__ void extrapolated_combineCicleSolutions_kernel(
//     double* x, double* temp, double* factor,
//     int start_circle_solver, PolarGrid* grid)
// {
//     int i_r = 1 + 2 * (blockIdx.x * blockDim.x + threadIdx.x);
//     int i_theta = blockIdx.y * blockDim.y + threadIdx.y;

//     if (i_r < grid->numberSmootherCircles()  && i_theta < grid->ntheta()) {
//         int matrix_index = i_r/2 * grid->ntheta();
//         x[matrix_index + i_theta] -= factor[i_r/2] * temp[matrix_index + i_theta];
//     }
// }


// __global__ void solve_circle_diagonals_kernel(
//     double* x, PolarGrid* grid, double* circle_diagonals)
// {
//     // int i_r = 2 + 2 * (blockIdx.x * blockDim.x + threadIdx.x);
//     // int i_theta = blockIdx.y * blockDim.y + threadIdx.y;

//     // if (i_r >= 2 && i_r < grid->numberSmootherCircles() && i_theta < grid->ntheta()) {
//     //     x[grid->index(i_r, i_theta)] /= circle_diagonals[i_r/2 * grid->ntheta() + i_theta];
//     // }
// }


// void ExtrapolatedSmootherTakeGPU::solveAsc_BlackCircle(GPU_Vector<double>& x, const GPU_Vector<double>& rhs, GPU_Vector<double>& temp)
// {
//     const PolarGrid& grid = level_.grid();

//     const int start_black_circles = (grid.numberSmootherCircles() % 2 == 0) ? 1 : 0;

//     /* Can run in parallel if needed. */
//     if(start_black_circles == 0){
//         // std::vector<double> host_solution(grid.ntheta());
//         // cudaMemcpy(host_solution.data(), x.data(), grid.ntheta() * sizeof(double), cudaMemcpyDeviceToHost);

//         // inner_boundary_mumps_solver_.job    = JOB_COMPUTE_SOLUTION;
//         // inner_boundary_mumps_solver_.nrhs   = 1; // single rhs vector
//         // inner_boundary_mumps_solver_.nz_rhs = grid.ntheta(); // non-zeros in rhs
//         // inner_boundary_mumps_solver_.rhs    = host_solution.data();
//         // inner_boundary_mumps_solver_.lrhs   = grid.ntheta(); // leading dimension of rhs
//         // dmumps_c(&inner_boundary_mumps_solver_);
//         // if (inner_boundary_mumps_solver_.info[0] != 0) {
//         //     std::cerr << "Error solving the system: " << inner_boundary_mumps_solver_.info[0] << std::endl;
//         // }
//         // cudaMemcpy(x.data(), host_solution.data(), grid.ntheta() * sizeof(double), cudaMemcpyHostToDevice);
//     }

//     if(start_black_circles == 1){
//         int start_black_circle_solver = (grid.numberSmootherCircles() % 2 == 0) ? 1 : 2;
//         int start = start_black_circle_solver * grid.ntheta();
//         int batch_count = (grid.numberSmootherCircles()-1) / 2;
//         int batch_stride = grid.ntheta();
//         /* cusparseDgtsv2StridedBatch could run in parallel if we use different pBuffer_. */
//         cusparseDgtsv2StridedBatch(
//             sparse_handle_, 
//             grid.ntheta(), 
//             circle_lower_diagonals_, 
//             circle_main_diagonals_, 
//             circle_upper_diagonals_, 
//             x.data() + start, 
//             batch_count, batch_stride, 
//             pBuffer_
//         );

//         cusparseDgtsv2StridedBatch(
//             sparse_handle_, 
//             grid.ntheta(), 
//             circle_lower_diagonals_, 
//             circle_main_diagonals_, 
//             circle_upper_diagonals_, 
//             temp.data() + start, 
//             batch_count, batch_stride, 
//             pBuffer_
//         );

//         int factor_blockSize = 256;
//         int factor_numBlocks = (batch_count + factor_blockSize - 1) / factor_blockSize;
//         extrapolated_computeCicleFactor_kernel<<<factor_numBlocks, factor_blockSize>>>(
//             x.data(), temp.data(), factor_,
//             circle_lower_diagonals_, circle_main_diagonals_, circle_upper_diagonals_,
//             sherman_morrison_gammas_,
//             start_black_circle_solver,
//             level_.device_grid()
//         );
//         cudaDeviceSynchronize();

//         dim3 blockDim(16, 16);
//         dim3 gridDim((batch_count + blockDim.x - 1) / blockDim.x, (grid.ntheta() + blockDim.y - 1) / blockDim.y);
//         extrapolated_combineCicleSolutions_kernel<<<gridDim, blockDim>>>(
//             x.data(), temp.data(), factor_,
//             start_black_circle_solver, level_.device_grid()
//         );
//         cudaDeviceSynchronize();
//     }
//     else{
//         dim3 blockDim(16, 16);
//         dim3 gridDim((grid.numberSmootherCircles()/2 + blockDim.x - 1) / blockDim.x, (grid.ntheta() + blockDim.y - 1) / blockDim.y);
//         solve_circle_diagonals_kernel<<<gridDim, blockDim>>>(
//             x.data(), level_.device_grid(), circle_diagonals_
//         );
//         cudaDeviceSynchronize();
//     }
// }





// void ExtrapolatedSmootherTakeGPU::solveAsc_WhiteCircle(GPU_Vector<double>& x, const GPU_Vector<double>& rhs, GPU_Vector<double>& temp)
// {
//     const PolarGrid& grid = level_.grid();

//     const int start_white_circles = (grid.numberSmootherCircles() % 2 == 0) ? 0 : 1;

//     /* Can run in parallel if needed. */
//     if(start_white_circles == 0){
//         // std::vector<double> host_solution(grid.ntheta());
//         // cudaMemcpy(host_solution.data(), x.data(), grid.ntheta() * sizeof(double), cudaMemcpyDeviceToHost);

//         // inner_boundary_mumps_solver_.job    = JOB_COMPUTE_SOLUTION;
//         // inner_boundary_mumps_solver_.nrhs   = 1; // single rhs vector
//         // inner_boundary_mumps_solver_.nz_rhs = grid.ntheta(); // non-zeros in rhs
//         // inner_boundary_mumps_solver_.rhs    = host_solution.data();
//         // inner_boundary_mumps_solver_.lrhs   = grid.ntheta(); // leading dimension of rhs
//         // dmumps_c(&inner_boundary_mumps_solver_);
//         // if (inner_boundary_mumps_solver_.info[0] != 0) {
//         //     std::cerr << "Error solving the system: " << inner_boundary_mumps_solver_.info[0] << std::endl;
//         // }
//         // cudaMemcpy(x.data(), host_solution.data(), grid.ntheta() * sizeof(double), cudaMemcpyHostToDevice);
//     }

//     if(start_white_circles == 1){
//         int start_white_circle_solver = (grid.numberSmootherCircles() % 2 == 0) ? 2 : 1;
//         int start = start_white_circle_solver * grid.ntheta();
//         int batch_count = grid.numberSmootherCircles() / 2;
//         int batch_stride = grid.ntheta();
//         /* cusparseDgtsv2StridedBatch could run in parallel if we use different pBuffer_. */
//         cusparseDgtsv2StridedBatch(
//             sparse_handle_, 
//             grid.ntheta(), 
//             circle_lower_diagonals_, 
//             circle_main_diagonals_, 
//             circle_upper_diagonals_, 
//             x.data() + start, 
//             batch_count, batch_stride, 
//             pBuffer_
//         );

//         cusparseDgtsv2StridedBatch(
//             sparse_handle_, 
//             grid.ntheta(), 
//             circle_lower_diagonals_, 
//             circle_main_diagonals_, 
//             circle_upper_diagonals_, 
//             temp.data() + start, 
//             batch_count, batch_stride, 
//             pBuffer_
//         );

//         int factor_blockSize = 256;
//         int factor_numBlocks = (batch_count + factor_blockSize - 1) / factor_blockSize;
//         extrapolated_computeCicleFactor_kernel<<<factor_numBlocks, factor_blockSize>>>(
//             x.data(), temp.data(), factor_,
//             circle_lower_diagonals_, circle_main_diagonals_, circle_upper_diagonals_,
//             sherman_morrison_gammas_,
//             start_white_circle_solver,
//             level_.device_grid()
//         );
//         cudaDeviceSynchronize();

//         dim3 blockDim(16, 16);
//         dim3 gridDim((batch_count + blockDim.x - 1) / blockDim.x, (grid.ntheta() + blockDim.y - 1) / blockDim.y);
//         extrapolated_combineCicleSolutions_kernel<<<gridDim, blockDim>>>(
//             x.data(), temp.data(), factor_,
//             start_white_circle_solver, level_.device_grid()
//         );
//         cudaDeviceSynchronize();
//     }
//     else{
//         dim3 blockDim(16, 16);
//         dim3 gridDim((grid.numberSmootherCircles()/2 + blockDim.x - 1) / blockDim.x, (grid.ntheta() + blockDim.y - 1) / blockDim.y);
//         solve_circle_diagonals_kernel<<<gridDim, blockDim>>>(
//             x.data(), level_.device_grid(), circle_diagonals_
//         );
//         cudaDeviceSynchronize();
//     }
// }