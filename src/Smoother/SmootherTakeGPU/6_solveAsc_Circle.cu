#include "../../../include/Smoother/SmootherTakeGPU/smoother.h"

#include <chrono>

__global__ void combineCicleFactor_kernel(
    double* x, double* temp, double* factor,
    double* circle_lower_diagonals, double* circle_main_diagonals, double* circle_upper_diagonals,
    double* sherman_morrison_gammas,
    int start_circle_solver, PolarGrid* grid)
{
    int i_r = start_circle_solver + 2 * (blockIdx.x * blockDim.x + threadIdx.x);
    if(i_r < grid->numberSmootherCircles()){
        int matrix_index = i_r * grid->ntheta();
        double alpha = circle_lower_diagonals[matrix_index + 0];
        double gamma = sherman_morrison_gammas[i_r];
        factor[i_r] = (x[matrix_index + 0] + x[matrix_index + grid->ntheta()-1] * alpha / gamma) /
            (1.0 + temp[matrix_index + 0] + temp[matrix_index + grid->ntheta()-1] * alpha / gamma);
    }
}


__global__ void combineCicleSolutions_kernel(
    double* x, double* temp, double* factor,
    int start_circle_solver, PolarGrid* grid)
{
    int i_r = start_circle_solver + 2 * (blockIdx.x * blockDim.x + threadIdx.x);
    int i_theta = blockIdx.y * blockDim.y + threadIdx.y;

    if (i_r < grid->numberSmootherCircles()  && i_theta < grid->ntheta()) {
        int matrix_index = i_r * grid->ntheta();
        x[matrix_index + i_theta] -= factor[i_r] * temp[matrix_index + i_theta];
    }
}



void SmootherTakeGPU::solveAsc_BlackCircle(GPU_Vector<double>& x, const GPU_Vector<double>& rhs, GPU_Vector<double>& temp)
{
    const PolarGrid& grid = level_.grid();

    const int start_black_circles = (grid.numberSmootherCircles() % 2 == 0) ? 1 : 0;

    /* Can run in parallel if needed. */
    if(start_black_circles == 0){
        /* This approach is very slow. */
        /* I wasnt able to find working analysis + solve methods in cuSparse/cuSolver. */
        /* Thus we use Mumps instead. */

        // // Source: x, Destination: temp
        // cudaMemcpy(temp.data(), x.data(), grid.ntheta() * sizeof(double), cudaMemcpyDeviceToDevice);

        // int singularity = 0; 
        // double tol = 1e-10;  
        // int reorder = 0; 

        // int n = grid.ntheta();
        // int nnzA = DirBC_Interior_ ? grid.ntheta() : 4 * grid.ntheta();

        // cusolverSpDcsrlsvchol(
        //     solver_handle_,
        //     n, nnzA, descrA_,
        //     csrValA_, csrRowPtrA_, csrColIndA_,
        //     temp.data(),
        //     tol, reorder,
        //     x.data(),
        //     &singularity
        // );

        /* Using Mumps */
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
        cudaMemcpy(x.data(), host_solution.data(), grid.ntheta() * sizeof(double), cudaMemcpyHostToDevice);
    }

    int start_black_circle_solver = (grid.numberSmootherCircles() % 2 == 0) ? 1 : 2;
    int start = start_black_circle_solver * grid.ntheta();
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
    combineCicleFactor_kernel<<<factor_numBlocks, factor_blockSize>>>(
        x.data(), temp.data(), factor_,
        circle_lower_diagonals_, circle_main_diagonals_, circle_upper_diagonals_,
        sherman_morrison_gammas_,
        start_black_circle_solver,
        level_.device_grid()
    );
    cudaDeviceSynchronize();

    dim3 blockDim(16, 16);
    dim3 gridDim((batch_count + blockDim.x - 1) / blockDim.x, (grid.ntheta() + blockDim.y - 1) / blockDim.y);
    combineCicleSolutions_kernel<<<gridDim, blockDim>>>(
        x.data(), temp.data(), factor_,
        start_black_circle_solver, level_.device_grid()
    );
    cudaDeviceSynchronize();
}



void SmootherTakeGPU::solveAsc_WhiteCircle(GPU_Vector<double>& x, const GPU_Vector<double>& rhs, GPU_Vector<double>& temp)
{
    const PolarGrid& grid = level_.grid();

    const int start_white_circles = (grid.numberSmootherCircles() % 2 == 0) ? 0 : 1;

    /* Can run in parallel if needed. */
    if(start_white_circles == 0){
        /* This approach is very slow. */
        /* I wasnt able to find working analysis + solve methods in cuSparse/cuSolver. */
        /* Thus we use Mumps instead. */

        // // Source: x, Destination: temp
        // cudaMemcpy(temp.data(), x.data(), grid.ntheta() * sizeof(double), cudaMemcpyDeviceToDevice);

        // int singularity = 0; 
        // double tol = 1e-10;  
        // int reorder = 0; 

        // int n = grid.ntheta();
        // int nnzA = DirBC_Interior_ ? grid.ntheta() : 4 * grid.ntheta();

        // cusolverSpDcsrlsvchol(
        //     solver_handle_,
        //     n, nnzA, descrA_,
        //     csrValA_, csrRowPtrA_, csrColIndA_,
        //     temp.data(),
        //     tol, reorder,
        //     x.data(),
        //     &singularity
        // );

        /* Using Mumps */
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
        cudaMemcpy(x.data(), host_solution.data(), grid.ntheta() * sizeof(double), cudaMemcpyHostToDevice);
    }

    int start_white_circle_solver = (grid.numberSmootherCircles() % 2 == 0) ? 2 : 1;
    int start = start_white_circle_solver * grid.ntheta();
    int batch_count = (grid.numberSmootherCircles()-1) / 2;
    int batch_stride = 2 * grid.ntheta();
    /* cusparseDgtsv2StridedBatch could run in parallel if we use s different pBuffer_. */
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
    combineCicleFactor_kernel<<<factor_numBlocks, factor_blockSize>>>(
        x.data(), temp.data(), factor_,
        circle_lower_diagonals_, circle_main_diagonals_, circle_upper_diagonals_,
        sherman_morrison_gammas_,
        start_white_circle_solver,
        level_.device_grid()
    );
    cudaDeviceSynchronize();

    dim3 blockDim(16, 16);
    dim3 gridDim((batch_count + blockDim.x - 1) / blockDim.x, (grid.ntheta() + blockDim.y - 1) / blockDim.y);
    combineCicleSolutions_kernel<<<gridDim, blockDim>>>(
        x.data(), temp.data(), factor_,
        start_white_circle_solver, level_.device_grid()
    );
    cudaDeviceSynchronize();
}