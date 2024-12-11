#include "../../include/GMGPolar/gmgpolar.h"

#include "../../include/LinearAlgebra/Vector/vector_operations.h"
#include "../../include/LinearAlgebra/Vector/gpu_vector_operations.h"

std::pair<double, double> GMGPolar::computeExactError(Level& level, const Vector<double>& solution,
                                                      Vector<double>& error)
{
    assert(exact_solution_ != nullptr);

    const PolarGrid& grid        = level.grid();
    const LevelCache& levelCache = level.levelCache();
    const auto& sin_theta_cache  = levelCache.sin_theta();
    const auto& cos_theta_cache  = levelCache.cos_theta();

    assert(solution.size() == error.size());
    assert(solution.size() == grid.numberOfNodes());

#pragma omp parallel
    {
#pragma omp for nowait
        for (int i_r = 0; i_r < grid.numberSmootherCircles(); i_r++) {
            double r = grid.radius(i_r);
            for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
                double theta                    = grid.theta(i_theta);
                double sin_theta                = sin_theta_cache[i_theta];
                double cos_theta                = cos_theta_cache[i_theta];
                error[grid.index(i_r, i_theta)] = exact_solution_->exact_solution(r, theta, sin_theta, cos_theta) -
                                                  solution[grid.index(i_r, i_theta)];
            }
        }
#pragma omp for nowait
        for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
            double theta     = grid.theta(i_theta);
            double sin_theta = sin_theta_cache[i_theta];
            double cos_theta = cos_theta_cache[i_theta];
            for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++) {
                double r                        = grid.radius(i_r);
                error[grid.index(i_r, i_theta)] = exact_solution_->exact_solution(r, theta, sin_theta, cos_theta) -
                                                  solution[grid.index(i_r, i_theta)];
            }
        }
    }

    double weighted_euclidean_error = l2_norm(error) / sqrt(grid.numberOfNodes());
    double infinity_error           = infinity_norm(error);

    return std::make_pair(weighted_euclidean_error, infinity_error);
}


__global__ void computeExactError_kernel(
    double* solution, double* error,
    PolarGrid* grid, ExactSolution* exact_solution,
    double* sin_theta_cache, double* cos_theta_cache) 
{

    int i_r = blockIdx.x * 14 + threadIdx.x - 1;
    int i_theta = blockIdx.y * 14 + threadIdx.y - 1;

    if (i_r < 0 || i_r >= grid->nr() || i_theta < 0 || i_theta >= grid->ntheta()) return;

    double r = grid->radius(i_r);
    double theta = grid->theta(i_theta);

    double sin_theta = sin_theta_cache[i_theta];
    double cos_theta = cos_theta_cache[i_theta];

    int index = grid->index(i_r, i_theta);

    error[index] = exact_solution->exact_solution(r, theta, sin_theta, cos_theta) - solution[index];
}


std::pair<double, double> GMGPolar::computeExactError(Level& level, const GPU_Vector<double>& solution,
                                                      GPU_Vector<double>& error)
{
    assert(exact_solution_ != nullptr);

    const PolarGrid& grid = level.grid();

    assert(solution.size() == error.size());
    assert(solution.size() == grid.numberOfNodes());

    const GPU_Vector<double>& sin_theta_cache = level.levelCache().GPU_sin_theta();
    const GPU_Vector<double>& cos_theta_cache = level.levelCache().GPU_cos_theta();

    ExactSolution* device_exact_solution;
    cudaMalloc(&device_exact_solution, sizeof(ExactSolution));
    cudaMemcpy(device_exact_solution, exact_solution_.get(), sizeof(ExactSolution), cudaMemcpyHostToDevice);

    dim3 threadsPerBlock(16, 16);
    dim3 numBlocks((grid.nr() + 14 - 1) / 14,
                   (grid.ntheta() + 14 - 1) / 14);

    computeExactError_kernel<<<numBlocks, threadsPerBlock>>>(
        solution.data(), error.data(), 
        level.device_grid(), device_exact_solution,
        sin_theta_cache.data(), cos_theta_cache.data()
    );

    cudaDeviceSynchronize();

    cudaFree(device_exact_solution);

    double weighted_euclidean_error = l2_norm(error) / sqrt(grid.numberOfNodes());
    double infinity_error           = infinity_norm(error);

    return std::make_pair(weighted_euclidean_error, infinity_error);

}