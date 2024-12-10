#include "../../../include/Smoother/SmootherTakeGPU/smoother.h"

__global__ void adjustAscCircle_ShermanMorrison_kernel(
    double* circle_lower_diagonals, double* circle_main_diagonals, double* circle_upper_diagonals,
    double* sherman_morrison_gammas,
    PolarGrid* grid)
{
    int i_r = blockIdx.x * blockDim.x + threadIdx.x;
    if(i_r < 0 || i_r >= grid->numberSmootherCircles()) return;

    int circle_m = grid->ntheta();

    if(i_r == 0){
        sherman_morrison_gammas[i_r] = 0.0; // Unused
    } 
    else if(i_r > 0 && i_r < grid->numberSmootherCircles())
    {
        /* gamma_ = -main_diagonal(0); */
        double gamma = -circle_main_diagonals[i_r * circle_m + 0];
        /* main_diagonal(0) -= gamma_; */
        circle_main_diagonals[i_r * circle_m + 0] -= gamma;
        /* main_diagonal(matrix_dimension_ - 1) -= cyclic_corner_element() * cyclic_corner_element() / gamma_; */
        double alpha = circle_lower_diagonals[i_r * circle_m + 0];
        double beta = circle_upper_diagonals[i_r * circle_m + circle_m-1];
        circle_main_diagonals[i_r * circle_m + circle_m-1] -= alpha * beta / gamma;
        /* We will need gamma later again. */
        sherman_morrison_gammas[i_r] = gamma;
    }
}

void SmootherTakeGPU::adjustAscCircle_ShermanMorrison()
{
    const PolarGrid& grid = level_.grid();

    int blockSize = 256; 
    int numBlocks = (grid.numberSmootherCircles() + blockSize - 1) / blockSize;

    adjustAscCircle_ShermanMorrison_kernel<<<numBlocks, blockSize>>>(
        circle_lower_diagonals_, circle_main_diagonals_, circle_upper_diagonals_,
        sherman_morrison_gammas_,
        level_.device_grid()
    );
    cudaDeviceSynchronize();
}