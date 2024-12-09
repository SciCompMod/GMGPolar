
#include "../../include/Interpolation/interpolation.h"

#include "../../include/LinearAlgebra/Vector/gpu_vector_operations.h"

__global__ void applyFMG_Interpolation_kernel(PolarGrid* fine_grid, double* result, PolarGrid* coarse_grid, double* x) {
    int i_r = blockIdx.x * blockDim.x + threadIdx.x;
    int i_theta = blockIdx.y * blockDim.y + threadIdx.y;

    if (i_r >= fine_grid->nr() || i_theta >= fine_grid->ntheta()) return;

    int i_r_coarse = i_r >> 1;
    int i_theta_coarse = i_theta >> 1;

    /* Case 1: On the boundary */
    if(i_r == 0 || i_r == fine_grid->nr() - 1){
        if(i_theta & 1){
            double k0 = coarse_grid->angularSpacing(i_theta_coarse-1);
            double k1 = fine_grid->angularSpacing(i_theta-1);
            double k2 = fine_grid->angularSpacing(i_theta);
            double k3 = coarse_grid->angularSpacing(i_theta_coarse+1);
            
            double w_theta0 = - k1/k0 * k2/(k0+k1+k2) * (k2+k3)/(k0+k1+k2+k3);
            double w_theta1 = (k0+k1)/k0 * k2/(k1+k2) * (k2+k3)/(k1+k2+k3);
            double w_theta2 = (k0+k1)/(k0+k1+k2) * k1/(k1+k2) * (k2+k3)/k3;
            double w_theta3 = - (k0+k1)/(k0+k1+k2+k3) * k1/(k1+k2+k3) * k2/k3;
            
            result[fine_grid->index(i_r, i_theta)] = (
                w_theta0 * x[coarse_grid->index(i_r_coarse, i_theta_coarse-1)] + /* (0, -3) */ \
                w_theta1 * x[coarse_grid->index(i_r_coarse, i_theta_coarse  )] + /* (0, -1) */ \
                w_theta2 * x[coarse_grid->index(i_r_coarse, i_theta_coarse+1)] + /* (0, +1) */ \
                w_theta3 * x[coarse_grid->index(i_r_coarse, i_theta_coarse+2)] /* (0, +3) */ \
            );
        }
        else{
            result[fine_grid->index(i_r, i_theta)] =
                x[coarse_grid->index(i_r_coarse, i_theta_coarse)]; /* center */
        }
    }
    /* Case 2: Next to the boundary */
    else if(i_r == 1 || i_r == fine_grid->nr() - 2){
        if(i_theta & 1){
            double k0 = coarse_grid->angularSpacing(i_theta_coarse-1);
            double k1 = fine_grid->angularSpacing(i_theta-1);
            double k2 = fine_grid->angularSpacing(i_theta);
            double k3 = coarse_grid->angularSpacing(i_theta_coarse+1);
            
            double w_theta0 = - k1/k0 * k2/(k0+k1+k2) * (k2+k3)/(k0+k1+k2+k3);
            double w_theta1 = (k0+k1)/k0 * k2/(k1+k2) * (k2+k3)/(k1+k2+k3);
            double w_theta2 = (k0+k1)/(k0+k1+k2) * k1/(k1+k2) * (k2+k3)/k3;
            double w_theta3 = - (k0+k1)/(k0+k1+k2+k3) * k1/(k1+k2+k3) * k2/k3;
            
            double left_value = (
                w_theta0 * x[coarse_grid->index(i_r_coarse, i_theta_coarse-1)] + /* (-1, -3) */
                w_theta1 * x[coarse_grid->index(i_r_coarse, i_theta_coarse  )] + /* (-1, -1) */
                w_theta2 * x[coarse_grid->index(i_r_coarse, i_theta_coarse+1)] + /* (-1, +1) */
                w_theta3 * x[coarse_grid->index(i_r_coarse, i_theta_coarse+2)]   /* (-1, +3) */
            );
            double right_value = (
                w_theta0 * x[coarse_grid->index(i_r_coarse+1, i_theta_coarse-1)] + /* (+1, -3) */
                w_theta1 * x[coarse_grid->index(i_r_coarse+1, i_theta_coarse  )] + /* (+1, -1) */
                w_theta2 * x[coarse_grid->index(i_r_coarse+1, i_theta_coarse+1)] + /* (+1, +1) */
                w_theta3 * x[coarse_grid->index(i_r_coarse+1, i_theta_coarse+2)]   /* (+1, +3) */
            );

            double h1 = fine_grid->radialSpacing(i_r-1);
            double h2 = fine_grid->radialSpacing(i_r);
            result[fine_grid->index(i_r, i_theta)] = (h1 * left_value + h2 * right_value) / (h1 + h2);
        }
        else{
            double h1 = fine_grid->radialSpacing(i_r-1);
            double h2 = fine_grid->radialSpacing(i_r);
            result[fine_grid->index(i_r, i_theta)] = (
                h1 * x[coarse_grid->index(i_r_coarse, i_theta_coarse)] + /* left */
                h2 * x[coarse_grid->index(i_r_coarse+1, i_theta_coarse)] /* right */
            ) / (h1 + h2);
        }
    }
    else{
        /* Case 3: In the interior */
        if(i_r & 1){
            if(i_theta & 1){
                double k0 = coarse_grid->angularSpacing(i_theta_coarse-1);
                double k1 = fine_grid->angularSpacing(i_theta-1);
                double k2 = fine_grid->angularSpacing(i_theta);
                double k3 = coarse_grid->angularSpacing(i_theta_coarse+1);
                
                double w_theta0 = - k1/k0 * k2/(k0+k1+k2) * (k2+k3)/(k0+k1+k2+k3);
                double w_theta1 = (k0+k1)/k0 * k2/(k1+k2) * (k2+k3)/(k1+k2+k3);
                double w_theta2 = (k0+k1)/(k0+k1+k2) * k1/(k1+k2) * (k2+k3)/k3;
                double w_theta3 = - (k0+k1)/(k0+k1+k2+k3) * k1/(k1+k2+k3) * k2/k3;
                
                double outer_left_value = (
                    w_theta0 * x[coarse_grid->index(i_r_coarse-1, i_theta_coarse-1)] + /* (-3, -3) */
                    w_theta1 * x[coarse_grid->index(i_r_coarse-1, i_theta_coarse  )] + /* (-3, -1) */
                    w_theta2 * x[coarse_grid->index(i_r_coarse-1, i_theta_coarse+1)] + /* (-3, +1) */
                    w_theta3 * x[coarse_grid->index(i_r_coarse-1, i_theta_coarse+2)]   /* (-3, +3) */
                );
                double inner_left_value = (
                    w_theta0 * x[coarse_grid->index(i_r_coarse, i_theta_coarse-1)] + /* (-1, -3) */
                    w_theta1 * x[coarse_grid->index(i_r_coarse, i_theta_coarse  )] + /* (-1, -1) */
                    w_theta2 * x[coarse_grid->index(i_r_coarse, i_theta_coarse+1)] + /* (-1, +1) */
                    w_theta3 * x[coarse_grid->index(i_r_coarse, i_theta_coarse+2)]   /* (-1, +3) */
                );
                double inner_right_value = (
                    w_theta0 * x[coarse_grid->index(i_r_coarse+1, i_theta_coarse-1)] + /* (+1, -3) */
                    w_theta1 * x[coarse_grid->index(i_r_coarse+1, i_theta_coarse  )] + /* (+1, -1) */
                    w_theta2 * x[coarse_grid->index(i_r_coarse+1, i_theta_coarse+1)] + /* (+1, +1) */
                    w_theta3 * x[coarse_grid->index(i_r_coarse+1, i_theta_coarse+2)]   /* (+1, +3) */
                );
                double outer_right_value = (
                    w_theta0 * x[coarse_grid->index(i_r_coarse+2, i_theta_coarse-1)] + /* (+3, -3) */
                    w_theta1 * x[coarse_grid->index(i_r_coarse+2, i_theta_coarse  )] + /* (+3, -1) */
                    w_theta2 * x[coarse_grid->index(i_r_coarse+2, i_theta_coarse+1)] + /* (+3, +1) */
                    w_theta3 * x[coarse_grid->index(i_r_coarse+2, i_theta_coarse+2)]   /* (+3, +3) */
                );
                
                double h0 = coarse_grid->radialSpacing(i_r_coarse-1);
                double h1 = fine_grid->radialSpacing(i_r-1);
                double h2 = fine_grid->radialSpacing(i_r);
                double h3 = coarse_grid->radialSpacing(i_r_coarse+1);
                
                double w_r0 = - h1/h0 * h2/(h0+h1+h2) * (h2+h3)/(h0+h1+h2+h3);
                double w_r1 = (h0+h1)/h0 * h2/(h1+h2) * (h2+h3)/(h1+h2+h3);
                double w_r2 = (h0+h1)/(h0+h1+h2) * h1/(h1+h2) * (h2+h3)/h3;
                double w_r3 = - (h0+h1)/(h0+h1+h2+h3) * h1/(h1+h2+h3) * h2/h3;
                
                result[fine_grid->index(i_r, i_theta)] = (
                    w_r0 * outer_left_value +
                    w_r1 * inner_left_value +
                    w_r2 * inner_right_value +
                    w_r3 * outer_right_value
                );
            }
            else{
                double h0 = coarse_grid->radialSpacing(i_r_coarse-1);
                double h1 = fine_grid->radialSpacing(i_r-1);
                double h2 = fine_grid->radialSpacing(i_r);
                double h3 = coarse_grid->radialSpacing(i_r_coarse+1);
                
                double w_r0 = - h1/h0 * h2/(h0+h1+h2) * (h2+h3)/(h0+h1+h2+h3);
                double w_r1 = (h0+h1)/h0 * h2/(h1+h2) * (h2+h3)/(h1+h2+h3);
                double w_r2 = (h0+h1)/(h0+h1+h2) * h1/(h1+h2) * (h2+h3)/h3;
                double w_r3 = - (h0+h1)/(h0+h1+h2+h3) * h1/(h1+h2+h3) * h2/h3;
                
                result[fine_grid->index(i_r, i_theta)] = (
                    w_r0 * x[coarse_grid->index(i_r_coarse-1, i_theta_coarse)] + /* (-3, 0) */
                    w_r1 * x[coarse_grid->index(i_r_coarse, i_theta_coarse  )] + /* (-1, 0) */
                    w_r2 * x[coarse_grid->index(i_r_coarse+1, i_theta_coarse)] + /* (+1, 0) */
                    w_r3 * x[coarse_grid->index(i_r_coarse+2, i_theta_coarse)] /* (+3, 0) */
                );
            }
        }
        else{
            if(i_theta & 1){
                double k0 = coarse_grid->angularSpacing(i_theta_coarse-1);
                double k1 = fine_grid->angularSpacing(i_theta-1);
                double k2 = fine_grid->angularSpacing(i_theta);
                double k3 = coarse_grid->angularSpacing(i_theta_coarse+1);
                
                double w_theta0 = - k1/k0 * k2/(k0+k1+k2) * (k2+k3)/(k0+k1+k2+k3);
                double w_theta1 = (k0+k1)/k0 * k2/(k1+k2) * (k2+k3)/(k1+k2+k3);
                double w_theta2 = (k0+k1)/(k0+k1+k2) * k1/(k1+k2) * (k2+k3)/k3;
                double w_theta3 = - (k0+k1)/(k0+k1+k2+k3) * k1/(k1+k2+k3) * k2/k3;
                
                result[fine_grid->index(i_r, i_theta)]  = (
                    w_theta0 * x[coarse_grid->index(i_r_coarse, i_theta_coarse-1)] + /* (0, -3) */
                    w_theta1 * x[coarse_grid->index(i_r_coarse, i_theta_coarse  )] + /* (0, -1) */
                    w_theta2 * x[coarse_grid->index(i_r_coarse, i_theta_coarse+1)] + /* (0, +1) */
                    w_theta3 * x[coarse_grid->index(i_r_coarse, i_theta_coarse+2)]   /* (0, +3) */
                );
            }
            else{
                result[fine_grid->index(i_r, i_theta)] = x[coarse_grid->index(i_r_coarse, i_theta_coarse)]; /* center */
            }
        }
    }
}

/* Remark: This injection is not scaled. */
void Interpolation::applyFMGInterpolation(const Level& fromLevel, const Level& toLevel, GPU_Vector<double>& result, const GPU_Vector<double>& x) const
{
    assert(toLevel.level() == fromLevel.level() - 1);

    const PolarGrid& coarse_grid = fromLevel.grid();
    const PolarGrid& fine_grid = toLevel.grid();

    assert(x.size() == coarse_grid.numberOfNodes());
    assert(result.size() == fine_grid.numberOfNodes());

    dim3 threadsPerBlock(16, 16);
    dim3 numBlocks((fine_grid.nr() + threadsPerBlock.x - 1) / threadsPerBlock.x,
                   (fine_grid.ntheta() + threadsPerBlock.y - 1) / threadsPerBlock.y);

    applyFMG_Interpolation_kernel<<<numBlocks, threadsPerBlock>>>(
        toLevel.device_grid(), result.data(), fromLevel.device_grid(), x.data());

    cudaDeviceSynchronize();
}