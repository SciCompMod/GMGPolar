// #include "../../include/Residual/ResidualTakeGPU/residual.h"

// __global__ void applyResidual_kernel(
//     double* result, double* rhs, double* x,
//     PolarGrid* grid, bool DirBC_Interior,
//     DomainGeometry* domain_geometry,
//     double* coeff_alpha_cache, double* coeff_beta_cache,
//     double* sin_theta_cache, double* cos_theta_cache) 
// {
//     int i_r = blockIdx.x * blockDim.x + threadIdx.x;
//     int i_theta = blockIdx.y * blockDim.y + threadIdx.y;

//     if (i_r >= grid->nr() || i_theta >= grid->ntheta()) return;

//     int local_i_r = threadIdx.x;
//     int local_i_theta = threadIdx.y;
//     int local_index = threadIdx.x + threadIdx.y * blockDim.x;

//     int local_i_r_corner = blockIdx.x * blockDim.x;
//     int local_i_theta_corner = blockIdx.y * blockDim.y;

//     __shared__ double s_detDF[18][18];
//     __shared__ double s_arr[18][18];
//     __shared__ double s_att[18][18];
//     __shared__ double s_art[18][18];
//     __shared__ double s_x[18][18];

//     int partition_size = (18 * 18) / (16 * 16);
//     int remainder = (18 * 18) % (16 * 16);

//     int shared_global_start = local_index * partition_size + min(local_index, remainder);
//     int shared_global_end = shared_global_start + partition_size + (local_index < remainder ? 1 : 0);

//     for (int shared_index = shared_global_start; shared_index < shared_global_end; shared_index++) {
//         int shared_i_r = shared_index % 18;
//         int shared_i_theta = shared_index / 18;

//         int i_r = local_i_r_corner + shared_i_r - 1;
//         int i_theta = local_i_theta_corner + shared_i_theta - 1;

//         if (i_r >= 0 && i_r < grid->nr() && i_theta >= 0 && i_theta < grid->ntheta()) {
//             double r = grid->radius(i_r);
//             double theta = grid->theta(i_theta);

//             double sin_theta = sin_theta_cache[i_theta];
//             double cos_theta = cos_theta_cache[i_theta];
            
//             double Jrr = domain_geometry->dFx_dr(r, theta, sin_theta, cos_theta);
//             double Jtr = domain_geometry->dFy_dr(r, theta, sin_theta, cos_theta);
//             double Jrt = domain_geometry->dFx_dt(r, theta, sin_theta, cos_theta);
//             double Jtt = domain_geometry->dFy_dt(r, theta, sin_theta, cos_theta);

//             double coeff_alpha = coeff_alpha_cache[i_r];

//             double detDF = Jrr * Jtt - Jrt * Jtr;
//             double arr = 0.5 * (Jtt * Jtt + Jrt * Jrt) * coeff_alpha / fabs(detDF);
//             double att = 0.5 * (Jtr * Jtr + Jrr * Jrr) * coeff_alpha / fabs(detDF);
//             double art = (- Jtt * Jtr - Jrt * Jrr) * coeff_alpha / fabs(detDF);

//             s_detDF[shared_i_r][shared_i_theta] = detDF;
//             s_arr[shared_i_r][shared_i_theta] = arr;
//             s_att[shared_i_r][shared_i_theta] = att;
//             s_art[shared_i_r][shared_i_theta] = art;

//             s_x[shared_i_r][shared_i_theta] = x[grid->index(i_r, i_theta)];
//         }
//     }

//     __syncthreads();

//     bool isOnInnerBoundary = (i_r == 0);
//     bool isOnOuterBoundary = (i_r == grid->nr() - 1);

//     double h1 = DirBC_Interior ? 
//         ((!isOnInnerBoundary) ? grid->radialSpacing(i_r - 1) : 0.0) :
//         ((!isOnInnerBoundary) ? grid->radialSpacing(i_r - 1) : 2.0 * grid->radius(0));
//     double h2 = (!isOnOuterBoundary) ? grid->radialSpacing(i_r - 1) : 0.0;
//     double k1 = grid->angularSpacing(i_theta - 1);                                                          
//     double k2 = grid->angularSpacing(i_theta);

//     int center_index = grid->index(i_r, i_theta);
//     int s_i_r = local_i_r + 1;
//     int s_i_theta = local_i_theta + 1;

//     if (!isOnInnerBoundary && !isOnOuterBoundary) {    

//         double coeff1 = 0.5 * (k1 + k2) / h1;
//         double coeff2 = 0.5 * (k1 + k2) / h2;
//         double coeff3 = 0.5 * (h1 + h2) / k1;
//         double coeff4 = 0.5 * (h1 + h2) / k2;

//         result[center_index] = rhs[center_index] - 

//             (0.25 * (h1 + h2) * (k1 + k2) * coeff_beta_cache[i_r] * fabs(s_detDF[s_i_r][s_i_theta]) * s_x[s_i_r][s_i_theta] /* beta_{i,j} */ 
                                                                                                               
//             - coeff1 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r-1][s_i_theta]) * (s_x[s_i_r-1][s_i_theta] - s_x[s_i_r][s_i_theta]) /* Left - Center: (Left) */          
//             - coeff2 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r+1][s_i_theta]) * (s_x[s_i_r+1][s_i_theta] - s_x[s_i_r][s_i_theta]) /* Right - Center: (Right) */      
//             - coeff3 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta-1]) * (s_x[s_i_r][s_i_theta-1] - s_x[s_i_r][s_i_theta]) /* Bottom - Center: (Bottom) */  
//             - coeff4 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta+1]) * (s_x[s_i_r][s_i_theta+1] - s_x[s_i_r][s_i_theta]) /* Top - Center: (Top) */              
                                                                                                                   
//             - 0.25 * (s_art[s_i_r-1][s_i_theta] + s_art[s_i_r][s_i_theta-1]) * s_x[s_i_r-1][s_i_theta-1] /* Bottom Left */                             
//             + 0.25 * (s_art[s_i_r+1][s_i_theta] + s_art[s_i_r][s_i_theta-1]) * s_x[s_i_r+1][s_i_theta-1] /* Bottom Right */                          
//             + 0.25 * (s_art[s_i_r-1][s_i_theta] + s_art[s_i_r][s_i_theta+1]) * s_x[s_i_r-1][s_i_theta+1] /* Top Left */                                      
//             - 0.25 * (s_art[s_i_r+1][s_i_theta] + s_art[s_i_r][s_i_theta+1]) * s_x[s_i_r+1][s_i_theta+1] /* Top Right */                                   
//         );
//     }
//     else if((isOnInnerBoundary && !DirBC_Interior)){

//         double coeff1 = 0.5 * (k1 + k2) / h1;
//         double coeff2 = 0.5 * (k1 + k2) / h2;
//         double coeff3 = 0.5 * (h1 + h2) / k1;
//         double coeff4 = 0.5 * (h1 + h2) / k2;

//         result[center_index] = rhs[center_index] - 

//             (0.25 * (h1 + h2) * (k1 + k2) * coeff_beta_cache[i_r] * fabs(s_detDF[s_i_r][s_i_theta]) * s_x[s_i_r][s_i_theta] /* beta_{i,j} */ 
                                                                                                               
//             - coeff1 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r-1][s_i_theta]) * (s_x[s_i_r-1][s_i_theta] - s_x[s_i_r][s_i_theta]) /* Left - Center: (Left) */          
//             - coeff2 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r+1][s_i_theta]) * (s_x[s_i_r+1][s_i_theta] - s_x[s_i_r][s_i_theta]) /* Right - Center: (Right) */      
//             - coeff3 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta-1]) * (s_x[s_i_r][s_i_theta-1] - s_x[s_i_r][s_i_theta]) /* Bottom - Center: (Bottom) */  
//             - coeff4 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta+1]) * (s_x[s_i_r][s_i_theta+1] - s_x[s_i_r][s_i_theta]) /* Top - Center: (Top) */              
                                                                                                                                              
//             + 0.25 * (s_art[s_i_r+1][s_i_theta] + s_art[s_i_r][s_i_theta-1]) * s_x[s_i_r+1][s_i_theta-1] /* Bottom Right */                                                              
//             - 0.25 * (s_art[s_i_r+1][s_i_theta] + s_art[s_i_r][s_i_theta+1]) * s_x[s_i_r+1][s_i_theta+1] /* Top Right */                                   
//         );

//     }
//     else if((isOnInnerBoundary && DirBC_Interior) || isOnOuterBoundary){
//         result[center_index] = rhs[center_index] - s_x[s_i_r][s_i_theta];
//         // result[center_index] = rhs[center_index] - x[center_index]; 
//     }
// }


// /* ------------------ */
// /* result = rhs - A*x */
// void ResidualTakeGPU::computeResidual(GPU_Vector<double>& result, const GPU_Vector<double>& rhs, const GPU_Vector<double>& x) const {

//     const PolarGrid& grid = level_.grid();

//     assert(result.size() == grid.numberOfNodes());
//     assert(rhs.size() == grid.numberOfNodes());
//     assert(x.size() == grid.numberOfNodes());

//     const GPU_Vector<double>& sin_theta_cache = level_.levelCache().GPU_sin_theta();
//     const GPU_Vector<double>& cos_theta_cache = level_.levelCache().GPU_cos_theta();

//     const GPU_Vector<double>& coeff_alpha_cache = level_.levelCache().GPU_coeff_alpha();
//     const GPU_Vector<double>& coeff_beta_cache = level_.levelCache().GPU_coeff_beta();

//     DomainGeometry* device_domain_geometry;
//     cudaMalloc(&device_domain_geometry, sizeof(DomainGeometry));
//     cudaMemcpy(device_domain_geometry, &domain_geometry_, sizeof(DomainGeometry), cudaMemcpyHostToDevice);

//     /* We use precomputed DensityProfileCoefficients values. */
//     // DensityProfileCoefficients* device_density_profile;
//     // cudaMalloc(&device_density_profile, sizeof(DensityProfileCoefficients));
//     // cudaMemcpy(device_density_profile, &density_profile_coefficients_, sizeof(DensityProfileCoefficients), cudaMemcpyHostToDevice);
    
//     dim3 threadsPerBlock(16, 16);
//     dim3 numBlocks((grid.nr() + threadsPerBlock.x - 1) / threadsPerBlock.x,
//                    (grid.ntheta() + threadsPerBlock.y - 1) / threadsPerBlock.y);

//     applyResidual_kernel<<<numBlocks, threadsPerBlock>>>(
//         result.data(), rhs.data(), x.data(), 
//         level_.device_grid(), DirBC_Interior_,
//         device_domain_geometry, 
//         coeff_alpha_cache.data(), coeff_beta_cache.data(), 
//         sin_theta_cache.data(), cos_theta_cache.data()
//     );

//     cudaDeviceSynchronize();

//     /* We use precomputed DensityProfileCoefficients values. */
//     cudaFree(device_domain_geometry);
//     // cudaFree(device_density_profile);
// }
























#include "../../include/Residual/ResidualTakeGPU/residual.h"

__global__ void applyResidual_kernel(
    double* result, double* rhs, double* x,
    PolarGrid* grid, bool DirBC_Interior,
    DomainGeometry* domain_geometry,
    double* coeff_alpha_cache, double* coeff_beta_cache,
    double* sin_theta_cache, double* cos_theta_cache) 
{

    int i_r = blockIdx.x * 14 + threadIdx.x - 1;
    int i_theta = blockIdx.y * 14 + threadIdx.y - 1;

    if(i_r == -1 && !DirBC_Interior){
        i_r = 0;
        i_theta += grid->ntheta() / 2;
    }
    i_theta = grid->wrapThetaIndex(i_theta);

    if (i_r < 0 || i_r >= grid->nr() || i_theta < 0 || i_theta >= grid->ntheta()) return;

    __shared__ double s_detDF[16][16];
    __shared__ double s_arr[16][16];
    __shared__ double s_att[16][16];
    __shared__ double s_art[16][16];
    __shared__ double s_x[16][16];

    int s_i_r = threadIdx.x;
    int s_i_theta = threadIdx.y;

    int center_index = grid->index(i_r, i_theta);
    s_x[s_i_r][s_i_theta] = x[center_index];

    double r = grid->radius(i_r);
    double theta = grid->theta(i_theta);

    double sin_theta = sin_theta_cache[i_theta];
    double cos_theta = cos_theta_cache[i_theta];
    
    double Jrr = domain_geometry->dFx_dr(r, theta, sin_theta, cos_theta);
    double Jtr = domain_geometry->dFy_dr(r, theta, sin_theta, cos_theta);
    double Jrt = domain_geometry->dFx_dt(r, theta, sin_theta, cos_theta);
    double Jtt = domain_geometry->dFy_dt(r, theta, sin_theta, cos_theta);

    double coeff_alpha = coeff_alpha_cache[i_r];

    double detDF = Jrr * Jtt - Jrt * Jtr;
    double arr = 0.5 * (Jtt * Jtt + Jrt * Jrt) * coeff_alpha / fabs(detDF);
    double att = 0.5 * (Jtr * Jtr + Jrr * Jrr) * coeff_alpha / fabs(detDF);
    double art = (- Jtt * Jtr - Jrt * Jrr) * coeff_alpha / fabs(detDF);

    s_detDF[s_i_r][s_i_theta] = detDF;
    s_arr[s_i_r][s_i_theta] = arr;
    s_att[s_i_r][s_i_theta] = att;
    s_art[s_i_r][s_i_theta] = art;

    __syncthreads();

    if (i_r < 0 || i_r >= grid->nr() || i_theta < 0 || i_theta >= grid->ntheta()) return;
    if (s_i_r == 0 || s_i_r == 15 || s_i_theta == 0 || s_i_theta == 15) return;

    bool isOnInnerBoundary = (i_r == 0);
    bool isOnOuterBoundary = (i_r == grid->nr() - 1);

    double h1 = DirBC_Interior ? 
        ((!isOnInnerBoundary) ? grid->radialSpacing(i_r - 1) : 0.0) :
        ((!isOnInnerBoundary) ? grid->radialSpacing(i_r - 1) : 2.0 * grid->radius(0));
    double h2 = (!isOnOuterBoundary) ? grid->radialSpacing(i_r) : 0.0;
    double k1 = grid->angularSpacing(i_theta - 1);                                                          
    double k2 = grid->angularSpacing(i_theta);

    if (!isOnInnerBoundary && !isOnOuterBoundary) {    

        double coeff1 = 0.5 * (k1 + k2) / h1;
        double coeff2 = 0.5 * (k1 + k2) / h2;
        double coeff3 = 0.5 * (h1 + h2) / k1;
        double coeff4 = 0.5 * (h1 + h2) / k2;

        result[center_index] = rhs[center_index] - 

            (0.25 * (h1 + h2) * (k1 + k2) * coeff_beta_cache[i_r] * fabs(s_detDF[s_i_r][s_i_theta]) * s_x[s_i_r][s_i_theta] /* beta_{i,j} */ 
                                                                                                               
            - coeff1 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r-1][s_i_theta]) * (s_x[s_i_r-1][s_i_theta] - s_x[s_i_r][s_i_theta]) /* Left - Center: (Left) */          
            - coeff2 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r+1][s_i_theta]) * (s_x[s_i_r+1][s_i_theta] - s_x[s_i_r][s_i_theta]) /* Right - Center: (Right) */      
            - coeff3 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta-1]) * (s_x[s_i_r][s_i_theta-1] - s_x[s_i_r][s_i_theta]) /* Bottom - Center: (Bottom) */  
            - coeff4 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta+1]) * (s_x[s_i_r][s_i_theta+1] - s_x[s_i_r][s_i_theta]) /* Top - Center: (Top) */              
                                                                                                                   
            - 0.25 * (s_art[s_i_r-1][s_i_theta] + s_art[s_i_r][s_i_theta-1]) * s_x[s_i_r-1][s_i_theta-1] /* Bottom Left */                             
            + 0.25 * (s_art[s_i_r+1][s_i_theta] + s_art[s_i_r][s_i_theta-1]) * s_x[s_i_r+1][s_i_theta-1] /* Bottom Right */                          
            + 0.25 * (s_art[s_i_r-1][s_i_theta] + s_art[s_i_r][s_i_theta+1]) * s_x[s_i_r-1][s_i_theta+1] /* Top Left */                                      
            - 0.25 * (s_art[s_i_r+1][s_i_theta] + s_art[s_i_r][s_i_theta+1]) * s_x[s_i_r+1][s_i_theta+1] /* Top Right */                                   
        );
    }
    else if((isOnInnerBoundary && !DirBC_Interior)){

        double coeff1 = 0.5 * (k1 + k2) / h1;
        double coeff2 = 0.5 * (k1 + k2) / h2;
        double coeff3 = 0.5 * (h1 + h2) / k1;
        double coeff4 = 0.5 * (h1 + h2) / k2;

        result[center_index] = rhs[center_index] - 

            (0.25 * (h1 + h2) * (k1 + k2) * coeff_beta_cache[i_r] * fabs(s_detDF[s_i_r][s_i_theta]) * s_x[s_i_r][s_i_theta] /* beta_{i,j} */ 
                                                                                                               
            - coeff1 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r-1][s_i_theta]) * (s_x[s_i_r-1][s_i_theta] - s_x[s_i_r][s_i_theta]) /* Left - Center: (Left) */          
            - coeff2 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r+1][s_i_theta]) * (s_x[s_i_r+1][s_i_theta] - s_x[s_i_r][s_i_theta]) /* Right - Center: (Right) */      
            - coeff3 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta-1]) * (s_x[s_i_r][s_i_theta-1] - s_x[s_i_r][s_i_theta]) /* Bottom - Center: (Bottom) */  
            - coeff4 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta+1]) * (s_x[s_i_r][s_i_theta+1] - s_x[s_i_r][s_i_theta]) /* Top - Center: (Top) */              
                                                                                                                                              
            + 0.25 * (s_art[s_i_r+1][s_i_theta] + s_art[s_i_r][s_i_theta-1]) * s_x[s_i_r+1][s_i_theta-1] /* Bottom Right */                                                              
            - 0.25 * (s_art[s_i_r+1][s_i_theta] + s_art[s_i_r][s_i_theta+1]) * s_x[s_i_r+1][s_i_theta+1] /* Top Right */                                   
        );

    }
    else if((isOnInnerBoundary && DirBC_Interior) || isOnOuterBoundary){
        result[center_index] = rhs[center_index] - s_x[s_i_r][s_i_theta];
    }
}


/* ------------------ */
/* result = rhs - A*x */
void ResidualTakeGPU::computeResidual(GPU_Vector<double>& result, const GPU_Vector<double>& rhs, const GPU_Vector<double>& x) const {

    const PolarGrid& grid = level_.grid();

    assert(result.size() == grid.numberOfNodes());
    assert(rhs.size() == grid.numberOfNodes());
    assert(x.size() == grid.numberOfNodes());

    const GPU_Vector<double>& sin_theta_cache = level_.levelCache().GPU_sin_theta();
    const GPU_Vector<double>& cos_theta_cache = level_.levelCache().GPU_cos_theta();

    const GPU_Vector<double>& coeff_alpha_cache = level_.levelCache().GPU_coeff_alpha();
    const GPU_Vector<double>& coeff_beta_cache = level_.levelCache().GPU_coeff_beta();

    DomainGeometry* device_domain_geometry;
    cudaMalloc(&device_domain_geometry, sizeof(DomainGeometry));
    cudaMemcpy(device_domain_geometry, &domain_geometry_, sizeof(DomainGeometry), cudaMemcpyHostToDevice);

    /* We use precomputed DensityProfileCoefficients values. */
    // DensityProfileCoefficients* device_density_profile;
    // cudaMalloc(&device_density_profile, sizeof(DensityProfileCoefficients));
    // cudaMemcpy(device_density_profile, &density_profile_coefficients_, sizeof(DensityProfileCoefficients), cudaMemcpyHostToDevice);
    
    dim3 threadsPerBlock(16, 16);
    dim3 numBlocks((grid.nr() + 14 - 1) / 14,
                   (grid.ntheta() + 14 - 1) / 14);

    applyResidual_kernel<<<numBlocks, threadsPerBlock>>>(
        result.data(), rhs.data(), x.data(), 
        level_.device_grid(), DirBC_Interior_,
        device_domain_geometry, 
        coeff_alpha_cache.data(), coeff_beta_cache.data(), 
        sin_theta_cache.data(), cos_theta_cache.data()
    );

    cudaDeviceSynchronize();

    /* We use precomputed DensityProfileCoefficients values. */
    cudaFree(device_domain_geometry);
    // cudaFree(device_density_profile);
}

