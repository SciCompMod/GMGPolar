#include "../../../include/ExtrapolatedSmoother/ExtrapolatedSmootherTakeGPU/extrapolatedSmoother.h"

__global__ void extrapolated_applyAscOrtho_Circle_kernel(
    double* x, double* rhs, double* temp,
    double* circle_lower_diagonals, double* circle_main_diagonals, double* circle_upper_diagonals,
    double* sherman_morrison_gammas,
    PolarGrid* grid, bool DirBC_Interior,
    int start_i_r,
    DomainGeometry* domain_geometry,
    double* coeff_alpha_cache, double* coeff_beta_cache,
    double* sin_theta_cache, double* cos_theta_cache) 
{
    /* The stencil is computed on a 14x14 grid. */
    /* We use a 16x16 halo block to compute the expensive values. */
    /* This minimizes threads beeing idle. */
    int i_r = blockIdx.x * 14 + threadIdx.x - 1;
    int i_theta = blockIdx.y * 14 + threadIdx.y - 1;

    /* Adjust for across origin and periodic boundary. */
    if(i_r == -1 && !DirBC_Interior){
        i_r = 0;
        i_theta += grid->ntheta() / 2;
    }
    i_theta = grid->wrapThetaIndex(i_theta);

    /* Node lies outside the circle section with halo. */
    /* The halo node i_r == grid->numberSmootherCircles() lies in the radial section. */
    if (i_r < 0 || i_r > grid->numberSmootherCircles() || i_theta < 0 || i_theta >= grid->ntheta()) return;

    /* Share expensive to compute values. */
    __shared__ double s_arr[16][16];
    __shared__ double s_att[16][16];
    __shared__ double s_art[16][16];
    __shared__ double s_x[16][16];

    /* Indexing into shared data */
    int s_i_r = threadIdx.x;
    int s_i_theta = threadIdx.y;

    /* Current node index */
    int center_index = grid->index(i_r, i_theta);
    s_x[s_i_r][s_i_theta] = x[center_index];

    /* Compute Jacobian on current node */
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

    /* Share data to nodes in local grid block. */
    s_arr[s_i_r][s_i_theta] = arr;
    s_att[s_i_r][s_i_theta] = att;
    s_art[s_i_r][s_i_theta] = art;

    __syncthreads();

    /* Node color and smoother color doesnt match. */
    if(i_r % 2 != start_i_r) return;
    /* Node lies outside of the circle section. */
    if(i_r < 0 || i_r >= grid->numberSmootherCircles() || i_theta < 0 || i_theta >= grid->ntheta()) return;
    /* Node lies on the halo. */
    if(s_i_r == 0 || s_i_r == 15 || s_i_theta == 0 || s_i_theta == 15) return;

    /* Prepare temp for the 2nd solution in the Shermann-Morrison formula. */
    if(i_r > 0){
        if(i_r & 1){
            int matrix_index = i_r * grid->ntheta();
            if(i_theta == 0){
                temp[matrix_index + i_theta] = sherman_morrison_gammas[i_r];
            }
            else if(i_theta > 0 && i_theta < grid->ntheta()-1){
                temp[matrix_index + i_theta] = 0.0;
            }
            else if(i_theta == grid->ntheta()-1){
                temp[matrix_index + i_theta] = circle_upper_diagonals[matrix_index + grid->ntheta() - 1];
            }
        }
    }

    /* Compute neighbor distances */
    bool isOnInnerBoundary = (i_r == 0);

    double h1 = DirBC_Interior ? 
        ((!isOnInnerBoundary) ? grid->radialSpacing(i_r - 1) : 0.0) :
        ((!isOnInnerBoundary) ? grid->radialSpacing(i_r - 1) : 2.0 * grid->radius(0));
    double h2 = grid->radialSpacing(i_r);
    double k1 = grid->angularSpacing(i_theta - 1);                                                          
    double k2 = grid->angularSpacing(i_theta);

    double coeff1 = (h1 != 0.0) ? 0.5 * (k1 + k2) / h1 : 0.0;
    double coeff2 = (h2 != 0.0) ? 0.5 * (k1 + k2) / h2 : 0.0;
    double coeff3 = (k1 != 0.0) ? 0.5 * (h1 + h2) / k1 : 0.0;
    double coeff4 = (k2 != 0.0) ? 0.5 * (h1 + h2) / k2 : 0.0;

    if (!isOnInnerBoundary) {
        if (i_r & 1) {                                                                                           
            /* i_r % 2 == 1 and i_theta % 2 == 1 */                                                              
            /* | X | O | X | */                                                                                  
            /* |   |   |   | */                                                                                  
            /* | O | Õ | O | */                                                                                  
            /* |   |   |   | */                                                                                  
            /* | X | O | X | */                                                                                  
            /* or */                                                                                             
            /* i_r % 2 == 1 and i_theta % 2 == 0 */                                                              
            /* | O | O | O | */                                                                                  
            /* |   |   |   | */                                                                                  
            /* | X | Õ | X | */                                                                                  
            /* |   |   |   | */                                                                                  
            /* | O | O | O | */        
            x[center_index] = rhs[center_index] - (
                - coeff1 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r-1][s_i_theta]) * s_x[s_i_r-1][s_i_theta] /* Left */  
                - coeff2 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r+1][s_i_theta]) * s_x[s_i_r+1][s_i_theta] /* Right */   

                - 0.25 * (s_art[s_i_r-1][s_i_theta] + s_art[s_i_r][s_i_theta-1]) * s_x[s_i_r-1][s_i_theta-1] /* Bottom Left */
                + 0.25 * (s_art[s_i_r+1][s_i_theta] + s_art[s_i_r][s_i_theta-1]) * s_x[s_i_r+1][s_i_theta-1] /* Bottom Right */
                + 0.25 * (s_art[s_i_r-1][s_i_theta] + s_art[s_i_r][s_i_theta+1]) * s_x[s_i_r-1][s_i_theta+1] /* Top Left */
                - 0.25 * (s_art[s_i_r+1][s_i_theta] + s_art[s_i_r][s_i_theta+1]) * s_x[s_i_r+1][s_i_theta+1] /* Top Right */
            );                                                
        }                                                                                                        
        else {                                                                                                   
            if (i_theta & 1) {                                                                                   
                /* i_r % 2 == 0 and i_theta % 2 == 1 */                                                          
                /* | O | X | O | */                                                                              
                /* |   |   |   | */                                                                              
                /* | O | Õ | O | */                                                                              
                /* |   |   |   | */                                                                              
                /* | O | X | O | */     
                x[center_index] = rhs[center_index] - (
                    - coeff1 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r-1][s_i_theta]) * s_x[s_i_r-1][s_i_theta] /* Left */  
                    - coeff2 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r+1][s_i_theta]) * s_x[s_i_r+1][s_i_theta] /* Right */   
                    - coeff3 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta-1]) * s_x[s_i_r][s_i_theta-1] /* Bottom */               
                    - coeff4 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta+1]) * s_x[s_i_r][s_i_theta+1] /* Top */       

                    - 0.25 * (s_art[s_i_r-1][s_i_theta] + s_art[s_i_r][s_i_theta-1]) * s_x[s_i_r-1][s_i_theta-1] /* Bottom Left */
                    + 0.25 * (s_art[s_i_r+1][s_i_theta] + s_art[s_i_r][s_i_theta-1]) * s_x[s_i_r+1][s_i_theta-1] /* Bottom Right */
                    + 0.25 * (s_art[s_i_r-1][s_i_theta] + s_art[s_i_r][s_i_theta+1]) * s_x[s_i_r-1][s_i_theta+1] /* Top Left */
                    - 0.25 * (s_art[s_i_r+1][s_i_theta] + s_art[s_i_r][s_i_theta+1]) * s_x[s_i_r+1][s_i_theta+1] /* Top Right */
                );                                                                             
            }                                                                                                    
            else {                                                                                               
                /* i_r % 2 == 0 and i_theta % 2 == 0 */                                                          
                /* | O | O | O | */                                                                              
                /* |   |   |   | */                                                                              
                /* | O | X̃ | O | */                                                                              
                /* |   |   |   | */                                                                              
                /* | O | O | O | */                                                                              
                x[center_index] = s_x[s_i_r][s_i_theta];                                                                      
            }                                                                                                    
        } 
    }
    else if(isOnInnerBoundary && !DirBC_Interior){
        if (i_theta & 1) {                                                                                                                                        
            /* i_theta % 2 == 1 */                              
            /* -| X | O | X | */                                
            /* -|   |   |   | */                                
            /* -| Õ | O | O | */                                
            /* -|   |   |   | */                                
            /* -| X | O | X | */      
            x[center_index] = rhs[center_index] - (
                - coeff2 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r+1][s_i_theta]) * s_x[s_i_r+1][s_i_theta] /* Right */   
                - coeff3 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta-1]) * s_x[s_i_r][s_i_theta-1] /* Bottom */               
                - coeff4 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta+1]) * s_x[s_i_r][s_i_theta+1] /* Top */       

                + 0.25 * (s_art[s_i_r+1][s_i_theta] + s_art[s_i_r][s_i_theta-1]) * s_x[s_i_r+1][s_i_theta-1] /* Bottom Right */
                - 0.25 * (s_art[s_i_r+1][s_i_theta] + s_art[s_i_r][s_i_theta+1]) * s_x[s_i_r+1][s_i_theta+1] /* Top Right */
            );                                            
        }     
        else{
            /* i_theta % 2 == 0 */                                                                                                                             
            /* -| O | O | O | */                                                                                                                               
            /* -|   |   |   | */                                                                                                                               
            /* -| X̃ | O | X | */                                                                                                                               
            /* -|   |   |   | */                                                                                                                               
            /* -| O | O | O | */                                                                                                                               
            x[center_index] = s_x[s_i_r][s_i_theta]; 
        }   
    }                                               
    else if(isOnInnerBoundary && DirBC_Interior){
        if (i_theta & 1) {                                      
            /* i_theta % 2 == 1 */                              
            /* || X | O | X | */                                
            /* ||   |   |   | */                                
            /* || Õ | O | O | */                                
            /* ||   |   |   | */                                
            /* || X | O | X | */                                
            x[center_index] = rhs[center_index];                         
        }                                                       
        else {                                                  
            /* i_theta % 2 == 0 */                              
            /* || O | O | O | */                                
            /* ||   |   |   | */                                
            /* || X̃ | O | X | */                                
            /* ||   |   |   | */                                
            /* || O | O | O | */                                
            x[center_index] = s_x[s_i_r][s_i_theta];                           
        }                 
    }
}



void ExtrapolatedSmootherTakeGPU::applyAscOrtho_BlackCircle(
    GPU_Vector<double>& x, const GPU_Vector<double>& rhs, GPU_Vector<double>& temp, 
    DomainGeometry* device_domain_geometry)
{

    const PolarGrid& grid = level_.grid();

    const GPU_Vector<double>& sin_theta_cache = level_.levelCache().GPU_sin_theta();
    const GPU_Vector<double>& cos_theta_cache = level_.levelCache().GPU_cos_theta();

    const GPU_Vector<double>& coeff_alpha_cache = level_.levelCache().GPU_coeff_alpha();
    const GPU_Vector<double>& coeff_beta_cache = level_.levelCache().GPU_coeff_beta();

    const int start_black_circles = (grid.numberSmootherCircles() % 2 == 0) ? 1 : 0;

    dim3 threadsPerBlock(16, 16);
    dim3 numBlocks((grid.numberSmootherCircles() + 14 - 1) / 14, (grid.ntheta() + 14 - 1) / 14);
    extrapolated_applyAscOrtho_Circle_kernel<<<numBlocks, threadsPerBlock>>>(
        x.data(), rhs.data(), temp.data(),
        circle_lower_diagonals_, circle_main_diagonals_, circle_upper_diagonals_,
        sherman_morrison_gammas_,
        level_.device_grid(), DirBC_Interior_,
        start_black_circles,
        device_domain_geometry, 
        coeff_alpha_cache.data(), coeff_beta_cache.data(), 
        sin_theta_cache.data(), cos_theta_cache.data()
    );
    cudaDeviceSynchronize();
}



void ExtrapolatedSmootherTakeGPU::applyAscOrtho_WhiteCircle(
    GPU_Vector<double>& x, const GPU_Vector<double>& rhs, GPU_Vector<double>& temp, 
    DomainGeometry* device_domain_geometry)
{
    const PolarGrid& grid = level_.grid();

    const GPU_Vector<double>& sin_theta_cache = level_.levelCache().GPU_sin_theta();
    const GPU_Vector<double>& cos_theta_cache = level_.levelCache().GPU_cos_theta();

    const GPU_Vector<double>& coeff_alpha_cache = level_.levelCache().GPU_coeff_alpha();
    const GPU_Vector<double>& coeff_beta_cache = level_.levelCache().GPU_coeff_beta();

    const int start_white_circles = (grid.numberSmootherCircles() % 2 == 0) ? 0 : 1;

    dim3 threadsPerBlock(16, 16);
    dim3 numBlocks((grid.numberSmootherCircles() + 14 - 1) / 14, (grid.ntheta() + 14 - 1) / 14);
    extrapolated_applyAscOrtho_Circle_kernel<<<numBlocks, threadsPerBlock>>>(
        x.data(), rhs.data(), temp.data(),
        circle_lower_diagonals_, circle_main_diagonals_, circle_upper_diagonals_,
        sherman_morrison_gammas_,
        level_.device_grid(), DirBC_Interior_,
        start_white_circles,
        device_domain_geometry, 
        coeff_alpha_cache.data(), coeff_beta_cache.data(), 
        sin_theta_cache.data(), cos_theta_cache.data()
    );
    cudaDeviceSynchronize();
}