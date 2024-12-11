#include "../../../include/ExtrapolatedSmoother/ExtrapolatedSmootherTakeGPU/extrapolatedSmoother.h"

__global__ void build_AscMatrices_kernel(
    double* circle_lower_diagonals, double* circle_main_diagonals, double* circle_upper_diagonals,
    double* radial_lower_diagonals, double* radial_main_diagonals, double* radial_upper_diagonals,
    int* d_inner_boundary_matrix_row_indices, 
    int* d_inner_boundary_matrix_column_indices,
    double* d_inner_boundary_matrix_values,
    PolarGrid* grid, bool DirBC_Interior,
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

    /* Node lies outside of the grid. */
    if (i_r < 0 || i_r >= grid->nr() || i_theta < 0 || i_theta >= grid->ntheta()) return;

    /* Share expensive to compute values. */
    __shared__ double s_detDF[16][16];
    __shared__ double s_arr[16][16];
    __shared__ double s_att[16][16];

    /* Indexing into shared data */
    int s_i_r = threadIdx.x;
    int s_i_theta = threadIdx.y;

    /* Current node index */
    int center_index = grid->index(i_r, i_theta);
    
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

    /* Share data to nodes in local grid block. */
    s_detDF[s_i_r][s_i_theta] = detDF;
    s_arr[s_i_r][s_i_theta] = arr;
    s_att[s_i_r][s_i_theta] = att;

    __syncthreads();

    /* Node lies outside of the grid. */
    if (i_r < 0 || i_r >= grid->nr() || i_theta < 0 || i_theta >= grid->ntheta()) return;
    /* Node lies on the halo. */
    if(s_i_r == 0 || s_i_r == 15 || s_i_theta == 0 || s_i_theta == 15) return;

    /* Compute neighbor distances */
    bool isOnInnerBoundary = (i_r == 0);
    bool isOnOuterBoundary = (i_r == grid->nr() - 1);

    double h1 = DirBC_Interior ? 
        ((!isOnInnerBoundary) ? grid->radialSpacing(i_r - 1) : 0.0) :
        ((!isOnInnerBoundary) ? grid->radialSpacing(i_r - 1) : 2.0 * grid->radius(0));
    double h2 = (!isOnOuterBoundary) ? grid->radialSpacing(i_r) : 0.0;
    double k1 = grid->angularSpacing(i_theta - 1);                                                          
    double k2 = grid->angularSpacing(i_theta);

    double coeff1 = (h1 != 0.0) ? 0.5 * (k1 + k2) / h1 : 0.0;
    double coeff2 = (h2 != 0.0) ? 0.5 * (k1 + k2) / h2 : 0.0;
    double coeff3 = (k1 != 0.0) ? 0.5 * (h1 + h2) / k1 : 0.0;
    double coeff4 = (k2 != 0.0) ? 0.5 * (h1 + h2) / k2 : 0.0;

    int i_theta_M1 = grid->wrapThetaIndex(i_theta - 1);
    int i_theta_P1 = grid->wrapThetaIndex(i_theta + 1);

    const int numberSmootherCircles = grid->numberSmootherCircles(); 
    const int lengthSmootherRadial  = grid->lengthSmootherRadial(); 

    int row, column;
    double value;

    int circle_m = grid->ntheta();
    int radial_m = grid->lengthSmootherRadial();

    /* ------------------------------------------ */
    /* Circle Section: Node in the inner boundary */
    /* ------------------------------------------ */     
    if(i_r == 0){
        if(DirBC_Interior){
            /* ------------------------------------------------ */
            /* Case 1: Dirichlet boundary on the inner boundary */
            /* ------------------------------------------------ */      
            d_inner_boundary_matrix_row_indices[i_theta] = i_theta + 1;
            d_inner_boundary_matrix_column_indices[i_theta] = i_theta + 1;
            d_inner_boundary_matrix_values[i_theta] = 1.0;
        }
        else{
            /* ------------------------------------------------------------- */
            /* Case 2: Across origin discretization on the interior boundary */
            /* ------------------------------------------------------------- */
            /* h1 gets replaced with 2 * R0. */
            /* (i_r-1,i_theta) gets replaced with (i_r, i_theta + (grid.ntheta()>>1)). */
            /* Some more adjustments from the changing the 9-point stencil to the artifical 7-point stencil. */

            int i_theta_AcrossOrigin = grid->wrapThetaIndex(i_theta + grid->ntheta() / 2);

            const int center_index = i_theta;
            const int left_index   = i_theta_AcrossOrigin;

            int center_nz_index;
            if (!DirBC_Interior) {
                if (i_theta % 2 == 0) {
                    center_nz_index = 3 * (i_theta / 2);
                }
                else {
                    center_nz_index = 3 * (i_theta / 2) + 1;
                }
            }
            else {
                center_nz_index = i_theta;
            }


            if (i_theta & 1) {           
                /* i_theta % 2 == 1 */   
                /* -| X | O | X | */     
                /* -|   |   |   | */     
                /* -| Õ | O | O | */     
                /* -|   |   |   | */     
                /* -| X | O | X | */

                /* Stencil Indexing */
                int StencilType_Center = 0;
                int StencilType_Left = 1;

                double center_value = (
                    + 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta_cache[i_r] * fabs(s_detDF[s_i_r][s_i_theta])
                    + coeff1 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r-1][s_i_theta]) 
                    + coeff2 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r+1][s_i_theta])
                    + coeff3 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta-1]) 
                    + coeff4 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta+1])
                );
                d_inner_boundary_matrix_row_indices[center_nz_index + StencilType_Center] = center_index + 1;
                d_inner_boundary_matrix_column_indices[center_nz_index + StencilType_Center] = center_index+ 1;
                d_inner_boundary_matrix_values[center_nz_index + StencilType_Center] = center_value;

                double left_value = -coeff1 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r-1][s_i_theta]);
                d_inner_boundary_matrix_row_indices[center_nz_index + StencilType_Left] = center_index + 1;
                d_inner_boundary_matrix_column_indices[center_nz_index + StencilType_Left] = left_index + 1;
                d_inner_boundary_matrix_values[center_nz_index + StencilType_Left] = left_value;
            }                            
            else {                       
                /* i_theta % 2 == 0 */   
                /* -| O | O | O | */     
                /* -|   |   |   | */     
                /* -| X̃ | O | X | */     
                /* -|   |   |   | */     
                /* -| O | O | O | */      
                int StencilType_Center = 0;     
                d_inner_boundary_matrix_row_indices[center_nz_index + StencilType_Center] = i_theta + 1;
                d_inner_boundary_matrix_column_indices[center_nz_index + StencilType_Center]  = i_theta + 1;
                d_inner_boundary_matrix_values[center_nz_index + StencilType_Center]  = 1.0;
            }         
        }
    }
    /* ------------------------------------------ */
    /* Node in the interior of the Circle Section */
    /* ------------------------------------------ */
    else if(i_r > 0 && i_r < numberSmootherCircles){
        int center_index = i_theta;                          
        int bottom_index = i_theta_M1;                       
        int top_index    = i_theta_P1;                       
        /* -------------------------- */                     
        /* Cyclic Tridiagonal Section */                     
        /* i_r % 2 == 1               */                     
        if (i_r & 1) {                                       
            /* i_theta % 2 == 1 */ /* i_theta % 2 == 0 */    
            /* | X | O | X | */ /* | O | O | O | */          
            /* |   |   |   | */ /* |   |   |   | */          
            /* | 0 | Õ | O | */ /* | X | Õ | X | */ 
            /* |   |   |   | */ /* |   |   |   | */          
            /* | X | O | X | */ /* | O | O | O | */          
                                                                
            /* Center: (Left, Right, Bottom, Top) */  
            row = center_index;
            column = center_index;
            value = (
                + 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta_cache[i_r] * fabs(s_detDF[s_i_r][s_i_theta])
                + coeff1 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r-1][s_i_theta]) 
                + coeff2 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r+1][s_i_theta])
                + coeff3 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta-1]) 
                + coeff4 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta+1])
            );
            if (row == column) circle_main_diagonals[i_r * circle_m + row] = value;
            else if (row == column + 1) circle_lower_diagonals[i_r * circle_m + row] = value;
            else if (row == column - 1) circle_upper_diagonals[i_r * circle_m + row] = value;
            else if (row == 0 && column == circle_m - 1) circle_lower_diagonals[i_r * circle_m + row] = value;
            else if (row == circle_m - 1 && column == 0) circle_upper_diagonals[i_r * circle_m + row] = value;

            /* Bottom */  
            row = center_index;
            column = bottom_index;
            value = - coeff3 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta-1]);
            if (row == column) circle_main_diagonals[i_r * circle_m + row] = value;
            else if (row == column + 1) circle_lower_diagonals[i_r * circle_m + row] = value;
            else if (row == column - 1) circle_upper_diagonals[i_r * circle_m + row] = value;
            else if (row == 0 && column == circle_m - 1) circle_lower_diagonals[i_r * circle_m + row] = value;
            else if (row == circle_m - 1 && column == 0) circle_upper_diagonals[i_r * circle_m + row] = value;

            /* Top */  
            row = center_index;
            column = top_index;
            value = - coeff4 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta+1]);
            if (row == column) circle_main_diagonals[i_r * circle_m + row] = value;
            else if (row == column + 1) circle_lower_diagonals[i_r * circle_m + row] = value;
            else if (row == column - 1) circle_upper_diagonals[i_r * circle_m + row] = value;
            else if (row == 0 && column == circle_m - 1) circle_lower_diagonals[i_r * circle_m + row] = value;
            else if (row == circle_m - 1 && column == 0) circle_upper_diagonals[i_r * circle_m + row] = value;
        }

        /* ---------------- */                               
        /* Diagonal Section */                               
        /* i_r % 2 == 0     */                               
        else {                                               
            /* i_theta % 2 == 1 */ /* i_theta % 2 == 0 */    
            /* | O | X | O | */ /* | O | O | O | */          
            /* |   |   |   | */ /* |   |   |   | */          
            /* | O | Õ | O | */ /* | O | X̃ | O | */ 
            /* |   |   |   | */ /* |   |   |   | */          
            /* | O | X | O | */ /* | O | O | O | */         

            if (i_theta & 1) { /* i_theta % 2 == 1 */        
                /* Center: (Left, Right, Bottom, Top) */     
                row = center_index;
                column = center_index;
                value = (
                    + 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta_cache[i_r] * fabs(s_detDF[s_i_r][s_i_theta])
                    + coeff1 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r-1][s_i_theta]) 
                    + coeff2 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r+1][s_i_theta])
                    + coeff3 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta-1]) 
                    + coeff4 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta+1])
                );
                circle_main_diagonals[i_r * circle_m + row] = value; 
            }                                                
            else { /* i_theta % 2 == 0 */                    
                /* Center: Coarse */                         
                row = center_index;
                column = center_index;
                value = 1.0;
                circle_main_diagonals[i_r * circle_m + row] = value;        
            }        
        }
    }
    /* --------------------------------------------- */
    /* Radial Section: Node next to circular section */
    /* --------------------------------------------- */
    else if(i_r == numberSmootherCircles){
        int center_index = i_r - numberSmootherCircles;
        int right_index = i_r - numberSmootherCircles + 1;
        if (i_theta & 1) {                                       
            /* i_theta % 2 == 1 and i_r % 2 == 1 */              
            /* | X | O | X || O   X   O   X  */                  
            /* |   |   |   || -------------- */                  
            /* | 0 | O | O || Õ   O   O   O  */                  
            /* |   |   |   || -------------- */                  
            /* | X | O | X || O   X   O   X  */                  
            /* or */                                             
            /* i_theta % 2 == 1 and i_r % 2 == 0 */              
            /* | O | X | O || X   O   X   O  */                  
            /* |   |   |   || -------------- */                  
            /* | 0 | O | O || Õ   O   O   O  */                  
            /* |   |   |   || -------------- */                  
            /* | O | X | O || X   O   X   O  */ 

            /* Center: (Left, Right, Bottom, Top) */  
            row = center_index;
            column = center_index;
            value = (
                + 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta_cache[i_r] * fabs(s_detDF[s_i_r][s_i_theta])
                + coeff1 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r-1][s_i_theta]) 
                + coeff2 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r+1][s_i_theta])
                + coeff3 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta-1]) 
                + coeff4 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta+1])
            );   
            if (row == column) radial_main_diagonals[i_theta * radial_m + row] = value;
            else if (row == column + 1) radial_lower_diagonals[i_theta * radial_m + row] = value;
            else if (row == column - 1) radial_upper_diagonals[i_theta * radial_m + row] = value;
            else if (row == 0 && column == radial_m - 1) radial_lower_diagonals[i_theta * radial_m + row] = value;
            else if (row == radial_m - 1 && column == 0) radial_upper_diagonals[i_theta * radial_m + row] = value;

            /* Right */  
            row = center_index;
            column = right_index;
            value = - coeff2 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r + 1][s_i_theta]);  
            if (row == column) radial_main_diagonals[i_theta * radial_m + row] = value;
            else if (row == column + 1) radial_lower_diagonals[i_theta * radial_m + row] = value;
            else if (row == column - 1) radial_upper_diagonals[i_theta * radial_m + row] = value;
            else if (row == 0 && column == radial_m - 1) radial_lower_diagonals[i_theta * radial_m + row] = value;
            else if (row == radial_m - 1 && column == 0) radial_upper_diagonals[i_theta * radial_m + row] = value;                                                           
        }
        else {                                                                                                                       
            if (i_r & 1) {                                       
                /* i_theta % 2 == 0 and i_r % 2 == 1 */          
                /* | O | O | O || O   O   O   O  */              
                /* |   |   |   || -------------- */              
                /* | X | O | X || Õ   X   O   X  */              
                /* |   |   |   || -------------- */              
                /* | O | O | O || O   O   O   O  */              
                                                                    
                /* Center: (Left, Right, Bottom, Top) */         
                row = center_index;
                column = center_index;
                value = (
                    + 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta_cache[i_r] * fabs(s_detDF[s_i_r][s_i_theta])
                    + coeff1 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r-1][s_i_theta]) 
                    + coeff2 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r+1][s_i_theta])
                    + coeff3 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta-1]) 
                    + coeff4 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta+1])
                );   
                radial_main_diagonals[i_theta * radial_m + row] = value;
            }                                                    
            else {                                               
                /* i_theta % 2 == 0 and i_r % 2 == 0 */          
                /* | O | O | O || O   O   O   O  */              
                /* |   |   |   || -------------- */              
                /* | O | X | O || X̃   O   X   O  */              
                /* |   |   |   || -------------- */              
                /* | O | O | O || O   O   O   O  */              
                /* Center: Coarse */                             
                row = center_index;
                column = center_index;
                value = 1.0;
                radial_main_diagonals[i_theta * radial_m + row] = value;            
            }                                                    
        }  
    }
    /* ------------------------------------------ */
    /* Node in the interior of the Radial Section */
    /* ------------------------------------------ */
    else if(i_r > numberSmootherCircles && i_r < grid->nr()-2){
        int center_index = i_r - numberSmootherCircles;
        int left_index   = i_r - numberSmootherCircles - 1;
        int right_index  = i_r - numberSmootherCircles + 1;

        /* ------------------- */                                                                         
        /* Tridiagonal Section */                                                                         
        /* i_theta % 2 == 1    */                                                                         
        if (i_theta & 1) {                                                                                
            /* i_r % 2 == 1 */ /* i_r % 2 == 0 */                                                         
            /* ---------- */ /* ---------- */                                                             
            /* X   O   X  */ /* O   X   O  */                                                             
            /* ---------- */ /* ---------- */                                                             
            /* O   Õ   O  */ /* O   Õ   O  */                                                    
            /* ---------- */ /* ---------- */                                                             
            /* X   O   X  */ /* O   X   O  */                                                             
            /* ---------- */ /* ---------- */

            /* Center: (Left, Right, Bottom, Top) */  
            row = center_index;
            column = center_index;
            value = (
                + 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta_cache[i_r] * fabs(s_detDF[s_i_r][s_i_theta])
                + coeff1 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r-1][s_i_theta]) 
                + coeff2 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r+1][s_i_theta])
                + coeff3 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta-1]) 
                + coeff4 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta+1])
            );   

            if (row == column) radial_main_diagonals[i_theta * radial_m + row] = value;
            else if (row == column + 1) radial_lower_diagonals[i_theta * radial_m + row] = value;
            else if (row == column - 1) radial_upper_diagonals[i_theta * radial_m + row] = value;
            else if (row == 0 && column == radial_m - 1) radial_lower_diagonals[i_theta * radial_m + row] = value;
            else if (row == radial_m - 1 && column == 0) radial_upper_diagonals[i_theta * radial_m + row] = value;

            /* Left */  
            row = center_index;
            column = left_index;
            value = -coeff1 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r-1][s_i_theta]);
            if (row == column) radial_main_diagonals[i_theta * radial_m + row] = value;
            else if (row == column + 1) radial_lower_diagonals[i_theta * radial_m + row] = value;
            else if (row == column - 1) radial_upper_diagonals[i_theta * radial_m + row] = value;
            else if (row == 0 && column == radial_m - 1) radial_lower_diagonals[i_theta * radial_m + row] = value;
            else if (row == radial_m - 1 && column == 0) radial_upper_diagonals[i_theta * radial_m + row] = value;

            /* Right */  
            row = center_index;
            column = right_index;
            value = - coeff2 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r+1][s_i_theta]); 
            if (row == column) radial_main_diagonals[i_theta * radial_m + row] = value;
            else if (row == column + 1) radial_lower_diagonals[i_theta * radial_m + row] = value;
            else if (row == column - 1) radial_upper_diagonals[i_theta * radial_m + row] = value;
            else if (row == 0 && column == radial_m - 1) radial_lower_diagonals[i_theta * radial_m + row] = value;
            else if (row == radial_m - 1 && column == 0) radial_upper_diagonals[i_theta * radial_m + row] = value;
     
        }
        /* ---------------- */                                                                            
        /* Diagonal Section */                                                                            
        /* i_theta % 2 == 0 */                                                                            
        else {                                                                                            
            /* i_r % 2 == 1 */ /* i_r % 2 == 0 */                                                         
            /* ---------- */ /* ---------- */                                                             
            /* O   O   O  */ /* O   O   O  */                                                             
            /* ---------- */ /* ---------- */                                                             
            /* X   Õ   X  */ /* O   X̃   O  */                                                    
            /* ---------- */ /* ---------- */                                                             
            /* O   O   O  */ /* O   O   O  */                                                             
            /* ---------- */ /* ---------- */                                                                                                                                                                              
            if (i_r & 1) { /* i_r % 2 == 1 */                                                             
                /* Center: (Left, Right, Bottom, Top) */                                                  
                row = center_index;
                column = center_index;
                value = (
                    + 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta_cache[i_r] * fabs(s_detDF[s_i_r][s_i_theta])
                    + coeff1 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r-1][s_i_theta]) 
                    + coeff2 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r+1][s_i_theta])
                    + coeff3 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta-1]) 
                    + coeff4 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta+1])
                );   
                radial_main_diagonals[i_theta * radial_m + row] = value;                                 
            }                                                                                             
            else { /* i_r % 2 == 0 */                                                                     
                /* Center: Coarse */                                                                      
                row = center_index;
                column = center_index;
                value = 1.0;
                radial_main_diagonals[i_theta * radial_m + row] = value;                                                               
            }                                                                                             
        }                
    }
    /* ------------------------------------------- */
    /* Radial Section: Node next to outer boundary */
    /* ------------------------------------------- */
    else if(i_r == grid->nr()-2){
        assert(i_r % 2 == 1);    

        int center_index = i_r - numberSmootherCircles;
        int left_index   = i_r - numberSmootherCircles - 1;
        int right_index   = i_r - numberSmootherCircles + 1;

        if (i_theta & 1) {                                                                                      
            /* i_theta % 2 == 1 */                                                                              
            /* ---------------|| */                                                                             
            /* O   X   O   X  || */                                                                             
            /* ---------------|| */                                                                             
            /* O   O   Õ   O  || */                                                                             
            /* ---------------|| */                                                                             
            /* O   X   O   X  || */                                                                             
            /* ---------------|| */

            /* Center: (Left, Right, Bottom, Top) */  
            row = center_index;
            column = center_index;
            value = (
                + 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta_cache[i_r] * fabs(s_detDF[s_i_r][s_i_theta])
                + coeff1 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r-1][s_i_theta]) 
                + coeff2 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r+1][s_i_theta])
                + coeff3 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta-1]) 
                + coeff4 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta+1])
            );   
            if (row == column) radial_main_diagonals[i_theta * radial_m + row] = value;
            else if (row == column + 1) radial_lower_diagonals[i_theta * radial_m + row] = value;
            else if (row == column - 1) radial_upper_diagonals[i_theta * radial_m + row] = value;
            else if (row == 0 && column == radial_m - 1) radial_lower_diagonals[i_theta * radial_m + row] = value;
            else if (row == radial_m - 1 && column == 0) radial_upper_diagonals[i_theta * radial_m + row] = value;

            /* Left */  
            row = center_index;
            column = left_index;
            value = -coeff1 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r-1][s_i_theta]);
            if (row == column) radial_main_diagonals[i_theta * radial_m + row] = value;
            else if (row == column + 1) radial_lower_diagonals[i_theta * radial_m + row] = value;
            else if (row == column - 1) radial_upper_diagonals[i_theta * radial_m + row] = value;
            else if (row == 0 && column == radial_m - 1) radial_lower_diagonals[i_theta * radial_m + row] = value;
            else if (row == radial_m - 1 && column == 0) radial_upper_diagonals[i_theta * radial_m + row] = value;

            /* Right */  
            row = center_index;
            column = right_index;
            value = 0.0;  /* Make tridiagonal matrix symmetric */  
            if (row == column) radial_main_diagonals[i_theta * radial_m + row] = value;
            else if (row == column + 1) radial_lower_diagonals[i_theta * radial_m + row] = value;
            else if (row == column - 1) radial_upper_diagonals[i_theta * radial_m + row] = value;
            else if (row == 0 && column == radial_m - 1) radial_lower_diagonals[i_theta * radial_m + row] = value;
            else if (row == radial_m - 1 && column == 0) radial_upper_diagonals[i_theta * radial_m + row] = value;
        }                                                                                                       
        else {                                                                                                  
            /* i_theta % 2 == 0 */                                                                              
            /* ---------------|| */                                                                             
            /* O   O   O   O  || */                                                                             
            /* ---------------|| */                                                                             
            /* O   X   Õ   X  || */                                                                             
            /* ---------------|| */                                                                             
            /* O   O   O   O  || */                                                                             
            /* ---------------|| */                                                                             
                                                                                                                
            /* Center: (Left, Right, Bottom, Top) */  
            row = center_index;
            column = center_index;
            value = (
                + 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta_cache[i_r] * fabs(s_detDF[s_i_r][s_i_theta])
                + coeff1 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r-1][s_i_theta]) 
                + coeff2 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r+1][s_i_theta])
                + coeff3 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta-1]) 
                + coeff4 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta+1])
            );   
            radial_main_diagonals[i_theta * radial_m + row] = value;                                                             
        }           
    }
    /* ------------------------------------------ */
    /* Radial Section: Node on the outer boundary */
    /* ------------------------------------------ */
    else if(i_r == grid->nr()-1){
        assert(!i_r % 2 == 0); 

        int center_index = i_r - numberSmootherCircles;
        int left_index   = i_r - numberSmootherCircles - 1;  

        if (i_theta & 1) {                                                                                      
            /* i_theta % 2 == 1 */                                                                              
            /* -----------|| */                                                                                 
            /* X   O   X  || */                                                                                 
            /* -----------|| */                                                                                 
            /* O   O   Õ  || */                                                                                 
            /* -----------|| */                                                                                 
            /* X   O   X  || */                                                                                 
            /* -----------|| */   
            row = center_index;
            column = center_index;
            value = 1.0;
            if (row == column) radial_main_diagonals[i_theta * radial_m + row] = value;
            else if (row == column + 1) radial_lower_diagonals[i_theta * radial_m + row] = value;
            else if (row == column - 1) radial_upper_diagonals[i_theta * radial_m + row] = value;
            else if (row == 0 && column == radial_m - 1) radial_lower_diagonals[i_theta * radial_m + row] = value;
            else if (row == radial_m - 1 && column == 0) radial_upper_diagonals[i_theta * radial_m + row] = value;

            row = center_index;
            column = left_index;
            value = 0.0;  /* Make tridiagonal matrix symmetric */   
            if (row == column) radial_main_diagonals[i_theta * radial_m + row] = value;
            else if (row == column + 1) radial_lower_diagonals[i_theta * radial_m + row] = value;
            else if (row == column - 1) radial_upper_diagonals[i_theta * radial_m + row] = value;
            else if (row == 0 && column == radial_m - 1) radial_lower_diagonals[i_theta * radial_m + row] = value;
            else if (row == radial_m - 1 && column == 0) radial_upper_diagonals[i_theta * radial_m + row] = value;
        }                                                                                                       
        else {                                                                                                  
            /* i_theta % 2 == 0 */                                                                              
            /* -----------|| */                                                                                 
            /* O   O   O  || */                                                                                 
            /* -----------|| */                                                                                 
            /* X   O   X̃  || */                                                                                 
            /* -----------|| */                                                                                 
            /* O   O   O  || */                                                                                 
            /* -----------|| */                                                                                 
            row = center_index;
            column = center_index;
            value = 1.0;
            radial_main_diagonals[i_theta * radial_m + row] = value;                                                  
        }          
    }
}



void ExtrapolatedSmootherTakeGPU::buildAscMatrices()
{
    const PolarGrid& grid = level_.grid();

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
    
    /* The stencil is computed on a 14x14 grid. */
    /* We use a 16x16 halo block to compute the expensive values. */
    /* This minimizes threads beeing idle. */
    dim3 threadsPerBlock(16, 16);
    dim3 numBlocks((grid.nr() + 14 - 1) / 14, (grid.ntheta() + 14 - 1) / 14);
    build_AscMatrices_kernel<<<numBlocks, threadsPerBlock>>>(
        circle_lower_diagonals_, circle_main_diagonals_, circle_upper_diagonals_,
        radial_lower_diagonals_, radial_main_diagonals_, radial_upper_diagonals_,
        d_inner_boundary_matrix_row_indices_,
        d_inner_boundary_matrix_column_indices_,
        d_inner_boundary_matrix_values_,
        level_.device_grid(), DirBC_Interior_,
        device_domain_geometry, 
        coeff_alpha_cache.data(), coeff_beta_cache.data(), 
        sin_theta_cache.data(), cos_theta_cache.data()
    );
    cudaDeviceSynchronize();

    int nnz = DirBC_Interior_ ? grid.ntheta() : grid.ntheta() / 2 + 2 * (grid.ntheta() / 2); 
    cudaMemcpy(inner_boundary_matrix_row_indices_.get(), d_inner_boundary_matrix_row_indices_, nnz * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(inner_boundary_matrix_column_indices_.get(), d_inner_boundary_matrix_column_indices_, nnz * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(inner_boundary_matrix_values_.get(), d_inner_boundary_matrix_values_, nnz * sizeof(double), cudaMemcpyDeviceToHost);

    /* We use precomputed DensityProfileCoefficients values. */
    cudaFree(device_domain_geometry);
    // cudaFree(device_density_profile);
}
