#include "../../../include/Smoother/SmootherTakeGPU/smoother.h"

__global__ void build_AscMatrices_kernel(
    double* circle_lower_diagonals, double* circle_main_diagonals, double* circle_upper_diagonals,
    double* radial_lower_diagonals, double* radial_main_diagonals, double* radial_upper_diagonals,
    double* csrValA, int* csrRowPtrA, int* csrColIndA,
    int* d_inner_boundary_matrix_row_indices, 
    int* d_inner_boundary_matrix_column_indices,
    double* d_inner_boundary_matrix_values,
    PolarGrid* grid, bool DirBC_Interior,
    DomainGeometry* domain_geometry,
    double* coeff_alpha_cache, double* coeff_beta_cache,
    double* sin_theta_cache, double* cos_theta_cache)
{
    /* The thread block covers a 14x14 region within a 16x16 shared memory block (1-cell halo). */
    const int global_i_r = blockIdx.x * 14 + threadIdx.x - 1;
    const int global_i_theta = blockIdx.y * 14 + threadIdx.y - 1;

    /* Adjust for across origin and periodic boundary. */
    int i_r = global_i_r;
    int i_theta = global_i_theta;
    if(i_r == -1 && !DirBC_Interior){
        i_r = 0;
        i_theta += grid->ntheta() / 2;
    }
    i_theta = grid->wrapThetaIndex(i_theta);

    /* Define bounds for valid global indices (domain + halo). */
    const int min_i_r = DirBC_Interior ? 0 : -1; 
    const int max_i_r = grid->nr() - 1;
    const int min_i_theta = -1; 
    const int max_i_theta = grid->ntheta();

    /* Exit if outside of the computational domain and halo region. */
    if (global_i_r < min_i_r || global_i_r > max_i_r || global_i_theta < min_i_theta || global_i_theta > max_i_theta) return;

    /* Allocate shared memory with padding for avoiding bank conflicts. */
    __shared__ double s_arr[16][16 + 1];
    __shared__ double s_att[16][16 + 1];

    /* Local (shared memory) thread indices. */
    const int s_i_r = threadIdx.x;
    const int s_i_theta = threadIdx.y;

    /* Compute Jacobian on current node */
    const double r = grid->radius(i_r);
    const double theta = grid->theta(i_theta);

    const double sin_theta = sin_theta_cache[i_theta];
    const double cos_theta = cos_theta_cache[i_theta];
    
    const double Jrr = domain_geometry->dFx_dr(r, theta, sin_theta, cos_theta);
    const double Jtr = domain_geometry->dFy_dr(r, theta, sin_theta, cos_theta);
    const double Jrt = domain_geometry->dFx_dt(r, theta, sin_theta, cos_theta);
    const double Jtt = domain_geometry->dFy_dt(r, theta, sin_theta, cos_theta);

    const double coeff_alpha = coeff_alpha_cache[i_r];

    const double detDF = Jrr * Jtt - Jrt * Jtr;
    const double arr = 0.5 * (Jtt * Jtt + Jrt * Jrt) * coeff_alpha / fabs(detDF);
    const double att = 0.5 * (Jtr * Jtr + Jrr * Jrr) * coeff_alpha / fabs(detDF);

    /* Share data to nodes in local grid block. */
    s_arr[s_i_r][s_i_theta] = arr;
    s_att[s_i_r][s_i_theta] = att;

    __syncthreads();

    /* Node lies outside of the grid. */
    if (global_i_r < 0 || global_i_r >= grid->nr() || global_i_theta < 0 || global_i_theta >= grid->ntheta()) return;
    /* Node lies on the halo. */
    if (s_i_r == 0 || s_i_r == 15 || s_i_theta == 0 || s_i_theta == 15) return;

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
            int n = grid->ntheta();
            int nnz = grid->ntheta();
            if (i_theta < n) {
                csrRowPtrA[i_theta] = i_theta; 
                csrColIndA[i_theta] = i_theta;
                csrValA[i_theta] = 1.0;
            }
            if (i_theta == 0) {
                csrRowPtrA[n] = n; 
            }

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
            const int bottom_index = i_theta_M1;
            const int top_index    = i_theta_P1;

            const int center_nz_index = 4 * i_theta;

            /* There are 4 non zero entries per row. */
            csrRowPtrA[i_theta] = 4 * i_theta; 
            if(i_theta == 0) csrRowPtrA[grid->ntheta()] = 4 * grid->ntheta(); 

            /* Stencil Indexing */
            int StencilType_Center = 0;
            int StencilType_Left = 1;
            int StencilType_Bottom = 2;
            int StencilType_Top = 3;

            double center_value = (
                + 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta_cache[i_r] * fabs(detDF)
                + coeff1 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r-1][s_i_theta]) 
                + coeff2 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r+1][s_i_theta])
                + coeff3 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta-1]) 
                + coeff4 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta+1])
            );
            csrColIndA[center_nz_index + StencilType_Center] = center_index;
            csrValA[center_nz_index + StencilType_Center] = center_value;
            d_inner_boundary_matrix_row_indices[center_nz_index + StencilType_Center] = center_index + 1;
            d_inner_boundary_matrix_column_indices[center_nz_index + StencilType_Center] = center_index+ 1;
            d_inner_boundary_matrix_values[center_nz_index + StencilType_Center] = center_value;

            double left_value = -coeff1 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r-1][s_i_theta]);
            csrColIndA[center_nz_index + StencilType_Left] = left_index;
            csrValA[center_nz_index + StencilType_Left] = left_value;
            d_inner_boundary_matrix_row_indices[center_nz_index + StencilType_Left] = center_index + 1;
            d_inner_boundary_matrix_column_indices[center_nz_index + StencilType_Left] = left_index + 1;
            d_inner_boundary_matrix_values[center_nz_index + StencilType_Left] = left_value;

            double bottom_value = -coeff3 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta-1]);
            csrColIndA[center_nz_index + StencilType_Bottom] = bottom_index;
            csrValA[center_nz_index + StencilType_Bottom] = bottom_value;
            d_inner_boundary_matrix_row_indices[center_nz_index + StencilType_Bottom] = center_index + 1;
            d_inner_boundary_matrix_column_indices[center_nz_index + StencilType_Bottom] = bottom_index + 1;
            d_inner_boundary_matrix_values[center_nz_index + StencilType_Bottom] = bottom_value;

            double top_value = -coeff4 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta+1]);      
            csrColIndA[center_nz_index + StencilType_Top] = top_index;
            csrValA[center_nz_index + StencilType_Top] = top_value;
            d_inner_boundary_matrix_row_indices[center_nz_index + StencilType_Top] = center_index + 1;
            d_inner_boundary_matrix_column_indices[center_nz_index + StencilType_Top] = top_index + 1;
            d_inner_boundary_matrix_values[center_nz_index + StencilType_Top] = top_value;
        }
    }
    /* ------------------------------------------ */
    /* Node in the interior of the Circle Section */
    /* ------------------------------------------ */
    else if(i_r > 0 && i_r < numberSmootherCircles){

        int center_index = i_theta;
        int bottom_index = i_theta_M1;
        int top_index = i_theta_P1;

        /* Center: (Left, Right, Bottom, Top) */  
        row = center_index;
        column = center_index;
        value = (
            + 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta_cache[i_r] * fabs(detDF)
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
    /* --------------------------------------------- */
    /* Radial Section: Node next to circular section */
    /* --------------------------------------------- */
    else if(i_r == numberSmootherCircles){

        int center_index = i_r - numberSmootherCircles;
        int right_index = i_r - numberSmootherCircles + 1;

        /* Center: (Left, Right, Bottom, Top) */  
        row = center_index;
        column = center_index;
        value = (
            + 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta_cache[i_r] * fabs(detDF)
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
    /* ------------------------------------------ */
    /* Node in the interior of the Radial Section */
    /* ------------------------------------------ */
    else if(i_r > numberSmootherCircles && i_r < grid->nr()-2){
        int center_index = i_r - numberSmootherCircles;
        int left_index   = i_r - numberSmootherCircles - 1;
        int right_index  = i_r - numberSmootherCircles + 1;

        /* Center: (Left, Right, Bottom, Top) */  
        row = center_index;
        column = center_index;
        value = (
            + 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta_cache[i_r] * fabs(detDF)
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
    /* ------------------------------------------- */
    /* Radial Section: Node next to outer boundary */
    /* ------------------------------------------- */
    else if(i_r == grid->nr()-2){
        int center_index = i_r - numberSmootherCircles;
        int left_index   = i_r - numberSmootherCircles - 1;
        int right_index   = i_r - numberSmootherCircles + 1;

        /* Center: (Left, Right, Bottom, Top) */  
        row = center_index;
        column = center_index;
        value = (
            + 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta_cache[i_r] * fabs(detDF)
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
        value = 0.0;
        if (row == column) radial_main_diagonals[i_theta * radial_m + row] = value;
        else if (row == column + 1) radial_lower_diagonals[i_theta * radial_m + row] = value;
        else if (row == column - 1) radial_upper_diagonals[i_theta * radial_m + row] = value;
        else if (row == 0 && column == radial_m - 1) radial_lower_diagonals[i_theta * radial_m + row] = value;
        else if (row == radial_m - 1 && column == 0) radial_upper_diagonals[i_theta * radial_m + row] = value;
    }
    /* ------------------------------------------ */
    /* Radial Section: Node on the outer boundary */
    /* ------------------------------------------ */
    else if(i_r == grid->nr()-1){

        int center_index = i_r - numberSmootherCircles;
        int left_index = i_r - numberSmootherCircles - 1;

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
        value = 0.0;
        if (row == column) radial_main_diagonals[i_theta * radial_m + row] = value;
        else if (row == column + 1) radial_lower_diagonals[i_theta * radial_m + row] = value;
        else if (row == column - 1) radial_upper_diagonals[i_theta * radial_m + row] = value;
        else if (row == 0 && column == radial_m - 1) radial_lower_diagonals[i_theta * radial_m + row] = value;
        else if (row == radial_m - 1 && column == 0) radial_upper_diagonals[i_theta * radial_m + row] = value;
    }
}



void SmootherTakeGPU::buildAscMatrices()
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
        csrValA_, csrRowPtrA_, csrColIndA_,
        d_inner_boundary_matrix_row_indices_,
        d_inner_boundary_matrix_column_indices_,
        d_inner_boundary_matrix_values_,
        level_.device_grid(), DirBC_Interior_,
        device_domain_geometry, 
        coeff_alpha_cache.data(), coeff_beta_cache.data(), 
        sin_theta_cache.data(), cos_theta_cache.data()
    );
    cudaDeviceSynchronize();

    int nnz = DirBC_Interior_ ? grid.ntheta() : 4 * grid.ntheta();
    cudaMemcpy(inner_boundary_matrix_row_indices_.get(), d_inner_boundary_matrix_row_indices_, nnz * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(inner_boundary_matrix_column_indices_.get(), d_inner_boundary_matrix_column_indices_, nnz * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(inner_boundary_matrix_values_.get(), d_inner_boundary_matrix_values_, nnz * sizeof(double), cudaMemcpyDeviceToHost);

    /* We use precomputed DensityProfileCoefficients values. */
    cudaFree(device_domain_geometry);
    // cudaFree(device_density_profile);
}
