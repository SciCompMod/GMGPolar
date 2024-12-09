#include "../../../include/Smoother/SmootherTakeGPU/smoother.h"


// __global__ void factorize_Asc_Circle_kernel(
//     GPU_SymmetricTridiagonalSolver<double>** circle_tridiagonal_solver,
//     PolarGrid* grid
// )
// {
//     int tid = threadIdx.x + blockIdx.x * blockDim.x;

//     if(tid <= 0 || tid >= grid->numberSmootherCircles()) return;

//     GPU_SymmetricTridiagonalSolver<double>* matrix = circle_tridiagonal_solver[tid];

//     matrix->factorize();
// }

// __global__ void factorize_Asc_Radial_kernel(
//     GPU_SymmetricTridiagonalSolver<double>** radial_tridiagonal_solver,
//     PolarGrid* grid
// )
// {
//     int tid = threadIdx.x + blockIdx.x * blockDim.x;

//     if(tid < 0 || tid >= grid->ntheta()) return;

//     GPU_SymmetricTridiagonalSolver<double>* matrix = radial_tridiagonal_solver[tid];

//     matrix->factorize();
// }

// __global__ void build_Asc_kernel(
//     GPU_SymmetricTridiagonalSolver<double>** circle_tridiagonal_solver,
//     GPU_SymmetricTridiagonalSolver<double>** radial_tridiagonal_solver,
//     PolarGrid* grid, bool DirBC_Interior,
//     DomainGeometry* domain_geometry,
//     double* coeff_alpha_cache, double* coeff_beta_cache,
//     double* sin_theta_cache, double* cos_theta_cache
// )
// {
//     int i_r = blockIdx.x * 14 + threadIdx.x - 1;
//     int i_theta = blockIdx.y * 14 + threadIdx.y - 1;

//     if(i_r == -1 && !DirBC_Interior){
//         i_r = 0;
//         i_theta += grid->ntheta() / 2;
//     }
//     i_theta = grid->wrapThetaIndex(i_theta);

//     if (i_r < 0 || i_r >= grid->nr() || i_theta < 0 || i_theta >= grid->ntheta()) return;

//     __shared__ double s_detDF[16][16];
//     __shared__ double s_arr[16][16];
//     __shared__ double s_att[16][16];

//     int s_i_r = threadIdx.x;
//     int s_i_theta = threadIdx.y;

//     int center_index = grid->index(i_r, i_theta);
 
//     double r = grid->radius(i_r);
//     double theta = grid->theta(i_theta);

//     double sin_theta = sin_theta_cache[i_theta];
//     double cos_theta = cos_theta_cache[i_theta];
    
//     double Jrr = domain_geometry->dFx_dr(r, theta, sin_theta, cos_theta);
//     double Jtr = domain_geometry->dFy_dr(r, theta, sin_theta, cos_theta);
//     double Jrt = domain_geometry->dFx_dt(r, theta, sin_theta, cos_theta);
//     double Jtt = domain_geometry->dFy_dt(r, theta, sin_theta, cos_theta);

//     double coeff_alpha = coeff_alpha_cache[i_r];

//     double detDF = Jrr * Jtt - Jrt * Jtr;
//     double arr = 0.5 * (Jtt * Jtt + Jrt * Jrt) * coeff_alpha / fabs(detDF);
//     double att = 0.5 * (Jtr * Jtr + Jrr * Jrr) * coeff_alpha / fabs(detDF);

//     s_detDF[s_i_r][s_i_theta] = detDF;
//     s_arr[s_i_r][s_i_theta] = arr;
//     s_att[s_i_r][s_i_theta] = att;

//     __syncthreads();

//     if (i_r < 0 || i_r >= grid->nr() || i_theta < 0 || i_theta >= grid->ntheta()) return;
//     if(s_i_r == 0 || s_i_r == 15 || s_i_theta == 0 || s_i_theta == 15) return;

//     bool isOnInnerBoundary = (i_r == 0);
//     bool isOnOuterBoundary = (i_r == grid->nr() - 1);

//     double h1 = DirBC_Interior ? 
//         ((!isOnInnerBoundary) ? grid->radialSpacing(i_r - 1) : 0.0) :
//         ((!isOnInnerBoundary) ? grid->radialSpacing(i_r - 1) : 2.0 * grid->radius(0));
//     double h2 = (!isOnOuterBoundary) ? grid->radialSpacing(i_r) : 0.0;
//     double k1 = grid->angularSpacing(i_theta - 1);                                                          
//     double k2 = grid->angularSpacing(i_theta);

//     double coeff1 = (h1 != 0.0) ? 0.5 * (k1 + k2) / h1 : 0.0;
//     double coeff2 = (h2 != 0.0) ? 0.5 * (k1 + k2) / h2 : 0.0;
//     double coeff3 = (k1 != 0.0) ? 0.5 * (h1 + h2) / k1 : 0.0;
//     double coeff4 = (k2 != 0.0) ? 0.5 * (h1 + h2) / k2 : 0.0;

//     int i_theta_M1 = grid->wrapThetaIndex(i_theta - 1);
//     int i_theta_P1 = grid->wrapThetaIndex(i_theta + 1);

//     const int numberSmootherCircles = grid->numberSmootherCircles(); 
//     const int lengthSmootherRadial  = grid->lengthSmootherRadial(); 

//     int row, column;
//     double value;

//     if(i_r == 0){
//         if(DirBC_Interior){

//         }
//         else{

//         }
//     }
//     /* ------------------------------------------ */
//     /* Node in the interior of the Circle Section */
//     /* ------------------------------------------ */
//     else if(i_r > 0 && i_r < numberSmootherCircles){

//         GPU_SymmetricTridiagonalSolver<double>* matrix = circle_tridiagonal_solver[i_r];

//         int center_index = i_theta;
//         int bottom_index = i_theta_M1;
//         int top_index = i_theta_P1;

//         /* Center: (Left, Right, Bottom, Top) */  
//         row = center_index;
//         column = center_index;
//         value = (
//             + 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta_cache[i_r] * fabs(s_detDF[s_i_r][s_i_theta])
//             + coeff1 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r-1][s_i_theta]) 
//             + coeff2 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r+1][s_i_theta])
//             + coeff3 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta-1]) 
//             + coeff4 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta+1])
//         );   
//         if (row == column)
//             matrix->main_diagonal(row) = value;
//         else if (row == column - 1)
//             matrix->sub_diagonal(row) = value;
//         else if (row == 0 && column == matrix->columns() - 1)
//             matrix->cyclic_corner_element() = value;   

//         /* Bottom */  
//         row = center_index;
//         column = bottom_index;
//         value = - coeff3 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta-1]);
//         if (row == column)
//             matrix->main_diagonal(row) = value;
//         else if (row == column - 1)
//             matrix->sub_diagonal(row) = value;
//         else if (row == 0 && column == matrix->columns() - 1)
//             matrix->cyclic_corner_element() = value;   

//         /* Top */  
//         row = center_index;
//         column = top_index;
//         value = - coeff4 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta+1]);
//         if (row == column)
//             matrix->main_diagonal(row) = value;
//         else if (row == column - 1)
//             matrix->sub_diagonal(row) = value;
//         else if (row == 0 && column == matrix->columns() - 1)
//             matrix->cyclic_corner_element() = value;

//         // matrix->cyclic_corner_element() = 8.0;

//         // printf("Cyclic corner element value: %f\n", matrix->cyclic_corner_element());
//     }
//     /* --------------------------------------------- */
//     /* Radial Section: Node next to circular section */
//     /* --------------------------------------------- */
//     else if(i_r == numberSmootherCircles){

//         GPU_SymmetricTridiagonalSolver<double>* matrix = radial_tridiagonal_solver[i_theta];

//         int center_index = i_r - numberSmootherCircles;
//         int right_index = i_r - numberSmootherCircles + 1;

//         /* Center: (Left, Right, Bottom, Top) */  
//         row = center_index;
//         column = center_index;
//         value = (
//             + 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta_cache[i_r] * fabs(s_detDF[s_i_r][s_i_theta])
//             + coeff1 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r-1][s_i_theta]) 
//             + coeff2 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r+1][s_i_theta])
//             + coeff3 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta-1]) 
//             + coeff4 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta+1])
//         );   
//         if (row == column)
//             matrix->main_diagonal(row) = value;
//         else if (row == column - 1)
//             matrix->sub_diagonal(row) = value;
//         else if (row == 0 && column == matrix->columns() - 1)
//             matrix->cyclic_corner_element() = value;   

//         /* Right */  
//         row = center_index;
//         column = right_index;
//         value = - coeff2 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r + 1][s_i_theta]);  
//         if (row == column)
//             matrix->main_diagonal(row) = value;
//         else if (row == column - 1)
//             matrix->sub_diagonal(row) = value;
//         else if (row == 0 && column == matrix->columns() - 1)
//             matrix->cyclic_corner_element() = value;
//     }
//     /* ------------------------------------------ */
//     /* Node in the interior of the Radial Section */
//     /* ------------------------------------------ */
//     else if(i_r > numberSmootherCircles && i_r < grid->nr()-2){

//         GPU_SymmetricTridiagonalSolver<double>* matrix = radial_tridiagonal_solver[i_theta];

//         int center_index = i_r - numberSmootherCircles;
//         int left_index   = i_r - numberSmootherCircles - 1;
//         int right_index  = i_r - numberSmootherCircles + 1;

//         /* Center: (Left, Right, Bottom, Top) */  
//         row = center_index;
//         column = center_index;
//         value = (
//             + 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta_cache[i_r] * fabs(s_detDF[s_i_r][s_i_theta])
//             + coeff1 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r-1][s_i_theta]) 
//             + coeff2 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r+1][s_i_theta])
//             + coeff3 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta-1]) 
//             + coeff4 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta+1])
//         );   
//         if (row == column)
//             matrix->main_diagonal(row) = value;
//         else if (row == column - 1)
//             matrix->sub_diagonal(row) = value;
//         else if (row == 0 && column == matrix->columns() - 1)
//             matrix->cyclic_corner_element() = value;   

//         /* Left */  
//         row = center_index;
//         column = left_index;
//         value = -coeff1 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r-1][s_i_theta]);
//         if (row == column)
//             matrix->main_diagonal(row) = value;
//         else if (row == column - 1)
//             matrix->sub_diagonal(row) = value;
//         else if (row == 0 && column == matrix->columns() - 1)
//             matrix->cyclic_corner_element() = value;   

//         /* Right */  
//         row = center_index;
//         column = right_index;
//         value = - coeff2 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r+1][s_i_theta]); 
//         if (row == column)
//             matrix->main_diagonal(row) = value;
//         else if (row == column - 1)
//             matrix->sub_diagonal(row) = value;
//         else if (row == 0 && column == matrix->columns() - 1)
//             matrix->cyclic_corner_element() = value;
//     }
//     /* ------------------------------------------- */
//     /* Radial Section: Node next to outer boundary */
//     /* ------------------------------------------- */
//     else if(i_r == grid->nr()-2){
//         GPU_SymmetricTridiagonalSolver<double>* matrix = radial_tridiagonal_solver[i_theta];

//         int center_index = i_r - numberSmootherCircles;
//         int left_index   = i_r - numberSmootherCircles - 1;

//         /* Center: (Left, Right, Bottom, Top) */  
//         row = center_index;
//         column = center_index;
//         value = (
//             + 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta_cache[i_r] * fabs(s_detDF[s_i_r][s_i_theta])
//             + coeff1 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r-1][s_i_theta]) 
//             + coeff2 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r+1][s_i_theta])
//             + coeff3 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta-1]) 
//             + coeff4 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta+1])
//         );   
//         if (row == column)
//             matrix->main_diagonal(row) = value;
//         else if (row == column - 1)
//             matrix->sub_diagonal(row) = value;
//         else if (row == 0 && column == matrix->columns() - 1)
//             matrix->cyclic_corner_element() = value;   

//         /* Left */  
//         row = center_index;
//         column = left_index;
//         value = -coeff1 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r-1][s_i_theta]);
//         if (row == column)
//             matrix->main_diagonal(row) = value;
//         else if (row == column - 1)
//             matrix->sub_diagonal(row) = value;
//         else if (row == 0 && column == matrix->columns() - 1)
//             matrix->cyclic_corner_element() = value;   
//     }
//     /* ------------------------------------------ */
//     /* Radial Section: Node on the outer boundary */
//     /* ------------------------------------------ */
//     else if(i_r == grid->nr()-1){
//         GPU_SymmetricTridiagonalSolver<double>* matrix = radial_tridiagonal_solver[i_theta];

//         int center_index = i_r - numberSmootherCircles;

//         row = center_index;
//         column = center_index;
//         value = 1.0;
//         if (row == column)
//             matrix->main_diagonal(row) = value;
//         else if (row == column - 1)
//             matrix->sub_diagonal(row) = value;
//         else if (row == 0 && column == matrix->columns() - 1)
//             matrix->cyclic_corner_element() = value;   
//     }
// }

// void SmootherTakeGPU::buildAscMatrices()
// {
//     /* -------------------------------------- */
//     /* Part 1: Allocate Asc Smoother matrices */
//     /* -------------------------------------- */

//     const PolarGrid& grid = level_.grid();

//     const int number_smoother_circles = grid.numberSmootherCircles();
//     const int length_smoother_radial  = grid.lengthSmootherRadial();

//     const int num_circle_nodes = grid.ntheta();
//     circle_tridiagonal_solver_.resize(number_smoother_circles);

//     const int num_radial_nodes = length_smoother_radial;
//     radial_tridiagonal_solver_.resize(grid.ntheta());

//     // ---------------- //
//     // Circular Section //
//     // ---------------- //
//     // Remark: circle_tridiagonal_solver_[0] is unitialized.
//     // Please use inner_boundary_circle_matrix_ instead!
//     for (int circle_Asc_index = 0; circle_Asc_index < number_smoother_circles; circle_Asc_index++) {

//         /* Inner boundary circle */
//         if (circle_Asc_index == 0) {
//             // // Although the matrix is symmetric, we need to store all its entries, so we disable the symmetry.
//             // const int nnz                 = getNonZeroCountCircleAsc(circle_Asc_index);
//             // inner_boundary_circle_matrix_ = SparseMatrix<double>(num_circle_nodes, num_circle_nodes, nnz);
//             // inner_boundary_circle_matrix_.is_symmetric(false);
//         }

//         /* Interior Circle Section */
//         else {
//             auto& solverMatrix = circle_tridiagonal_solver_[circle_Asc_index];
//             solverMatrix       = GPU_SymmetricTridiagonalSolver<double>(num_circle_nodes);
//             solverMatrix.is_cyclic(true);
//         }
//     }

//     // -------------- //
//     // Radial Section //
//     // -------------- //
//     for (int radial_Asc_index = 0; radial_Asc_index < grid.ntheta(); radial_Asc_index++) {
//         auto& solverMatrix = radial_tridiagonal_solver_[radial_Asc_index];
//         solverMatrix       = GPU_SymmetricTridiagonalSolver<double>(num_radial_nodes);
//         solverMatrix.is_cyclic(false);
//     }
    
//     /* ---------------------------------- */
//     /* Part 2: Fill Asc Smoother matrices */
//     /* ---------------------------------- */

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
    
//     cudaMalloc(&d_circle_tridiagonal_solver_, circle_tridiagonal_solver_.size() * sizeof(GPU_SymmetricTridiagonalSolver<double>*));
//     cudaMalloc(&d_radial_tridiagonal_solver_, radial_tridiagonal_solver_.size() * sizeof(GPU_SymmetricTridiagonalSolver<double>*));

//     for (int i = 0; i < circle_tridiagonal_solver_.size(); ++i) {
//         GPU_SymmetricTridiagonalSolver<double>* device_solver;
//         cudaMalloc(&device_solver, sizeof(GPU_SymmetricTridiagonalSolver<double>));
//         cudaMemcpy(device_solver, &circle_tridiagonal_solver_[i], sizeof(GPU_SymmetricTridiagonalSolver<double>), cudaMemcpyHostToDevice);
//         cudaMemcpy(d_circle_tridiagonal_solver_ + i, &device_solver, sizeof(GPU_SymmetricTridiagonalSolver<double>*), cudaMemcpyHostToDevice);
//     }

//     for (int i = 0; i < radial_tridiagonal_solver_.size(); ++i) {
//         GPU_SymmetricTridiagonalSolver<double>* device_solver;
//         cudaMalloc(&device_solver, sizeof(GPU_SymmetricTridiagonalSolver<double>));
//         cudaMemcpy(device_solver, &radial_tridiagonal_solver_[i], sizeof(GPU_SymmetricTridiagonalSolver<double>), cudaMemcpyHostToDevice);
//         cudaMemcpy(d_radial_tridiagonal_solver_ + i, &device_solver, sizeof(GPU_SymmetricTridiagonalSolver<double>*), cudaMemcpyHostToDevice);
//     }

//     dim3 threadsPerBlock(16, 16);
//     dim3 numBlocks((grid.nr() + 14 - 1) / 14, (grid.ntheta() + 14 - 1) / 14);
//     build_Asc_kernel<<<numBlocks, threadsPerBlock>>>(
//         d_circle_tridiagonal_solver_, d_radial_tridiagonal_solver_, 
//         level_.device_grid(), DirBC_Interior_,
//         device_domain_geometry, 
//         coeff_alpha_cache.data(), coeff_beta_cache.data(), 
//         sin_theta_cache.data(), cos_theta_cache.data()
//     );
//     cudaDeviceSynchronize();

//     factorize_Asc_Circle_kernel<<<1, grid.numberSmootherCircles()>>>(
//         d_circle_tridiagonal_solver_, level_.device_grid()
//     );
//     factorize_Asc_Radial_kernel<<<1, grid.ntheta()>>>(
//         d_radial_tridiagonal_solver_, level_.device_grid()
//     );
//     cudaDeviceSynchronize();
// }