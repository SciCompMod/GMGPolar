#include "../../../include/Smoother/SmootherTakeGPU/smoother.h"

#include "../../../include/LinearAlgebra/Vector/gpu_vector_operations.h"

#include <cmath>



// __global__ void solve_Asc_Circle_kernel(
//     double* x, double* rhs, double* temp,
//     GPU_SymmetricTridiagonalSolver<double>** circle_tridiagonal_solver,
//     PolarGrid* grid, bool DirBC_Interior,
//     int start_i_r
// ) 
// {
//     int i_r = 2 * (threadIdx.x + blockIdx.x * blockDim.x) + start_i_r;

//     if(i_r < 0 || i_r >= grid->numberSmootherCircles()) return;

//     const int start = grid->index(i_r, 0);
//     const int end   = start + grid->ntheta();

//     if(i_r == 0){

//     }
//     else{
//         circle_tridiagonal_solver[i_r]->solveInPlace(x + start, temp + start);
//     }
// }



// __global__ void solve_Asc_Radial_kernel(
//     double* x, double* rhs, double* temp,
//     GPU_SymmetricTridiagonalSolver<double>** radial_tridiagonal_solver,
//     PolarGrid* grid, bool DirBC_Interior,
//     int start_i_theta
// ) 
// {
//     int i_theta = 2 * (threadIdx.x + blockIdx.x * blockDim.x) + start_i_theta;

//     if(i_theta <= 0 || i_theta >= grid->ntheta()) return;

//     const int start = grid->index(grid->numberSmootherCircles(), i_theta);
//     const int end   = start + grid->lengthSmootherRadial();

//     radial_tridiagonal_solver[i_theta]->solveInPlace(temp + start);
    
//     for (int i = start; i < end; i++)
//     {
//         x[i] = std::move(temp[i]);
//     }
// }




// __global__ void applyAscOrtho_Circle_kernel(
//     double* x, double* rhs, double* temp,
//     PolarGrid* grid, bool DirBC_Interior,
//     int start_i_r,
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

//     if (i_r < 0 || i_r > grid->numberSmootherCircles() || i_theta < 0 || i_theta >= grid->ntheta()) return;

//     __shared__ double s_arr[16][16];
//     __shared__ double s_art[16][16];
//     __shared__ double s_x[16][16];

//     int s_i_r = threadIdx.x;
//     int s_i_theta = threadIdx.y;

//     int center_index = grid->index(i_r, i_theta);
//     s_x[s_i_r][s_i_theta] = x[center_index];

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
//     double art = (- Jtt * Jtr - Jrt * Jrr) * coeff_alpha / fabs(detDF);

//     s_arr[s_i_r][s_i_theta] = arr;
//     s_art[s_i_r][s_i_theta] = art;

//     __syncthreads();

//     if(i_r % 2 != start_i_r) return;

//     if (i_r < 0 || i_r >= grid->numberSmootherCircles() || i_theta < 0 || i_theta >= grid->ntheta()) return;
//     if(s_i_r == 0 || s_i_r == 15 || s_i_theta == 0 || s_i_theta == 15) return;

//     bool isOnInnerBoundary = (i_r == 0);

//     double h1 = DirBC_Interior ? 
//         ((!isOnInnerBoundary) ? grid->radialSpacing(i_r - 1) : 0.0) :
//         ((!isOnInnerBoundary) ? grid->radialSpacing(i_r - 1) : 2.0 * grid->radius(0));
//     double h2 = grid->radialSpacing(i_r);
//     double k1 = grid->angularSpacing(i_theta - 1);                                                          
//     double k2 = grid->angularSpacing(i_theta);

//     if (!isOnInnerBoundary) {   

//         double coeff1 = 0.5 * (k1 + k2) / h1;
//         double coeff2 = 0.5 * (k1 + k2) / h2;

//         temp[center_index] = rhs[center_index] - (
//             - coeff1 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r-1][s_i_theta]) * s_x[s_i_r-1][s_i_theta] /* Left */  
//             - coeff2 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r+1][s_i_theta]) * s_x[s_i_r+1][s_i_theta] /* Right */   

//             - 0.25 * (s_art[s_i_r-1][s_i_theta] + s_art[s_i_r][s_i_theta-1]) * s_x[s_i_r-1][s_i_theta-1] /* Bottom Left */
//             + 0.25 * (s_art[s_i_r+1][s_i_theta] + s_art[s_i_r][s_i_theta-1]) * s_x[s_i_r+1][s_i_theta-1] /* Bottom Right */
//             + 0.25 * (s_art[s_i_r-1][s_i_theta] + s_art[s_i_r][s_i_theta+1]) * s_x[s_i_r-1][s_i_theta+1] /* Top Left */
//             - 0.25 * (s_art[s_i_r+1][s_i_theta] + s_art[s_i_r][s_i_theta+1]) * s_x[s_i_r+1][s_i_theta+1] /* Top Right */
//         );
//     }
//     else if(isOnInnerBoundary && !DirBC_Interior){

//         double coeff2 = 0.5 * (k1 + k2) / h2;

//         temp[center_index] = rhs[center_index] - (
//             - coeff2 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r+1][s_i_theta]) * s_x[s_i_r+1][s_i_theta] /* Right */

//             + 0.25 * (s_art[s_i_r+1][s_i_theta] + s_art[s_i_r][s_i_theta-1]) *  s_x[s_i_r+1][s_i_theta-1] /* Bottom Right */
//             - 0.25 * (s_art[s_i_r+1][s_i_theta] + s_art[s_i_r][s_i_theta+1]) * s_x[s_i_r+1][s_i_theta+1] /* Top Right */  
//         );
//     }
//     else if(isOnInnerBoundary && DirBC_Interior){
//         temp[center_index] = rhs[center_index];           
//     }
// }

// __global__ void applyAscOrtho_Radial_kernel(
//     double* x, double* rhs, double* temp,
//     PolarGrid* grid, bool DirBC_Interior,
//     int start_i_theta,
//     DomainGeometry* domain_geometry,
//     double* coeff_alpha_cache, double* coeff_beta_cache,
//     double* sin_theta_cache, double* cos_theta_cache
// ) 
// {
//     int i_r = grid->numberSmootherCircles() + blockIdx.x * 14 + threadIdx.x - 1;
//     int i_theta = blockIdx.y * 14 + threadIdx.y - 1;

//     if(i_r == -1 && !DirBC_Interior){
//         i_r = 0;
//         i_theta += grid->ntheta() / 2;
//     }
//     i_theta = grid->wrapThetaIndex(i_theta);

//     if (i_r < grid->numberSmootherCircles() - 1 || i_r >= grid->nr() || i_theta < 0 || i_theta >= grid->ntheta()) return;

//     __shared__ double s_arr[16][16];
//     __shared__ double s_att[16][16];
//     __shared__ double s_art[16][16];
//     __shared__ double s_x[16][16];

//     int s_i_r = threadIdx.x;
//     int s_i_theta = threadIdx.y;

//     int center_index = grid->index(i_r, i_theta);
//     s_x[s_i_r][s_i_theta] = x[center_index];

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
//     double art = (- Jtt * Jtr - Jrt * Jrr) * coeff_alpha / fabs(detDF);

//     s_arr[s_i_r][s_i_theta] = arr;
//     s_att[s_i_r][s_i_theta] = att;
//     s_art[s_i_r][s_i_theta] = art;

//     __syncthreads();

//     if(i_theta % 2 != start_i_theta) return;

//     if(i_r < grid->numberSmootherCircles() || i_r >= grid->nr() || i_theta < 0 || i_theta >= grid->ntheta()) return;
//     if(s_i_r == 0 || s_i_r == 15 || s_i_theta == 0 || s_i_theta == 15) return;

//     bool isOnOuterBoundary = (i_r == grid->nr()-1);
//     bool isNextToCircleSection = (i_r == grid->numberSmootherCircles());

//     double h1 = grid->radialSpacing(i_r-1);
//     double h2 = ((!isOnOuterBoundary) ? grid->radialSpacing(i_r) : 0.0);
//     double k1 = grid->angularSpacing(i_theta - 1);                                                          
//     double k2 = grid->angularSpacing(i_theta);

//     if(i_r >= grid->numberSmootherCircles() && i_r < grid->nr() - 1){
//         double coeff3 = 0.5 * (h1 + h2) / k1;
//         double coeff4 = 0.5 * (h1 + h2) / k2;
//         temp[center_index] = rhs[center_index] - (
//             - coeff3 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta-1]) * s_x[s_i_r][s_i_theta-1] /* Bottom */
//             - coeff4 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta+1]) * s_x[s_i_r][s_i_theta+1] /* Top */
//             - 0.25 * (s_art[s_i_r-1][s_i_theta] + s_art[s_i_r][s_i_theta-1]) * s_x[s_i_r-1][s_i_theta-1] /* Bottom Left */
//             + 0.25 * (s_art[s_i_r+1][s_i_theta] + s_art[s_i_r][s_i_theta-1]) * s_x[s_i_r+1][s_i_theta-1] /* Bottom Right */
//             + 0.25 * (s_art[s_i_r-1][s_i_theta] + s_art[s_i_r][s_i_theta+1]) * s_x[s_i_r-1][s_i_theta+1] /* Top Left */
//             - 0.25 * (s_art[s_i_r+1][s_i_theta] + s_art[s_i_r][s_i_theta+1]) * s_x[s_i_r+1][s_i_theta+1] /* Top Right */
//         );
//     }
//     if (i_r == grid->numberSmootherCircles()){
//         double coeff1 = 0.5 * (k1 + k2) / h1;
//         temp[center_index] -=  (
//             - coeff1 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r-1][s_i_theta]) * s_x[s_i_r-1][s_i_theta] /* Left */ 
//         );
//     }
//     if (i_r == grid->nr() - 2){
//         double coeff2 = 0.5 * (k1 + k2) / h2;
//         temp[center_index] -= (
//             /* "Right" is part of the radial Asc smoother matrices, */ 
//             /* but is shifted over to the rhs to make the radial Asc smoother matrices symmetric. */ 
//             /* Note that the circle Asc smoother matrices are symmetric by default. */
//             - coeff2 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r+1][s_i_theta]) * rhs[grid->index(i_r + 1, i_theta)] /* Right */ 
//         );        
//     }
//     if (i_r == grid->nr() - 1){
//         temp[center_index] = rhs[center_index];
//     }
// }

// void SmootherTakeGPU::smoothingInPlace(GPU_Vector<double>& x, const GPU_Vector<double>& rhs, GPU_Vector<double>& temp)
// {
//     const PolarGrid& grid = level_.grid();

//     assert(x.size() == grid.numberOfNodes());
//     assert(rhs.size() == grid.numberOfNodes());
//     assert(temp.size() == grid.numberOfNodes());

//     const GPU_Vector<double>& sin_theta_cache = level_.levelCache().GPU_sin_theta();
//     const GPU_Vector<double>& cos_theta_cache = level_.levelCache().GPU_cos_theta();

//     const GPU_Vector<double>& coeff_alpha_cache = level_.levelCache().GPU_coeff_alpha();
//     const GPU_Vector<double>& coeff_beta_cache = level_.levelCache().GPU_coeff_beta();

//     dim3 threadsPerBlock(1,1);
//     dim3 numBlocks(1,1);

//     DomainGeometry* device_domain_geometry;
//     cudaMalloc(&device_domain_geometry, sizeof(DomainGeometry));
//     cudaMemcpy(device_domain_geometry, &domain_geometry_, sizeof(DomainGeometry), cudaMemcpyHostToDevice);

//     /* We use precomputed DensityProfileCoefficients values. */
//     // DensityProfileCoefficients* device_density_profile;
//     // cudaMalloc(&device_density_profile, sizeof(DensityProfileCoefficients));
//     // cudaMemcpy(device_density_profile, &density_profile_coefficients_, sizeof(DensityProfileCoefficients), cudaMemcpyHostToDevice);
    

//     // threadsPerBlock = dim3(16, 16);
//     // numBlocks = dim3((grid.numberSmootherCircles() + 14 - 1) / 14, (grid.ntheta() + 14 - 1) / 14);
//     const int start_black_circles = (grid.numberSmootherCircles() % 2 == 0) ? 1 : 0;
//     // applyAscOrtho_Circle_kernel<<<numBlocks, threadsPerBlock>>>(
//     //     x.data(), rhs.data(), temp.data(), 
//     //     level_.device_grid(), DirBC_Interior_,
//     //     start_black_circles,
//     //     device_domain_geometry, 
//     //     coeff_alpha_cache.data(), coeff_beta_cache.data(), 
//     //     sin_theta_cache.data(), cos_theta_cache.data()
//     // );
//     // cudaDeviceSynchronize();


//     solve_Asc_Circle_kernel<<<1, (grid.numberSmootherCircles()+1)/2>>>(
//         x.data(), rhs.data(), temp.data(),
//         d_circle_tridiagonal_solver_,
//         level_.device_grid(), DirBC_Interior_,
//         start_black_circles
//     );
//     cudaDeviceSynchronize();


//     // threadsPerBlock = dim3(16, 16);
//     // numBlocks = dim3((grid.numberSmootherCircles() + 14 - 1) / 14, (grid.ntheta() + 14 - 1) / 14);
//     const int start_white_circles = (grid.numberSmootherCircles() % 2 == 0) ? 0 : 1;
//     // applyAscOrtho_Circle_kernel<<<numBlocks, threadsPerBlock>>>(
//     //     x.data(), rhs.data(), temp.data(), 
//     //     level_.device_grid(), DirBC_Interior_,
//     //     start_white_circles,
//     //     device_domain_geometry, 
//     //     coeff_alpha_cache.data(), coeff_beta_cache.data(), 
//     //     sin_theta_cache.data(), cos_theta_cache.data()
//     // );
//     // cudaDeviceSynchronize();

//     solve_Asc_Circle_kernel<<<1, (grid.numberSmootherCircles())/2>>>(
//         x.data(), rhs.data(), temp.data(),
//         d_circle_tridiagonal_solver_,
//         level_.device_grid(), DirBC_Interior_,
//         start_white_circles
//     );
//     cudaDeviceSynchronize();


//     // threadsPerBlock = dim3(16, 16);
//     // numBlocks = dim3((grid.lengthSmootherRadial() + 14 - 1) / 14, (grid.ntheta() + 14 - 1) / 14);
//     // const int start_black_radials = 0;
//     // applyAscOrtho_Radial_kernel<<<numBlocks, threadsPerBlock>>>(
//     //     x.data(), rhs.data(), temp.data(), 
//     //     level_.device_grid(), DirBC_Interior_,
//     //     start_black_radials,
//     //     device_domain_geometry, 
//     //     coeff_alpha_cache.data(), coeff_beta_cache.data(), 
//     //     sin_theta_cache.data(), cos_theta_cache.data()
//     // );
//     // cudaDeviceSynchronize();

//     // threadsPerBlock = dim3(16, 16);
//     // numBlocks = dim3((grid.lengthSmootherRadial() + 14 - 1) / 14, (grid.ntheta() + 14 - 1) / 14);
//     // const int start_white_radials = 1;
//     // applyAscOrtho_Radial_kernel<<<numBlocks, threadsPerBlock>>>(
//     //     x.data(), rhs.data(), temp.data(), 
//     //     level_.device_grid(), DirBC_Interior_,
//     //     start_white_radials,
//     //     device_domain_geometry, 
//     //     coeff_alpha_cache.data(), coeff_beta_cache.data(), 
//     //     sin_theta_cache.data(), cos_theta_cache.data()
//     // );
//     // cudaDeviceSynchronize();


//     /* We use precomputed DensityProfileCoefficients values. */
//     cudaFree(device_domain_geometry);
//     // cudaFree(device_density_profile);

//     // x = temp;

//     // std::cout<<"DEVICE"<<std::endl;
//     // std::cout<<x<<std::endl;
//     // std::cout<<temp<<std::endl;
// }




















































// #include "../../../include/Smoother/SmootherTakeGPU/smoother.h"

// #include "../../../include/LinearAlgebra/Vector/gpu_vector_operations.h"

// #include <cmath>

// __global__ void applyAscOrtho_BlackCircle_kernel(
// ) 
// {
//     if (threadIdx.x == 0 && blockIdx.x == 0 && threadIdx.y == 0 && blockIdx.y == 0) {
//         printf("applyAscOrtho_BlackCircle_kernel\n");
//     }
// }

// __global__ void applyAscOrtho_WhiteCircle_kernel(
// ) 
// {
//     unsigned long long waitTime = 50000000; // Example: 1 second (in microseconds)
//     unsigned long long startClock = clock64(); // Get the start time using CUDA clock

//     if (! (threadIdx.x == 0 && blockIdx.x == 0 && threadIdx.y == 0 && blockIdx.y == 0)) return;
//     while (clock64() - startClock < waitTime) {
//         // printf(" WAIT applyAscOrtho_WhiteCircle_kernel\n");
//         // Do nothing, just wait
//     }
//     double dwa = sin(323.0) * sin(threadIdx.x * 4.0);
//     for (int i = 0; i < 10000000; i++)
//     {
//        dwa += sin(323.0) * sin(i * 4.0);
//     }
    

//     if (threadIdx.x == 0 && blockIdx.x == 0 && threadIdx.y == 0 && blockIdx.y == 0) {
//         printf("applyAscOrtho_WhiteCircle_kernel\n");
//         printf("Thread (%d, %d) - dwa: %f\n", threadIdx.x, threadIdx.y, dwa);
//     }
// }

// __global__ void applyAscOrtho_BlackRadial_kernel(
// ) 
// {
//         if (threadIdx.x == 0 && blockIdx.x == 0 && threadIdx.y == 0 && blockIdx.y == 0) {
//         printf("applyAscOrtho_BlackRadial_kernel\n");
//     }
// }

// __global__ void applyAscOrtho_WhiteRadial_kernel(
// ) 
// {
//     if (threadIdx.x == 0 && blockIdx.x == 0 && threadIdx.y == 0 && blockIdx.y == 0) {
//         printf("applyAscOrtho_WhiteRadial_kernel\n");
//     }
// }

// __global__ void solveAsc_BlackCircle_kernel(
// ) 
// {
//         if (threadIdx.x == 0 && blockIdx.x == 0 && threadIdx.y == 0 && blockIdx.y == 0) {
//         printf("solveAsc_BlackCircle_kernel\n");
//     }
// }

// __global__ void solveAsc_WhiteCircle_kernel(
// ) 
// {
//         if (threadIdx.x == 0 && blockIdx.x == 0 && threadIdx.y == 0 && blockIdx.y == 0) {
//         printf("solveAsc_WhiteCircle_kernel\n");
//     }
// }

// __global__ void solveAsc_BlackRadial_kernel(
// ) 
// {
//         if (threadIdx.x == 0 && blockIdx.x == 0 && threadIdx.y == 0 && blockIdx.y == 0) {
//         printf("solveAsc_BlackRadial_kernel\n");
//     }
// }

// __global__ void solveAsc_WhiteRadial_kernel(
// ) 
// {
//         if (threadIdx.x == 0 && blockIdx.x == 0 && threadIdx.y == 0 && blockIdx.y == 0) {
//         printf("solveAsc_WhiteRadial_kernel\n");
//     }
// }

// __global__ void solveAsc_InteriorBlackCircle_kernel(
// ) 
// {
//         if (threadIdx.x == 0 && blockIdx.x == 0 && threadIdx.y == 0 && blockIdx.y == 0) {
//         printf("solveAsc_InteriorBlackCircle_kernel\n");
//     }
// }

// __global__ void solveAsc_InteriorWhiteCircle_kernel(
// ) 
// {
//         if (threadIdx.x == 0 && blockIdx.x == 0 && threadIdx.y == 0 && blockIdx.y == 0) {
//         printf("solveAsc_InteriorWhiteCircle_kernel\n");
//     }
// }











// void SmootherTakeGPU::smoothingInPlace(GPU_Vector<double>& x, const GPU_Vector<double>& rhs, GPU_Vector<double>& temp)
// {
//     const PolarGrid& grids = level_.grid();

//     assert(x.size() == grids.numberOfNodes());
//     assert(rhs.size() == grids.numberOfNodes());
//     assert(temp.size() == grids.numberOfNodes());

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
    
//     dim3 grid(16, 16);
//     dim3 block((grids.nr() + 14 - 1) / 14, (grids.ntheta() + 14 - 1) / 14);

//     // cudaStream_t stream1, stream2, stream3;
//     // cudaStreamCreate(&stream1);
//     // cudaStreamCreate(&stream2);
//     // cudaStreamCreate(&stream3);

//     // cudaEvent_t event1, event2, event3, event4, event5, event6, event7, event8, event9, event10;
//     // cudaEventCreate(&event1);
//     // cudaEventCreate(&event2);
//     // cudaEventCreate(&event3);
//     // cudaEventCreate(&event4);
//     // cudaEventCreate(&event5);
//     // cudaEventCreate(&event6);
//     // cudaEventCreate(&event7);
//     // cudaEventCreate(&event8);
//     // cudaEventCreate(&event9);
//     // cudaEventCreate(&event10);

//     // Event1: applyAscOrtho_BlackCircle_kernel. 
//     // Event2: solveAsc_BlackCircle_kernel. Waits for Event1.
//     // Event3: solveAsc_InteriorBlackCircle_kernel. Waits for Event1.
//     // Event4: solveAsc_WhiteCircle_kernel. Waits for Event2 and Event3.
//     // Event5: solveAsc_WhiteCircle_kernel. Waits for Event4.
//     // Event6: solveAsc_InteriorWhiteCircle_kernel. Waits for Event5.
//     // Event7: applyAscOrtho_BlackRadial_kernel. Waits for Event2 and Event3.
//     // Event8: solveAsc_BlackRadial_kernel. Waits for Event7.
//     // Event9: applyAscOrtho_WhiteRadial_kernel. Waits for Event8.
//     // Event10: solveAsc_WhiteRadial_kernel. Waits for Event9.

//     // /* Black Circle */
//     // applyAscOrtho_BlackCircle_kernel<<<numBlocks, threadsPerBlock, 0, stream1>>>(
//     //     x.data(), rhs.data(), temp.data(), 
//     //     level_.device_grid(), DirBC_Interior_,
//     //     device_domain_geometry, 
//     //     coeff_alpha_cache.data(), coeff_beta_cache.data(), 
//     //     sin_theta_cache.data(), cos_theta_cache.data()
//     // );
//     // cudaEventRecord(event1, stream1);

//     // cudaStreamWaitEvent(stream1, event1, 0);
//     // solveAsc_BlackCircle_kernel<<<numBlocks, threadsPerBlock, 0, stream1>>>(
//     //     x.data(), rhs.data(), temp.data(), 
//     //     level_.device_grid(), DirBC_Interior_,
//     //     device_domain_geometry, 
//     //     coeff_alpha_cache.data(), coeff_beta_cache.data(), 
//     //     sin_theta_cache.data(), cos_theta_cache.data()
//     // );
//     // cudaEventRecord(event2, stream1);

//     // cudaStreamWaitEvent(stream2, event1, 0);
//     // solveAsc_InteriorBlackCircle_kernel<<<numBlocks, threadsPerBlock, 0, stream2>>>(
//     //     x.data(), rhs.data(), temp.data(), 
//     //     level_.device_grid(), DirBC_Interior_,
//     //     device_domain_geometry, 
//     //     coeff_alpha_cache.data(), coeff_beta_cache.data(), 
//     //     sin_theta_cache.data(), cos_theta_cache.data()
//     // );
//     // cudaEventRecord(event3, stream2);

//     // /* White Circle */
//     // cudaStreamWaitEvent(stream1, event2, 0);
//     // cudaStreamWaitEvent(stream1, event3, 0);
//     // applyAscOrtho_WhiteCircle_kernel<<<numBlocks, threadsPerBlock, 0, stream1>>>(
//     //     x.data(), rhs.data(), temp.data(), 
//     //     level_.device_grid(), DirBC_Interior_,
//     //     device_domain_geometry, 
//     //     coeff_alpha_cache.data(), coeff_beta_cache.data(), 
//     //     sin_theta_cache.data(), cos_theta_cache.data()
//     // );
//     // cudaEventRecord(event4, stream1);

//     // cudaStreamWaitEvent(stream1, event4, 0);
//     // solveAsc_WhiteCircle_kernel<<<numBlocks, threadsPerBlock, 0, stream1>>>(
//     //     x.data(), rhs.data(), temp.data(), 
//     //     level_.device_grid(), DirBC_Interior_,
//     //     device_domain_geometry, 
//     //     coeff_alpha_cache.data(), coeff_beta_cache.data(), 
//     //     sin_theta_cache.data(), cos_theta_cache.data()
//     // );
//     // cudaEventRecord(event5, stream1);

//     // cudaStreamWaitEvent(stream2, event4, 0);
//     // solveAsc_InteriorWhiteCircle_kernel<<<numBlocks, threadsPerBlock, 0, stream2>>>(
//     //     x.data(), rhs.data(), temp.data(), 
//     //     level_.device_grid(), DirBC_Interior_,
//     //     device_domain_geometry, 
//     //     coeff_alpha_cache.data(), coeff_beta_cache.data(), 
//     //     sin_theta_cache.data(), cos_theta_cache.data()
//     // );
//     // cudaEventRecord(event6, stream2);

//     // /* Black Radial */
//     // cudaStreamWaitEvent(stream3, event2, 0);
//     // cudaStreamWaitEvent(stream3, event3, 0);
//     // applyAscOrtho_BlackRadial_kernel<<<numBlocks, threadsPerBlock, 0, stream3>>>(
//     //     x.data(), rhs.data(), temp.data(), 
//     //     level_.device_grid(), DirBC_Interior_,
//     //     device_domain_geometry, 
//     //     coeff_alpha_cache.data(), coeff_beta_cache.data(), 
//     //     sin_theta_cache.data(), cos_theta_cache.data()
//     // );
//     // cudaEventRecord(event7, stream3);

//     // cudaStreamWaitEvent(stream3, event8, 0);
//     // solveAsc_BlackRadial_kernel<<<numBlocks, threadsPerBlock, 0, stream3>>>(
//     //     x.data(), rhs.data(), temp.data(), 
//     //     level_.device_grid(), DirBC_Interior_,
//     //     device_domain_geometry, 
//     //     coeff_alpha_cache.data(), coeff_beta_cache.data(), 
//     //     sin_theta_cache.data(), cos_theta_cache.data()
//     // );
//     // cudaEventRecord(event8, stream3);

//     // /* White Radial */
//     // cudaStreamWaitEvent(stream3, event8, 0);
//     // applyAscOrtho_WhiteRadial_kernel<<<numBlocks, threadsPerBlock, 0, stream3>>>(
//     //     x.data(), rhs.data(), temp.data(), 
//     //     level_.device_grid(), DirBC_Interior_,
//     //     device_domain_geometry, 
//     //     coeff_alpha_cache.data(), coeff_beta_cache.data(), 
//     //     sin_theta_cache.data(), cos_theta_cache.data()
//     // );
//     // cudaEventRecord(event8, stream3);

//     // cudaStreamWaitEvent(stream3, event9, 0);
//     // solveAsc_WhiteRadial_kernel<<<numBlocks, threadsPerBlock, 0, stream3>>>(
//     //     x.data(), rhs.data(), temp.data(), 
//     //     level_.device_grid(), DirBC_Interior_,
//     //     device_domain_geometry, 
//     //     coeff_alpha_cache.data(), coeff_beta_cache.data(), 
//     //     sin_theta_cache.data(), cos_theta_cache.data()
//     // );
//     // cudaEventRecord(event10, stream3);







// //     cudaStream_t streams[10];
// //     cudaEvent_t events[10];

// //     // Create streams
// //     for (int i = 0; i < 10; i++) {
// //         cudaStreamCreate(&streams[i]);
// //     }

// //     // Create events
// //     for (int i = 0; i < 10; i++) {
// //         cudaEventCreate(&events[i]);
// //     }

// //    dim3 gridSize(16, 16);
// //     dim3 blockSize((grid.nr() + 14 - 1) / 14, (grid.ntheta() + 14 - 1) / 14);

// //     // Launch Event1: applyAscOrtho_BlackCircle_kernel
// //     applyAscOrtho_BlackCircle_kernel<<<gridSize, blockSize, 0, streams[0]>>>();
// //     cudaEventRecord(events[0], streams[0]);

// //     // Launch Event2: solveAsc_BlackCircle_kernel, waits for Event1
// //     solveAsc_BlackCircle_kernel<<<gridSize, blockSize, 0, streams[1]>>>();
// //     cudaStreamWaitEvent(streams[1], events[0], 0);  // Event2 waits for Event1
// //     cudaEventRecord(events[1], streams[1]);

// //     // Launch Event3: solveAsc_InteriorBlackCircle_kernel, waits for Event1
// //     solveAsc_InteriorBlackCircle_kernel<<<gridSize, blockSize, 0, streams[2]>>>();
// //     cudaStreamWaitEvent(streams[2], events[0], 0);  // Event3 waits for Event1
// //     cudaEventRecord(events[2], streams[2]);

// //     // Launch Event4: applyAscOrtho_WhiteCircle_kernel, waits for Event2 and Event3
// //     applyAscOrtho_WhiteCircle_kernel<<<gridSize, blockSize, 0, streams[3]>>>();
// //     cudaStreamWaitEvent(streams[3], events[1], 0);  // Event4 waits for Event2
// //     cudaStreamWaitEvent(streams[3], events[2], 0);  // Event4 waits for Event3
// //     cudaEventRecord(events[3], streams[3]);

// //     // Launch Event5: solveAsc_WhiteCircle_kernel, waits for Event4
// //     solveAsc_WhiteCircle_kernel<<<gridSize, blockSize, 0, streams[4]>>>();
// //     cudaStreamWaitEvent(streams[4], events[3], 0);  // Event5 waits for Event4
// //     cudaEventRecord(events[4], streams[4]);

// //     // Launch Event6: solveAsc_InteriorWhiteCircle_kernel, waits for Event5
// //     solveAsc_InteriorWhiteCircle_kernel<<<gridSize, blockSize, 0, streams[5]>>>();
// //     cudaStreamWaitEvent(streams[5], events[4], 0);  // Event6 waits for Event5
// //     cudaEventRecord(events[5], streams[5]);

// //     // Launch Event7: applyAscOrtho_BlackRadial_kernel, waits for Event2 and Event3
// //     applyAscOrtho_BlackRadial_kernel<<<gridSize, blockSize, 0, streams[6]>>>();
// //     cudaStreamWaitEvent(streams[6], events[1], 0);  // Event7 waits for Event2
// //     cudaStreamWaitEvent(streams[6], events[2], 0);  // Event7 waits for Event3
// //     cudaEventRecord(events[6], streams[6]);

// //     // Launch Event8: solveAsc_BlackRadial_kernel, waits for Event7
// //     solveAsc_BlackRadial_kernel<<<gridSize, blockSize, 0, streams[7]>>>();
// //     cudaStreamWaitEvent(streams[7], events[6], 0);  // Event8 waits for Event7
// //     cudaEventRecord(events[7], streams[7]);

// //     // Launch Event9: applyAscOrtho_WhiteRadial_kernel, waits for Event8
// //     applyAscOrtho_WhiteRadial_kernel<<<gridSize, blockSize, 0, streams[8]>>>();
// //     cudaStreamWaitEvent(streams[8], events[7], 0);  // Event9 waits for Event8
// //     cudaEventRecord(events[8], streams[8]);

// //     // Launch Event10: solveAsc_WhiteRadial_kernel, waits for Event9
// //     solveAsc_WhiteRadial_kernel<<<gridSize, blockSize, 0, streams[9]>>>();
// //     cudaStreamWaitEvent(streams[9], events[8], 0);  // Event10 waits for Event9
// //     cudaEventRecord(events[9], streams[9]);

// //     // Synchronize all streams (optional)
// //     for (int i = 0; i < 10; i++) {
// //         cudaStreamSynchronize(streams[i]);
// //     }

// //     // Clean up
// //     for (int i = 0; i < 10; i++) {
// //         cudaStreamDestroy(streams[i]);
// //         cudaEventDestroy(events[i]);
// //     }




// //    dim3 grid(16, 16);
// //     dim3 block((1000 + 14 - 1) / 14, (1000 + 14 - 1) / 14);



// // Create CUDA events
// cudaEvent_t event1, event2, event3, event4, event5, event6, event7, event8, event9, event10;
// cudaEventCreate(&event1);
// cudaEventCreate(&event2);
// cudaEventCreate(&event3);
// cudaEventCreate(&event4);
// cudaEventCreate(&event5);
// cudaEventCreate(&event6);
// cudaEventCreate(&event7);
// cudaEventCreate(&event8);
// cudaEventCreate(&event9);
// cudaEventCreate(&event10);

// // Create CUDA streams
// cudaStream_t stream1, stream2, stream3, stream4, stream5, stream6, stream7, stream8, stream9, stream10;
// cudaStreamCreate(&stream1);
// cudaStreamCreate(&stream2);
// cudaStreamCreate(&stream3);
// cudaStreamCreate(&stream4);
// cudaStreamCreate(&stream5);
// cudaStreamCreate(&stream6);
// cudaStreamCreate(&stream7);
// cudaStreamCreate(&stream8);
// cudaStreamCreate(&stream9);
// cudaStreamCreate(&stream10);

// // Launch Event1 kernel in stream1 (No dependencies)
// applyAscOrtho_BlackCircle_kernel<<<grid, block, 0, stream1>>>();
// cudaEventRecord(event1, stream1);

// // Launch Event2 kernel (depends on Event1) in stream2
// cudaStreamWaitEvent(stream2, event1, 0);
// solveAsc_BlackCircle_kernel<<<grid, block, 0, stream2>>>();
// cudaEventRecord(event2, stream2);

// // Launch Event3 kernel (depends on Event1) in stream3
// cudaStreamWaitEvent(stream3, event1, 0);
// solveAsc_InteriorBlackCircle_kernel<<<grid, block, 0, stream3>>>();
// cudaEventRecord(event3, stream3);

// // Launch Event4 kernel (depends on Event2 and Event3) in stream4
// cudaStreamWaitEvent(stream4, event2, 0);
// cudaStreamWaitEvent(stream4, event3, 0);
// applyAscOrtho_WhiteCircle_kernel<<<grid, block, 0, stream4>>>();
// cudaEventRecord(event4, stream4);

// // Launch Event5 kernel (depends on Event4) in stream5
// cudaStreamWaitEvent(stream5, event4, 0);
// solveAsc_WhiteCircle_kernel<<<grid, block, 0, stream5>>>();
// cudaEventRecord(event5, stream5);

// // Launch Event6 kernel (depends on Event5) in stream6
// cudaStreamWaitEvent(stream6, event5, 0);
// solveAsc_InteriorWhiteCircle_kernel<<<grid, block, 0, stream6>>>();
// cudaEventRecord(event6, stream6);

// // Launch Event7 kernel (depends on Event2 and Event3) in stream7
// cudaStreamWaitEvent(stream7, event2, 0);
// cudaStreamWaitEvent(stream7, event3, 0);
// applyAscOrtho_BlackRadial_kernel<<<grid, block, 0, stream7>>>();
// cudaEventRecord(event7, stream7);

// // Launch Event8 kernel (depends on Event7) in stream8
// cudaStreamWaitEvent(stream8, event7, 0);
// solveAsc_BlackRadial_kernel<<<grid, block, 0, stream8>>>();
// cudaEventRecord(event8, stream8);

// // Launch Event9 kernel (depends on Event8) in stream9
// cudaStreamWaitEvent(stream9, event8, 0);
// applyAscOrtho_WhiteRadial_kernel<<<grid, block, 0, stream9>>>();
// cudaEventRecord(event9, stream9);

// // Launch Event10 kernel (depends on Event9) in stream10
// cudaStreamWaitEvent(stream10, event9, 0);
// solveAsc_WhiteRadial_kernel<<<grid, block, 0, stream10>>>();
// cudaEventRecord(event10, stream10);

// // Wait for the final event to complete
// cudaEventSynchronize(event10);

// // Clean up resources
// cudaEventDestroy(event1);
// cudaEventDestroy(event2);
// cudaEventDestroy(event3);
// cudaEventDestroy(event4);
// cudaEventDestroy(event5);
// cudaEventDestroy(event6);
// cudaEventDestroy(event7);
// cudaEventDestroy(event8);
// cudaEventDestroy(event9);
// cudaEventDestroy(event10);

// cudaStreamDestroy(stream1);
// cudaStreamDestroy(stream2);
// cudaStreamDestroy(stream3);
// cudaStreamDestroy(stream4);
// cudaStreamDestroy(stream5);
// cudaStreamDestroy(stream6);
// cudaStreamDestroy(stream7);
// cudaStreamDestroy(stream8);
// cudaStreamDestroy(stream9);
// cudaStreamDestroy(stream10);





//     cudaDeviceSynchronize();

//     /* We use precomputed DensityProfileCoefficients values. */
//     cudaFree(device_domain_geometry);
//     // cudaFree(device_density_profile);
// }