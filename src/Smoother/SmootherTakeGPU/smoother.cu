#include "../../../include/Smoother/SmootherTakeGPU/smoother.h"

#include <chrono>

SmootherTakeGPU::SmootherTakeGPU(const Level& level, const DomainGeometry& domain_geometry,
                   const DensityProfileCoefficients& density_profile_coefficients, bool DirBC_Interior)
    : level_(level)
    , domain_geometry_(domain_geometry)
    , density_profile_coefficients_(density_profile_coefficients)
    , DirBC_Interior_(DirBC_Interior)
    , circle_main_diagonals_(nullptr)
    , circle_lower_diagonals_(nullptr)
    , circle_upper_diagonals_(nullptr)
    , radial_main_diagonals_(nullptr)
    , radial_lower_diagonals_(nullptr)
    , radial_upper_diagonals_(nullptr)
    , pBuffer_(nullptr)
{
    const PolarGrid& grid = level.grid();

    int nr = grid.nr();
    int ntheta = grid.ntheta();
    int number_smoother_circles = grid.numberSmootherCircles();
    int length_smoother_radial = grid.lengthSmootherRadial();

    int circle_batch_count = number_smoother_circles;
    int circle_m = ntheta;
    cudaMalloc(&circle_lower_diagonals_, circle_m * circle_batch_count * sizeof(double));
    cudaMalloc(&circle_main_diagonals_, circle_m * circle_batch_count * sizeof(double));
    cudaMalloc(&circle_upper_diagonals_, circle_m * circle_batch_count * sizeof(double));
    cudaMemset(circle_lower_diagonals_, 0, circle_m * circle_batch_count * sizeof(double));
    cudaMemset(circle_main_diagonals_, 0, circle_m * circle_batch_count * sizeof(double));
    cudaMemset(circle_upper_diagonals_, 0, circle_m * circle_batch_count * sizeof(double));
    cudaMalloc(&sherman_morrison_gammas_, circle_batch_count * sizeof(double));

    int radial_batch_count = ntheta;
    int radial_m = length_smoother_radial;
    cudaMalloc(&radial_lower_diagonals_, radial_m * radial_batch_count * sizeof(double));
    cudaMalloc(&radial_main_diagonals_, radial_m * radial_batch_count * sizeof(double));
    cudaMalloc(&radial_upper_diagonals_, radial_m * radial_batch_count * sizeof(double));
    cudaMemset(radial_lower_diagonals_, 0, radial_m * radial_batch_count * sizeof(double));
    cudaMemset(radial_main_diagonals_, 0, radial_m * radial_batch_count * sizeof(double));
    cudaMemset(radial_upper_diagonals_, 0, radial_m * radial_batch_count * sizeof(double));
    

    cusparseCreate(&handle_);

    size_t black_circle_pBufferSizeInBytes;
    cusparseDgtsv2StridedBatch_bufferSizeExt(
        handle_, circle_m, 
        circle_lower_diagonals_, circle_main_diagonals_, circle_upper_diagonals_, 
        nullptr, 
        (circle_batch_count + 1) / 2, 
        2 * circle_m, 
        &black_circle_pBufferSizeInBytes
    );
    size_t white_circle_pBufferSizeInBytes;
    cusparseDgtsv2StridedBatch_bufferSizeExt(
        handle_, circle_m, 
        circle_lower_diagonals_, circle_main_diagonals_, circle_upper_diagonals_, 
        nullptr, 
        (circle_batch_count) / 2, 
        2 * circle_m, 
        &white_circle_pBufferSizeInBytes
    );

    size_t black_radial_pBufferSizeInBytes;
    cusparseDgtsv2StridedBatch_bufferSizeExt(
        handle_, radial_m, 
        radial_lower_diagonals_, radial_main_diagonals_, radial_upper_diagonals_, 
        nullptr, 
        (radial_batch_count) / 2, 
        2 * radial_m, 
        &black_radial_pBufferSizeInBytes
    );
    size_t white_radial_pBufferSizeInBytes;
    cusparseDgtsv2StridedBatch_bufferSizeExt(
        handle_, radial_m, 
        radial_lower_diagonals_, radial_main_diagonals_, radial_upper_diagonals_, 
        nullptr, 
        (radial_batch_count) / 2, 
        2 * radial_m, 
        &white_radial_pBufferSizeInBytes
    );

    size_t max_pBufferSizeInBytes = std::max({
        black_circle_pBufferSizeInBytes, 
        white_circle_pBufferSizeInBytes, 
        black_radial_pBufferSizeInBytes, 
        white_radial_pBufferSizeInBytes
    });

    cudaMalloc(&pBuffer_, max_pBufferSizeInBytes);

    buildAscMatrices();
    shermannMorrisonAdjustment();
}


SmootherTakeGPU::~SmootherTakeGPU() {

    // Destroy cusparse handle
    cusparseDestroy(handle_);

    // Free allocated device memory for circle diagonals
    if (circle_lower_diagonals_) {
        cudaFree(circle_lower_diagonals_);
        circle_lower_diagonals_ = nullptr;
    }
    if (circle_main_diagonals_) {
        cudaFree(circle_main_diagonals_);
        circle_main_diagonals_ = nullptr;
    }
    if (circle_upper_diagonals_) {
        cudaFree(circle_upper_diagonals_);
        circle_upper_diagonals_ = nullptr;
    }
    if (sherman_morrison_gammas_) {
        cudaFree(sherman_morrison_gammas_);
        sherman_morrison_gammas_ = nullptr;
    }

    // Free allocated device memory for radial diagonals
    if (radial_lower_diagonals_) {
        cudaFree(radial_lower_diagonals_);
        radial_lower_diagonals_ = nullptr;
    }
    if (radial_main_diagonals_) {
        cudaFree(radial_main_diagonals_);
        radial_main_diagonals_ = nullptr;
    }
    if (radial_upper_diagonals_) {
        cudaFree(radial_upper_diagonals_);
        radial_upper_diagonals_ = nullptr;
    }

    // Free the buffer memory
    if (pBuffer_) {
        cudaFree(pBuffer_);
        pBuffer_ = nullptr;
    }
}



__global__ void build_Asc_kernel(
    double* circle_lower_diagonals, double* circle_main_diagonals, double* circle_upper_diagonals,
    double* radial_lower_diagonals, double* radial_main_diagonals, double* radial_upper_diagonals,
    PolarGrid* grid, bool DirBC_Interior,
    DomainGeometry* domain_geometry,
    double* coeff_alpha_cache, double* coeff_beta_cache,
    double* sin_theta_cache, double* cos_theta_cache
)
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

    int s_i_r = threadIdx.x;
    int s_i_theta = threadIdx.y;

    int center_index = grid->index(i_r, i_theta);
 
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

    s_detDF[s_i_r][s_i_theta] = detDF;
    s_arr[s_i_r][s_i_theta] = arr;
    s_att[s_i_r][s_i_theta] = att;

    __syncthreads();

    if (i_r < 0 || i_r >= grid->nr() || i_theta < 0 || i_theta >= grid->ntheta()) return;
    if(s_i_r == 0 || s_i_r == 15 || s_i_theta == 0 || s_i_theta == 15) return;

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

    if(i_r == 0){
        if(DirBC_Interior){

        }
        else{

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
    /* ------------------------------------------- */
    /* Radial Section: Node next to outer boundary */
    /* ------------------------------------------- */
    else if(i_r == grid->nr()-2){
        int center_index = i_r - numberSmootherCircles;
        int left_index   = i_r - numberSmootherCircles - 1;

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
    }
    /* ------------------------------------------ */
    /* Radial Section: Node on the outer boundary */
    /* ------------------------------------------ */
    else if(i_r == grid->nr()-1){

        int center_index = i_r - numberSmootherCircles;

        row = center_index;
        column = center_index;
        value = 1.0;
        if (row == column) radial_main_diagonals[i_theta * radial_m + row] = value;
        else if (row == column + 1) radial_lower_diagonals[i_theta * radial_m + row] = value;
        else if (row == column - 1) radial_upper_diagonals[i_theta * radial_m + row] = value;
        else if (row == 0 && column == radial_m - 1) radial_lower_diagonals[i_theta * radial_m + row] = value;
        else if (row == radial_m - 1 && column == 0) radial_upper_diagonals[i_theta * radial_m + row] = value;
    }
}


void SmootherTakeGPU::buildAscMatrices()
{
    std::cout<<"Start Building GPU Asc Matrices!"<<std::endl;
    std::cout<<""<<std::endl;

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
    
    dim3 threadsPerBlock(16, 16);
    dim3 numBlocks((grid.nr() + 14 - 1) / 14, (grid.ntheta() + 14 - 1) / 14);
    build_Asc_kernel<<<numBlocks, threadsPerBlock>>>(
        circle_lower_diagonals_, circle_main_diagonals_, circle_upper_diagonals_,
        radial_lower_diagonals_, radial_main_diagonals_, radial_upper_diagonals_,
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



__global__ void shermannMorrisonAdjustment_kernel(
    double* circle_lower_diagonals, double* circle_main_diagonals, double* circle_upper_diagonals,
    double* sherman_morrison_gammas,
    PolarGrid* grid
)
{
    int i_r = blockIdx.x * blockDim.x + threadIdx.x;

    if(i_r >= grid->numberSmootherCircles()) return;

    int circle_m = grid->ntheta();

    if(i_r == 0){
        sherman_morrison_gammas[i_r] = 0.0; // Unused
    } 
    else if(i_r > 0 && i_r < grid->numberSmootherCircles())
    {
        double local_gamma = - circle_main_diagonals[i_r * circle_m + 0];
        sherman_morrison_gammas[i_r] = local_gamma;

        circle_main_diagonals[i_r * circle_m + 0] -= local_gamma;

        circle_main_diagonals[i_r * circle_m + circle_m-1] -= 
            circle_lower_diagonals[i_r * circle_m + 0] * circle_upper_diagonals[i_r * circle_m + circle_m-1] / local_gamma;
    }
}



void SmootherTakeGPU::shermannMorrisonAdjustment()
{
    const PolarGrid& grid = level_.grid();

    int blockSize = 256; 
    int numBlocks = (grid.numberSmootherCircles() + blockSize - 1) / blockSize;

    shermannMorrisonAdjustment_kernel<<<numBlocks, blockSize>>>(
        circle_lower_diagonals_, circle_main_diagonals_, circle_upper_diagonals_,
        sherman_morrison_gammas_,
        level_.device_grid()
    );
    cudaDeviceSynchronize();
}













__global__ void applyAscOrtho_Circle_kernel(
    double* x, double* rhs, double* temp,
    double* circle_lower_diagonals, double* circle_main_diagonals, double* circle_upper_diagonals,
    double* sherman_morrison_gammas,
    PolarGrid* grid, bool DirBC_Interior,
    int start_i_r,
    DomainGeometry* domain_geometry,
    double* coeff_alpha_cache, double* coeff_beta_cache,
    double* sin_theta_cache, double* cos_theta_cache
) 
{
    int i_r = blockIdx.x * 14 + threadIdx.x - 1;
    int i_theta = blockIdx.y * 14 + threadIdx.y - 1;

    if(i_r == -1 && !DirBC_Interior){
        i_r = 0;
        i_theta += grid->ntheta() / 2;
    }
    i_theta = grid->wrapThetaIndex(i_theta);

    if (i_r < 0 || i_r > grid->numberSmootherCircles() || i_theta < 0 || i_theta >= grid->ntheta()) return;

    __shared__ double s_arr[16][16];
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
    double art = (- Jtt * Jtr - Jrt * Jrr) * coeff_alpha / fabs(detDF);

    s_arr[s_i_r][s_i_theta] = arr;
    s_art[s_i_r][s_i_theta] = art;

    __syncthreads();

    if(i_r % 2 != start_i_r) return;

    if (i_r < 0 || i_r >= grid->numberSmootherCircles() || i_theta < 0 || i_theta >= grid->ntheta()) return;
    if(s_i_r == 0 || s_i_r == 15 || s_i_theta == 0 || s_i_theta == 15) return;

    if(i_r > 0){
        int circle_m = grid->ntheta();
        if(i_theta == 0){
            temp[i_r * circle_m + i_theta] = sherman_morrison_gammas[i_r];
        }
        else if(i_theta > 0 && i_theta < grid->ntheta()-1){
            temp[i_r * circle_m + i_theta] = 0.0;
        }
        else if(i_theta == grid->ntheta()-1){
            temp[i_r * circle_m + i_theta] = circle_upper_diagonals[i_r * circle_m + circle_m-1];
        }
    }

    bool isOnInnerBoundary = (i_r == 0);

    double h1 = DirBC_Interior ? 
        ((!isOnInnerBoundary) ? grid->radialSpacing(i_r - 1) : 0.0) :
        ((!isOnInnerBoundary) ? grid->radialSpacing(i_r - 1) : 2.0 * grid->radius(0));
    double h2 = grid->radialSpacing(i_r);
    double k1 = grid->angularSpacing(i_theta - 1);                                                          
    double k2 = grid->angularSpacing(i_theta);

    if (!isOnInnerBoundary) {   

        double coeff1 = 0.5 * (k1 + k2) / h1;
        double coeff2 = 0.5 * (k1 + k2) / h2;

        x[center_index] = rhs[center_index] - (
            - coeff1 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r-1][s_i_theta]) * s_x[s_i_r-1][s_i_theta] /* Left */  
            - coeff2 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r+1][s_i_theta]) * s_x[s_i_r+1][s_i_theta] /* Right */   

            - 0.25 * (s_art[s_i_r-1][s_i_theta] + s_art[s_i_r][s_i_theta-1]) * s_x[s_i_r-1][s_i_theta-1] /* Bottom Left */
            + 0.25 * (s_art[s_i_r+1][s_i_theta] + s_art[s_i_r][s_i_theta-1]) * s_x[s_i_r+1][s_i_theta-1] /* Bottom Right */
            + 0.25 * (s_art[s_i_r-1][s_i_theta] + s_art[s_i_r][s_i_theta+1]) * s_x[s_i_r-1][s_i_theta+1] /* Top Left */
            - 0.25 * (s_art[s_i_r+1][s_i_theta] + s_art[s_i_r][s_i_theta+1]) * s_x[s_i_r+1][s_i_theta+1] /* Top Right */
        );
    }
    else if(isOnInnerBoundary && !DirBC_Interior){

        double coeff2 = 0.5 * (k1 + k2) / h2;

        x[center_index] = rhs[center_index] - (
            - coeff2 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r+1][s_i_theta]) * s_x[s_i_r+1][s_i_theta] /* Right */

            + 0.25 * (s_art[s_i_r+1][s_i_theta] + s_art[s_i_r][s_i_theta-1]) *  s_x[s_i_r+1][s_i_theta-1] /* Bottom Right */
            - 0.25 * (s_art[s_i_r+1][s_i_theta] + s_art[s_i_r][s_i_theta+1]) * s_x[s_i_r+1][s_i_theta+1] /* Top Right */  
        );
    }
    else if(isOnInnerBoundary && DirBC_Interior){
        x[center_index] = rhs[center_index];           
    }
}

__global__ void applyAscOrtho_Radial_kernel(
    double* x, double* rhs,
    PolarGrid* grid, bool DirBC_Interior,
    int start_i_theta,
    DomainGeometry* domain_geometry,
    double* coeff_alpha_cache, double* coeff_beta_cache,
    double* sin_theta_cache, double* cos_theta_cache
) 
{
    int i_r = grid->numberSmootherCircles() + blockIdx.x * 14 + threadIdx.x - 1;
    int i_theta = blockIdx.y * 14 + threadIdx.y - 1;

    if(i_r == -1 && !DirBC_Interior){
        i_r = 0;
        i_theta += grid->ntheta() / 2;
    }
    i_theta = grid->wrapThetaIndex(i_theta);

    if (i_r < grid->numberSmootherCircles() - 1 || i_r >= grid->nr() || i_theta < 0 || i_theta >= grid->ntheta()) return;

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

    s_arr[s_i_r][s_i_theta] = arr;
    s_att[s_i_r][s_i_theta] = att;
    s_art[s_i_r][s_i_theta] = art;

    __syncthreads();

    if(i_theta % 2 != start_i_theta) return;

    if(i_r < grid->numberSmootherCircles() || i_r >= grid->nr() || i_theta < 0 || i_theta >= grid->ntheta()) return;
    if(s_i_r == 0 || s_i_r == 15 || s_i_theta == 0 || s_i_theta == 15) return;

    bool isOnOuterBoundary = (i_r == grid->nr()-1);
    bool isNextToCircleSection = (i_r == grid->numberSmootherCircles());

    double h1 = grid->radialSpacing(i_r-1);
    double h2 = ((!isOnOuterBoundary) ? grid->radialSpacing(i_r) : 0.0);
    double k1 = grid->angularSpacing(i_theta - 1);                                                          
    double k2 = grid->angularSpacing(i_theta);

    if(i_r >= grid->numberSmootherCircles() && i_r < grid->nr() - 1){
        double coeff3 = 0.5 * (h1 + h2) / k1;
        double coeff4 = 0.5 * (h1 + h2) / k2;
        x[center_index] = rhs[center_index] - (
            - coeff3 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta-1]) * s_x[s_i_r][s_i_theta-1] /* Bottom */
            - coeff4 * (s_att[s_i_r][s_i_theta] + s_att[s_i_r][s_i_theta+1]) * s_x[s_i_r][s_i_theta+1] /* Top */
            - 0.25 * (s_art[s_i_r-1][s_i_theta] + s_art[s_i_r][s_i_theta-1]) * s_x[s_i_r-1][s_i_theta-1] /* Bottom Left */
            + 0.25 * (s_art[s_i_r+1][s_i_theta] + s_art[s_i_r][s_i_theta-1]) * s_x[s_i_r+1][s_i_theta-1] /* Bottom Right */
            + 0.25 * (s_art[s_i_r-1][s_i_theta] + s_art[s_i_r][s_i_theta+1]) * s_x[s_i_r-1][s_i_theta+1] /* Top Left */
            - 0.25 * (s_art[s_i_r+1][s_i_theta] + s_art[s_i_r][s_i_theta+1]) * s_x[s_i_r+1][s_i_theta+1] /* Top Right */
        );
    }
    if (i_r == grid->numberSmootherCircles()){
        double coeff1 = 0.5 * (k1 + k2) / h1;
        x[center_index] -=  (
            - coeff1 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r-1][s_i_theta]) * s_x[s_i_r-1][s_i_theta] /* Left */ 
        );
    }
    if (i_r == grid->nr() - 2){
        double coeff2 = 0.5 * (k1 + k2) / h2;
        x[center_index] -= (
            /* "Right" is part of the radial Asc smoother matrices, */ 
            /* but is shifted over to the rhs to make the radial Asc smoother matrices symmetric. */ 
            /* Note that the circle Asc smoother matrices are symmetric by default. */
            - coeff2 * (s_arr[s_i_r][s_i_theta] + s_arr[s_i_r+1][s_i_theta]) * rhs[grid->index(i_r + 1, i_theta)] /* Right */ 
        );        
    }
    if (i_r == grid->nr() - 1){
        x[center_index] = rhs[center_index];
    }
}






__global__ void combine_circle_solutions(
    double* x, double* temp,
    double* circle_lower_diagonals, double* circle_main_diagonals, double* circle_upper_diagonals,
    double* sherman_morrison_gammas,
    int start_solver_i_r,
    PolarGrid* grid
) 
{
    int i_r = start_solver_i_r + 2 * (blockIdx.x * 16 + threadIdx.x);
    int i_theta = blockIdx.y * 16 + threadIdx.y;

    if(i_r < 0 || i_r >= grid->numberSmootherCircles() || i_theta < 0 || i_theta >= grid->ntheta()) return;

    if(i_r > 0){
        int circle_m = grid->ntheta();
        const double dot_product_x_v = x[circle_m * i_r + 0] + circle_lower_diagonals[circle_m * i_r + 0] / sherman_morrison_gammas[i_r] * x[circle_m * i_r + circle_m - 1];
        const double dot_product_u_v = temp[circle_m * i_r + 0] + circle_lower_diagonals[circle_m * i_r + 0] / sherman_morrison_gammas[i_r]  * temp[circle_m * i_r + circle_m - 1];
        const double factor          = dot_product_x_v / (1.0 + dot_product_u_v);

        x[circle_m * i_r + i_theta] -= factor * temp[circle_m * i_r + i_theta];

    }

}




void SmootherTakeGPU::smoothingInPlace(GPU_Vector<double>& x, const GPU_Vector<double>& rhs, GPU_Vector<double>& temp)
{
    std::cout<<"Start GPU Smoother!"<<std::endl;
    std::cout<<""<<std::endl;

    const PolarGrid& grid = level_.grid();

    assert(x.size() == grid.numberOfNodes());
    assert(rhs.size() == grid.numberOfNodes());
    assert(temp.size() == grid.numberOfNodes());

    const GPU_Vector<double>& sin_theta_cache = level_.levelCache().GPU_sin_theta();
    const GPU_Vector<double>& cos_theta_cache = level_.levelCache().GPU_cos_theta();

    const GPU_Vector<double>& coeff_alpha_cache = level_.levelCache().GPU_coeff_alpha();
    const GPU_Vector<double>& coeff_beta_cache = level_.levelCache().GPU_coeff_beta();

    dim3 threadsPerBlock(1,1);
    dim3 numBlocks(1,1);

    DomainGeometry* device_domain_geometry;
    cudaMalloc(&device_domain_geometry, sizeof(DomainGeometry));
    cudaMemcpy(device_domain_geometry, &domain_geometry_, sizeof(DomainGeometry), cudaMemcpyHostToDevice);

    /* We use precomputed DensityProfileCoefficients values. */
    // DensityProfileCoefficients* device_density_profile;
    // cudaMalloc(&device_density_profile, sizeof(DensityProfileCoefficients));
    // cudaMemcpy(device_density_profile, &density_profile_coefficients_, sizeof(DensityProfileCoefficients), cudaMemcpyHostToDevice);
    
    threadsPerBlock = dim3(16, 16);
    numBlocks = dim3((grid.numberSmootherCircles() + 14 - 1) / 14, (grid.ntheta() + 14 - 1) / 14);
    const int start_black_circles = (grid.numberSmootherCircles() % 2 == 0) ? 1 : 0;
    applyAscOrtho_Circle_kernel<<<numBlocks, threadsPerBlock>>>(
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

    int start_solver_black_i_r = (grid.numberSmootherCircles() % 2 == 0) ? 1 : 2;
    cusparseDgtsv2StridedBatch(
        handle_, 
        grid.ntheta(), 
        circle_lower_diagonals_ + start_solver_black_i_r * grid.ntheta(), 
        circle_main_diagonals_ + start_solver_black_i_r * grid.ntheta(), 
        circle_upper_diagonals_ + start_solver_black_i_r * grid.ntheta(), 
        x.data() + start_solver_black_i_r * grid.ntheta(), 
        (grid.numberSmootherCircles()) / 2, 
        2 * grid.ntheta(), 
        pBuffer_
    );

    cusparseDgtsv2StridedBatch(
        handle_,
        grid.ntheta(), 
        circle_lower_diagonals_ + start_solver_black_i_r * grid.ntheta(), 
        circle_main_diagonals_ + start_solver_black_i_r * grid.ntheta(), 
        circle_upper_diagonals_ + start_solver_black_i_r * grid.ntheta(), 
        temp.data() + start_solver_black_i_r * grid.ntheta(), 
        (grid.numberSmootherCircles()) / 2, 
        2 * grid.ntheta(), 
        pBuffer_
    );

    threadsPerBlock = dim3(16, 16);
    numBlocks = dim3(((grid.numberSmootherCircles())/2 + 16 - 1) / 16, (grid.ntheta() + 16 - 1) / 16);
    combine_circle_solutions<<<numBlocks, threadsPerBlock>>>(
        x.data(), temp.data(),
        circle_lower_diagonals_, circle_main_diagonals_, circle_upper_diagonals_,
        sherman_morrison_gammas_,
        start_solver_black_i_r,
        level_.device_grid()
    );
    cudaDeviceSynchronize();


    threadsPerBlock = dim3(16, 16);
    numBlocks = dim3((grid.numberSmootherCircles() + 14 - 1) / 14, (grid.ntheta() + 14 - 1) / 14);
    const int start_white_circles = (grid.numberSmootherCircles() % 2 == 0) ? 0 : 1;
    applyAscOrtho_Circle_kernel<<<numBlocks, threadsPerBlock>>>(
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


   int start_solver_white_i_r = (grid.numberSmootherCircles() % 2 == 0) ? 2 : 1;
    cusparseDgtsv2StridedBatch(
        handle_, 
        grid.ntheta(), 
        circle_lower_diagonals_ + start_solver_white_i_r * grid.ntheta(), 
        circle_main_diagonals_ + start_solver_white_i_r * grid.ntheta(), 
        circle_upper_diagonals_ + start_solver_white_i_r * grid.ntheta(), 
        x.data() + start_solver_white_i_r * grid.ntheta(), 
        (grid.numberSmootherCircles()) / 2, 
        2 * grid.ntheta(), 
        pBuffer_
    );

    cusparseDgtsv2StridedBatch(
        handle_,
        grid.ntheta(), 
        circle_lower_diagonals_ + start_solver_white_i_r * grid.ntheta(), 
        circle_main_diagonals_ + start_solver_white_i_r * grid.ntheta(), 
        circle_upper_diagonals_ + start_solver_white_i_r * grid.ntheta(), 
        temp.data() + start_solver_white_i_r * grid.ntheta(), 
        (grid.numberSmootherCircles()) / 2, 
        2 * grid.ntheta(), 
        pBuffer_
    );

    threadsPerBlock = dim3(16, 16);
    numBlocks = dim3(((grid.numberSmootherCircles())/2 + 16 - 1) / 16, (grid.ntheta() + 16 - 1) / 16);
    combine_circle_solutions<<<numBlocks, threadsPerBlock>>>(
        x.data(), temp.data(),
        circle_lower_diagonals_, circle_main_diagonals_, circle_upper_diagonals_,
        sherman_morrison_gammas_,
        start_solver_white_i_r,
        level_.device_grid()
    );
    cudaDeviceSynchronize();




    threadsPerBlock = dim3(16, 16);
    numBlocks = dim3((grid.lengthSmootherRadial() + 14 - 1) / 14, (grid.ntheta() + 14 - 1) / 14);
    const int start_black_radials = 0;
    applyAscOrtho_Radial_kernel<<<numBlocks, threadsPerBlock>>>(
        x.data(), rhs.data(),
        level_.device_grid(), DirBC_Interior_,
        start_black_radials,
        device_domain_geometry, 
        coeff_alpha_cache.data(), coeff_beta_cache.data(), 
        sin_theta_cache.data(), cos_theta_cache.data()
    );
    cudaDeviceSynchronize();


    cusparseDgtsv2StridedBatch(
        handle_, grid.lengthSmootherRadial(), radial_lower_diagonals_, radial_main_diagonals_, radial_upper_diagonals_, x.data() + grid.numberCircularSmootherNodes(), grid.ntheta() / 2, 2 * grid.lengthSmootherRadial(), pBuffer_
    );

    threadsPerBlock = dim3(16, 16);
    numBlocks = dim3((grid.lengthSmootherRadial() + 14 - 1) / 14, (grid.ntheta() + 14 - 1) / 14);
    const int start_white_radials = 1;
    applyAscOrtho_Radial_kernel<<<numBlocks, threadsPerBlock>>>(
        x.data(), rhs.data(),
        level_.device_grid(), DirBC_Interior_,
        start_white_radials,
        device_domain_geometry, 
        coeff_alpha_cache.data(), coeff_beta_cache.data(), 
        sin_theta_cache.data(), cos_theta_cache.data()
    );
    cudaDeviceSynchronize();

    cusparseDgtsv2StridedBatch(
        handle_, grid.lengthSmootherRadial(), radial_lower_diagonals_ + grid.lengthSmootherRadial(), radial_main_diagonals_ + grid.lengthSmootherRadial(), radial_upper_diagonals_ + grid.lengthSmootherRadial(), x.data() + grid.numberCircularSmootherNodes() + grid.lengthSmootherRadial(), grid.ntheta() / 2, 2 * grid.lengthSmootherRadial(), pBuffer_
    );

    /* We use precomputed DensityProfileCoefficients values. */
    cudaFree(device_domain_geometry);
    // cudaFree(device_density_profile);

    std::cout<<x<<std::endl;
}







