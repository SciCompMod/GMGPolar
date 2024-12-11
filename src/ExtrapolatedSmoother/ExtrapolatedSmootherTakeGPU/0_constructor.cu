#include "../../../include/ExtrapolatedSmoother/ExtrapolatedSmootherTakeGPU/extrapolatedSmoother.h"

/* Constructor */

ExtrapolatedSmootherTakeGPU::ExtrapolatedSmootherTakeGPU(const Level& level, const DomainGeometry& domain_geometry,
                   const DensityProfileCoefficients& density_profile_coefficients, bool DirBC_Interior)
    /* Constructor Members */
    : level_(level)
    , domain_geometry_(domain_geometry)
    , density_profile_coefficients_(density_profile_coefficients)
    , DirBC_Interior_(DirBC_Interior)
    /* Circle Tridiagonal Matrices */
    , circle_main_diagonals_(nullptr)
    , circle_lower_diagonals_(nullptr)
    , circle_upper_diagonals_(nullptr)
    , sherman_morrison_gammas_(nullptr)
    /* Radial Tridiagonal Matrices */
    , radial_main_diagonals_(nullptr)
    , radial_lower_diagonals_(nullptr)
    , radial_upper_diagonals_(nullptr)
    /* Tridiagonal Solver Buffer */
    , pBuffer_(nullptr)
    /* Inner Boundary Mumps COO Matrix */
    , inner_boundary_matrix_row_indices_(nullptr)
    , inner_boundary_matrix_column_indices_(nullptr)
    , inner_boundary_matrix_values_(nullptr)
    , d_inner_boundary_matrix_row_indices_(nullptr)
    , d_inner_boundary_matrix_column_indices_(nullptr)
    , d_inner_boundary_matrix_values_(nullptr)
{
    const PolarGrid& grid = level.grid();

    int nr = grid.nr();
    int ntheta = grid.ntheta();
    int number_smoother_circles = grid.numberSmootherCircles();
    int length_smoother_radial = grid.lengthSmootherRadial();

    int circle_batch_count = number_smoother_circles;
    int circle_m = ntheta;
    /* Cyclic Tridiagonal Circle Matrices */
    cudaMalloc(&circle_lower_diagonals_, circle_m * circle_batch_count * sizeof(double));
    cudaMemset(circle_lower_diagonals_, 0, circle_m * circle_batch_count * sizeof(double));
    cudaMalloc(&circle_main_diagonals_, circle_m * circle_batch_count * sizeof(double));
    cudaMemset(circle_main_diagonals_, 0, circle_m * circle_batch_count * sizeof(double));
    cudaMalloc(&circle_upper_diagonals_, circle_m * circle_batch_count * sizeof(double));
    cudaMemset(circle_upper_diagonals_, 0, circle_m * circle_batch_count * sizeof(double));
    /* Cuda does not supply a cyclic tridiagonal solver. */
    /* Thus we use the Sherman–Morrison formula to reduce the problem to a simple tridiagonal problem with two right hand sides. */
    cudaMalloc(&sherman_morrison_gammas_, circle_batch_count * sizeof(double));
    cudaMemset(sherman_morrison_gammas_, 0, circle_batch_count * sizeof(double));
    cudaMalloc(&factor_, circle_batch_count * sizeof(double));
    cudaMemset(factor_, 0, circle_batch_count * sizeof(double));
    /* Remark: The 1st cylic tridiagonal matrix on the interior boundary is unused. */

    int radial_batch_count = ntheta;
    int radial_m = length_smoother_radial;
    /* Tridiagonal Radial Matrices */
    cudaMalloc(&radial_lower_diagonals_, radial_m * radial_batch_count * sizeof(double));
    cudaMemset(radial_lower_diagonals_, 0, radial_m * radial_batch_count * sizeof(double));
    cudaMalloc(&radial_main_diagonals_, radial_m * radial_batch_count * sizeof(double));
    cudaMemset(radial_main_diagonals_, 0, radial_m * radial_batch_count * sizeof(double));
    cudaMalloc(&radial_upper_diagonals_, radial_m * radial_batch_count * sizeof(double));
    cudaMemset(radial_upper_diagonals_, 0, radial_m * radial_batch_count * sizeof(double));

    /* Tridiaginal Cuda Solver */
    cusparseCreate(&sparse_handle_);

    /* Allocate Tridiagonal Solver Buffer. */
    size_t bufferSizeInBytes_Circle;
    cusparseDgtsv2StridedBatch_bufferSizeExt(
        sparse_handle_, ntheta, nullptr, nullptr, nullptr, nullptr, 
        number_smoother_circles / 2, 2 * ntheta, &bufferSizeInBytes_Circle);

    size_t bufferSizeInBytes_Radial;
    cusparseDgtsv2StridedBatch_bufferSizeExt(
        sparse_handle_, length_smoother_radial, nullptr, nullptr, nullptr, nullptr, 
        ntheta / 2, 2 * length_smoother_radial, &bufferSizeInBytes_Radial);

    /* The Tridiagonal solvers require 5 * batch_count * m * sizeof(double) bytes. */
    /* Alternatively use four different pBuffers for each tridiagonal batch solver. */
    size_t max_pBufferSizeInBytes = std::max(bufferSizeInBytes_Circle, bufferSizeInBytes_Radial);
    cudaMalloc(&pBuffer_, max_pBufferSizeInBytes);

    /* Inner Boundary Mumps COO Matrix */
    inner_boundary_matrix_nnz_ = DirBC_Interior_ ? grid.ntheta() : grid.ntheta() / 2 + 2 * (grid.ntheta() / 2); 

    inner_boundary_matrix_row_indices_ = std::make_unique<int[]>(inner_boundary_matrix_nnz_);
    inner_boundary_matrix_column_indices_ = std::make_unique<int[]>(inner_boundary_matrix_nnz_);
    inner_boundary_matrix_values_ = std::make_unique<double[]>(inner_boundary_matrix_nnz_);

    cudaMalloc(&d_inner_boundary_matrix_row_indices_, inner_boundary_matrix_nnz_ * sizeof(int));
    cudaMemset(d_inner_boundary_matrix_row_indices_, 0, inner_boundary_matrix_nnz_* sizeof(int));
    cudaMalloc(&d_inner_boundary_matrix_column_indices_, inner_boundary_matrix_nnz_ * sizeof(int));
    cudaMemset(d_inner_boundary_matrix_column_indices_, 0, inner_boundary_matrix_nnz_* sizeof(int));
    cudaMalloc(&d_inner_boundary_matrix_values_, inner_boundary_matrix_nnz_ * sizeof(double));
    cudaMemset(d_inner_boundary_matrix_values_, 0, inner_boundary_matrix_nnz_* sizeof(double));

    /* Build Smoother Matrices we have allocated. */
    buildAscMatrices();

    initializeMumps();

    /* The cyclic tridiagonal Matrices need to be adjusted to a system of a non-cyclic tridiagonal matrices. */
    adjustAscCircle_ShermanMorrison(); 
}



/* Destructor */

ExtrapolatedSmootherTakeGPU::~ExtrapolatedSmootherTakeGPU() {
    /* Cyclic Tridiagonal Circle Matrices */
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
    /* Cuda does not supply a cyclic tridiagonal solver. */
    /* Thus we use the Sherman–Morrison formula to reduce the problem to a simple tridiagonal problem with two right hand sides. */
    if (sherman_morrison_gammas_) {
        cudaFree(sherman_morrison_gammas_);
        sherman_morrison_gammas_ = nullptr;
    }
    if (factor_) {
        cudaFree(factor_);
        factor_ = nullptr;
    }

    /* Tridiagonal Radial Matrices */
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

    /* Tridiaginal Cuda Solver */
    cusparseDestroy(sparse_handle_);

    /* Free Tridiagonal Solver Buffer. */
    if (pBuffer_) {
        cudaFree(pBuffer_);
        pBuffer_ = nullptr;
    }

    /* Inner Boundary Mumps COO Matrix */
    if (d_inner_boundary_matrix_row_indices_) {
        cudaFree(d_inner_boundary_matrix_row_indices_);
        d_inner_boundary_matrix_row_indices_ = nullptr;
    }
    if (d_inner_boundary_matrix_column_indices_) {
        cudaFree(d_inner_boundary_matrix_column_indices_);
        d_inner_boundary_matrix_column_indices_ = nullptr;
    }
    if (d_inner_boundary_matrix_values_) {
        cudaFree(d_inner_boundary_matrix_values_);
        d_inner_boundary_matrix_values_ = nullptr;
    }

    finalizeMumpsSolver();
}