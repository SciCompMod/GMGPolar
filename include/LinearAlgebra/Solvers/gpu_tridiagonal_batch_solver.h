// #include <cuda_runtime.h>
// #include <cusparse_v2.h>
// #include <stdexcept>

// class GPU_SymmetricTridiagonalBatchSolver {
// public:
//     GPU_SymmetricTridiagonalBatchSolver(int m, int batchCount, bool isCyclic)
//         : m(m)
//         , batchCount(batchCount)
//         , batchStride(m)
//         , isCyclic(isCyclic)
//         , d_dl(nullptr)
//         , d_d(nullptr)
//         , d_du(nullptr)
//         , pBuffer(nullptr) 
//     {
//         if (m < 3) {
//             throw std::invalid_argument("Matrix size (m) must be at least 3.");
//         }

//         // Initialize cuSPARSE handle
//         cusparseStatus_t status = cusparseCreate(&handle);
//         if (status != CUSPARSE_STATUS_SUCCESS) {
//             throw std::runtime_error("Failed to create cuSPARSE handle.");
//         }

//         allocateMemory();
//         determineBufferSize();
//     }

//     ~GPU_SymmetricTridiagonalBatchSolver() {
//         freeMemory();
//         cusparseDestroy(handle);
//     }

//     // Allocate GPU memory for matrix diagonals
//     void allocateMemory() {
//         cudaMalloc(&d_dl, batchCount * batchStride * sizeof(double));
//         cudaMalloc(&d_d, batchCount * batchStride * sizeof(double));
//         cudaMalloc(&d_du, batchCount * batchStride * sizeof(double));
//     }

//     // Deallocate GPU memory
//     void freeMemory() {
//         cudaFree(d_dl);
//         cudaFree(d_d);
//         cudaFree(d_du);
//         cudaFree(pBuffer);
//     }

//     // Determine buffer size required for cuSPARSE solver
//     void determineBufferSize() {
//         // bufferSize = 5 * batchCount * m * sizeof(double)
//         cusparseStatus_t status = cusparseDgtsv2StridedBatch_bufferSizeExt(
//             handle, m, d_dl, d_d, d_du, nullptr, batchCount, batchStride, &bufferSize
//         );

//         if (status != CUSPARSE_STATUS_SUCCESS) {
//             throw std::runtime_error("Failed to determine buffer size.");
//         }

//         cudaMalloc(&pBuffer, bufferSize);
//     }

//     // Solve the batch of tridiagonal systems
//     void solve(double* d_x) {
//         cusparseStatus_t status = cusparseDgtsv2StridedBatch(
//             handle, m, d_dl, d_d, d_du, d_x, batchCount, batchStride, pBuffer
//         );

//         if (status != CUSPARSE_STATUS_SUCCESS) {
//             throw std::runtime_error("Failed to solve the system.");
//         }
//     }

//     // Device-side accessors for sub-diagonal
//     __device__ __forceinline__ double& sub_diagonal(int batch_index, int index) {
//         assert(0 <= index && index < m - 1);
//         return d_dl[batch_index * batchStride + index];
//     }

//     // Device-side accessors for main diagonal
//     __device__ __forceinline__ double& main_diagonal(int batch_index, int index) {
//         return d_d[batch_index * batchStride + index];
//     }

//     // Device-side accessors for upper diagonal
//     __device__ __forceinline__ double& upper_diagonal(int batch_index, int index) {
//         assert(1 <= index && index < m);
//         return d_du[batch_index * batchStride + index];
//     }


//     // Device-side accessors for sub-diagonal
//     __device__ __forceinline__ double& sub_cyclic_corner_element(int batch_index, int index) {
//         assert(0 <= index && index < m - 1);
//         return d_dl[batch_index * batchStride + index];
//     }

// private:
//     int m;                  // Size of each tridiagonal system
//     int batchCount;         // Number of systems
//     int batchStride;        // Stride between systems
//     double* d_dl;           // Lower diagonal (on GPU)
//     double* d_d;            // Main diagonal (on GPU)
//     double* d_du;           // Upper diagonal (on GPU)
//     void* pBuffer;          // Temporary buffer (on GPU)
//     cusparseHandle_t handle; // cuSPARSE handle
//     size_t bufferSize;      // Buffer size for cuSPARSE
// };
