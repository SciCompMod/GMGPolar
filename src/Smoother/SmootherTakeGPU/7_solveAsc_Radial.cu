#include "../../../include/Smoother/SmootherTakeGPU/smoother.h"


void SmootherTakeGPU::solveAsc_BlackRadial(GPU_Vector<double>& x, const GPU_Vector<double>& rhs)
{
    const PolarGrid& grid = level_.grid();

    const int start_black_radial_solver = 0;
    int start = start_black_radial_solver * grid.lengthSmootherRadial();
    int batch_count = grid.ntheta() / 2;
    int batch_stride = 2 * grid.lengthSmootherRadial();
    cusparseDgtsv2StridedBatch(
        sparse_handle_, 
        grid.lengthSmootherRadial(), 
        radial_lower_diagonals_ + start,
        radial_main_diagonals_ + start,
        radial_upper_diagonals_ + start, 
        x.data() + grid.numberCircularSmootherNodes() + start, 
        batch_count, batch_stride, pBuffer_
    );
}


void SmootherTakeGPU::solveAsc_WhiteRadial(GPU_Vector<double>& x, const GPU_Vector<double>& rhs)
{
    const PolarGrid& grid = level_.grid();

    const int start_white_radial_solver = 1;
    int start = start_white_radial_solver * grid.lengthSmootherRadial();
    int batch_count = grid.ntheta() / 2;
    int batch_stride = 2 * grid.lengthSmootherRadial();
    cusparseDgtsv2StridedBatch(
        sparse_handle_, 
        grid.lengthSmootherRadial(), 
        radial_lower_diagonals_ + start,
        radial_main_diagonals_ + start,
        radial_upper_diagonals_ + start, 
        x.data() + grid.numberCircularSmootherNodes() + start, 
        batch_count, batch_stride, pBuffer_
    );
}