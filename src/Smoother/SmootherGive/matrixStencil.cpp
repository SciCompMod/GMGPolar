#include "../../../include/Smoother/SmootherGive/smootherGive.h"

#include "../../../include/Smoother/SmootherTake/smootherTake.h"

const Stencil& SmootherGive::getStencil(int i_r) const
{
    if (i_r == 0) {
        return DirBC_Interior_ ? stencil_DB_ : circle_stencil_across_origin_;
    }
    else {
        throw std::out_of_range("getStencil: Only i_r = 0 implemented.");
    }
}

int SmootherGive::getNonZeroCountCircleAsc(int i_r) const
{
    if (i_r == 0) {
        const int size_stencil_inner_boundary = DirBC_Interior_ ? 1 : 4;
        return size_stencil_inner_boundary * grid_.ntheta();
    }
    else {
        throw std::out_of_range("getNonZeroCountCircleAsc: Only i_r = 0 implemented.");
    }
}

int SmootherGive::getCircleAscIndex(int i_r, int i_theta) const
{
    if (i_r == 0) {
        const int size_stencil_inner_boundary = DirBC_Interior_ ? 1 : 4;
        return size_stencil_inner_boundary * i_theta;
    }
    else {
        throw std::out_of_range("getCircleAscIndex: Only i_r = 0 implemented.");
    }
}
