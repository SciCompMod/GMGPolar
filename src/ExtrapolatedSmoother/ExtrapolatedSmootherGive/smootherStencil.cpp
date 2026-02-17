#include "../../../include/ExtrapolatedSmoother/ExtrapolatedSmootherGive/extrapolatedSmootherGive.h"

const Stencil& ExtrapolatedSmootherGive::getStencil(int i_r, int i_theta) const
{
    assert(0 <= i_r && i_r < grid_.nr());

    assert((grid_.ntheta() / 2) % 2 == 0);

    if (i_r == 0) {
        if (i_theta % 2 == 0) {
            return stencil_center_;
        }
        else {
            if (!DirBC_Interior_) {
                return stencil_center_left_;
            }
            else {
                return stencil_center_;
            }
        }
    }

    throw std::out_of_range("getStencil: Only i_r = 0 implemented.");
}

int ExtrapolatedSmootherGive::getNonZeroCountCircleAsc(const int i_r) const
{
    assert(i_r >= 0 && i_r < grid_.numberSmootherCircles());

    assert((grid_.ntheta() / 2) % 2 == 0);

    if (i_r == 0) {
        if (!DirBC_Interior_) {
            return grid_.ntheta() / 2 + 2 * (grid_.ntheta() / 2);
        }
        else {
            return grid_.ntheta();
        }
    }

    throw std::out_of_range("nnz_circle_Asc: Only i_r = 0 implemented.");
}

int ExtrapolatedSmootherGive::getCircleAscIndex(const int i_r, const int i_theta) const
{
    assert(i_r >= 0 && i_r < grid_.numberSmootherCircles());

    assert((grid_.ntheta() / 2) % 2 == 0);

    if (i_r == 0) {
        if (!DirBC_Interior_) {
            if (i_theta % 2 == 0) {
                return 3 * (i_theta / 2);
            }
            else {
                return 3 * (i_theta / 2) + 1;
            }
        }
        else {
            return i_theta;
        }
    }

    throw std::out_of_range("ptr_nz_index_circle_Asc: Only i_r = 0 implemented.");
}
