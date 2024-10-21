#include "../../../include/ExtrapolatedSmoother/ExtrapolatedSmootherTake/extrapolatedSmootherTake.h"

const Stencil& ExtrapolatedSmootherTake::getStencil(int i_r, int i_theta) const {
    assert(0 <= i_r && i_r < grid_.nr());

    assert((grid_.ntheta() / 2) % 2 == 0);

    if(i_r == 0){
        static const Stencil stencil_Center = 
            {-1, -1, -1,
            -1,  0, -1,
            -1, -1, -1};

        static const Stencil stencil_Center_Left = 
            {-1, -1, -1,
            1,  0, -1,
            -1, -1, -1};

        if (i_theta % 2 == 0) {
            return stencil_Center;
        } else {
            if (!DirBC_Interior_) {
                return stencil_Center_Left;
            } else {
                return stencil_Center;
            }
        }
    }

    throw std::out_of_range("getStencil: Only i_r = 0 implemented.");
}

int ExtrapolatedSmootherTake::getNonZeroCountCircleAsc(const int i_r) const {
    assert(i_r >= 0 && i_r < grid_.numberSmootherCircles());

    assert((grid_.ntheta() / 2) % 2 == 0);

    if(i_r == 0){
        if (!DirBC_Interior_) {
            return grid_.ntheta()/2 + 2*(grid_.ntheta()/2);
        } else {
            return grid_.ntheta();
        }

    }

    throw std::out_of_range("nnz_circle_Asc: Only i_r = 0 implemented.");
}


int ExtrapolatedSmootherTake::getCircleAscIndex(const int i_r, const int i_theta) const {
    assert(i_r >= 0 && i_r < grid_.numberSmootherCircles());

    assert((grid_.ntheta() / 2) % 2 == 0);

    if(i_r == 0){
        if (!DirBC_Interior_) {
            if(i_theta % 2 == 0){
                return 3 * (i_theta/2);
            }
            else{
                return 3 * (i_theta/2) + 1;
            }
        }
        else{
            return i_theta;
        }
    }

    throw std::out_of_range("ptr_nz_index_circle_Asc: Only i_r = 0 implemented.");
}



int ExtrapolatedSmootherTake::getNonZeroCountRadialAsc(const int i_theta) const {
    throw std::out_of_range("ExtrapolatedSmoother: nnz_radial_Asc not implemented.");
}



int ExtrapolatedSmootherTake::getRadialAscIndex(const int i_r, const int i_theta) const {
    throw std::out_of_range("ExtrapolatedSmoother: ptr_nz_index_radial_Asc not implemented.");
}
