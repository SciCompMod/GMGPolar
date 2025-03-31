#include "../../../include/DirectSolver/DirectSolverTakeCustomLU/directSolverTakeCustomLU.h"

int DirectSolverTakeCustomLU::getStencilSize(int global_index) const
{
    int i_r, i_theta;
    grid_.multiIndex(global_index, i_r, i_theta);

    const int size_stencil_inner_boundary      = DirBC_Interior_ ? 1 : 7;
    const int size_stencil_next_inner_boundary = DirBC_Interior_ ? 9 : 9;
    const int size_stencil_interior            = 9;
    const int size_stencil_next_outer_boundary = 9;
    const int size_stencil_outer_boundary      = 1;

    if ((i_r > 1 && i_r < grid_.nr() - 2) || (i_r == 1 && !DirBC_Interior_)) {
        return size_stencil_interior;
    }
    else if (i_r == 0 && !DirBC_Interior_) {
        return size_stencil_inner_boundary;
    }
    else if ((i_r == 0 && DirBC_Interior_) || i_r == grid_.nr() - 1) {
        return size_stencil_outer_boundary;
    }
    else if (i_r == 1 && DirBC_Interior_) {
        return size_stencil_next_inner_boundary;
    }
    else if (i_r == grid_.nr() - 2) {
        return size_stencil_next_outer_boundary;
    }
    throw std::out_of_range("Invalid index for stencil");
}

const Stencil& DirectSolverTakeCustomLU::getStencil(int i_r) const
{
    assert(0 <= i_r && i_r < grid_.nr());
    assert(grid_.nr() >= 4);

    if ((i_r > 1 && i_r < grid_.nr() - 2) || (i_r == 1 && !DirBC_Interior_)) {
        return stencil_interior_;
    }
    else if (i_r == 0 && !DirBC_Interior_) {
        return stencil_across_origin_;
    }
    else if ((i_r == 0 && DirBC_Interior_) || i_r == grid_.nr() - 1) {
        return stencil_DB_;
    }
    else if (i_r == 1 && DirBC_Interior_) {
        return stencil_next_inner_DB_;
    }
    else if (i_r == grid_.nr() - 2) {
        return stencil_next_outer_DB_;
    }
    throw std::out_of_range("Invalid index for stencil");
}

int DirectSolverTakeCustomLU::getNonZeroCountSolverMatrix() const
{
    const int size_stencil_inner_boundary      = DirBC_Interior_ ? 1 : 7;
    const int size_stencil_next_inner_boundary = DirBC_Interior_ ? 9 : 9;
    const int size_stencil_interior            = 9;
    const int size_stencil_next_outer_boundary = 9;
    const int size_stencil_outer_boundary      = 1;

    assert(grid_.nr() >= 4);

    return grid_.ntheta() *
           (size_stencil_inner_boundary + size_stencil_next_inner_boundary + (grid_.nr() - 4) * size_stencil_interior +
            size_stencil_next_outer_boundary + size_stencil_outer_boundary);
}
