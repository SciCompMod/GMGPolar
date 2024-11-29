#include "../../../include/Smoother/SmootherTake/smootherTake.h"

const Stencil& SmootherTake::getStencil(int i_r) const
{
    assert(0 <= i_r && i_r < grid_.nr());

    assert(grid_.numberSmootherCircles() >= 2);
    assert(grid_.lengthSmootherRadial() >= 3);

    const int numberSmootherCircles = grid_.numberSmootherCircles();

    if (i_r < numberSmootherCircles) {
        /* Circle Section */
        if (i_r > 0 && i_r < numberSmootherCircles) {
            return circle_stencil_interior_;
        }
        else if (i_r == 0) {
            return DirBC_Interior_ ? stencil_DB_ : circle_stencil_across_origin_;
        }
    }
    else {
        /* Radial Section */
        if (i_r > numberSmootherCircles && i_r < grid_.nr() - 2) {
            return radial_stencil_interior_;
        }
        else if (i_r == numberSmootherCircles) {
            return radial_stencil_next_circular_smoothing_;
        }
        else if (i_r == grid_.nr() - 1) {
            return stencil_DB_;
        }
        else if (i_r == grid_.nr() - 2) {
            return radial_stencil_next_outer_DB_;
        }
    }
    throw std::out_of_range("Invalid index for stencil");
}

int SmootherTake::getNonZeroCountCircleAsc(const int i_r) const
{
    assert(i_r >= 0 && i_r < grid_.numberSmootherCircles());

    const int numberSmootherCircles = grid_.numberSmootherCircles();

    const int size_stencil_inner_boundary = DirBC_Interior_ ? 1 : 4;
    const int size_stencil_interior = 3;

    if (i_r > 0) {
        return size_stencil_interior * grid_.ntheta();
    }
    else if (i_r == 0) {
        return size_stencil_inner_boundary * grid_.ntheta();
    }
    throw std::out_of_range("Invalid index for nnz_circle_Asc");
}

int SmootherTake::getCircleAscIndex(const int i_r, const int i_theta) const
{
    assert(i_r >= 0 && i_r < grid_.numberSmootherCircles());

    const int numberSmootherCircles = grid_.numberSmootherCircles();

    const int size_stencil_inner_boundary = DirBC_Interior_ ? 1 : 4;
    const int size_stencil_interior = 3;

    if (i_r > 0) {
        return size_stencil_interior * i_theta;
    }
    else {
        return size_stencil_inner_boundary * i_theta;
    }
}

int SmootherTake::getNonZeroCountRadialAsc(const int i_theta) const
{
    assert(i_theta >= 0 && i_theta < grid_.ntheta());

    const int size_stencil_next_circluar_smoothing = 2;
    const int size_stencil_interior = 3;
    const int size_stencil_next_outer_boundary = 2;
    const int size_stencil_outer_boundary = 1;

    assert(grid_.lengthSmootherRadial() >= 3);

    return size_stencil_next_circluar_smoothing + (grid_.lengthSmootherRadial() - 3) * size_stencil_interior + size_stencil_next_outer_boundary +
           size_stencil_outer_boundary;
}

int SmootherTake::getRadialAscIndex(const int i_r, const int i_theta) const
{
    assert(i_theta >= 0 && i_theta < grid_.ntheta());

    const int size_stencil_next_circluar_smoothing = 2;
    const int size_stencil_interior = 3;
    const int size_stencil_next_outer_boundary = 2;
    const int size_stencil_outer_boundary = 1;

    assert(grid_.lengthSmootherRadial() >= 3);
    assert(grid_.numberSmootherCircles() >= 2);

    const int numberSmootherCircles = grid_.numberSmootherCircles();

    if (i_r > numberSmootherCircles && i_r < grid_.nr() - 2) {
        return size_stencil_next_circluar_smoothing + (i_r - numberSmootherCircles - 1) * size_stencil_interior;
    }
    else if (i_r == numberSmootherCircles) {
        return 0;
    }
    else if (i_r == grid_.nr() - 2) {
        return size_stencil_next_circluar_smoothing + (grid_.lengthSmootherRadial() - 3) * size_stencil_interior;
    }
    else if (i_r == grid_.nr() - 1) {
        return size_stencil_next_circluar_smoothing + (grid_.lengthSmootherRadial() - 3) * size_stencil_interior + size_stencil_next_outer_boundary;
    }
    throw std::out_of_range("Invalid index for stencil");
}
