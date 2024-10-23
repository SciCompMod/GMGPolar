#include "../../../include/DirectSolver/DirectSolverTake/directSolverTake.h"

const Stencil& DirectSolverTake::getStencil(int i_r) const
{
    assert(0 <= i_r && i_r < grid_.nr());
    assert(grid_.nr() >= 4);

    if ((i_r > 1 && i_r < grid_.nr() - 2) || (i_r == 1 && !DirBC_Interior_))
    {
        return stencil_interior_;
    }
    else if (i_r == 0 && !DirBC_Interior_)
    {
        return stencil_across_origin_;
    }
    else if ((i_r == 0 && DirBC_Interior_) || i_r == grid_.nr() - 1)
    {
        return stencil_DB_;
    }
    else if (i_r == 1 && DirBC_Interior_)
    {
        return stencil_next_inner_DB_;
    }
    else if (i_r == grid_.nr() - 2)
    {
        return stencil_next_outer_DB_;
    }
    throw std::out_of_range("Invalid index for stencil");
}

// clang-format off
int DirectSolverTake::getNonZeroCountSolverMatrix() const {
    const int size_stencil_inner_boundary = DirBC_Interior_ ? 1 : 7 ;
    const int size_stencil_next_inner_boundary = DirBC_Interior_ ? 6 : 9;
    const int size_stencil_interior = 9;
    const int size_stencil_next_outer_boundary = 6;
    const int size_stencil_outer_boundary = 1;

    assert(grid_.nr() >= 4);

    return grid_.ntheta() * ( 
        size_stencil_inner_boundary +
        size_stencil_next_inner_boundary +
        (grid_.nr()-4) * size_stencil_interior +
        size_stencil_next_outer_boundary +
        size_stencil_outer_boundary
    ); 
}

/* ----------------------------------------------------------------- */
/* If the indexing is not smoother-based, please adjust the indexing */
int DirectSolverTake::getSolverMatrixIndex(const int i_r, const int i_theta) const {
    const int size_stencil_inner_boundary = DirBC_Interior_ ? 1 : 7 ;
    const int size_stencil_next_inner_boundary = DirBC_Interior_ ? 6 : 9;
    const int size_stencil_interior = 9;
    const int size_stencil_next_outer_boundary = 6;
    const int size_stencil_outer_boundary = 1;

    assert(grid_.nr() >= 4);
    assert(grid_.numberSmootherCircles() >= 2);
    assert(grid_.lengthSmootherRadial() >= 2);

    if(1 < i_r && i_r < grid_.numberSmootherCircles()) {
        // Interior: Circle index section
        const int prior_inner_boundary_nodes = grid_.ntheta();
        const int prior_next_inner_boundary_nodes = grid_.ntheta();
        const int prior_interior_nodes = (i_r-2) * grid_.ntheta() + i_theta;
        return
            size_stencil_inner_boundary * prior_inner_boundary_nodes +
            size_stencil_next_inner_boundary * prior_next_inner_boundary_nodes +
            size_stencil_interior * prior_interior_nodes;
    }
    if(i_r >= grid_.numberSmootherCircles() && i_r < grid_.nr() - 2) {
        // Interior: Radial index section
        const int prior_inner_boundary_nodes = grid_.ntheta();
        const int prior_next_inner_boundary_nodes = grid_.ntheta();
        const int prior_interior_nodes = 
            grid_.ntheta() * (grid_.numberSmootherCircles()-2) +
            i_theta * (grid_.lengthSmootherRadial()-2) + 
            i_r - grid_.numberSmootherCircles();
        const int prior_next_outer_boundary_nodes = i_theta;
        const int prior_outer_boundary_nodes = i_theta;
        return
            size_stencil_inner_boundary * prior_inner_boundary_nodes +
            size_stencil_next_inner_boundary * prior_next_inner_boundary_nodes +
            size_stencil_interior * prior_interior_nodes + 
            size_stencil_next_outer_boundary * prior_next_outer_boundary_nodes +
            size_stencil_outer_boundary * prior_outer_boundary_nodes;
    }
    else if(i_r == 0) {
        // Inner boundary
        const int prior_inner_boundary_nodes = i_theta;
        return
            size_stencil_inner_boundary * prior_inner_boundary_nodes;
    }
    else if(i_r == 1) {
        // Next to inner boundary
        const int prior_inner_boundary_nodes = grid_.ntheta();
        const int prior_next_inner_boundary_nodes = i_theta;
        return
            size_stencil_inner_boundary * prior_inner_boundary_nodes +
            size_stencil_next_inner_boundary  * prior_next_inner_boundary_nodes;
    }
    else if(i_r == grid_.nr()-2) {
        // Next to outer boundary
        const int prior_inner_boundary_nodes = grid_.ntheta();
        const int prior_next_inner_boundary_nodes = grid_.ntheta();
        const int prior_interior_nodes = 
            grid_.ntheta() * (grid_.numberSmootherCircles()-2) +
            (i_theta + 1) * (grid_.lengthSmootherRadial()-2);
        const int prior_next_outer_boundary_nodes = i_theta;
        const int prior_outer_boundary_nodes = i_theta;
        return
            size_stencil_inner_boundary * prior_inner_boundary_nodes +
            size_stencil_next_inner_boundary * prior_next_inner_boundary_nodes +
            size_stencil_interior * prior_interior_nodes + 
            size_stencil_next_outer_boundary * prior_next_outer_boundary_nodes +
            size_stencil_outer_boundary * prior_outer_boundary_nodes;
    }
    else if(i_r == grid_.nr()-1) {
        // Outer boundary
        const int prior_inner_boundary_nodes = grid_.ntheta();
        const int prior_next_inner_boundary_nodes = grid_.ntheta();
        const int prior_interior_nodes = 
            grid_.ntheta() * (grid_.numberSmootherCircles()-2) +
            (i_theta + 1) * (grid_.lengthSmootherRadial()-2);
        const int prior_next_outer_boundary_nodes = i_theta + 1;
        const int prior_outer_boundary_nodes = i_theta;
        return
            size_stencil_inner_boundary * prior_inner_boundary_nodes +
            size_stencil_next_inner_boundary * prior_next_inner_boundary_nodes +
            size_stencil_interior * prior_interior_nodes + 
            size_stencil_next_outer_boundary * prior_next_outer_boundary_nodes +
            size_stencil_outer_boundary * prior_outer_boundary_nodes;
    }
    throw std::out_of_range("Invalid index for stencil");
}
