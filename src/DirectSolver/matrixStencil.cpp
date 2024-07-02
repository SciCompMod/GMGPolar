#include "../../include/DirectSolver/directSolver.h"

const Stencil& DirectSolver::get_stencil(int i_r) const {
    assert(0 <= i_r && i_r < grid_.nr());

    static const Stencil stencil_interior = 
        {7, 4, 8,
        1, 0, 2,
        5, 3, 6};

    static const Stencil stencil_across_origin = 
        {-1, 4, 6,
        1, 0, 2,
        -1, 3, 5};

    static const Stencil stencil_DB = 
        {-1, -1, -1,
        -1,  0, -1,
        -1, -1, -1};

    static const Stencil stencil_next_inner_DB = 
        {-1, 3, 5,
        -1, 0, 1,
        -1, 2, 4};

    static const Stencil stencil_next_outer_DB = 
        {5, 3, -1,
        1, 0, -1,
        4, 2, -1};

    assert(grid_.nr() >= 4);

    if ((i_r > 1 && i_r < grid_.nr() - 2) || (i_r == 1 && !DirBC_Interior_)) {
        return stencil_interior;
    } 
    else if(i_r == 0 && !DirBC_Interior_) {
        return stencil_across_origin;
    }
    else if((i_r == 0 && DirBC_Interior_) || i_r == grid_.nr() - 1) {
        return stencil_DB;
    }
    else if(i_r == 1 && DirBC_Interior_) {
        return stencil_next_inner_DB;
    }
    else if(i_r == grid_.nr() - 2) {
        return stencil_next_outer_DB;
    }
    throw std::out_of_range("Invalid index for stencil");
}

int DirectSolver::nnz_matrixA() const {
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

int DirectSolver::ptr_nz_index_matrixA(const int i_r, const int i_theta) const {
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
