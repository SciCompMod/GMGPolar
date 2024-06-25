#include "../../include/Smoother/smoother.h"

const Stencil& Smoother::get_stencil(int i_r) const {

    assert(0 <= i_r && i_r < grid_.nr());

    assert(grid_.numberSmootherCircles() >= 2);
    assert(grid_.lengthSmootherRadial() >= 3);

    static const Stencil stencil_DB = 
        {-1, -1, -1,
        -1,  0, -1,
        -1, -1, -1};

    /* Circle Stencils */

    static const Stencil circle_stencil_interior = 
        {-1, 2, -1,
        -1,  0, -1,
        -1, 1, -1};

    static const Stencil circle_stencil_across_origin = 
        {-1, 3, -1,
        1,  0, -1,
        -1, 2, -1};

    /* Radial Stencils */

    static const Stencil radial_stencil_interior = 
        {-1, -1, -1,
        1, 0, 2,
        -1, -1, -1};

    static const Stencil radial_stencil_next_outer_DB = 
        {-1, -1, -1,
        1, 0, -1,
        -1, -1, -1};

    static const Stencil radial_stencil_next_circular_smoothing = 
        {-1, -1, -1,
        -1, 0, 1,
        -1, -1, -1};

    const int numberSmootherCircles = grid_.numberSmootherCircles();

    if(i_r < numberSmootherCircles){
        /* Circle Section */
        if(i_r > 0 && i_r < numberSmootherCircles){
            return circle_stencil_interior;
        } 
        else if(i_r == 0){
            return DirBC_Interior_? stencil_DB : circle_stencil_across_origin;
        }

    } else{
        /* Radial Section */
        if(i_r > numberSmootherCircles && i_r < grid_.nr()-2){
            return radial_stencil_interior;
        }
        else if(i_r == numberSmootherCircles){
            return radial_stencil_next_circular_smoothing;
        }
        else if(i_r == grid_.nr()-1){
            return stencil_DB;
        }
        else if(i_r == grid_.nr()-2){
            return radial_stencil_next_outer_DB;
        }
    }
    throw std::out_of_range("Invalid index for stencil");
}

int Smoother::nnz_circle_Asc(const int i_r) const {
    assert(i_r >= 0 && i_r < grid_.numberSmootherCircles());

    const int numberSmootherCircles = grid_.numberSmootherCircles();

    const int size_stencil_inner_boundary = DirBC_Interior_ ? 1 : 4;
    const int size_stencil_interior = 3;

    if(i_r > 0){
        return size_stencil_interior * grid_.ntheta();
    }
    else if(i_r == 0){
        return size_stencil_inner_boundary * grid_.ntheta();
    }
    throw std::out_of_range("Invalid index for nnz_circle_Asc");
}


int Smoother::ptr_nz_index_circle_Asc(const int i_r, const int i_theta) const {
    assert(i_r >= 0 && i_r < grid_.numberSmootherCircles());

    const int numberSmootherCircles = grid_.numberSmootherCircles();

    const int size_stencil_inner_boundary = DirBC_Interior_ ? 1 : 4;
    const int size_stencil_interior = 3;

    if(i_r > 0){
        return size_stencil_interior * i_theta;
    } else{
        return size_stencil_inner_boundary * i_theta;
    }
}



int Smoother::nnz_radial_Asc(const int i_theta) const {
    assert(i_theta >= 0 && i_theta < grid_.ntheta());

    const int size_stencil_next_circluar_smoothing = 2;
    const int size_stencil_interior = 3;
    const int size_stencil_next_outer_boundary = 2;
    const int size_stencil_outer_boundary = 1;

    assert(grid_.lengthSmootherRadial() >= 3);

    return 
        size_stencil_next_circluar_smoothing +
        (grid_.lengthSmootherRadial()-3) * size_stencil_interior +
        size_stencil_next_outer_boundary + 
        size_stencil_outer_boundary;
}



int Smoother::ptr_nz_index_radial_Asc(const int i_r, const int i_theta) const {
    assert(i_theta >= 0 && i_theta < grid_.ntheta());

    const int size_stencil_next_circluar_smoothing = 2;
    const int size_stencil_interior = 3;
    const int size_stencil_next_outer_boundary = 2;
    const int size_stencil_outer_boundary = 1;

    assert(grid_.lengthSmootherRadial() >= 3);
    assert(grid_.numberSmootherCircles() >= 2);

    const int numberSmootherCircles = grid_.numberSmootherCircles();

    if(i_r > numberSmootherCircles && i_r < grid_.nr()-2){
        return 
            size_stencil_next_circluar_smoothing + 
            (i_r-numberSmootherCircles-1) * size_stencil_interior;
    }
    else if(i_r == numberSmootherCircles){
        return 0;
    }
    else if(i_r == grid_.nr()-2){
        return 
            size_stencil_next_circluar_smoothing +
            (grid_.lengthSmootherRadial()-3) * size_stencil_interior;
    }
    else if(i_r == grid_.nr()-1){
        return 
            size_stencil_next_circluar_smoothing + 
            (grid_.lengthSmootherRadial()-3) * size_stencil_interior + 
            size_stencil_next_outer_boundary;
    }
    throw std::out_of_range("Invalid index for stencil");
}
