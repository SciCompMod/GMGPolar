#include "../../include/Operator/stencil.h"

Stencil::Stencil(std::initializer_list<int> init) : values{} {
    std::copy(init.begin(), init.end(), values.begin());
}

int Stencil::operator[](StencilType type) const {
    return values[static_cast<size_t>(type)];
}


namespace matrixA_Stencil {
    const Stencil& get_stencil(const PolarGrid& grid, int i_r, int DirBC_Interior) {
        static const Stencil stencil_DB = 
            {-1, -1, -1,
            -1,  0, -1,
            -1, -1, -1};
            
        static const Stencil stencil_across = 
            {-1, 4, 6,
            1, 0, 2,
            -1, 3, 5};
            
        static const Stencil stencil_interior = 
            {7, 4, 8,
            1, 0, 2,
            5, 3, 6};

        if (i_r > 0 && i_r < grid.nr() - 1) {
            return stencil_interior;
        } else if (i_r == grid.nr() - 1 || (i_r == 0 && DirBC_Interior)) {
            return stencil_DB;
        } else if (i_r == 0 && !DirBC_Interior) {
            return stencil_across;
        }
        
        throw std::out_of_range("Invalid index for stencil");
    }


    int nnz_matrixA(const PolarGrid& grid, const int DirBC_Interior){
        const int inner_circle_stencil = DirBC_Interior ? 1 : 7;
        const int interior_stencil = 9;
        const int outer_circle_stencil = 1;
        return inner_circle_stencil * grid.ntheta() + 
            interior_stencil * (grid.nr()-2) * grid.ntheta() +
            outer_circle_stencil * grid.ntheta();
    }

    int ptr_nz_index_matrixA(const PolarGrid& grid, const int i_r, const int i_theta, const int DirBC_Interior){
        const int inner_circle_stencil = DirBC_Interior ? 1 : 7;
        const int interior_stencil = 9;
        const int outer_circle_stencil = 1;
        if(0 < i_r && i_r < grid.numberSmootherCircles()){
            // Interior. Circle index section
            const int number_prior_inner_boundary_nodes = grid.ntheta();
            const int number_prior_interior_circle_nodes = (i_r-1) * grid.ntheta() + i_theta;
            return inner_circle_stencil * number_prior_inner_boundary_nodes + 
                interior_stencil * number_prior_interior_circle_nodes;
        } else if(grid.numberSmootherCircles() <= i_r && i_r < grid.nr()-1){
            // Interior. Radial index section
            const int number_prior_inner_boundary_nodes = grid.ntheta();
            const int number_prior_interior_circle_nodes = grid.ntheta() * (grid.numberSmootherCircles()-1);
            const int number_prior_interior_radial_nodes = i_theta * (grid.lengthSmootherRadial()-1) + i_r - grid.numberSmootherCircles();
            const int number_prior_outer_boundary_nodes = i_theta;
            return inner_circle_stencil * number_prior_inner_boundary_nodes + 
                interior_stencil * number_prior_interior_circle_nodes + 
                interior_stencil * number_prior_interior_radial_nodes + 
                outer_circle_stencil * number_prior_outer_boundary_nodes;
        } else if(i_r == 0){
            // Boundary on the interior, inner circle.
            const int number_prior_inner_boundary_nodes = i_theta;
            return inner_circle_stencil * number_prior_inner_boundary_nodes;
        } else{
            // Boundary on the outside, outer circle.
            const int number_prior_inner_boundary_nodes = grid.ntheta();
            const int number_prior_interior_circle_nodes = grid.ntheta() * (grid.numberSmootherCircles()-1);
            const int number_prior_interior_radial_nodes = (i_theta+1) * (grid.lengthSmootherRadial()-1);
            const int number_prior_outer_boundary_nodes = i_theta;
            return inner_circle_stencil * number_prior_inner_boundary_nodes + 
                interior_stencil * number_prior_interior_circle_nodes + 
                interior_stencil * number_prior_interior_radial_nodes + 
                outer_circle_stencil * number_prior_outer_boundary_nodes;
        }
    }
}



namespace symmetric_matrixA_Stencil {
    const Stencil& get_stencil(const PolarGrid& grid, int i_r, int DirBC_Interior) {
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
            
        static const Stencil stencil_across = 
            {-1, 4, 6,
            1, 0, 2,
            -1, 3, 5};
            
        static const Stencil stencil_interior = 
            {7, 4, 8,
            1, 0, 2,
            5, 3, 6};

        if ((i_r > 1 && i_r < grid.nr() - 2) || i_r == 1 && !DirBC_Interior) {
            return stencil_interior;
        } else if (i_r == grid.nr() - 1 || (i_r == 0 && DirBC_Interior)) {
            return stencil_DB;
        } else if (i_r == 0 && !DirBC_Interior) {
            return stencil_across;
        } else if(i_r == 1 && DirBC_Interior){
            return stencil_next_inner_DB;
        } else if(i_r == grid.nr() - 2){
            return stencil_next_outer_DB;
        }
        
        throw std::out_of_range("Invalid index for stencil");
    }


    int nnz_matrixA(const PolarGrid& grid, const int DirBC_Interior){
        const int inner_circle_stencil = DirBC_Interior ? 1 : 7;
        const int next_inner_circle_stencil = DirBC_Interior ? 6 : 9;
        const int interior_stencil = 9;
        const int next_outer_circle_stencil = 6;
        const int outer_circle_stencil = 1;
        return inner_circle_stencil * grid.ntheta() +
            next_inner_circle_stencil * grid.ntheta() +
            interior_stencil * (grid.nr()-4) * grid.ntheta() +
            next_outer_circle_stencil * grid.ntheta() +
            outer_circle_stencil * grid.ntheta();
    }


    int ptr_nz_index_matrixA(const PolarGrid& grid, const int i_r, const int i_theta, const int DirBC_Interior){
        const int inner_circle_stencil = DirBC_Interior ? 1 : 7;
        const int next_inner_circle_stencil = DirBC_Interior ? 6 : 9;
        const int interior_stencil = 9;
        const int next_outer_circle_stencil = 6;
        const int outer_circle_stencil = 1;
        if(1 < i_r && i_r < grid.numberSmootherCircles()){
            // Interior. Circle index section
            const int number_prior_inner_boundary_nodes = grid.ntheta();
            const int number_prior_next_inner_boundary_nodes = grid.ntheta();
            const int number_prior_interior_circle_nodes = (i_r-2) * grid.ntheta() + i_theta;
            return inner_circle_stencil * number_prior_inner_boundary_nodes + 
                next_inner_circle_stencil * number_prior_next_inner_boundary_nodes +
                interior_stencil * number_prior_interior_circle_nodes;
        } else if(grid.numberSmootherCircles() <= i_r && i_r < grid.nr()-2){
            // Interior. Radial index section
            const int number_prior_inner_boundary_nodes = grid.ntheta();
            const int number_prior_next_inner_boundary_nodes = grid.ntheta();
            const int number_prior_interior_circle_nodes = grid.ntheta() * (grid.numberSmootherCircles()-2);
            const int number_prior_interior_radial_nodes = i_theta * (grid.lengthSmootherRadial()-2) + i_r - grid.numberSmootherCircles();
            const int number_prior_next_outer_boundary_nodes = i_theta;
            const int number_prior_outer_boundary_nodes = i_theta;
            return inner_circle_stencil * number_prior_inner_boundary_nodes + 
                next_inner_circle_stencil * number_prior_next_inner_boundary_nodes +
                interior_stencil * number_prior_interior_circle_nodes + 
                interior_stencil * number_prior_interior_radial_nodes + 
                next_outer_circle_stencil * number_prior_next_outer_boundary_nodes +
                outer_circle_stencil * number_prior_outer_boundary_nodes;
        } else if(i_r == 0){
            // Boundary on the interior, inner circle.
            const int number_prior_inner_boundary_nodes = i_theta;
            return inner_circle_stencil * number_prior_inner_boundary_nodes;
        } else if(i_r == 1){
            // Node next to interior boundary
            const int number_prior_inner_boundary_nodes = grid.ntheta();
            const int number_prior_next_inner_boundary_nodes = i_theta;
            return  inner_circle_stencil * number_prior_inner_boundary_nodes + 
                next_inner_circle_stencil * number_prior_next_inner_boundary_nodes;
        } else if(i_r == grid.nr() - 2){
            // Node next to outer boundary
            const int number_prior_inner_boundary_nodes = grid.ntheta();
            const int number_prior_next_inner_boundary_nodes = grid.ntheta();
            const int number_prior_interior_circle_nodes = grid.ntheta() * (grid.numberSmootherCircles()-2);
            const int number_prior_interior_radial_nodes = (i_theta+1) * (grid.lengthSmootherRadial()-2);
            const int number_prior_next_outer_boundary_nodes = i_theta;
            const int number_prior_outer_boundary_nodes = i_theta;
            return inner_circle_stencil * number_prior_inner_boundary_nodes + 
                next_inner_circle_stencil * number_prior_next_inner_boundary_nodes +
                interior_stencil * number_prior_interior_circle_nodes + 
                interior_stencil * number_prior_interior_radial_nodes + 
                next_outer_circle_stencil * number_prior_next_outer_boundary_nodes +
                outer_circle_stencil * number_prior_outer_boundary_nodes;
        } else{ // i_r == grid.nr() - 1
            // Boundary on the outside, outer circle.
            const int number_prior_inner_boundary_nodes = grid.ntheta();
            const int number_prior_next_inner_boundary_nodes = grid.ntheta();
            const int number_prior_interior_circle_nodes = grid.ntheta() * (grid.numberSmootherCircles()-2);
            const int number_prior_interior_radial_nodes = (i_theta+1) * (grid.lengthSmootherRadial()-2);
            const int number_prior_next_outer_boundary_nodes = (i_theta + 1);
            const int number_prior_outer_boundary_nodes = i_theta;
            return inner_circle_stencil * number_prior_inner_boundary_nodes + 
                next_inner_circle_stencil * number_prior_next_inner_boundary_nodes +
                interior_stencil * number_prior_interior_circle_nodes + 
                interior_stencil * number_prior_interior_radial_nodes + 
                next_outer_circle_stencil * number_prior_next_outer_boundary_nodes +
                outer_circle_stencil * number_prior_outer_boundary_nodes;
        }
    }
}