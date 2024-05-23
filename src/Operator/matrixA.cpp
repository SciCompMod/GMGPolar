#include "../../include/Operator/operator.h"

#include "../../include/common/constants.h"

void Operator::applyATake0(const Level& onLevel, Vector<scalar_t>& result, const Vector<scalar_t>& x) const{
    const PolarGrid& grid = onLevel.grid();
    const ExactFunctions& exactFuncs = onLevel.exactFunctions();

    assert(x.size() == grid.number_of_nodes());
    assert(result.size() == grid.number_of_nodes());

    #pragma omp parallel for
    for(int index = 0; index < grid.number_of_nodes(); index ++){
        double arr, arr_left, arr_right, arr_bottom, arr_top;
        double att, att_left, att_right, att_bottom, att_top;
        double art, art_left, art_right, art_bottom, art_top;
        double detDFinv, detDFinv_left, detDFinv_right, detDFinv_bottom, detDFinv_top;

        double coeff_alpha, coeff_beta;
        double sin_theta, cos_theta;

        MultiIndex node = grid.multiindex(index);
        Point coords = grid.polar_coordinates(node);

        coeff_alpha = (*alpha_)(coords[0]);
        coeff_beta = (*beta_)(coords[0]);

        sin_theta = sin(coords[1]); cos_theta = cos(coords[1]);

        arr_att_art(coords[0], coords[1], sin_theta, cos_theta, coeff_alpha, arr, att, art, detDFinv);

        std::array<std::pair<scalar_t,scalar_t>, space_dimension> neighbor_distance;
        grid.adjacent_neighbor_distances(node, neighbor_distance);

        std::array<std::pair<int,int>, space_dimension> neighbors;
        grid.adjacent_neighbors_of(node, neighbors);

        if(neighbors[0].second == -1){
            // Dirichlet Boundary
            result[index] = x[index];
        } else{

            /* Gather arr, art, att values from adjacent neighbors */
            if(neighbors[0].first != -1){
                MultiIndex left_node = grid.multiindex(neighbors[0].first);
                Point left_coords = grid.polar_coordinates(left_node); 
                double coeff_alpha_left = (*alpha_)(left_coords[0]);
                sin_theta = sin(left_coords[1]); cos_theta = cos(left_coords[1]);
                arr_att_art(left_coords[0], left_coords[1], sin_theta, cos_theta, coeff_alpha_left, arr_left, att_left, art_left, detDFinv_left);     
            }else{
                MultiIndex across_origin_node(0, (node[1] + grid.ntheta() / 2) % grid.ntheta());
                Point across_origin_coords = grid.polar_coordinates(across_origin_node); 
                scalar_t coeff_alpha_left = (*alpha_)(across_origin_coords[0]);
                sin_theta = sin(across_origin_coords[1]); cos_theta = cos(across_origin_coords[1]);
                arr_att_art(across_origin_coords[0], across_origin_coords[1], sin_theta, cos_theta, coeff_alpha_left, arr_left, att_left, art_left, detDFinv_left);  
            }

            // Right
            if(neighbors[0].second != -1){
                MultiIndex right_node = grid.multiindex(neighbors[0].second);
                Point right_coords = grid.polar_coordinates(right_node); 
                scalar_t coeff_alpha_right = exactFuncs.coeffs1(right_coords[0], Rmax_);
                sin_theta = sin(right_coords[1]); cos_theta = cos(right_coords[1]);
                arr_att_art(right_coords[0], right_coords[1], sin_theta, cos_theta, coeff_alpha_right, arr_right, att_right, art_right, detDFinv_right);     
            }
            // Bottom
            if(neighbors[1].first != -1){
                MultiIndex bottom_node = grid.multiindex(neighbors[1].first);
                Point bottom_coords = grid.polar_coordinates(bottom_node); 
                scalar_t coeff_alpha_bottom = (*alpha_)(bottom_coords[0]);
                sin_theta = sin(bottom_coords[1]); cos_theta = cos(bottom_coords[1]);
                arr_att_art(bottom_coords[0], bottom_coords[1], sin_theta, cos_theta, coeff_alpha_bottom, arr_bottom, att_bottom, art_bottom, detDFinv_bottom);   
            }
            // Top
            if(neighbors[1].second != -1){
                MultiIndex top_node = grid.multiindex(neighbors[1].second);
                Point top_coords = grid.polar_coordinates(top_node); 
                scalar_t coeff_alpha_top = (*alpha_)(top_coords[0]);
                sin_theta = sin(top_coords[1]); cos_theta = cos(top_coords[1]);
                arr_att_art(top_coords[0], top_coords[1], sin_theta, cos_theta, coeff_alpha_top, arr_top, att_top, art_top, detDFinv_top);   
            }

            scalar_t h1 = neighbor_distance[0].first;
            scalar_t h2 = neighbor_distance[0].second;
            scalar_t k1 = neighbor_distance[1].first;
            scalar_t k2 = neighbor_distance[1].second;

            // beta_{s,t} / f_{s,t} //
            scalar_t value = (h1+h2)*(k1+k2) * coeff_beta * fabs(detDFinv) / 4  * x[index];

            // ------ //
            // Center //
            // ------ //
            // Left
            if(neighbors[0].first != -1){
                value += 0.5 * (k1 + k2) / h1 * (arr + arr_left) * x[index];
            }else{
                MultiIndex across_origin_node(0, (node[1] + grid.ntheta() / 2) % grid.ntheta());
                value += 0.5 * (k1 + k2) / (2 * grid.radius(0)) * (arr + arr_left) * x[index];
            }
            // Right
            if(neighbors[0].second != -1){
                value += 0.5 * (k1 + k2) / h2 * (arr + arr_right) * x[index];
            }
            // Bottom
            if(neighbors[1].first != -1){
                value += 0.5 * (h1 + h2) / k1 * (att + att_bottom) * x[index];
            }
            // Top
            if(neighbors[1].second != -1){
                value += 0.5 * (h1 + h2) / k2 * (att + att_top) * x[index];
            }

            // ---- //
            // Left //
            if(neighbors[0].first != -1){
                value += -0.5 * (k1 + k2) / h1 * (arr + arr_left) * x[neighbors[0].first];
            }else{
                MultiIndex across_origin_node(0, (node[1] + grid.ntheta() / 2) % grid.ntheta());
                value += -0.5 * (k1 + k2) / (2 * grid.radius(0)) * (arr + arr_left) * x[grid.index(across_origin_node)];
            }
            // ----- //
            // Right //
            if(neighbors[0].second != -1){
                value += -0.5 * (k1 + k2) / h2 * (arr + arr_right) * x[neighbors[0].second];
            }
            // ------ //
            // Bottom //
            if(neighbors[1].first != -1){
                value += -0.5 * (h1 + h2) / k1 * (att + att_bottom) * x[neighbors[1].first];
            }
            // --- //
            // Top //
            if(neighbors[1].second != -1){
                value += -0.5 * (h1 + h2) / k2 * (att + att_top) * x[neighbors[1].second];
            }


            grid.diagonal_neighbors_of(node, neighbors);

            // Bottom Left
            if(neighbors[0].first != -1){
                value += -0.25 * (art_bottom + art_left) * x[neighbors[0].first];
            }
            // Bottom Right
            if(neighbors[0].second != -1){
                value += 0.25 * (art_bottom + art_right) * x[neighbors[0].second];
            }
            // Top Left
            if(neighbors[1].first != -1){
                value += 0.25 * (art_top + art_left) * x[neighbors[1].first];
            }
            // Top Right
            if(neighbors[1].second != -1){
                value += -0.25 * (art_top + art_right) * x[neighbors[1].second];
            }

            result[index] = value;
        }
    }
}



#define NODE_APPLY_A_GIVE(i_r, i_theta, grid, result, x, arr, att, art, coeff_beta, detDFinv) \
do { \
    /* Node in the interior */ \
    if (i_r > 0 && i_r < grid.nr() - 2) { \
        scalar_t h1 = grid.r_dist(i_r-1); \
        scalar_t h2 = grid.r_dist(i_r); \
        scalar_t k1 = grid.theta_dist(i_theta-1); \
        scalar_t k2 = grid.theta_dist(i_theta); \
        \
        double coeff1 = 0.5*(k1+k2)/h1; \
        double coeff2 = 0.5*(k1+k2)/h2; \
        double coeff3 = 0.5*(h1+h2)/k1; \
        double coeff4 = 0.5*(h1+h2)/k2; \
        /* Fill result(i,j) */ \
        result[grid.index(i_r,i_theta)] += \
            0.25 * (h1+h2)*(k1+k2) * coeff_beta * fabs(detDFinv) * x[grid.index(i_r,i_theta)] /* beta_{i,j} */ \
            - coeff1 * arr * x[grid.index(i_r-1,i_theta)] /* Left */ \
            - coeff2 * arr * x[grid.index(i_r+1,i_theta)] /* Right */ \
            - coeff3 * att * x[grid.index(i_r,i_theta-1)] /* Bottom */ \
            - coeff4 * att * x[grid.index(i_r,i_theta+1)] /* Top */ \
            /* Center: (Left, Right, Bottom, Top) */ \
            + ((coeff1 + coeff2) * arr + (coeff3 + coeff4) * att) * x[grid.index(i_r,i_theta)]; \
        /* Fill result(i-1,j) */ \
        result[grid.index(i_r-1,i_theta)] += \
            - coeff1 * arr * x[grid.index(i_r,i_theta)] /* Right */ \
            + coeff1 * arr * x[grid.index(i_r-1,i_theta)] /* Center: (Right) */ \
            - 0.25 * art * x[grid.index(i_r,i_theta+1)] /* Top Right */ \
            + 0.25 * art * x[grid.index(i_r,i_theta-1)]; /* Bottom Right */ \
        /* Fill result(i+1,j) */ \
        result[grid.index(i_r+1,i_theta)] += \
            - coeff2 * arr * x[grid.index(i_r,i_theta)] /* Left */ \
            + coeff2 * arr * x[grid.index(i_r+1,i_theta)] /* Center: (Left) */ \
            + 0.25 * art * x[grid.index(i_r,i_theta+1)] /* Top Left */ \
            - 0.25 * art * x[grid.index(i_r,i_theta-1)]; /* Bottom Left */ \
        /* Fill result(i,j-1) */ \
        result[grid.index(i_r,i_theta-1)] += \
            - coeff3 * att * x[grid.index(i_r,i_theta)] /* Top */ \
            + coeff3 * att * x[grid.index(i_r,i_theta-1)] /* Center: (Top) */ \
            - 0.25 * art * x[grid.index(i_r+1,i_theta)] /* Top Right */ \
            + 0.25 * art * x[grid.index(i_r-1,i_theta)]; /* Top Left */ \
        /* Fill result(i,j+1) */ \
        result[grid.index(i_r,i_theta+1)] += \
            - coeff4 * att * x[grid.index(i_r,i_theta)] /* Bottom */ \
            + coeff4 * att * x[grid.index(i_r,i_theta+1)] /* Center: (Bottom) */ \
            + 0.25 * art * x[grid.index(i_r+1,i_theta)] /* Bottom Right */ \
            - 0.25 * art * x[grid.index(i_r-1,i_theta)]; /* Bottom Left */ \
        \
    /* Node in the inner circle */ \
    } else if (i_r == 0) { \
        /* h1 gets replaced with 2 * R0. */ \
        /* (i_r-1,i_theta) gets replaced with (i_r, i_theta + (grid.ntheta()>>1)). */ \
        /* Some more adjustments from the changing the 9-point stencil to the artifical 7-point stencil. */ \
        scalar_t h1 = 2 * grid.radius(0); \
        scalar_t h2 = grid.r_dist(i_r); \
        scalar_t k1 = grid.theta_dist(i_theta-1); \
        scalar_t k2 = grid.theta_dist(i_theta); \
        \
        double coeff1 = 0.5*(k1+k2)/h1; \
        double coeff2 = 0.5*(k1+k2)/h2; \
        double coeff3 = 0.5*(h1+h2)/k1; \
        double coeff4 = 0.5*(h1+h2)/k2; \
        /* Fill result(i,j) */ \
        result[grid.index(i_r,i_theta)] += \
            0.25 * (h1+h2)*(k1+k2) * coeff_beta * fabs(detDFinv) * x[grid.index(i_r,i_theta)] /* beta_{i,j} */ \
            - coeff1 * arr * x[grid.index(i_r, i_theta + (grid.ntheta()>>1))] /* Left */ \
            - coeff2 * arr * x[grid.index(i_r+1,i_theta)] /* Right */ \
            - coeff3 * att * x[grid.index(i_r,i_theta-1)] /* Bottom */ \
            - coeff4 * att * x[grid.index(i_r,i_theta+1)] /* Top */ \
            /* Center: (Left, Right, Bottom, Top) */ \
            + ((coeff1 + coeff2) * arr + (coeff3 + coeff4) * att) * x[grid.index(i_r,i_theta)]; \
        /* Fill result(i-1,j) */ \
        result[grid.index(i_r, i_theta + (grid.ntheta()>>1))] += \
            - coeff1 * arr * x[grid.index(i_r,i_theta)] /* Right */ \
            + coeff1 * arr * x[grid.index(i_r, i_theta + (grid.ntheta()>>1))]; /* Center: (Right) */ \
        /*  + 0.25 * art * x[grid.index(i_r,i_theta+1)] // Top Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */ \
        /*  - 0.25 * art * x[grid.index(i_r,i_theta-1)]; // Bottom Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */ \
        /* Fill result(i+1,j) */ \
        result[grid.index(i_r+1,i_theta)] += \
            - coeff2 * arr * x[grid.index(i_r,i_theta)] /* Left */ \
            + coeff2 * arr * x[grid.index(i_r+1,i_theta)] /* Center: (Left) */ \
            + 0.25 * art * x[grid.index(i_r,i_theta+1)] /* Top Left */ \
            - 0.25 * art * x[grid.index(i_r,i_theta-1)]; /* Bottom Left */ \
        /* Fill result(i,j-1) */ \
        result[grid.index(i_r,i_theta-1)] += \
            - coeff3 * att * x[grid.index(i_r,i_theta)] /* Top */ \
            + coeff3 * att * x[grid.index(i_r,i_theta-1)] /* Center: (Top) */ \
            - 0.25 * art * x[grid.index(i_r+1,i_theta)]; /* Top Right */ \
        /*  + 0.25 * art * x[grid.index(i_r-1,i_theta)]; // Top Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */ \
        /* Fill result(i,j+1) */ \
        result[grid.index(i_r,i_theta+1)] += \
            - coeff4 * att * x[grid.index(i_r,i_theta)] /* Bottom */ \
            + coeff4 * att * x[grid.index(i_r,i_theta+1)] /* Center: (Bottom) */ \
            + 0.25 * art * x[grid.index(i_r+1,i_theta)]; /* Bottom Right */ \
        /*  - 0.25 * art * x[grid.index(i_r-1,i_theta)]; // Bottom Left: REMOVED DUE TO ARTIFICAL 7 POINT STENCIL */ \
        \
    /* Node in the 2nd outer circle */ \
    } else if (i_r == grid.nr() - 2) { \
        scalar_t h1 = grid.r_dist(i_r-1); \
        scalar_t h2 = grid.r_dist(i_r); \
        scalar_t k1 = grid.theta_dist(i_theta-1); \
        scalar_t k2 = grid.theta_dist(i_theta); \
        \
        double coeff1 = 0.5*(k1+k2)/h1; \
        double coeff2 = 0.5*(k1+k2)/h2; \
        double coeff3 = 0.5*(h1+h2)/k1; \
        double coeff4 = 0.5*(h1+h2)/k2; \
        /* Fill result(i,j) */ \
        result[grid.index(i_r,i_theta)] += \
            (h1+h2)*(k1+k2) * coeff_beta * fabs(detDFinv) / 4 * x[grid.index(i_r,i_theta)] /* beta_{i,j} */ \
            - coeff1 * arr * x[grid.index(i_r-1,i_theta)] /* Left */ \
            - coeff2 * arr * x[grid.index(i_r+1,i_theta)] /* Right */ \
            - coeff3 * att * x[grid.index(i_r,i_theta-1)] /* Bottom */ \
            - coeff4 * att * x[grid.index(i_r,i_theta+1)] /* Top */ \
            + ((coeff1 + coeff2) * arr + (coeff3 + coeff4) * att) * x[grid.index(i_r,i_theta)]; /* Center: (Left, Right, Bottom, Top) */ \
        /* Fill result(i-1,j) */ \
        result[grid.index(i_r-1,i_theta)] += \
            - coeff1 * arr * x[grid.index(i_r,i_theta)] /* Right */ \
            + coeff1 * arr * x[grid.index(i_r-1,i_theta)] /* Center: (Right) */ \
            - 0.25 * art * x[grid.index(i_r,i_theta+1)] /* Top Right */ \
            + 0.25 * art * x[grid.index(i_r,i_theta-1)]; /* Bottom Right */ \
        /* Don't give to the dirichlet boundary! */ \
        /* Fill result(i+1,j) */ \
        /* result[grid.index(i_r+1,i_theta)] += */ \
        /*     - coeff2 * arr * x[grid.index(i_r,i_theta)] // Left */ \
        /*     + coeff2 * arr * x[grid.index(i_r+1,i_theta)] // Center: (Left) */ \
        /*     + 0.25 * art * x[grid.index(i_r,i_theta+1)] // Top Left */ \
        /*     - 0.25 * art * x[grid.index(i_r,i_theta-1)]; // Bottom Left */ \
        /* Fill result(i,j-1) */ \
        result[grid.index(i_r,i_theta-1)] += \
            - coeff3 * att * x[grid.index(i_r,i_theta)] /* Top */ \
            + coeff3 * att * x[grid.index(i_r,i_theta-1)] /* Center: (Top) */ \
            - 0.25 * art * x[grid.index(i_r+1,i_theta)] /* Top Right */ \
            + 0.25 * art * x[grid.index(i_r-1,i_theta)]; /* Top Left */ \
        /* Fill result(i,j+1) */ \
        result[grid.index(i_r,i_theta+1)] += \
            - coeff4 * att * x[grid.index(i_r,i_theta)] /* Bottom */ \
            + coeff4 * att * x[grid.index(i_r,i_theta+1)] /* Center: (Bottom) */ \
            + 0.25 * art * x[grid.index(i_r+1,i_theta)] /* Bottom Right */ \
            - 0.25 * art * x[grid.index(i_r-1,i_theta)]; /* Bottom Left */ \
        \
    /* Node in the outer circle */ \
    } else if (i_r == grid.nr() - 1) { \
        /* Dirichlet boundary */ \
        result[grid.index(i_r,i_theta)] = x[grid.index(i_r,i_theta)]; \
        /* Give value to the interior nodes! */ \
        scalar_t h1 = grid.r_dist(i_r-1); \
        scalar_t k1 = grid.theta_dist(i_theta-1); \
        scalar_t k2 = grid.theta_dist(i_theta); \
        \
        double coeff1 = 0.5*(k1+k2)/h1; \
        /* Fill result(i-1,j) */ \
        result[grid.index(i_r-1,i_theta)] += \
            - coeff1 * arr * x[grid.index(i_r,i_theta)] /* Right */ \
            + coeff1 * arr * x[grid.index(i_r-1,i_theta)] /* Center: (Right) */ \
            - 0.25 * art * x[grid.index(i_r,i_theta+1)] /* Top Right */ \
            + 0.25 * art * x[grid.index(i_r,i_theta-1)]; /* Bottom Right */ \
        \
    } \
} while(0)


void Operator::applyAGive(const Level& onLevel, Vector<scalar_t>& result, const Vector<scalar_t>& x) const{
    const PolarGrid& grid = onLevel.grid();
    const ExactFunctions& exactFuncs = onLevel.exactFunctions();

    assert(x.size() == grid.number_of_nodes());
    assert(result.size() == grid.number_of_nodes());

    assign(result, 0.0);

    // ---------------------- //
    // OpenMP Parallelization //
    // ---------------------- //
    const int numThreads = omp_get_max_threads();
    omp_set_num_threads(numThreads);
    const int minimalChunkSize = 4;
    const int zone = 2;

    // Distribute Tasks to each thread
    TaskDistribution CircleSmootherTasks(grid.numberSmootherCircles(), minimalChunkSize, numThreads);
    TaskDistribution RadialSmootherTasks(grid.ntheta(), minimalChunkSize, numThreads);

    // The acroos origin neighbor nodes need to be in the same TaskDistribution to prevent race condition.
    assert(grid.numberSmootherCircles() >= 1);

    #pragma omp parallel num_threads(numThreads)
    {   
        double r, theta;
        double sin_theta, cos_theta;
        double arr, att, art;
        double coeff_alpha, coeff_beta;
        double detDFinv;

        const int threadID = omp_get_thread_num();
        // ---------------------------------------------------------- //
        // Take care of the separation strips of the circular smoother //
        // ---------------------------------------------------------- //
        const int i_r_start = CircleSmootherTasks.getStart(threadID);
        const int i_r_end = CircleSmootherTasks.getEnd(threadID);
        const int i_r_separate = std::min(i_r_end - i_r_start, zone);

        // For loop matches circular access pattern
        for (int i_r = i_r_end - i_r_separate; i_r < i_r_end; i_r++){
            r = grid.radius(i_r);
            coeff_alpha = (*alpha_)(r);
            coeff_beta = (*beta_)(r);
            for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++){
                theta = grid.theta(i_theta);
                sin_theta = sin_theta_[i_theta];
                cos_theta = cos_theta_[i_theta];
                arr_att_art(r, theta, sin_theta, cos_theta, coeff_alpha, arr, att, art, detDFinv);
                NODE_APPLY_A_GIVE(i_r, i_theta, grid, result, x, arr, att, art, coeff_beta, detDFinv);
            }
        }
        
        #pragma omp barrier

        // -------------------------------------------------------- //
        // Take care of the separation strips of the radial smoother //
        // -------------------------------------------------------- //
        const int i_theta_start = RadialSmootherTasks.getStart(threadID);
        const int i_theta_end = RadialSmootherTasks.getEnd(threadID);
        const int i_theta_seperate = std::min(i_theta_end-i_theta_start, zone);

        // For loop matches radial access pattern
        for (int i_theta = i_theta_start; i_theta < i_theta_start + i_theta_seperate; i_theta++){
            theta = grid.theta(i_theta);
            sin_theta = sin_theta_[i_theta];
            cos_theta = cos_theta_[i_theta];
            for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++){
                r = grid.radius(i_r);
                coeff_alpha = (*alpha_)(r);
                coeff_beta = (*beta_)(r);
                arr_att_art(r, theta, sin_theta, cos_theta, coeff_alpha, arr, att, art, detDFinv);
                NODE_APPLY_A_GIVE(i_r, i_theta, grid, result, x, arr, att, art, coeff_beta, detDFinv);
            }
        }

        #pragma omp barrier

        // ------------------------------------------ //
        // Take care of the circular smoother section //
        // ------------------------------------------ //
        // For loop matches circular access pattern
        for (int i_r = i_r_start; i_r < i_r_end - i_r_separate; i_r++){
            r = grid.radius(i_r);
            coeff_alpha = (*alpha_)(r);
            coeff_beta = (*beta_)(r);
            for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++){
                theta = grid.theta(i_theta);
                sin_theta = sin_theta_[i_theta];
                cos_theta = cos_theta_[i_theta];
                arr_att_art(r, theta, sin_theta, cos_theta, coeff_alpha, arr, att, art, detDFinv);
                NODE_APPLY_A_GIVE(i_r, i_theta, grid, result, x, arr, att, art, coeff_beta, detDFinv);
            }
        }

        // ---------------------------------------- //
        // Take care of the radial smoother section //
        // ---------------------------------------- //
        // For loop matches radial access pattern
        for (int i_theta = i_theta_start + i_theta_seperate; i_theta < i_theta_end; i_theta++){
            theta = grid.theta(i_theta);
            sin_theta = sin_theta_[i_theta];
            cos_theta = cos_theta_[i_theta];
            for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++){
                r = grid.radius(i_r);
                coeff_alpha = (*alpha_)(r);
                coeff_beta = (*beta_)(r);
                arr_att_art(r, theta, sin_theta, cos_theta, coeff_alpha, arr, att, art, detDFinv);
                NODE_APPLY_A_GIVE(i_r, i_theta, grid, result, x, arr, att, art, coeff_beta, detDFinv);
            }
        }
    }
}






#define GIVE_CIRCLE(i_r) \
do { \
    r = grid.radius(i_r); \
    coeff_alpha = (*alpha_)(r); \
    coeff_beta = (*beta_)(r); \
    for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++){ \
        theta = grid.theta(i_theta); \
        sin_theta = sin_theta_[i_theta]; \
        cos_theta = cos_theta_[i_theta]; \
        arr_att_art(r, theta, sin_theta, cos_theta, coeff_alpha, arr, att, art, detDFinv); \
        NODE_APPLY_A_GIVE(i_r, i_theta, grid, result, x, arr, att, art, coeff_beta, detDFinv); \
    } \
} while(0)


#define GIVE_RADIAL(i_theta) \
do { \
    theta = grid.theta(i_theta); \
    sin_theta = sin_theta_[i_theta]; \
    cos_theta = cos_theta_[i_theta]; \
    for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++){ \
        r = grid.radius(i_r); \
        coeff_alpha = (*alpha_)(r); \
        coeff_beta = (*beta_)(r); \
        arr_att_art(r, theta, sin_theta, cos_theta, coeff_alpha, arr, att, art, detDFinv); \
        NODE_APPLY_A_GIVE(i_r, i_theta, grid, result, x, arr, att, art, coeff_beta, detDFinv); \
    } \
} while(0)



void Operator::applyAGiveMutex(Level& onLevel, Vector<scalar_t>& result, const Vector<scalar_t>& x) {
    const PolarGrid& grid = onLevel.grid();
    const ExactFunctions& exactFuncs = onLevel.exactFunctions();

    assert(x.size() == grid.number_of_nodes());
    assert(result.size() == grid.number_of_nodes());

    assign(result, 0.0);

    auto start_time = std::chrono::high_resolution_clock::now();

    #pragma omp parallel
    {
        #pragma omp for nowait
        for(int i = 0; i < numberMod1Circles; i++){
            Mod1CircleDepCounter[i] = 2;
        }
        #pragma omp for nowait
        for(int i = 0; i < numberMod2Circles; i++){
            Mod2CircleDepCounter[i] = 2;
        }
        #pragma omp for nowait
        for(int i = 0; i < numberDiv3Radials; i++){
            Mod1RadialDepCounter[i] = 2;
        }
        #pragma omp for
        for(int i = 0; i < numberDiv3Radials; i++){
            Mod2RadialDepCounter[i] = 2;
        }
    }

    if(additional_circle_task == 1){}
    else if(additional_circle_task == 2){
        Mod1CircleDepCounter.back() = 1;
    }
    else if(additional_circle_task == 0){
        Mod1CircleDepCounter.back() = 1;
        Mod2CircleDepCounter.back() = 1;
    }

    int circle_index;

    double r, theta;
    double sin_theta, cos_theta;
    double arr, att, art;
    double coeff_alpha, coeff_beta;
    double detDFinv;

    std::cout<<numberMod0Circles<<std::endl;

    int threads_used = std::min(numberMod0Circles + numberDiv3Radials, omp_get_max_threads());

    #pragma omp parallel num_threads(threads_used) private(circle_index, r, theta, sin_theta, cos_theta, arr, att, art, coeff_alpha, coeff_beta, detDFinv)
    #pragma omp single
    {
        for (int localMod0CircleIndex = 0; localMod0CircleIndex < numberMod0Circles; localMod0CircleIndex++){
            // Mod 0 Circles
            #pragma omp task
            {
                // std::cout<<"Start: Mod 0 Circle "<<localMod0CircleIndex<< ".\n";
                circle_index = 3 * localMod0CircleIndex;
                GIVE_CIRCLE(circle_index);
                // std::this_thread::sleep_for(std::chrono::milliseconds(200*localMod0CircleIndex));
                // // Work finished
                // std::cout<<circle_index<< " Finished: Mod 0 Circle "<<localMod0CircleIndex<< ".\n";

                // Update left Mod 1 dependency 
                int leftlocalMod1CircleIndex = -1;
                if(additional_circle_task != 1 || localMod0CircleIndex < numberMod0Circles-1) leftlocalMod1CircleIndex = localMod0CircleIndex;
                if(leftlocalMod1CircleIndex != -1){
                    Mod1CircleMutexes[leftlocalMod1CircleIndex].lock();
                    Mod1CircleDepCounter[leftlocalMod1CircleIndex] --;
                    int leftlocalMod1Circle_resolved = Mod1CircleDepCounter[leftlocalMod1CircleIndex];
                    Mod1CircleMutexes[leftlocalMod1CircleIndex].unlock();
                    if(leftlocalMod1Circle_resolved == 0){
                        // Mod 1 Circles
                        #pragma omp task
                        {
                            // std::cout<<"Start: Mod 1 Circle "<<leftlocalMod1CircleIndex<< ".\n";
                            circle_index = 3 * (leftlocalMod1CircleIndex) + 1;
                            GIVE_CIRCLE(circle_index);
                            // std::this_thread::sleep_for(std::chrono::milliseconds(100*leftlocalMod1CircleIndex));
                            // // // Work finished
                            // std::cout<<circle_index << " Finished: Mod 1 Circle "<<leftlocalMod1CircleIndex<< ".\n";

                            if(leftlocalMod1CircleIndex == 0){
                                std::cout<<"Start Radials!"<<std::endl;

                                for (int localMod0RadialIndex = 0; localMod0RadialIndex < numberDiv3Radials; localMod0RadialIndex++){
                                    // Mod 0 Radials
                                    #pragma omp task
                                    {
                                        if(localMod0RadialIndex == 0){
                                            GIVE_RADIAL(0);
                                            if(additional_radial_task >= 1){
                                                GIVE_RADIAL(1);
                                            }
                                        }
                                        else{
                                            int radial_index = 3 * localMod0RadialIndex + additional_radial_task;
                                            GIVE_RADIAL(radial_index);
                                        }

                                        int leftlocalMod1RadialIndex = (localMod0RadialIndex + 2) % numberDiv3Radials;
                                        Mod1RadialMutexes[leftlocalMod1RadialIndex].lock();
                                        Mod1RadialDepCounter[leftlocalMod1RadialIndex] --;
                                        int leftlocalMod1Radial_resolved = Mod1RadialDepCounter[leftlocalMod1RadialIndex];
                                        Mod1RadialMutexes[leftlocalMod1RadialIndex].unlock();

                                        if(leftlocalMod1Radial_resolved == 0){
                                            // Mod 1 Radials
                                            #pragma omp task
                                            {
                                                if(leftlocalMod1RadialIndex == 0){
                                                    if(additional_radial_task == 0){
                                                        GIVE_RADIAL(1);
                                                    }
                                                    else if(additional_radial_task == 1){
                                                        GIVE_RADIAL(2);
                                                    }
                                                    else if(additional_radial_task == 2){
                                                        GIVE_RADIAL(2);
                                                        GIVE_RADIAL(3);
                                                    }
                                                }
                                                else{
                                                    int radial_index = 3 * leftlocalMod1RadialIndex + additional_radial_task + 1;
                                                    GIVE_RADIAL(radial_index);
                                                }

                                                int leftlocalMod2RadialIndex = (leftlocalMod1RadialIndex + 2) % numberDiv3Radials;
                                                Mod2RadialMutexes[leftlocalMod2RadialIndex].lock();
                                                Mod2RadialDepCounter[leftlocalMod2RadialIndex] --;
                                                int leftlocalMod2Radial_resolved = Mod2RadialDepCounter[leftlocalMod2RadialIndex];
                                                Mod2RadialMutexes[leftlocalMod2RadialIndex].unlock();

                                                if(leftlocalMod2Radial_resolved == 0){
                                                    // Mod 2 Radials
                                                    #pragma omp task
                                                    {
                                                        int radial_index = 3 * leftlocalMod2RadialIndex + additional_radial_task + 2;
                                                        GIVE_RADIAL(radial_index);
                                                    }
                                                }

                                                int rightlocalMod2RadialIndex = (leftlocalMod1RadialIndex + 2) % numberDiv3Radials;
                                                Mod2RadialMutexes[rightlocalMod2RadialIndex].lock();
                                                Mod2RadialDepCounter[rightlocalMod2RadialIndex] --;
                                                int rightlocalMod2Radial_resolved = Mod2RadialDepCounter[rightlocalMod2RadialIndex];
                                                
                                                Mod2RadialMutexes[rightlocalMod2RadialIndex].unlock();

                                                if(rightlocalMod2Radial_resolved == 0){
                                                    // Mod 2 Radials
                                                    #pragma omp task
                                                    {
                                                        int radial_index = 3 * rightlocalMod2RadialIndex + additional_radial_task + 2;
                                                        GIVE_RADIAL(radial_index);
                                                    }
                                                }
                                            }
                                        }

                                        int rightlocalMod1RadialIndex = localMod0RadialIndex;
                                        Mod1RadialMutexes[rightlocalMod1RadialIndex].lock();
                                        Mod1RadialDepCounter[rightlocalMod1RadialIndex] --;
                                        int rightlocalMod1Radial_resolved = Mod1RadialDepCounter[rightlocalMod1RadialIndex];
                                        
                                        Mod1RadialMutexes[rightlocalMod1RadialIndex].unlock();

                                        if(rightlocalMod1Radial_resolved == 0){
                                        // Mod 1 Radials
                                            #pragma omp task
                                            {
                                                if(rightlocalMod1RadialIndex == 0){
                                                    if(additional_radial_task == 0){
                                                        GIVE_RADIAL(1);
                                                    }
                                                    else if(additional_radial_task == 1){
                                                        GIVE_RADIAL(2);
                                                    }
                                                    else if(additional_radial_task == 2){
                                                        GIVE_RADIAL(2);
                                                        GIVE_RADIAL(3);
                                                    }
                                                }
                                                else{
                                                    int radial_index = 3 * rightlocalMod1RadialIndex + additional_radial_task + 1;
                                                    GIVE_RADIAL(radial_index);
                                                }   

                                                int leftlocalMod2RadialIndex = (rightlocalMod1RadialIndex + 2) % numberDiv3Radials;
                                                Mod2RadialMutexes[leftlocalMod2RadialIndex].lock();
                                                Mod2RadialDepCounter[leftlocalMod2RadialIndex] --;
                                                int leftlocalMod2Radial_resolved = Mod2RadialDepCounter[leftlocalMod2RadialIndex];
                                                
                                                Mod2RadialMutexes[leftlocalMod2RadialIndex].unlock();

                                                if(leftlocalMod2Radial_resolved == 0){
                                                    // Mod 2 Radials
                                                    #pragma omp task
                                                    {
                                                        int radial_index = 3 * leftlocalMod2RadialIndex + additional_radial_task + 2;
                                                        GIVE_RADIAL(radial_index);
                                                    }
                                                }

                                                int rightlocalMod2RadialIndex = (rightlocalMod1RadialIndex + 2) % numberDiv3Radials;
                                                Mod2RadialMutexes[rightlocalMod2RadialIndex].lock();
                                                Mod2RadialDepCounter[rightlocalMod2RadialIndex] --;
                                                int rightlocalMod2Radial_resolved = Mod2RadialDepCounter[rightlocalMod2RadialIndex];
                                                
                                                Mod2RadialMutexes[rightlocalMod2RadialIndex == 0].unlock();

                                                if(rightlocalMod2Radial_resolved){
                                                    // Mod 2 Radials
                                                    #pragma omp task
                                                    {
                                                        int radial_index = 3 * rightlocalMod2RadialIndex + additional_radial_task + 2;
                                                        GIVE_RADIAL(radial_index);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            } // END RADIALS

                            // Update left Mod 2 dependency 
                            int leftlocalMod2CircleIndex = -1;
                            if(additional_circle_task == 0 || leftlocalMod1CircleIndex < numberMod1Circles-1) leftlocalMod2CircleIndex = leftlocalMod1CircleIndex;
                            if(leftlocalMod2CircleIndex != -1){
                                //std::cout<< "T "<<additional_circle_task<< ", " << leftlocalMod1CircleIndex << ", "<<numberMod1Circles<<std::endl;
                                // std::cout<<"S1: "<<leftlocalMod2CircleIndex<<std::endl;
                                Mod2CircleMutexes[leftlocalMod2CircleIndex].lock();
                                Mod2CircleDepCounter[leftlocalMod2CircleIndex] --;
                                int leftlocalMod2Circle_resolved = Mod2CircleDepCounter[leftlocalMod2CircleIndex];
                                
                                Mod2CircleMutexes[leftlocalMod2CircleIndex].unlock();
                                if(leftlocalMod2Circle_resolved == 0){
                                    // Mod 2 Circle Task
                                    #pragma omp task
                                    {
                                        // std::cout<<"Start: Mod 2 Circle "<<leftlocalMod2CircleIndex<< ".\n";
                                        circle_index = 3 * (leftlocalMod2CircleIndex)+2;
                                        GIVE_CIRCLE(circle_index);
                                        // std::this_thread::sleep_for(std::chrono::milliseconds(10*leftlocalMod2CircleIndex));
                                        // // // Work finished
                                        //std::cout<<circle_index<< " Finished: Mod 2 Circle "<<leftlocalMod2CircleIndex<< ".\n";
                                        // std::cout<<3 * (leftlocalMod2CircleIndex)+2<<std::endl;

                                    }
                                }
                            }
                            // Update right Mod 2 dependency 
                            int rightlocalMod2CircleIndex = -1;
                            if(leftlocalMod1CircleIndex > 0) rightlocalMod2CircleIndex = leftlocalMod1CircleIndex-1;
                            if(rightlocalMod2CircleIndex != -1){
                                // std::cout<<"S2 "<<rightlocalMod2CircleIndex<<std::endl;
                                Mod2CircleMutexes[rightlocalMod2CircleIndex].lock();
                                Mod2CircleDepCounter[rightlocalMod2CircleIndex]--;
                                int rightlocalMod2Circle_resolved = Mod2CircleDepCounter[rightlocalMod2CircleIndex];
                                
                                Mod2CircleMutexes[rightlocalMod2CircleIndex].unlock();
                                if(rightlocalMod2Circle_resolved == 0){
                                    // Mod 1 Circle Task
                                    #pragma omp task
                                    {
                                        // std::cout<<"Start: Mod 2 Circle "<<rightlocalMod2CircleIndex<< ".\n";
                                        circle_index = 3 * (rightlocalMod2CircleIndex) + 2;
                                        GIVE_CIRCLE(circle_index);
                                        // std::this_thread::sleep_for(std::chrono::milliseconds(10*rightlocalMod2CircleIndex));
                                        // // // Work finished
                                        //std::cout<<circle_index<<" Finished: Mod 2 Circle "<<rightlocalMod2CircleIndex<< ".\n";

                                        // std::cout<<3 * (rightlocalMod2CircleIndex)+2<<std::endl;
                                    }
                                }
                            }
                        }
                    }
                }
                // Update right Mod 1 dependency 
                int rightlocalMod1CircleIndex = -1;
                if(localMod0CircleIndex > 0) rightlocalMod1CircleIndex = localMod0CircleIndex-1;
                if(rightlocalMod1CircleIndex != -1){
                    Mod1CircleMutexes[rightlocalMod1CircleIndex].lock();
                    Mod1CircleDepCounter[rightlocalMod1CircleIndex] --;
                    int rightlocalMod1Circle_resolved = Mod1CircleDepCounter[rightlocalMod1CircleIndex];
                    
                    Mod1CircleMutexes[rightlocalMod1CircleIndex].unlock();
                    if(rightlocalMod1Circle_resolved == 0){
                        // Mod 1 Circle Task
                        #pragma omp task
                        {
                            // std::cout<<"Start: Mod 1 Circle "<<rightlocalMod1CircleIndex<< ".\n";
                            circle_index = 3 * (rightlocalMod1CircleIndex)+1;
                            GIVE_CIRCLE(circle_index);
                            // std::this_thread::sleep_for(std::chrono::milliseconds(100 * rightlocalMod1CircleIndex));
                            // // // Work finished
                            //std::cout<<circle_index<<" Finished: Mod 1 Circle "<<rightlocalMod1CircleIndex<< ".\n";

                            // std::cout<<3 * (rightlocalMod1CircleIndex)+1<<std::endl;
                            
                            if(rightlocalMod1CircleIndex == 0){
                                std::cout<<"Start Radials!"<<std::endl;

                                for (int localMod0RadialIndex = 0; localMod0RadialIndex < numberDiv3Radials; localMod0RadialIndex++){
                                    // Mod 0 Radials
                                    #pragma omp task
                                    {
                                        if(localMod0RadialIndex == 0){
                                            GIVE_RADIAL(0);
                                            if(additional_radial_task >= 1){
                                                GIVE_RADIAL(1);
                                            }
                                        }
                                        else{
                                            int radial_index = 3 * localMod0RadialIndex + additional_radial_task;
                                            GIVE_RADIAL(radial_index);
                                        }

                                        int leftlocalMod1RadialIndex = (localMod0RadialIndex + 2) % numberDiv3Radials;
                                        Mod1RadialMutexes[leftlocalMod1RadialIndex].lock();
                                        Mod1RadialDepCounter[leftlocalMod1RadialIndex]--;
                                        int leftlocalMod1Radial_resolved = Mod1RadialDepCounter[leftlocalMod1RadialIndex];
                                        
                                        Mod1RadialMutexes[leftlocalMod1RadialIndex].unlock();

                                        if(leftlocalMod1Radial_resolved == 0){
                                            // Mod 1 Radials
                                            #pragma omp task
                                            {
                                                if(leftlocalMod1RadialIndex == 0){
                                                    if(additional_radial_task == 0){
                                                        GIVE_RADIAL(1);
                                                    }
                                                    else if(additional_radial_task == 1){
                                                        GIVE_RADIAL(2);
                                                    }
                                                    else if(additional_radial_task == 2){
                                                        GIVE_RADIAL(2);
                                                        GIVE_RADIAL(3);
                                                    }
                                                }
                                                else{
                                                    int radial_index = 3 * leftlocalMod1RadialIndex + additional_radial_task + 1;
                                                    GIVE_RADIAL(radial_index);
                                                }

                                                int leftlocalMod2RadialIndex = (leftlocalMod1RadialIndex + 2) % numberDiv3Radials;
                                                Mod2RadialMutexes[leftlocalMod2RadialIndex].lock();
                                                Mod2RadialDepCounter[leftlocalMod2RadialIndex]--;
                                                int leftlocalMod2Radial_resolved = Mod2RadialDepCounter[leftlocalMod2RadialIndex];
                                                
                                                Mod2RadialMutexes[leftlocalMod2RadialIndex].unlock();

                                                if(leftlocalMod2Radial_resolved == 0){
                                                    // Mod 2 Radials
                                                    #pragma omp task
                                                    {
                                                        int radial_index = 3 * leftlocalMod2RadialIndex + additional_radial_task + 2;
                                                        GIVE_RADIAL(radial_index);
                                                    }
                                                }

                                                int rightlocalMod2RadialIndex = (leftlocalMod1RadialIndex + 2) % numberDiv3Radials;
                                                Mod2RadialMutexes[rightlocalMod2RadialIndex].lock();
                                                Mod2RadialDepCounter[rightlocalMod2RadialIndex] --;
                                                int rightlocalMod2Radial_resolved = Mod2RadialDepCounter[rightlocalMod2RadialIndex];
                                                
                                                Mod2RadialMutexes[rightlocalMod2RadialIndex].unlock();

                                                if(rightlocalMod2Radial_resolved == 0){
                                                    // Mod 2 Radials
                                                    #pragma omp task
                                                    {
                                                        int radial_index = 3 * rightlocalMod2RadialIndex + additional_radial_task + 2;
                                                        GIVE_RADIAL(radial_index);
                                                    }
                                                }
                                            }
                                        }

                                        int rightlocalMod1RadialIndex = localMod0RadialIndex;
                                        Mod1RadialMutexes[rightlocalMod1RadialIndex].lock();
                                        Mod1RadialDepCounter[rightlocalMod1RadialIndex] --;
                                        int rightlocalMod1Radial_resolved = Mod1RadialDepCounter[rightlocalMod1RadialIndex];
                                        
                                        Mod1RadialMutexes[rightlocalMod1RadialIndex].unlock();

                                        if(rightlocalMod1Radial_resolved == 0){
                                        // Mod 1 Radials
                                            #pragma omp task
                                            {
                                                if(rightlocalMod1RadialIndex == 0){
                                                    if(additional_radial_task == 0){
                                                        GIVE_RADIAL(1);
                                                    }
                                                    else if(additional_radial_task == 1){
                                                        GIVE_RADIAL(2);
                                                    }
                                                    else if(additional_radial_task == 2){
                                                        GIVE_RADIAL(2);
                                                        GIVE_RADIAL(3);
                                                    }
                                                }
                                                else{
                                                    int radial_index = 3 * rightlocalMod1RadialIndex + additional_radial_task + 1;
                                                    GIVE_RADIAL(radial_index);
                                                }   

                                                int leftlocalMod2RadialIndex = (rightlocalMod1RadialIndex + 2) % numberDiv3Radials;
                                                Mod2RadialMutexes[leftlocalMod2RadialIndex].lock();
                                                Mod2RadialDepCounter[leftlocalMod2RadialIndex]--;
                                                int leftlocalMod2Radial_resolved = Mod2RadialDepCounter[leftlocalMod2RadialIndex];
                                                
                                                Mod2RadialMutexes[leftlocalMod2RadialIndex].unlock();

                                                if(leftlocalMod2Radial_resolved == 0){
                                                    // Mod 2 Radials
                                                    #pragma omp task
                                                    {
                                                        int radial_index = 3 * leftlocalMod2RadialIndex + additional_radial_task + 2;
                                                        GIVE_RADIAL(radial_index);
                                                    }
                                                }

                                                int rightlocalMod2RadialIndex = (rightlocalMod1RadialIndex + 2) % numberDiv3Radials;
                                                Mod2RadialMutexes[rightlocalMod2RadialIndex].lock();
                                                Mod2RadialDepCounter[rightlocalMod2RadialIndex]--;
                                                int rightlocalMod2Radial_resolved = Mod2RadialDepCounter[rightlocalMod2RadialIndex];
                                                
                                                Mod2RadialMutexes[rightlocalMod2RadialIndex].unlock();

                                                if(rightlocalMod2Radial_resolved == 0){
                                                    // Mod 2 Radials
                                                    #pragma omp task
                                                    {
                                                        int radial_index = 3 * rightlocalMod2RadialIndex + additional_radial_task + 2;
                                                        GIVE_RADIAL(radial_index);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            } // END RADIALS

                            // Update left Mod 2 dependency 
                            int leftlocalMod2CircleIndex = -1;
                            if(additional_circle_task == 0 || rightlocalMod1CircleIndex < numberMod1Circles-1) leftlocalMod2CircleIndex = rightlocalMod1CircleIndex;
                            if(leftlocalMod2CircleIndex != -1){
                                Mod2CircleMutexes[leftlocalMod2CircleIndex].lock();
                                Mod2CircleDepCounter[leftlocalMod2CircleIndex]--;
                                int leftlocalMod2Circle_resolved = Mod2CircleDepCounter[leftlocalMod2CircleIndex];
                                
                                Mod2CircleMutexes[leftlocalMod2CircleIndex].unlock();
                                if(leftlocalMod2Circle_resolved == 0){
                                    // Mod 2 Circle Task
                                    #pragma omp task
                                    {
                                        // std::cout<<"Start: Mod 2 Circle "<<leftlocalMod2CircleIndex<< ".\n";

                                        circle_index = 3 * (leftlocalMod2CircleIndex)+2;
                                        GIVE_CIRCLE(circle_index);
                                        // std::this_thread::sleep_for(std::chrono::milliseconds(10 * leftlocalMod2CircleIndex));
                                        // // // Work finished
                                        //std::cout<<circle_index<<" Finished: Mod 2 Circle "<<leftlocalMod2CircleIndex<< ".\n";

                                        // std::cout<<3 * (leftlocalMod2CircleIndex)+2<<std::endl;

                            
                                    }
                                }
                            }
                            // Update right Mod 2 dependency 
                            int rightlocalMod2CircleIndex = -1;
                            if(rightlocalMod1CircleIndex > 0) rightlocalMod2CircleIndex = rightlocalMod1CircleIndex-1;
                            if(rightlocalMod2CircleIndex != -1){
                                Mod2CircleMutexes[rightlocalMod2CircleIndex].lock();
                                Mod2CircleDepCounter[rightlocalMod2CircleIndex]--;
                                int rightlocalMod2Circle_resolved = Mod2CircleDepCounter[rightlocalMod2CircleIndex];
                                
                                Mod2CircleMutexes[rightlocalMod2CircleIndex].unlock();
                                if(rightlocalMod2Circle_resolved == 0){
                                    // Mod 1 Circle Task
                                    #pragma omp task
                                    {
                                        // std::cout<<"Start: Mod 2 Circle "<<rightlocalMod2CircleIndex<< ".\n";
                                        circle_index = 3 * (rightlocalMod2CircleIndex) + 2;
                                        GIVE_CIRCLE(circle_index);
                                        // std::this_thread::sleep_for(std::chrono::milliseconds(10 * rightlocalMod2CircleIndex));
                                        // // // Work finished
                                        //std::cout<<circle_index<<" Finished: Mod 2 Circle "<<rightlocalMod2CircleIndex<< ".\n";

                                        // std::cout<<3 * (rightlocalMod2CircleIndex)+2<<std::endl;
                                    
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();

    auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "Mutex time: " << duration_ms.count() << " milliseconds" << std::endl;
}






void Operator::applyAGiveTasks(const Level& onLevel, Vector<scalar_t>& result, const Vector<scalar_t>& x) const{
    const PolarGrid& grid = onLevel.grid();
    const ExactFunctions& exactFuncs = onLevel.exactFunctions();

    assert(x.size() == grid.number_of_nodes());
    assert(result.size() == grid.number_of_nodes());

    assign(result, 0.0);

    // ---------------------- //
    // OpenMP Parallelization //
    // ---------------------- //
    const int numThreads = omp_get_max_threads();
    omp_set_num_threads(numThreads);

    const int S1 = grid.numberSmootherCircles();
    const int S2 = grid.ntheta();
    const int T = S1 + S2;

    int* dep = new int[T];

    const int S2_wait = S2 % 3;
    const int S2_start = std::max(S1-2, 0);

    #pragma omp parallel shared(dep) num_threads(numThreads)
    {
        double r, theta;
        double sin_theta, cos_theta;
        double arr, att, art;
        double coeff_alpha, coeff_beta;
        double detDFinv;
        #pragma omp single
        {
            for(int i_r = S1 - 1; i_r >= 0; i_r -= 3) {
                #pragma omp task shared(dep) firstprivate(i_r) depend(out: dep[i_r])
                {
                    r = grid.radius(i_r);
                    coeff_alpha = (*alpha_)(r);
                    coeff_beta = (*beta_)(r);
                    for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++){
                        theta = grid.theta(i_theta);
                        sin_theta = sin_theta_[i_theta];
                        cos_theta = cos_theta_[i_theta];
                        arr_att_art(r, theta, sin_theta, cos_theta, coeff_alpha, arr, att, art, detDFinv);
                        NODE_APPLY_A_GIVE(i_r, i_theta, grid, result, x, arr, att, art, coeff_beta, detDFinv);
                    }
                }
            }

            for(int i_r = S1 - 2; i_r >= 0; i_r -= 3) {
                #pragma omp task shared(dep) firstprivate(i_r) depend(in: dep[i_r-2], dep[i_r+1]) depend(out: dep[i_r])
                {
                    r = grid.radius(i_r);
                    coeff_alpha = (*alpha_)(r);
                    coeff_beta = (*beta_)(r);
                    for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++){
                        theta = grid.theta(i_theta);
                        sin_theta = sin_theta_[i_theta];
                        cos_theta = cos_theta_[i_theta];
                        arr_att_art(r, theta, sin_theta, cos_theta, coeff_alpha, arr, att, art, detDFinv);
                        NODE_APPLY_A_GIVE(i_r, i_theta, grid, result, x, arr, att, art, coeff_beta, detDFinv);
                    }
                }
            }

            for(int i_r = S1 - 3; i_r >= 0; i_r -= 3) {
                #pragma omp task shared(dep) firstprivate(i_r) depend(in: dep[i_r-2], dep[i_r+1]) depend(out: dep[i_r])
                {
                    r = grid.radius(i_r);
                    coeff_alpha = (*alpha_)(r);
                    coeff_beta = (*beta_)(r);
                    for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++){
                        theta = grid.theta(i_theta);
                        sin_theta = sin_theta_[i_theta];
                        cos_theta = cos_theta_[i_theta];
                        arr_att_art(r, theta, sin_theta, cos_theta, coeff_alpha, arr, att, art, detDFinv);
                        NODE_APPLY_A_GIVE(i_r, i_theta, grid, result, x, arr, att, art, coeff_beta, detDFinv);
                    }
                }
            }

            for(int i_theta = 0; i_theta < S2 - S2_wait; i_theta += 3) {
                #pragma omp task shared(dep) firstprivate(i_theta) depend(out: dep[S1+i_theta]) depend(in: dep[S2_start])
                {
                    theta = grid.theta(i_theta);
                    sin_theta = sin_theta_[i_theta];
                    cos_theta = cos_theta_[i_theta];
                    for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++){
                        r = grid.radius(i_r);
                        coeff_alpha = (*alpha_)(r);
                        coeff_beta = (*beta_)(r);
                        arr_att_art(r, theta, sin_theta, cos_theta, coeff_alpha, arr, att, art, detDFinv);
                        NODE_APPLY_A_GIVE(i_r, i_theta, grid, result, x, arr, att, art, coeff_beta, detDFinv);
                    }
                }
            }

            for(int i_theta = 1; i_theta < S2 - S2_wait; i_theta += 3) {
                #pragma omp task firstprivate(i_theta) depend(out: dep[S1+i_theta]) depend(in: dep[S2_start], dep[S1+i_theta-1], dep[S1+i_theta+2])
                {
                    theta = grid.theta(i_theta);
                    sin_theta = sin_theta_[i_theta];
                    cos_theta = cos_theta_[i_theta];
                    for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++){
                        r = grid.radius(i_r);
                        coeff_alpha = (*alpha_)(r);
                        coeff_beta = (*beta_)(r);
                        arr_att_art(r, theta, sin_theta, cos_theta, coeff_alpha, arr, att, art, detDFinv);
                        NODE_APPLY_A_GIVE(i_r, i_theta, grid, result, x, arr, att, art, coeff_beta, detDFinv);
                    }
                }
            }

            for(int i_theta = 2; i_theta < S2 - S2_wait; i_theta += 3) {
                #pragma omp task shared(dep) firstprivate(i_theta) depend(out: dep[S1+i_theta]) depend(in: dep[S2_start], dep[S1+i_theta-1], dep[S1+i_theta+2])
                {
                    theta = grid.theta(i_theta);
                    sin_theta = sin_theta_[i_theta];
                    cos_theta = cos_theta_[i_theta];
                    for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++){
                        r = grid.radius(i_r);
                        coeff_alpha = (*alpha_)(r);
                        coeff_beta = (*beta_)(r);
                        arr_att_art(r, theta, sin_theta, cos_theta, coeff_alpha, arr, att, art, detDFinv);
                        NODE_APPLY_A_GIVE(i_r, i_theta, grid, result, x, arr, att, art, coeff_beta, detDFinv);
                    }
                }
            }

            if( S2_wait >= 1 ){
                int i_theta = S2 - S2_wait;
                #pragma omp task shared(dep) firstprivate(i_theta) depend(out: dep[S1+i_theta]) depend(in: dep[S2_start],  dep[S1+i_theta-1])
                {
                    theta = grid.theta(i_theta);
                    sin_theta = sin_theta_[i_theta];
                    cos_theta = cos_theta_[i_theta];
                    for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++){
                        r = grid.radius(i_r);
                        coeff_alpha = (*alpha_)(r);
                        coeff_beta = (*beta_)(r);
                        arr_att_art(r, theta, sin_theta, cos_theta, coeff_alpha, arr, att, art, detDFinv);
                        NODE_APPLY_A_GIVE(i_r, i_theta, grid, result, x, arr, att, art, coeff_beta, detDFinv);
                    }
                }
            }

            if(S2_wait >= 2){
                int i_theta = S2 - S2_wait + 1;
                #pragma omp task shared(dep) firstprivate(i_theta) depend(out: dep[S1+i_theta]) depend(in: dep[S2_start],  dep[S1+i_theta-1])
                {
                    theta = grid.theta(i_theta);
                    sin_theta = sin_theta_[i_theta];
                    cos_theta = cos_theta_[i_theta];
                    for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++){
                        r = grid.radius(i_r);
                        coeff_alpha = (*alpha_)(r);
                        coeff_beta = (*beta_)(r);
                        arr_att_art(r, theta, sin_theta, cos_theta, coeff_alpha, arr, att, art, detDFinv);
                        NODE_APPLY_A_GIVE(i_r, i_theta, grid, result, x, arr, att, art, coeff_beta, detDFinv);
                    }
                }
            }
        }
    }
    delete[] dep;
}


























                // if (i_r > 0 && i_r < grid.nr() - 2) {
                //     scalar_t h1 = grid.r_dist(i_r-1); 
                //     scalar_t h2 = grid.r_dist(i_r);
                //     scalar_t k1 = grid.theta_dist(i_theta-1);
                //     scalar_t k2 = grid.theta_dist(i_theta);

                //     double coeff1 = 0.5*(k1+k2)/h1;
                //     double coeff2 = 0.5*(k1+k2)/h2;
                //     double coeff3 = 0.5*(h1+h2)/k1;
                //     double coeff4 = 0.5*(h1+h2)/k2;
                
                //     result[grid.index(i_r,i_theta)] += 
                //         (h1+h2)*(k1+k2) * coeff_beta * fabs(detDFinv) / 4 // f_{i,j}
                //         - coeff1 * arr * x[grid.index(i_r-1,i_theta)] // Left
                //         - coeff2 * arr * x[grid.index(i_r+1,i_theta)] // Right
                //         - coeff3 * att * x[grid.index(i_r,i_theta-1)] // Bottom
                //         - coeff4 * att * x[grid.index(i_r,i_theta+1)] // Top
                //         + ((coeff1 + coeff2) * arr + (coeff3 + coeff4) * att) * x[grid.index(i_r,i_theta)]; // Center: (Left, Right, Bottom, Top)

                //     result[grid.index(i_r-1,i_theta)] += 
                //         - coeff2 * arr * x[grid.index(i_r,i_theta)] // Right
                //         + coeff2 * arr * x[grid.index(i_r-1,i_theta)] // Center: (Right)
                //         - art * x[grid.index(i_r,i_theta+1)] / 4 // Top Right
                //         + art * x[grid.index(i_r,i_theta-1)] / 4; // Bottom Right

                //     result[grid.index(i_r+1,i_theta)] +=
                //         - coeff1 * arr * x[grid.index(i_r,i_theta)] // Left
                //         + coeff1 * arr * x[grid.index(i_r+1,i_theta)] // Center: (Left)
                //         + art * x[grid.index(i_r,i_theta+1)] / 4 // Top Left
                //         - art * x[grid.index(i_r,i_theta-1)] / 4; // Bottom Left

                //     result[grid.index(i_r,i_theta-1)] +=
                //         - coeff4 * att * x[grid.index(i_r,i_theta)] // Top
                //         + coeff4 * att * x[grid.index(i_r,i_theta-1)] // Center: (Top)
                //         - art * x[grid.index(i_r+1,i_theta)] / 4 // Top Right
                //         + art * x[grid.index(i_r-1,i_theta)] / 4; // Top Left

                //     result[grid.index(i_r,i_theta+1)] +=
                //         - coeff3 * att * x[grid.index(i_r,i_theta)] // Bottom
                //         + coeff3 * att * x[grid.index(i_r,i_theta+1)] // Center: (Bottom)
                //         + art * x[grid.index(i_r+1,i_theta)] / 4 // Bottom Right
                //         - art * x[grid.index(i_r-1,i_theta)] / 4; // Bottom Left 

                // } else if (i_r == 0) {
                //     // h1 gets replaced with 2 * R0 and
                //     // (i_r-1,i_theta) gets replaced with (i_r, i_theta + (grid.ntheta()>>1))
                //     scalar_t h1 = 2 * grid.radius(0); 
                //     scalar_t h2 = grid.r_dist(i_r);
                //     scalar_t k1 = grid.theta_dist(i_theta-1);
                //     scalar_t k2 = grid.theta_dist(i_theta);

                //     double coeff1 = 0.5*(k1+k2)/h1;
                //     double coeff2 = 0.5*(k1+k2)/h2;
                //     double coeff3 = 0.5*(h1+h2)/k1;
                //     double coeff4 = 0.5*(h1+h2)/k2;

                //     result[grid.index(i_r,i_theta)] += 
                //         (h1+h2)*(k1+k2) * coeff_beta * fabs(detDFinv) / 4 // f_{i,j}
                //         - coeff1 * arr * x[grid.index(i_r,i_theta+(grid.ntheta()>>1))] // Left
                //         - coeff2 * arr * x[grid.index(i_r+1,i_theta)] // Right
                //         - coeff3 * att * x[grid.index(i_r,i_theta-1)] // Bottom
                //         - coeff4 * att * x[grid.index(i_r,i_theta+1)] // Top
                //         + ((coeff1 + coeff2) * arr + (coeff3 + coeff4) * att) * x[grid.index(i_r,i_theta)]; // Center: (Left, Right, Bottom, Top)

                //     result[grid.index(i_r,i_theta+(grid.ntheta()>>1))] += 
                //         - coeff2 * arr * x[grid.index(i_r,i_theta)] // Right
                //         + coeff2 * arr * x[grid.index(i_r,i_theta+(grid.ntheta()>>1))] // Center: (Right)
                //         - art * x[grid.index(i_r,i_theta+1)] / 4 // Top Right
                //         + art * x[grid.index(i_r,i_theta-1)] / 4; // Bottom Right

                //     result[grid.index(i_r+1,i_theta)] +=
                //         - coeff1 * arr * x[grid.index(i_r,i_theta)] // Left
                //         + coeff1 * arr * x[grid.index(i_r+1,i_theta)] // Center: (Left)
                //         + art * x[grid.index(i_r,i_theta+1)] / 4 // Top Left
                //         - art * x[grid.index(i_r,i_theta-1)] / 4; // Bottom Left

                //     result[grid.index(i_r,i_theta-1)] +=
                //         - coeff4 * att * x[grid.index(i_r,i_theta)] // Top
                //         + coeff4 * att * x[grid.index(i_r,i_theta-1)] // Center: (Top)
                //         - art * x[grid.index(i_r+1,i_theta)] / 4 // Top Right
                //         + art * x[grid.index(i_r,i_theta+(grid.ntheta()>>1))] / 4; // Top Left

                //     result[grid.index(i_r,i_theta+1)] +=
                //         - coeff3 * att * x[grid.index(i_r,i_theta)] // Bottom
                //         + coeff3 * att * x[grid.index(i_r,i_theta+1)] // Center: (Bottom)
                //         + art * x[grid.index(i_r+1,i_theta)] / 4 // Bottom Right
                //         - art * x[grid.index(i_r,i_theta+(grid.ntheta()>>1))] / 4; // Bottom Left    

                // } else if (i_r == grid.nr() - 2) {
                //     scalar_t h1 = grid.r_dist(i_r-1); 
                //     scalar_t h2 = grid.r_dist(i_r);
                //     scalar_t k1 = grid.theta_dist(i_theta-1);
                //     scalar_t k2 = grid.theta_dist(i_theta);

                //     double coeff1 = 0.5*(k1+k2)/h1;
                //     double coeff2 = 0.5*(k1+k2)/h2;
                //     double coeff3 = 0.5*(h1+h2)/k1;
                //     double coeff4 = 0.5*(h1+h2)/k2;
                
                //     result[grid.index(i_r,i_theta)] += 
                //         (h1+h2)*(k1+k2) * coeff_beta * fabs(detDFinv) / 4 // f_{i,j}
                //         - coeff1 * arr * x[grid.index(i_r-1,i_theta)] // Left
                //         - coeff2 * arr * x[grid.index(i_r+1,i_theta)] // Right
                //         - coeff3 * att * x[grid.index(i_r,i_theta-1)] // Bottom
                //         - coeff4 * att * x[grid.index(i_r,i_theta+1)] // Top
                //         + ((coeff1 + coeff2) * arr + (coeff3 + coeff4) * att) * x[grid.index(i_r,i_theta)]; // Center: (Left, Right, Bottom, Top)

                //     result[grid.index(i_r-1,i_theta)] += 
                //         - coeff2 * arr * x[grid.index(i_r,i_theta)] // Right
                //         + coeff2 * arr * x[grid.index(i_r-1,i_theta)] // Center: (Right)
                //         - art * x[grid.index(i_r,i_theta+1)] / 4 // Top Right
                //         + art * x[grid.index(i_r,i_theta-1)] / 4; // Bottom Right

                //     /* Don't write to the dirichlet boundary part! */
                //     // result[grid.index(i_r+1,i_theta)] +=
                //     //     - coeff1 * arr * x[grid.index(i_r,i_theta)] // Left
                //     //     + coeff1 * arr * x[grid.index(i_r+1,i_theta)] // Center: (Left)
                //     //     + art * x[grid.index(i_r,i_theta+1)] / 4 // Top Left
                //     //     - art * x[grid.index(i_r,i_theta-1)] / 4; // Bottom Left

                //     result[grid.index(i_r,i_theta-1)] +=
                //         - coeff4 * att * x[grid.index(i_r,i_theta)] // Top
                //         + coeff4 * att * x[grid.index(i_r,i_theta-1)] // Center: (Top)
                //         - art * x[grid.index(i_r+1,i_theta)] / 4 // Top Right
                //         + art * x[grid.index(i_r-1,i_theta)] / 4; // Top Left

                //     result[grid.index(i_r,i_theta+1)] +=
                //         - coeff3 * att * x[grid.index(i_r,i_theta)] // Bottom
                //         + coeff3 * att * x[grid.index(i_r,i_theta+1)] // Center: (Bottom)
                //         + art * x[grid.index(i_r+1,i_theta)] / 4 // Bottom Right
                //         - art * x[grid.index(i_r-1,i_theta)] / 4; // Bottom Left 
                // } else if (i_r == grid.nr() - 1) {
                //     // Dirichlet boundary
                //     result[grid.index(i_r,i_theta)] = x[grid.index(i_r,i_theta)];
                // }





































// void Operator::applyATake(const Level& onLevel, Vector<scalar_t>& result, const Vector<scalar_t>& x) const{
//     const PolarGrid& grid = onLevel.grid();
//     const ExactFunctions& exactFuncs = onLevel.exactFunctions();

//     assert(x.size() == grid.number_of_nodes());
//     assert(result.size() == grid.number_of_nodes());

//     const int numThreads = omp_get_max_threads();
//     omp_set_num_threads(numThreads);

//     #pragma omp parallel
//     {
//         double r, theta;
//         double arr, att, art;
//         double coeff_alpha, coeff_beta;
//         double Jrr, Jrt, Jtr, Jtt;
//         double detDFinv;
//         // double coeff, beta;
//         // Circular Smoother section
//         // For loop matches circular access pattern
//         #pragma omp for nowait
//         for (int i_r = 0; i_r < grid.numberSmootherCircles(); i_r++){
//             r = grid.radius(i_r);
//             coeff_alpha = exactFuncs.coeffs1(r, Rmax_);
//             coeff_beta = exactFuncs.coeffs2(r, Rmax_);
//             for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++){
//                 theta = grid.theta(i_theta);

//                 Jrr = exactFuncs.J_rr(r, theta, kappa_eps_, delta_e_, Rmax_);
//                 Jrt = exactFuncs.J_rt(r, theta, kappa_eps_, delta_e_, Rmax_);
//                 Jtr = exactFuncs.J_tr(r, theta, kappa_eps_, delta_e_, Rmax_);
//                 Jtt = exactFuncs.J_tt(r, theta, kappa_eps_, delta_e_, Rmax_);
//                 detDFinv = fabs(Jrr * Jtt - Jrt * Jtr);

//                 arr = 0.5 * (Jtt * Jtt + Jrt * Jrt) * coeff_alpha / detDFinv;
//                 art = -0.25 * (Jtt * Jtr + Jrt * Jrr) * coeff_alpha / detDFinv;
//                 att = 0.5 * (Jtr * Jtr + Jrr * Jrr) * coeff_alpha / detDFinv;

//                 // arr_att_art_beta(exactFuncs, r, theta, arr, att, art, coeff_beta);
//             }
//         }

//         // Radial smoother section
//         // For loop matches radial access pattern
//         #pragma omp for
//         for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++){
//             theta = grid.theta(i_theta);

//             for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++){
//                 r = grid.radius(i_r);

//                 // arr_att_art_beta(exactFuncs, r, theta, arr, att, art, coeff_beta);
//                 coeff_alpha = exactFuncs.coeffs1(r, Rmax_);
//                 coeff_beta = exactFuncs.coeffs2(r, Rmax_);

//                 Jrr = exactFuncs.J_rr(r, theta, kappa_eps_, delta_e_, Rmax_);
//                 Jrt = exactFuncs.J_rt(r, theta, kappa_eps_, delta_e_, Rmax_);
//                 Jtr = exactFuncs.J_tr(r, theta, kappa_eps_, delta_e_, Rmax_);
//                 Jtt = exactFuncs.J_tt(r, theta, kappa_eps_, delta_e_, Rmax_);
//                 detDFinv = fabs(Jrr * Jtt - Jrt * Jtr);

//                 arr = 0.5 * (Jtt * Jtt + Jrt * Jrt) * coeff_alpha / detDFinv;
//                 art = -0.25 * (Jtt * Jtr + Jrt * Jrr) * coeff_alpha / detDFinv;
//                 att = 0.5 * (Jtr * Jtr + Jrr * Jrr) * coeff_alpha / detDFinv;

//                 // arr_att_art(exactFuncs, r, theta, arr, att, art);
//             }
//         }
//     }
// }



// void Operator::applyATake(const Level& onLevel, Vector<scalar_t>& result, const Vector<scalar_t>& x) const{
//     const PolarGrid& grid = onLevel.grid();
//     const ExactFunctions& exactFuncs = onLevel.exactFunctions();

//     assert(x.size() == grid.number_of_nodes());
//     assert(result.size() == grid.number_of_nodes());

//     const int numThreads = omp_get_max_threads();
//     omp_set_num_threads(numThreads);

//     #pragma omp parallel
//     {
//         double r, theta;
//         double arr, att, art;
//         double coeff_alpha, coeff_beta;
//         double Jrr, Jrt, Jtr, Jtt;
//         double detDFinv;
//         // double coeff, beta;
//         // Circular Smoother section
//         // For loop matches circular access pattern
//         #pragma omp for nowait
//         for (int i_r = 0; i_r < grid.numberSmootherCircles(); i_r++){
//             r = grid.radius(i_r);
//             coeff_alpha = exactFuncs.coeffs1(r, Rmax_);
//             coeff_beta = exactFuncs.coeffs2(r, Rmax_);
//             for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++){
//                 theta = grid.theta(i_theta);

//                 Jrr = exactFuncs.J_rr(r, theta, kappa_eps_, delta_e_, Rmax_);
//                 Jrt = exactFuncs.J_rt(r, theta, kappa_eps_, delta_e_, Rmax_);
//                 Jtr = exactFuncs.J_tr(r, theta, kappa_eps_, delta_e_, Rmax_);
//                 Jtt = exactFuncs.J_tt(r, theta, kappa_eps_, delta_e_, Rmax_);
//                 detDFinv = fabs(Jrr * Jtt - Jrt * Jtr);

//                 arr = 0.5 * (Jtt * Jtt + Jrt * Jrt) * coeff_alpha / detDFinv;
//                 art = -0.25 * (Jtt * Jtr + Jrt * Jrr) * coeff_alpha / detDFinv;
//                 att = 0.5 * (Jtr * Jtr + Jrr * Jrr) * coeff_alpha / detDFinv;

//                 // arr_att_art_beta(exactFuncs, r, theta, arr, att, art, coeff_beta);
//             }
//         }

//         // Radial smoother section
//         // For loop matches radial access pattern
//         #pragma omp for
//         for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++){
//             theta = grid.theta(i_theta);

//             for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++){
//                 r = grid.radius(i_r);

//                 // arr_att_art_beta(exactFuncs, r, theta, arr, att, art, coeff_beta);
//                 coeff_alpha = exactFuncs.coeffs1(r, Rmax_);
//                 coeff_beta = exactFuncs.coeffs2(r, Rmax_);

//                 Jrr = exactFuncs.J_rr(r, theta, kappa_eps_, delta_e_, Rmax_);
//                 Jrt = exactFuncs.J_rt(r, theta, kappa_eps_, delta_e_, Rmax_);
//                 Jtr = exactFuncs.J_tr(r, theta, kappa_eps_, delta_e_, Rmax_);
//                 Jtt = exactFuncs.J_tt(r, theta, kappa_eps_, delta_e_, Rmax_);
//                 detDFinv = fabs(Jrr * Jtt - Jrt * Jtr);

//                 arr = 0.5 * (Jtt * Jtt + Jrt * Jrt) * coeff_alpha / detDFinv;
//                 art = -0.25 * (Jtt * Jtr + Jrt * Jrr) * coeff_alpha / detDFinv;
//                 att = 0.5 * (Jtr * Jtr + Jrr * Jrr) * coeff_alpha / detDFinv;

//                 // arr_att_art(exactFuncs, r, theta, arr, att, art);
//             }
//         }
//     }
// }




// double Operator::detDFinv(const ExactFunctions& exactFuncs, double r, double theta) const
// {
//     double Jrr = exactFuncs.J_rr(r, theta, kappa_eps_, delta_e_, Rmax_);
//     double Jrt = exactFuncs.J_rt(r, theta, kappa_eps_, delta_e_, Rmax_);
//     double Jtr = exactFuncs.J_tr(r, theta, kappa_eps_, delta_e_, Rmax_);
//     double Jtt = exactFuncs.J_tt(r, theta, kappa_eps_, delta_e_, Rmax_);
    
//     double detDFinv_r = Jrr * Jtt - Jrt * Jtr;

//     return detDFinv_r;
// }

// double Operator::arr(const ExactFunctions& exactFuncs, double r, double theta) const
// {
//     double Jrt = exactFuncs.J_rt(r, theta, kappa_eps_, delta_e_, Rmax_);
//     double Jtt = exactFuncs.J_tt(r, theta, kappa_eps_, delta_e_, Rmax_);
    
//     double detDFinv_r = detDFinv(exactFuncs, r, theta);
//     double alpha_coeff_r = exactFuncs.coeffs1(r, Rmax_);

//     double arr_r = 0.5 * (Jtt * Jtt + Jrt * Jrt) * alpha_coeff_r / fabs(detDFinv_r);

//     return arr_r;
// }

// double Operator::art(const ExactFunctions& exactFuncs, double r, double theta) const
// {
//     double Jrr = exactFuncs.J_rr(r, theta, kappa_eps_, delta_e_, Rmax_);
//     double Jrt = exactFuncs.J_rt(r, theta, kappa_eps_, delta_e_, Rmax_);
//     double Jtr = exactFuncs.J_tr(r, theta, kappa_eps_, delta_e_, Rmax_);
//     double Jtt = exactFuncs.J_tt(r, theta, kappa_eps_, delta_e_, Rmax_);

//     double detDFinv_r = detDFinv(exactFuncs, r, theta);
//     double alpha_coeff_r = exactFuncs.coeffs1(r, Rmax_);

//     double art_r = -0.25 * (Jtt * Jtr + Jrt * Jrr) * alpha_coeff_r / fabs(detDFinv_r);

//     return art_r;
// }

// double Operator::att(const ExactFunctions& exactFuncs, double r, double theta) const
// {
//     double Jrr = exactFuncs.J_rr(r, theta, kappa_eps_, delta_e_, Rmax_);
//     double Jtr = exactFuncs.J_tr(r, theta, kappa_eps_, delta_e_, Rmax_);

//     double detDFinv_r = detDFinv(exactFuncs, r, theta);
//     double alpha_coeff_r = exactFuncs.coeffs1(r, Rmax_);
    
//     double att_r = 0.5 * (Jtr * Jtr + Jrr * Jrr) * alpha_coeff_r / fabs(detDFinv_r);

//     return att_r;
// }

// inline void Operator::arr_att_art_beta(const ExactFunctions& exactFuncs, double r, double theta, double& arr, double& att, double& art, double& coeff_beta) const
// {
//     double Jrr = exactFuncs.J_rr(r, theta, kappa_eps_, delta_e_, Rmax_);
//     double Jrt = exactFuncs.J_rt(r, theta, kappa_eps_, delta_e_, Rmax_);
//     double Jtr = exactFuncs.J_tr(r, theta, kappa_eps_, delta_e_, Rmax_);
//     double Jtt = exactFuncs.J_tt(r, theta, kappa_eps_, delta_e_, Rmax_);
    
//     double detDFinv_r = Jrr * Jtt - Jrt * Jtr;

//     double coeff_alpha = exactFuncs.coeffs1(r, Rmax_);

//     arr = 0.5 * (Jtt * Jtt + Jrt * Jrt) * coeff_alpha / fabs(detDFinv_r);
//     art = -0.25 * (Jtt * Jtr + Jrt * Jrr) * coeff_alpha / fabs(detDFinv_r);
//     att = 0.5 * (Jtr * Jtr + Jrr * Jrr) * coeff_alpha / fabs(detDFinv_r);
//     coeff_beta = exactFuncs.coeffs2(r, Rmax_);
// }

// inline void Operator::arr_att_art(const ExactFunctions& exactFuncs, double r, double theta, double& arr, double& att, double& art) const
// {
//     double Jrr = exactFuncs.J_rr(r, theta, kappa_eps_, delta_e_, Rmax_);
//     double Jrt = exactFuncs.J_rt(r, theta, kappa_eps_, delta_e_, Rmax_);
//     double Jtr = exactFuncs.J_tr(r, theta, kappa_eps_, delta_e_, Rmax_);
//     double Jtt = exactFuncs.J_tt(r, theta, kappa_eps_, delta_e_, Rmax_);
    
//     double detDFinv_r = Jrr * Jtt - Jrt * Jtr;

//     double coeff_alpha = exactFuncs.coeffs1(r, Rmax_);

//     arr = 0.5 * (Jtt * Jtt + Jrt * Jrt) * coeff_alpha / fabs(detDFinv_r);
//     art = -0.25 * (Jtt * Jtr + Jrt * Jrr) * coeff_alpha / fabs(detDFinv_r);
//     att = 0.5 * (Jtr * Jtr + Jrr * Jrr) * coeff_alpha / fabs(detDFinv_r);
// }




























// void Operator::applyATake0(const Level& onLevel, Vector<scalar_t>& result, const Vector<scalar_t>& x) const{
//     const PolarGrid& grid = onLevel.grid();
//     const ExactFunctions& exactFuncs = onLevel.exactFunctions();

//     assert(x.size() == grid.number_of_nodes());
//     assert(result.size() == grid.number_of_nodes());

//     double arr, att, art, coeff_beta;

//     #pragma omp parallel for
//     for (int index = 0; index < grid.number_of_nodes(); index++) {
//         std::array<std::pair<scalar_t,scalar_t>, space_dimension> neighbor_distance;
//         MultiIndex node = grid.multiindex(index);
//         Point polar_coords = grid.polar_coordinates(node);

//         // arr_att_art_beta(exactFuncs, polar_coords[0], polar_coords[1], arr, att, art, coeff_beta);

//     }
// }



// void Operator::applyATake(const Level& onLevel, Vector<scalar_t>& result, const Vector<scalar_t>& x) const{
//     const PolarGrid& grid = onLevel.grid();
//     const ExactFunctions& exactFuncs = onLevel.exactFunctions();

//     assert(x.size() == grid.number_of_nodes());
//     assert(result.size() == grid.number_of_nodes());

//     const int numThreads = omp_get_max_threads();
//     omp_set_num_threads(numThreads);

//     #pragma omp parallel
//     {
//         double r, theta;
//         double arr, att, art;
//         double coeff_alpha, coeff_beta;
//         double Jrr, Jrt, Jtr, Jtt;
//         double detDFinv;
//         // double coeff, beta;
//         // Circular Smoother section
//         // For loop matches circular access pattern
//         #pragma omp for nowait
//         for (int i_r = 0; i_r < grid.numberSmootherCircles(); i_r++){
//             r = grid.radius(i_r);
//             coeff_alpha = exactFuncs.coeffs1(r, Rmax_);
//             coeff_beta = exactFuncs.coeffs2(r, Rmax_);
//             for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++){
//                 theta = grid.theta(i_theta);
//                 arr_att_art(exactFuncs, r, theta, i_theta, coeff_alpha, arr, art, att, detDFinv);

//             }
//         }

//         // Radial smoother section
//         // For loop matches radial access pattern
//         #pragma omp for
//         for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++){
//             theta = grid.theta(i_theta);
//             for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++){
//                 r = grid.radius(i_r);
//                 coeff_alpha = exactFuncs.coeffs1(r, Rmax_);
//                 coeff_beta = exactFuncs.coeffs2(r, Rmax_);
//                 arr_att_art(exactFuncs, r, theta, i_theta, coeff_alpha, arr, art, att, detDFinv);


//             }
//         }
//     }
// }