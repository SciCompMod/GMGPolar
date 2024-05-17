#include "../../include/Operator/operator.h"

void Operator::arr_att_art(const ExactFunctions& exactFuncs, double r, double theta, int i_theta, double coeff_alpha, double& arr, double& att, double& art, double& detDFinv) const{
    // With precomputed sinus and cosin values
    double Jrr = exactFuncs.J_rr(r, theta, kappa_eps_, delta_e_, Rmax_, sin_theta_[i_theta], cos_theta_[i_theta]);
    double Jrt = exactFuncs.J_rt(r, theta, kappa_eps_, delta_e_, Rmax_, sin_theta_[i_theta], cos_theta_[i_theta]);
    double Jtr = exactFuncs.J_tr(r, theta, kappa_eps_, delta_e_, Rmax_, sin_theta_[i_theta], cos_theta_[i_theta]);
    double Jtt = exactFuncs.J_tt(r, theta, kappa_eps_, delta_e_, Rmax_, sin_theta_[i_theta], cos_theta_[i_theta]);

    // Evaluate sinus and cosin values as needed
    // double Jrr = exactFuncs.J_rr(r, theta, kappa_eps_, delta_e_, Rmax_);
    // double Jrt = exactFuncs.J_rt(r, theta, kappa_eps_, delta_e_, Rmax_);
    // double Jtr = exactFuncs.J_tr(r, theta, kappa_eps_, delta_e_, Rmax_);
    // double Jtt = exactFuncs.J_tt(r, theta, kappa_eps_, delta_e_, Rmax_);
    
    detDFinv = Jrr * Jtt - Jrt * Jtr;
    arr = 0.5 * (Jtt * Jtt + Jrt * Jrt) * coeff_alpha / fabs(detDFinv);
    art = -0.25 * (Jtt * Jtr + Jrt * Jrr) * coeff_alpha / fabs(detDFinv);
    att = 0.5 * (Jtr * Jtr + Jrr * Jrr) * coeff_alpha / fabs(detDFinv);
}


void Operator::applyATake0(const Level& onLevel, Vector<scalar_t>& result, const Vector<scalar_t>& x) const{
    const PolarGrid& grid = onLevel.grid();
    const ExactFunctions& exactFuncs = onLevel.exactFunctions();

    assert(x.size() == grid.number_of_nodes());
    assert(result.size() == grid.number_of_nodes());

    double arr, arr_left, arr_right, arr_bottom, arr_top;
    double att, att_left, att_right, att_bottom, att_top;
    double art, art_left, art_right, art_bottom, art_top;
    double detDFinv, detDFinv_left, detDFinv_right, detDFinv_bottom, detDFinv_top;

    double coeff_alpha, coeff_beta;

    #pragma omp parallel for
    for(int index = 0; index < grid.number_of_nodes(); index ++){
        MultiIndex node = grid.multiindex(index);
        Point coords = grid.polar_coordinates(node);

        coeff_alpha = exactFuncs.coeffs1(coords[0], Rmax_);
        coeff_beta = exactFuncs.coeffs2(coords[0], Rmax_);

        arr_att_art(exactFuncs, coords[0], coords[1], node[1], coeff_alpha, arr, att, art, detDFinv);

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
                double coeff_alpha_left = exactFuncs.coeffs1(left_coords[0], Rmax_);
                arr_att_art(exactFuncs, left_coords[0], left_coords[1], left_node[1], coeff_alpha_left, arr_left, att_left, art_left, detDFinv_left);     
            }else{
                MultiIndex across_origin_node(0, (node[1] + grid.ntheta() / 2) % grid.ntheta());
                Point across_origin_coords = grid.polar_coordinates(across_origin_node); 

                std::cout<<across_origin_coords[0]<<", "<<across_origin_coords[1]<<std::endl;

                scalar_t coeff_alpha_left = exactFuncs.coeffs1(across_origin_coords[0], Rmax_);
                arr_att_art(exactFuncs, across_origin_coords[0], across_origin_coords[1], across_origin_node[1], coeff_alpha_left, arr_left, att_left, art_left, detDFinv_left);  
            }

            // Right
            if(neighbors[0].second != -1){
                MultiIndex right_node = grid.multiindex(neighbors[0].second);
                Point right_coords = grid.polar_coordinates(right_node); 
                scalar_t coeff_alpha_right = exactFuncs.coeffs1(right_coords[0], Rmax_);
                arr_att_art(exactFuncs, right_coords[0], right_coords[1], right_node[1], coeff_alpha_right, arr_right, att_right, art_right, detDFinv_right);     
            }
            // Bottom
            if(neighbors[1].first != -1){
                MultiIndex bottom_node = grid.multiindex(neighbors[1].first);
                Point bottom_coords = grid.polar_coordinates(bottom_node); 
                scalar_t coeff_alpha_bottom = exactFuncs.coeffs1(bottom_coords[0], Rmax_);
                arr_att_art(exactFuncs, bottom_coords[0], bottom_coords[1], bottom_node[1], coeff_alpha_bottom, arr_bottom, att_bottom, art_bottom, detDFinv_bottom);   
            }
            // Top
            if(neighbors[1].second != -1){
                MultiIndex top_node = grid.multiindex(neighbors[1].second);
                Point top_coords = grid.polar_coordinates(top_node); 
                scalar_t coeff_alpha_top = exactFuncs.coeffs1(top_coords[0], Rmax_);
                arr_att_art(exactFuncs, top_coords[0], top_coords[1], top_node[1], coeff_alpha_top, arr_top, att_top, art_top, detDFinv_top);   
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
        /* Don't give to the dirichlet boundary part! */ \
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
        /* Give values to the interior nodes! */ \
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
        double arr, att, art;
        double coeff_alpha, coeff_beta;
        double detDFinv;

        const int threadID = omp_get_thread_num();
        // ---------------------------------------------------------- //
        // Take care of the speration strips of the circular smoother //
        // ---------------------------------------------------------- //
        const int i_r_start = CircleSmootherTasks.getStart(threadID);
        const int i_r_end = CircleSmootherTasks.getEnd(threadID);
        const int i_r_separate = std::min(i_r_end - i_r_start, zone);

        // For loop matches circular access pattern
        for (int i_r = i_r_end - i_r_separate; i_r < i_r_end; i_r++){
            r = grid.radius(i_r);
            coeff_alpha = exactFuncs.coeffs1(r, Rmax_);
            coeff_beta = exactFuncs.coeffs2(r, Rmax_);
            for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++){
                theta = grid.theta(i_theta);
                arr_att_art(exactFuncs, r, theta, i_theta, coeff_alpha, arr, att, art, detDFinv);
                NODE_APPLY_A_GIVE(i_r, i_theta, grid, result, x, arr, att, art, coeff_beta, detDFinv);
            }
        }
        
        #pragma omp barrier

        // -------------------------------------------------------- //
        // Take care of the speration strips of the radial smoother //
        // -------------------------------------------------------- //
        const int i_theta_start = RadialSmootherTasks.getStart(threadID);
        const int i_theta_end = RadialSmootherTasks.getEnd(threadID);
        const int i_theta_seperate = std::min(i_theta_end-i_theta_start, zone);

        // For loop matches radial access pattern
        for (int i_theta = i_theta_start; i_theta < i_theta_start + i_theta_seperate; i_theta++){
            theta = grid.theta(i_theta);
            for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++){
                r = grid.radius(i_r);
                coeff_alpha = exactFuncs.coeffs1(r, Rmax_);
                coeff_beta = exactFuncs.coeffs2(r, Rmax_);
                arr_att_art(exactFuncs, r, theta, i_theta, coeff_alpha, arr, att, art, detDFinv);
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
            coeff_alpha = exactFuncs.coeffs1(r, Rmax_);
            coeff_beta = exactFuncs.coeffs2(r, Rmax_);
            for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++){
                theta = grid.theta(i_theta);
                arr_att_art(exactFuncs, r, theta, i_theta, coeff_alpha, arr, att, art, detDFinv);


                // if (i_r > 0 && i_r < grid.nr() - 2) {
                //     // scalar_t h1 = grid.r_dist(i_r-1); 
                //     // scalar_t h2 = grid.r_dist(i_r);
                //     // scalar_t k1 = grid.theta_dist(i_theta-1);
                //     // scalar_t k2 = grid.theta_dist(i_theta);

                //     // double coeff1 = 0.5*(k1+k2)/h1;
                //     // double coeff2 = 0.5*(k1+k2)/h2;
                //     // double coeff3 = 0.5*(h1+h2)/k1;
                //     // double coeff4 = 0.5*(h1+h2)/k2;

                //     // if(i_r == 4 && i_theta ==4){
                //     //     std::cout<<"START: "<<std::endl;
                //     //     std::cout<<" Beta: "<<0.25*(h1+h2)*(k1+k2) * coeff_beta * fabs(detDFinv) * x[grid.index(i_r,i_theta)]<<
                //     //         "\n Left "<< - coeff1 * arr * x[grid.index(i_r-1,i_theta)]<<", "<<- coeff2 * arr * x[grid.index(i_r+1,i_theta)]<<
                //     //         "\n Right "<< - coeff2 * arr * x[grid.index(i_r+1,i_theta)]<<
                //     //         "\n Bottom "<< - coeff3 * att * x[grid.index(i_r,i_theta-1)]<<
                //     //         "\n Top "<< - coeff4 * att * x[grid.index(i_r,i_theta+1)]<<
                //     //         "\n Center (all) " << ((coeff1 + 0) * arr + (0 + 0) * att) * x[grid.index(i_r,i_theta)]<<
                //     //         ", "<<((0 + coeff2) * arr + (0 + 0) * att) * x[grid.index(i_r,i_theta)]<<
                //     //         ", "<<((0 + 0) * arr + (coeff3 + 0) * att) * x[grid.index(i_r,i_theta)]<<
                //     //         ", "<<((0 + 0) * arr + (0 + coeff4) * att) * x[grid.index(i_r,i_theta)]<<  std::endl;
                //     // }
                
                //     // result[grid.index(i_r,i_theta)] += 
                //     //     0.25*(h1+h2)*(k1+k2) * coeff_beta * fabs(detDFinv) * x[grid.index(i_r,i_theta)] // beta_{i,j}
                //     //     - coeff1 * arr * x[grid.index(i_r-1,i_theta)] // Left
                //     //     - coeff2 * arr * x[grid.index(i_r+1,i_theta)] // Right
                //     //     - coeff3 * att * x[grid.index(i_r,i_theta-1)] // Bottom
                //     //     - coeff4 * att * x[grid.index(i_r,i_theta+1)] // Top
                //     //     + ((coeff1 + coeff2) * arr + (coeff3 + coeff4) * att) * x[grid.index(i_r,i_theta)]; // Center: (Left, Right, Bottom, Top)


                //     // if(i_r == 5 && i_theta ==4){
                //     //         std::cout<<"\n Right "<<  - coeff2 * arr * x[grid.index(i_r,i_theta)] <<
                //     //         "\n Center: (Right) "<< coeff2 * arr * x[grid.index(i_r-1,i_theta)] <<
                //     //         "\n Top Right "<<- 0.25 * art * x[grid.index(i_r,i_theta+1)] <<
                //     //         "\n Bottom Right " <<0.25 * art * x[grid.index(i_r,i_theta-1)]<<std::endl;
                //     // }

                //     // result[grid.index(i_r-1,i_theta)] += 
                //     //     - coeff2 * arr * x[grid.index(i_r,i_theta)] // Right
                //     //     + coeff2 * arr * x[grid.index(i_r-1,i_theta)] // Center: (Right)
                //     //     - 0.25 * art * x[grid.index(i_r,i_theta+1)] // Top Right
                //     //     + 0.25 * art * x[grid.index(i_r,i_theta-1)]; // Bottom Right

                //     // if(i_r == 3 && i_theta ==4){
                //     //         std::cout<<"\n Left "<< - coeff1 * arr * x[grid.index(i_r,i_theta)] <<
                //     //         "\n Center: (Left) "<<coeff1 * arr * x[grid.index(i_r+1,i_theta)]<<
                //     //         "\n Top Left "<< 0.25 * art * x[grid.index(i_r,i_theta+1)] <<
                //     //         "\n Bottom Left " << - 0.25 * art * x[grid.index(i_r,i_theta-1)]<<std::endl;
                //     // }

                //     // result[grid.index(i_r+1,i_theta)] +=
                //     //     - coeff1 * arr * x[grid.index(i_r,i_theta)] // Left
                //     //     + coeff1 * arr * x[grid.index(i_r+1,i_theta)] // Center: (Left)
                //     //     + 0.25 * art * x[grid.index(i_r,i_theta+1)] // Top Left
                //     //     - 0.25 * art * x[grid.index(i_r,i_theta-1)]; // Bottom Left

                //     // if(i_r == 4 && i_theta ==5){
                //     //         std::cout<<"\n Top "<<  - coeff4 * att * x[grid.index(i_r,i_theta)] <<
                //     //         "\n Center: (Top) "<< + coeff4 * att * x[grid.index(i_r,i_theta-1)] <<
                //     //         "\n Top Right "<<- 0.25 * art * x[grid.index(i_r+1,i_theta)] <<
                //     //         "\n Top Left " <<0.25 * art * x[grid.index(i_r-1,i_theta)]<<std::endl;
                //     // }

                //     // result[grid.index(i_r,i_theta-1)] +=
                //     //     - coeff4 * att * x[grid.index(i_r,i_theta)] // Top
                //     //     + coeff4 * att * x[grid.index(i_r,i_theta-1)] // Center: (Top)
                //     //     - 0.25 * art * x[grid.index(i_r+1,i_theta)] // Top Right
                //     //     + 0.25 * art * x[grid.index(i_r-1,i_theta)]; // Top Left

                //     // if(i_r == 4 && i_theta ==3){
                //     //         std::cout<<"\n Bottom "<<  - coeff3 * att * x[grid.index(i_r,i_theta)] <<
                //     //         "\n Center: (Bottom) "<< coeff3 * att * x[grid.index(i_r,i_theta+1)] <<
                //     //         "\n Bottom Right "<<0.25 * art * x[grid.index(i_r+1,i_theta)]  <<
                //     //         "\n Bottom Left " << - 0.25 * art * x[grid.index(i_r-1,i_theta)]<<std::endl;
                //     // }

                //     // result[grid.index(i_r,i_theta+1)] +=
                //     //     - coeff3 * att * x[grid.index(i_r,i_theta)] // Bottom
                //     //     + coeff3 * att * x[grid.index(i_r,i_theta+1)] // Center: (Bottom)
                //     //     + 0.25 * art * x[grid.index(i_r+1,i_theta)] // Bottom Right
                //     //     - 0.25 * art * x[grid.index(i_r-1,i_theta)]; // Bottom Left 

                       

                // } else if (i_r == 0) {
                //     // // h1 gets replaced with 2 * R0 and
                //     // // (i_r-1,i_theta) gets replaced with (i_r, i_theta + (grid.ntheta()>>1))
                //     // scalar_t h1 = 2 * grid.radius(0); 
                //     // scalar_t h2 = grid.r_dist(i_r);
                //     // scalar_t k1 = grid.theta_dist(i_theta-1);
                //     // scalar_t k2 = grid.theta_dist(i_theta);

                //     // double coeff1 = 0.5*(k1+k2)/h1;
                //     // double coeff2 = 0.5*(k1+k2)/h2;
                //     // double coeff3 = 0.5*(h1+h2)/k1;
                //     // double coeff4 = 0.5*(h1+h2)/k2;

                //     // result[grid.index(i_r,i_theta)] += 
                //     //     (h1+h2)*(k1+k2) * coeff_beta * fabs(detDFinv) / 4 // f_{i,j}
                //     //     - coeff1 * arr * x[grid.index(i_r,i_theta+(grid.ntheta()>>1))] // Left
                //     //     - coeff2 * arr * x[grid.index(i_r+1,i_theta)] // Right
                //     //     - coeff3 * att * x[grid.index(i_r,i_theta-1)] // Bottom
                //     //     - coeff4 * att * x[grid.index(i_r,i_theta+1)] // Top
                //     //     + ((coeff1 + coeff2) * arr + (coeff3 + coeff4) * att) * x[grid.index(i_r,i_theta)]; // Center: (Left, Right, Bottom, Top)

                //     // result[grid.index(i_r,i_theta+(grid.ntheta()>>1))] += 
                //     //     - coeff2 * arr * x[grid.index(i_r,i_theta)] // Right
                //     //     + coeff2 * arr * x[grid.index(i_r,i_theta+(grid.ntheta()>>1))] // Center: (Right)
                //     //     - art * x[grid.index(i_r,i_theta+1)] / 4 // Top Right
                //     //     + art * x[grid.index(i_r,i_theta-1)] / 4; // Bottom Right

                //     // result[grid.index(i_r+1,i_theta)] +=
                //     //     - coeff1 * arr * x[grid.index(i_r,i_theta)] // Left
                //     //     + coeff1 * arr * x[grid.index(i_r+1,i_theta)] // Center: (Left)
                //     //     + art * x[grid.index(i_r,i_theta+1)] / 4 // Top Left
                //     //     - art * x[grid.index(i_r,i_theta-1)] / 4; // Bottom Left

                //     // result[grid.index(i_r,i_theta-1)] +=
                //     //     - coeff4 * att * x[grid.index(i_r,i_theta)] // Top
                //     //     + coeff4 * att * x[grid.index(i_r,i_theta-1)] // Center: (Top)
                //     //     - art * x[grid.index(i_r+1,i_theta)] / 4 // Top Right
                //     //     + art * x[grid.index(i_r,i_theta+(grid.ntheta()>>1))] / 4; // Top Left

                //     // result[grid.index(i_r,i_theta+1)] +=
                //     //     - coeff3 * att * x[grid.index(i_r,i_theta)] // Bottom
                //     //     + coeff3 * att * x[grid.index(i_r,i_theta+1)] // Center: (Bottom)
                //     //     + art * x[grid.index(i_r+1,i_theta)] / 4 // Bottom Right
                //     //     - art * x[grid.index(i_r,i_theta+(grid.ntheta()>>1))] / 4; // Bottom Left    

                // } else if (i_r == grid.nr() - 2) {
                //     // scalar_t h1 = grid.r_dist(i_r-1); 
                //     // scalar_t h2 = grid.r_dist(i_r);
                //     // scalar_t k1 = grid.theta_dist(i_theta-1);
                //     // scalar_t k2 = grid.theta_dist(i_theta);

                //     // double coeff1 = 0.5*(k1+k2)/h1;
                //     // double coeff2 = 0.5*(k1+k2)/h2;
                //     // double coeff3 = 0.5*(h1+h2)/k1;
                //     // double coeff4 = 0.5*(h1+h2)/k2;
                
                //     // result[grid.index(i_r,i_theta)] += 
                //     //     (h1+h2)*(k1+k2) * coeff_beta * fabs(detDFinv) / 4 // f_{i,j}
                //     //     - coeff1 * arr * x[grid.index(i_r-1,i_theta)] // Left
                //     //     - coeff2 * arr * x[grid.index(i_r+1,i_theta)] // Right
                //     //     - coeff3 * att * x[grid.index(i_r,i_theta-1)] // Bottom
                //     //     - coeff4 * att * x[grid.index(i_r,i_theta+1)] // Top
                //     //     + ((coeff1 + coeff2) * arr + (coeff3 + coeff4) * att) * x[grid.index(i_r,i_theta)]; // Center: (Left, Right, Bottom, Top)

                //     // result[grid.index(i_r-1,i_theta)] += 
                //     //     - coeff2 * arr * x[grid.index(i_r,i_theta)] // Right
                //     //     + coeff2 * arr * x[grid.index(i_r-1,i_theta)] // Center: (Right)
                //     //     - art * x[grid.index(i_r,i_theta+1)] / 4 // Top Right
                //     //     + art * x[grid.index(i_r,i_theta-1)] / 4; // Bottom Right

                //     // /* Don't write to the dirichlet boundary part! */
                //     // // result[grid.index(i_r+1,i_theta)] +=
                //     // //     - coeff1 * arr * x[grid.index(i_r,i_theta)] // Left
                //     // //     + coeff1 * arr * x[grid.index(i_r+1,i_theta)] // Center: (Left)
                //     // //     + art * x[grid.index(i_r,i_theta+1)] / 4 // Top Left
                //     // //     - art * x[grid.index(i_r,i_theta-1)] / 4; // Bottom Left

                //     // result[grid.index(i_r,i_theta-1)] +=
                //     //     - coeff4 * att * x[grid.index(i_r,i_theta)] // Top
                //     //     + coeff4 * att * x[grid.index(i_r,i_theta-1)] // Center: (Top)
                //     //     - art * x[grid.index(i_r+1,i_theta)] / 4 // Top Right
                //     //     + art * x[grid.index(i_r-1,i_theta)] / 4; // Top Left

                //     // result[grid.index(i_r,i_theta+1)] +=
                //     //     - coeff3 * att * x[grid.index(i_r,i_theta)] // Bottom
                //     //     + coeff3 * att * x[grid.index(i_r,i_theta+1)] // Center: (Bottom)
                //     //     + art * x[grid.index(i_r+1,i_theta)] / 4 // Bottom Right
                //     //     - art * x[grid.index(i_r-1,i_theta)] / 4; // Bottom Left 
                // } else if (i_r == grid.nr() - 1) {
                //     // Dirichlet boundary
                //     // result[grid.index(i_r,i_theta)] = x[grid.index(i_r,i_theta)];
                // }







                NODE_APPLY_A_GIVE(i_r, i_theta, grid, result, x, arr, att, art, coeff_beta, detDFinv);
            }
        }

        // ---------------------------------------- //
        // Take care of the radial smoother section //
        // ---------------------------------------- //
        // For loop matches radial access pattern
        for (int i_theta = i_theta_start + i_theta_seperate; i_theta < i_theta_end; i_theta++){
            theta = grid.theta(i_theta);
            for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++){
                r = grid.radius(i_r);
                coeff_alpha = exactFuncs.coeffs1(r, Rmax_);
                coeff_beta = exactFuncs.coeffs2(r, Rmax_);
                arr_att_art(exactFuncs, r, theta, i_theta, coeff_alpha, arr, att, art, detDFinv);
                NODE_APPLY_A_GIVE(i_r, i_theta, grid, result, x, arr, att, art, coeff_beta, detDFinv);
            }
        }
    }
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

    #pragma omp parallel shared(dep)
    {
        double r, theta;
        double arr, att, art;
        double coeff_alpha, coeff_beta;
        double detDFinv;
        #pragma omp single
        {
            for(int i_r = S1 - 1; i_r >= 0; i_r -= 3) {
                #pragma omp task shared(dep) firstprivate(i_r) depend(out: dep[i_r])
                {
                    r = grid.radius(i_r);
                    coeff_alpha = exactFuncs.coeffs1(r, Rmax_);
                    coeff_beta = exactFuncs.coeffs2(r, Rmax_);
                    for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++){
                        theta = grid.theta(i_theta);
                        arr_att_art(exactFuncs, r, theta, i_theta, coeff_alpha, arr, att, art, detDFinv);
                        NODE_APPLY_A_GIVE(i_r, i_theta, grid, result, x, arr, att, art, coeff_beta, detDFinv);
                    }
                }
            }

            for(int i_r = S1 - 2; i_r >= 0; i_r -= 3) {
                #pragma omp task shared(dep) firstprivate(i_r) depend(in: dep[i_r-2], dep[i_r+1]) depend(out: dep[i_r])
                {
                    r = grid.radius(i_r);
                    coeff_alpha = exactFuncs.coeffs1(r, Rmax_);
                    coeff_beta = exactFuncs.coeffs2(r, Rmax_);
                    for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++){
                        theta = grid.theta(i_theta);
                        arr_att_art(exactFuncs, r, theta, i_theta, coeff_alpha, arr, att, art, detDFinv);
                        NODE_APPLY_A_GIVE(i_r, i_theta, grid, result, x, arr, att, art, coeff_beta, detDFinv);
                    }
                }
            }

            for(int i_r = S1 - 3; i_r >= 0; i_r -= 3) {
                #pragma omp task shared(dep) firstprivate(i_r) depend(in: dep[i_r-2], dep[i_r+1]) depend(out: dep[i_r])
                {
                    r = grid.radius(i_r);
                    coeff_alpha = exactFuncs.coeffs1(r, Rmax_);
                    coeff_beta = exactFuncs.coeffs2(r, Rmax_);
                    for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++){
                        theta = grid.theta(i_theta);
                        arr_att_art(exactFuncs, r, theta, i_theta, coeff_alpha, arr, att, art, detDFinv);
                        NODE_APPLY_A_GIVE(i_r, i_theta, grid, result, x, arr, att, art, coeff_beta, detDFinv);
                    }
                }
            }

            for(int i_theta = 0; i_theta < S2 - S2_wait; i_theta += 3) {
                #pragma omp task shared(dep) firstprivate(i_theta) depend(out: dep[S1+i_theta]) depend(in: dep[S2_start])
                {
                    theta = grid.theta(i_theta);
                    for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++){
                        r = grid.radius(i_r);
                        coeff_alpha = exactFuncs.coeffs1(r, Rmax_);
                        coeff_beta = exactFuncs.coeffs2(r, Rmax_);
                        arr_att_art(exactFuncs, r, theta, i_theta, coeff_alpha, arr, att, art, detDFinv);
                        NODE_APPLY_A_GIVE(i_r, i_theta, grid, result, x, arr, att, art, coeff_beta, detDFinv);
                    }
                }
            }

            for(int i_theta = 1; i_theta < S2 - S2_wait; i_theta += 3) {
                #pragma omp task firstprivate(i_theta) depend(out: dep[S1+i_theta]) depend(in: dep[S2_start], dep[S1+i_theta-1], dep[S1+i_theta+2])
                {
                    theta = grid.theta(i_theta);
                    for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++){
                        r = grid.radius(i_r);
                        coeff_alpha = exactFuncs.coeffs1(r, Rmax_);
                        coeff_beta = exactFuncs.coeffs2(r, Rmax_);
                        arr_att_art(exactFuncs, r, theta, i_theta, coeff_alpha, arr, att, art, detDFinv);
                        NODE_APPLY_A_GIVE(i_r, i_theta, grid, result, x, arr, att, art, coeff_beta, detDFinv);
                    }
                }
            }

            for(int i_theta = 2; i_theta < S2 - S2_wait; i_theta += 3) {
                #pragma omp task shared(dep) firstprivate(i_theta) depend(out: dep[S1+i_theta]) depend(in: dep[S2_start], dep[S1+i_theta-1], dep[S1+i_theta+2])
                {
                    theta = grid.theta(i_theta);
                    for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++){
                        r = grid.radius(i_r);
                        coeff_alpha = exactFuncs.coeffs1(r, Rmax_);
                        coeff_beta = exactFuncs.coeffs2(r, Rmax_);
                        arr_att_art(exactFuncs, r, theta, i_theta, coeff_alpha, arr, att, art, detDFinv);
                        NODE_APPLY_A_GIVE(i_r, i_theta, grid, result, x, arr, att, art, coeff_beta, detDFinv);
                    }
                }
            }

            if( S2_wait >= 1 ){
                int i_theta = S2 - S2_wait;
                #pragma omp task shared(dep) firstprivate(i_theta) depend(out: dep[S1+i_theta]) depend(in: dep[S2_start],  dep[S1+i_theta-1])
                {
                    theta = grid.theta(i_theta);
                    for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++){
                        r = grid.radius(i_r);
                        coeff_alpha = exactFuncs.coeffs1(r, Rmax_);
                        coeff_beta = exactFuncs.coeffs2(r, Rmax_);
                        arr_att_art(exactFuncs, r, theta, i_theta, coeff_alpha, arr, att, art, detDFinv);
                        NODE_APPLY_A_GIVE(i_r, i_theta, grid, result, x, arr, att, art, coeff_beta, detDFinv);
                    }
                }
            }

            if( S2_wait >= 2 ){
                int i_theta = S2 - S2_wait + 1;
                #pragma omp task shared(dep) firstprivate(i_theta) depend(out: dep[S1+i_theta]) depend(in: dep[S2_start],  dep[S1+i_theta-1])
                {
                    theta = grid.theta(i_theta);
                    for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++){
                        r = grid.radius(i_r);
                        coeff_alpha = exactFuncs.coeffs1(r, Rmax_);
                        coeff_beta = exactFuncs.coeffs2(r, Rmax_);
                        arr_att_art(exactFuncs, r, theta, i_theta, coeff_alpha, arr, att, art, detDFinv);
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