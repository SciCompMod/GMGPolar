#include "../../include/Interpolation/interpolation.h"

// RestrictionTake is faster than the RestrictionGive!

// For the restriction we use R = P^T.
// The restriction for an siotropic mesh reduces to
//           |1  2  1|
// R = 1/4 * |2  4  2| = P^T
//           |1  2  1|

void Interpolation::applyRestrictionTake0(const Level& fromLevel, const Level& toLevel, Vector<scalar_t>& result, const Vector<scalar_t>& x) const{
    assert(toLevel.level() == fromLevel.level() + 1);

    const PolarGrid& fineGrid = fromLevel.grid();
    const PolarGrid& coarseGrid = toLevel.grid();

    assert(x.size() == fineGrid.number_of_nodes());
    assert(result.size() == coarseGrid.number_of_nodes());

    #pragma omp parallel for
    for(int index = 0; index < coarseGrid.number_of_nodes(); index ++){
        MultiIndex coarse_node = coarseGrid.multiindex(index);
        MultiIndex fine_node(coarse_node[0]*2, coarse_node[1]*2);

        std::array<std::pair<scalar_t,scalar_t>, space_dimension> neighbor_distance;
        std::array<std::pair<int,int>, space_dimension> neighbors;

        // Center
        scalar_t value = x[fineGrid.index(fine_node)];

        fineGrid.adjacent_neighbors_of(fine_node, neighbors);

        // Left
        if(neighbors[0].first != -1){
            fineGrid.adjacent_neighbor_distances(fineGrid.multiindex(neighbors[0].first), neighbor_distance);
            value += neighbor_distance[0].second* x[neighbors[0].first] / (neighbor_distance[0].first + neighbor_distance[0].second);
        }

        // Right
        if(neighbors[0].second != -1){
            fineGrid.adjacent_neighbor_distances(fineGrid.multiindex(neighbors[0].second), neighbor_distance);
            value += neighbor_distance[0].first* x[neighbors[0].second] / (neighbor_distance[0].first + neighbor_distance[0].second);
        }
        
        // Bottom
        if(neighbors[1].first != -1){
            fineGrid.adjacent_neighbor_distances(fineGrid.multiindex(neighbors[1].first), neighbor_distance);
            value += neighbor_distance[1].second* x[neighbors[1].first] / (neighbor_distance[1].first + neighbor_distance[1].second);
        }
        
        // Top
        if(neighbors[1].second != -1){
            fineGrid.adjacent_neighbor_distances(fineGrid.multiindex(neighbors[1].second), neighbor_distance);
            value += neighbor_distance[1].first* x[neighbors[1].second] / (neighbor_distance[1].first + neighbor_distance[1].second);
        }

        fineGrid.diagonal_neighbors_of(fine_node, neighbors);

        // Bottom Left
        if(neighbors[0].first != -1){
            fineGrid.adjacent_neighbor_distances(fineGrid.multiindex(neighbors[0].first), neighbor_distance);
            value += neighbor_distance[0].second*neighbor_distance[1].second * x[neighbors[0].first] / ((neighbor_distance[0].first+neighbor_distance[0].second)*(neighbor_distance[1].first+neighbor_distance[1].second));
        }

        // Bottom Right
        if(neighbors[0].second != -1){
            fineGrid.adjacent_neighbor_distances(fineGrid.multiindex(neighbors[0].second), neighbor_distance);
            value += neighbor_distance[0].first*neighbor_distance[1].second* x[neighbors[0].second] / ((neighbor_distance[0].first + neighbor_distance[0].second)* (neighbor_distance[1].first + neighbor_distance[1].second));
        }
        
        // Top Left
        if(neighbors[1].first != -1){
            fineGrid.adjacent_neighbor_distances(fineGrid.multiindex(neighbors[1].first), neighbor_distance);
            value += neighbor_distance[0].second * neighbor_distance[1].first* x[neighbors[1].first] / ((neighbor_distance[0].first + neighbor_distance[0].second)* (neighbor_distance[1].first + neighbor_distance[1].second));
        }
        
        // Top Right
        if(neighbors[1].second != -1){
            fineGrid.adjacent_neighbor_distances(fineGrid.multiindex(neighbors[1].second), neighbor_distance);
            value +=  neighbor_distance[0].first * neighbor_distance[1].first* x[neighbors[1].second] / ((neighbor_distance[0].first + neighbor_distance[0].second)* (neighbor_distance[1].first + neighbor_distance[1].second));
        }

        result[index] = value;
    }
}

// ------------------------------------------ //
// Optimized version of applyRestrictionTake0 //
// ------------------------------------------ //

void Interpolation::applyRestrictionTake(const Level& fromLevel, const Level& toLevel, Vector<scalar_t>& result, const Vector<scalar_t>& x) const{
    assert(toLevel.level() == fromLevel.level() + 1);

    const PolarGrid& fineGrid = fromLevel.grid();
    const PolarGrid& coarseGrid = toLevel.grid();

    assert(x.size() == fineGrid.number_of_nodes());
    assert(result.size() == coarseGrid.number_of_nodes());

    // ---------------------- //
    // OpenMP Parallelization //
    // ---------------------- //
    const int numThreads = omp_get_max_threads();
    omp_set_num_threads(numThreads);
    const int minimalChunkSize = 1;

    // Distribute Tasks to each thread
    TaskDistribution CircleSmootherTasks(coarseGrid.numberSmootherCircles(), minimalChunkSize, numThreads);
    TaskDistribution RadialSmootherTasks(coarseGrid.ntheta(), minimalChunkSize, numThreads);

    #pragma omp parallel num_threads(numThreads)
    {   
        const int threadID = omp_get_thread_num();
        // ------------------------------------------ //
        // Take care of the circular Smoother section //
        // ------------------------------------------ //
        const int i_coarse_r_start = CircleSmootherTasks.getStart(threadID);
        const int i_coarse_r_end = CircleSmootherTasks.getEnd(threadID);

        const int innerSmootherCircle = 0;
        const int outerSmootherCircle = coarseGrid.numberSmootherCircles();

        // For loop matches circular access pattern
        for (int i_r_coarse = i_coarse_r_start; i_r_coarse < i_coarse_r_end; i_r_coarse++){
            int i_r = i_r_coarse << 1;
            for (int i_theta_coarse = 0; i_theta_coarse < coarseGrid.ntheta(); i_theta_coarse++){
                int i_theta = i_theta_coarse << 1;
                if(i_r_coarse > innerSmootherCircle && i_r_coarse < outerSmootherCircle - 1){
                    scalar_t h1 = fineGrid.r_dist(i_r-2);
                    scalar_t h2 = fineGrid.r_dist(i_r-1);
                    scalar_t h3 = fineGrid.r_dist(i_r);
                    scalar_t h4 = fineGrid.r_dist(i_r+1);

                    scalar_t k1 = fineGrid.theta_dist(i_theta-2);
                    scalar_t k2 = fineGrid.theta_dist(i_theta-1);
                    scalar_t k3 = fineGrid.theta_dist(i_theta);
                    scalar_t k4 = fineGrid.theta_dist(i_theta+1);

                    result[coarseGrid.index(i_r_coarse,i_theta_coarse)] =   
                        // Center
                        x[fineGrid.index(i_r,i_theta)] +
                        // Left, Right, Bottom, Top
                        h2 * x[fineGrid.index(i_r-1, i_theta)] / (h1+h2) +
                        h3 * x[fineGrid.index(i_r+1, i_theta)] / (h3+h4) +
                        k2 * x[fineGrid.index(i_r, i_theta-1)] / (k1+k2) +
                        k3 * x[fineGrid.index(i_r, i_theta+1)] / (k3+k4) +
                        // Bottom Left, Bottom Right, Top Left, Top Right
                        h2*k2 * x[fineGrid.index(i_r-1, i_theta-1)] / ((h1+h2)*(k1+k2)) +
                        h3*k2 * x[fineGrid.index(i_r+1, i_theta-1)] / ((h3+h4)*(k1+k2)) +
                        h2*k3 * x[fineGrid.index(i_r-1, i_theta+1)] / ((h1+h2)*(k3+k4)) +
                        h3*k3 * x[fineGrid.index(i_r+1, i_theta+1)] / ((h3+h4)*(k3+k4));
                } else{
                    /* First and Last Circle have to be checked for domain boundary */
                    // Middle Part
                    scalar_t k1 = fineGrid.theta_dist(i_theta-2);
                    scalar_t k2 = fineGrid.theta_dist(i_theta-1);
                    scalar_t k3 = fineGrid.theta_dist(i_theta);
                    scalar_t k4 = fineGrid.theta_dist(i_theta+1);
                    // Center, Bottom, Top
                    scalar_t value = x[fineGrid.index(i_r,i_theta)] +
                        k2 * x[fineGrid.index(i_r, i_theta-1)] / (k1+k2) +
                        k3 * x[fineGrid.index(i_r, i_theta+1)] / (k3+k4);

                    if(i_r_coarse > 0){
                        // Left Part
                        scalar_t h1 = fineGrid.r_dist(i_r-2);
                        scalar_t h2 = fineGrid.r_dist(i_r-1);
                        // Left, Bottom Left, Bottom Right
                        value += h2 * x[fineGrid.index(i_r-1, i_theta)] / (h1+h2) +
                            h2*k2 * x[fineGrid.index(i_r-1, i_theta-1)] / ((h1+h2)*(k1+k2)) +
                            h2*k3 * x[fineGrid.index(i_r-1, i_theta+1)] / ((h1+h2)*(k3+k4));                       
                    } 
                    if(i_r_coarse < coarseGrid.nr() - 1){
                       // Right Part
                        scalar_t h3 = fineGrid.r_dist(i_r);
                        scalar_t h4 = fineGrid.r_dist(i_r+1);
                        // Right, Bottom Right, Top Right
                        value += h3 * x[fineGrid.index(i_r+1, i_theta)] / (h3+h4) +
                            h3*k2 * x[fineGrid.index(i_r+1, i_theta-1)] / ((h3+h4)*(k1+k2)) +
                            h3*k3 * x[fineGrid.index(i_r+1, i_theta+1)] / ((h3+h4)*(k3+k4));
                    }
                    result[coarseGrid.index(i_r_coarse,i_theta_coarse)] = value;
                }
            }
        }

        // ---------------------------------------- //
        // Take care of the radial smoother section //
        // ---------------------------------------- //
        const int i_coarse_theta_start = RadialSmootherTasks.getStart(threadID);
        const int i_coarse_theta_end = RadialSmootherTasks.getEnd(threadID);

        const int i_r_coarse_start = coarseGrid.numberSmootherCircles();
        const int i_r_coarse_end = coarseGrid.nr();

        const int firstRadialSmootherNodes = coarseGrid.numberSmootherCircles();
        const int lastRadialSmootherNodes = coarseGrid.nr();

        // For loop matches circular access pattern
        for (int i_theta_coarse = i_coarse_theta_start; i_theta_coarse < i_coarse_theta_end; i_theta_coarse++){
            int i_theta = i_theta_coarse << 1;
            for (int i_r_coarse = i_r_coarse_start; i_r_coarse < i_r_coarse_end; i_r_coarse++){
                int i_r = i_r_coarse << 1;
                if(i_r_coarse > firstRadialSmootherNodes && i_r_coarse < lastRadialSmootherNodes - 1){
                    scalar_t h1 = fineGrid.r_dist(i_r-2);
                    scalar_t h2 = fineGrid.r_dist(i_r-1);
                    scalar_t h3 = fineGrid.r_dist(i_r);
                    scalar_t h4 = fineGrid.r_dist(i_r+1);

                    scalar_t k1 = fineGrid.theta_dist(i_theta-2);
                    scalar_t k2 = fineGrid.theta_dist(i_theta-1);
                    scalar_t k3 = fineGrid.theta_dist(i_theta);
                    scalar_t k4 = fineGrid.theta_dist(i_theta+1);

                    result[coarseGrid.index(i_r_coarse,i_theta_coarse)] =   
                        // Center
                        x[fineGrid.index(i_r,i_theta)] +
                        // Left, Right, Bottom, Top
                        h2 * x[fineGrid.index(i_r-1, i_theta)] / (h1+h2) +
                        h3 * x[fineGrid.index(i_r+1, i_theta)] / (h3+h4) +
                        k2 * x[fineGrid.index(i_r, i_theta-1)] / (k1+k2) +
                        k3 * x[fineGrid.index(i_r, i_theta+1)] / (k3+k4) +
                        // Bottom Left, Bottom Right, Top Left, Top Right
                        h2*k2 * x[fineGrid.index(i_r-1, i_theta-1)] / ((h1+h2)*(k1+k2)) +
                        h3*k2 * x[fineGrid.index(i_r+1, i_theta-1)] / ((h3+h4)*(k1+k2)) +
                        h2*k3 * x[fineGrid.index(i_r-1, i_theta+1)] / ((h1+h2)*(k3+k4)) +
                        h3*k3 * x[fineGrid.index(i_r+1, i_theta+1)] / ((h3+h4)*(k3+k4));
                } else{
                    /* First and Last radial node have to be checked for domain boundary */
                    // Middle Part
                    scalar_t k1 = fineGrid.theta_dist(i_theta-2);
                    scalar_t k2 = fineGrid.theta_dist(i_theta-1);
                    scalar_t k3 = fineGrid.theta_dist(i_theta);
                    scalar_t k4 = fineGrid.theta_dist(i_theta+1);
                    // Center, Bottom, Top
                    scalar_t value = x[fineGrid.index(i_r,i_theta)] +
                        k2 * x[fineGrid.index(i_r, i_theta-1)] / (k1+k2) +
                        k3 * x[fineGrid.index(i_r, i_theta+1)] / (k3+k4);
                    if(i_r_coarse > 0){
                        // Left Part
                        scalar_t h1 = fineGrid.r_dist(i_r-2);
                        scalar_t h2 = fineGrid.r_dist(i_r-1);
                        // Left, Bottom Left, Bottom Right
                        value += h2 * x[fineGrid.index(i_r-1, i_theta)] / (h1+h2) +
                            h2*k2 * x[fineGrid.index(i_r-1, i_theta-1)] / ((h1+h2)*(k1+k2)) +
                            h2*k3 * x[fineGrid.index(i_r-1, i_theta+1)] / ((h1+h2)*(k3+k4));               
                    } 
                    if(i_r_coarse < coarseGrid.nr() - 1){
                       // Right Part
                        scalar_t h3 = fineGrid.r_dist(i_r);
                        scalar_t h4 = fineGrid.r_dist(i_r+1);
                        // Right, Bottom Right, Top Right
                        value += h3 * x[fineGrid.index(i_r+1, i_theta)] / (h3+h4) +
                            h3*k2 * x[fineGrid.index(i_r+1, i_theta-1)] / ((h3+h4)*(k1+k2)) +
                            h3*k3 * x[fineGrid.index(i_r+1, i_theta+1)] / ((h3+h4)*(k3+k4));
                    }
                    result[coarseGrid.index(i_r_coarse,i_theta_coarse)] = value;
                }
            }
        }
    }
}


// ----------------------------------------------- //
// Task dependency version of applyRestrictionTake //
// ----------------------------------------------- //

void Interpolation::applyRestrictionTakeTasks(const Level& fromLevel, const Level& toLevel, Vector<scalar_t>& result, const Vector<scalar_t>& x) const{
    assert(toLevel.level() == fromLevel.level() + 1);

    const PolarGrid& fineGrid = fromLevel.grid();
    const PolarGrid& coarseGrid = toLevel.grid();

    assert(x.size() == fineGrid.number_of_nodes());
    assert(result.size() == coarseGrid.number_of_nodes());

    // ---------------------- //
    // OpenMP Parallelization //
    // ---------------------- //
    const int numThreads = omp_get_max_threads();
    omp_set_num_threads(numThreads);

    #pragma omp parallel num_threads(numThreads)
    {   
        #pragma omp single
        {
            // ------------------------------------------ //
            // Take care of the circular Smoother section //
            // ------------------------------------------ //
            const int innerSmootherCircle = 0;
            const int outerSmootherCircle = coarseGrid.numberSmootherCircles();
            for (int i_r_coarse = 0; i_r_coarse < coarseGrid.numberSmootherCircles(); i_r_coarse++){
                #pragma omp task firstprivate(i_r_coarse)
                {
                    int i_r = i_r_coarse << 1;
                    for (int i_theta_coarse = 0; i_theta_coarse < coarseGrid.ntheta(); i_theta_coarse++){
                        int i_theta = i_theta_coarse << 1;
                        if(i_r_coarse > innerSmootherCircle && i_r_coarse < outerSmootherCircle - 1){
                            scalar_t h1 = fineGrid.r_dist(i_r-2);
                            scalar_t h2 = fineGrid.r_dist(i_r-1);
                            scalar_t h3 = fineGrid.r_dist(i_r);
                            scalar_t h4 = fineGrid.r_dist(i_r+1);

                            scalar_t k1 = fineGrid.theta_dist(i_theta-2);
                            scalar_t k2 = fineGrid.theta_dist(i_theta-1);
                            scalar_t k3 = fineGrid.theta_dist(i_theta);
                            scalar_t k4 = fineGrid.theta_dist(i_theta+1);

                            result[coarseGrid.index(i_r_coarse,i_theta_coarse)] =   
                                // Center
                                x[fineGrid.index(i_r,i_theta)] +
                                // Left, Right, Bottom, Top
                                h2 * x[fineGrid.index(i_r-1, i_theta)] / (h1+h2) +
                                h3 * x[fineGrid.index(i_r+1, i_theta)] / (h3+h4) +
                                k2 * x[fineGrid.index(i_r, i_theta-1)] / (k1+k2) +
                                k3 * x[fineGrid.index(i_r, i_theta+1)] / (k3+k4) +
                                // Bottom Left, Bottom Right, Top Left, Top Right
                                h2*k2 * x[fineGrid.index(i_r-1, i_theta-1)] / ((h1+h2)*(k1+k2)) +
                                h3*k2 * x[fineGrid.index(i_r+1, i_theta-1)] / ((h3+h4)*(k1+k2)) +
                                h2*k3 * x[fineGrid.index(i_r-1, i_theta+1)] / ((h1+h2)*(k3+k4)) +
                                h3*k3 * x[fineGrid.index(i_r+1, i_theta+1)] / ((h3+h4)*(k3+k4));
                        } else{
                            /* First and Last Circle have to be checked for domain boundary */
                            // Middle Part
                            scalar_t k1 = fineGrid.theta_dist(i_theta-2);
                            scalar_t k2 = fineGrid.theta_dist(i_theta-1);
                            scalar_t k3 = fineGrid.theta_dist(i_theta);
                            scalar_t k4 = fineGrid.theta_dist(i_theta+1);
                            // Center, Bottom, Top
                            scalar_t value = x[fineGrid.index(i_r,i_theta)] +
                                k2 * x[fineGrid.index(i_r, i_theta-1)] / (k1+k2) +
                                k3 * x[fineGrid.index(i_r, i_theta+1)] / (k3+k4);

                            if(i_r_coarse > 0){
                                // Left Part
                                scalar_t h1 = fineGrid.r_dist(i_r-2);
                                scalar_t h2 = fineGrid.r_dist(i_r-1);
                                // Left, Bottom Left, Bottom Right
                                value += h2 * x[fineGrid.index(i_r-1, i_theta)] / (h1+h2) +
                                    h2*k2 * x[fineGrid.index(i_r-1, i_theta-1)] / ((h1+h2)*(k1+k2)) +
                                    h2*k3 * x[fineGrid.index(i_r-1, i_theta+1)] / ((h1+h2)*(k3+k4));                       
                            } 
                            if(i_r_coarse < coarseGrid.nr() - 1){
                            // Right Part
                                scalar_t h3 = fineGrid.r_dist(i_r);
                                scalar_t h4 = fineGrid.r_dist(i_r+1);
                                // Right, Bottom Right, Top Right
                                value += h3 * x[fineGrid.index(i_r+1, i_theta)] / (h3+h4) +
                                    h3*k2 * x[fineGrid.index(i_r+1, i_theta-1)] / ((h3+h4)*(k1+k2)) +
                                    h3*k3 * x[fineGrid.index(i_r+1, i_theta+1)] / ((h3+h4)*(k3+k4));
                            }
                            result[coarseGrid.index(i_r_coarse,i_theta_coarse)] = value;
                        }
                    }
                }
            }

            // ---------------------------------------- //
            // Take care of the radial smoother section //
            // ---------------------------------------- //
            const int firstRadialSmootherNodes = coarseGrid.numberSmootherCircles();
            const int lastRadialSmootherNodes = coarseGrid.nr();
            // For loop matches circular access pattern
            for (int i_theta_coarse = 0; i_theta_coarse < coarseGrid.ntheta(); i_theta_coarse++){
                #pragma omp task firstprivate(i_theta_coarse)
                {
                    int i_theta = i_theta_coarse << 1;
                    for (int i_r_coarse = coarseGrid.numberSmootherCircles(); i_r_coarse < coarseGrid.nr(); i_r_coarse++){
                        int i_r = i_r_coarse << 1;
                        if(i_r_coarse > firstRadialSmootherNodes && i_r_coarse < lastRadialSmootherNodes - 1){
                            scalar_t h1 = fineGrid.r_dist(i_r-2);
                            scalar_t h2 = fineGrid.r_dist(i_r-1);
                            scalar_t h3 = fineGrid.r_dist(i_r);
                            scalar_t h4 = fineGrid.r_dist(i_r+1);

                            scalar_t k1 = fineGrid.theta_dist(i_theta-2);
                            scalar_t k2 = fineGrid.theta_dist(i_theta-1);
                            scalar_t k3 = fineGrid.theta_dist(i_theta);
                            scalar_t k4 = fineGrid.theta_dist(i_theta+1);

                            result[coarseGrid.index(i_r_coarse,i_theta_coarse)] =   
                                // Center
                                x[fineGrid.index(i_r,i_theta)] +
                                // Left, Right, Bottom, Top
                                h2 * x[fineGrid.index(i_r-1, i_theta)] / (h1+h2) +
                                h3 * x[fineGrid.index(i_r+1, i_theta)] / (h3+h4) +
                                k2 * x[fineGrid.index(i_r, i_theta-1)] / (k1+k2) +
                                k3 * x[fineGrid.index(i_r, i_theta+1)] / (k3+k4) +
                                // Bottom Left, Bottom Right, Top Left, Top Right
                                h2*k2 * x[fineGrid.index(i_r-1, i_theta-1)] / ((h1+h2)*(k1+k2)) +
                                h3*k2 * x[fineGrid.index(i_r+1, i_theta-1)] / ((h3+h4)*(k1+k2)) +
                                h2*k3 * x[fineGrid.index(i_r-1, i_theta+1)] / ((h1+h2)*(k3+k4)) +
                                h3*k3 * x[fineGrid.index(i_r+1, i_theta+1)] / ((h3+h4)*(k3+k4));
                        } else{
                            /* First and Last radial node have to be checked for domain boundary */
                            // Middle Part
                            scalar_t k1 = fineGrid.theta_dist(i_theta-2);
                            scalar_t k2 = fineGrid.theta_dist(i_theta-1);
                            scalar_t k3 = fineGrid.theta_dist(i_theta);
                            scalar_t k4 = fineGrid.theta_dist(i_theta+1);
                            // Center, Bottom, Top
                            scalar_t value = x[fineGrid.index(i_r,i_theta)] +
                                k2 * x[fineGrid.index(i_r, i_theta-1)] / (k1+k2) +
                                k3 * x[fineGrid.index(i_r, i_theta+1)] / (k3+k4);
                            if(i_r_coarse > 0){
                                // Left Part
                                scalar_t h1 = fineGrid.r_dist(i_r-2);
                                scalar_t h2 = fineGrid.r_dist(i_r-1);
                                // Left, Bottom Left, Bottom Right
                                value += h2 * x[fineGrid.index(i_r-1, i_theta)] / (h1+h2) +
                                    h2*k2 * x[fineGrid.index(i_r-1, i_theta-1)] / ((h1+h2)*(k1+k2)) +
                                    h2*k3 * x[fineGrid.index(i_r-1, i_theta+1)] / ((h1+h2)*(k3+k4));               
                            } 
                            if(i_r_coarse < coarseGrid.nr() - 1){
                            // Right Part
                                scalar_t h3 = fineGrid.r_dist(i_r);
                                scalar_t h4 = fineGrid.r_dist(i_r+1);
                                // Right, Bottom Right, Top Right
                                value += h3 * x[fineGrid.index(i_r+1, i_theta)] / (h3+h4) +
                                    h3*k2 * x[fineGrid.index(i_r+1, i_theta-1)] / ((h3+h4)*(k1+k2)) +
                                    h3*k3 * x[fineGrid.index(i_r+1, i_theta+1)] / ((h3+h4)*(k3+k4));
                            }
                            result[coarseGrid.index(i_r_coarse,i_theta_coarse)] = value;
                        }
                    }
                }
            }
        }
    }
}