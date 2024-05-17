#include "../../include/Interpolation/interpolation.h"

// RestrictionTake is faster than the RestrictionGive!

// For the restriction we use R = P^T.
// The restriction for an siotropic mesh reduces to
//           |1  2  1|
// R = 1/4 * |2  4  2| = P^T
//           |1  2  1|

void Interpolation::applyRestrictionGive0(const Level& fromLevel, const Level& toLevel, Vector<scalar_t>& result, const Vector<scalar_t>& x) const{
    assert(toLevel.level() == fromLevel.level() + 1);

    const PolarGrid& fineGrid = fromLevel.grid();
    const PolarGrid& coarseGrid = toLevel.grid();

    assert(x.size() == fineGrid.number_of_nodes());
    assert(result.size() == coarseGrid.number_of_nodes());
    
    #pragma omp declare reduction(vec_float_plus : Vector<scalar_t> : \
                                std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<scalar_t>())) \
                        initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))

    assign(result, 0.0);

    #pragma omp parallel for reduction(vec_float_plus : result)
    for(int index = 0; index < fineGrid.number_of_nodes(); index ++){
        std::array<std::pair<scalar_t,scalar_t>, space_dimension> neighbor_distance;
        std::array<std::pair<int,int>, space_dimension> neighbors;

        MultiIndex fine_node = fineGrid.multiindex(index);

        // Fine node appears in coarse grid
        if(fine_node[0] % 2 == 0 && fine_node[1] % 2 == 0){
            // Input x needs a fine grid index: x[FINE_INDEX]
            // Result needs a coarse grid index: result[COARSE_INDEX]
            MultiIndex coarse_node(fine_node[0] / 2, fine_node[1] / 2);
            result[coarseGrid.index(coarse_node)] += x[index];
        }

        // Fine node between coarse nodes in theta direction
        if(fine_node[0] % 2 == 0 && fine_node[1] % 2 == 1){
            fineGrid.adjacent_neighbor_distances(fine_node, neighbor_distance);
            scalar_t k1 = neighbor_distance[1].first;
            scalar_t k2 = neighbor_distance[1].second;

            fineGrid.adjacent_neighbors_of(fine_node, neighbors);

            MultiIndex bottom_coarse_node(fine_node[0] / 2, fine_node[1] / 2);
            MultiIndex top_coarse_node(fine_node[0] / 2, (fine_node[1] / 2 + 1) % coarseGrid.ntheta());

            result[coarseGrid.index(bottom_coarse_node)] += k1 *x[index] / (k1 + k2);
            result[coarseGrid.index(top_coarse_node)] += k2 *x[index]/ (k1+k2);
        }


        // Fine node between coarse nodes in radial direction
        if(fine_node[0] % 2 == 1 && fine_node[1] % 2 == 0){
            fineGrid.adjacent_neighbor_distances(fine_node, neighbor_distance);
            scalar_t h1 = neighbor_distance[0].first;
            scalar_t h2 = neighbor_distance[0].second;

            fineGrid.adjacent_neighbors_of(fine_node, neighbors);

            MultiIndex left_coarse_node(fine_node[0] / 2, fine_node[1] / 2);
            MultiIndex right_coarse_node(fine_node[0] / 2 + 1, fine_node[1] / 2);

            result[coarseGrid.index(left_coarse_node)] +=  (h1 *x[index] ) / (h1+h2);;
            result[coarseGrid.index(right_coarse_node)] +=  (h2 *x[index] ) / (h1+h2);;
        }


        //Fine node in the center of four coarse nodes
        if(fine_node[0] % 2 == 1 && fine_node[1] % 2 == 1){

            fineGrid.adjacent_neighbor_distances(fine_node, neighbor_distance);
            scalar_t h1 = neighbor_distance[0].first;
            scalar_t h2 = neighbor_distance[0].second;
            scalar_t k1 = neighbor_distance[1].first;
            scalar_t k2 = neighbor_distance[1].second;

            fineGrid.adjacent_neighbors_of(fine_node, neighbors);

            MultiIndex bottom_left_coarse_node(fine_node[0] / 2, fine_node[1] / 2);
            MultiIndex bottom_right_coarse_node(fine_node[0] / 2 + 1, fine_node[1] / 2);
            MultiIndex top_left_coarse_node(fine_node[0] / 2, (fine_node[1] / 2 + 1) % coarseGrid.ntheta());
            MultiIndex top_right_node(fine_node[0] / 2 + 1, (fine_node[1] / 2 + 1) % coarseGrid.ntheta());

            result[coarseGrid.index(bottom_left_coarse_node)] += h1*k1 *x[index]  / ((h1 + h2) * (k1 + k2));
            result[coarseGrid.index(bottom_right_coarse_node)] += h2*k1*x[index]  / ((h1 + h2) * (k1 + k2));
            result[coarseGrid.index(top_left_coarse_node)] += h1*k2 *x[index] / ((h1 + h2) * (k1 + k2));
            result[coarseGrid.index(top_right_node)] += h2*k2*x[index] / ((h1 + h2) * (k1 + k2)) ;
        }
    }
}

// ------------------------------------------ //
// Optimized version of applyRestrictionGive0 //
// ------------------------------------------ //

#define FINE_NODE_GIVE_RESTRICTION(i_r, i_theta, i_r_coarse, i_theta_coarse, fineGrid, coarseGrid, result, x) \
do { \
    if(i_r & 1) { \
        if(i_theta & 1){ \
            /* Fine node in the center of four coarse nodes */ \
            scalar_t h1 = fineGrid.r_dist(i_r-1); \
            scalar_t h2 = fineGrid.r_dist(i_r); \
            scalar_t k1 = fineGrid.theta_dist(i_theta-1); \
            scalar_t k2 = fineGrid.theta_dist(i_theta); \
            scalar_t divisor = (h1+h2) * (k1+k2); \
            scalar_t value = x[fineGrid.index(i_r, i_theta)]; \
            /* Bottom Left, Bottom Right, Top Left, Top Right */ \
            result[coarseGrid.index(i_r_coarse,   i_theta_coarse)] += h1*k1 * value / divisor; \
            result[coarseGrid.index(i_r_coarse+1, i_theta_coarse)] += h2*k1 * value / divisor; \
            result[coarseGrid.index(i_r_coarse,   i_theta_coarse+1)] += h1*k2 * value / divisor; \
            result[coarseGrid.index(i_r_coarse+1, i_theta_coarse+1)] += h2*k2 * value / divisor; \
        } \
        else{ \
            /* Fine node between coarse nodes in radial direction */ \
            scalar_t h1 = fineGrid.r_dist(i_r-1); \
            scalar_t h2 = fineGrid.r_dist(i_r); \
            scalar_t divisor = (h1+h2); \
            scalar_t value = x[fineGrid.index(i_r, i_theta)]; \
            /* Left; Right */ \
            result[coarseGrid.index(i_r_coarse,   i_theta_coarse)] += h1 * value / divisor; \
            result[coarseGrid.index(i_r_coarse+1, i_theta_coarse)] += h2 * value / divisor; \
        } \
    } \
    else{ \
        if(i_theta & 1){ \
            /* Fine node between coarse nodes in theta direction */ \
            scalar_t k1 = fineGrid.theta_dist(i_theta-1); \
            scalar_t k2 = fineGrid.theta_dist(i_theta); \
            scalar_t divisor = (k1+k2); \
            scalar_t value = x[fineGrid.index(i_r, i_theta)]; \
            /* Bottom, Top */ \
            result[coarseGrid.index(i_r_coarse, i_theta_coarse)] += k1 * value / divisor; \
            result[coarseGrid.index(i_r_coarse, i_theta_coarse+1)] += k2 * value / divisor; \
        } \
        else{ \
            /* Fine node appears in coarse grid */ \
            result[coarseGrid.index(i_r_coarse, i_theta_coarse)] += x[fineGrid.index(i_r, i_theta)]; \
        } \
    } \
} while(0)

void Interpolation::applyRestrictionGive(const Level& fromLevel, const Level& toLevel, Vector<scalar_t>& result, const Vector<scalar_t>& x) const{
    assert(toLevel.level() == fromLevel.level() + 1);

    const PolarGrid& fineGrid = fromLevel.grid();
    const PolarGrid& coarseGrid = toLevel.grid();

    assert(x.size() == fineGrid.number_of_nodes());
    assert(result.size() == coarseGrid.number_of_nodes());

    assign(result, 0.0);

    // ---------------------- //
    // OpenMP Parallelization //
    // ---------------------- //
    const int numThreads = omp_get_max_threads();
    omp_set_num_threads(numThreads);
    const int minimalChunkSize = 4; // Could be reduced to 3 by adding more barriers.
    const int zone = 2;

    // Distribute Tasks to each thread
    TaskDistribution CircleSmootherTasks(fineGrid.numberSmootherCircles(), minimalChunkSize, numThreads);
    TaskDistribution RadialSmootherTasks(fineGrid.ntheta(), minimalChunkSize, numThreads);

    #pragma omp parallel num_threads(numThreads)
    {   
        const int threadID = omp_get_thread_num();
        // ---------------------------------------------------------- //
        // Take care of the speration strips of the circular smoother //
        // ---------------------------------------------------------- //
        const int i_r_start = CircleSmootherTasks.getStart(threadID);
        const int i_r_end = CircleSmootherTasks.getEnd(threadID);
        const int i_r_separate = std::min(i_r_end - i_r_start, zone);

        // For loop matches circular access pattern
        for (int i_r = i_r_end - i_r_separate; i_r < i_r_end; i_r++){
            int i_r_coarse = i_r >> 1;
            for (int i_theta = 0; i_theta < fineGrid.ntheta(); i_theta++){
                int i_theta_coarse = i_theta >> 1;
                FINE_NODE_GIVE_RESTRICTION(i_r, i_theta, i_r_coarse, i_theta_coarse, fineGrid, coarseGrid, result, x);
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
            int i_theta_coarse = i_theta >> 1;
            for (int i_r = fineGrid.numberSmootherCircles(); i_r < fineGrid.nr(); i_r++){
                int i_r_coarse = i_r >> 1;
                FINE_NODE_GIVE_RESTRICTION(i_r, i_theta, i_r_coarse, i_theta_coarse, fineGrid, coarseGrid, result, x);
            }
        }

        #pragma omp barrier

        // ------------------------------------------ //
        // Take care of the circular smoother section //
        // ------------------------------------------ //
        // For loop matches circular access pattern
        for (int i_r = i_r_start; i_r < i_r_end - i_r_separate; i_r++){
            int i_r_coarse = i_r >> 1;
            for (int i_theta = 0; i_theta < fineGrid.ntheta(); i_theta++){
                int i_theta_coarse = i_theta >> 1;
                FINE_NODE_GIVE_RESTRICTION(i_r, i_theta, i_r_coarse, i_theta_coarse, fineGrid, coarseGrid, result, x);
            }
        }

        // ---------------------------------------- //
        // Take care of the radial smoother section //
        // ---------------------------------------- //
        // For loop matches radial access pattern
        for (int i_theta = i_theta_start + i_theta_seperate; i_theta < i_theta_end; i_theta++){
            int i_theta_coarse = i_theta >> 1;
            for (int i_r = fineGrid.numberSmootherCircles(); i_r < fineGrid.nr(); i_r++){
                int i_r_coarse = i_r >> 1;
                FINE_NODE_GIVE_RESTRICTION(i_r, i_theta, i_r_coarse, i_theta_coarse, fineGrid, coarseGrid, result, x);
            }
        }
    }
}

// ----------------------------------------------- //
// Task dependency version of applyRestrictionGive //
// ----------------------------------------------- //

// The cost of task creation (when including dependencies) increases proportionally with the number of threads utilized.
// Given the abundance of tasks with low computational complexity, the overhead of task creation with dependencies
// significantly surpasses the time required for actual computation.
// It is therefore not recommended to use tasks with dependencies for our usage.
void Interpolation::applyRestrictionGiveTasks(const Level& fromLevel, const Level& toLevel, Vector<scalar_t>& result, const Vector<scalar_t>& x) const{
    assert(toLevel.level() == fromLevel.level() + 1);

    const PolarGrid& fineGrid = fromLevel.grid();
    const PolarGrid& coarseGrid = toLevel.grid();

    assert(x.size() == fineGrid.number_of_nodes());
    assert(result.size() == coarseGrid.number_of_nodes());

    assign(result, 0.0);

    // ---------------------- //
    // OpenMP Parallelization //
    // ---------------------- //
    const int numThreads = omp_get_max_threads();
    omp_set_num_threads(numThreads);

    const int S1 = fineGrid.numberSmootherCircles();
    const int S2 = fineGrid.ntheta();
    const int T = S1 + S2;

    int* dep = new int[T];

    const int S2_wait = S2 % 3;
    const int S2_start = std::max(S1-2, 0);

    #pragma omp parallel shared(dep)
    {
        #pragma omp single
        {
            for(int i_r = S1 - 1; i_r >= 0; i_r -= 3) {
                #pragma omp task shared(dep) firstprivate(i_r) depend(out: dep[i_r])
                {
                    int i_r_coarse = i_r >> 1;
                    for (int i_theta = 0; i_theta < fineGrid.ntheta(); i_theta++){
                        int i_theta_coarse = i_theta >> 1;
                        FINE_NODE_GIVE_RESTRICTION(i_r, i_theta, i_r_coarse, i_theta_coarse, fineGrid, coarseGrid, result, x);
                    }
                }
            }

            for(int i_r = S1 - 2; i_r >= 0; i_r -= 3) {
                #pragma omp task shared(dep) firstprivate(i_r) depend(in: dep[i_r-2], dep[i_r+1]) depend(out: dep[i_r])
                {
                    int i_r_coarse = i_r >> 1;
                    for (int i_theta = 0; i_theta < fineGrid.ntheta(); i_theta++){
                        int i_theta_coarse = i_theta >> 1;
                        FINE_NODE_GIVE_RESTRICTION(i_r, i_theta, i_r_coarse, i_theta_coarse, fineGrid, coarseGrid, result, x);
                    }
                }
            }

            for(int i_r = S1 - 3; i_r >= 0; i_r -= 3) {
                #pragma omp task shared(dep) firstprivate(i_r) depend(in: dep[i_r-2], dep[i_r+1]) depend(out: dep[i_r])
                {
                    int i_r_coarse = i_r >> 1;
                    for (int i_theta = 0; i_theta < fineGrid.ntheta(); i_theta++){
                        int i_theta_coarse = i_theta >> 1;
                        FINE_NODE_GIVE_RESTRICTION(i_r, i_theta, i_r_coarse, i_theta_coarse, fineGrid, coarseGrid, result, x);
                    }
                }
            }

            for(int i_theta = 0; i_theta < S2 - S2_wait; i_theta += 3) {
                #pragma omp task shared(dep) firstprivate(i_theta) depend(out: dep[S1+i_theta]) depend(in: dep[S2_start])
                {
                    int i_theta_coarse = i_theta >> 1;
                    for (int i_r = fineGrid.numberSmootherCircles(); i_r < fineGrid.nr(); i_r++){
                        int i_r_coarse = i_r >> 1;
                        FINE_NODE_GIVE_RESTRICTION(i_r, i_theta, i_r_coarse, i_theta_coarse, fineGrid, coarseGrid, result, x);
                    }
                }
            }

            for(int i_theta = 1; i_theta < S2 - S2_wait; i_theta += 3) {
                #pragma omp task firstprivate(i_theta) depend(out: dep[S1+i_theta]) depend(in: dep[S2_start], dep[S1+i_theta-1], dep[S1+i_theta+2])
                {
                    int i_theta_coarse = i_theta >> 1;
                    for (int i_r = fineGrid.numberSmootherCircles(); i_r < fineGrid.nr(); i_r++){
                        int i_r_coarse = i_r >> 1;
                        FINE_NODE_GIVE_RESTRICTION(i_r, i_theta, i_r_coarse, i_theta_coarse, fineGrid, coarseGrid, result, x);
                    }
                }
            }

            for(int i_theta = 2; i_theta < S2 - S2_wait; i_theta += 3) {
                #pragma omp task shared(dep) firstprivate(i_theta) depend(out: dep[S1+i_theta]) depend(in: dep[S2_start], dep[S1+i_theta-1], dep[S1+i_theta+2])
                {
                    int i_theta_coarse = i_theta >> 1;
                    for (int i_r = fineGrid.numberSmootherCircles(); i_r < fineGrid.nr(); i_r++){
                        int i_r_coarse = i_r >> 1;
                        FINE_NODE_GIVE_RESTRICTION(i_r, i_theta, i_r_coarse, i_theta_coarse, fineGrid, coarseGrid, result, x);
                    }
                }
            }

            if( S2_wait >= 1 ){
                int i_theta = S2 - S2_wait;
                #pragma omp task shared(dep) firstprivate(i_theta) depend(out: dep[S1+i_theta]) depend(in: dep[S2_start],  dep[S1+i_theta-1])
                {
                    int i_theta_coarse = i_theta >> 1;
                    for (int i_r = fineGrid.numberSmootherCircles(); i_r < fineGrid.nr(); i_r++){
                        int i_r_coarse = i_r >> 1;
                        FINE_NODE_GIVE_RESTRICTION(i_r, i_theta, i_r_coarse, i_theta_coarse, fineGrid, coarseGrid, result, x);
                    }
                }
            }

            if( S2_wait >= 2 ){
                int i_theta = S2 - S2_wait + 1;
                #pragma omp task shared(dep) firstprivate(i_theta) depend(out: dep[S1+i_theta]) depend(in: dep[S2_start],  dep[S1+i_theta-1])
                {
                    int i_theta_coarse = i_theta >> 1;
                    for (int i_r = fineGrid.numberSmootherCircles(); i_r < fineGrid.nr(); i_r++){
                        int i_r_coarse = i_r >> 1;
                        FINE_NODE_GIVE_RESTRICTION(i_r, i_theta, i_r_coarse, i_theta_coarse, fineGrid, coarseGrid, result, x);
                    }
                }
            }
        }
    }
    delete[] dep;
}


