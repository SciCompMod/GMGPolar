#include "../../include/Interpolation/interpolation.h"

// RestrictionTake is faster than the RestrictionGive!

// For the restriction we use R = P^T.
// The restriction for an siotropic mesh reduces to
//           |1  2  1|
// R = 1/4 * |2  4  2| = P^T
//           |1  2  1|

void Interpolation::applyRestrictionGive0(const Level& fromLevel, const Level& toLevel, Vector<double>& result, const Vector<double>& x) const{
    assert(toLevel.level() == fromLevel.level() + 1);

    const PolarGrid& fineGrid = fromLevel.grid();
    const PolarGrid& coarseGrid = toLevel.grid();

    assert(x.size() == fineGrid.number_of_nodes());
    assert(result.size() == coarseGrid.number_of_nodes());
    
    assign(result, 0.0);

    for(int index = 0; index < fineGrid.number_of_nodes(); index ++){
        std::array<std::pair<double,double>, space_dimension> neighbor_distance;
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
            double k1 = neighbor_distance[1].first;
            double k2 = neighbor_distance[1].second;

            fineGrid.adjacent_neighbors_of(fine_node, neighbors);

            MultiIndex bottom_coarse_node(fine_node[0] / 2, fine_node[1] / 2);
            MultiIndex top_coarse_node(fine_node[0] / 2, (fine_node[1] / 2 + 1) % coarseGrid.ntheta());

            result[coarseGrid.index(bottom_coarse_node)] += k1 *x[index] / (k1 + k2);
            result[coarseGrid.index(top_coarse_node)] += k2 *x[index]/ (k1+k2);
        }


        // Fine node between coarse nodes in radial direction
        if(fine_node[0] % 2 == 1 && fine_node[1] % 2 == 0){
            fineGrid.adjacent_neighbor_distances(fine_node, neighbor_distance);
            double h1 = neighbor_distance[0].first;
            double h2 = neighbor_distance[0].second;

            fineGrid.adjacent_neighbors_of(fine_node, neighbors);

            MultiIndex left_coarse_node(fine_node[0] / 2, fine_node[1] / 2);
            MultiIndex right_coarse_node(fine_node[0] / 2 + 1, fine_node[1] / 2);

            result[coarseGrid.index(left_coarse_node)] +=  (h1 *x[index] ) / (h1+h2);;
            result[coarseGrid.index(right_coarse_node)] +=  (h2 *x[index] ) / (h1+h2);;
        }


        //Fine node in the center of four coarse nodes
        if(fine_node[0] % 2 == 1 && fine_node[1] % 2 == 1){

            fineGrid.adjacent_neighbor_distances(fine_node, neighbor_distance);
            double h1 = neighbor_distance[0].first;
            double h2 = neighbor_distance[0].second;
            double k1 = neighbor_distance[1].first;
            double k2 = neighbor_distance[1].second;

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
            double h1 = fineGrid.r_dist(i_r-1); \
            double h2 = fineGrid.r_dist(i_r); \
            double k1 = fineGrid.theta_dist(i_theta-1); \
            double k2 = fineGrid.theta_dist(i_theta); \
            double divisor = (h1+h2) * (k1+k2); \
            double value = x[fineGrid.index(i_r, i_theta)]; \
            /* Bottom Left, Bottom Right, Top Left, Top Right */ \
            result[coarseGrid.index(i_r_coarse,   i_theta_coarse)] += h1*k1 * value / divisor; \
            result[coarseGrid.index(i_r_coarse+1, i_theta_coarse)] += h2*k1 * value / divisor; \
            result[coarseGrid.index(i_r_coarse,   i_theta_coarse+1)] += h1*k2 * value / divisor; \
            result[coarseGrid.index(i_r_coarse+1, i_theta_coarse+1)] += h2*k2 * value / divisor; \
        } \
        else{ \
            /* Fine node between coarse nodes in radial direction */ \
            double h1 = fineGrid.r_dist(i_r-1); \
            double h2 = fineGrid.r_dist(i_r); \
            double divisor = (h1+h2); \
            double value = x[fineGrid.index(i_r, i_theta)]; \
            /* Left; Right */ \
            result[coarseGrid.index(i_r_coarse,   i_theta_coarse)] += h1 * value / divisor; \
            result[coarseGrid.index(i_r_coarse+1, i_theta_coarse)] += h2 * value / divisor; \
        } \
    } \
    else{ \
        if(i_theta & 1){ \
            /* Fine node between coarse nodes in theta direction */ \
            double k1 = fineGrid.theta_dist(i_theta-1); \
            double k2 = fineGrid.theta_dist(i_theta); \
            double divisor = (k1+k2); \
            double value = x[fineGrid.index(i_r, i_theta)]; \
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



#define CIRCLE_SECTION_GIVE_RESTRICTION(i_r) \
do { \
    int i_r_coarse = i_r >> 1; \
    for (int i_theta = 0; i_theta < fineGrid.ntheta(); i_theta++){ \
        int i_theta_coarse = i_theta >> 1; \
        FINE_NODE_GIVE_RESTRICTION(i_r, i_theta, i_r_coarse, i_theta_coarse, fineGrid, coarseGrid, result, x); \
    } \
} while(0)


#define RADIAL_SECTION_GIVE_RESTRICTION(i_theta) \
do { \
    int i_theta_coarse = i_theta >> 1; \
    for (int i_r = fineGrid.numberSmootherCircles(); i_r < fineGrid.nr(); i_r++){ \
        int i_r_coarse = i_r >> 1; \
        FINE_NODE_GIVE_RESTRICTION(i_r, i_theta, i_r_coarse, i_theta_coarse, fineGrid, coarseGrid, result, x); \
    } \
} while(0)


// ----------------------------------------------- //
// Task dependency version of applyRestrictionGive //
// ----------------------------------------------- //

void Interpolation::applyRestrictionGiveTasks(const Level& fromLevel, const Level& toLevel, Vector<double>& result, const Vector<double>& x) const{
    assert(toLevel.level() == fromLevel.level() + 1);

    omp_set_num_threads(maxOpenMPThreads_);

    const PolarGrid& fineGrid = fromLevel.grid();
    const PolarGrid& coarseGrid = toLevel.grid();

    assert(x.size() == fineGrid.number_of_nodes());
    assert(result.size() == coarseGrid.number_of_nodes());

    assign(result, 0.0);

    const int numCircleTasks = fineGrid.numberSmootherCircles();
    const int additionalRadialTasks = fineGrid.ntheta() % 3;
    const int numRadialTasks = fineGrid.ntheta() - additionalRadialTasks;

    assert(numCircleTasks >= 2);
    assert(numRadialTasks >= 3 && numRadialTasks % 3 == 0);

    /* Make sure to deallocate at the end */
    int* dep = new int [numCircleTasks + numRadialTasks];

    const int openMPTaskThreads = taskingThreads_[fromLevel.level()];
    omp_set_num_threads(openMPTaskThreads);
    #pragma omp parallel num_threads(openMPTaskThreads) /* Outside variable are shared by default */
    {
        #pragma omp single
        {
            /* ------------ */
            /* Circle Tasks */
            /* ------------ */

            /* Mod 0 Circles */
            for(int circle_task = 0; circle_task < numCircleTasks; circle_task += 3) {
                #pragma omp task \
                    depend(out: dep[circle_task])
                {
                    int i_r = fineGrid.numberSmootherCircles() - circle_task - 1;    
                    CIRCLE_SECTION_GIVE_RESTRICTION(i_r);
                }
            }
            /* Mod 2 Circles */
            for(int circle_task = 1; circle_task < numCircleTasks; circle_task += 3) {
                #pragma omp task \
                    depend(out: dep[circle_task]) \
                    depend(in: dep[circle_task-1], dep[circle_task+2])   
                {
                    int i_r = fineGrid.numberSmootherCircles() - circle_task - 1;    
                    CIRCLE_SECTION_GIVE_RESTRICTION(i_r);
                }
            }
            /* Mod 2 Circles */
            for(int circle_task = 2; circle_task < numCircleTasks; circle_task += 3) {
                #pragma omp task \
                    depend(out: dep[circle_task]) \
                    depend(in: dep[circle_task-1], dep[circle_task+2])   
                {
                    int i_r = fineGrid.numberSmootherCircles() - circle_task - 1;    
                    CIRCLE_SECTION_GIVE_RESTRICTION(i_r);
                }
            }

            /* ------------ */
            /* Radial Tasks */
            /* ------------ */

            /* Mod 0 Radials */
            for(int radial_task = 0; radial_task < numRadialTasks; radial_task += 3) {
                #pragma omp task \
                    depend(out: dep[numCircleTasks+radial_task]) \
                    depend(in: dep[1])
                {
                    if(radial_task > 0){
                        int i_theta = radial_task + additionalRadialTasks;    
                        RADIAL_SECTION_GIVE_RESTRICTION(i_theta);
                    } else{
                        if(additionalRadialTasks == 0){
                            RADIAL_SECTION_GIVE_RESTRICTION(0);
                        } 
                        else if(additionalRadialTasks >= 1){
                            RADIAL_SECTION_GIVE_RESTRICTION(0);
                            RADIAL_SECTION_GIVE_RESTRICTION(1);
                        }
                    }
                }
            }
            /* Mod 1 Radials */
            for(int radial_task = 1; radial_task < numRadialTasks; radial_task += 3) {
                #pragma omp task \
                    depend(out: dep[numCircleTasks + radial_task]) \
                    depend(in: \
                        dep[numCircleTasks + radial_task-1], \
                        dep[numCircleTasks + (radial_task+2) % numRadialTasks])   
                {
                    if(radial_task > 1){
                        int i_theta = radial_task + additionalRadialTasks;    
                        RADIAL_SECTION_GIVE_RESTRICTION(i_theta);
                    } else {
                        if(additionalRadialTasks == 0){
                            RADIAL_SECTION_GIVE_RESTRICTION(1);
                        } 
                        else if(additionalRadialTasks == 1){
                            RADIAL_SECTION_GIVE_RESTRICTION(2);
                        }
                        else if(additionalRadialTasks == 2){
                            RADIAL_SECTION_GIVE_RESTRICTION(2);
                            RADIAL_SECTION_GIVE_RESTRICTION(3);
                        }
                    }
                }
            }
            /* Mod 2 Radials */
            for(int radial_task = 2; radial_task < numRadialTasks; radial_task += 3) {
                #pragma omp task \
                    depend(out: dep[numCircleTasks + radial_task]) \
                    depend(in: \
                        dep[numCircleTasks + radial_task-1], \
                        dep[numCircleTasks + (radial_task+2) % numRadialTasks])   
                {
                    int i_theta = radial_task + additionalRadialTasks;    
                    RADIAL_SECTION_GIVE_RESTRICTION(i_theta);
                }
            }
        }
    }
    omp_set_num_threads(maxOpenMPThreads_);

    delete[] dep;
}










void Interpolation::applyRestrictionGive(const Level& fromLevel, const Level& toLevel, Vector<double>& result, const Vector<double>& x) const{
    // assert(toLevel.level() == fromLevel.level() + 1);

    // const PolarGrid& fineGrid = fromLevel.grid();
    // const PolarGrid& coarseGrid = toLevel.grid();

    // assert(x.size() == fineGrid.number_of_nodes());
    // assert(result.size() == coarseGrid.number_of_nodes());

    // assign(result, 0.0);

    // // ---------------------- //
    // // OpenMP Parallelization //
    // // ---------------------- //
    // const int numThreads = omp_get_max_threads();
    // omp_set_num_threads(numThreads);
    // const int minimalChunkSize = 4;
    // const int zone = 2;

    // // Distribute Tasks to each thread
    // TaskDistribution CircleSmootherTasks(fineGrid.numberSmootherCircles(), minimalChunkSize, numThreads);
    // TaskDistribution RadialSmootherTasks(fineGrid.ntheta(), minimalChunkSize, numThreads);

    // #pragma omp parallel num_threads(numThreads)
    // {   
    //     const int threadID = omp_get_thread_num();
    //     // ---------------------------------------------------------- //
    //     // Take care of the speration strips of the circular smoother //
    //     // ---------------------------------------------------------- //
    //     const int i_r_start = CircleSmootherTasks.getStart(threadID);
    //     const int i_r_end = CircleSmootherTasks.getEnd(threadID);
    //     const int i_r_separate = std::min(i_r_end - i_r_start, zone);

    //     // For loop matches circular access pattern
    //     for (int i_r = i_r_end - i_r_separate; i_r < i_r_end; i_r++){
    //         int i_r_coarse = i_r >> 1;
    //         for (int i_theta = 0; i_theta < fineGrid.ntheta(); i_theta++){
    //             int i_theta_coarse = i_theta >> 1;
    //             FINE_NODE_GIVE_RESTRICTION(i_r, i_theta, i_r_coarse, i_theta_coarse, fineGrid, coarseGrid, result, x);
    //         }
    //     }
        
    //     #pragma omp barrier

    //     // -------------------------------------------------------- //
    //     // Take care of the speration strips of the radial smoother //
    //     // -------------------------------------------------------- //
    //     const int i_theta_start = RadialSmootherTasks.getStart(threadID);
    //     const int i_theta_end = RadialSmootherTasks.getEnd(threadID);
    //     const int i_theta_seperate = std::min(i_theta_end-i_theta_start, zone);

    //     // For loop matches radial access pattern
    //     for (int i_theta = i_theta_start; i_theta < i_theta_start + i_theta_seperate; i_theta++){
    //         int i_theta_coarse = i_theta >> 1;
    //         for (int i_r = fineGrid.numberSmootherCircles(); i_r < fineGrid.nr(); i_r++){
    //             int i_r_coarse = i_r >> 1;
    //             FINE_NODE_GIVE_RESTRICTION(i_r, i_theta, i_r_coarse, i_theta_coarse, fineGrid, coarseGrid, result, x);
    //         }
    //     }

    //     #pragma omp barrier

    //     // ------------------------------------------ //
    //     // Take care of the circular smoother section //
    //     // ------------------------------------------ //
    //     // For loop matches circular access pattern
    //     for (int i_r = i_r_start; i_r < i_r_end - i_r_separate; i_r++){
    //         int i_r_coarse = i_r >> 1;
    //         for (int i_theta = 0; i_theta < fineGrid.ntheta(); i_theta++){
    //             int i_theta_coarse = i_theta >> 1;
    //             FINE_NODE_GIVE_RESTRICTION(i_r, i_theta, i_r_coarse, i_theta_coarse, fineGrid, coarseGrid, result, x);
    //         }
    //     }

    //     // ---------------------------------------- //
    //     // Take care of the radial smoother section //
    //     // ---------------------------------------- //
    //     // For loop matches radial access pattern
    //     for (int i_theta = i_theta_start + i_theta_seperate; i_theta < i_theta_end; i_theta++){
    //         int i_theta_coarse = i_theta >> 1;
    //         for (int i_r = fineGrid.numberSmootherCircles(); i_r < fineGrid.nr(); i_r++){
    //             int i_r_coarse = i_r >> 1;
    //             FINE_NODE_GIVE_RESTRICTION(i_r, i_theta, i_r_coarse, i_theta_coarse, fineGrid, coarseGrid, result, x);
    //         }
    //     }
    // }
}
