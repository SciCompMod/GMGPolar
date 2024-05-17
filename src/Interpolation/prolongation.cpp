#include "../../include/Interpolation/interpolation.h"

// We use the anisotropic bilinear interpolation stencil.
// For an isotropic mesh, this stencil reduces to
//           |1  2  1|
// P = 1/4 * |2  4  2|
//           |1  2  1|

void Interpolation::applyProlongation0(const Level& fromLevel, const Level& toLevel, Vector<scalar_t>& result, const Vector<scalar_t>& x) const{
    assert(toLevel.level() == fromLevel.level() - 1);

    const PolarGrid& coarseGrid = fromLevel.grid();
    const PolarGrid& fineGrid = toLevel.grid();

    assert(x.size() == coarseGrid.number_of_nodes());
    assert(result.size() == fineGrid.number_of_nodes());

    #pragma omp parallel for
    for (int index = 0; index < fineGrid.number_of_nodes(); index++) {
        std::array<std::pair<scalar_t,scalar_t>, space_dimension> neighbor_distance;

        MultiIndex fine_node = fineGrid.multiindex(index);
        MultiIndex coarse_node(fine_node[0] / 2, fine_node[1] / 2); // Nearest lower left coarse node in the fine grid.

        if(fine_node[0] % 2 == 0 && fine_node[1] % 2 == 0){
            // Fine node appears in coarse grid
            result[index] = x[coarseGrid.index(coarse_node)];
        }

        if(fine_node[0] % 2 == 0 && fine_node[1] % 2 == 1){
            // Fine node between two coarse nodes in theta direction
            // X
            // |
            // O
            // |
            // X
            fineGrid.adjacent_neighbor_distances(fine_node, neighbor_distance);

            scalar_t k1 = neighbor_distance[1].first;
            scalar_t k2 = neighbor_distance[1].second;

            MultiIndex bottomNeighbor(coarse_node[0], coarse_node[1]);
            MultiIndex topNeighbor(coarse_node[0], (coarse_node[1] + 1) % coarseGrid.ntheta());

            result[index] = 
                (k1 * x[coarseGrid.index(bottomNeighbor)] + 
                k2 * x[coarseGrid.index(topNeighbor)]) / (k1 + k2);  
        }

        if(fine_node[0] % 2 == 1 && fine_node[1] % 2 == 0){
            // Fine node between two coarse nodes in radial direction
            // X -- O -- X
            fineGrid.adjacent_neighbor_distances(fine_node, neighbor_distance);

            scalar_t h1 = neighbor_distance[0].first;
            scalar_t h2 = neighbor_distance[0].second;

            MultiIndex leftNeighbor(coarse_node[0], coarse_node[1]);
            MultiIndex rightNeighbor(coarse_node[0] + 1, coarse_node[1]);

            result[index] = 
                (h1 * x[coarseGrid.index(leftNeighbor)] + 
                h2 * x[coarseGrid.index(rightNeighbor)]) / (h1 + h2); 
        }

        if(fine_node[0] % 2 == 1 && fine_node[1] % 2 == 1){   
            // Interpolates a fine node value based on four neighboring coarse nodes.
            // Fine node lies in the center of four coarse nodes forming a cross shape:
            //
            //           X     X
            /*            \   /                         */
            //              O <-- Fine Node (i_r, i_theta)
            /*            /   \                         */
            //           X     X
            //
            fineGrid.adjacent_neighbor_distances(fine_node, neighbor_distance);

            scalar_t h1 = neighbor_distance[0].first;
            scalar_t h2 = neighbor_distance[0].second;
            scalar_t k1 = neighbor_distance[1].first;
            scalar_t k2 = neighbor_distance[1].second;
            
            MultiIndex bottom_left_neighbor(coarse_node[0], coarse_node[1]);
            MultiIndex bottom_right_neighbor(coarse_node[0] + 1, coarse_node[1]);
            MultiIndex top_left_neighbor(coarse_node[0], (coarse_node[1] + 1) % coarseGrid.ntheta());
            MultiIndex top_right_neighbor(coarse_node[0] + 1, (coarse_node[1] + 1) % coarseGrid.ntheta());

            result[index] = 
                (h1 * k1 * x[coarseGrid.index(bottom_left_neighbor)] +
                h2 * k1 * x[coarseGrid.index(bottom_right_neighbor)] +
                h1 * k2 * x[coarseGrid.index(top_left_neighbor)] +
                h2 * k2 * x[coarseGrid.index(top_right_neighbor)])
                / ((h1 + h2) * (k1 + k2));
        }
    }
}


// --------------------------------------- //
// Optimized version of applyProlongation0 //
// --------------------------------------- //

#define FINE_NODE_TAKE_PROLONGATION() \
do { \
    if(i_r & 1) { \
        if(i_theta & 1) { \
            /* i_r % 2 == 1, i_theta % 2 == 1 */ \
            /* Fine node in the center of four coarse nodes */ \
            scalar_t h1 = fineGrid.r_dist(i_r-1); \
            scalar_t h2 = fineGrid.r_dist(i_r); \
            scalar_t k1 = fineGrid.theta_dist(i_theta-1); \
            scalar_t k2 = fineGrid.theta_dist(i_theta); \
            scalar_t divisor = (h1+h2) * (k1+k2); \
            scalar_t value = ( \
                h1*k1*x[coarseGrid.index(i_r_coarse,   i_theta_coarse)] +  /* Bottom left */ \
                h2*k1*x[coarseGrid.index(i_r_coarse+1, i_theta_coarse)] +  /* Bottom right */ \
                h1*k2*x[coarseGrid.index(i_r_coarse,   i_theta_coarse+1)] + /* Top left */ \
                h2*k2*x[coarseGrid.index(i_r_coarse+1, i_theta_coarse+1)]   /* Top right */ \
            ); \
            result[fineGrid.index(i_r, i_theta)] = value / divisor; \
        } \
        else { \
            /* i_r % 2 == 1, i_theta % 2 == 0 */ \
            /* Fine node between coarse nodes in radial direction */ \
            scalar_t h1 = fineGrid.r_dist(i_r-1); \
            scalar_t h2 = fineGrid.r_dist(i_r); \
            scalar_t divisor = (h1+h2); \
            scalar_t value = ( \
                h1*x[coarseGrid.index(i_r_coarse,   i_theta_coarse)] +  /* left */ \
                h2*x[coarseGrid.index(i_r_coarse+1, i_theta_coarse)]    /* right */ \
            ); \
            result[fineGrid.index(i_r, i_theta)] = value / divisor; \
        } \
    } \
    else { \
        if(i_theta & 1) { \
            /* i_r % 2 == 0, i_theta % 2 == 1 */ \
            /* Fine node between coarse nodes in theta direction */ \
            scalar_t k1 = fineGrid.theta_dist(i_theta-1); \
            scalar_t k2 = fineGrid.theta_dist(i_theta); \
            scalar_t divisor = (k1+k2); \
            scalar_t value = ( \
                k1*x[coarseGrid.index(i_r_coarse, i_theta_coarse)] +  /* bottom */ \
                k2*x[coarseGrid.index(i_r_coarse, i_theta_coarse+1)]   /* top */ \
            ); \
            result[fineGrid.index(i_r, i_theta)] = value / divisor; \
        } \
        else { \
            /* i_r % 2 == 0, i_theta % 2 == 0 */ \
            /* Fine node appears in coarse grid */ \
            result[fineGrid.index(i_r, i_theta)] = x[coarseGrid.index(i_r_coarse, i_theta_coarse)]; /* center */ \
        } \
    } \
} while(0)

void Interpolation::applyProlongation(const Level& fromLevel, const Level& toLevel, Vector<scalar_t>& result, const Vector<scalar_t>& x) const{
    assert(toLevel.level() == fromLevel.level() - 1);

    const PolarGrid& coarseGrid = fromLevel.grid();
    const PolarGrid& fineGrid = toLevel.grid();

    assert(x.size() == coarseGrid.number_of_nodes());
    assert(result.size() == fineGrid.number_of_nodes());

    // ---------------------- //
    // OpenMP Parallelization //
    // ---------------------- //
    const int numThreads = omp_get_max_threads();
    omp_set_num_threads(numThreads);

    #pragma omp parallel
    {
        // Circular Smoother section
        // For loop matches circular access pattern
        #pragma omp for nowait
        for (int i_r = 0; i_r < fineGrid.numberSmootherCircles(); i_r++){
            int i_r_coarse = i_r >> 1;
            for (int i_theta = 0; i_theta < fineGrid.ntheta(); i_theta++){
                int i_theta_coarse = i_theta >> 1;
                FINE_NODE_TAKE_PROLONGATION();
            }
        }

        // Radial smoother section
        // For loop matches radial access pattern
        #pragma omp for
        for (int i_theta = 0; i_theta < fineGrid.ntheta(); i_theta++){
            int i_theta_coarse = i_theta >> 1;
            for (int i_r = fineGrid.numberSmootherCircles(); i_r < fineGrid.nr(); i_r++){
                int i_r_coarse = i_r >> 1;
                FINE_NODE_TAKE_PROLONGATION();
            }
        }
    }
}








// --------------------------------------- //
// Optimized version of applyProlongation0 // With custom parallelization
// --------------------------------------- //

// void Interpolation::applyProlongation(const Level& fromLevel, const Level& toLevel, Vector<scalar_t>& result, const Vector<scalar_t>& x) const{
//     assert(toLevel.level() == fromLevel.level() - 1);

//     const PolarGrid& coarseGrid = fromLevel.grid();
//     const PolarGrid& fineGrid = toLevel.grid();

//     assert(x.size() == coarseGrid.number_of_nodes());
//     assert(result.size() == fineGrid.number_of_nodes());

//     // ---------------------- //
//     // OpenMP Parallelization //
//     // ---------------------- //
//     const int numThreads = omp_get_max_threads();
//     omp_set_num_threads(numThreads);
//     const int minimalChunkSize = 1;

//     // Distribute Tasks to each thread
//     TaskDistribution CircleSmootherTasks(fineGrid.numberSmootherCircles(), minimalChunkSize, numThreads);
//     TaskDistribution RadialSmootherTasks(fineGrid.ntheta(), minimalChunkSize, numThreads);

//     #pragma omp parallel num_threads(numThreads)
//     {   
//         int threadID = omp_get_thread_num();
//         // ------------------------------------------ //
//         // Take care of the circular Smoother section //
//         // ------------------------------------------ //
//         const int i_r_start = CircleSmootherTasks.getStart(threadID);
//         const int i_r_end = CircleSmootherTasks.getEnd(threadID);
//         // For loop matches circular access pattern
//         for (int i_r = i_r_start; i_r < i_r_end; i_r++){
//             int i_r_coarse = i_r >> 1;
//             for (int i_theta = 0; i_theta < fineGrid.ntheta(); i_theta++){
//                 int i_theta_coarse = i_theta >> 1;
//                 FINE_NODE_TAKE_PROLONGATION();
//             }
//         }
//         // ---------------------------------------- //
//         // Take care of the radial smoother section //
//         // ---------------------------------------- //
//         threadID = numThreads - threadID - 1; // Distribute work more evenly
//         const int i_theta_start = RadialSmootherTasks.getStart(threadID);
//         const int i_theta_end = RadialSmootherTasks.getEnd(threadID);
//         // For loop matches radial access pattern
//         for (int i_theta = i_theta_start; i_theta < i_theta_end; i_theta++){
//             int i_theta_coarse = i_theta >> 1;
//             for (int i_r = fineGrid.numberSmootherCircles(); i_r < fineGrid.nr(); i_r++){
//                 int i_r_coarse = i_r >> 1;
//                 FINE_NODE_TAKE_PROLONGATION();
//             }
//         }
//     }
// }
