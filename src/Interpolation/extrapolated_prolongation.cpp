#include "../../include/Interpolation/interpolation.h"

void Interpolation::applyExtrapolatedProlongation0(
    const Level& fromLevel, const Level& toLevel, Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> result,
    const Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> x) const
{
    assert(toLevel.level_depth() == fromLevel.level_depth() - 1);

    const PolarGrid& coarseGrid = fromLevel.grid();
    const PolarGrid& fineGrid   = toLevel.grid();

    assert(x.size() == coarseGrid.numberOfNodes());
    assert(result.size() == fineGrid.numberOfNodes());

#pragma omp parallel for num_threads(threads_per_level_[toLevel.level_depth()])
    for (int index = 0; index < fineGrid.numberOfNodes(); index++) {
        std::array<std::pair<double, double>, space_dimension> neighbor_distance;

        MultiIndex fine_node = fineGrid.multiIndex(index);
        MultiIndex coarse_node(fine_node[0] / 2, fine_node[1] / 2); // Nearest lower left coarse node in the fine grid.

        if (fine_node[0] % 2 == 0 && fine_node[1] % 2 == 0) {
            // Fine node appears in coarse grid
            result[index] = x[coarseGrid.index(coarse_node)];
        }

        if (fine_node[0] % 2 == 0 && fine_node[1] % 2 == 1) {
            // Fine node between two coarse nodes in theta direction
            // X
            // |
            // O
            // |
            // X
            MultiIndex bottomNeighbor(coarse_node[0], coarse_node[1]);
            MultiIndex topNeighbor(coarse_node[0], (coarse_node[1] + 1) % coarseGrid.ntheta());
            result[index] = 0.5 * (x[coarseGrid.index(bottomNeighbor)] + x[coarseGrid.index(topNeighbor)]);
        }

        if (fine_node[0] % 2 == 1 && fine_node[1] % 2 == 0) {
            // Fine node between two coarse nodes in radial direction
            // X -- O -- X
            MultiIndex leftNeighbor(coarse_node[0], coarse_node[1]);
            MultiIndex rightNeighbor(coarse_node[0] + 1, coarse_node[1]);
            result[index] = 0.5 * (x[coarseGrid.index(leftNeighbor)] + x[coarseGrid.index(rightNeighbor)]);
        }

        if (fine_node[0] % 2 == 1 && fine_node[1] % 2 == 1) {
            // Interpolates a fine node value based on two neighboring coarse nodes.
            // Fine node lies in the center of four coarse nodes forming a cross shape:
            //
            //           X
            /*            \                             */
            //              O <-- Fine Node (i_r, i_theta)
            /*                \                         */
            //                 X
            //
            MultiIndex bottom_right_neighbor(coarse_node[0] + 1, coarse_node[1]);
            MultiIndex top_left_neighbor(coarse_node[0], (coarse_node[1] + 1) % coarseGrid.ntheta());
            result[index] = 0.5 * (x[coarseGrid.index(bottom_right_neighbor)] + x[coarseGrid.index(top_left_neighbor)]);
        }
    }
}

// --------------------------------------- //
// Optimized version of applyProlongation0 //
// --------------------------------------- //

#define FINE_NODE_EXTRAPOLATED_PROLONGATION()                                                                          \
    do {                                                                                                               \
        if (i_r & 1) {                                                                                                 \
            if (i_theta & 1) {                                                                                         \
                /* i_r % 2 == 1, i_theta % 2 == 1 */                                                                   \
                /* Fine node in the center of four coarse nodes */                                                     \
                result[fineGrid.index(i_r, i_theta)] =                                                                 \
                    0.5 * (x[coarseGrid.index(i_r_coarse + 1, i_theta_coarse)] + /* Bottom right */                    \
                           x[coarseGrid.index(i_r_coarse, i_theta_coarse + 1)] /* Top left */                          \
                          );                                                                                           \
            }                                                                                                          \
            else {                                                                                                     \
                /* i_r % 2 == 1, i_theta % 2 == 0 */                                                                   \
                /* Fine node between coarse nodes in radial direction */                                               \
                result[fineGrid.index(i_r, i_theta)] =                                                                 \
                    0.5 * (x[coarseGrid.index(i_r_coarse, i_theta_coarse)] + /* left */                                \
                           x[coarseGrid.index(i_r_coarse + 1, i_theta_coarse)] /* right */                             \
                          );                                                                                           \
            }                                                                                                          \
        }                                                                                                              \
        else {                                                                                                         \
            if (i_theta & 1) {                                                                                         \
                /* i_r % 2 == 0, i_theta % 2 == 1 */                                                                   \
                /* Fine node between coarse nodes in theta direction */                                                \
                result[fineGrid.index(i_r, i_theta)] =                                                                 \
                    0.5 * (x[coarseGrid.index(i_r_coarse, i_theta_coarse)] + /* bottom */                              \
                           x[coarseGrid.index(i_r_coarse, i_theta_coarse + 1)] /* top */                               \
                          );                                                                                           \
            }                                                                                                          \
            else {                                                                                                     \
                /* i_r % 2 == 0, i_theta % 2 == 0 */                                                                   \
                /* Fine node appears in coarse grid */                                                                 \
                result[fineGrid.index(i_r, i_theta)] = x[coarseGrid.index(i_r_coarse, i_theta_coarse)]; /* center */   \
            }                                                                                                          \
        }                                                                                                              \
    } while (0)

void Interpolation::applyExtrapolatedProlongation(
    const Level& fromLevel, const Level& toLevel, Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> result,
    const Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> x) const
{
    assert(toLevel.level_depth() == fromLevel.level_depth() - 1);

    const PolarGrid& coarseGrid = fromLevel.grid();
    const PolarGrid& fineGrid   = toLevel.grid();

    assert(x.size() == coarseGrid.numberOfNodes());
    assert(result.size() == fineGrid.numberOfNodes());

#pragma omp parallel num_threads(threads_per_level_[toLevel.level_depth()]) if (fineGrid.numberOfNodes() > 10'000)
    {
/* Circluar Indexing Section */
/* For loop matches circular access pattern */
#pragma omp for nowait
        for (int i_r = 0; i_r < fineGrid.numberSmootherCircles(); i_r++) {
            int i_r_coarse = i_r / 2;
            for (int i_theta = 0; i_theta < fineGrid.ntheta(); i_theta++) {
                int i_theta_coarse = i_theta / 2;
                FINE_NODE_EXTRAPOLATED_PROLONGATION();
            }
        }

/* Radial Indexing Section */
/* For loop matches radial access pattern */
#pragma omp for nowait
        for (int i_theta = 0; i_theta < fineGrid.ntheta(); i_theta++) {
            int i_theta_coarse = i_theta / 2;
            for (int i_r = fineGrid.numberSmootherCircles(); i_r < fineGrid.nr(); i_r++) {
                int i_r_coarse = i_r / 2;
                FINE_NODE_EXTRAPOLATED_PROLONGATION();
            }
        }
    }
}
