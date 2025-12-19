#include "../../include/Interpolation/interpolation.h"

// We use the anisotropic bilinear interpolation stencil.
// For an isotropic mesh, this stencil reduces to
//           |1  2  1|
// P = 1/4 * |2  4  2|
//           |1  2  1|

void Interpolation::applyProlongation0(const PolarGrid& coarse_grid, const PolarGrid& fine_grid,
                                       Vector<double> fine_result, ConstVector<double> coarse_values) const
{
    assert(coarse_values.size() == static_cast<uint>(coarse_grid.numberOfNodes()));
    assert(fine_result.size() == static_cast<uint>(fine_grid.numberOfNodes()));

#pragma omp parallel for num_threads(max_omp_threads_)
    for (int index = 0; index < fine_grid.numberOfNodes(); index++) {
        std::array<std::pair<double, double>, space_dimension> neighbor_distance;

        MultiIndex fine_node = fine_grid.multiIndex(index);
        MultiIndex coarse_node(fine_node[0] / 2, fine_node[1] / 2); // Nearest lower left coarse node in the fine grid.

        if (fine_node[0] % 2 == 0 && fine_node[1] % 2 == 0) {
            // Fine node appears in coarse grid
            fine_result[index] = coarse_values[coarse_grid.index(coarse_node)];
        }

        if (fine_node[0] % 2 == 0 && fine_node[1] % 2 == 1) {
            // Fine node between two coarse nodes in theta direction
            // X
            // |
            // O
            // |
            // X
            fine_grid.adjacentNeighborDistances(fine_node, neighbor_distance);

            double k1 = neighbor_distance[1].first;
            double k2 = neighbor_distance[1].second;

            MultiIndex bottomNeighbor(coarse_node[0], coarse_node[1]);
            MultiIndex topNeighbor(coarse_node[0], (coarse_node[1] + 1) % coarse_grid.ntheta());

            fine_result[index] = (k1 * coarse_values[coarse_grid.index(bottomNeighbor)] +
                                  k2 * coarse_values[coarse_grid.index(topNeighbor)]) /
                                 (k1 + k2);
        }

        if (fine_node[0] % 2 == 1 && fine_node[1] % 2 == 0) {
            // Fine node between two coarse nodes in radial direction
            // X -- O -- X
            fine_grid.adjacentNeighborDistances(fine_node, neighbor_distance);

            double h1 = neighbor_distance[0].first;
            double h2 = neighbor_distance[0].second;

            MultiIndex leftNeighbor(coarse_node[0], coarse_node[1]);
            MultiIndex rightNeighbor(coarse_node[0] + 1, coarse_node[1]);

            fine_result[index] = (h1 * coarse_values[coarse_grid.index(leftNeighbor)] +
                                  h2 * coarse_values[coarse_grid.index(rightNeighbor)]) /
                                 (h1 + h2);
        }

        if (fine_node[0] % 2 == 1 && fine_node[1] % 2 == 1) {
            // Interpolates a fine node value based on four neighboring coarse nodes.
            // Fine node lies in the center of four coarse nodes forming a cross shape:
            //
            //           X     X
            /*            \   /                         */
            //              O <-- Fine Node (i_r, i_theta)
            /*            /   \                         */
            //           X     X
            //
            fine_grid.adjacentNeighborDistances(fine_node, neighbor_distance);

            double h1 = neighbor_distance[0].first;
            double h2 = neighbor_distance[0].second;
            double k1 = neighbor_distance[1].first;
            double k2 = neighbor_distance[1].second;

            MultiIndex bottom_left_neighbor(coarse_node[0], coarse_node[1]);
            MultiIndex bottom_right_neighbor(coarse_node[0] + 1, coarse_node[1]);
            MultiIndex top_left_neighbor(coarse_node[0], (coarse_node[1] + 1) % coarse_grid.ntheta());
            MultiIndex top_right_neighbor(coarse_node[0] + 1, (coarse_node[1] + 1) % coarse_grid.ntheta());

            fine_result[index] = (h1 * k1 * coarse_values[coarse_grid.index(bottom_left_neighbor)] +
                                  h2 * k1 * coarse_values[coarse_grid.index(bottom_right_neighbor)] +
                                  h1 * k2 * coarse_values[coarse_grid.index(top_left_neighbor)] +
                                  h2 * k2 * coarse_values[coarse_grid.index(top_right_neighbor)]) /
                                 ((h1 + h2) * (k1 + k2));
        }
    }
}

// --------------------------------------- //
// Optimized version of applyProlongation0 //
// --------------------------------------- //

#define FINE_NODE_PROLONGATION()                                                                                       \
    do {                                                                                                               \
        if (i_r & 1) {                                                                                                 \
            if (i_theta & 1) {                                                                                         \
                /* i_r % 2 == 1, i_theta % 2 == 1 */                                                                   \
                /* Fine node in the center of four coarse nodes */                                                     \
                double h1             = fine_grid.radialSpacing(i_r - 1);                                              \
                double h2             = fine_grid.radialSpacing(i_r);                                                  \
                double k1             = fine_grid.angularSpacing(i_theta - 1);                                         \
                double k2             = fine_grid.angularSpacing(i_theta);                                             \
                int i_theta_coarse_P1 = coarse_grid.wrapThetaIndex(i_theta_coarse + 1);                                \
                double divisor        = (h1 + h2) * (k1 + k2);                                                         \
                double value =                                                                                         \
                    (h1 * k1 * coarse_values[coarse_grid.index(i_r_coarse, i_theta_coarse)] + /* Bottom left */        \
                     h2 * k1 * coarse_values[coarse_grid.index(i_r_coarse + 1, i_theta_coarse)] + /* Bottom right */   \
                     h1 * k2 * coarse_values[coarse_grid.index(i_r_coarse, i_theta_coarse_P1)] + /* Top left */        \
                     h2 * k2 * coarse_values[coarse_grid.index(i_r_coarse + 1, i_theta_coarse_P1)] /* Top right */     \
                    );                                                                                                 \
                fine_result[fine_grid.index(i_r, i_theta)] = value / divisor;                                          \
            }                                                                                                          \
            else {                                                                                                     \
                /* i_r % 2 == 1, i_theta % 2 == 0 */                                                                   \
                /* Fine node between coarse nodes in radial direction */                                               \
                double h1      = fine_grid.radialSpacing(i_r - 1);                                                     \
                double h2      = fine_grid.radialSpacing(i_r);                                                         \
                double divisor = (h1 + h2);                                                                            \
                double value   = (h1 * coarse_values[coarse_grid.index(i_r_coarse, i_theta_coarse)] + /* left */       \
                                h2 * coarse_values[coarse_grid.index(i_r_coarse + 1, i_theta_coarse)] /* right */    \
                );                                                                                                   \
                fine_result[fine_grid.index(i_r, i_theta)] = value / divisor;                                          \
            }                                                                                                          \
        }                                                                                                              \
        else {                                                                                                         \
            if (i_theta & 1) {                                                                                         \
                /* i_r % 2 == 0, i_theta % 2 == 1 */                                                                   \
                /* Fine node between coarse nodes in theta direction */                                                \
                double k1             = fine_grid.angularSpacing(i_theta - 1);                                         \
                double k2             = fine_grid.angularSpacing(i_theta);                                             \
                int i_theta_coarse_P1 = coarse_grid.wrapThetaIndex(i_theta_coarse + 1);                                \
                double divisor        = (k1 + k2);                                                                     \
                double value = (k1 * coarse_values[coarse_grid.index(i_r_coarse, i_theta_coarse)] + /* bottom */       \
                                k2 * coarse_values[coarse_grid.index(i_r_coarse, i_theta_coarse_P1)] /* top */         \
                );                                                                                                     \
                fine_result[fine_grid.index(i_r, i_theta)] = value / divisor;                                          \
            }                                                                                                          \
            else {                                                                                                     \
                /* i_r % 2 == 0, i_theta % 2 == 0 */                                                                   \
                /* Fine node appears in coarse grid */                                                                 \
                fine_result[fine_grid.index(i_r, i_theta)] =                                                           \
                    coarse_values[coarse_grid.index(i_r_coarse, i_theta_coarse)]; /* center */                         \
            }                                                                                                          \
        }                                                                                                              \
    } while (0)

void Interpolation::applyProlongation(const PolarGrid& coarse_grid, const PolarGrid& fine_grid,
                                      Vector<double> fine_result, ConstVector<double> coarse_values) const
{
    assert(coarse_values.size() == static_cast<uint>(coarse_grid.numberOfNodes()));
    assert(fine_result.size() == static_cast<uint>(fine_grid.numberOfNodes()));

#pragma omp parallel num_threads(max_omp_threads_)
    {
/* Circluar Indexing Section */
/* For loop matches circular access pattern */
#pragma omp for nowait
        for (int i_r = 0; i_r < fine_grid.numberSmootherCircles(); i_r++) {
            int i_r_coarse = i_r / 2;
            for (int i_theta = 0; i_theta < fine_grid.ntheta(); i_theta++) {
                int i_theta_coarse = i_theta / 2;
                FINE_NODE_PROLONGATION();
            }
        }

/* Radial Indexing Section */
/* For loop matches radial access pattern */
#pragma omp for nowait
        for (int i_theta = 0; i_theta < fine_grid.ntheta(); i_theta++) {
            int i_theta_coarse = i_theta / 2;
            for (int i_r = fine_grid.numberSmootherCircles(); i_r < fine_grid.nr(); i_r++) {
                int i_r_coarse = i_r / 2;
                FINE_NODE_PROLONGATION();
            }
        }
    }
}
