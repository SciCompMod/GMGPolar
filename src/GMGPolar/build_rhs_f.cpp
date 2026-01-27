#include "../../include/GMGPolar/gmgpolar.h"

void IGMGPolar::build_rhs_f(const Level& level, Vector<double> rhs_f, const BoundaryConditions& boundary_conditions,
                            const SourceTerm& source_term)
{
    const PolarGrid& grid = level.grid();
    assert(std::ssize(rhs_f) == grid.numberOfNodes());

#pragma omp parallel
    {
// ----------------------------------------- //
// Store rhs values (circular index section) //
// ----------------------------------------- //
#pragma omp for nowait
        for (int i_r = 0; i_r < grid.numberSmootherCircles(); i_r++) {
            double r = grid.radius(i_r);
            for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
                double theta = grid.theta(i_theta);

                if ((0 < i_r && i_r < grid.nr() - 1) || (i_r == 0 && !DirBC_Interior_)) {
                    rhs_f[grid.index(i_r, i_theta)] = source_term(i_r, i_theta);
                }
                else if (i_r == 0 && DirBC_Interior_) {
                    rhs_f[grid.index(i_r, i_theta)] = boundary_conditions.u_D_Interior(r, theta);
                }
                else if (i_r == grid.nr() - 1) {
                    rhs_f[grid.index(i_r, i_theta)] = boundary_conditions.u_D(r, theta);
                }
            }
        }

// --------------------------------------- //
// Store rhs values (radial index section) //
// --------------------------------------- //
#pragma omp for
        for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
            double theta = grid.theta(i_theta);

            for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++) {
                double r = grid.radius(i_r);
                if ((0 < i_r && i_r < grid.nr() - 1) || (i_r == 0 && !DirBC_Interior_)) {
                    rhs_f[grid.index(i_r, i_theta)] = source_term(i_r, i_theta);
                }
                else if (i_r == 0 && DirBC_Interior_) {
                    rhs_f[grid.index(i_r, i_theta)] = boundary_conditions.u_D_Interior(r, theta);
                }
                else if (i_r == grid.nr() - 1) {
                    rhs_f[grid.index(i_r, i_theta)] = boundary_conditions.u_D(r, theta);
                }
            }
        }
    }
}
