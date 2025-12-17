#include "../../include/GMGPolar/gmgpolar.h"

template <DomainGeometryConcept DomainGeometry>
void GMGPolar<DomainGeometry>::discretize_rhs_f(const Level& level, Vector<double> rhs_f)
{
    const PolarGrid& grid = level.grid();
    assert(rhs_f.size() == static_cast<uint>(grid.numberOfNodes()));

    if (level.levelCache().cacheDomainGeometry()) {
        /* DomainGeometry is cached */
        const auto& detDF_cache = level.levelCache().detDF();
#pragma omp parallel
        {
// ---------------------------------------------- //
// Discretize rhs values (circular index section) //
// ---------------------------------------------- //
#pragma omp for nowait
            for (int i_r = 0; i_r < grid.numberSmootherCircles(); i_r++) {
                double r = grid.radius(i_r);
                for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
                    double theta = grid.theta(i_theta);
                    if ((0 < i_r && i_r < grid.nr() - 1) || (i_r == 0 && !DirBC_Interior_)) {
                        double h1          = (i_r == 0) ? 2.0 * grid.radius(0) : grid.radialSpacing(i_r - 1);
                        double h2          = grid.radialSpacing(i_r);
                        double k1          = grid.angularSpacing(i_theta - 1);
                        double k2          = grid.angularSpacing(i_theta);
                        const double detDF = detDF_cache[grid.index(i_r, i_theta)];
                        rhs_f[grid.index(i_r, i_theta)] *= 0.25 * (h1 + h2) * (k1 + k2) * fabs(detDF);
                    }
                    else if (i_r == 0 && DirBC_Interior_) {
                        rhs_f[grid.index(i_r, i_theta)] *= 1.0;
                    }
                    else if (i_r == grid.nr() - 1) {
                        rhs_f[grid.index(i_r, i_theta)] *= 1.0;
                    }
                }
            }

// -------------------------------------------- //
// Discretize rhs values (radial index section) //
// -------------------------------------------- //
#pragma omp for nowait
            for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
                double theta = grid.theta(i_theta);
                for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++) {
                    double r = grid.radius(i_r);
                    if ((0 < i_r && i_r < grid.nr() - 1) || (i_r == 0 && !DirBC_Interior_)) {
                        double h1          = (i_r == 0) ? 2.0 * grid.radius(0) : grid.radialSpacing(i_r - 1);
                        double h2          = grid.radialSpacing(i_r);
                        double k1          = grid.angularSpacing(i_theta - 1);
                        double k2          = grid.angularSpacing(i_theta);
                        const double detDF = detDF_cache[grid.index(i_r, i_theta)];
                        rhs_f[grid.index(i_r, i_theta)] *= 0.25 * (h1 + h2) * (k1 + k2) * fabs(detDF);
                    }
                    else if (i_r == 0 && DirBC_Interior_) {
                        rhs_f[grid.index(i_r, i_theta)] *= 1.0;
                    }
                    else if (i_r == grid.nr() - 1) {
                        rhs_f[grid.index(i_r, i_theta)] *= 1.0;
                    }
                }
            }
        }
    }
    else {
        /* DomainGeometry is not cached */

#pragma omp parallel
        {
// ---------------------------------------------- //
// Discretize rhs values (circular index section) //
// ---------------------------------------------- //
#pragma omp for nowait
            for (int i_r = 0; i_r < grid.numberSmootherCircles(); i_r++) {
                double r = grid.radius(i_r);
                for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
                    double theta = grid.theta(i_theta);

                    if ((0 < i_r && i_r < grid.nr() - 1) || (i_r == 0 && !DirBC_Interior_)) {
                        double h1 = (i_r == 0) ? 2.0 * grid.radius(0) : grid.radialSpacing(i_r - 1);
                        double h2 = grid.radialSpacing(i_r);
                        double k1 = grid.angularSpacing(i_theta - 1);
                        double k2 = grid.angularSpacing(i_theta);
                        /* Calculate the elements of the Jacobian matrix for the transformation mapping */
                        /* The Jacobian matrix is: */
                        /* [Jrr, Jrt] */
                        /* [Jtr, Jtt] */
                        double Jrr = domain_geometry_.dFx_dr(r, theta);
                        double Jtr = domain_geometry_.dFy_dr(r, theta);
                        double Jrt = domain_geometry_.dFx_dt(r, theta);
                        double Jtt = domain_geometry_.dFy_dt(r, theta);
                        /* Compute the determinant of the Jacobian matrix */
                        double detDF = Jrr * Jtt - Jrt * Jtr;
                        rhs_f[grid.index(i_r, i_theta)] *= 0.25 * (h1 + h2) * (k1 + k2) * fabs(detDF);
                    }
                    else if (i_r == 0 && DirBC_Interior_) {
                        rhs_f[grid.index(i_r, i_theta)] *= 1.0;
                    }
                    else if (i_r == grid.nr() - 1) {
                        rhs_f[grid.index(i_r, i_theta)] *= 1.0;
                    }
                }
            }

// -------------------------------------------- //
// Discretize rhs values (radial index section) //
// -------------------------------------------- //
#pragma omp for nowait
            for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
                double theta = grid.theta(i_theta);

                for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++) {
                    double r = grid.radius(i_r);
                    if ((0 < i_r && i_r < grid.nr() - 1) || (i_r == 0 && !DirBC_Interior_)) {
                        double h1 = (i_r == 0) ? 2.0 * grid.radius(0) : grid.radialSpacing(i_r - 1);
                        double h2 = grid.radialSpacing(i_r);
                        double k1 = grid.angularSpacing(i_theta - 1);
                        double k2 = grid.angularSpacing(i_theta);
                        /* Calculate the elements of the Jacobian matrix for the transformation mapping */
                        /* The Jacobian matrix is: */
                        /* [Jrr, Jrt] */
                        /* [Jtr, Jtt] */
                        double Jrr = domain_geometry_.dFx_dr(r, theta);
                        double Jtr = domain_geometry_.dFy_dr(r, theta);
                        double Jrt = domain_geometry_.dFx_dt(r, theta);
                        double Jtt = domain_geometry_.dFy_dt(r, theta);
                        /* Compute the determinant of the Jacobian matrix */
                        double detDF = Jrr * Jtt - Jrt * Jtr;
                        rhs_f[grid.index(i_r, i_theta)] *= 0.25 * (h1 + h2) * (k1 + k2) * fabs(detDF);
                    }
                    else if (i_r == 0 && DirBC_Interior_) {
                        rhs_f[grid.index(i_r, i_theta)] *= 1.0;
                    }
                    else if (i_r == grid.nr() - 1) {
                        rhs_f[grid.index(i_r, i_theta)] *= 1.0;
                    }
                }
            }
        }
    }
}
