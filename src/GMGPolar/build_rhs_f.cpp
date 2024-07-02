#include "../../include/GMGPolar/gmgpolar.h"

void GMGPolar::build_rhs_f_discretization(const Level& level, Vector<double>& rhs_f_discretization){
    const PolarGrid& grid = level.grid();
    assert(rhs_f_discretization.size() == grid.number_of_nodes());
    
    const auto& sin_theta_cache = level.levelCache().sin_theta();
    const auto& cos_theta_cache = level.levelCache().cos_theta();

    const int maxThreads = omp_get_max_threads();
    #pragma omp parallel num_threads(maxThreads)
    {
        // ----------------------------------------- //
        // Store rhs values (circular index section) //
        // ----------------------------------------- //
        #pragma omp for nowait
        for (int i_r = 0; i_r < grid.numberSmootherCircles(); i_r++){
            double r = grid.radius(i_r);
            for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++){
                double theta = grid.theta(i_theta);
                double sin_theta = sin_theta_cache[i_theta];
                double cos_theta = cos_theta_cache[i_theta];
                if(0 < i_r && i_r < grid.nr()-1 || (i_r == 0 && !DirBC_Interior)){
                    double h1 = (i_r == 0) ? 2.0 * grid.radius(0) : grid.r_dist(i_r-1);
                    double h2 = grid.r_dist(i_r);
                    double k1 = grid.theta_dist(i_theta-1);
                    double k2 = grid.theta_dist(i_theta);
                    /* Calculate the elements of the Jacobian matrix for the transformation mapping */
                    /* The Jacobian matrix is: */
                    /* [Jrr, Jrt] */
                    /* [Jtr, Jtt] */
                    const double Jrr = domain_geometry_.dFx_dr(r, theta, sin_theta, cos_theta);
                    const double Jtr = domain_geometry_.dFy_dr(r, theta, sin_theta, cos_theta);
                    const double Jrt = domain_geometry_.dFx_dt(r, theta, sin_theta, cos_theta);
                    const double Jtt = domain_geometry_.dFy_dt(r, theta, sin_theta, cos_theta);
                    /* Compute the determinant of the Jacobian matrix */
                    const double detDF = Jrr * Jtt - Jrt * Jtr;
                    rhs_f_discretization[grid.index(i_r,i_theta)] =
                        0.25 * (h1+h2)*(k1+k2) * system_parameters_.rhs_f(r, theta, sin_theta, cos_theta) * fabs(detDF);
                }
                else if(i_r == 0 && DirBC_Interior){
                    rhs_f_discretization[grid.index(i_r,i_theta)] =
                        system_parameters_.u_D_Interior(r, theta, sin_theta, cos_theta);
                }
                else if(i_r = grid.nr()-1){
                    rhs_f_discretization[grid.index(i_r,i_theta)] =
                        system_parameters_.u_D(r, theta, sin_theta, cos_theta);
                }
            }
        }

        // --------------------------------------- //
        // Store rhs values (radial index section) //
        // --------------------------------------- //
        #pragma omp for nowait
        for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++){
            double theta = grid.theta(i_theta);
            double sin_theta = sin_theta_cache[i_theta];
            double cos_theta = cos_theta_cache[i_theta];
            for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++){
                double r = grid.radius(i_r);
                if(0 < i_r && i_r < grid.nr()-1 || (i_r == 0 && !DirBC_Interior)){
                    double h1 = (i_r == 0) ? 2.0 * grid.radius(0) : grid.r_dist(i_r-1);
                    double h2 = grid.r_dist(i_r);
                    double k1 = grid.theta_dist(i_theta-1);
                    double k2 = grid.theta_dist(i_theta);
                    /* Calculate the elements of the Jacobian matrix for the transformation mapping */
                    /* The Jacobian matrix is: */
                    /* [Jrr, Jrt] */
                    /* [Jtr, Jtt] */
                    const double Jrr = domain_geometry_.dFx_dr(r, theta, sin_theta, cos_theta);
                    const double Jtr = domain_geometry_.dFy_dr(r, theta, sin_theta, cos_theta);
                    const double Jrt = domain_geometry_.dFx_dt(r, theta, sin_theta, cos_theta);
                    const double Jtt = domain_geometry_.dFy_dt(r, theta, sin_theta, cos_theta);
                    /* Compute the determinant of the Jacobian matrix */
                    const double detDF = Jrr * Jtt - Jrt * Jtr;
                    rhs_f_discretization[grid.index(i_r,i_theta)] =
                        0.25 * (h1+h2)*(k1+k2) * system_parameters_.rhs_f(r, theta, sin_theta, cos_theta) * fabs(detDF);
                }
                else if(i_r == 0 && DirBC_Interior){
                    rhs_f_discretization[grid.index(i_r,i_theta)] =
                        system_parameters_.u_D_Interior(r, theta, sin_theta, cos_theta);
                }
                else if(i_r = grid.nr()-1){
                    rhs_f_discretization[grid.index(i_r,i_theta)] =
                        system_parameters_.u_D(r, theta, sin_theta, cos_theta);
                } 
            }
        }
    }
}
