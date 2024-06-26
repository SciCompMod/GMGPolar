#include "../../include/Level/level.h"

LevelCache::LevelCache(const PolarGrid& grid, const DomainGeometry& domain_geometry, const SystemParameters& system_parameters, const bool DirBC_Interior) : 
    sin_theta_(grid.ntheta()), 
    cos_theta_(grid.ntheta()),
    rhs_(grid.number_of_nodes())
{
    const int numThreads = omp_get_max_threads();
    omp_set_num_threads(numThreads);

    #pragma omp parallel num_threads(numThreads)
    {   
        // ------------------------------ //
        // Store sinus and cosinus values //
        // ------------------------------ //
        #pragma omp for 
        for (int i = 0; i < grid.ntheta(); i++) {
            sin_theta_[i] = sin(grid.theta(i));
            cos_theta_[i] = cos(grid.theta(i));
        }

        #pragma omp barrier

        // ---------------------------------------------------------- //
        // Store rhs values (circular index section) //
        // ---------------------------------------------------------- //
        #pragma omp for nowait
        for (int i_r = 0; i_r < grid.numberSmootherCircles(); i_r++){
            double r = grid.radius(i_r);
            for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++){
                double theta = grid.theta(i_theta);
                double sin_theta = sin_theta_[i_theta];
                double cos_theta = cos_theta_[i_theta];
                if(0 < i_r && i_r < grid.nr() || (i_r == 0 && !DirBC_Interior)){
                    double h1 = (i_r == 0) ? 2.0 * grid.radius(0) : grid.r_dist(i_r-1);
                    double h2 = grid.r_dist(i_r);
                    double k1 = grid.theta_dist(i_theta-1);
                    double k2 = grid.theta_dist(i_theta);
                    /* Calculate the elements of the Jacobian matrix for the transformation mapping */
                    /* The Jacobian matrix is: */
                    /* [Jrr, Jrt] */
                    /* [Jtr, Jtt] */
                    const double Jrr = domain_geometry.dFx_dr(r, theta, sin_theta, cos_theta);
                    const double Jtr = domain_geometry.dFy_dr(r, theta, sin_theta, cos_theta);
                    const double Jrt = domain_geometry.dFx_dt(r, theta, sin_theta, cos_theta);
                    const double Jtt = domain_geometry.dFy_dt(r, theta, sin_theta, cos_theta);
                    /* Compute the determinant of the Jacobian matrix */
                    const double detDF = Jrr * Jtt - Jrt * Jtr;
                    rhs_[grid.index(i_r,i_theta)] =
                        0.25 * (h1+h2)*(k1+k2) * system_parameters.rhs_f(r, theta, sin_theta, cos_theta) * fabs(detDF);
                }
                else if(i_r == 0 && DirBC_Interior){
                    rhs_[grid.index(i_r,i_theta)] =
                        system_parameters.u_D_Interior(r, theta, sin_theta, cos_theta);
                }
                else if(i_r = grid.nr()){
                    rhs_[grid.index(i_r,i_theta)] =
                        system_parameters.u_D(r, theta, sin_theta, cos_theta);
                }
            }
        }

        // -------------------------------------------------------- //
        // Store rhs values (radial index section) //
        // -------------------------------------------------------- //
        #pragma omp for nowait
        for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++){
            double theta = grid.theta(i_theta);
            double sin_theta = sin_theta_[i_theta];
            double cos_theta = cos_theta_[i_theta];
            for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++){
                double r = grid.radius(i_r);
                if(0 < i_r && i_r < grid.nr() || (i_r == 0 && !DirBC_Interior)){
                    double h1 = (i_r == 0) ? 2.0 * grid.radius(0) : grid.r_dist(i_r-1);
                    double h2 = grid.r_dist(i_r);
                    double k1 = grid.theta_dist(i_theta-1);
                    double k2 = grid.theta_dist(i_theta);
                    /* Calculate the elements of the Jacobian matrix for the transformation mapping */
                    /* The Jacobian matrix is: */
                    /* [Jrr, Jrt] */
                    /* [Jtr, Jtt] */
                    const double Jrr = domain_geometry.dFx_dr(r, theta, sin_theta, cos_theta);
                    const double Jtr = domain_geometry.dFy_dr(r, theta, sin_theta, cos_theta);
                    const double Jrt = domain_geometry.dFx_dt(r, theta, sin_theta, cos_theta);
                    const double Jtt = domain_geometry.dFy_dt(r, theta, sin_theta, cos_theta);
                    /* Compute the determinant of the Jacobian matrix */
                    const double detDF = Jrr * Jtt - Jrt * Jtr;
                    rhs_[grid.index(i_r,i_theta)] =
                        0.25 * (h1+h2)*(k1+k2) * system_parameters.rhs_f(r, theta, sin_theta, cos_theta) * fabs(detDF);
                }
                else if(i_r == 0 && DirBC_Interior){
                    rhs_[grid.index(i_r,i_theta)] =
                        system_parameters.u_D_Interior(r, theta, sin_theta, cos_theta);
                }
                else if(i_r = grid.nr()){
                    rhs_[grid.index(i_r,i_theta)] =
                        system_parameters.u_D(r, theta, sin_theta, cos_theta);
                } 
            }
        }
    }    
}



LevelCache::LevelCache(const Level& previous_level, const PolarGrid& current_grid, const bool DirBC_Interior) : 
    sin_theta_(current_grid.ntheta()), 
    cos_theta_(current_grid.ntheta()),
    rhs_(current_grid.number_of_nodes())
{
    const auto& previous_level_cache = previous_level.levelCache();
    const auto& previous_grid = previous_level.grid();

    const auto& previous_rhs = previous_level_cache.rhs();

    const int numThreads = omp_get_max_threads();
    omp_set_num_threads(numThreads);

    const int minimalChunkSize = 32;
    TaskDistribution boundary_tasks(current_grid.ntheta(), minimalChunkSize, numThreads);

    #pragma omp parallel num_threads(numThreads)
    {   
        // ------------------------------ //
        // Store sinus and cosinus values //
        // ------------------------------ //
        #pragma omp for 
        for (int i = 0; i < current_grid.ntheta(); i++) {
            sin_theta_[i] = sin(current_grid.theta(i));
            cos_theta_[i] = cos(current_grid.theta(i));
        }

        #pragma omp barrier

        // ------------------------------------------- //
        // Store rhs_f values (circular index section) //
        // ----------------------------------------- //
        #pragma omp for nowait
        for (int i_r = 0; i_r < previous_grid.numberSmootherCircles(); i_r += 2){
            for (int i_theta = 0; i_theta < previous_grid.nr(); i_theta += 2){
                rhs_[current_grid.index(i_r>>1,i_theta>>1)] = 
                    previous_rhs[previous_grid.index(i_r,i_theta)];
            }
        }

        // ----------------------------------------- //
        // Store rhs_f values (radial index section) //
        // ----------------------------------------- //
        #pragma omp for nowait
        for (int i_theta = 0; i_theta < previous_grid.ntheta(); i_theta += 2){
            for (int i_r = previous_grid.numberSmootherCircles(); i_r < previous_grid.nr(); i_r += 2){
                rhs_[current_grid.index(i_r>>1,i_theta>>1)] = 
                    previous_rhs[previous_grid.index(i_r,i_theta)];
            }
        }
    }
}


const std::vector<double>& LevelCache::sin_theta() const {
    return sin_theta_;
}
const std::vector<double>& LevelCache::cos_theta() const {
    return cos_theta_;
}

const Vector<double>& LevelCache::rhs() const {
    return rhs_;
}

Vector<double>& LevelCache::rhs() {
    return rhs_;
}
