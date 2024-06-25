#include "../../include/Level/level.h"

LevelCache::LevelCache(const PolarGrid& grid, const SystemParameters& system_parameters, const bool DirBC_Interior) : 
    sin_theta_(grid.ntheta()), 
    cos_theta_(grid.ntheta()),
    rhs_f_(grid.number_of_nodes()),
    u_D_Interior_(DirBC_Interior ? grid.ntheta() : 0),
    u_D_(grid.ntheta())
{
    const int numThreads = omp_get_max_threads();
    omp_set_num_threads(numThreads);

    const int minimalChunkSize = 32;
    TaskDistribution boundary_tasks(grid.ntheta(), minimalChunkSize, numThreads);

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

        // --------------------------------- //
        // Store u_D and u_D_Interior values //
        // --------------------------------- //

        const int threadID = omp_get_thread_num();
        const int local_start = boundary_tasks.getStart(threadID);
        const int local_end = boundary_tasks.getEnd(threadID);

        double R0 = grid.radius(0);
        double Rmax = grid.radius(grid.nr()-1);

        for (int i_theta = local_start; i_theta < local_end; i_theta++){
            double theta = grid.theta(i_theta);
            double sin_theta = sin_theta_[i_theta];
            double cos_theta = cos_theta_[i_theta];
            R0 = grid.radius(0);
            Rmax = grid.radius(grid.nr()-1);
            if (DirBC_Interior) {
                u_D_Interior_[i_theta] = system_parameters.u_D_Interior(R0, theta, sin_theta, cos_theta);
            }
            u_D_[i_theta] = system_parameters.u_D(Rmax, theta, sin_theta, cos_theta);
        }

        // ------------------------------------------- //
        // Store rhs_f values (circular index section) //
        // ------------------------------------------- //
        #pragma omp for nowait
        for (int i_r = 0; i_r < grid.numberSmootherCircles(); i_r++){
            double r = grid.radius(i_r);
            for (int i_theta = 0; i_theta < grid.nr(); i_theta++){
                double theta = grid.theta(i_theta);
                double sin_theta = sin_theta_[i_theta];
                double cos_theta = cos_theta_[i_theta];
                rhs_f_[grid.index(i_r,i_theta)] = system_parameters.rhs_f(r, theta, sin_theta, cos_theta);
            }
        }

        // ----------------------------------------- //
        // Store rhs_f values (radial index section) //
        // ----------------------------------------- //
        #pragma omp for nowait
        for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++){
            double theta = grid.theta(i_theta);
            double sin_theta = sin_theta_[i_theta];
            double cos_theta = cos_theta_[i_theta];
            for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++){
                double r = grid.radius(i_r);
                rhs_f_[grid.index(i_r,i_theta)] = system_parameters.rhs_f(r, theta, sin_theta, cos_theta);
            }
        }
    }    
}



LevelCache::LevelCache(const Level& previous_level, const PolarGrid& current_grid, const bool DirBC_Interior) : 
    sin_theta_(current_grid.ntheta()), 
    cos_theta_(current_grid.ntheta()),
    rhs_f_(current_grid.number_of_nodes()),
    u_D_Interior_(DirBC_Interior ? current_grid.ntheta() : 0),
    u_D_(current_grid.ntheta())
{
    const auto& previous_level_cache = previous_level.levelCache();
    const auto& previous_grid = previous_level.grid();

    const auto& previous_rhs_f = previous_level_cache.rhs_f();
    const auto& previous_u_D_Interior = previous_level_cache.u_D_Interior();
    const auto& previous_u_D = previous_level_cache.u_D();

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
                rhs_f_[current_grid.index(i_r>>1,i_theta>>1)] = previous_rhs_f[previous_grid.index(i_r,i_theta)];
            }
        }

        // ----------------------------------------- //
        // Store rhs_f values (radial index section) //
        // ----------------------------------------- //
        #pragma omp for nowait
        for (int i_theta = 0; i_theta < previous_grid.ntheta(); i_theta += 2){
            for (int i_r = previous_grid.numberSmootherCircles(); i_r < previous_grid.nr(); i_r += 2){
                rhs_f_[current_grid.index(i_r>>1,i_theta>>1)] = previous_rhs_f[previous_grid.index(i_r,i_theta)];
            }
        }

        // --------------------------------- //
        // Store u_D and u_D_Interior values //
        // --------------------------------- //

        const int threadID = omp_get_thread_num();
        const int local_start = boundary_tasks.getStart(threadID);
        const int local_end = boundary_tasks.getEnd(threadID);

        for (int i_theta = local_start; i_theta < local_end; i_theta++){
            if (DirBC_Interior) {
                u_D_Interior_[i_theta] = previous_u_D_Interior[2*i_theta];
            }
            u_D_[i_theta] = previous_u_D[2*i_theta];
        }
    }
}


const std::vector<double>& LevelCache::sin_theta() const {
    return sin_theta_;
}
const std::vector<double>& LevelCache::cos_theta() const {
    return cos_theta_;
}

const Vector<double>& LevelCache::rhs_f() const {
    return rhs_f_;
}
const std::vector<double>& LevelCache::u_D() const {
    return u_D_;
}
const std::vector<double>& LevelCache::u_D_Interior() const {
    return u_D_Interior_;
}
