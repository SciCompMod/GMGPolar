#include "../../include/Level/level.h"

LevelCache::LevelCache(const PolarGrid& grid) : 
    sin_theta_(grid.ntheta()), 
    cos_theta_(grid.ntheta())
{
    const int maxThreads = omp_get_max_threads();
    const int minimalChunkSize = 64;
    TaskDistribution sin_cos_partition(grid.ntheta(), minimalChunkSize, maxThreads);

    #pragma omp parallel num_threads(maxThreads)
    {   
        // ------------------------------ //
        // Store sinus and cosinus values //
        // ------------------------------ //
        const int threadID = omp_get_thread_num();
        const int local_start = sin_cos_partition.getStart(threadID);
        const int local_end = sin_cos_partition.getEnd(threadID);
        for (int i = local_start; i < local_end; ++i) {
            const double theta = grid.theta(i);
            sin_theta_[i] = sin(theta);
            cos_theta_[i] = cos(theta);
        }
    }
}



LevelCache::LevelCache(const Level& previous_level, const PolarGrid& current_grid) : 
    sin_theta_(current_grid.ntheta()), 
    cos_theta_(current_grid.ntheta())
{
    /* const auto& previous_level_cache = previous_level.levelCache(); */
    /* const auto& previous_grid = previous_level.grid(); */

    const int maxThreads = omp_get_max_threads();
    const int minimalChunkSize = 64;
    TaskDistribution sin_cos_partition(current_grid.ntheta(), minimalChunkSize, maxThreads);

    #pragma omp parallel num_threads(maxThreads)
    {   
        // ------------------------------ //
        // Store sinus and cosinus values //
        // ------------------------------ //
        const int threadID = omp_get_thread_num();
        const int local_start = sin_cos_partition.getStart(threadID);
        const int local_end = sin_cos_partition.getEnd(threadID);
        for (int i = local_start; i < local_end; ++i) {
            const double theta = current_grid.theta(i);
            sin_theta_[i] = sin(theta);
            cos_theta_[i] = cos(theta);
        }
    }
}


const std::vector<double>& LevelCache::sin_theta() const {
    return sin_theta_;
}
const std::vector<double>& LevelCache::cos_theta() const {
    return cos_theta_;
}
