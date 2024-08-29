#include "../../include/Level/level.h"

LevelCache::LevelCache(const PolarGrid& grid, const DensityProfileCoefficients& density_profile_coefficients) : 
    sin_theta_(grid.ntheta()), 
    cos_theta_(grid.ntheta()),
    coeff_alpha_(grid.nr()),
    coeff_beta_(grid.nr())
{
    for (std::size_t i_theta = 0; i_theta < grid.ntheta(); i_theta++){
        const double theta = grid.theta(i_theta);
        sin_theta_[i_theta] = sin(theta);
        cos_theta_[i_theta] = cos(theta);
    }
    for (std::size_t i_r = 0; i_r < grid.nr(); i_r++){
        const double r = grid.radius(i_r);
        coeff_alpha_[i_r] = density_profile_coefficients.alpha(r);
        coeff_beta_[i_r] = density_profile_coefficients.beta(r);
    }
}

LevelCache::LevelCache(const Level& previous_level, const PolarGrid& current_grid) : 
    sin_theta_(current_grid.ntheta()), 
    cos_theta_(current_grid.ntheta()),
    coeff_alpha_(current_grid.nr()),
    coeff_beta_(current_grid.nr())
{
    const auto& previous_level_cache = previous_level.levelCache();

    for (std::size_t i_theta = 0; i_theta < current_grid.ntheta(); i_theta++){
        const double theta = current_grid.theta(i_theta);
        sin_theta_[i_theta] = sin(theta);
        cos_theta_[i_theta] = cos(theta);
    }

    for (std::size_t i_r = 0; i_r < current_grid.nr(); i_r++){
        coeff_alpha_[i_r] = previous_level_cache.coeff_alpha()[2*i_r];
        coeff_beta_[i_r] = previous_level_cache.coeff_beta()[2*i_r];
    }
}

const std::vector<double>& LevelCache::sin_theta() const {
    return sin_theta_;
}
const std::vector<double>& LevelCache::cos_theta() const {
    return cos_theta_;
}

const std::vector<double>& LevelCache::coeff_alpha() const {
    return coeff_alpha_;
}
const std::vector<double>& LevelCache::coeff_beta() const {
    return coeff_beta_;
}
