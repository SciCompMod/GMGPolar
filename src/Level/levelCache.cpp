#include "../../include/Level/level.h"

LevelCache::LevelCache(const PolarGrid& grid, const DensityProfileCoefficients& density_profile_coefficients,
                       const DomainGeometry& domain_geometry, const bool cache_density_profile_coefficients,
                       const bool cache_domain_geometry)
    : domain_geometry_(domain_geometry)
    , density_profile_coefficients_(density_profile_coefficients)
    , sin_theta_(grid.ntheta())
    , cos_theta_(grid.ntheta())
    , cache_density_profile_coefficients_(cache_density_profile_coefficients)
    // If the domain geometry is cached, we don't need to cache the alpha coefficient
    , coeff_alpha_((cache_density_profile_coefficients && !cache_domain_geometry) ? grid.numberOfNodes() : 0)
    , coeff_beta_(cache_density_profile_coefficients ? grid.numberOfNodes() : 0)
    , cache_domain_geometry_(cache_domain_geometry)
    , arr_(cache_domain_geometry ? grid.numberOfNodes() : 0)
    , att_(cache_domain_geometry ? grid.numberOfNodes() : 0)
    , art_(cache_domain_geometry ? grid.numberOfNodes() : 0)
    , detDF_(cache_domain_geometry ? grid.numberOfNodes() : 0)
{
#pragma omp parallel for
    for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
        const double theta  = grid.theta(i_theta);
        sin_theta_[i_theta] = sin(theta);
        cos_theta_[i_theta] = cos(theta);
    }

    if (cache_density_profile_coefficients_) {
#pragma omp parallel for
        for (int i_r = 0; i_r < grid.nr(); i_r++) {
            const double r = grid.radius(i_r);
            for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
                const double theta = grid.theta(i_theta);
                const double index = grid.index(i_r, i_theta);
                if (!cache_domain_geometry_) {
                    coeff_alpha_[index] = density_profile_coefficients.alpha(r, theta);
                }
                coeff_beta_[index] = density_profile_coefficients.beta(r, theta);
            }
        }
    }

    if (cache_domain_geometry_) {
#pragma omp parallel for
        for (int i_r = 0; i_r < grid.numberSmootherCircles(); i_r++) {
            const double r = grid.radius(i_r);
            for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
                const double theta     = grid.theta(i_theta);
                const double sin_theta = sin_theta_[i_theta];
                const double cos_theta = cos_theta_[i_theta];
                const double index     = grid.index(i_r, i_theta);

                double coeff_alpha = density_profile_coefficients.alpha(r, theta);

                double arr, att, art, detDF;
                compute_jacobian_elements(domain_geometry_, r, theta, sin_theta, cos_theta, coeff_alpha, arr, att, art,
                                          detDF);
                detDF_[index] = detDF;
                arr_[index]   = arr;
                att_[index]   = att;
                art_[index]   = art;
            }
        }

#pragma omp parallel for
        for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
            const double theta     = grid.theta(i_theta);
            const double sin_theta = sin_theta_[i_theta];
            const double cos_theta = cos_theta_[i_theta];
            for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++) {
                const double r     = grid.radius(i_r);
                const double index = grid.index(i_r, i_theta);

                double coeff_alpha;
                if (cache_density_profile_coefficients_ && !cache_domain_geometry_) {
                    coeff_alpha = coeff_alpha_[index];
                }
                else {
                    coeff_alpha = density_profile_coefficients.alpha(r, theta);
                }

                double arr, att, art, detDF;
                compute_jacobian_elements(domain_geometry_, r, theta, sin_theta, cos_theta, coeff_alpha, arr, att, art,
                                          detDF);
                detDF_[index] = detDF;
                arr_[index]   = arr;
                att_[index]   = att;
                art_[index]   = art;
            }
        }
    }
}

LevelCache::LevelCache(const Level& previous_level, const PolarGrid& current_grid)
    : domain_geometry_(previous_level.levelCache().domainGeometry())
    , density_profile_coefficients_(previous_level.levelCache().densityProfileCoefficients())
    , sin_theta_(current_grid.ntheta())
    , cos_theta_(current_grid.ntheta())
    , cache_density_profile_coefficients_(previous_level.levelCache().cacheDensityProfileCoefficients())
    , coeff_alpha_(previous_level.levelCache().coeff_alpha().size() > 0 ? current_grid.numberOfNodes() : 0)
    , coeff_beta_(previous_level.levelCache().coeff_beta().size() > 0 ? current_grid.numberOfNodes() : 0)
    , cache_domain_geometry_(previous_level.levelCache().cacheDomainGeometry())
    , arr_(previous_level.levelCache().arr().size() > 0 ? current_grid.numberOfNodes() : 0)
    , att_(previous_level.levelCache().att().size() > 0 ? current_grid.numberOfNodes() : 0)
    , art_(previous_level.levelCache().art().size() > 0 ? current_grid.numberOfNodes() : 0)
    , detDF_(previous_level.levelCache().detDF().size() > 0 ? current_grid.numberOfNodes() : 0)
{
    const auto& previous_level_cache = previous_level.levelCache();

    for (int i_theta = 0; i_theta < current_grid.ntheta(); i_theta++) {
        const double theta  = current_grid.theta(i_theta);
        sin_theta_[i_theta] = previous_level_cache.sin_theta()[2 * i_theta];
        cos_theta_[i_theta] = previous_level_cache.cos_theta()[2 * i_theta];
    }

    if (previous_level_cache.cacheDensityProfileCoefficients()) {
#pragma omp parallel for
        for (int i_r = 0; i_r < current_grid.nr(); i_r++) {
            for (int i_theta = 0; i_theta < current_grid.ntheta(); i_theta++) {
                const int current_index  = current_grid.index(i_r, i_theta);
                const int previous_index = previous_level.grid().index(2 * i_r, 2 * i_theta);

                if (!previous_level_cache.cacheDomainGeometry()) {
                    coeff_alpha_[current_index] = previous_level_cache.coeff_alpha()[previous_index];
                }
                coeff_beta_[current_index] = previous_level_cache.coeff_beta()[previous_index];
            }
        }
    }

    if (previous_level_cache.cacheDomainGeometry()) {
#pragma omp parallel for
        for (int i_r = 0; i_r < current_grid.numberSmootherCircles(); i_r++) {
            for (int i_theta = 0; i_theta < current_grid.ntheta(); i_theta++) {
                const int current_index  = current_grid.index(i_r, i_theta);
                const int previous_index = previous_level.grid().index(2 * i_r, 2 * i_theta);
                arr_[current_index]      = previous_level_cache.arr()[previous_index];
                att_[current_index]      = previous_level_cache.att()[previous_index];
                art_[current_index]      = previous_level_cache.art()[previous_index];
                detDF_[current_index]    = previous_level_cache.detDF()[previous_index];
            }
        }
#pragma omp parallel for
        for (int i_theta = 0; i_theta < current_grid.ntheta(); i_theta++) {
            for (int i_r = current_grid.numberSmootherCircles(); i_r < current_grid.nr(); i_r++) {
                const int current_index  = current_grid.index(i_r, i_theta);
                const int previous_index = previous_level.grid().index(2 * i_r, 2 * i_theta);
                arr_[current_index]      = previous_level_cache.arr()[previous_index];
                att_[current_index]      = previous_level_cache.att()[previous_index];
                art_[current_index]      = previous_level_cache.art()[previous_index];
                detDF_[current_index]    = previous_level_cache.detDF()[previous_index];
            }
        }
    }
}

const DensityProfileCoefficients& LevelCache::densityProfileCoefficients() const
{
    return density_profile_coefficients_;
}
const DomainGeometry& LevelCache::domainGeometry() const
{
    return domain_geometry_;
}

const std::vector<double>& LevelCache::sin_theta() const
{
    return sin_theta_;
}
const std::vector<double>& LevelCache::cos_theta() const
{
    return cos_theta_;
}

bool LevelCache::cacheDensityProfileCoefficients() const
{
    return cache_density_profile_coefficients_;
}
const Vector<double>& LevelCache::coeff_alpha() const
{
    return coeff_alpha_;
}
const Vector<double>& LevelCache::coeff_beta() const
{
    return coeff_beta_;
}

bool LevelCache::cacheDomainGeometry() const
{
    return cache_domain_geometry_;
}
const Vector<double>& LevelCache::arr() const
{
    return arr_;
}
const Vector<double>& LevelCache::att() const
{
    return att_;
}
const Vector<double>& LevelCache::art() const
{
    return art_;
}
const Vector<double>& LevelCache::detDF() const
{
    return detDF_;
}