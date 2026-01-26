#include "../../include/Level/level.h"

LevelCache::LevelCache(const PolarGrid& grid, const DensityProfileCoefficients& density_profile_coefficients,
                       const DomainGeometry& domain_geometry, const bool cache_density_profile_coefficients,
                       const bool cache_domain_geometry)
    : domain_geometry_(domain_geometry)
    , density_profile_coefficients_(density_profile_coefficients)
    , cache_density_profile_coefficients_(cache_density_profile_coefficients)
    // If the domain geometry is cached, we don't need to cache the alpha coefficient
    , coeff_alpha_("coeff_alpha",
                   (cache_density_profile_coefficients && !cache_domain_geometry) ? grid.numberOfNodes() : 0)
    , coeff_beta_("coeff_beta", cache_density_profile_coefficients ? grid.numberOfNodes() : 0)
    , cache_domain_geometry_(cache_domain_geometry)
    , arr_("arr", cache_domain_geometry ? grid.numberOfNodes() : 0)
    , att_("att", cache_domain_geometry ? grid.numberOfNodes() : 0)
    , art_("art", cache_domain_geometry ? grid.numberOfNodes() : 0)
    , detDF_("detDF", cache_domain_geometry ? grid.numberOfNodes() : 0)
{
    // Pre-compute and store alpha/beta coefficients at all grid nodes to avoid
    // repeated expensive evaluations during runtime computations
    if (cache_density_profile_coefficients_) {
#pragma omp parallel for
        for (int i_r = 0; i_r < grid.nr(); i_r++) {
            const double r = grid.radius(i_r);
            for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
                const double theta = grid.theta(i_theta);
                const int index    = grid.index(i_r, i_theta);
                if (!cache_domain_geometry_) {
                    coeff_alpha_(index) = density_profile_coefficients.alpha(r, theta);
                }
                coeff_beta_(index) = density_profile_coefficients.beta(r, theta);
            }
        }
    }

    // Pre-compute and store Jacobian matrix elements (arr, att, art, detDF) at all grid nodes
    // to avoid repeated coordinate transformation calculations during domain operations
    if (cache_domain_geometry_) {
        // We split the loops into two regions to better respect the
        // access patterns of the smoother and improve cache locality

        // Circular Indexing Section
#pragma omp parallel for
        for (int i_r = 0; i_r < grid.numberSmootherCircles(); i_r++) {
            const double r = grid.radius(i_r);
            for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
                const double theta       = grid.theta(i_theta);
                const int index          = grid.index(i_r, i_theta);
                const double coeff_alpha = density_profile_coefficients.alpha(r, theta);

                double arr, att, art, detDF;
                compute_jacobian_elements(domain_geometry_, r, theta, coeff_alpha, arr, att, art, detDF);
                detDF_(index) = detDF;
                arr_(index)   = arr;
                att_(index)   = att;
                art_(index)   = art;
            }
        }
        // Radial Indexing Section
#pragma omp parallel for
        for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
            const double theta = grid.theta(i_theta);
            for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++) {
                const double r           = grid.radius(i_r);
                const int index          = grid.index(i_r, i_theta);
                const double coeff_alpha = density_profile_coefficients.alpha(r, theta);

                double arr, att, art, detDF;
                compute_jacobian_elements(domain_geometry_, r, theta, coeff_alpha, arr, att, art, detDF);
                detDF_(index) = detDF;
                arr_(index)   = arr;
                att_(index)   = att;
                art_(index)   = art;
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

bool LevelCache::cacheDensityProfileCoefficients() const
{
    return cache_density_profile_coefficients_;
}
ConstVector<double> LevelCache::coeff_alpha() const
{
    return coeff_alpha_;
}
ConstVector<double> LevelCache::coeff_beta() const
{
    return coeff_beta_;
}

bool LevelCache::cacheDomainGeometry() const
{
    return cache_domain_geometry_;
}
ConstVector<double> LevelCache::arr() const
{
    return arr_;
}
ConstVector<double> LevelCache::att() const
{
    return att_;
}
ConstVector<double> LevelCache::art() const
{
    return art_;
}
ConstVector<double> LevelCache::detDF() const
{
    return detDF_;
}
