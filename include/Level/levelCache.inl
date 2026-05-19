#pragma once

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
LevelCache<DomainGeometry, DensityProfileCoefficients>::LevelCache(
    const PolarGrid<DefaultMemorySpace>& grid, const DensityProfileCoefficients& density_profile_coefficients,
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
	    Kokkos::parallel_for(
						"Cache density profile coefficients",
Kokkos::MDRangePolicy<Kokkos::DefaultHostExecutionSpace, Kokkos::Rank<2>>( // Rank of the index space
                {0, 0}, // Starting point of the index space
                {grid.nr(), grid.ntheta()} // Ending point of the index space
                ),
            // Kokkos lambda function to execute for each point in the index space
            KOKKOS_LAMBDA(const int i_r, const int i_theta) {
            const double r = grid.radius(i_r);
                const double theta = grid.theta(i_theta);
                const int index    = grid.index(i_r, i_theta);
                if (!cache_domain_geometry_) {
                    coeff_alpha_(index) = density_profile_coefficients.alpha(r, theta);
                }
                coeff_beta_(index) = density_profile_coefficients.beta(r, theta);
        });
    }

    // Pre-compute and store Jacobian matrix elements (arr, att, art, detDF) at all grid nodes
    // to avoid repeated coordinate transformation calculations during domain operations
    if (cache_domain_geometry_) {
        // We split the loops into two regions to better respect the
        // access patterns of the smoother and improve cache locality

        // Circular Indexing Section
	    Kokkos::parallel_for(
						"Cache domain geometry (circular indexing)",
Kokkos::MDRangePolicy<Kokkos::DefaultHostExecutionSpace, Kokkos::Rank<2>>( // Rank of the index space
                {0, 0}, // Starting point of the index space
                {grid.numberSmootherCircles(), grid.ntheta()} // Ending point of the index space
                ),
            // Kokkos lambda function to execute for each point in the index space
            KOKKOS_LAMBDA(const int i_r, const int i_theta) {
            const double r = grid.radius(i_r);
                const double theta       = grid.theta(i_theta);
                const int index          = grid.index(i_r, i_theta);
                const double coeff_alpha = density_profile_coefficients.alpha(r, theta);

                double arr, att, art, detDF;
                compute_jacobian_elements(domain_geometry_, r, theta, coeff_alpha, arr, att, art, detDF);
                detDF_(index) = detDF;
                arr_(index)   = arr;
                att_(index)   = att;
                art_(index)   = art;
        });
        // Radial Indexing Section
	    Kokkos::parallel_for(
						"Cache domain geometry (radial indexing)",
Kokkos::MDRangePolicy<Kokkos::DefaultHostExecutionSpace, Kokkos::Rank<2>>( // Rank of the index space
                {0, grid.numberSmootherCircles()}, // Starting point of the index space
                {grid.ntheta(), grid.nr()} // Ending point of the index space
                ),
            // Kokkos lambda function to execute for each point in the index space
            KOKKOS_LAMBDA(const int i_theta, const int i_r) {
            const double theta = grid.theta(i_theta);
                const double r           = grid.radius(i_r);
                const int index          = grid.index(i_r, i_theta);
                const double coeff_alpha = density_profile_coefficients.alpha(r, theta);

                double arr, att, art, detDF;
                compute_jacobian_elements(domain_geometry_, r, theta, coeff_alpha, arr, att, art, detDF);
                detDF_(index) = detDF;
                arr_(index)   = arr;
                att_(index)   = att;
                art_(index)   = art;
        });
    }
    Kokkos::fence();
}

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
const DensityProfileCoefficients&
LevelCache<DomainGeometry, DensityProfileCoefficients>::densityProfileCoefficients() const
{
    return density_profile_coefficients_;
}

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
const DomainGeometry& LevelCache<DomainGeometry, DensityProfileCoefficients>::domainGeometry() const
{
    return domain_geometry_;
}

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
bool LevelCache<DomainGeometry, DensityProfileCoefficients>::cacheDensityProfileCoefficients() const
{
    return cache_density_profile_coefficients_;
}

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
ConstVector<double> LevelCache<DomainGeometry, DensityProfileCoefficients>::coeff_alpha() const
{
    return coeff_alpha_;
}

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
ConstVector<double> LevelCache<DomainGeometry, DensityProfileCoefficients>::coeff_beta() const
{
    return coeff_beta_;
}

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
bool LevelCache<DomainGeometry, DensityProfileCoefficients>::cacheDomainGeometry() const
{
    return cache_domain_geometry_;
}

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
ConstVector<double> LevelCache<DomainGeometry, DensityProfileCoefficients>::arr() const
{
    return arr_;
}

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
ConstVector<double> LevelCache<DomainGeometry, DensityProfileCoefficients>::att() const
{
    return att_;
}

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
ConstVector<double> LevelCache<DomainGeometry, DensityProfileCoefficients>::art() const
{
    return art_;
}

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
ConstVector<double> LevelCache<DomainGeometry, DensityProfileCoefficients>::detDF() const
{
    return detDF_;
}
