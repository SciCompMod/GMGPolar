#pragma once

namespace residual_take
{

static inline void node_apply_a_take(const int i_r, const int i_theta, const PolarGrid& grid, bool DirBC_Interior,
                                     Vector<double>& result, ConstVector<double>& x, ConstVector<double>& arr,
                                     ConstVector<double>& att, ConstVector<double>& art, ConstVector<double>& detDF,
                                     ConstVector<double>& coeff_beta)
{
    const int center = grid.index(i_r, i_theta);

    if ((i_r == 0 && DirBC_Interior) || (i_r == grid.nr() - 1)) {
        result[center] = x[center];
        return;
    }

    // Across origin: h1 gets replaced with 2 * R0.
    const double h1 = (i_r == 0) ? 2.0 * grid.radius(0) : grid.radialSpacing(i_r - 1);
    const double h2 = grid.radialSpacing(i_r);
    const double k1 = grid.angularSpacing(i_theta - 1);
    const double k2 = grid.angularSpacing(i_theta);

    const double coeff1 = 0.5 * (k1 + k2) / h1;
    const double coeff2 = 0.5 * (k1 + k2) / h2;
    const double coeff3 = 0.5 * (h1 + h2) / k1;
    const double coeff4 = 0.5 * (h1 + h2) / k2;
    const double coeff5 = 0.25 * (h1 + h2) * (k1 + k2);

    const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
    const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

    // Across origin: (i_r - 1, i_theta) gets replaced with (i_r, i_theta + (grid.ntheta() / 2)).
    const int left =
        (i_r > 0) ? grid.index(i_r - 1, i_theta) : grid.index(i_r, grid.wrapThetaIndex(i_theta + grid.ntheta() / 2));
    const int right  = grid.index(i_r + 1, i_theta);
    const int bottom = grid.index(i_r, i_theta_M1);
    const int top    = grid.index(i_r, i_theta_P1);

    double value = 0.0;

    value += coeff5 * coeff_beta[center] * std::fabs(detDF[center]) * x[center]; /* beta_{i,j} */

    value += coeff1 * (arr[center] + arr[left]) * (x[center] - x[left]); /* Center: (Left) - Left */
    value += coeff2 * (arr[center] + arr[right]) * (x[center] - x[right]); /* Center: (Right) - Right */
    value += coeff3 * (att[center] + att[bottom]) * (x[center] - x[bottom]); /* Center: (Bottom) - Bottom */
    value += coeff4 * (att[center] + att[top]) * (x[center] - x[top]); /* Center: (Top) - Top */

    // Across origin: Reduce 9-point stencil to the artifical 7-point stencil.
    if (i_r > 0) {
        const int bottom_left = grid.index(i_r - 1, i_theta_M1);
        const int top_left    = grid.index(i_r - 1, i_theta_P1);
        value += -0.25 * (art[left] + art[bottom]) * x[bottom_left]; /* Bottom Left */
        value += +0.25 * (art[left] + art[top]) * x[top_left]; /* Top Left */
    }

    const int bottom_right = grid.index(i_r + 1, i_theta_M1);
    const int top_right    = grid.index(i_r + 1, i_theta_P1);
    value += +0.25 * (art[right] + art[bottom]) * x[bottom_right]; /* Bottom Right */
    value += -0.25 * (art[right] + art[top]) * x[top_right]; /* Top Right */

    result[center] = value;
}

} // namespace residual_take

template <class LevelCacheType>
void ResidualTake<LevelCacheType>::applySystemOperator(Vector<double> result, ConstVector<double> x) const
{
    assert(result.size() == x.size());

    using residual_take::node_apply_a_take;

    const PolarGrid& grid     = Residual<LevelCacheType>::grid_;
    const bool DirBC_Interior = Residual<LevelCacheType>::DirBC_Interior_;

    assert(Residual<LevelCacheType>::level_cache_.cacheDensityProfileCoefficients());
    assert(Residual<LevelCacheType>::level_cache_.cacheDomainGeometry());

    ConstVector<double> arr        = Residual<LevelCacheType>::level_cache_.arr();
    ConstVector<double> att        = Residual<LevelCacheType>::level_cache_.att();
    ConstVector<double> art        = Residual<LevelCacheType>::level_cache_.art();
    ConstVector<double> detDF      = Residual<LevelCacheType>::level_cache_.detDF();
    ConstVector<double> coeff_beta = Residual<LevelCacheType>::level_cache_.coeff_beta();

    /* We split the loops into two regions to better respect the */
    /* access patterns of the smoother and improve cache locality. */

    // The For loop matches circular access pattern */
    Kokkos::parallel_for(
        "Residual Take: Apply System Operator (Circular)",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>( // Rank of the index space
            {0, 0}, // Starting point of the index space
            {grid.numberSmootherCircles(), grid.ntheta()} // Ending point of the index space
            ),
        // Kokkos lambda function to execute for each point in the index space
        KOKKOS_LAMBDA(const int i_r, const int i_theta) {
            node_apply_a_take(i_r, i_theta, grid, DirBC_Interior, result, x, arr, att, art, detDF, coeff_beta);
        });

    /* For loop matches radial access pattern */
    Kokkos::parallel_for(
        "Residual Take: Apply System Operator (Radial)",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>( // Rank of the index space
            {0, grid.numberSmootherCircles()}, // Starting point of the index space
            {grid.ntheta(), grid.nr()} // Ending point of the index space
            ),
        // Kokkos lambda function to execute for each point in the index space
        KOKKOS_LAMBDA(const int i_theta, const int i_r) {
            node_apply_a_take(i_r, i_theta, grid, DirBC_Interior, result, x, arr, att, art, detDF, coeff_beta);
        });

    Kokkos::fence();
}