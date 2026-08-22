#pragma once

namespace residual_take
{

static KOKKOS_INLINE_FUNCTION void
applySystemMatrixTakeInterior(const int i_r, const int i_theta, const PolarGrid& grid, Vector<double>& result,
                              ConstVector<double>& x, ConstVector<double>& arr, ConstVector<double>& att,
                              ConstVector<double>& art, ConstVector<double>& detDF, ConstVector<double>& coeff_beta)
{
    KOKKOS_ASSERT(0 < i_r && i_r < grid.nr() - 1);

    const int i_theta_M1 = grid.wrapThetaIndex(i_theta - 1);
    const int i_theta_P1 = grid.wrapThetaIndex(i_theta + 1);

    const int bottom_left  = grid.index(i_r - 1, i_theta_M1);
    const int bottom       = grid.index(i_r, i_theta_M1);
    const int bottom_right = grid.index(i_r + 1, i_theta_M1);
    const int left         = grid.index(i_r - 1, i_theta);
    const int center       = grid.index(i_r, i_theta);
    const int right        = grid.index(i_r + 1, i_theta);
    const int top_left     = grid.index(i_r - 1, i_theta_P1);
    const int top          = grid.index(i_r, i_theta_P1);
    const int top_right    = grid.index(i_r + 1, i_theta_P1);

    const double h1 = grid.radialSpacing(i_r - 1);
    const double h2 = grid.radialSpacing(i_r);
    const double k1 = grid.angularSpacing(i_theta_M1);
    const double k2 = grid.angularSpacing(i_theta);

    const double coeff1 = 0.5 * (k1 + k2) / h1;
    const double coeff2 = 0.5 * (k1 + k2) / h2;
    const double coeff3 = 0.5 * (h1 + h2) / k1;
    const double coeff4 = 0.5 * (h1 + h2) / k2;
    const double coeff5 = 0.25 * (h1 + h2) * (k1 + k2);

    result[center] = (+coeff5 * coeff_beta[center] * Kokkos::fabs(detDF[center]) * x[center] /* beta_{i,j} */

                      + coeff1 * (arr[center] + arr[left]) * (x[center] - x[left]) /* Center: (Left) - Left */
                      + coeff2 * (arr[center] + arr[right]) * (x[center] - x[right]) /* Center: (Right) - Right */
                      + coeff3 * (att[center] + att[bottom]) * (x[center] - x[bottom]) /* Center: (Bottom) - Bottom */
                      + coeff4 * (att[center] + att[top]) * (x[center] - x[top]) /* Center: (Top) - Top */

                      - (art[left] + art[bottom]) * x[bottom_left] /* Bottom Left */
                      + (art[left] + art[top]) * x[top_left] /* Top Left */
                      + (art[right] + art[bottom]) * x[bottom_right] /* Bottom Right */
                      - (art[right] + art[top]) * x[top_right] /* Top Right */
    );
}

static KOKKOS_INLINE_FUNCTION void applySystemMatrixTakeBoundary(const int i_r, const int i_theta,
                                                                 const PolarGrid& grid, bool DirBC_Interior,
                                                                 Vector<double>& result, ConstVector<double>& x,
                                                                 ConstVector<double>& arr, ConstVector<double>& att,
                                                                 ConstVector<double>& art, ConstVector<double>& detDF,
                                                                 ConstVector<double>& coeff_beta)
{
    KOKKOS_ASSERT(i_r == 0 || i_r == grid.nr() - 1);

    if ((i_r == 0 && DirBC_Interior) || (i_r == grid.nr() - 1)) {
        const int center = grid.index(i_r, i_theta);
        result[center]   = x[center];
        return;
    }

    const int i_theta_M1     = grid.wrapThetaIndex(i_theta - 1);
    const int i_theta_P1     = grid.wrapThetaIndex(i_theta + 1);
    const int i_theta_Across = grid.wrapThetaIndex(i_theta + grid.ntheta() / 2);

    const int left         = grid.index(i_r, i_theta_Across);
    const int bottom       = grid.index(i_r, i_theta_M1);
    const int center       = grid.index(i_r, i_theta);
    const int top          = grid.index(i_r, i_theta_P1);
    const int bottom_right = grid.index(i_r + 1, i_theta_M1);
    const int right        = grid.index(i_r + 1, i_theta);
    const int top_right    = grid.index(i_r + 1, i_theta_P1);

    // Across origin: h1 gets replaced with 2 * R0.
    const double h1 = 2.0 * grid.radius(0);
    const double h2 = grid.radialSpacing(i_r);
    const double k1 = grid.angularSpacing(i_theta - 1);
    const double k2 = grid.angularSpacing(i_theta);

    const double coeff1 = 0.5 * (k1 + k2) / h1;
    const double coeff2 = 0.5 * (k1 + k2) / h2;
    const double coeff3 = 0.5 * (h1 + h2) / k1;
    const double coeff4 = 0.5 * (h1 + h2) / k2;
    const double coeff5 = 0.25 * (h1 + h2) * (k1 + k2);

    result[center] = (+coeff5 * coeff_beta[center] * Kokkos::fabs(detDF[center]) * x[center] /* beta_{i,j} */

                      + coeff1 * (arr[center] + arr[left]) * (x[center] - x[left]) /* Center: (Left) - Left */
                      + coeff2 * (arr[center] + arr[right]) * (x[center] - x[right]) /* Center: (Right) - Right */
                      + coeff3 * (att[center] + att[bottom]) * (x[center] - x[bottom]) /* Center: (Bottom) - Bottom */
                      + coeff4 * (att[center] + att[top]) * (x[center] - x[top]) /* Center: (Top) - Top */

                      + (art[right] + art[bottom]) * x[bottom_right] /* Bottom Right */
                      - (art[right] + art[top]) * x[top_right] /* Top Right */
    );
}

} // namespace residual_take

template <class LevelCacheType>
void ResidualTake<LevelCacheType>::applySystemOperator(Vector<double> result, ConstVector<double> x) const
{
    using residual_take::applySystemMatrixTakeBoundary;
    using residual_take::applySystemMatrixTakeInterior;

    assert(result.size() == x.size());

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
        "Residual Take: Apply System Operator Boundary (Circular)",
        Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>( // Rank of the index space
            {0, 0}, // Starting point of the index space
            {1, grid.ntheta()} // Ending point of the index space
            ),
        // Kokkos lambda function to execute for each point in the index space
        KOKKOS_LAMBDA(const int i_r, const int i_theta) {
            applySystemMatrixTakeBoundary(i_r, i_theta, grid, DirBC_Interior, result, x, arr, att, art, detDF,
                                          coeff_beta);
        });

    Kokkos::parallel_for(
        "Residual Take: Apply System Operator Interior (Circular)",
        Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>( // Rank of the index space
            {1, 0}, // Starting point of the index space
            {grid.numberSmootherCircles(), grid.ntheta()} // Ending point of the index space
            ),
        // Kokkos lambda function to execute for each point in the index space
        KOKKOS_LAMBDA(const int i_r, const int i_theta) {
            applySystemMatrixTakeInterior(i_r, i_theta, grid, result, x, arr, att, art, detDF, coeff_beta);
        });

    /* For loop matches radial access pattern */
    Kokkos::parallel_for(
        "Residual Take: Apply System Operator Interior (Radial)",
        Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>( // Rank of the index space
            {0, grid.numberSmootherCircles()}, // Starting point of the index space
            {grid.ntheta(), grid.nr() - 1} // Ending point of the index space
            ),
        // Kokkos lambda function to execute for each point in the index space
        KOKKOS_LAMBDA(const int i_theta, const int i_r) {
            applySystemMatrixTakeInterior(i_r, i_theta, grid, result, x, arr, att, art, detDF, coeff_beta);
        });

    Kokkos::parallel_for(
        "Residual Take: Apply System Operator Boundary (Radial)",
        Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>( // Rank of the index space
            {0, grid.nr() - 1}, // Starting point of the index space
            {grid.ntheta(), grid.nr()} // Ending point of the index space
            ),
        // Kokkos lambda function to execute for each point in the index space
        KOKKOS_LAMBDA(const int i_theta, const int i_r) {
            applySystemMatrixTakeBoundary(i_r, i_theta, grid, DirBC_Interior, result, x, arr, att, art, detDF,
                                          coeff_beta);
        });

    Kokkos::fence();
}
