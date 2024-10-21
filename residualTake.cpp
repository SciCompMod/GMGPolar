// #include "../../include/Residual/residual.h"

// Residual::Residual(
//     const PolarGrid& grid, const LevelCache& level_cache, 
//     const DomainGeometry& domain_geometry, const DensityProfileCoefficients& density_profile_coefficients,
//     bool DirBC_Interior, int num_omp_threads
// ) :
//     grid_(grid),
//     level_cache_(level_cache),
//     domain_geometry_(domain_geometry),
//     density_profile_coefficients_(density_profile_coefficients),
//     DirBC_Interior_(DirBC_Interior),
//     num_omp_threads_(num_omp_threads)
// {}


#define NODE_APPLY_A_TAKE(i_r, i_theta, \
    grid, DirBC_Interior, \
    result, x, factor, \
    arr, att, art, detDF, coeff_beta) \
do { \
    /* -------------------- */ \
    /* Node in the interior */ \
    /* -------------------- */ \
    if (i_r > 0 && i_r < grid.nr() - 1) { \
        double h1 = grid.radialSpacing(i_r-1); \
        double h2 = grid.radialSpacing(i_r); \
        double k1 = grid.angularSpacing(i_theta-1); \
        double k2 = grid.angularSpacing(i_theta); \
        \
        double coeff1 = 0.5*(k1+k2)/h1; \
        double coeff2 = 0.5*(k1+k2)/h2; \
        double coeff3 = 0.5*(h1+h2)/k1; \
        double coeff4 = 0.5*(h1+h2)/k2; \
        \
        result[grid.index(i_r,i_theta)] += factor * ( \
            0.25 * (h1+h2)*(k1+k2) * coeff_beta[i_r] * fabs(detDF[grid.index(i_r,i_theta)]) * x[grid.index(i_r,i_theta)] /* beta_{i,j} */ \
            \
            - coeff1 * (arr[grid.index(i_r,i_theta)] + arr[grid.index(i_r-1,i_theta)]) * x[grid.index(i_r-1,i_theta)] /* Left */ \
            - coeff2 * (arr[grid.index(i_r,i_theta)] + arr[grid.index(i_r+1,i_theta)]) * x[grid.index(i_r+1,i_theta)] /* Right */ \
            - coeff3 * (att[grid.index(i_r,i_theta)] + att[grid.index(i_r,i_theta-1)]) * x[grid.index(i_r,i_theta-1)] /* Bottom */ \
            - coeff4 * (att[grid.index(i_r,i_theta)] + att[grid.index(i_r,i_theta+1)]) * x[grid.index(i_r,i_theta+1)] /* Top */ \
            \
            - 0.25 * (art[grid.index(i_r-1,i_theta)] + art[grid.index(i_r,i_theta-1)]) * x[grid.index(i_r-1,i_theta-1)] /* Bottom Left */ \
            + 0.25 * (art[grid.index(i_r+1,i_theta)] + art[grid.index(i_r,i_theta-1)]) * x[grid.index(i_r+1,i_theta-1)] /* Bottom Right */ \
            + 0.25 * (art[grid.index(i_r-1,i_theta)] + art[grid.index(i_r,i_theta+1)]) * x[grid.index(i_r-1,i_theta+1)] /* Top Left */ \
            - 0.25 * (art[grid.index(i_r+1,i_theta)] + art[grid.index(i_r,i_theta+1)]) * x[grid.index(i_r+1,i_theta+1)] /* Top Right */ \
            \
            + ( \
                + coeff1 * (arr[grid.index(i_r,i_theta)] + arr[grid.index(i_r-1,i_theta)]) /* Center: (Left) */ \
                + coeff2 * (arr[grid.index(i_r,i_theta)] + arr[grid.index(i_r+1,i_theta)]) /* Center: (Right) */ \
                + coeff3 * (att[grid.index(i_r,i_theta)] + att[grid.index(i_r,i_theta-1)]) /* Center: (Bottom) */ \
                + coeff4 * (att[grid.index(i_r,i_theta)] + att[grid.index(i_r,i_theta+1)]) /* Center: (Top) */ \
            ) * x[grid.index(i_r,i_theta)] \
        ); \
    /* -------------------------- */ \
    /* Node on the inner boundary */ \
    /* -------------------------- */ \
    } else if (i_r == 0) { \
        /* ------------------------------------------------ */ \
        /* Case 1: Dirichlet boundary on the inner boundary */ \
        /* ------------------------------------------------ */ \
        if(DirBC_Interior){ \
            result[grid.index(i_r,i_theta)] += factor * x[grid.index(i_r,i_theta)]; \
        } else{ \
            /* ------------------------------------------------------------- */ \
            /* Case 2: Across origin discretization on the interior boundary */ \
            /* ------------------------------------------------------------- */ \
            /* h1 gets replaced with 2 * R0. */ \
            /* (i_r-1,i_theta) gets replaced with (i_r, i_theta + (grid.ntheta()>>1)). */ \
            /* Some more adjustments from the changing the 9-point stencil to the artifical 7-point stencil. */ \
            double h1 = 2.0 * grid.radius(0); \
            double h2 = grid.radialSpacing(i_r); \
            double k1 = grid.angularSpacing(i_theta-1); \
            double k2 = grid.angularSpacing(i_theta); \
            double coeff1 = 0.5*(k1+k2)/h1; \
            double coeff2 = 0.5*(k1+k2)/h2; \
            double coeff3 = 0.5*(h1+h2)/k1; \
            double coeff4 = 0.5*(h1+h2)/k2; \
            \
            result[grid.index(i_r,i_theta)] += factor * ( \
                0.25 * (h1+h2)*(k1+k2) * coeff_beta[i_r] * fabs(detDF[grid.index(i_r,i_theta)]) * x[grid.index(i_r,i_theta)] /* beta_{i,j} */ \
                \
                - coeff1 * (arr[grid.index(i_r,i_theta)] + arr[grid.index(i_r-1,i_theta)]) * (x[grid.index(i_r-1,i_theta)] - x[grid.index(i_r,i_theta)]) /* Left - Center: (Left) */ \
                - coeff2 * (arr[grid.index(i_r,i_theta)] + arr[grid.index(i_r+1,i_theta)]) * (x[grid.index(i_r+1,i_theta)] - x[grid.index(i_r,i_theta)]) /* Right - Center: (Right) */ \
                - coeff3 * (att[grid.index(i_r,i_theta)] + att[grid.index(i_r,i_theta-1)]) * (x[grid.index(i_r,i_theta-1)] - x[grid.index(i_r,i_theta)]) /* Bottom - Center: (Bottom) */ \
                - coeff4 * (att[grid.index(i_r,i_theta)] + att[grid.index(i_r,i_theta+1)]) * (x[grid.index(i_r,i_theta+1)] - x[grid.index(i_r,i_theta)]) /* Top - Center: (Top) */ \
                \
                - 0.25 * (art[grid.index(i_r-1,i_theta)] + art[grid.index(i_r,i_theta-1)]) * x[grid.index(i_r-1,i_theta-1)] /* Bottom Left */ \
                + 0.25 * (art[grid.index(i_r+1,i_theta)] + art[grid.index(i_r,i_theta-1)]) * x[grid.index(i_r+1,i_theta-1)] /* Bottom Right */ \
                + 0.25 * (art[grid.index(i_r-1,i_theta)] + art[grid.index(i_r,i_theta+1)]) * x[grid.index(i_r-1,i_theta+1)] /* Top Left */ \
                - 0.25 * (art[grid.index(i_r+1,i_theta)] + art[grid.index(i_r,i_theta+1)]) * x[grid.index(i_r+1,i_theta+1)] /* Top Right */ \
            ); \
        } \
    /* ----------------------------- */ \
    /* Node on to the outer boundary */ \
    /* ----------------------------- */ \
    } else if (i_r == grid.nr() - 1) { \
        /* Fill result of (i,j) */ \
        result[grid.index(i_r,i_theta)] += factor * x[grid.index(i_r,i_theta)]; \
    } \
} while(0)


// void Residual::applyAGiveCircleSection(const int i_r, Vector<double>& result, const Vector<double>& x, const double& factor) const 
// {    
//     const auto& arr = level_cache_.arr();
//     const auto& att = level_cache_.att();
//     const auto& art = level_cache_.art();
//     const auto& detDF = level_cache_.detDF();
//     const auto& coeff_beta = level_cache_.coeff_beta();

//     for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++){
//         NODE_APPLY_A_TAKE(i_r, i_theta,
//             grid_, DirBC_Interior_,
//             result, x, factor,
//             arr, att, art, detDF, coeff_beta);
//     }
// }

// void Residual::applyAGiveRadialSection(const int i_theta, Vector<double>& result, const Vector<double>& x, const double& factor) const 
// {
//     const auto& arr = level_cache_.arr();
//     const auto& att = level_cache_.att();
//     const auto& art = level_cache_.art();
//     const auto& detDF = level_cache_.detDF();
//     const auto& coeff_beta = level_cache_.coeff_beta();

//     for (int i_r = grid_.numberSmootherCircles(); i_r < grid_.nr(); i_r++){
//         NODE_APPLY_A_TAKE(i_r, i_theta,
//             grid_, DirBC_Interior_,
//             result, x, factor,
//             arr, att, art, detDF, coeff_beta);
//     }
// }


// /* ------------------ */
// /* result = rhs - A*x */
// void Residual::computeResidual(Vector<double>& result, const Vector<double>& rhs, const Vector<double>& x) const {
//     assert(result.size() == x.size());

//     bool use_simple_parallelism = true; // Fastest: true

//     omp_set_num_threads(num_omp_threads_);

//     result = rhs;

//     const double factor = -1.0;

//     if(omp_get_max_threads() == 1){
//         /* Single-threaded execution */
//         for(int i_r = 0; i_r < grid_.numberSmootherCircles(); i_r++) {
//             applyAGiveCircleSection(i_r, result, x, factor);
//         }
//         for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
//             applyAGiveRadialSection(i_theta, result, x, factor);
//         }
//     }
//     else{

//         #pragma omp parallel
//         {
//             #pragma omp for nowait
//             for(int i_r = 0; i_r < grid_.numberSmootherCircles(); i_r++) {
//                 applyAGiveCircleSection(i_r, result, x, factor);
//             }
//             #pragma omp for nowait
//             for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
//                 applyAGiveRadialSection(i_theta, result, x, factor);
//             }
//         }
//     }
// }