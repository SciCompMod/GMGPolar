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


// void Residual::applyAGiveCircleSection(const int i_r, Vector<double>& result, const Vector<double>& x, const double& factor) const 
// {    

// }

// void Residual::applyAGiveRadialSection(const int i_theta, Vector<double>& result, const Vector<double>& x, const double& factor) const 
// {
// }


// /* ------------------ */
// /* result = rhs - A*x */
// void Residual::computeResidual(Vector<double>& result, const Vector<double>& rhs, const Vector<double>& x) const {
//     assert(result.size() == x.size());

//     const auto& arr = level_cache_.arr();
//     const auto& att = level_cache_.att();
//     const auto& art = level_cache_.art();
//     const auto& detDF = level_cache_.detDF();
//     const auto& coeff_beta = level_cache_.coeff_beta();

//     auto applyInteriorStencil = [&](int i_r, int i_theta, int center) {
//         const double h1 = grid_.radialSpacing(i_r-1);
//         const double h2 = grid_.radialSpacing(i_r);
//         const double k1 = grid_.angularSpacing(i_theta-1);
//         const double k2 = grid_.angularSpacing(i_theta);

//         const double coeff1 = 0.5 * (k1 + k2) / h1;
//         const double coeff2 = 0.5 * (k1 + k2) / h2;
//         const double coeff3 = 0.5 * (h1 + h2) / k1;
//         const double coeff4 = 0.5 * (h1 + h2) / k2;

//         const int i_theta_M1 = grid_.wrapThetaIndex(i_theta-1);
//         const int i_theta_P1 = grid_.wrapThetaIndex(i_theta+1);

//         const int bottom_left = grid_.index(i_r - 1, i_theta_M1);
//         const int left = grid_.index(i_r - 1, i_theta);
//         const int top_left = grid_.index(i_r - 1, i_theta_P1);
//         const int bottom = grid_.index(i_r, i_theta_M1);
//         const int top = grid_.index(i_r, i_theta_P1);
//         const int bottom_right = grid_.index(i_r + 1, i_theta_M1);
//         const int right = grid_.index(i_r + 1, i_theta);
//         const int top_right = grid_.index(i_r + 1, i_theta_P1);

//         result[center] = rhs[center] - (
//             0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[i_r] * fabs(detDF[center]) * x[center] - 
//             coeff1 * (arr[center] + arr[left]) * (x[left] - x[center]) - 
//             coeff2 * (arr[center] + arr[right]) * (x[right] - x[center]) - 
//             coeff3 * (att[center] + att[bottom]) * (x[bottom] - x[center]) - 
//             coeff4 * (att[center] + att[top]) * (x[top] - x[center]) - 
//             0.25 * (art[left] + art[bottom]) * x[bottom_left] + 
//             0.25 * (art[right] + art[bottom]) * x[bottom_right] + 
//             0.25 * (art[left] + art[top]) * x[top_left] - 
//             0.25 * (art[right] + art[top]) * x[top_right] + 
//         );
//     };

//     auto applyInnerBoundaryCondition = [&](int i_r, int i_theta, int center) {
//         if (DirBC_Interior_) {
//             // Dirichlet boundary on the inner boundary
//             result[center] = rhs[center] - x[center];
//         } else {
//             // Across origin discretization on the interior boundary
//             double h1 = 2.0 * grid_.radius(0);
//             double h2 = grid_.radialSpacing(i_r);
//             double k1 = grid_.angularSpacing(i_theta - 1);
//             double k2 = grid_.angularSpacing(i_theta);

//             double coeff1 = 0.5 * (k1 + k2) / h1;
//             double coeff2 = 0.5 * (k1 + k2) / h2;
//             double coeff3 = 0.5 * (h1 + h2) / k1;
//             double coeff4 = 0.5 * (h1 + h2) / k2;

//             const int i_theta_M1 = grid_.wrapThetaIndex(i_theta - 1);
//             const int i_theta_P1 = grid_.wrapThetaIndex(i_theta + 1);
//             const int i_theta_wrap = grid_.wrapThetaIndex(i_theta + grid_.ntheta() / 2);

//             const int left = grid_.index(i_r, i_theta_wrap);
//             const int bottom = grid_.index(i_r, i_theta_M1);
//             const int top = grid_.index(i_r, i_theta_P1);
//             const int bottom_right = grid_.index(i_r + 1, i_theta_M1);
//             const int right = grid_.index(i_r + 1, i_theta);
//             const int top_right = grid_.index(i_r + 1, i_theta_P1);

//             result[center] = rhs[center] - (
//                 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[i_r] * fabs(detDF[center]) * x[center] - 
//                 coeff1 * (arr[center] + arr[left]) * x[left] - 
//                 coeff2 * (arr[center] + arr[right]) * x[right] - 
//                 coeff3 * (att[center] + att[bottom]) * x[bottom] - 
//                 coeff4 * (att[center] + att[top]) * x[top] + 
//                 0.25 * (art[right] + art[bottom]) * x[bottom_right] - 
//                 0.25 * (art[right] + art[top]) * x[top_right] + 
//                 (coeff1 * (arr[center] + arr[left]) + coeff2 * (arr[center] + arr[right]) + 
//                  coeff3 * (att[center] + att[bottom]) + coeff4 * (att[center] + att[top])) * x[center]
//             );
//         }
//     };

//     auto applyOuterBoundaryCondition = [&](int center) {
//         result[center] = rhs[center] - x[center];
//     };

//     #pragma omp parallel
//     {
//         #pragma omp for nowait
//         for(int i_r = 0; i_r < grid_.numberSmootherCircles(); i_r++) {
//             for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
//                 const int center = grid_.index(i_r, i_theta);
//                 if (i_r > 0 && i_r < grid_.nr() - 1) {
//                     applyInteriorStencil(i_r, i_theta, center);
//                 } else if (i_r == 0) {
//                     applyInnerBoundaryCondition(i_r, i_theta, center);
//                 } else if (i_r == grid_.nr() - 1) {
//                     applyOuterBoundaryCondition(center);
//                 }
//             }
//         }

//         #pragma omp for nowait
//         for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
//             for (int i_r = grid_.numberSmootherCircles(); i_r < grid_.nr(); i_r++) {
//                 const int center = grid_.index(i_r, i_theta);
//                 if (i_r > 0 && i_r < grid_.nr() - 1) {
//                     applyInteriorStencil(i_r, i_theta, center);
//                 } else if (i_r == 0) {
//                     applyInnerBoundaryCondition(i_r, i_theta, center);
//                 } else if (i_r == grid_.nr() - 1) {
//                     applyOuterBoundaryCondition(center);
//                 }
//             }
//         }
//     }
// }





// /* ------------------ */
// /* result = rhs - A*x */
// void Residual::computeResidual(Vector<double>& result, const Vector<double>& rhs, const Vector<double>& x) const {
//     assert(result.size() == x.size());

//     const auto& arr = level_cache_.arr();
//     const auto& att = level_cache_.att();
//     const auto& art = level_cache_.art();
//     const auto& detDF = level_cache_.detDF();
//     const auto& coeff_beta = level_cache_.coeff_beta();

//     auto applyInteriorStencil = [&](int i_r, int i_theta, int center) {
//         const double h1 = grid_.radialSpacing(i_r-1);
//         const double h2 = grid_.radialSpacing(i_r);
//         const double k1 = grid_.angularSpacing(i_theta-1);
//         const double k2 = grid_.angularSpacing(i_theta);

//         const double coeff1 = 0.5 * (k1 + k2) / h1;
//         const double coeff2 = 0.5 * (k1 + k2) / h2;
//         const double coeff3 = 0.5 * (h1 + h2) / k1;
//         const double coeff4 = 0.5 * (h1 + h2) / k2;

//         const int i_theta_M1 = grid_.wrapThetaIndex(i_theta-1);
//         const int i_theta_P1 = grid_.wrapThetaIndex(i_theta+1);

//         const int bottom_left = grid_.index(i_r - 1, i_theta_M1);
//         const int left = grid_.index(i_r - 1, i_theta);
//         const int top_left = grid_.index(i_r - 1, i_theta_P1);
//         const int bottom = grid_.index(i_r, i_theta_M1);
//         const int top = grid_.index(i_r, i_theta_P1);
//         const int bottom_right = grid_.index(i_r + 1, i_theta_M1);
//         const int right = grid_.index(i_r + 1, i_theta);
//         const int top_right = grid_.index(i_r + 1, i_theta_P1);

//         result[center] = rhs[center] - (
//             0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[i_r] * fabs(detDF[center]) * x[center] - 
//             coeff1 * (arr[center] + arr[left]) * x[left] - 
//             coeff2 * (arr[center] + arr[right]) * x[right] - 
//             coeff3 * (att[center] + att[bottom]) * x[bottom] - 
//             coeff4 * (att[center] + att[top]) * x[top] - 
//             0.25 * (art[left] + art[bottom]) * x[bottom_left] + 
//             0.25 * (art[right] + art[bottom]) * x[bottom_right] + 
//             0.25 * (art[left] + art[top]) * x[top_left] - 
//             0.25 * (art[right] + art[top]) * x[top_right] + 
//             (coeff1 * (arr[center] + arr[left]) + coeff2 * (arr[center] + arr[right]) + 
//              coeff3 * (att[center] + att[bottom]) + coeff4 * (att[center] + att[top])) * x[center]
//         );
//     };

//     auto applyInnerBoundaryCondition = [&](int i_r, int i_theta, int center) {
//         if (DirBC_Interior_) {
//             // Dirichlet boundary on the inner boundary
//             result[center] = rhs[center] - x[center];
//         } else {
//             // Across origin discretization on the interior boundary
//             double h1 = 2.0 * grid_.radius(0);
//             double h2 = grid_.radialSpacing(i_r);
//             double k1 = grid_.angularSpacing(i_theta - 1);
//             double k2 = grid_.angularSpacing(i_theta);

//             double coeff1 = 0.5 * (k1 + k2) / h1;
//             double coeff2 = 0.5 * (k1 + k2) / h2;
//             double coeff3 = 0.5 * (h1 + h2) / k1;
//             double coeff4 = 0.5 * (h1 + h2) / k2;

//             const int i_theta_M1 = grid_.wrapThetaIndex(i_theta - 1);
//             const int i_theta_P1 = grid_.wrapThetaIndex(i_theta + 1);
//             const int i_theta_wrap = grid_.wrapThetaIndex(i_theta + grid_.ntheta() / 2);

//             const int left = grid_.index(i_r, i_theta_wrap);
//             const int bottom = grid_.index(i_r, i_theta_M1);
//             const int top = grid_.index(i_r, i_theta_P1);
//             const int bottom_right = grid_.index(i_r + 1, i_theta_M1);
//             const int right = grid_.index(i_r + 1, i_theta);
//             const int top_right = grid_.index(i_r + 1, i_theta_P1);

//             result[center] = rhs[center] - (
//                 0.25 * (h1 + h2) * (k1 + k2) * coeff_beta[i_r] * fabs(detDF[center]) * x[center] - 
//                 coeff1 * (arr[center] + arr[left]) * x[left] - 
//                 coeff2 * (arr[center] + arr[right]) * x[right] - 
//                 coeff3 * (att[center] + att[bottom]) * x[bottom] - 
//                 coeff4 * (att[center] + att[top]) * x[top] + 
//                 0.25 * (art[right] + art[bottom]) * x[bottom_right] - 
//                 0.25 * (art[right] + art[top]) * x[top_right] + 
//                 (coeff1 * (arr[center] + arr[left]) + coeff2 * (arr[center] + arr[right]) + 
//                  coeff3 * (att[center] + att[bottom]) + coeff4 * (att[center] + att[top])) * x[center]
//             );
//         }
//     };

//     auto applyOuterBoundaryCondition = [&](int center) {
//         result[center] = rhs[center] - x[center];
//     };

//     #pragma omp parallel
//     {
//         #pragma omp for nowait
//         for(int i_r = 0; i_r < grid_.numberSmootherCircles(); i_r++) {
//             for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
//                 const int center = grid_.index(i_r, i_theta);
//                 if (i_r > 0 && i_r < grid_.nr() - 1) {
//                     applyInteriorStencil(i_r, i_theta, center);
//                 } else if (i_r == 0) {
//                     applyInnerBoundaryCondition(i_r, i_theta, center);
//                 } else if (i_r == grid_.nr() - 1) {
//                     applyOuterBoundaryCondition(center);
//                 }
//             }
//         }

//         #pragma omp for nowait
//         for (int i_theta = 0; i_theta < grid_.ntheta(); i_theta++) {
//             for (int i_r = grid_.numberSmootherCircles(); i_r < grid_.nr(); i_r++) {
//                 const int center = grid_.index(i_r, i_theta);
//                 if (i_r > 0 && i_r < grid_.nr() - 1) {
//                     applyInteriorStencil(i_r, i_theta, center);
//                 } else if (i_r == 0) {
//                     applyInnerBoundaryCondition(i_r, i_theta, center);
//                 } else if (i_r == grid_.nr() - 1) {
//                     applyOuterBoundaryCondition(center);
//                 }
//             }
//         }
//     }
// }