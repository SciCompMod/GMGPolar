#include "../../include/Level/level.h"

LevelCache::LevelCache(
    const ProcessingType processing_type,
    const PolarGrid& grid,
    const DensityProfileCoefficients& density_profile_coefficients,
    const DomainGeometry& domain_geometry
) : 
    processing_type_(processing_type),

    sin_theta_(processing_type != ProcessingType::GPU ? grid.ntheta() : 0), 
    cos_theta_(processing_type != ProcessingType::GPU ? grid.ntheta() : 0),

    coeff_alpha_(processing_type != ProcessingType::GPU ? grid.nr() : 0),
    coeff_beta_(processing_type != ProcessingType::GPU ? grid.nr() : 0),

    arr_(processing_type != ProcessingType::GPU ? grid.numberOfNodes() : 0),
    att_(processing_type != ProcessingType::GPU ? grid.numberOfNodes() : 0),
    art_(processing_type != ProcessingType::GPU ? grid.numberOfNodes() : 0),
    detDF_(processing_type != ProcessingType::GPU ? grid.numberOfNodes() : 0),

    gpu_sin_theta_(processing_type == ProcessingType::GPU ? grid.ntheta() : 0), 
    gpu_cos_theta_(processing_type == ProcessingType::GPU ? grid.ntheta() : 0),

    gpu_coeff_alpha_(processing_type == ProcessingType::GPU ? grid.nr() : 0), 
    gpu_coeff_beta_(processing_type == ProcessingType::GPU ? grid.nr() : 0)
{
    if(processing_type_ != ProcessingType::GPU){
        #pragma omp parallel for
        for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++){
            const double theta = grid.theta(i_theta);
            sin_theta_[i_theta] = sin(theta);
            cos_theta_[i_theta] = cos(theta);
        }
        #pragma omp parallel for
        for (int i_r = 0; i_r < grid.nr(); i_r++){
            const double r = grid.radius(i_r);
            coeff_alpha_[i_r] = density_profile_coefficients.alpha(r);
            coeff_beta_[i_r] = density_profile_coefficients.beta(r);
        }

        #pragma omp parallel for
        for (int i_r = 0; i_r < grid.numberSmootherCircles(); i_r++) {
            const double r = grid.radius(i_r);
            double coeff_alpha = coeff_alpha_[i_r];;
            for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
                const double theta = grid.theta(i_theta);
                const double sin_theta = sin_theta_[i_theta];
                const double cos_theta = cos_theta_[i_theta];
                const double index = grid.index(i_r, i_theta);
                /* Calculate the elements of the Jacobian matrix for the transformation mapping */
                /* The Jacobian matrix is: */
                /* [Jrr, Jrt] */
                /* [Jtr, Jtt] */
                const double Jrr = domain_geometry.dFx_dr(r, theta, sin_theta, cos_theta);
                const double Jtr = domain_geometry.dFy_dr(r, theta, sin_theta, cos_theta);
                const double Jrt = domain_geometry.dFx_dt(r, theta, sin_theta, cos_theta);
                const double Jtt = domain_geometry.dFy_dt(r, theta, sin_theta, cos_theta);
                /* Compute the determinant of the Jacobian matrix */
                detDF_[index] = Jrr * Jtt - Jrt * Jtr;
                /* Compute the elements of the symmetric matrix: */
                /* 0.5 * alpha * DF^{-1} * DF^{-T} * |det(DF)| */
                /* which is represented by: */ \
                /* [arr, 0.5*art] */ \
                /* [0.5*atr, att] */ \
                arr_[index]  = 0.5 * (Jtt * Jtt + Jrt * Jrt) * coeff_alpha / fabs(detDF_[index] );
                att_[index]  = 0.5 * (Jtr * Jtr + Jrr * Jrr) * coeff_alpha / fabs(detDF_[index] );
                art_[index]  = (- Jtt * Jtr - Jrt * Jrr) * coeff_alpha / fabs(detDF_[index] );
                /* Note that the inverse Jacobian matrix DF^{-1} is: */
                /* 1.0 / det(DF) *   */
                /* [Jtt, -Jrt] */
                /* [-Jtr, Jrr] */
            }
        }

        #pragma omp parallel for
        for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
            const double theta = grid.theta(i_theta);
            const double sin_theta = sin_theta_[i_theta];
            const double cos_theta = cos_theta_[i_theta];
            for (int i_r = grid.numberSmootherCircles(); i_r < grid.nr(); i_r++) {
                const double r = grid.radius(i_r);
                const double index = grid.index(i_r, i_theta);
                double coeff_alpha = coeff_alpha_[i_r];
                /* Calculate the elements of the Jacobian matrix for the transformation mapping */
                /* The Jacobian matrix is: */
                /* [Jrr, Jrt] */
                /* [Jtr, Jtt] */
                const double Jrr = domain_geometry.dFx_dr(r, theta, sin_theta, cos_theta);
                const double Jtr = domain_geometry.dFy_dr(r, theta, sin_theta, cos_theta);
                const double Jrt = domain_geometry.dFx_dt(r, theta, sin_theta, cos_theta);
                const double Jtt = domain_geometry.dFy_dt(r, theta, sin_theta, cos_theta);
                /* Compute the determinant of the Jacobian matrix */
                detDF_[index] = Jrr * Jtt - Jrt * Jtr;
                /* Compute the elements of the symmetric matrix: */
                /* 0.5 * alpha * DF^{-1} * DF^{-T} * |det(DF)| */
                /* which is represented by: */ \
                /* [arr, 0.5*art] */ \
                /* [0.5*atr, att] */ \
                arr_[index]  = 0.5 * (Jtt * Jtt + Jrt * Jrt) * coeff_alpha / fabs(detDF_[index] );
                att_[index]  = 0.5 * (Jtr * Jtr + Jrr * Jrr) * coeff_alpha / fabs(detDF_[index] );
                art_[index]  = (- Jtt * Jtr - Jrt * Jrr) * coeff_alpha / fabs(detDF_[index] );
                /* Note that the inverse Jacobian matrix DF^{-1} is: */
                /* 1.0 / det(DF) *   */
                /* [Jtt, -Jrt] */
                /* [-Jtr, Jrr] */
            }
        }
    }
    else if(processing_type_ == ProcessingType::GPU){
        Vector<double> sin_theta_host(grid.ntheta());
        Vector<double> cos_theta_host(grid.ntheta());
        #pragma omp parallel for
        for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++){
            const double theta = grid.theta(i_theta);
            sin_theta_host[i_theta] = sin(theta);
            cos_theta_host[i_theta] = cos(theta);
        }
        copyHostToDevice(sin_theta_host, gpu_sin_theta_);
        copyHostToDevice(cos_theta_host, gpu_cos_theta_);

        Vector<double> coeff_alpha_host(grid.nr());
        Vector<double> coeff_beta_host(grid.nr());
        for (int i_r = 0; i_r < grid.nr(); i_r++){
            const double r = grid.radius(i_r);
            coeff_alpha_host[i_r] = density_profile_coefficients.alpha(r);
            coeff_beta_host[i_r] = density_profile_coefficients.beta(r);
        }
        copyHostToDevice(coeff_alpha_host, gpu_coeff_alpha_);
        copyHostToDevice(coeff_beta_host, gpu_coeff_beta_);
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

const Vector<double>& LevelCache::arr() const {
    return arr_;
}
const Vector<double>& LevelCache::att() const {
    return att_;
}
const Vector<double>& LevelCache::art() const {
    return art_;
}
const Vector<double>& LevelCache::detDF() const {
    return detDF_;
}

const GPU_Vector<double>& LevelCache::GPU_sin_theta() const {
    return gpu_sin_theta_;
}
const GPU_Vector<double>& LevelCache::GPU_cos_theta() const {
    return gpu_cos_theta_;
}

const GPU_Vector<double>& LevelCache::GPU_coeff_alpha() const {
    return gpu_coeff_alpha_;
}
const GPU_Vector<double>& LevelCache::GPU_coeff_beta() const {
    return gpu_coeff_beta_;
}