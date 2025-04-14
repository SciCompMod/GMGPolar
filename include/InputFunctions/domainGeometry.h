#pragma once

#include <cuda_runtime.h>
#include <cmath>
#include <cassert>

class DomainGeometry {
public:
    DomainGeometry();
#ifdef GEOM_SHAFRANOV     
    explicit DomainGeometry(const double& Rmax, const double& elongation_kappa, const double& shift_delta);
#else
    explicit DomainGeometry(const double& Rmax, const double& inverse_aspect_ratio_epsilon, const double& ellipticity_e);
#endif

    ~DomainGeometry() = default;

#ifdef GEOM_SHAFRANOV    
    /* Invertible mapping F: \Omega_{Ref} -> \Omega. */
    __host__ __device__ __forceinline__ 
    double Fx(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
        return (1.0 - elongation_kappa) * (r/Rmax) * cos_theta - shift_delta * (r/Rmax) * (r/Rmax);
    }
    __host__ __device__ __forceinline__ 
    double Fy(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
        return (1.0 + elongation_kappa) * (r/Rmax) * sin_theta;
    }

    /* Jacobian matrix of the mapping F. */
    __host__ __device__ __forceinline__ 
    double dFx_dr(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
        return ( (Rmax - elongation_kappa * Rmax) * cos_theta - 2.0 * shift_delta * r ) / ( Rmax * Rmax );
    }
    __host__ __device__ __forceinline__ 
    double dFy_dr(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
        return (elongation_kappa + 1.0) * sin_theta / Rmax;
    }
    __host__ __device__ __forceinline__ 
    double dFx_dt(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
        return ((elongation_kappa - 1.0) * r * sin_theta) / Rmax;
    }
    __host__ __device__ __forceinline__ 
    double dFy_dt(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
        return ((elongation_kappa + 1.0) * r * cos_theta) / Rmax;
    }
#else
    /* Invertible mapping F: \Omega_{Ref} -> \Omega. */
    __host__ __device__ __forceinline__ 
    double Fx(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
        double temp = sqrt(1.0 + inverse_aspect_ratio_epsilon * (inverse_aspect_ratio_epsilon + 2.0 * (r / Rmax) * cos_theta));
        return (1.0 - temp) / inverse_aspect_ratio_epsilon;
    }
    __host__ __device__ __forceinline__ 
    double Fy(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
        double temp = sqrt(1.0 + inverse_aspect_ratio_epsilon * (inverse_aspect_ratio_epsilon + 2.0 * (r / Rmax) * cos_theta));
        return ellipticity_e * factor_xi * (r / Rmax) * sin_theta / (1.0 + (1.0 - temp));
    }

    /* Jacobian matrix of the mapping F. */
    __host__ __device__ __forceinline__ 
    double dFx_dr(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
        double temp = sqrt(1.0 + inverse_aspect_ratio_epsilon * (2.0 * (r / Rmax) * cos_theta + inverse_aspect_ratio_epsilon));
        return -cos_theta / (Rmax * temp);
    }
    __host__ __device__ __forceinline__ 
    double dFy_dr(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
        double temp = sqrt(1.0 + inverse_aspect_ratio_epsilon * (2.0 * (r / Rmax) * cos_theta + inverse_aspect_ratio_epsilon));
        return (ellipticity_e * factor_xi * sin_theta) / (Rmax * (2.0 - temp)) +
            (ellipticity_e * factor_xi * inverse_aspect_ratio_epsilon * r * sin_theta * cos_theta) /
            (Rmax * Rmax * temp * (2.0 - temp) * (2.0 - temp));
    }
    __host__ __device__ __forceinline__ 
    double dFx_dt(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
        double temp = sqrt(1.0 + inverse_aspect_ratio_epsilon * (2.0 * (r / Rmax) * cos_theta + inverse_aspect_ratio_epsilon));
        return (r / Rmax) * sin_theta / temp;
    }
    __host__ __device__ __forceinline__ 
    double dFy_dt(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
        double temp = sqrt(1.0 + inverse_aspect_ratio_epsilon * (2.0 * (r / Rmax) * cos_theta + inverse_aspect_ratio_epsilon));
        return (ellipticity_e * factor_xi * (r / Rmax) * cos_theta) / (2.0 - temp) -
            (ellipticity_e * factor_xi * inverse_aspect_ratio_epsilon * (r / Rmax) * (r / Rmax) * sin_theta * sin_theta) /
            (temp * (2.0 - temp) * (2.0 - temp));
    }
#endif

private:
    const double Rmax = 1.3;
    
#ifdef GEOM_SHAFRANOV
    // Shafranov
    const double elongation_kappa = 0.3;
    const double shift_delta      = 0.2; 
#else
    // Czarny
    const double inverse_aspect_ratio_epsilon = 0.3;
    const double ellipticity_e = 1.4; 
    double factor_xi;
#endif


    void initializeGeometry();
};
