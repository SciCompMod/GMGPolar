#pragma once

#include <cuda_runtime.h>
#include <cassert>
#include <cmath>

class BoundaryConditions
{
public:
    BoundaryConditions();
    explicit BoundaryConditions(
        const double& Rmax, const double& inverse_aspect_ratio_epsilon, const double& ellipticity_e);

    ~BoundaryConditions() = default;

    __host__ __device__ __forceinline__ 
    double u_D(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
        return 0.4096 * pow((r/Rmax), 6.0) * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta);
    }

    /* Only used if DirBC_Interior = true */
    __host__ __device__ __forceinline__ 
    double u_D_Interior(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
        return 0.4096 * pow((r/Rmax), 6.0) * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta);
    }

private:
    const double Rmax = 1.3;
    const double inverse_aspect_ratio_epsilon = 0.3;
    const double ellipticity_e = 1.4;

    void initializeGeometry();
    double factor_xi;
};
