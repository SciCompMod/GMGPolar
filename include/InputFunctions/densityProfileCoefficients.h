#pragma once

#include <cuda_runtime.h>
#include <cmath>

class DensityProfileCoefficients
{
public:
    DensityProfileCoefficients() = default;
    explicit DensityProfileCoefficients(const double& Rmax, const double& alpha);

    ~DensityProfileCoefficients() = default;

    __host__ __device__ __forceinline__ 
    double alpha(const double& r) const {
        return 0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r/Rmax) - 11.1111111111111);
    }

    __host__ __device__ __forceinline__ 
    double beta(const double& r) const {
        return pow((0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r/Rmax) - 11.1111111111111)), (double)((-1)));
    }

    // Only used in custom mesh generation -> refinement_radius
    double getAlphaJump() const;

private:
    const double Rmax = 1.3;
    const double alpha_jump = 0.66 * 1.3;
};
