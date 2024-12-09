#pragma once

#include <cuda_runtime.h>
#include <cmath>
#include <cassert>

class ExactSolution {
public:
    ExactSolution();
    ~ExactSolution() = default;
    explicit ExactSolution(const double& Rmax, const double& inverse_aspect_ratio_epsilon, const double& ellipticity_e);

    __host__ __device__ 
    double exact_solution(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
        return 0.4096 * pow((r/Rmax), 6.0) * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta);
    }
    
private:
    const double Rmax = 1.3;
    const double inverse_aspect_ratio_epsilon = 0.3;
    const double ellipticity_e = 1.4;

    double factor_xi;
    void initializeGeometry();
};


