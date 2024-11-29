#pragma once

#include <cmath>

#include "../sourceTerm.h"

class CartesianR2_ZoniShiftedGyro_ShafranovGeometry : public SourceTerm { 
public:
    CartesianR2_ZoniShiftedGyro_ShafranovGeometry() = default;
    explicit CartesianR2_ZoniShiftedGyro_ShafranovGeometry(const double& Rmax, const double& elongation_kappa, const double& shift_delta);
    virtual ~CartesianR2_ZoniShiftedGyro_ShafranovGeometry() = default;

    double rhs_f(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const override;

private:
    const double Rmax = 1.3;
    const double elongation_kappa = 0.3;
    const double shift_delta = 0.2;
};
