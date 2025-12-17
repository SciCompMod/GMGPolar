#pragma once

#include <cmath>

#include "../sourceTerm.h"

class CartesianR2_ZoniShiftedGyro_ShafranovGeometry : public SourceTerm
{
public:

    explicit CartesianR2_ZoniShiftedGyro_ShafranovGeometry(PolarGrid const& grid, double Rmax, double elongation_kappa, double shift_delta);
    virtual ~CartesianR2_ZoniShiftedGyro_ShafranovGeometry() = default;

    double operator()(double r, double theta) const override;

private:
    PolarGrid const& grid_;
    const double Rmax             = 1.3;
    const double elongation_kappa = 0.3;
    const double shift_delta      = 0.2;
};
