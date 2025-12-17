#pragma once

#include <cmath>

#include "../sourceTerm.h"

class CartesianR6_ZoniShiftedGyro_ShafranovGeometry : public SourceTerm
{
public:
    explicit CartesianR6_ZoniShiftedGyro_ShafranovGeometry(PolarGrid const& grid, double Rmax, double elongation_kappa,
                                                           double shift_delta);
    virtual ~CartesianR6_ZoniShiftedGyro_ShafranovGeometry() = default;

    double operator()(int i_r, int i_theta) const override;

private:
    PolarGrid const& grid_;
    const double Rmax             = 1.3;
    const double elongation_kappa = 0.3;
    const double shift_delta      = 0.2;
};
