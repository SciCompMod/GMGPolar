#pragma once

#include <cmath>

#include "../sourceTerm.h"

class PolarR6_ZoniShifted_ShafranovGeometry : public SourceTerm
{
public:
    explicit PolarR6_ZoniShifted_ShafranovGeometry(PolarGrid const& grid, double Rmax, double elongation_kappa,
                                                   double shift_delta);
    virtual ~PolarR6_ZoniShifted_ShafranovGeometry() = default;

    double operator()(int i_r, int i_theta) const override;

private:
    PolarGrid const& grid_;
    const double Rmax             = 1.3;
    const double elongation_kappa = 0.3;
    const double shift_delta      = 0.2;
};
