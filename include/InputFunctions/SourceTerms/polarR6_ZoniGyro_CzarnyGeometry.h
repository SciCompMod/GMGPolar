#pragma once

#include <cmath>

#include "../sourceTerm.h"

#include "../../PolarGrid/polargrid.h"

class PolarR6_ZoniGyro_CzarnyGeometry : public SourceTerm
{
public:
    explicit PolarR6_ZoniGyro_CzarnyGeometry(PolarGrid const& grid, double Rmax, double inverse_aspect_ratio_epsilon,
                                             double ellipticity_e);
    virtual ~PolarR6_ZoniGyro_CzarnyGeometry() = default;

    double operator()(std::size_t i_r, std::size_t i_theta) const override;

private:
    PolarGrid const& grid_;
    const double Rmax                         = 1.3;
    const double inverse_aspect_ratio_epsilon = 0.3;
    const double ellipticity_e                = 1.4;

    void initializeGeometry();
    double factor_xi;
};
