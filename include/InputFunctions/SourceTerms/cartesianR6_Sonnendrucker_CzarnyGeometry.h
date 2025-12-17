#pragma once

#include <cmath>

#include "../sourceTerm.h"

#include "../../PolarGrid/polargrid.h"

class CartesianR6_Sonnendrucker_CzarnyGeometry : public SourceTerm
{
public:
    explicit CartesianR6_Sonnendrucker_CzarnyGeometry(PolarGrid const& grid, double Rmax,
                                                      double inverse_aspect_ratio_epsilon, double ellipticity_e);
    virtual ~CartesianR6_Sonnendrucker_CzarnyGeometry() = default;

    double operator()(std::size_t i_r, std::size_t i_theta) const override;

private:
    PolarGrid const& grid_;
    const double Rmax                         = 1.3;
    const double inverse_aspect_ratio_epsilon = 0.3;
    const double ellipticity_e                = 1.4;

    void initializeGeometry();
    double factor_xi;
};
