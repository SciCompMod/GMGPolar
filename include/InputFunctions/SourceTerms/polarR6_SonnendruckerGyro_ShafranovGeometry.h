#pragma once

#include <cmath>

#include "../sourceTerm.h"

class PolarR6_SonnendruckerGyro_ShafranovGeometry : public SourceTerm
{
public:

    explicit PolarR6_SonnendruckerGyro_ShafranovGeometry(PolarGrid const& grid, double Rmax, double elongation_kappa, double shift_delta);
    virtual ~PolarR6_SonnendruckerGyro_ShafranovGeometry() = default;

    double operator()(double r, double theta) const override;

private:
    PolarGrid const& grid_;
    const double Rmax             = 1.3;
    const double elongation_kappa = 0.3;
    const double shift_delta      = 0.2;
};
