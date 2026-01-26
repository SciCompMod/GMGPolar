#pragma once

#include <cmath>

#include "../sourceTerm.h"

#include "../../PolarGrid/polargrid.h"

class CartesianR2_SonnendruckerGyro_CircularGeometry : public SourceTerm
{
public:
    explicit CartesianR2_SonnendruckerGyro_CircularGeometry(PolarGrid const& grid, double Rmax);
    virtual ~CartesianR2_SonnendruckerGyro_CircularGeometry() = default;

    double operator()(std::size_t i_r, std::size_t i_theta) const override;

private:
    PolarGrid const& grid_;
    const double Rmax = 1.3;
};
