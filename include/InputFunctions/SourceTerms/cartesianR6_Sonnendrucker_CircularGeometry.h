#pragma once

#include <cmath>

#include "../sourceTerm.h"

#include "../../PolarGrid/polargrid.h"

class CartesianR6_Sonnendrucker_CircularGeometry : public SourceTerm
{
public:
    explicit CartesianR6_Sonnendrucker_CircularGeometry(PolarGrid const& grid, double Rmax);
    virtual ~CartesianR6_Sonnendrucker_CircularGeometry() = default;

    double operator()(std::size_t i_r, std::size_t i_theta) const override;

private:
    PolarGrid const& grid_;
    const double Rmax = 1.3;
};
