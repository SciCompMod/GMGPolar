#pragma once

#include <cmath>

#include "../sourceTerm.h"

class CartesianR6_Sonnendrucker_CircularGeometry : public SourceTerm
{
public:

    explicit CartesianR6_Sonnendrucker_CircularGeometry(PolarGrid const& grid, double Rmax);
    virtual ~CartesianR6_Sonnendrucker_CircularGeometry() = default;

    double rhs_f(double r, double theta) const override;

private:
    PolarGrid const& grid_;
    const double Rmax = 1.3;
};
