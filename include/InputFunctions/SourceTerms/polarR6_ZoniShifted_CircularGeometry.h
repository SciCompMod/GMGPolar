#pragma once

#include <cmath>

#include "../sourceTerm.h"

class PolarR6_ZoniShifted_CircularGeometry : public SourceTerm
{
public:

    explicit PolarR6_ZoniShifted_CircularGeometry(PolarGrid const& grid, double Rmax);
    virtual ~PolarR6_ZoniShifted_CircularGeometry() = default;

    double operator()(double r, double theta) const override;

private:
    PolarGrid const& grid_;
    const double Rmax = 1.3;
};
