#pragma once

#include <cmath>

#include "../sourceTerm.h"

class PolarR6_ZoniShifted_CircularGeometry : public SourceTerm
{
public:
    PolarR6_ZoniShifted_CircularGeometry() = default;
    explicit PolarR6_ZoniShifted_CircularGeometry(const double& Rmax);
    virtual ~PolarR6_ZoniShifted_CircularGeometry() = default;

    double rhs_f(const double& r, const double& theta) const override;

private:
    const double Rmax = 1.3;
};
