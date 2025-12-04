#pragma once

#include <cmath>

#include "../sourceTerm.h"

class Refined_ZoniShiftedGyro_CircularGeometry : public SourceTerm
{
public:
    Refined_ZoniShiftedGyro_CircularGeometry() = default;
    explicit Refined_ZoniShiftedGyro_CircularGeometry(const double& Rmax);
    virtual ~Refined_ZoniShiftedGyro_CircularGeometry() = default;

    double rhs_f(const double& r, const double& theta) const override;

private:
    const double Rmax = 1.3;
};
