#pragma once

#include <cmath>

#include "../sourceTerm.h"

class Refined_ZoniShiftedGyro_CulhamGeometry : public SourceTerm
{
public:

    explicit Refined_ZoniShiftedGyro_CulhamGeometry(PolarGrid const& grid, double Rmax);
    virtual ~Refined_ZoniShiftedGyro_CulhamGeometry() = default;

    double rhs_f(double r, double theta) const override;

private:
    PolarGrid const& grid_;
    const double Rmax = 1.3;
};
