#pragma once

#include <cmath>

#include "../sourceTerm.h"

class Refined_ZoniShiftedGyro_CircularGeometry : public SourceTerm
{
public:
    explicit Refined_ZoniShiftedGyro_CircularGeometry(PolarGrid const& grid, double Rmax);
    virtual ~Refined_ZoniShiftedGyro_CircularGeometry() = default;

    double operator()(int i_r, int i_theta) const override;

private:
    PolarGrid const& grid_;
    const double Rmax = 1.3;
};
