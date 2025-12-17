#pragma once

#include <cmath>

#include "../sourceTerm.h"

#include "../../PolarGrid/polargrid.h"

class PolarR6_ZoniShifted_CircularGeometry : public SourceTerm
{
public:
    explicit PolarR6_ZoniShifted_CircularGeometry(PolarGrid const& grid, double Rmax);
    virtual ~PolarR6_ZoniShifted_CircularGeometry() = default;

    double operator()(int i_r, int i_theta) const override;

private:
    PolarGrid const& grid_;
    const double Rmax = 1.3;
};
