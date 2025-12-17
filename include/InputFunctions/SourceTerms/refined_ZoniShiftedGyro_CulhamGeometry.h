#pragma once

#include <cmath>

#include "../sourceTerm.h"

#include "../../PolarGrid/polargrid.h"

class Refined_ZoniShiftedGyro_CulhamGeometry : public SourceTerm
{
public:
    explicit Refined_ZoniShiftedGyro_CulhamGeometry(PolarGrid const& grid, double Rmax);
    virtual ~Refined_ZoniShiftedGyro_CulhamGeometry() = default;

    double operator()(int i_r, int i_theta) const override;

private:
    PolarGrid const& grid_;
    const double Rmax = 1.3;
};
