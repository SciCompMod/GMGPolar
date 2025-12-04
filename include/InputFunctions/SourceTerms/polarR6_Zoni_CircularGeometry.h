#pragma once

#include <cmath>

#include "../sourceTerm.h"

class PolarR6_Zoni_CircularGeometry : public SourceTerm
{
public:
    PolarR6_Zoni_CircularGeometry() = default;
    explicit PolarR6_Zoni_CircularGeometry(double Rmax);
    virtual ~PolarR6_Zoni_CircularGeometry() = default;

    double rhs_f(double r, double theta) const override;

private:
    const double Rmax = 1.3;
};
