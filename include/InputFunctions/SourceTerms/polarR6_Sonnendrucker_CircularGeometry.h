#pragma once

#include <cmath>

#include "../sourceTerm.h"

class PolarR6_Sonnendrucker_CircularGeometry : public SourceTerm { 
public:
    PolarR6_Sonnendrucker_CircularGeometry() = default;
    explicit PolarR6_Sonnendrucker_CircularGeometry(const double& Rmax);
    virtual ~PolarR6_Sonnendrucker_CircularGeometry() = default;

    double rhs_f(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const override;

private:
    const double Rmax = 1.3;
};
