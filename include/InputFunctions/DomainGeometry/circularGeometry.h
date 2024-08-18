#pragma once

#include <cmath>

#include "../domainGeometry.h"

class CircularGeometry : public DomainGeometry {
public:
    CircularGeometry() = default;
    explicit CircularGeometry(const double& Rmax);

    virtual ~CircularGeometry() = default;

    double Fx(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const override;
    double Fy(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const override;
    double dFx_dr(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const override;
    double dFy_dr(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const override;
    double dFx_dt(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const override;
    double dFy_dt(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const override;

private:
    const double Rmax = 1.3;
};

#include "circularGeometry.inl"