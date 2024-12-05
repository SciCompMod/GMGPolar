#pragma once

#include <cmath>

#include "../domainGeometry.h"

/* Triangular Shape and Ellipticity */

class CzarnyGeometry : public DomainGeometry
{
public:
    explicit CzarnyGeometry();
    explicit CzarnyGeometry(const double& Rmax, const double& inverse_aspect_ratio_epsilon,
                            const double& ellipticity_e);

    virtual ~CzarnyGeometry() = default;

    double Fx(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const override;
    double Fy(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const override;
    double dFx_dr(const double& r, const double& theta, const double& sin_theta,
                  const double& cos_theta) const override;
    double dFy_dr(const double& r, const double& theta, const double& sin_theta,
                  const double& cos_theta) const override;
    double dFx_dt(const double& r, const double& theta, const double& sin_theta,
                  const double& cos_theta) const override;
    double dFy_dt(const double& r, const double& theta, const double& sin_theta,
                  const double& cos_theta) const override;

private:
    const double Rmax                         = 1.3;
    const double inverse_aspect_ratio_epsilon = 0.3;
    const double ellipticity_e                = 1.4;

    void initializeGeometry();
    double factor_xi;
};

#include "czarnyGeometry.inl"
