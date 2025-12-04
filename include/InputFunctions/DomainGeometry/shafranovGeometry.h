#pragma once

#include <cmath>

#include "../domainGeometry.h"

/* Strechted Ellipse with a Shafranov Shift */

class ShafranovGeometry : public DomainGeometry
{
public:
    ShafranovGeometry() = default;
    explicit ShafranovGeometry(const double& Rmax, const double& elongation_kappa, const double& shift_delta);

    virtual ~ShafranovGeometry() = default;

    double Fx(const double& r, const double& theta) const override;
    double Fy(const double& r, const double& theta) const override;
    double dFx_dr(const double& r, const double& theta) const override;
    double dFy_dr(const double& r, const double& theta) const override;
    double dFx_dt(const double& r, const double& theta) const override;
    double dFy_dt(const double& r, const double& theta) const override;

private:
    const double Rmax             = 1.3;
    const double elongation_kappa = 0.3;
    const double shift_delta      = 0.2;
};

#include "shafranovGeometry.inl"