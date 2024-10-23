#pragma once

#include <cmath>
#include <memory>
#include <array>
#include <cstdint>

#include "../domainGeometry.h"

class CulhamGeometry : public DomainGeometry {
public:
    CulhamGeometry();
    explicit CulhamGeometry(const double& Rmax);

    virtual ~CulhamGeometry() = default;

    double Fx(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const override;
    double Fy(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const override;
    double dFx_dr(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const override;
    double dFy_dr(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const override;
    double dFx_dt(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const override;
    double dFy_dt(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const override;

private:
    const double Rmax = 1.3;

    void initializeGeometry();


};

// #include "culhamGeometry2.inl"