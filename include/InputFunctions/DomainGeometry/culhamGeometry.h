#pragma once

#include <cmath>
#include <memory>
#include <array>
#include <cstdint>

#include "../domainGeometry.h"

class CulhamGeometry : public DomainGeometry
{
public:
    CulhamGeometry();
    explicit CulhamGeometry(const double& Rmax);

    virtual ~CulhamGeometry() = default;

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
    const double Rmax = 1.3;

    void initializeGeometry();

    double my_sum(std::array<double, 1001>& f, int64_t start_idx, int64_t end_idx) const;
    double q(double rr) const;
    double dq(double rr) const;
    double p(double rr) const;
    double dp(double rr) const;
    double dg(double rr, double g) const;
    double double_deriv(double rr, double c, double g, double dg, double val, double d_val) const;
    double g(double rr) const;
    double Delta(double rr) const;
    double Delta_prime(double rr) const;
    double E(double rr) const;
    double T(double rr) const;
    double E_prime(double rr) const;
    double T_prime(double rr) const;
    double P(double rr) const;
    double dP(double rr) const;

    double rr;
    double dr;
    double dr_h;
    std::array<double, 1001> r;
    int64_t i;
    double dg_1;
    double dE_1;
    double dT_1;
    double ddE_1;
    double ddT_1;
    double r2;
    double g_2;
    double dg_2;
    double E_2;
    double T_2;
    double dE_2;
    double dT_2;
    double ddE_2;
    double ddT_2;
    double g_3;
    double dg_3;
    double E_3;
    double T_3;
    double dE_3;
    double dT_3;
    double ddE_3;
    double ddT_3;
    double g_4;
    double dg_4;
    double E_4;
    double T_4;
    double dE_4;
    double dT_4;
    double ddE_4;
    double ddT_4;
    double current_Ea;
    double current_Ta;
    std::array<double, 1001> f;
    std::array<double, 1001> integ_contents;
    double integral;
    double current_Delta_a;
    size_t i_0001;
    std::array<double, 1001> g_array;
    std::array<double, 1001> Delta_array;
    std::array<double, 1001> Delta_prime_array;
    std::array<double, 1001> E_array;
    std::array<double, 1001> T_array;
    std::array<double, 1001> E_prime_array;
    std::array<double, 1001> T_prime_array;
};

#include "culhamGeometry.inl"