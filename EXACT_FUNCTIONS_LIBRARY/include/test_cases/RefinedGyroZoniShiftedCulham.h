#pragma once

#include <stdlib.h>
#include <array>
#include <stdint.h>
#include <vector>
#include <math.h>
#include <stdbool.h>
#include "../ExactFunctions/exactfunctions.h"

class RefinedGyroZoniShiftedCulham: public ExactFunctions
{
public:
    RefinedGyroZoniShiftedCulham();

    double x(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const override;
    double y(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const override;
    double x(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, double sin_theta, double cos_theta) const override;
    double y(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, double sin_theta, double cos_theta) const override;

    double J_rr(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const override;
    double J_rt(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const override;
    double J_tr(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const override;
    double J_tt(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const override;
    double J_rr(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, double sin_theta, double cos_theta) const override;
    double J_rt(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, double sin_theta, double cos_theta) const override;
    double J_tr(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, double sin_theta, double cos_theta) const override;
    double J_tt(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, double sin_theta, double cos_theta) const override;

    double coeffs1(double r, double Rmax) const override;
    double coeffs2(double r, double Rmax) const override;

    double rho_glob(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const override;
    double rho_pole(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const override;
    double rho_glob(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, double sin_theta, double cos_theta) const override;
    double rho_pole(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, double sin_theta, double cos_theta) const override;

    double phi_exact(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const override;
    double phi_exact(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, double sin_theta, double cos_theta) const override;

private:
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

    double J_xr(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const;
    double J_xq(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const;
    double J_yr(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const;
    double J_yq(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const;
    double J_xr(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, double sin_theta, double cos_theta) const;
    double J_xq(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, double sin_theta, double cos_theta) const;
    double J_yr(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, double sin_theta, double cos_theta) const;
    double J_yq(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, double sin_theta, double cos_theta) const;

    std::array<double, 1001> g_array;
    std::array<double, 1001> Delta_array;
    std::array<double, 1001> Delta_prime_array;
    std::array<double, 1001> E_array;
    std::array<double, 1001> T_array;
    std::array<double, 1001> E_prime_array;
    std::array<double, 1001> T_prime_array;
};
