#pragma once

#include <stdlib.h>
#include <array>
#include <stdint.h>
#include <vector>
#include "exact_funcs.h"

class PolarR6GyroZoniShiftedCulham: public ExactFuncs
{
public:
    PolarR6GyroZoniShiftedCulham();
    double x(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const override;
    void x(std::vector<double> const& r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol) const override;
    void x(double r, std::vector<double> const& theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const override;
    double y(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const override;
    void y(std::vector<double> const& r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol) const override;
    void y(double r, std::vector<double> const& theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const override;
    double J_rr(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const override;
    void J_rr(std::vector<double> const& r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol) const override;
    void J_rr(double r, std::vector<double> const& theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const override;
    double J_rt(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const override;
    void J_rt(std::vector<double> const& r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol) const override;
    void J_rt(double r, std::vector<double> const& theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const override;
    double J_tr(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const override;
    void J_tr(std::vector<double> const& r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol) const override;
    void J_tr(double r, std::vector<double> const& theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const override;
    double J_tt(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const override;
    void J_tt(std::vector<double> const& r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol) const override;
    void J_tt(double r, std::vector<double> const& theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const override;
    double coeffs1(double r, double Rmax) const override;
    void coeffs1(std::vector<double> const& r, double Rmax, std::vector<double>& sol) const override;
    double coeffs2(double r, double Rmax) const override;
    void coeffs2(std::vector<double> const& r, double Rmax, std::vector<double>& sol) const override;
    double rho_glob(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const override;
    void rho_glob(std::vector<double> const& r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol) const override;
    void rho_glob(double r, std::vector<double> const& theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const override;
    double rho_pole(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const override;
    void rho_pole(std::vector<double> const& r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol) const override;
    void rho_pole(double r, std::vector<double> const& theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const override;
    double phi_exact(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const override;
    void phi_exact(std::vector<double> const& r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol) const override;
    void phi_exact(double r, std::vector<double> const& theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const override;
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
    void J_xr(std::vector<double> const& r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol) const;
    void J_xr(double r, std::vector<double> const& theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const;
    double J_xq(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const;
    void J_xq(std::vector<double> const& r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol) const;
    void J_xq(double r, std::vector<double> const& theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const;
    double J_yr(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const;
    void J_yr(std::vector<double> const& r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol) const;
    void J_yr(double r, std::vector<double> const& theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const;
    double J_yq(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const;
    void J_yq(std::vector<double> const& r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol) const;
    void J_yq(double r, std::vector<double> const& theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const;
    std::array<double, 1001> g_array;
    std::array<double, 1001> Delta_array;
    std::array<double, 1001> Delta_prime_array;
    std::array<double, 1001> E_array;
    std::array<double, 1001> T_array;
    std::array<double, 1001> E_prime_array;
    std::array<double, 1001> T_prime_array;
};


