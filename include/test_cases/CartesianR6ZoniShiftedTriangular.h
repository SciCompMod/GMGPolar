#pragma once

#include <stdlib.h>
#include <vector>
#include "exact_funcs.h"


class CartesianR6ZoniShiftedTriangular: public ExactFuncs
{
public:
    double x(double r, double theta, double map2_epsilon, double map2_e, double Rmax) const override;
    void x(std::vector<double> const& r, double theta, double map2_epsilon, double map2_e, double Rmax, std::vector<double>& sol) const override;
    void x(double r, std::vector<double> const& theta, double map2_epsilon, double map2_e, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const override;
    double y(double r, double theta, double map2_epsilon, double map2_e, double Rmax) const override;
    void y(std::vector<double> const& r, double theta, double map2_epsilon, double map2_e, double Rmax, std::vector<double>& sol) const override;
    void y(double r, std::vector<double> const& theta, double map2_epsilon, double map2_e, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const override;
    double J_rr(double r, double theta, double map2_epsilon, double map2_e, double Rmax) const override;
    void J_rr(std::vector<double> const& r, double theta, double map2_epsilon, double map2_e, double Rmax, std::vector<double>& sol) const override;
    void J_rr(double r, std::vector<double> const& theta, double map2_epsilon, double map2_e, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const override;
    double J_rt(double r, double theta, double map2_epsilon, double map2_e, double Rmax) const override;
    void J_rt(std::vector<double> const& r, double theta, double map2_epsilon, double map2_e, double Rmax, std::vector<double>& sol) const override;
    void J_rt(double r, std::vector<double> const& theta, double map2_epsilon, double map2_e, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const override;
    double J_tr(double r, double theta, double map2_epsilon, double map2_e, double Rmax) const override;
    void J_tr(std::vector<double> const& r, double theta, double map2_epsilon, double map2_e, double Rmax, std::vector<double>& sol) const override;
    void J_tr(double r, std::vector<double> const& theta, double map2_epsilon, double map2_e, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const override;
    double J_tt(double r, double theta, double map2_epsilon, double map2_e, double Rmax) const override;
    void J_tt(std::vector<double> const& r, double theta, double map2_epsilon, double map2_e, double Rmax, std::vector<double>& sol) const override;
    void J_tt(double r, std::vector<double> const& theta, double map2_epsilon, double map2_e, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const override;
    double rho_glob(double r, double theta, double map2_epsilon, double map2_e, double Rmax) const override;
    void rho_glob(std::vector<double> const& r, double theta, double map2_epsilon, double map2_e, double Rmax, std::vector<double>& sol) const override;
    void rho_glob(double r, std::vector<double> const& theta, double map2_epsilon, double map2_e, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const override;
    double rho_pole(double r, double theta, double map2_epsilon, double map2_e, double Rmax) const override;
    void rho_pole(std::vector<double> const& r, double theta, double map2_epsilon, double map2_e, double Rmax, std::vector<double>& sol) const override;
    void rho_pole(double r, std::vector<double> const& theta, double map2_epsilon, double map2_e, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const override;
    double coeffs1(double r, double Rmax) const override;
    void coeffs1(std::vector<double> const& r, double Rmax, std::vector<double>& sol) const override;
    double coeffs2(double r, double Rmax) const override;
    void coeffs2(std::vector<double> const& r, double Rmax, std::vector<double>& sol) const override;
    double phi_exact(double r, double theta, double map2_epsilon, double map2_e, double Rmax) const override;
    void phi_exact(std::vector<double> const& r, double theta, double map2_epsilon, double map2_e, double Rmax, std::vector<double>& sol) const override;
    void phi_exact(double r, std::vector<double> const& theta, double map2_epsilon, double map2_e, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const override;
private:
    double J_xs(double r, double theta, double map2_epsilon, double map2_e, double Rmax) const;
    void J_xs(std::vector<double> const& r, double theta, double map2_epsilon, double map2_e, double Rmax, std::vector<double>& sol) const;
    void J_xs(double r, std::vector<double> const& theta, double map2_epsilon, double map2_e, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const;
    double J_xt(double r, double theta, double map2_epsilon, double map2_e, double Rmax) const;
    void J_xt(std::vector<double> const& r, double theta, double map2_epsilon, double map2_e, double Rmax, std::vector<double>& sol) const;
    void J_xt(double r, std::vector<double> const& theta, double map2_epsilon, double map2_e, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const;
    double J_ys(double r, double theta, double map2_epsilon, double map2_e, double Rmax) const;
    void J_ys(std::vector<double> const& r, double theta, double map2_epsilon, double map2_e, double Rmax, std::vector<double>& sol) const;
    void J_ys(double r, std::vector<double> const& theta, double map2_epsilon, double map2_e, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const;
    double J_yt(double r, double theta, double map2_epsilon, double map2_e, double Rmax) const;
    void J_yt(std::vector<double> const& r, double theta, double map2_epsilon, double map2_e, double Rmax, std::vector<double>& sol) const;
    void J_yt(double r, std::vector<double> const& theta, double map2_epsilon, double map2_e, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const;
};


