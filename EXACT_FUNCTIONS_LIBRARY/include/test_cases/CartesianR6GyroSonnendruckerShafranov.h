#pragma once

#include <stdlib.h>
#include <vector>
#include <cmath>
#include "../ExactFunctions/exactfunctions.h"

class CartesianR6GyroSonnendruckerShafranov: public ExactFunctions
{
public:
    double x(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const override;
    double y(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const override;
    double x(double r, double theta, double map1_kappa, double map1_delta, double Rmax, double sin_theta, double cos_theta) const override;
    double y(double r, double theta, double map1_kappa, double map1_delta, double Rmax, double sin_theta, double cos_theta) const override;
 
    double J_rr(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const override;
    double J_rt(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const override;
    double J_tr(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const override;
    double J_tt(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const override;
    double J_rr(double r, double theta, double map1_kappa, double map1_delta, double Rmax, double sin_theta, double cos_theta) const override;
    double J_rt(double r, double theta, double map1_kappa, double map1_delta, double Rmax, double sin_theta, double cos_theta) const override;
    double J_tr(double r, double theta, double map1_kappa, double map1_delta, double Rmax, double sin_theta, double cos_theta) const override;
    double J_tt(double r, double theta, double map1_kappa, double map1_delta, double Rmax, double sin_theta, double cos_theta) const override;

    double coeffs1(double r, double Rmax) const override;
    double coeffs2(double r, double Rmax) const override;

    double rho_glob(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const override;
    double rho_pole(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const override;
    double rho_glob(double r, double theta, double map1_kappa, double map1_delta, double Rmax, double sin_theta, double cos_theta) const override;
    double rho_pole(double r, double theta, double map1_kappa, double map1_delta, double Rmax, double sin_theta, double cos_theta) const override;

    double phi_exact(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const override;
    double phi_exact(double r, double theta, double map1_kappa, double map1_delta, double Rmax, double sin_theta, double cos_theta) const override;

private:
    double J_xs(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const;
    double J_xt(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const;
    double J_ys(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const;
    double J_yt(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const;
    double J_xs(double r, double theta, double map1_kappa, double map1_delta, double Rmax, double sin_theta, double cos_theta) const;
    double J_xt(double r, double theta, double map1_kappa, double map1_delta, double Rmax, double sin_theta, double cos_theta) const;
    double J_ys(double r, double theta, double map1_kappa, double map1_delta, double Rmax, double sin_theta, double cos_theta) const;
    double J_yt(double r, double theta, double map1_kappa, double map1_delta, double Rmax, double sin_theta, double cos_theta) const;
};

