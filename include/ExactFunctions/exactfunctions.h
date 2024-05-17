#pragma once

class ExactFunctions
{
public:
    virtual ~ExactFunctions()
    {
    }
    virtual double x(double r, double theta, double kappa_eps, double delta_e, double Rmax) const = 0;
    virtual double y(double r, double theta, double kappa_eps, double delta_e, double Rmax) const = 0;
    virtual double x(double r, double theta, double kappa_eps, double delta_e, double Rmax, double sin_theta, double cos_theta) const = 0;
    virtual double y(double r, double theta, double kappa_eps, double delta_e, double Rmax, double sin_theta, double cos_theta) const = 0;

    virtual double J_rr(double r, double theta, double kappa_eps, double delta_e, double Rmax) const = 0;
    virtual double J_rt(double r, double theta, double kappa_eps, double delta_e, double Rmax) const = 0;
    virtual double J_tr(double r, double theta, double kappa_eps, double delta_e, double Rmax) const = 0;
    virtual double J_tt(double r, double theta, double kappa_eps, double delta_e, double Rmax) const = 0;
    virtual double J_rr(double r, double theta, double kappa_eps, double delta_e, double Rmax, double sin_theta, double cos_theta) const = 0;
    virtual double J_rt(double r, double theta, double kappa_eps, double delta_e, double Rmax, double sin_theta, double cos_theta) const = 0;
    virtual double J_tr(double r, double theta, double kappa_eps, double delta_e, double Rmax, double sin_theta, double cos_theta) const = 0;
    virtual double J_tt(double r, double theta, double kappa_eps, double delta_e, double Rmax, double sin_theta, double cos_theta) const = 0;
    
    virtual double coeffs1(double r, double Rmax) const = 0;
    virtual double coeffs2(double r, double Rmax) const = 0;

    virtual double rho_glob(double r, double theta, double kappa_eps, double delta_e, double Rmax) const = 0;        
    virtual double rho_pole(double r, double theta, double kappa_eps, double delta_e, double Rmax) const = 0;
    virtual double rho_glob(double r, double theta, double kappa_eps, double delta_e, double Rmax, double sin_theta, double cos_theta) const = 0;        
    virtual double rho_pole(double r, double theta, double kappa_eps, double delta_e, double Rmax, double sin_theta, double cos_theta) const = 0;

    virtual double phi_exact(double r, double theta, double kappa_eps, double delta_e, double Rmax) const = 0;
    virtual double phi_exact(double r, double theta, double kappa_eps, double delta_e, double Rmax, double sin_theta, double cos_theta) const = 0;
};
