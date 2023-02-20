#include "CartesianR2PoissonShafranov.h"
#include <stdlib.h>
#include <math.h>


/*........................................*/
double CartesianR2PoissonShafranov::x(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return (-map1_delta) * ((r/Rmax) * (r/Rmax)) - map1_kappa * (r/Rmax) * cos(theta) + (r/Rmax) * cos(theta);
}
/*........................................*/
void CartesianR2PoissonShafranov::x(std::vector<double> const& r, double theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (-map1_delta) * ((r[i]/Rmax) * (r[i]/Rmax)) - map1_kappa * (r[i]/Rmax) * cos(theta) + (r[i]/Rmax) * cos(theta);
    }
}
/*........................................*/
void CartesianR2PoissonShafranov::x(double r, std::vector<double> const& theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (-map1_delta) * ((r/Rmax) * (r/Rmax)) - map1_kappa * (r/Rmax) * cos_theta[i] + (r/Rmax) * cos_theta[i];
    }
}
/*........................................*/
double CartesianR2PoissonShafranov::y(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return map1_kappa * (r/Rmax) * sin(theta) + (r/Rmax) * sin(theta);
}
/*........................................*/
void CartesianR2PoissonShafranov::y(std::vector<double> const& r, double theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = map1_kappa * (r[i]/Rmax) * sin(theta) + (r[i]/Rmax) * sin(theta);
    }
}
/*........................................*/
void CartesianR2PoissonShafranov::y(double r, std::vector<double> const& theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = map1_kappa * (r/Rmax) * sin_theta[i] + (r/Rmax) * sin_theta[i];
    }
}
/*........................................*/
double CartesianR2PoissonShafranov::J_rr(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta))/Rmax;
}
/*........................................*/
void CartesianR2PoissonShafranov::J_rr(std::vector<double> const& r, double theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = ((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta))/Rmax;
    }
}
/*........................................*/
void CartesianR2PoissonShafranov::J_rr(double r, std::vector<double> const& theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i])/Rmax;
    }
}
/*........................................*/
double CartesianR2PoissonShafranov::J_rt(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return (r/Rmax) * (map1_kappa * sin(theta) - sin(theta));
}
/*........................................*/
void CartesianR2PoissonShafranov::J_rt(std::vector<double> const& r, double theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (r[i]/Rmax) * (map1_kappa * sin(theta) - sin(theta));
    }
}
/*........................................*/
void CartesianR2PoissonShafranov::J_rt(double r, std::vector<double> const& theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (r/Rmax) * (map1_kappa * sin_theta[i] - sin_theta[i]);
    }
}
/*........................................*/
double CartesianR2PoissonShafranov::J_tr(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return ((map1_kappa + 1.0) * sin(theta))/Rmax;
}
/*........................................*/
void CartesianR2PoissonShafranov::J_tr(std::vector<double> const& r, double theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = ((map1_kappa + 1.0) * sin(theta))/Rmax;
    }
}
/*........................................*/
void CartesianR2PoissonShafranov::J_tr(double r, std::vector<double> const& theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = ((map1_kappa + 1.0) * sin_theta[i])/Rmax;
    }
}
/*........................................*/
double CartesianR2PoissonShafranov::J_tt(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return (r/Rmax) * (map1_kappa * cos(theta) + cos(theta));
}
/*........................................*/
void CartesianR2PoissonShafranov::J_tt(std::vector<double> const& r, double theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (r[i]/Rmax) * (map1_kappa * cos(theta) + cos(theta));
    }
}
/*........................................*/
void CartesianR2PoissonShafranov::J_tt(double r, std::vector<double> const& theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (r/Rmax) * (map1_kappa * cos_theta[i] + cos_theta[i]);
    }
}
/*........................................*/
double CartesianR2PoissonShafranov::J_xs(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return (-cos(theta)) / (2.0 * map1_delta * (r/Rmax) * cos(theta) + map1_kappa * pow(sin(theta), 2.0) + map1_kappa * pow(cos(theta), 2.0) - pow(sin(theta), 2.0) - pow(cos(theta), 2.0));
}
/*........................................*/
void CartesianR2PoissonShafranov::J_xs(std::vector<double> const& r, double theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (-cos(theta)) / (2.0 * map1_delta * (r[i]/Rmax) * cos(theta) + map1_kappa * pow(sin(theta), 2.0) + map1_kappa * pow(cos(theta), 2.0) - pow(sin(theta), 2.0) - pow(cos(theta), 2.0));
    }
}
/*........................................*/
void CartesianR2PoissonShafranov::J_xs(double r, std::vector<double> const& theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (-cos_theta[i]) / (2.0 * map1_delta * (r/Rmax) * cos_theta[i] + map1_kappa * pow(sin_theta[i], 2.0) + map1_kappa * pow(cos_theta[i], 2.0) - pow(sin_theta[i], 2.0) - pow(cos_theta[i], 2.0));
    }
}
/*........................................*/
double CartesianR2PoissonShafranov::J_xt(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return (map1_kappa * sin(theta) - sin(theta)) / (2.0 * map1_delta * map1_kappa * (r/Rmax) * cos(theta) + 2.0 * map1_delta * (r/Rmax) * cos(theta) + map1_kappa * map1_kappa * pow(sin(theta), 2.0) + map1_kappa * map1_kappa * pow(cos(theta), 2.0) - pow(sin(theta), 2.0) - pow(cos(theta), 2.0));
}
/*........................................*/
void CartesianR2PoissonShafranov::J_xt(std::vector<double> const& r, double theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (map1_kappa * sin(theta) - sin(theta)) / (2.0 * map1_delta * map1_kappa * (r[i]/Rmax) * cos(theta) + 2.0 * map1_delta * (r[i]/Rmax) * cos(theta) + map1_kappa * map1_kappa * pow(sin(theta), 2.0) + map1_kappa * map1_kappa * pow(cos(theta), 2.0) - pow(sin(theta), 2.0) - pow(cos(theta), 2.0));
    }
}
/*........................................*/
void CartesianR2PoissonShafranov::J_xt(double r, std::vector<double> const& theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (map1_kappa * sin_theta[i] - sin_theta[i]) / (2.0 * map1_delta * map1_kappa * (r/Rmax) * cos_theta[i] + 2.0 * map1_delta * (r/Rmax) * cos_theta[i] + map1_kappa * map1_kappa * pow(sin_theta[i], 2.0) + map1_kappa * map1_kappa * pow(cos_theta[i], 2.0) - pow(sin_theta[i], 2.0) - pow(cos_theta[i], 2.0));
    }
}
/*........................................*/
double CartesianR2PoissonShafranov::J_ys(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return sin(theta) / (2.0 * map1_delta * ((r/Rmax) * (r/Rmax)) * cos(theta) + map1_kappa * (r/Rmax) * pow(sin(theta), 2.0) + map1_kappa * (r/Rmax) * pow(cos(theta), 2.0) - (r/Rmax) * pow(sin(theta), 2.0) - (r/Rmax) * pow(cos(theta), 2.0));
}
/*........................................*/
void CartesianR2PoissonShafranov::J_ys(std::vector<double> const& r, double theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = sin(theta) / (2.0 * map1_delta * ((r[i]/Rmax) * (r[i]/Rmax)) * cos(theta) + map1_kappa * (r[i]/Rmax) * pow(sin(theta), 2.0) + map1_kappa * (r[i]/Rmax) * pow(cos(theta), 2.0) - (r[i]/Rmax) * pow(sin(theta), 2.0) - (r[i]/Rmax) * pow(cos(theta), 2.0));
    }
}
/*........................................*/
void CartesianR2PoissonShafranov::J_ys(double r, std::vector<double> const& theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = sin_theta[i] / (2.0 * map1_delta * ((r/Rmax) * (r/Rmax)) * cos_theta[i] + map1_kappa * (r/Rmax) * pow(sin_theta[i], 2.0) + map1_kappa * (r/Rmax) * pow(cos_theta[i], 2.0) - (r/Rmax) * pow(sin_theta[i], 2.0) - (r/Rmax) * pow(cos_theta[i], 2.0));
    }
}
/*........................................*/
double CartesianR2PoissonShafranov::J_yt(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return (2.0 * map1_delta * (r/Rmax) + map1_kappa * cos(theta) - cos(theta)) / (2.0 * map1_delta * map1_kappa * ((r/Rmax) * (r/Rmax)) * cos(theta) + 2.0 * map1_delta * ((r/Rmax) * (r/Rmax)) * cos(theta) + map1_kappa * map1_kappa * (r/Rmax) * pow(sin(theta), 2.0) + map1_kappa * map1_kappa * (r/Rmax) * pow(cos(theta), 2.0) - (r/Rmax) * pow(sin(theta), 2.0) - (r/Rmax) * pow(cos(theta), 2.0));
}
/*........................................*/
void CartesianR2PoissonShafranov::J_yt(std::vector<double> const& r, double theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (2.0 * map1_delta * (r[i]/Rmax) + map1_kappa * cos(theta) - cos(theta)) / (2.0 * map1_delta * map1_kappa * ((r[i]/Rmax) * (r[i]/Rmax)) * cos(theta) + 2.0 * map1_delta * ((r[i]/Rmax) * (r[i]/Rmax)) * cos(theta) + map1_kappa * map1_kappa * (r[i]/Rmax) * pow(sin(theta), 2.0) + map1_kappa * map1_kappa * (r[i]/Rmax) * pow(cos(theta), 2.0) - (r[i]/Rmax) * pow(sin(theta), 2.0) - (r[i]/Rmax) * pow(cos(theta), 2.0));
    }
}
/*........................................*/
void CartesianR2PoissonShafranov::J_yt(double r, std::vector<double> const& theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (2.0 * map1_delta * (r/Rmax) + map1_kappa * cos_theta[i] - cos_theta[i]) / (2.0 * map1_delta * map1_kappa * ((r/Rmax) * (r/Rmax)) * cos_theta[i] + 2.0 * map1_delta * ((r/Rmax) * (r/Rmax)) * cos_theta[i] + map1_kappa * map1_kappa * (r/Rmax) * pow(sin_theta[i], 2.0) + map1_kappa * map1_kappa * (r/Rmax) * pow(cos_theta[i], 2.0) - (r/Rmax) * pow(sin_theta[i], 2.0) - (r/Rmax) * pow(cos_theta[i], 2.0));
    }
}
/*........................................*/
double CartesianR2PoissonShafranov::rho_glob(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return 0.0;
}
/*........................................*/
void CartesianR2PoissonShafranov::rho_glob(std::vector<double> const& r, double theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = 0.0;
    }
}
/*........................................*/
void CartesianR2PoissonShafranov::rho_glob(double r, std::vector<double> const& theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = 0.0;
    }
}
/*........................................*/
double CartesianR2PoissonShafranov::rho_pole(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return 0.0;
}
/*........................................*/
void CartesianR2PoissonShafranov::rho_pole(std::vector<double> const& r, double theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = 0.0;
    }
}
/*........................................*/
void CartesianR2PoissonShafranov::rho_pole(double r, std::vector<double> const& theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = 0.0;
    }
}
/*........................................*/
double CartesianR2PoissonShafranov::coeffs1(double r, double Rmax) const
{
    return 0.0;
}
/*........................................*/
void CartesianR2PoissonShafranov::coeffs1(std::vector<double> const& r, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = 0.0;
    }
}
/*........................................*/
double CartesianR2PoissonShafranov::coeffs2(double r, double Rmax) const
{
    return 0.0;
}
/*........................................*/
void CartesianR2PoissonShafranov::coeffs2(std::vector<double> const& r, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = 0.0;
    }
}
/*........................................*/
double CartesianR2PoissonShafranov::phi_exact(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return (1.0 - (r/Rmax) * (r/Rmax)) * sin(M_PI * (2.0 * map1_kappa * (r/Rmax) * sin(theta) + 2.0 * (r/Rmax) * sin(theta))) * cos(M_PI * ((-2.0) * map1_delta * ((r/Rmax) * (r/Rmax)) - 2.0 * map1_kappa * (r/Rmax) * cos(theta) + 2.0 * (r/Rmax) * cos(theta)));
}
/*........................................*/
void CartesianR2PoissonShafranov::phi_exact(std::vector<double> const& r, double theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (1.0 - (r[i]/Rmax) * (r[i]/Rmax)) * sin(M_PI * (2.0 * map1_kappa * (r[i]/Rmax) * sin(theta) + 2.0 * (r[i]/Rmax) * sin(theta))) * cos(M_PI * ((-2.0) * map1_delta * ((r[i]/Rmax) * (r[i]/Rmax)) - 2.0 * map1_kappa * (r[i]/Rmax) * cos(theta) + 2.0 * (r[i]/Rmax) * cos(theta)));
    }
}
/*........................................*/
void CartesianR2PoissonShafranov::phi_exact(double r, std::vector<double> const& theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (1.0 - (r/Rmax) * (r/Rmax)) * sin(M_PI * (2.0 * map1_kappa * (r/Rmax) * sin_theta[i] + 2.0 * (r/Rmax) * sin_theta[i])) * cos(M_PI * ((-2.0) * map1_delta * ((r/Rmax) * (r/Rmax)) - 2.0 * map1_kappa * (r/Rmax) * cos_theta[i] + 2.0 * (r/Rmax) * cos_theta[i]));
    }
}
/*........................................*/
