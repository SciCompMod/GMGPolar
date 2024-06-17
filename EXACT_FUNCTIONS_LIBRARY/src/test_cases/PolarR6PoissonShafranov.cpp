// PolarR6 simulates solution (22) of Bourne et al. https://doi.org/10.1016/j.jcp.2023.112249
#include "../../include/test_cases/PolarR6PoissonShafranov.h"

double PolarR6PoissonShafranov::x(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return (-map1_delta) * ((r/Rmax) * (r/Rmax)) - map1_kappa * (r/Rmax) * cos(theta) + (r/Rmax) * cos(theta);
}

double PolarR6PoissonShafranov::y(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return map1_kappa * (r/Rmax) * sin(theta) + (r/Rmax) * sin(theta);
}

double PolarR6PoissonShafranov::J_rr(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta))/Rmax;
}

double PolarR6PoissonShafranov::J_rt(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return (r/Rmax) * (map1_kappa * sin(theta) - sin(theta));
}

double PolarR6PoissonShafranov::J_tr(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return ((map1_kappa + 1.0) * sin(theta))/Rmax;
}

double PolarR6PoissonShafranov::J_tt(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return (r/Rmax) * (map1_kappa * cos(theta) + cos(theta));
}

double PolarR6PoissonShafranov::J_xs(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return (-cos(theta)) / (2.0 * map1_delta * (r/Rmax) * cos(theta) + map1_kappa * pow(sin(theta), 2.0) + map1_kappa * pow(cos(theta), 2.0) - pow(sin(theta), 2.0) - pow(cos(theta), 2.0));
}

double PolarR6PoissonShafranov::J_xt(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return (map1_kappa * sin(theta) - sin(theta)) / (2.0 * map1_delta * map1_kappa * (r/Rmax) * cos(theta) + 2.0 * map1_delta * (r/Rmax) * cos(theta) + map1_kappa * map1_kappa * pow(sin(theta), 2.0) + map1_kappa * map1_kappa * pow(cos(theta), 2.0) - pow(sin(theta), 2.0) - pow(cos(theta), 2.0));
}

double PolarR6PoissonShafranov::J_ys(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return sin(theta) / (2.0 * map1_delta * ((r/Rmax) * (r/Rmax)) * cos(theta) + map1_kappa * (r/Rmax) * pow(sin(theta), 2.0) + map1_kappa * (r/Rmax) * pow(cos(theta), 2.0) - (r/Rmax) * pow(sin(theta), 2.0) - (r/Rmax) * pow(cos(theta), 2.0));
}

double PolarR6PoissonShafranov::J_yt(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return (2.0 * map1_delta * (r/Rmax) + map1_kappa * cos(theta) - cos(theta)) / (2.0 * map1_delta * map1_kappa * ((r/Rmax) * (r/Rmax)) * cos(theta) + 2.0 * map1_delta * ((r/Rmax) * (r/Rmax)) * cos(theta) + map1_kappa * map1_kappa * (r/Rmax) * pow(sin(theta), 2.0) + map1_kappa * map1_kappa * (r/Rmax) * pow(cos(theta), 2.0) - (r/Rmax) * pow(sin(theta), 2.0) - (r/Rmax) * pow(cos(theta), 2.0));
}

double PolarR6PoissonShafranov::rho_glob(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return 0.0;
}

double PolarR6PoissonShafranov::rho_pole(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return 0.0;
}

double PolarR6PoissonShafranov::coeffs1(double r, double Rmax) const
{
    return 0.0;
}

double PolarR6PoissonShafranov::coeffs2(double r, double Rmax) const
{
    return 0.0;
}

double PolarR6PoissonShafranov::phi_exact(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return 0.4096 * pow((r/Rmax), 6.0) * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta);
}



double PolarR6PoissonShafranov::x(double r, double theta, double map1_kappa, double map1_delta, double Rmax, double sin_theta, double cos_theta) const
{
    return (-map1_delta) * ((r/Rmax) * (r/Rmax)) - map1_kappa * (r/Rmax) * cos_theta + (r/Rmax) * cos_theta;
}

double PolarR6PoissonShafranov::y(double r, double theta, double map1_kappa, double map1_delta, double Rmax, double sin_theta, double cos_theta) const
{
    return map1_kappa * (r/Rmax) * sin_theta + (r/Rmax) * sin_theta;
}

double PolarR6PoissonShafranov::J_rr(double r, double theta, double map1_kappa, double map1_delta, double Rmax, double sin_theta, double cos_theta) const
{
    return ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta)/Rmax;
}

double PolarR6PoissonShafranov::J_rt(double r, double theta, double map1_kappa, double map1_delta, double Rmax, double sin_theta, double cos_theta) const
{
    return (r/Rmax) * (map1_kappa * sin_theta - sin_theta);
}

double PolarR6PoissonShafranov::J_tr(double r, double theta, double map1_kappa, double map1_delta, double Rmax, double sin_theta, double cos_theta) const
{
    return ((map1_kappa + 1.0) * sin_theta)/Rmax;
}

double PolarR6PoissonShafranov::J_tt(double r, double theta, double map1_kappa, double map1_delta, double Rmax, double sin_theta, double cos_theta) const
{
    return (r/Rmax) * (map1_kappa * cos_theta + cos_theta);
}

double PolarR6PoissonShafranov::J_xs(double r, double theta, double map1_kappa, double map1_delta, double Rmax, double sin_theta, double cos_theta) const
{
    return (-cos_theta) / (2.0 * map1_delta * (r/Rmax) * cos_theta + map1_kappa * pow(sin_theta, 2.0) + map1_kappa * pow(cos_theta, 2.0) - pow(sin_theta, 2.0) - pow(cos_theta, 2.0));
}

double PolarR6PoissonShafranov::J_xt(double r, double theta, double map1_kappa, double map1_delta, double Rmax, double sin_theta, double cos_theta) const
{
    return (map1_kappa * sin_theta - sin_theta) / (2.0 * map1_delta * map1_kappa * (r/Rmax) * cos_theta + 2.0 * map1_delta * (r/Rmax) * cos_theta + map1_kappa * map1_kappa * pow(sin_theta, 2.0) + map1_kappa * map1_kappa * pow(cos_theta, 2.0) - pow(sin_theta, 2.0) - pow(cos_theta, 2.0));
}

double PolarR6PoissonShafranov::J_ys(double r, double theta, double map1_kappa, double map1_delta, double Rmax, double sin_theta, double cos_theta) const
{
    return sin_theta / (2.0 * map1_delta * ((r/Rmax) * (r/Rmax)) * cos_theta + map1_kappa * (r/Rmax) * pow(sin_theta, 2.0) + map1_kappa * (r/Rmax) * pow(cos_theta, 2.0) - (r/Rmax) * pow(sin_theta, 2.0) - (r/Rmax) * pow(cos_theta, 2.0));
}

double PolarR6PoissonShafranov::J_yt(double r, double theta, double map1_kappa, double map1_delta, double Rmax, double sin_theta, double cos_theta) const
{
    return (2.0 * map1_delta * (r/Rmax) + map1_kappa * cos_theta - cos_theta) / (2.0 * map1_delta * map1_kappa * ((r/Rmax) * (r/Rmax)) * cos_theta + 2.0 * map1_delta * ((r/Rmax) * (r/Rmax)) * cos_theta + map1_kappa * map1_kappa * (r/Rmax) * pow(sin_theta, 2.0) + map1_kappa * map1_kappa * (r/Rmax) * pow(cos_theta, 2.0) - (r/Rmax) * pow(sin_theta, 2.0) - (r/Rmax) * pow(cos_theta, 2.0));
}

double PolarR6PoissonShafranov::rho_glob(double r, double theta, double map1_kappa, double map1_delta, double Rmax, double sin_theta, double cos_theta) const
{
    return 0.0;
}

double PolarR6PoissonShafranov::rho_pole(double r, double theta, double map1_kappa, double map1_delta, double Rmax, double sin_theta, double cos_theta) const
{
    return 0.0;
}

double PolarR6PoissonShafranov::phi_exact(double r, double theta, double map1_kappa, double map1_delta, double Rmax, double sin_theta, double cos_theta) const
{
    return 0.4096 * pow((r/Rmax), 6.0) * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta);
}
