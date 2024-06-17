// PolarR6 simulates solution (22) of Bourne et al. https://doi.org/10.1016/j.jcp.2023.112249
#include "../../include/test_cases/PolarR6PoissonCircular.h"

double PolarR6PoissonCircular::x(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (r/Rmax) * cos(theta);
}

double PolarR6PoissonCircular::y(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (r/Rmax) * sin(theta);
}

double PolarR6PoissonCircular::J_rr(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (cos(theta))/Rmax;
}

double PolarR6PoissonCircular::J_rt(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (-(r/Rmax)) * sin(theta);
}

double PolarR6PoissonCircular::J_tr(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (sin(theta))/Rmax;
}

double PolarR6PoissonCircular::J_tt(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (r/Rmax) * cos(theta);
}

double PolarR6PoissonCircular::J_xs(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (-pow(sin(theta), 2.0)) / cos(theta) + pow(cos(theta), (double)((-1)));
}

double PolarR6PoissonCircular::J_xt(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return sin(theta);
}

double PolarR6PoissonCircular::J_ys(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (-sin(theta)) / (r/Rmax);
}

double PolarR6PoissonCircular::J_yt(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return cos(theta) / (r/Rmax);
}

double PolarR6PoissonCircular::rho_glob(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return 0.0;
}

double PolarR6PoissonCircular::rho_pole(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return 0.0;
}

double PolarR6PoissonCircular::coeffs1(double r, double Rmax) const
{
    return 0.0;
}

double PolarR6PoissonCircular::coeffs2(double r, double Rmax) const
{
    return 0.0;
}

double PolarR6PoissonCircular::phi_exact(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return 0.4096 * pow((r/Rmax), 6.0) * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta);
}



double PolarR6PoissonCircular::x(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return (r/Rmax) * cos_theta;
}

double PolarR6PoissonCircular::y(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return (r/Rmax) * sin_theta;
}

double PolarR6PoissonCircular::J_rr(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return (cos_theta)/Rmax;
}

double PolarR6PoissonCircular::J_rt(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return (-(r/Rmax)) * sin_theta;
}

double PolarR6PoissonCircular::J_tr(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return (sin_theta)/Rmax;
}

double PolarR6PoissonCircular::J_tt(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return (r/Rmax) * cos_theta;
}

double PolarR6PoissonCircular::J_xs(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return (-pow(sin_theta, 2.0)) / cos_theta + pow(cos_theta, (double)((-1)));
}

double PolarR6PoissonCircular::J_xt(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return sin_theta;
}

double PolarR6PoissonCircular::J_ys(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return (-sin_theta) / (r/Rmax);
}

double PolarR6PoissonCircular::J_yt(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return cos_theta / (r/Rmax);
}

double PolarR6PoissonCircular::rho_glob(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return 0.0;
}

double PolarR6PoissonCircular::rho_pole(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return 0.0;
}

double PolarR6PoissonCircular::phi_exact(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return 0.4096 * pow((r/Rmax), 6.0) * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta);
}
