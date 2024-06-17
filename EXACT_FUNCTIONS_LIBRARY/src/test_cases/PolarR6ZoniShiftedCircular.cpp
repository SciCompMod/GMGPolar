// PolarR6 simulates solution (22) of Bourne et al. https://doi.org/10.1016/j.jcp.2023.112249
#include "../../include/test_cases/PolarR6ZoniShiftedCircular.h"

double PolarR6ZoniShiftedCircular::x(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (r/Rmax) * cos(theta);
}

double PolarR6ZoniShiftedCircular::y(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (r/Rmax) * sin(theta);
}

double PolarR6ZoniShiftedCircular::J_rr(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (cos(theta))/Rmax;
}

double PolarR6ZoniShiftedCircular::J_rt(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (-(r/Rmax)) * sin(theta);
}

double PolarR6ZoniShiftedCircular::J_tr(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (sin(theta))/Rmax;
}

double PolarR6ZoniShiftedCircular::J_tt(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (r/Rmax) * cos(theta);
}

double PolarR6ZoniShiftedCircular::J_xs(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (-pow(sin(theta), 2.0)) / cos(theta) + pow(cos(theta), (double)((-1)));
}

double PolarR6ZoniShiftedCircular::J_xt(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return sin(theta);
}

double PolarR6ZoniShiftedCircular::J_ys(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (-sin(theta)) / (r/Rmax);
}

double PolarR6ZoniShiftedCircular::J_yt(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return cos(theta) / (r/Rmax);
}

double PolarR6ZoniShiftedCircular::rho_glob(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (-pow((r/Rmax), 4.0)) * ((r/Rmax) * (12.288 * (r/Rmax) * pow(((r/Rmax) - 1.0), 4.0) * cos(11.0 * theta) + 17.2032 * pow(((r/Rmax) - 1.0), 5.0) * cos(11.0 * theta)) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) + (r/Rmax) * (2.4576 * (r/Rmax) * pow(((r/Rmax) - 1.0), 5.0) * cos(11.0 * theta) + 2.4576 * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta)) * (20.0 * pow(tanh(20.0 * (r/Rmax) - 14.0), 2.0) - 20.0) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) - 49.5616 * pow(((r/Rmax) - 1.0), 6.0) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) * cos(11.0 * theta) + 6.0 * (2.4576 * (r/Rmax) * pow(((r/Rmax) - 1.0), 5.0) * cos(11.0 * theta) + 2.4576 * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta)) * exp(-tanh(20.0 * (r/Rmax) - 14.0)));
}

double PolarR6ZoniShiftedCircular::rho_pole(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return 0.0;
}

double PolarR6ZoniShiftedCircular::coeffs1(double r, double Rmax) const
{ // With Rmax=1, equals alpha(r) from equation (18) of Bourne et al. https://doi.org/10.1016/j.jcp.2023.112249
    return exp(-tanh(20.0 * (r/Rmax) - 14.0));
}

double PolarR6ZoniShiftedCircular::coeffs2(double r, double Rmax) const
{
    return 0.0;
}

double PolarR6ZoniShiftedCircular::phi_exact(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return 0.4096 * pow((r/Rmax), 6.0) * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta);
}



double PolarR6ZoniShiftedCircular::x(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return (r/Rmax) * cos_theta;
}

double PolarR6ZoniShiftedCircular::y(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return (r/Rmax) * sin_theta;
}

double PolarR6ZoniShiftedCircular::J_rr(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return (cos_theta)/Rmax;
}

double PolarR6ZoniShiftedCircular::J_rt(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return (-(r/Rmax)) * sin_theta;
}

double PolarR6ZoniShiftedCircular::J_tr(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return (sin_theta)/Rmax;
}

double PolarR6ZoniShiftedCircular::J_tt(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return (r/Rmax) * cos_theta;
}

double PolarR6ZoniShiftedCircular::J_xs(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return (-pow(sin_theta, 2.0)) / cos_theta + pow(cos_theta, (double)((-1)));
}

double PolarR6ZoniShiftedCircular::J_xt(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return sin_theta;
}

double PolarR6ZoniShiftedCircular::J_ys(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return (-sin_theta) / (r/Rmax);
}

double PolarR6ZoniShiftedCircular::J_yt(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return cos_theta / (r/Rmax);
}

double PolarR6ZoniShiftedCircular::rho_glob(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return (-pow((r/Rmax), 4.0)) * ((r/Rmax) * (12.288 * (r/Rmax) * pow(((r/Rmax) - 1.0), 4.0) * cos(11.0 * theta) + 17.2032 * pow(((r/Rmax) - 1.0), 5.0) * cos(11.0 * theta)) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) + (r/Rmax) * (2.4576 * (r/Rmax) * pow(((r/Rmax) - 1.0), 5.0) * cos(11.0 * theta) + 2.4576 * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta)) * (20.0 * pow(tanh(20.0 * (r/Rmax) - 14.0), 2.0) - 20.0) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) - 49.5616 * pow(((r/Rmax) - 1.0), 6.0) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) * cos(11.0 * theta) + 6.0 * (2.4576 * (r/Rmax) * pow(((r/Rmax) - 1.0), 5.0) * cos(11.0 * theta) + 2.4576 * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta)) * exp(-tanh(20.0 * (r/Rmax) - 14.0)));
}

double PolarR6ZoniShiftedCircular::rho_pole(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return 0.0;
}

double PolarR6ZoniShiftedCircular::phi_exact(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return 0.4096 * pow((r/Rmax), 6.0) * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta);
}