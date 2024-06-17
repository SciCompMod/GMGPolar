#include "../../include/test_cases/CartesianR2ZoniCircular.h"

double CartesianR2ZoniCircular::x(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (r/Rmax) * cos(theta);
}

double CartesianR2ZoniCircular::y(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (r/Rmax) * sin(theta);
}

double CartesianR2ZoniCircular::J_rr(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (cos(theta))/Rmax;
}

double CartesianR2ZoniCircular::J_rt(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (-(r/Rmax)) * sin(theta);
}

double CartesianR2ZoniCircular::J_tr(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (sin(theta))/Rmax;
}

double CartesianR2ZoniCircular::J_tt(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (r/Rmax) * cos(theta);
}

double CartesianR2ZoniCircular::J_xs(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (-pow(sin(theta), 2.0)) / cos(theta) + pow(cos(theta), (double)((-1)));
}

double CartesianR2ZoniCircular::J_xt(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return sin(theta);
}

double CartesianR2ZoniCircular::J_ys(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (-sin(theta)) / (r/Rmax);
}

double CartesianR2ZoniCircular::J_yt(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return cos(theta) / (r/Rmax);
}

double CartesianR2ZoniCircular::rho_glob(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (-((r/Rmax) * (10.0 * pow(tanh(10.0 * (r/Rmax) - 5.0), 2.0) - 10.0) * ((-2.0) * (r/Rmax) * sin(2.0 * M_PI * (r/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r/Rmax) * cos(theta)) + 2.0 * M_PI * (1.0 - (r/Rmax) * (r/Rmax)) * sin(theta) * cos(2.0 * M_PI * (r/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r/Rmax) * cos(theta)) - 2.0 * M_PI * (1.0 - (r/Rmax) * (r/Rmax)) * sin(2.0 * M_PI * (r/Rmax) * sin(theta)) * sin(2.0 * M_PI * (r/Rmax) * cos(theta)) * cos(theta)) * exp(-tanh(10.0 * (r/Rmax) - 5.0)) + (r/Rmax) * ((-8.0) * M_PI * (r/Rmax) * sin(theta) * cos(2.0 * M_PI * (r/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r/Rmax) * cos(theta)) + 8.0 * M_PI * (r/Rmax) * sin(2.0 * M_PI * (r/Rmax) * sin(theta)) * sin(2.0 * M_PI * (r/Rmax) * cos(theta)) * cos(theta) - 4.0 * (M_PI * M_PI) * (1.0 - (r/Rmax) * (r/Rmax)) * pow(sin(theta), 2.0) * sin(2.0 * M_PI * (r/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r/Rmax) * cos(theta)) - 8.0 * (M_PI * M_PI) * (1.0 - (r/Rmax) * (r/Rmax)) * sin(theta) * sin(2.0 * M_PI * (r/Rmax) * cos(theta)) * cos(theta) * cos(2.0 * M_PI * (r/Rmax) * sin(theta)) - 4.0 * (M_PI * M_PI) * (1.0 - (r/Rmax) * (r/Rmax)) * sin(2.0 * M_PI * (r/Rmax) * sin(theta)) * pow(cos(theta), 2.0) * cos(2.0 * M_PI * (r/Rmax) * cos(theta)) - 2.0 * sin(2.0 * M_PI * (r/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r/Rmax) * cos(theta))) * exp(-tanh(10.0 * (r/Rmax) - 5.0)) + ((-2.0) * (r/Rmax) * sin(2.0 * M_PI * (r/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r/Rmax) * cos(theta)) + 2.0 * M_PI * (1.0 - (r/Rmax) * (r/Rmax)) * sin(theta) * cos(2.0 * M_PI * (r/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r/Rmax) * cos(theta)) - 2.0 * M_PI * (1.0 - (r/Rmax) * (r/Rmax)) * sin(2.0 * M_PI * (r/Rmax) * sin(theta)) * sin(2.0 * M_PI * (r/Rmax) * cos(theta)) * cos(theta)) * exp(-tanh(10.0 * (r/Rmax) - 5.0)) + ((-4.0) * (M_PI * M_PI) * (r/Rmax) * (1.0 - (r/Rmax) * (r/Rmax)) * pow(sin(theta), 2.0) * sin(2.0 * M_PI * (r/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r/Rmax) * cos(theta)) + 8.0 * (M_PI * M_PI) * (r/Rmax) * (1.0 - (r/Rmax) * (r/Rmax)) * sin(theta) * sin(2.0 * M_PI * (r/Rmax) * cos(theta)) * cos(theta) * cos(2.0 * M_PI * (r/Rmax) * sin(theta)) - 4.0 * (M_PI * M_PI) * (r/Rmax) * (1.0 - (r/Rmax) * (r/Rmax)) * sin(2.0 * M_PI * (r/Rmax) * sin(theta)) * pow(cos(theta), 2.0) * cos(2.0 * M_PI * (r/Rmax) * cos(theta)) - 2.0 * M_PI * (1.0 - (r/Rmax) * (r/Rmax)) * sin(theta) * cos(2.0 * M_PI * (r/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r/Rmax) * cos(theta)) + 2.0 * M_PI * (1.0 - (r/Rmax) * (r/Rmax)) * sin(2.0 * M_PI * (r/Rmax) * sin(theta)) * sin(2.0 * M_PI * (r/Rmax) * cos(theta)) * cos(theta)) * exp(-tanh(10.0 * (r/Rmax) - 5.0)))) / (r/Rmax);
}

double CartesianR2ZoniCircular::rho_pole(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return 0.0;
}

double CartesianR2ZoniCircular::coeffs1(double r, double Rmax) const
{
    return exp(-tanh(10.0 * (r/Rmax) - 5.0));
}

double CartesianR2ZoniCircular::coeffs2(double r, double Rmax) const
{
    return 0.0;
}

double CartesianR2ZoniCircular::phi_exact(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (1.0 - (r/Rmax) * (r/Rmax)) * sin(2.0 * M_PI * (r/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r/Rmax) * cos(theta));
}



double CartesianR2ZoniCircular::x(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return (r/Rmax) * cos_theta;
}

double CartesianR2ZoniCircular::y(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return (r/Rmax) * sin_theta;
}

double CartesianR2ZoniCircular::J_rr(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return (cos_theta)/Rmax;
}

double CartesianR2ZoniCircular::J_rt(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return (-(r/Rmax)) * sin_theta;
}

double CartesianR2ZoniCircular::J_tr(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return (sin_theta)/Rmax;
}

double CartesianR2ZoniCircular::J_tt(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return (r/Rmax) * cos_theta;
}

double CartesianR2ZoniCircular::J_xs(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return (-pow(sin_theta, 2.0)) / cos_theta + pow(cos_theta, (double)((-1)));
}

double CartesianR2ZoniCircular::J_xt(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return sin_theta;
}

double CartesianR2ZoniCircular::J_ys(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return (-sin_theta) / (r/Rmax);
}

double CartesianR2ZoniCircular::J_yt(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return cos_theta / (r/Rmax);
}

double CartesianR2ZoniCircular::rho_glob(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return (-((r/Rmax) * (10.0 * pow(tanh(10.0 * (r/Rmax) - 5.0), 2.0) - 10.0) * ((-2.0) * (r/Rmax) * sin(2.0 * M_PI * (r/Rmax) * sin_theta) * cos(2.0 * M_PI * (r/Rmax) * cos_theta) + 2.0 * M_PI * (1.0 - (r/Rmax) * (r/Rmax)) * sin_theta * cos(2.0 * M_PI * (r/Rmax) * sin_theta) * cos(2.0 * M_PI * (r/Rmax) * cos_theta) - 2.0 * M_PI * (1.0 - (r/Rmax) * (r/Rmax)) * sin(2.0 * M_PI * (r/Rmax) * sin_theta) * sin(2.0 * M_PI * (r/Rmax) * cos_theta) * cos_theta) * exp(-tanh(10.0 * (r/Rmax) - 5.0)) + (r/Rmax) * ((-8.0) * M_PI * (r/Rmax) * sin_theta * cos(2.0 * M_PI * (r/Rmax) * sin_theta) * cos(2.0 * M_PI * (r/Rmax) * cos_theta) + 8.0 * M_PI * (r/Rmax) * sin(2.0 * M_PI * (r/Rmax) * sin_theta) * sin(2.0 * M_PI * (r/Rmax) * cos_theta) * cos_theta - 4.0 * (M_PI * M_PI) * (1.0 - (r/Rmax) * (r/Rmax)) * pow(sin_theta, 2.0) * sin(2.0 * M_PI * (r/Rmax) * sin_theta) * cos(2.0 * M_PI * (r/Rmax) * cos_theta) - 8.0 * (M_PI * M_PI) * (1.0 - (r/Rmax) * (r/Rmax)) * sin_theta * sin(2.0 * M_PI * (r/Rmax) * cos_theta) * cos_theta * cos(2.0 * M_PI * (r/Rmax) * sin_theta) - 4.0 * (M_PI * M_PI) * (1.0 - (r/Rmax) * (r/Rmax)) * sin(2.0 * M_PI * (r/Rmax) * sin_theta) * pow(cos_theta, 2.0) * cos(2.0 * M_PI * (r/Rmax) * cos_theta) - 2.0 * sin(2.0 * M_PI * (r/Rmax) * sin_theta) * cos(2.0 * M_PI * (r/Rmax) * cos_theta)) * exp(-tanh(10.0 * (r/Rmax) - 5.0)) + ((-2.0) * (r/Rmax) * sin(2.0 * M_PI * (r/Rmax) * sin_theta) * cos(2.0 * M_PI * (r/Rmax) * cos_theta) + 2.0 * M_PI * (1.0 - (r/Rmax) * (r/Rmax)) * sin_theta * cos(2.0 * M_PI * (r/Rmax) * sin_theta) * cos(2.0 * M_PI * (r/Rmax) * cos_theta) - 2.0 * M_PI * (1.0 - (r/Rmax) * (r/Rmax)) * sin(2.0 * M_PI * (r/Rmax) * sin_theta) * sin(2.0 * M_PI * (r/Rmax) * cos_theta) * cos_theta) * exp(-tanh(10.0 * (r/Rmax) - 5.0)) + ((-4.0) * (M_PI * M_PI) * (r/Rmax) * (1.0 - (r/Rmax) * (r/Rmax)) * pow(sin_theta, 2.0) * sin(2.0 * M_PI * (r/Rmax) * sin_theta) * cos(2.0 * M_PI * (r/Rmax) * cos_theta) + 8.0 * (M_PI * M_PI) * (r/Rmax) * (1.0 - (r/Rmax) * (r/Rmax)) * sin_theta * sin(2.0 * M_PI * (r/Rmax) * cos_theta) * cos_theta * cos(2.0 * M_PI * (r/Rmax) * sin_theta) - 4.0 * (M_PI * M_PI) * (r/Rmax) * (1.0 - (r/Rmax) * (r/Rmax)) * sin(2.0 * M_PI * (r/Rmax) * sin_theta) * pow(cos_theta, 2.0) * cos(2.0 * M_PI * (r/Rmax) * cos_theta) - 2.0 * M_PI * (1.0 - (r/Rmax) * (r/Rmax)) * sin_theta * cos(2.0 * M_PI * (r/Rmax) * sin_theta) * cos(2.0 * M_PI * (r/Rmax) * cos_theta) + 2.0 * M_PI * (1.0 - (r/Rmax) * (r/Rmax)) * sin(2.0 * M_PI * (r/Rmax) * sin_theta) * sin(2.0 * M_PI * (r/Rmax) * cos_theta) * cos_theta) * exp(-tanh(10.0 * (r/Rmax) - 5.0)))) / (r/Rmax);
}

double CartesianR2ZoniCircular::rho_pole(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return 0.0;
}

double CartesianR2ZoniCircular::phi_exact(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return (1.0 - (r/Rmax) * (r/Rmax)) * sin(2.0 * M_PI * (r/Rmax) * sin_theta) * cos(2.0 * M_PI * (r/Rmax) * cos_theta);
}
