// PolarR6 simulates solution (22) of Bourne et al. https://doi.org/10.1016/j.jcp.2023.112249
#include "../../include/test_cases/PolarR6GyroZoniShiftedShafranov.h"

double PolarR6GyroZoniShiftedShafranov::x(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return (-map1_delta) * ((r/Rmax) * (r/Rmax)) - map1_kappa * (r/Rmax) * cos(theta) + (r/Rmax) * cos(theta);
}

double PolarR6GyroZoniShiftedShafranov::y(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return map1_kappa * (r/Rmax) * sin(theta) + (r/Rmax) * sin(theta);
}

double PolarR6GyroZoniShiftedShafranov::J_rr(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta))/Rmax;
}

double PolarR6GyroZoniShiftedShafranov::J_rt(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return (r/Rmax) * (map1_kappa * sin(theta) - sin(theta));
}

double PolarR6GyroZoniShiftedShafranov::J_tr(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return ((map1_kappa + 1.0) * sin(theta))/Rmax;
}

double PolarR6GyroZoniShiftedShafranov::J_tt(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return (r/Rmax) * (map1_kappa * cos(theta) + cos(theta));
}

double PolarR6GyroZoniShiftedShafranov::J_xs(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return (-cos(theta)) / (2.0 * map1_delta * (r/Rmax) * cos(theta) + map1_kappa * pow(sin(theta), 2.0) + map1_kappa * pow(cos(theta), 2.0) - pow(sin(theta), 2.0) - pow(cos(theta), 2.0));
}

double PolarR6GyroZoniShiftedShafranov::J_xt(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return (map1_kappa * sin(theta) - sin(theta)) / (2.0 * map1_delta * map1_kappa * (r/Rmax) * cos(theta) + 2.0 * map1_delta * (r/Rmax) * cos(theta) + map1_kappa * map1_kappa * pow(sin(theta), 2.0) + map1_kappa * map1_kappa * pow(cos(theta), 2.0) - pow(sin(theta), 2.0) - pow(cos(theta), 2.0));
}

double PolarR6GyroZoniShiftedShafranov::J_ys(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return sin(theta) / (2.0 * map1_delta * ((r/Rmax) * (r/Rmax)) * cos(theta) + map1_kappa * (r/Rmax) * pow(sin(theta), 2.0) + map1_kappa * (r/Rmax) * pow(cos(theta), 2.0) - (r/Rmax) * pow(sin(theta), 2.0) - (r/Rmax) * pow(cos(theta), 2.0));
}

double PolarR6GyroZoniShiftedShafranov::J_yt(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return (2.0 * map1_delta * (r/Rmax) + map1_kappa * cos(theta) - cos(theta)) / (2.0 * map1_delta * map1_kappa * ((r/Rmax) * (r/Rmax)) * cos(theta) + 2.0 * map1_delta * ((r/Rmax) * (r/Rmax)) * cos(theta) + map1_kappa * map1_kappa * (r/Rmax) * pow(sin(theta), 2.0) + map1_kappa * map1_kappa * (r/Rmax) * pow(cos(theta), 2.0) - (r/Rmax) * pow(sin(theta), 2.0) - (r/Rmax) * pow(cos(theta), 2.0));
}

double PolarR6GyroZoniShiftedShafranov::rho_glob(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return 0.4096 * pow((r/Rmax), 6.0) * pow(((r/Rmax) - 1.0), 6.0) * exp(tanh(20.0 * (r/Rmax) - 14.0)) * cos(11.0 * theta) - pow((r/Rmax), 4.0) * ((-9.0112) * map1_delta * (r/Rmax) * (map1_kappa - 1.0) * pow(((r/Rmax) - 1.0), 6.0) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) * sin(theta) * sin(11.0 * theta) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)) + 4.5056 * (r/Rmax) * pow(((r/Rmax) - 1.0), 6.0) * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)) * (20.0 * pow(tanh(20.0 * (r/Rmax) - 14.0), 2.0) - 20.0) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) * sin(11.0 * theta) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)) + 4.5056 * (r/Rmax) * pow(((r/Rmax) - 1.0), 6.0) * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)) * ((-2.0) * map1_delta * (map1_kappa - 1.0) * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)) * sin(theta) + 2.0 * map1_delta * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0)) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) * sin(11.0 * theta) / pow(((pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)), (3.0 / 2.0)) + 27.0336 * (r/Rmax) * pow(((r/Rmax) - 1.0), 5.0) * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) * sin(11.0 * theta) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)) + (r/Rmax) * (12.288 * (r/Rmax) * pow(((r/Rmax) - 1.0), 4.0) * cos(11.0 * theta) + 17.2032 * pow(((r/Rmax) - 1.0), 5.0) * cos(11.0 * theta)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)) + (r/Rmax) * (2.4576 * (r/Rmax) * pow(((r/Rmax) - 1.0), 5.0) * cos(11.0 * theta) + 2.4576 * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta)) * (20.0 * pow(tanh(20.0 * (r/Rmax) - 14.0), 2.0) - 20.0) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)) + (r/Rmax) * (2.4576 * (r/Rmax) * pow(((r/Rmax) - 1.0), 5.0) * cos(11.0 * theta) + 2.4576 * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta)) * ((-2.0) * map1_delta * (map1_kappa - 1.0) * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)) * sin(theta) + 2.0 * map1_delta * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) / pow(((pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)), (3.0 / 2.0)) - 49.5616 * pow(((r/Rmax) - 1.0), 6.0) * (pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) * cos(11.0 * theta) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)) - 4.5056 * pow(((r/Rmax) - 1.0), 6.0) * (pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (4.0 * map1_kappa * (pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * sin(theta) * cos(theta) - 1.0 / 2.0 * (pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta) + (2.0 * map1_kappa - 2.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) + 1.0 / 2.0 * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)) * (2.0 * pow((map1_kappa - 1.0), 2.0) * pow(sin(theta), 2.0) + 2.0 * (map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * cos(theta) + 2.0 * pow((map1_kappa + 1.0), 2.0) * cos(2.0 * theta))) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) * sin(11.0 * theta) / pow(((pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)), (3.0 / 2.0)) - 4.5056 * pow(((r/Rmax) - 1.0), 6.0) * (pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta) + (2.0 * map1_kappa - 2.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta)) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) * sin(11.0 * theta) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)) + 27.0336 * pow(((r/Rmax) - 1.0), 6.0) * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) * sin(11.0 * theta) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)) - ((-27.0336) * (r/Rmax) * pow(((r/Rmax) - 1.0), 5.0) * sin(11.0 * theta) - 27.0336 * pow(((r/Rmax) - 1.0), 6.0) * sin(11.0 * theta)) * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)) - (2.4576 * (r/Rmax) * pow(((r/Rmax) - 1.0), 5.0) * cos(11.0 * theta) + 2.4576 * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta)) * (pow((map1_kappa - 1.0), 2.0) * pow(sin(theta), 2.0) + (map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * cos(theta) + pow((map1_kappa + 1.0), 2.0) * cos(2.0 * theta)) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)) + 6.0 * (2.4576 * (r/Rmax) * pow(((r/Rmax) - 1.0), 5.0) * cos(11.0 * theta) + 2.4576 * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)) - (2.4576 * (r/Rmax) * pow(((r/Rmax) - 1.0), 5.0) * cos(11.0 * theta) + 2.4576 * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta)) * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)) * (4.0 * map1_kappa * (pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * sin(theta) * cos(theta) - 1.0 / 2.0 * (pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta) + (2.0 * map1_kappa - 2.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) + 1.0 / 2.0 * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)) * (2.0 * pow((map1_kappa - 1.0), 2.0) * pow(sin(theta), 2.0) + 2.0 * (map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * cos(theta) + 2.0 * pow((map1_kappa + 1.0), 2.0) * cos(2.0 * theta))) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) / pow(((pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)), (3.0 / 2.0))) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0));
}

double PolarR6GyroZoniShiftedShafranov::rho_pole(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return 0.0;
}

double PolarR6GyroZoniShiftedShafranov::coeffs1(double r, double Rmax) const
{ // With Rmax=1, equals alpha(r) from equation (18) of Bourne et al. https://doi.org/10.1016/j.jcp.2023.112249
    return exp(-tanh(20.0 * (r/Rmax) - 14.0));
}

double PolarR6GyroZoniShiftedShafranov::coeffs2(double r, double Rmax) const
{
    return exp(tanh(20.0 * (r/Rmax) - 14.0));
}

double PolarR6GyroZoniShiftedShafranov::phi_exact(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return 0.4096 * pow((r/Rmax), 6.0) * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta);
}



double PolarR6GyroZoniShiftedShafranov::x(double r, double theta, double map1_kappa, double map1_delta, double Rmax, double sin_theta, double cos_theta) const
{
    return (-map1_delta) * ((r/Rmax) * (r/Rmax)) - map1_kappa * (r/Rmax) * cos_theta + (r/Rmax) * cos_theta;
}

double PolarR6GyroZoniShiftedShafranov::y(double r, double theta, double map1_kappa, double map1_delta, double Rmax, double sin_theta, double cos_theta) const
{
    return map1_kappa * (r/Rmax) * sin_theta + (r/Rmax) * sin_theta;
}

double PolarR6GyroZoniShiftedShafranov::J_rr(double r, double theta, double map1_kappa, double map1_delta, double Rmax, double sin_theta, double cos_theta) const
{
    return ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta)/Rmax;
}

double PolarR6GyroZoniShiftedShafranov::J_rt(double r, double theta, double map1_kappa, double map1_delta, double Rmax, double sin_theta, double cos_theta) const
{
    return (r/Rmax) * (map1_kappa * sin_theta - sin_theta);
}

double PolarR6GyroZoniShiftedShafranov::J_tr(double r, double theta, double map1_kappa, double map1_delta, double Rmax, double sin_theta, double cos_theta) const
{
    return ((map1_kappa + 1.0) * sin_theta)/Rmax;
}

double PolarR6GyroZoniShiftedShafranov::J_tt(double r, double theta, double map1_kappa, double map1_delta, double Rmax, double sin_theta, double cos_theta) const
{
    return (r/Rmax) * (map1_kappa * cos_theta + cos_theta);
}

double PolarR6GyroZoniShiftedShafranov::J_xs(double r, double theta, double map1_kappa, double map1_delta, double Rmax, double sin_theta, double cos_theta) const
{
    return (-cos_theta) / (2.0 * map1_delta * (r/Rmax) * cos_theta + map1_kappa * pow(sin_theta, 2.0) + map1_kappa * pow(cos_theta, 2.0) - pow(sin_theta, 2.0) - pow(cos_theta, 2.0));
}

double PolarR6GyroZoniShiftedShafranov::J_xt(double r, double theta, double map1_kappa, double map1_delta, double Rmax, double sin_theta, double cos_theta) const
{
    return (map1_kappa * sin_theta - sin_theta) / (2.0 * map1_delta * map1_kappa * (r/Rmax) * cos_theta + 2.0 * map1_delta * (r/Rmax) * cos_theta + map1_kappa * map1_kappa * pow(sin_theta, 2.0) + map1_kappa * map1_kappa * pow(cos_theta, 2.0) - pow(sin_theta, 2.0) - pow(cos_theta, 2.0));
}

double PolarR6GyroZoniShiftedShafranov::J_ys(double r, double theta, double map1_kappa, double map1_delta, double Rmax, double sin_theta, double cos_theta) const
{
    return sin_theta / (2.0 * map1_delta * ((r/Rmax) * (r/Rmax)) * cos_theta + map1_kappa * (r/Rmax) * pow(sin_theta, 2.0) + map1_kappa * (r/Rmax) * pow(cos_theta, 2.0) - (r/Rmax) * pow(sin_theta, 2.0) - (r/Rmax) * pow(cos_theta, 2.0));
}

double PolarR6GyroZoniShiftedShafranov::J_yt(double r, double theta, double map1_kappa, double map1_delta, double Rmax, double sin_theta, double cos_theta) const
{
    return (2.0 * map1_delta * (r/Rmax) + map1_kappa * cos_theta - cos_theta) / (2.0 * map1_delta * map1_kappa * ((r/Rmax) * (r/Rmax)) * cos_theta + 2.0 * map1_delta * ((r/Rmax) * (r/Rmax)) * cos_theta + map1_kappa * map1_kappa * (r/Rmax) * pow(sin_theta, 2.0) + map1_kappa * map1_kappa * (r/Rmax) * pow(cos_theta, 2.0) - (r/Rmax) * pow(sin_theta, 2.0) - (r/Rmax) * pow(cos_theta, 2.0));
}

double PolarR6GyroZoniShiftedShafranov::rho_glob(double r, double theta, double map1_kappa, double map1_delta, double Rmax, double sin_theta, double cos_theta) const
{
    return 0.4096 * pow((r/Rmax), 6.0) * pow(((r/Rmax) - 1.0), 6.0) * exp(tanh(20.0 * (r/Rmax) - 14.0)) * cos(11.0 * theta) - pow((r/Rmax), 4.0) * ((-9.0112) * map1_delta * (r/Rmax) * (map1_kappa - 1.0) * pow(((r/Rmax) - 1.0), 6.0) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) * sin_theta * sin(11.0 * theta) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin_theta, 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta) * sin_theta + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)) + 4.5056 * (r/Rmax) * pow(((r/Rmax) - 1.0), 6.0) * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta) * sin_theta + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)) * (20.0 * pow(tanh(20.0 * (r/Rmax) - 14.0), 2.0) - 20.0) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) * sin(11.0 * theta) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin_theta, 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta) * sin_theta + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)) + 4.5056 * (r/Rmax) * pow(((r/Rmax) - 1.0), 6.0) * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta) * sin_theta + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)) * ((-2.0) * map1_delta * (map1_kappa - 1.0) * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta) * sin_theta + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)) * sin_theta + 2.0 * map1_delta * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin_theta, 2.0) + 2.0 * map1_kappa + 1.0)) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) * sin(11.0 * theta) / pow(((pow((map1_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin_theta, 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta) * sin_theta + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)), (3.0 / 2.0)) + 27.0336 * (r/Rmax) * pow(((r/Rmax) - 1.0), 5.0) * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta) * sin_theta + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) * sin(11.0 * theta) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin_theta, 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta) * sin_theta + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)) + (r/Rmax) * (12.288 * (r/Rmax) * pow(((r/Rmax) - 1.0), 4.0) * cos(11.0 * theta) + 17.2032 * pow(((r/Rmax) - 1.0), 5.0) * cos(11.0 * theta)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin_theta, 2.0) + 2.0 * map1_kappa + 1.0) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin_theta, 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta) * sin_theta + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)) + (r/Rmax) * (2.4576 * (r/Rmax) * pow(((r/Rmax) - 1.0), 5.0) * cos(11.0 * theta) + 2.4576 * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta)) * (20.0 * pow(tanh(20.0 * (r/Rmax) - 14.0), 2.0) - 20.0) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin_theta, 2.0) + 2.0 * map1_kappa + 1.0) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin_theta, 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta) * sin_theta + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)) + (r/Rmax) * (2.4576 * (r/Rmax) * pow(((r/Rmax) - 1.0), 5.0) * cos(11.0 * theta) + 2.4576 * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta)) * ((-2.0) * map1_delta * (map1_kappa - 1.0) * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta) * sin_theta + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)) * sin_theta + 2.0 * map1_delta * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin_theta, 2.0) + 2.0 * map1_kappa + 1.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin_theta, 2.0) + 2.0 * map1_kappa + 1.0) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) / pow(((pow((map1_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin_theta, 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta) * sin_theta + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)), (3.0 / 2.0)) - 49.5616 * pow(((r/Rmax) - 1.0), 6.0) * (pow((map1_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta), 2.0)) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) * cos(11.0 * theta) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin_theta, 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta) * sin_theta + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)) - 4.5056 * pow(((r/Rmax) - 1.0), 6.0) * (pow((map1_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta), 2.0)) * (4.0 * map1_kappa * (pow((map1_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta), 2.0)) * sin_theta * cos_theta - 1.0 / 2.0 * (pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta) + (2.0 * map1_kappa - 2.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta) * sin_theta) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin_theta, 2.0) + 2.0 * map1_kappa + 1.0) + 1.0 / 2.0 * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta) * sin_theta + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)) * (2.0 * pow((map1_kappa - 1.0), 2.0) * pow(sin_theta, 2.0) + 2.0 * (map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta) * cos_theta + 2.0 * pow((map1_kappa + 1.0), 2.0) * cos(2.0 * theta))) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) * sin(11.0 * theta) / pow(((pow((map1_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin_theta, 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta) * sin_theta + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)), (3.0 / 2.0)) - 4.5056 * pow(((r/Rmax) - 1.0), 6.0) * (pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta) + (2.0 * map1_kappa - 2.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta) * sin_theta) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) * sin(11.0 * theta) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin_theta, 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta) * sin_theta + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)) + 27.0336 * pow(((r/Rmax) - 1.0), 6.0) * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta) * sin_theta + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) * sin(11.0 * theta) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin_theta, 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta) * sin_theta + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)) - ((-27.0336) * (r/Rmax) * pow(((r/Rmax) - 1.0), 5.0) * sin(11.0 * theta) - 27.0336 * pow(((r/Rmax) - 1.0), 6.0) * sin(11.0 * theta)) * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta) * sin_theta + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin_theta, 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta) * sin_theta + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)) - (2.4576 * (r/Rmax) * pow(((r/Rmax) - 1.0), 5.0) * cos(11.0 * theta) + 2.4576 * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta)) * (pow((map1_kappa - 1.0), 2.0) * pow(sin_theta, 2.0) + (map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta) * cos_theta + pow((map1_kappa + 1.0), 2.0) * cos(2.0 * theta)) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin_theta, 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta) * sin_theta + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)) + 6.0 * (2.4576 * (r/Rmax) * pow(((r/Rmax) - 1.0), 5.0) * cos(11.0 * theta) + 2.4576 * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin_theta, 2.0) + 2.0 * map1_kappa + 1.0) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin_theta, 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta) * sin_theta + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)) - (2.4576 * (r/Rmax) * pow(((r/Rmax) - 1.0), 5.0) * cos(11.0 * theta) + 2.4576 * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta)) * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta) * sin_theta + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)) * (4.0 * map1_kappa * (pow((map1_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta), 2.0)) * sin_theta * cos_theta - 1.0 / 2.0 * (pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta) + (2.0 * map1_kappa - 2.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta) * sin_theta) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin_theta, 2.0) + 2.0 * map1_kappa + 1.0) + 1.0 / 2.0 * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta) * sin_theta + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)) * (2.0 * pow((map1_kappa - 1.0), 2.0) * pow(sin_theta, 2.0) + 2.0 * (map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta) * cos_theta + 2.0 * pow((map1_kappa + 1.0), 2.0) * cos(2.0 * theta))) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) / pow(((pow((map1_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin_theta, 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta) * sin_theta + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)), (3.0 / 2.0))) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin_theta, 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta + cos_theta) * sin_theta + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0));
}

double PolarR6GyroZoniShiftedShafranov::rho_pole(double r, double theta, double map1_kappa, double map1_delta, double Rmax, double sin_theta, double cos_theta) const
{
    return 0.0;
}

double PolarR6GyroZoniShiftedShafranov::phi_exact(double r, double theta, double map1_kappa, double map1_delta, double Rmax, double sin_theta, double cos_theta) const
{
    return 0.4096 * pow((r/Rmax), 6.0) * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta);
}
