// PolarR6 simulates solution (22) of Bourne et al. https://doi.org/10.1016/j.jcp.2023.112249
#include "../../include/test_cases/PolarR6GyroSonnendruckerCircular.h"

double PolarR6GyroSonnendruckerCircular::x(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (r/Rmax) * cos(theta);
}

double PolarR6GyroSonnendruckerCircular::y(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (r/Rmax) * sin(theta);
}

double PolarR6GyroSonnendruckerCircular::J_rr(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (cos(theta))/Rmax;
}

double PolarR6GyroSonnendruckerCircular::J_rt(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (-(r/Rmax)) * sin(theta);
}

double PolarR6GyroSonnendruckerCircular::J_tr(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (sin(theta))/Rmax;
}

double PolarR6GyroSonnendruckerCircular::J_tt(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (r/Rmax) * cos(theta);
}

double PolarR6GyroSonnendruckerCircular::J_xs(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (-pow(sin(theta), 2.0)) / cos(theta) + pow(cos(theta), (double)((-1)));
}

double PolarR6GyroSonnendruckerCircular::J_xt(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return sin(theta);
}

double PolarR6GyroSonnendruckerCircular::J_ys(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (-sin(theta)) / (r/Rmax);
}

double PolarR6GyroSonnendruckerCircular::J_yt(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return cos(theta) / (r/Rmax);
}

double PolarR6GyroSonnendruckerCircular::rho_glob(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return 0.4096 * pow((r/Rmax), 6.0) * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta) / (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r/Rmax) - 11.1111111111111)) - pow((r/Rmax), 4.0) * ((r/Rmax) * (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r/Rmax) - 11.1111111111111)) * (12.288 * (r/Rmax) * pow(((r/Rmax) - 1.0), 4.0) * cos(11.0 * theta) + 17.2032 * pow(((r/Rmax) - 1.0), 5.0) * cos(11.0 * theta)) - 5.03290747193186 * (r/Rmax) * (2.4576 * (r/Rmax) * pow(((r/Rmax) - 1.0), 5.0) * cos(11.0 * theta) + 2.4576 * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta)) / (208.641975308642 * pow(((r/Rmax) - 0.769230769230769), 2.0) + 1.0) - 49.5616 * (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r/Rmax) - 11.1111111111111)) * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta) + 6.0 * (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r/Rmax) - 11.1111111111111)) * (2.4576 * (r/Rmax) * pow(((r/Rmax) - 1.0), 5.0) * cos(11.0 * theta) + 2.4576 * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta)));
}

double PolarR6GyroSonnendruckerCircular::rho_pole(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return 0.0;
}

double PolarR6GyroSonnendruckerCircular::coeffs1(double r, double Rmax) const
{ // see Kuehn et al. https://doi.org/10.1007/s10915-022-01802-1, Eq. (2.3) with Rmax=1.3
    return 0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r/Rmax) - 11.1111111111111);
}

double PolarR6GyroSonnendruckerCircular::coeffs2(double r, double Rmax) const
{
    return pow((0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r/Rmax) - 11.1111111111111)), (double)((-1)));
}

double PolarR6GyroSonnendruckerCircular::phi_exact(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return 0.4096 * pow((r/Rmax), 6.0) * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta);
}



double PolarR6GyroSonnendruckerCircular::x(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return (r/Rmax) * cos_theta;
}

double PolarR6GyroSonnendruckerCircular::y(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return (r/Rmax) * sin_theta;
}

double PolarR6GyroSonnendruckerCircular::J_rr(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return (cos_theta)/Rmax;
}

double PolarR6GyroSonnendruckerCircular::J_rt(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return (-(r/Rmax)) * sin_theta;
}

double PolarR6GyroSonnendruckerCircular::J_tr(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return (sin_theta)/Rmax;
}

double PolarR6GyroSonnendruckerCircular::J_tt(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return (r/Rmax) * cos_theta;
}

double PolarR6GyroSonnendruckerCircular::J_xs(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return (-pow(sin_theta, 2.0)) / cos_theta + pow(cos_theta, (double)((-1)));
}

double PolarR6GyroSonnendruckerCircular::J_xt(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return sin_theta;
}

double PolarR6GyroSonnendruckerCircular::J_ys(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return (-sin_theta) / (r/Rmax);
}

double PolarR6GyroSonnendruckerCircular::J_yt(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return cos_theta / (r/Rmax);
}

double PolarR6GyroSonnendruckerCircular::rho_glob(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return 0.4096 * pow((r/Rmax), 6.0) * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta) / (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r/Rmax) - 11.1111111111111)) - pow((r/Rmax), 4.0) * ((r/Rmax) * (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r/Rmax) - 11.1111111111111)) * (12.288 * (r/Rmax) * pow(((r/Rmax) - 1.0), 4.0) * cos(11.0 * theta) + 17.2032 * pow(((r/Rmax) - 1.0), 5.0) * cos(11.0 * theta)) - 5.03290747193186 * (r/Rmax) * (2.4576 * (r/Rmax) * pow(((r/Rmax) - 1.0), 5.0) * cos(11.0 * theta) + 2.4576 * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta)) / (208.641975308642 * pow(((r/Rmax) - 0.769230769230769), 2.0) + 1.0) - 49.5616 * (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r/Rmax) - 11.1111111111111)) * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta) + 6.0 * (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r/Rmax) - 11.1111111111111)) * (2.4576 * (r/Rmax) * pow(((r/Rmax) - 1.0), 5.0) * cos(11.0 * theta) + 2.4576 * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta)));
}

double PolarR6GyroSonnendruckerCircular::rho_pole(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return 0.0;
}

double PolarR6GyroSonnendruckerCircular::phi_exact(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return 0.4096 * pow((r/Rmax), 6.0) * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta);
}

