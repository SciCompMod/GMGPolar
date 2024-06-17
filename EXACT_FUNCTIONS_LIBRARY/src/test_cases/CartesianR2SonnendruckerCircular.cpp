#include "../../include/test_cases/CartesianR2SonnendruckerCircular.h"

double CartesianR2SonnendruckerCircular::x(double r, double theta, double unused_1, double unused_2, double Rmax) const{
    return (r/Rmax) * cos(theta);
}
double CartesianR2SonnendruckerCircular::x(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const{
    return (r/Rmax) * cos_theta;
}

double CartesianR2SonnendruckerCircular::y(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (r/Rmax) * sin(theta);
}
double CartesianR2SonnendruckerCircular::y(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return (r/Rmax) * sin_theta;
}

double CartesianR2SonnendruckerCircular::J_rr(double r, double theta, double unused_1, double unused_2, double Rmax) const{
    return (cos(theta))/Rmax;
}
double CartesianR2SonnendruckerCircular::J_rr(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const{
    return (cos_theta)/Rmax;
}


double CartesianR2SonnendruckerCircular::J_rt(double r, double theta, double unused_1, double unused_2, double Rmax) const{
    return (-(r/Rmax)) * sin(theta);
}
double CartesianR2SonnendruckerCircular::J_rt(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const{
    return (-(r/Rmax)) * sin_theta;
}


double CartesianR2SonnendruckerCircular::J_tr(double r, double theta, double unused_1, double unused_2, double Rmax) const{
    return (sin(theta))/Rmax;
}
double CartesianR2SonnendruckerCircular::J_tr(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const{
    return (sin_theta)/Rmax;
}


double CartesianR2SonnendruckerCircular::J_tt(double r, double theta, double unused_1, double unused_2, double Rmax) const{
    return (r/Rmax) * cos(theta);
}
double CartesianR2SonnendruckerCircular::J_tt(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const{
    return (r/Rmax) * cos_theta;
}


double CartesianR2SonnendruckerCircular::coeffs1(double r, double Rmax) const
{ // see Kuehn et al. https://doi.org/10.1007/s10915-022-01802-1, Eq. (2.3) with Rmax=1.3
    return 0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r/Rmax) - 11.1111111111111);
}
double CartesianR2SonnendruckerCircular::coeffs2(double r, double Rmax) const
{
    return 0.0;
}


double CartesianR2SonnendruckerCircular::J_xs(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (-pow(sin(theta), 2.0)) / cos(theta) + pow(cos(theta), (double)((-1)));
}
double CartesianR2SonnendruckerCircular::J_xs(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return (-pow(sin_theta, 2.0)) / cos_theta + pow(cos_theta, (double)((-1)));
}


double CartesianR2SonnendruckerCircular::J_xt(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return sin(theta);
}
double CartesianR2SonnendruckerCircular::J_xt(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return sin_theta;
}


double CartesianR2SonnendruckerCircular::J_ys(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (-sin(theta)) / (r/Rmax);
}
double CartesianR2SonnendruckerCircular::J_ys(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return (-sin_theta) / (r/Rmax);
}


double CartesianR2SonnendruckerCircular::J_yt(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return cos(theta) / (r/Rmax);
}
double CartesianR2SonnendruckerCircular::J_yt(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return cos_theta / (r/Rmax);
}


double CartesianR2SonnendruckerCircular::rho_glob(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (-((r/Rmax) * (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r/Rmax) - 11.1111111111111)) * ((-8.0) * M_PI * (r/Rmax) * sin(theta) * cos(2.0 * M_PI * (r/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r/Rmax) * cos(theta)) + 8.0 * M_PI * (r/Rmax) * sin(2.0 * M_PI * (r/Rmax) * sin(theta)) * sin(2.0 * M_PI * (r/Rmax) * cos(theta)) * cos(theta) - 4.0 * (M_PI * M_PI) * (1.0 - (r/Rmax) * (r/Rmax)) * pow(sin(theta), 2.0) * sin(2.0 * M_PI * (r/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r/Rmax) * cos(theta)) - 8.0 * (M_PI * M_PI) * (1.0 - (r/Rmax) * (r/Rmax)) * sin(theta) * sin(2.0 * M_PI * (r/Rmax) * cos(theta)) * cos(theta) * cos(2.0 * M_PI * (r/Rmax) * sin(theta)) - 4.0 * (M_PI * M_PI) * (1.0 - (r/Rmax) * (r/Rmax)) * sin(2.0 * M_PI * (r/Rmax) * sin(theta)) * pow(cos(theta), 2.0) * cos(2.0 * M_PI * (r/Rmax) * cos(theta)) - 2.0 * sin(2.0 * M_PI * (r/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r/Rmax) * cos(theta))) - 5.03290747193186 * (r/Rmax) * ((-2.0) * (r/Rmax) * sin(2.0 * M_PI * (r/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r/Rmax) * cos(theta)) + 2.0 * M_PI * (1.0 - (r/Rmax) * (r/Rmax)) * sin(theta) * cos(2.0 * M_PI * (r/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r/Rmax) * cos(theta)) - 2.0 * M_PI * (1.0 - (r/Rmax) * (r/Rmax)) * sin(2.0 * M_PI * (r/Rmax) * sin(theta)) * sin(2.0 * M_PI * (r/Rmax) * cos(theta)) * cos(theta)) / (208.641975308642 * pow(((r/Rmax) - 0.769230769230769), 2.0) + 1.0) + (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r/Rmax) - 11.1111111111111)) * ((-2.0) * (r/Rmax) * sin(2.0 * M_PI * (r/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r/Rmax) * cos(theta)) + 2.0 * M_PI * (1.0 - (r/Rmax) * (r/Rmax)) * sin(theta) * cos(2.0 * M_PI * (r/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r/Rmax) * cos(theta)) - 2.0 * M_PI * (1.0 - (r/Rmax) * (r/Rmax)) * sin(2.0 * M_PI * (r/Rmax) * sin(theta)) * sin(2.0 * M_PI * (r/Rmax) * cos(theta)) * cos(theta)) + (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r/Rmax) - 11.1111111111111)) * ((-4.0) * (M_PI * M_PI) * (r/Rmax) * (1.0 - (r/Rmax) * (r/Rmax)) * pow(sin(theta), 2.0) * sin(2.0 * M_PI * (r/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r/Rmax) * cos(theta)) + 8.0 * (M_PI * M_PI) * (r/Rmax) * (1.0 - (r/Rmax) * (r/Rmax)) * sin(theta) * sin(2.0 * M_PI * (r/Rmax) * cos(theta)) * cos(theta) * cos(2.0 * M_PI * (r/Rmax) * sin(theta)) - 4.0 * (M_PI * M_PI) * (r/Rmax) * (1.0 - (r/Rmax) * (r/Rmax)) * sin(2.0 * M_PI * (r/Rmax) * sin(theta)) * pow(cos(theta), 2.0) * cos(2.0 * M_PI * (r/Rmax) * cos(theta)) - 2.0 * M_PI * (1.0 - (r/Rmax) * (r/Rmax)) * sin(theta) * cos(2.0 * M_PI * (r/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r/Rmax) * cos(theta)) + 2.0 * M_PI * (1.0 - (r/Rmax) * (r/Rmax)) * sin(2.0 * M_PI * (r/Rmax) * sin(theta)) * sin(2.0 * M_PI * (r/Rmax) * cos(theta)) * cos(theta)))) / (r/Rmax);
}
double CartesianR2SonnendruckerCircular::rho_glob(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return (-((r/Rmax) * (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r/Rmax) - 11.1111111111111)) * ((-8.0) * M_PI * (r/Rmax) * sin_theta * cos(2.0 * M_PI * (r/Rmax) * sin_theta) * cos(2.0 * M_PI * (r/Rmax) * cos_theta) + 8.0 * M_PI * (r/Rmax) * sin(2.0 * M_PI * (r/Rmax) * sin_theta) * sin(2.0 * M_PI * (r/Rmax) * cos_theta) * cos_theta - 4.0 * (M_PI * M_PI) * (1.0 - (r/Rmax) * (r/Rmax)) * pow(sin_theta, 2.0) * sin(2.0 * M_PI * (r/Rmax) * sin_theta) * cos(2.0 * M_PI * (r/Rmax) * cos_theta) - 8.0 * (M_PI * M_PI) * (1.0 - (r/Rmax) * (r/Rmax)) * sin_theta * sin(2.0 * M_PI * (r/Rmax) * cos_theta) * cos_theta * cos(2.0 * M_PI * (r/Rmax) * sin_theta) - 4.0 * (M_PI * M_PI) * (1.0 - (r/Rmax) * (r/Rmax)) * sin(2.0 * M_PI * (r/Rmax) * sin_theta) * pow(cos_theta, 2.0) * cos(2.0 * M_PI * (r/Rmax) * cos_theta) - 2.0 * sin(2.0 * M_PI * (r/Rmax) * sin_theta) * cos(2.0 * M_PI * (r/Rmax) * cos_theta)) - 5.03290747193186 * (r/Rmax) * ((-2.0) * (r/Rmax) * sin(2.0 * M_PI * (r/Rmax) * sin_theta) * cos(2.0 * M_PI * (r/Rmax) * cos_theta) + 2.0 * M_PI * (1.0 - (r/Rmax) * (r/Rmax)) * sin_theta * cos(2.0 * M_PI * (r/Rmax) * sin_theta) * cos(2.0 * M_PI * (r/Rmax) * cos_theta) - 2.0 * M_PI * (1.0 - (r/Rmax) * (r/Rmax)) * sin(2.0 * M_PI * (r/Rmax) * sin_theta) * sin(2.0 * M_PI * (r/Rmax) * cos_theta) * cos_theta) / (208.641975308642 * pow(((r/Rmax) - 0.769230769230769), 2.0) + 1.0) + (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r/Rmax) - 11.1111111111111)) * ((-2.0) * (r/Rmax) * sin(2.0 * M_PI * (r/Rmax) * sin_theta) * cos(2.0 * M_PI * (r/Rmax) * cos_theta) + 2.0 * M_PI * (1.0 - (r/Rmax) * (r/Rmax)) * sin_theta * cos(2.0 * M_PI * (r/Rmax) * sin_theta) * cos(2.0 * M_PI * (r/Rmax) * cos_theta) - 2.0 * M_PI * (1.0 - (r/Rmax) * (r/Rmax)) * sin(2.0 * M_PI * (r/Rmax) * sin_theta) * sin(2.0 * M_PI * (r/Rmax) * cos_theta) * cos_theta) + (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r/Rmax) - 11.1111111111111)) * ((-4.0) * (M_PI * M_PI) * (r/Rmax) * (1.0 - (r/Rmax) * (r/Rmax)) * pow(sin_theta, 2.0) * sin(2.0 * M_PI * (r/Rmax) * sin_theta) * cos(2.0 * M_PI * (r/Rmax) * cos_theta) + 8.0 * (M_PI * M_PI) * (r/Rmax) * (1.0 - (r/Rmax) * (r/Rmax)) * sin_theta * sin(2.0 * M_PI * (r/Rmax) * cos_theta) * cos_theta * cos(2.0 * M_PI * (r/Rmax) * sin_theta) - 4.0 * (M_PI * M_PI) * (r/Rmax) * (1.0 - (r/Rmax) * (r/Rmax)) * sin(2.0 * M_PI * (r/Rmax) * sin_theta) * pow(cos_theta, 2.0) * cos(2.0 * M_PI * (r/Rmax) * cos_theta) - 2.0 * M_PI * (1.0 - (r/Rmax) * (r/Rmax)) * sin_theta * cos(2.0 * M_PI * (r/Rmax) * sin_theta) * cos(2.0 * M_PI * (r/Rmax) * cos_theta) + 2.0 * M_PI * (1.0 - (r/Rmax) * (r/Rmax)) * sin(2.0 * M_PI * (r/Rmax) * sin_theta) * sin(2.0 * M_PI * (r/Rmax) * cos_theta) * cos_theta))) / (r/Rmax);
}


double CartesianR2SonnendruckerCircular::rho_pole(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return 0.0;
}
double CartesianR2SonnendruckerCircular::rho_pole(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return 0.0;
}


double CartesianR2SonnendruckerCircular::phi_exact(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (1.0 - (r/Rmax) * (r/Rmax)) * sin(2.0 * M_PI * (r/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r/Rmax) * cos(theta));
}
double CartesianR2SonnendruckerCircular::phi_exact(double r, double theta, double unused_1, double unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return (1.0 - (r/Rmax) * (r/Rmax)) * sin(2.0 * M_PI * (r/Rmax) * sin_theta) * cos(2.0 * M_PI * (r/Rmax) * cos_theta);
}
