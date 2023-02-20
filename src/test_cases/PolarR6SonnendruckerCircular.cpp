#include "PolarR6SonnendruckerCircular.h"
#include <stdlib.h>
#include <math.h>


/*........................................*/
double PolarR6SonnendruckerCircular::x(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (r/Rmax) * cos(theta);
}
/*........................................*/
void PolarR6SonnendruckerCircular::x(std::vector<double> const& r, double theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (r[i]/Rmax) * cos(theta);
    }
}
/*........................................*/
void PolarR6SonnendruckerCircular::x(double r, std::vector<double> const& theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (r/Rmax) * cos_theta[i];
    }
}
/*........................................*/
double PolarR6SonnendruckerCircular::y(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (r/Rmax) * sin(theta);
}
/*........................................*/
void PolarR6SonnendruckerCircular::y(std::vector<double> const& r, double theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (r[i]/Rmax) * sin(theta);
    }
}
/*........................................*/
void PolarR6SonnendruckerCircular::y(double r, std::vector<double> const& theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (r/Rmax) * sin_theta[i];
    }
}
/*........................................*/
double PolarR6SonnendruckerCircular::J_rr(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (cos(theta))/Rmax;
}
/*........................................*/
void PolarR6SonnendruckerCircular::J_rr(std::vector<double> const& r, double theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (cos(theta))/Rmax;
    }
}
/*........................................*/
void PolarR6SonnendruckerCircular::J_rr(double r, std::vector<double> const& theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (cos_theta[i])/Rmax;
    }
}
/*........................................*/
double PolarR6SonnendruckerCircular::J_rt(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (-(r/Rmax)) * sin(theta);
}
/*........................................*/
void PolarR6SonnendruckerCircular::J_rt(std::vector<double> const& r, double theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (-(r[i]/Rmax)) * sin(theta);
    }
}
/*........................................*/
void PolarR6SonnendruckerCircular::J_rt(double r, std::vector<double> const& theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (-(r/Rmax)) * sin_theta[i];
    }
}
/*........................................*/
double PolarR6SonnendruckerCircular::J_tr(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (sin(theta))/Rmax;
}
/*........................................*/
void PolarR6SonnendruckerCircular::J_tr(std::vector<double> const& r, double theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (sin(theta))/Rmax;
    }
}
/*........................................*/
void PolarR6SonnendruckerCircular::J_tr(double r, std::vector<double> const& theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (sin_theta[i])/Rmax;
    }
}
/*........................................*/
double PolarR6SonnendruckerCircular::J_tt(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (r/Rmax) * cos(theta);
}
/*........................................*/
void PolarR6SonnendruckerCircular::J_tt(std::vector<double> const& r, double theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (r[i]/Rmax) * cos(theta);
    }
}
/*........................................*/
void PolarR6SonnendruckerCircular::J_tt(double r, std::vector<double> const& theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (r/Rmax) * cos_theta[i];
    }
}
/*........................................*/
double PolarR6SonnendruckerCircular::J_xs(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (-pow(sin(theta), 2.0)) / cos(theta) + pow(cos(theta), (double)((-1)));
}
/*........................................*/
void PolarR6SonnendruckerCircular::J_xs(std::vector<double> const& r, double theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (-pow(sin(theta), 2.0)) / cos(theta) + pow(cos(theta), (double)((-1)));
    }
}
/*........................................*/
void PolarR6SonnendruckerCircular::J_xs(double r, std::vector<double> const& theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (-pow(sin_theta[i], 2.0)) / cos_theta[i] + pow(cos_theta[i], (double)((-1)));
    }
}
/*........................................*/
double PolarR6SonnendruckerCircular::J_xt(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return sin(theta);
}
/*........................................*/
void PolarR6SonnendruckerCircular::J_xt(std::vector<double> const& r, double theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = sin(theta);
    }
}
/*........................................*/
void PolarR6SonnendruckerCircular::J_xt(double r, std::vector<double> const& theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = sin_theta[i];
    }
}
/*........................................*/
double PolarR6SonnendruckerCircular::J_ys(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (-sin(theta)) / (r/Rmax);
}
/*........................................*/
void PolarR6SonnendruckerCircular::J_ys(std::vector<double> const& r, double theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (-sin(theta)) / (r[i]/Rmax);
    }
}
/*........................................*/
void PolarR6SonnendruckerCircular::J_ys(double r, std::vector<double> const& theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (-sin_theta[i]) / (r/Rmax);
    }
}
/*........................................*/
double PolarR6SonnendruckerCircular::J_yt(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return cos(theta) / (r/Rmax);
}
/*........................................*/
void PolarR6SonnendruckerCircular::J_yt(std::vector<double> const& r, double theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = cos(theta) / (r[i]/Rmax);
    }
}
/*........................................*/
void PolarR6SonnendruckerCircular::J_yt(double r, std::vector<double> const& theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = cos_theta[i] / (r/Rmax);
    }
}
/*........................................*/
double PolarR6SonnendruckerCircular::rho_glob(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (-pow((r/Rmax), 4.0)) * ((r/Rmax) * (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r/Rmax) - 11.1111111111111)) * (12.288 * (r/Rmax) * pow(((r/Rmax) - 1.0), 4.0) * cos(11.0 * theta) + 17.2032 * pow(((r/Rmax) - 1.0), 5.0) * cos(11.0 * theta)) - 5.03290747193186 * (r/Rmax) * (2.4576 * (r/Rmax) * pow(((r/Rmax) - 1.0), 5.0) * cos(11.0 * theta) + 2.4576 * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta)) / (208.641975308642 * pow(((r/Rmax) - 0.769230769230769), 2.0) + 1.0) - 49.5616 * (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r/Rmax) - 11.1111111111111)) * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta) + 6.0 * (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r/Rmax) - 11.1111111111111)) * (2.4576 * (r/Rmax) * pow(((r/Rmax) - 1.0), 5.0) * cos(11.0 * theta) + 2.4576 * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta)));
}
/*........................................*/
void PolarR6SonnendruckerCircular::rho_glob(std::vector<double> const& r, double theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (-pow((r[i]/Rmax), 4.0)) * ((r[i]/Rmax) * (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r[i]/Rmax) - 11.1111111111111)) * (12.288 * (r[i]/Rmax) * pow(((r[i]/Rmax) - 1.0), 4.0) * cos(11.0 * theta) + 17.2032 * pow(((r[i]/Rmax) - 1.0), 5.0) * cos(11.0 * theta)) - 5.03290747193186 * (r[i]/Rmax) * (2.4576 * (r[i]/Rmax) * pow(((r[i]/Rmax) - 1.0), 5.0) * cos(11.0 * theta) + 2.4576 * pow(((r[i]/Rmax) - 1.0), 6.0) * cos(11.0 * theta)) / (208.641975308642 * pow(((r[i]/Rmax) - 0.769230769230769), 2.0) + 1.0) - 49.5616 * (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r[i]/Rmax) - 11.1111111111111)) * pow(((r[i]/Rmax) - 1.0), 6.0) * cos(11.0 * theta) + 6.0 * (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r[i]/Rmax) - 11.1111111111111)) * (2.4576 * (r[i]/Rmax) * pow(((r[i]/Rmax) - 1.0), 5.0) * cos(11.0 * theta) + 2.4576 * pow(((r[i]/Rmax) - 1.0), 6.0) * cos(11.0 * theta)));
    }
}
/*........................................*/
void PolarR6SonnendruckerCircular::rho_glob(double r, std::vector<double> const& theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (-pow((r/Rmax), 4.0)) * ((r/Rmax) * (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r/Rmax) - 11.1111111111111)) * (12.288 * (r/Rmax) * pow(((r/Rmax) - 1.0), 4.0) * cos(11.0 * theta[i]) + 17.2032 * pow(((r/Rmax) - 1.0), 5.0) * cos(11.0 * theta[i])) - 5.03290747193186 * (r/Rmax) * (2.4576 * (r/Rmax) * pow(((r/Rmax) - 1.0), 5.0) * cos(11.0 * theta[i]) + 2.4576 * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta[i])) / (208.641975308642 * pow(((r/Rmax) - 0.769230769230769), 2.0) + 1.0) - 49.5616 * (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r/Rmax) - 11.1111111111111)) * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta[i]) + 6.0 * (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r/Rmax) - 11.1111111111111)) * (2.4576 * (r/Rmax) * pow(((r/Rmax) - 1.0), 5.0) * cos(11.0 * theta[i]) + 2.4576 * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta[i])));
    }
}
/*........................................*/
double PolarR6SonnendruckerCircular::rho_pole(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return 0.0;
}
/*........................................*/
void PolarR6SonnendruckerCircular::rho_pole(std::vector<double> const& r, double theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = 0.0;
    }
}
/*........................................*/
void PolarR6SonnendruckerCircular::rho_pole(double r, std::vector<double> const& theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = 0.0;
    }
}
/*........................................*/
double PolarR6SonnendruckerCircular::coeffs1(double r, double Rmax) const
{
    return 0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r/Rmax) - 11.1111111111111);
}
/*........................................*/
void PolarR6SonnendruckerCircular::coeffs1(std::vector<double> const& r, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = 0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r[i]/Rmax) - 11.1111111111111);
    }
}
/*........................................*/
double PolarR6SonnendruckerCircular::coeffs2(double r, double Rmax) const
{
    return 0.0;
}
/*........................................*/
void PolarR6SonnendruckerCircular::coeffs2(std::vector<double> const& r, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = 0.0;
    }
}
/*........................................*/
double PolarR6SonnendruckerCircular::phi_exact(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return 0.4096 * pow((r/Rmax), 6.0) * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta);
}
/*........................................*/
void PolarR6SonnendruckerCircular::phi_exact(std::vector<double> const& r, double theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = 0.4096 * pow((r[i]/Rmax), 6.0) * pow(((r[i]/Rmax) - 1.0), 6.0) * cos(11.0 * theta);
    }
}
/*........................................*/
void PolarR6SonnendruckerCircular::phi_exact(double r, std::vector<double> const& theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = 0.4096 * pow((r/Rmax), 6.0) * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta[i]);
    }
}
/*........................................*/
