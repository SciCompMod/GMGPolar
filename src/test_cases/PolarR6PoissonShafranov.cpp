#include "PolarR6PoissonShafranov.h"
#include <stdlib.h>
#include <math.h>


/*........................................*/
double PolarR6PoissonShafranov::x(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return (-map1_delta) * ((r/Rmax) * (r/Rmax)) - map1_kappa * (r/Rmax) * cos(theta) + (r/Rmax) * cos(theta);
}
/*........................................*/
void PolarR6PoissonShafranov::x(std::vector<double> const& r, double theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (-map1_delta) * ((r[i]/Rmax) * (r[i]/Rmax)) - map1_kappa * (r[i]/Rmax) * cos(theta) + (r[i]/Rmax) * cos(theta);
    }
}
/*........................................*/
void PolarR6PoissonShafranov::x(double r, std::vector<double> const& theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (-map1_delta) * ((r/Rmax) * (r/Rmax)) - map1_kappa * (r/Rmax) * cos_theta[i] + (r/Rmax) * cos_theta[i];
    }
}
/*........................................*/
double PolarR6PoissonShafranov::y(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return map1_kappa * (r/Rmax) * sin(theta) + (r/Rmax) * sin(theta);
}
/*........................................*/
void PolarR6PoissonShafranov::y(std::vector<double> const& r, double theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = map1_kappa * (r[i]/Rmax) * sin(theta) + (r[i]/Rmax) * sin(theta);
    }
}
/*........................................*/
void PolarR6PoissonShafranov::y(double r, std::vector<double> const& theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = map1_kappa * (r/Rmax) * sin_theta[i] + (r/Rmax) * sin_theta[i];
    }
}
/*........................................*/
double PolarR6PoissonShafranov::J_rr(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta))/Rmax;
}
/*........................................*/
void PolarR6PoissonShafranov::J_rr(std::vector<double> const& r, double theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = ((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta))/Rmax;
    }
}
/*........................................*/
void PolarR6PoissonShafranov::J_rr(double r, std::vector<double> const& theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i])/Rmax;
    }
}
/*........................................*/
double PolarR6PoissonShafranov::J_rt(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return (r/Rmax) * (map1_kappa * sin(theta) - sin(theta));
}
/*........................................*/
void PolarR6PoissonShafranov::J_rt(std::vector<double> const& r, double theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (r[i]/Rmax) * (map1_kappa * sin(theta) - sin(theta));
    }
}
/*........................................*/
void PolarR6PoissonShafranov::J_rt(double r, std::vector<double> const& theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (r/Rmax) * (map1_kappa * sin_theta[i] - sin_theta[i]);
    }
}
/*........................................*/
double PolarR6PoissonShafranov::J_tr(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return ((map1_kappa + 1.0) * sin(theta))/Rmax;
}
/*........................................*/
void PolarR6PoissonShafranov::J_tr(std::vector<double> const& r, double theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = ((map1_kappa + 1.0) * sin(theta))/Rmax;
    }
}
/*........................................*/
void PolarR6PoissonShafranov::J_tr(double r, std::vector<double> const& theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = ((map1_kappa + 1.0) * sin_theta[i])/Rmax;
    }
}
/*........................................*/
double PolarR6PoissonShafranov::J_tt(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return (r/Rmax) * (map1_kappa * cos(theta) + cos(theta));
}
/*........................................*/
void PolarR6PoissonShafranov::J_tt(std::vector<double> const& r, double theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (r[i]/Rmax) * (map1_kappa * cos(theta) + cos(theta));
    }
}
/*........................................*/
void PolarR6PoissonShafranov::J_tt(double r, std::vector<double> const& theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (r/Rmax) * (map1_kappa * cos_theta[i] + cos_theta[i]);
    }
}
/*........................................*/
double PolarR6PoissonShafranov::J_xs(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return (-cos(theta)) / (2.0 * map1_delta * (r/Rmax) * cos(theta) + map1_kappa * pow(sin(theta), 2.0) + map1_kappa * pow(cos(theta), 2.0) - pow(sin(theta), 2.0) - pow(cos(theta), 2.0));
}
/*........................................*/
void PolarR6PoissonShafranov::J_xs(std::vector<double> const& r, double theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (-cos(theta)) / (2.0 * map1_delta * (r[i]/Rmax) * cos(theta) + map1_kappa * pow(sin(theta), 2.0) + map1_kappa * pow(cos(theta), 2.0) - pow(sin(theta), 2.0) - pow(cos(theta), 2.0));
    }
}
/*........................................*/
void PolarR6PoissonShafranov::J_xs(double r, std::vector<double> const& theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (-cos_theta[i]) / (2.0 * map1_delta * (r/Rmax) * cos_theta[i] + map1_kappa * pow(sin_theta[i], 2.0) + map1_kappa * pow(cos_theta[i], 2.0) - pow(sin_theta[i], 2.0) - pow(cos_theta[i], 2.0));
    }
}
/*........................................*/
double PolarR6PoissonShafranov::J_xt(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return (map1_kappa * sin(theta) - sin(theta)) / (2.0 * map1_delta * map1_kappa * (r/Rmax) * cos(theta) + 2.0 * map1_delta * (r/Rmax) * cos(theta) + map1_kappa * map1_kappa * pow(sin(theta), 2.0) + map1_kappa * map1_kappa * pow(cos(theta), 2.0) - pow(sin(theta), 2.0) - pow(cos(theta), 2.0));
}
/*........................................*/
void PolarR6PoissonShafranov::J_xt(std::vector<double> const& r, double theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (map1_kappa * sin(theta) - sin(theta)) / (2.0 * map1_delta * map1_kappa * (r[i]/Rmax) * cos(theta) + 2.0 * map1_delta * (r[i]/Rmax) * cos(theta) + map1_kappa * map1_kappa * pow(sin(theta), 2.0) + map1_kappa * map1_kappa * pow(cos(theta), 2.0) - pow(sin(theta), 2.0) - pow(cos(theta), 2.0));
    }
}
/*........................................*/
void PolarR6PoissonShafranov::J_xt(double r, std::vector<double> const& theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (map1_kappa * sin_theta[i] - sin_theta[i]) / (2.0 * map1_delta * map1_kappa * (r/Rmax) * cos_theta[i] + 2.0 * map1_delta * (r/Rmax) * cos_theta[i] + map1_kappa * map1_kappa * pow(sin_theta[i], 2.0) + map1_kappa * map1_kappa * pow(cos_theta[i], 2.0) - pow(sin_theta[i], 2.0) - pow(cos_theta[i], 2.0));
    }
}
/*........................................*/
double PolarR6PoissonShafranov::J_ys(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return sin(theta) / (2.0 * map1_delta * ((r/Rmax) * (r/Rmax)) * cos(theta) + map1_kappa * (r/Rmax) * pow(sin(theta), 2.0) + map1_kappa * (r/Rmax) * pow(cos(theta), 2.0) - (r/Rmax) * pow(sin(theta), 2.0) - (r/Rmax) * pow(cos(theta), 2.0));
}
/*........................................*/
void PolarR6PoissonShafranov::J_ys(std::vector<double> const& r, double theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = sin(theta) / (2.0 * map1_delta * ((r[i]/Rmax) * (r[i]/Rmax)) * cos(theta) + map1_kappa * (r[i]/Rmax) * pow(sin(theta), 2.0) + map1_kappa * (r[i]/Rmax) * pow(cos(theta), 2.0) - (r[i]/Rmax) * pow(sin(theta), 2.0) - (r[i]/Rmax) * pow(cos(theta), 2.0));
    }
}
/*........................................*/
void PolarR6PoissonShafranov::J_ys(double r, std::vector<double> const& theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = sin_theta[i] / (2.0 * map1_delta * ((r/Rmax) * (r/Rmax)) * cos_theta[i] + map1_kappa * (r/Rmax) * pow(sin_theta[i], 2.0) + map1_kappa * (r/Rmax) * pow(cos_theta[i], 2.0) - (r/Rmax) * pow(sin_theta[i], 2.0) - (r/Rmax) * pow(cos_theta[i], 2.0));
    }
}
/*........................................*/
double PolarR6PoissonShafranov::J_yt(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return (2.0 * map1_delta * (r/Rmax) + map1_kappa * cos(theta) - cos(theta)) / (2.0 * map1_delta * map1_kappa * ((r/Rmax) * (r/Rmax)) * cos(theta) + 2.0 * map1_delta * ((r/Rmax) * (r/Rmax)) * cos(theta) + map1_kappa * map1_kappa * (r/Rmax) * pow(sin(theta), 2.0) + map1_kappa * map1_kappa * (r/Rmax) * pow(cos(theta), 2.0) - (r/Rmax) * pow(sin(theta), 2.0) - (r/Rmax) * pow(cos(theta), 2.0));
}
/*........................................*/
void PolarR6PoissonShafranov::J_yt(std::vector<double> const& r, double theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (2.0 * map1_delta * (r[i]/Rmax) + map1_kappa * cos(theta) - cos(theta)) / (2.0 * map1_delta * map1_kappa * ((r[i]/Rmax) * (r[i]/Rmax)) * cos(theta) + 2.0 * map1_delta * ((r[i]/Rmax) * (r[i]/Rmax)) * cos(theta) + map1_kappa * map1_kappa * (r[i]/Rmax) * pow(sin(theta), 2.0) + map1_kappa * map1_kappa * (r[i]/Rmax) * pow(cos(theta), 2.0) - (r[i]/Rmax) * pow(sin(theta), 2.0) - (r[i]/Rmax) * pow(cos(theta), 2.0));
    }
}
/*........................................*/
void PolarR6PoissonShafranov::J_yt(double r, std::vector<double> const& theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (2.0 * map1_delta * (r/Rmax) + map1_kappa * cos_theta[i] - cos_theta[i]) / (2.0 * map1_delta * map1_kappa * ((r/Rmax) * (r/Rmax)) * cos_theta[i] + 2.0 * map1_delta * ((r/Rmax) * (r/Rmax)) * cos_theta[i] + map1_kappa * map1_kappa * (r/Rmax) * pow(sin_theta[i], 2.0) + map1_kappa * map1_kappa * (r/Rmax) * pow(cos_theta[i], 2.0) - (r/Rmax) * pow(sin_theta[i], 2.0) - (r/Rmax) * pow(cos_theta[i], 2.0));
    }
}
/*........................................*/
double PolarR6PoissonShafranov::rho_glob(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return (-pow((r/Rmax), 4.0)) * ((-9.0112) * map1_delta * (r/Rmax) * (map1_kappa - 1.0) * pow(((r/Rmax) - 1.0), 6.0) * sin(theta) * sin(11.0 * theta) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)) + 4.5056 * (r/Rmax) * pow(((r/Rmax) - 1.0), 6.0) * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)) * ((-2.0) * map1_delta * (map1_kappa - 1.0) * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)) * sin(theta) + 2.0 * map1_delta * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0)) * sin(11.0 * theta) / pow(((pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)), (3.0 / 2.0)) + 27.0336 * (r/Rmax) * pow(((r/Rmax) - 1.0), 5.0) * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)) * sin(11.0 * theta) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)) + 1.0 * (r/Rmax) * (12.288 * (r/Rmax) * pow(((r/Rmax) - 1.0), 4.0) * cos(11.0 * theta) + 17.2032 * pow(((r/Rmax) - 1.0), 5.0) * cos(11.0 * theta)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)) + 1.0 * (r/Rmax) * (2.4576 * (r/Rmax) * pow(((r/Rmax) - 1.0), 5.0) * cos(11.0 * theta) + 2.4576 * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta)) * ((-2.0) * map1_delta * (map1_kappa - 1.0) * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)) * sin(theta) + 2.0 * map1_delta * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) / pow(((pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)), (3.0 / 2.0)) - 49.5616 * pow(((r/Rmax) - 1.0), 6.0) * (pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * cos(11.0 * theta) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)) - 4.5056 * pow(((r/Rmax) - 1.0), 6.0) * (pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (4.0 * map1_kappa * (pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * sin(theta) * cos(theta) - 1.0 / 2.0 * (pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta) + (2.0 * map1_kappa - 2.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) + 1.0 / 2.0 * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)) * (2.0 * pow((map1_kappa - 1.0), 2.0) * pow(sin(theta), 2.0) + 2.0 * (map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * cos(theta) + 2.0 * pow((map1_kappa + 1.0), 2.0) * cos(2.0 * theta))) * sin(11.0 * theta) / pow(((pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)), (3.0 / 2.0)) - 4.5056 * pow(((r/Rmax) - 1.0), 6.0) * (pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta) + (2.0 * map1_kappa - 2.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta)) * sin(11.0 * theta) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)) + 27.0336 * pow(((r/Rmax) - 1.0), 6.0) * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)) * sin(11.0 * theta) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)) - 1.0 * ((-27.0336) * (r/Rmax) * pow(((r/Rmax) - 1.0), 5.0) * sin(11.0 * theta) - 27.0336 * pow(((r/Rmax) - 1.0), 6.0) * sin(11.0 * theta)) * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)) - 1.0 * (2.4576 * (r/Rmax) * pow(((r/Rmax) - 1.0), 5.0) * cos(11.0 * theta) + 2.4576 * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta)) * (pow((map1_kappa - 1.0), 2.0) * pow(sin(theta), 2.0) + (map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * cos(theta) + pow((map1_kappa + 1.0), 2.0) * cos(2.0 * theta)) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)) + 6.0 * (2.4576 * (r/Rmax) * pow(((r/Rmax) - 1.0), 5.0) * cos(11.0 * theta) + 2.4576 * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)) - 1.0 * (2.4576 * (r/Rmax) * pow(((r/Rmax) - 1.0), 5.0) * cos(11.0 * theta) + 2.4576 * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta)) * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)) * (4.0 * map1_kappa * (pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * sin(theta) * cos(theta) - 1.0 / 2.0 * (pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta) + (2.0 * map1_kappa - 2.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) + 1.0 / 2.0 * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)) * (2.0 * pow((map1_kappa - 1.0), 2.0) * pow(sin(theta), 2.0) + 2.0 * (map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * cos(theta) + 2.0 * pow((map1_kappa + 1.0), 2.0) * cos(2.0 * theta))) / pow(((pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)), (3.0 / 2.0))) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0));
}
/*........................................*/
void PolarR6PoissonShafranov::rho_glob(std::vector<double> const& r, double theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (-pow((r[i]/Rmax), 4.0)) * ((-9.0112) * map1_delta * (r[i]/Rmax) * (map1_kappa - 1.0) * pow(((r[i]/Rmax) - 1.0), 6.0) * sin(theta) * sin(11.0 * theta) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)) + 4.5056 * (r[i]/Rmax) * pow(((r[i]/Rmax) - 1.0), 6.0) * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)) * ((-2.0) * map1_delta * (map1_kappa - 1.0) * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)) * sin(theta) + 2.0 * map1_delta * ((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0)) * sin(11.0 * theta) / pow(((pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)), (3.0 / 2.0)) + 27.0336 * (r[i]/Rmax) * pow(((r[i]/Rmax) - 1.0), 5.0) * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)) * sin(11.0 * theta) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)) + 1.0 * (r[i]/Rmax) * (12.288 * (r[i]/Rmax) * pow(((r[i]/Rmax) - 1.0), 4.0) * cos(11.0 * theta) + 17.2032 * pow(((r[i]/Rmax) - 1.0), 5.0) * cos(11.0 * theta)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)) + 1.0 * (r[i]/Rmax) * (2.4576 * (r[i]/Rmax) * pow(((r[i]/Rmax) - 1.0), 5.0) * cos(11.0 * theta) + 2.4576 * pow(((r[i]/Rmax) - 1.0), 6.0) * cos(11.0 * theta)) * ((-2.0) * map1_delta * (map1_kappa - 1.0) * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)) * sin(theta) + 2.0 * map1_delta * ((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) / pow(((pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)), (3.0 / 2.0)) - 49.5616 * pow(((r[i]/Rmax) - 1.0), 6.0) * (pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * cos(11.0 * theta) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)) - 4.5056 * pow(((r[i]/Rmax) - 1.0), 6.0) * (pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (4.0 * map1_kappa * (pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * sin(theta) * cos(theta) - 1.0 / 2.0 * (pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta) + (2.0 * map1_kappa - 2.0) * ((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) + 1.0 / 2.0 * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)) * (2.0 * pow((map1_kappa - 1.0), 2.0) * pow(sin(theta), 2.0) + 2.0 * (map1_kappa - 1.0) * ((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)) * cos(theta) + 2.0 * pow((map1_kappa + 1.0), 2.0) * cos(2.0 * theta))) * sin(11.0 * theta) / pow(((pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)), (3.0 / 2.0)) - 4.5056 * pow(((r[i]/Rmax) - 1.0), 6.0) * (pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta) + (2.0 * map1_kappa - 2.0) * ((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta)) * sin(11.0 * theta) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)) + 27.0336 * pow(((r[i]/Rmax) - 1.0), 6.0) * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)) * sin(11.0 * theta) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)) - 1.0 * ((-27.0336) * (r[i]/Rmax) * pow(((r[i]/Rmax) - 1.0), 5.0) * sin(11.0 * theta) - 27.0336 * pow(((r[i]/Rmax) - 1.0), 6.0) * sin(11.0 * theta)) * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)) - 1.0 * (2.4576 * (r[i]/Rmax) * pow(((r[i]/Rmax) - 1.0), 5.0) * cos(11.0 * theta) + 2.4576 * pow(((r[i]/Rmax) - 1.0), 6.0) * cos(11.0 * theta)) * (pow((map1_kappa - 1.0), 2.0) * pow(sin(theta), 2.0) + (map1_kappa - 1.0) * ((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)) * cos(theta) + pow((map1_kappa + 1.0), 2.0) * cos(2.0 * theta)) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)) + 6.0 * (2.4576 * (r[i]/Rmax) * pow(((r[i]/Rmax) - 1.0), 5.0) * cos(11.0 * theta) + 2.4576 * pow(((r[i]/Rmax) - 1.0), 6.0) * cos(11.0 * theta)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)) - 1.0 * (2.4576 * (r[i]/Rmax) * pow(((r[i]/Rmax) - 1.0), 5.0) * cos(11.0 * theta) + 2.4576 * pow(((r[i]/Rmax) - 1.0), 6.0) * cos(11.0 * theta)) * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)) * (4.0 * map1_kappa * (pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * sin(theta) * cos(theta) - 1.0 / 2.0 * (pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta) + (2.0 * map1_kappa - 2.0) * ((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) + 1.0 / 2.0 * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)) * (2.0 * pow((map1_kappa - 1.0), 2.0) * pow(sin(theta), 2.0) + 2.0 * (map1_kappa - 1.0) * ((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)) * cos(theta) + 2.0 * pow((map1_kappa + 1.0), 2.0) * cos(2.0 * theta))) / pow(((pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0)), (3.0 / 2.0))) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin(theta), 2.0) + pow(((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin(theta), 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r[i]/Rmax) - map1_kappa * cos(theta) + cos(theta)) * sin(theta) + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta)), 2.0));
    }
}
/*........................................*/
void PolarR6PoissonShafranov::rho_glob(double r, std::vector<double> const& theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (-pow((r/Rmax), 4.0)) * ((-9.0112) * map1_delta * (r/Rmax) * (map1_kappa - 1.0) * pow(((r/Rmax) - 1.0), 6.0) * sin_theta[i] * sin(11.0 * theta[i]) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin_theta[i], 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin_theta[i], 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]) * sin_theta[i] + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta[i])), 2.0)) + 4.5056 * (r/Rmax) * pow(((r/Rmax) - 1.0), 6.0) * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]) * sin_theta[i] + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta[i])) * ((-2.0) * map1_delta * (map1_kappa - 1.0) * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]) * sin_theta[i] + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta[i])) * sin_theta[i] + 2.0 * map1_delta * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin_theta[i], 2.0) + 2.0 * map1_kappa + 1.0)) * sin(11.0 * theta[i]) / pow(((pow((map1_kappa + 1.0), 2.0) * pow(sin_theta[i], 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin_theta[i], 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]) * sin_theta[i] + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta[i])), 2.0)), (3.0 / 2.0)) + 27.0336 * (r/Rmax) * pow(((r/Rmax) - 1.0), 5.0) * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]) * sin_theta[i] + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta[i])) * sin(11.0 * theta[i]) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin_theta[i], 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin_theta[i], 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]) * sin_theta[i] + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta[i])), 2.0)) + 1.0 * (r/Rmax) * (12.288 * (r/Rmax) * pow(((r/Rmax) - 1.0), 4.0) * cos(11.0 * theta[i]) + 17.2032 * pow(((r/Rmax) - 1.0), 5.0) * cos(11.0 * theta[i])) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin_theta[i], 2.0) + 2.0 * map1_kappa + 1.0) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin_theta[i], 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin_theta[i], 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]) * sin_theta[i] + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta[i])), 2.0)) + 1.0 * (r/Rmax) * (2.4576 * (r/Rmax) * pow(((r/Rmax) - 1.0), 5.0) * cos(11.0 * theta[i]) + 2.4576 * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta[i])) * ((-2.0) * map1_delta * (map1_kappa - 1.0) * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]) * sin_theta[i] + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta[i])) * sin_theta[i] + 2.0 * map1_delta * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin_theta[i], 2.0) + 2.0 * map1_kappa + 1.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin_theta[i], 2.0) + 2.0 * map1_kappa + 1.0) / pow(((pow((map1_kappa + 1.0), 2.0) * pow(sin_theta[i], 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin_theta[i], 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]) * sin_theta[i] + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta[i])), 2.0)), (3.0 / 2.0)) - 49.5616 * pow(((r/Rmax) - 1.0), 6.0) * (pow((map1_kappa + 1.0), 2.0) * pow(sin_theta[i], 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]), 2.0)) * cos(11.0 * theta[i]) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin_theta[i], 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin_theta[i], 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]) * sin_theta[i] + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta[i])), 2.0)) - 4.5056 * pow(((r/Rmax) - 1.0), 6.0) * (pow((map1_kappa + 1.0), 2.0) * pow(sin_theta[i], 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]), 2.0)) * (4.0 * map1_kappa * (pow((map1_kappa + 1.0), 2.0) * pow(sin_theta[i], 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]), 2.0)) * sin_theta[i] * cos_theta[i] - 1.0 / 2.0 * (pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta[i]) + (2.0 * map1_kappa - 2.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]) * sin_theta[i]) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin_theta[i], 2.0) + 2.0 * map1_kappa + 1.0) + 1.0 / 2.0 * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]) * sin_theta[i] + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta[i])) * (2.0 * pow((map1_kappa - 1.0), 2.0) * pow(sin_theta[i], 2.0) + 2.0 * (map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]) * cos_theta[i] + 2.0 * pow((map1_kappa + 1.0), 2.0) * cos(2.0 * theta[i]))) * sin(11.0 * theta[i]) / pow(((pow((map1_kappa + 1.0), 2.0) * pow(sin_theta[i], 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin_theta[i], 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]) * sin_theta[i] + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta[i])), 2.0)), (3.0 / 2.0)) - 4.5056 * pow(((r/Rmax) - 1.0), 6.0) * (pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta[i]) + (2.0 * map1_kappa - 2.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]) * sin_theta[i]) * sin(11.0 * theta[i]) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin_theta[i], 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin_theta[i], 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]) * sin_theta[i] + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta[i])), 2.0)) + 27.0336 * pow(((r/Rmax) - 1.0), 6.0) * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]) * sin_theta[i] + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta[i])) * sin(11.0 * theta[i]) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin_theta[i], 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin_theta[i], 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]) * sin_theta[i] + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta[i])), 2.0)) - 1.0 * ((-27.0336) * (r/Rmax) * pow(((r/Rmax) - 1.0), 5.0) * sin(11.0 * theta[i]) - 27.0336 * pow(((r/Rmax) - 1.0), 6.0) * sin(11.0 * theta[i])) * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]) * sin_theta[i] + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta[i])) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin_theta[i], 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin_theta[i], 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]) * sin_theta[i] + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta[i])), 2.0)) - 1.0 * (2.4576 * (r/Rmax) * pow(((r/Rmax) - 1.0), 5.0) * cos(11.0 * theta[i]) + 2.4576 * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta[i])) * (pow((map1_kappa - 1.0), 2.0) * pow(sin_theta[i], 2.0) + (map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]) * cos_theta[i] + pow((map1_kappa + 1.0), 2.0) * cos(2.0 * theta[i])) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin_theta[i], 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin_theta[i], 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]) * sin_theta[i] + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta[i])), 2.0)) + 6.0 * (2.4576 * (r/Rmax) * pow(((r/Rmax) - 1.0), 5.0) * cos(11.0 * theta[i]) + 2.4576 * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta[i])) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin_theta[i], 2.0) + 2.0 * map1_kappa + 1.0) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin_theta[i], 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin_theta[i], 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]) * sin_theta[i] + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta[i])), 2.0)) - 1.0 * (2.4576 * (r/Rmax) * pow(((r/Rmax) - 1.0), 5.0) * cos(11.0 * theta[i]) + 2.4576 * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta[i])) * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]) * sin_theta[i] + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta[i])) * (4.0 * map1_kappa * (pow((map1_kappa + 1.0), 2.0) * pow(sin_theta[i], 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]), 2.0)) * sin_theta[i] * cos_theta[i] - 1.0 / 2.0 * (pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta[i]) + (2.0 * map1_kappa - 2.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]) * sin_theta[i]) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin_theta[i], 2.0) + 2.0 * map1_kappa + 1.0) + 1.0 / 2.0 * ((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]) * sin_theta[i] + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta[i])) * (2.0 * pow((map1_kappa - 1.0), 2.0) * pow(sin_theta[i], 2.0) + 2.0 * (map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]) * cos_theta[i] + 2.0 * pow((map1_kappa + 1.0), 2.0) * cos(2.0 * theta[i]))) / pow(((pow((map1_kappa + 1.0), 2.0) * pow(sin_theta[i], 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin_theta[i], 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]) * sin_theta[i] + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta[i])), 2.0)), (3.0 / 2.0))) / sqrt((pow((map1_kappa + 1.0), 2.0) * pow(sin_theta[i], 2.0) + pow(((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]), 2.0)) * (map1_kappa * map1_kappa - 4.0 * map1_kappa * pow(sin_theta[i], 2.0) + 2.0 * map1_kappa + 1.0) - pow(((map1_kappa - 1.0) * ((-2.0) * map1_delta * (r/Rmax) - map1_kappa * cos_theta[i] + cos_theta[i]) * sin_theta[i] + 1.0 / 2.0 * pow((map1_kappa + 1.0), 2.0) * sin(2.0 * theta[i])), 2.0));
    }
}
/*........................................*/
double PolarR6PoissonShafranov::rho_pole(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return 0.0;
}
/*........................................*/
void PolarR6PoissonShafranov::rho_pole(std::vector<double> const& r, double theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = 0.0;
    }
}
/*........................................*/
void PolarR6PoissonShafranov::rho_pole(double r, std::vector<double> const& theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = 0.0;
    }
}
/*........................................*/
double PolarR6PoissonShafranov::coeffs1(double r, double Rmax) const
{
    return 1.0;
}
/*........................................*/
void PolarR6PoissonShafranov::coeffs1(std::vector<double> const& r, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = 1.0;
    }
}
/*........................................*/
double PolarR6PoissonShafranov::coeffs2(double r, double Rmax) const
{
    return 0.0;
}
/*........................................*/
void PolarR6PoissonShafranov::coeffs2(std::vector<double> const& r, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = 0.0;
    }
}
/*........................................*/
double PolarR6PoissonShafranov::phi_exact(double r, double theta, double map1_kappa, double map1_delta, double Rmax) const
{
    return 0.4096 * pow((r/Rmax), 6.0) * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta);
}
/*........................................*/
void PolarR6PoissonShafranov::phi_exact(std::vector<double> const& r, double theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = 0.4096 * pow((r[i]/Rmax), 6.0) * pow(((r[i]/Rmax) - 1.0), 6.0) * cos(11.0 * theta);
    }
}
/*........................................*/
void PolarR6PoissonShafranov::phi_exact(double r, std::vector<double> const& theta, double map1_kappa, double map1_delta, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = 0.4096 * pow((r/Rmax), 6.0) * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta[i]);
    }
}
/*........................................*/
