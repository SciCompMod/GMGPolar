// CartesianR6 simulates solution (23) of Bourne et al. https://doi.org/10.1016/j.jcp.2023.112249
#include "CartesianR6ZoniShiftedCircular.h"
#include <stdlib.h>
#include <math.h>


/*........................................*/
double CartesianR6ZoniShiftedCircular::x(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (r/Rmax) * cos(theta);
}
/*........................................*/
void CartesianR6ZoniShiftedCircular::x(std::vector<double> const& r, double theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (r[i]/Rmax) * cos(theta);
    }
}
/*........................................*/
void CartesianR6ZoniShiftedCircular::x(double r, std::vector<double> const& theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (r/Rmax) * cos_theta[i];
    }
}
/*........................................*/
double CartesianR6ZoniShiftedCircular::y(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (r/Rmax) * sin(theta);
}
/*........................................*/
void CartesianR6ZoniShiftedCircular::y(std::vector<double> const& r, double theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (r[i]/Rmax) * sin(theta);
    }
}
/*........................................*/
void CartesianR6ZoniShiftedCircular::y(double r, std::vector<double> const& theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (r/Rmax) * sin_theta[i];
    }
}
/*........................................*/
double CartesianR6ZoniShiftedCircular::J_rr(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (cos(theta))/Rmax;
}
/*........................................*/
void CartesianR6ZoniShiftedCircular::J_rr(std::vector<double> const& r, double theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (cos(theta))/Rmax;
    }
}
/*........................................*/
void CartesianR6ZoniShiftedCircular::J_rr(double r, std::vector<double> const& theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (cos_theta[i])/Rmax;
    }
}
/*........................................*/
double CartesianR6ZoniShiftedCircular::J_rt(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (-(r/Rmax)) * sin(theta);
}
/*........................................*/
void CartesianR6ZoniShiftedCircular::J_rt(std::vector<double> const& r, double theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (-(r[i]/Rmax)) * sin(theta);
    }
}
/*........................................*/
void CartesianR6ZoniShiftedCircular::J_rt(double r, std::vector<double> const& theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (-(r/Rmax)) * sin_theta[i];
    }
}
/*........................................*/
double CartesianR6ZoniShiftedCircular::J_tr(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (sin(theta))/Rmax;
}
/*........................................*/
void CartesianR6ZoniShiftedCircular::J_tr(std::vector<double> const& r, double theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (sin(theta))/Rmax;
    }
}
/*........................................*/
void CartesianR6ZoniShiftedCircular::J_tr(double r, std::vector<double> const& theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (sin_theta[i])/Rmax;
    }
}
/*........................................*/
double CartesianR6ZoniShiftedCircular::J_tt(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (r/Rmax) * cos(theta);
}
/*........................................*/
void CartesianR6ZoniShiftedCircular::J_tt(std::vector<double> const& r, double theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (r[i]/Rmax) * cos(theta);
    }
}
/*........................................*/
void CartesianR6ZoniShiftedCircular::J_tt(double r, std::vector<double> const& theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (r/Rmax) * cos_theta[i];
    }
}
/*........................................*/
double CartesianR6ZoniShiftedCircular::J_xs(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (-pow(sin(theta), 2.0)) / cos(theta) + pow(cos(theta), (double)((-1)));
}
/*........................................*/
void CartesianR6ZoniShiftedCircular::J_xs(std::vector<double> const& r, double theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (-pow(sin(theta), 2.0)) / cos(theta) + pow(cos(theta), (double)((-1)));
    }
}
/*........................................*/
void CartesianR6ZoniShiftedCircular::J_xs(double r, std::vector<double> const& theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (-pow(sin_theta[i], 2.0)) / cos_theta[i] + pow(cos_theta[i], (double)((-1)));
    }
}
/*........................................*/
double CartesianR6ZoniShiftedCircular::J_xt(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return sin(theta);
}
/*........................................*/
void CartesianR6ZoniShiftedCircular::J_xt(std::vector<double> const& r, double theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = sin(theta);
    }
}
/*........................................*/
void CartesianR6ZoniShiftedCircular::J_xt(double r, std::vector<double> const& theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = sin_theta[i];
    }
}
/*........................................*/
double CartesianR6ZoniShiftedCircular::J_ys(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (-sin(theta)) / (r/Rmax);
}
/*........................................*/
void CartesianR6ZoniShiftedCircular::J_ys(std::vector<double> const& r, double theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (-sin(theta)) / (r[i]/Rmax);
    }
}
/*........................................*/
void CartesianR6ZoniShiftedCircular::J_ys(double r, std::vector<double> const& theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (-sin_theta[i]) / (r/Rmax);
    }
}
/*........................................*/
double CartesianR6ZoniShiftedCircular::J_yt(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return cos(theta) / (r/Rmax);
}
/*........................................*/
void CartesianR6ZoniShiftedCircular::J_yt(std::vector<double> const& r, double theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = cos(theta) / (r[i]/Rmax);
    }
}
/*........................................*/
void CartesianR6ZoniShiftedCircular::J_yt(double r, std::vector<double> const& theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = cos_theta[i] / (r/Rmax);
    }
}
/*........................................*/
double CartesianR6ZoniShiftedCircular::rho_glob(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (-((r/Rmax) * (20.0 * pow(tanh(20.0 * (r/Rmax) - 14.0), 2.0) - 20.0) * (0.8192 * M_PI * pow(((r/Rmax) - 1.0), 6.0) * pow(((r/Rmax) + 1.0), 6.0) * sin(theta) * cos(2.0 * M_PI * (r/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r/Rmax) * cos(theta)) - 0.8192 * M_PI * pow(((r/Rmax) - 1.0), 6.0) * pow(((r/Rmax) + 1.0), 6.0) * sin(2.0 * M_PI * (r/Rmax) * sin(theta)) * sin(2.0 * M_PI * (r/Rmax) * cos(theta)) * cos(theta) + 2.4576 * pow(((r/Rmax) - 1.0), 6.0) * pow(((r/Rmax) + 1.0), 5.0) * sin(2.0 * M_PI * (r/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r/Rmax) * cos(theta)) + 2.4576 * pow(((r/Rmax) - 1.0), 5.0) * pow(((r/Rmax) + 1.0), 6.0) * sin(2.0 * M_PI * (r/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r/Rmax) * cos(theta))) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) + (r/Rmax) * ((-1.6384) * (M_PI * M_PI) * pow(((r/Rmax) - 1.0), 6.0) * pow(((r/Rmax) + 1.0), 6.0) * pow(sin(theta), 2.0) * sin(2.0 * M_PI * (r/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r/Rmax) * cos(theta)) - 3.2768 * (M_PI * M_PI) * pow(((r/Rmax) - 1.0), 6.0) * pow(((r/Rmax) + 1.0), 6.0) * sin(theta) * sin(2.0 * M_PI * (r/Rmax) * cos(theta)) * cos(theta) * cos(2.0 * M_PI * (r/Rmax) * sin(theta)) - 1.6384 * (M_PI * M_PI) * pow(((r/Rmax) - 1.0), 6.0) * pow(((r/Rmax) + 1.0), 6.0) * sin(2.0 * M_PI * (r/Rmax) * sin(theta)) * pow(cos(theta), 2.0) * cos(2.0 * M_PI * (r/Rmax) * cos(theta)) + 9.8304 * M_PI * pow(((r/Rmax) - 1.0), 6.0) * pow(((r/Rmax) + 1.0), 5.0) * sin(theta) * cos(2.0 * M_PI * (r/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r/Rmax) * cos(theta)) - 9.8304 * M_PI * pow(((r/Rmax) - 1.0), 6.0) * pow(((r/Rmax) + 1.0), 5.0) * sin(2.0 * M_PI * (r/Rmax) * sin(theta)) * sin(2.0 * M_PI * (r/Rmax) * cos(theta)) * cos(theta) + 12.288 * pow(((r/Rmax) - 1.0), 6.0) * pow(((r/Rmax) + 1.0), 4.0) * sin(2.0 * M_PI * (r/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r/Rmax) * cos(theta)) + 9.8304 * M_PI * pow(((r/Rmax) - 1.0), 5.0) * pow(((r/Rmax) + 1.0), 6.0) * sin(theta) * cos(2.0 * M_PI * (r/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r/Rmax) * cos(theta)) - 9.8304 * M_PI * pow(((r/Rmax) - 1.0), 5.0) * pow(((r/Rmax) + 1.0), 6.0) * sin(2.0 * M_PI * (r/Rmax) * sin(theta)) * sin(2.0 * M_PI * (r/Rmax) * cos(theta)) * cos(theta) + 29.4912 * pow(((r/Rmax) - 1.0), 5.0) * pow(((r/Rmax) + 1.0), 5.0) * sin(2.0 * M_PI * (r/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r/Rmax) * cos(theta)) + 12.288 * pow(((r/Rmax) - 1.0), 4.0) * pow(((r/Rmax) + 1.0), 6.0) * sin(2.0 * M_PI * (r/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r/Rmax) * cos(theta))) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) + (0.8192 * M_PI * pow(((r/Rmax) - 1.0), 6.0) * pow(((r/Rmax) + 1.0), 6.0) * sin(theta) * cos(2.0 * M_PI * (r/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r/Rmax) * cos(theta)) - 0.8192 * M_PI * pow(((r/Rmax) - 1.0), 6.0) * pow(((r/Rmax) + 1.0), 6.0) * sin(2.0 * M_PI * (r/Rmax) * sin(theta)) * sin(2.0 * M_PI * (r/Rmax) * cos(theta)) * cos(theta) + 2.4576 * pow(((r/Rmax) - 1.0), 6.0) * pow(((r/Rmax) + 1.0), 5.0) * sin(2.0 * M_PI * (r/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r/Rmax) * cos(theta)) + 2.4576 * pow(((r/Rmax) - 1.0), 5.0) * pow(((r/Rmax) + 1.0), 6.0) * sin(2.0 * M_PI * (r/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r/Rmax) * cos(theta))) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) + ((-1.6384) * (M_PI * M_PI) * (r/Rmax) * pow(((r/Rmax) - 1.0), 6.0) * pow(((r/Rmax) + 1.0), 6.0) * pow(sin(theta), 2.0) * sin(2.0 * M_PI * (r/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r/Rmax) * cos(theta)) + 3.2768 * (M_PI * M_PI) * (r/Rmax) * pow(((r/Rmax) - 1.0), 6.0) * pow(((r/Rmax) + 1.0), 6.0) * sin(theta) * sin(2.0 * M_PI * (r/Rmax) * cos(theta)) * cos(theta) * cos(2.0 * M_PI * (r/Rmax) * sin(theta)) - 1.6384 * (M_PI * M_PI) * (r/Rmax) * pow(((r/Rmax) - 1.0), 6.0) * pow(((r/Rmax) + 1.0), 6.0) * sin(2.0 * M_PI * (r/Rmax) * sin(theta)) * pow(cos(theta), 2.0) * cos(2.0 * M_PI * (r/Rmax) * cos(theta)) - 0.8192 * M_PI * pow(((r/Rmax) - 1.0), 6.0) * pow(((r/Rmax) + 1.0), 6.0) * sin(theta) * cos(2.0 * M_PI * (r/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r/Rmax) * cos(theta)) + 0.8192 * M_PI * pow(((r/Rmax) - 1.0), 6.0) * pow(((r/Rmax) + 1.0), 6.0) * sin(2.0 * M_PI * (r/Rmax) * sin(theta)) * sin(2.0 * M_PI * (r/Rmax) * cos(theta)) * cos(theta)) * exp(-tanh(20.0 * (r/Rmax) - 14.0)))) / (r/Rmax);
}
/*........................................*/
void CartesianR6ZoniShiftedCircular::rho_glob(std::vector<double> const& r, double theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (-((r[i]/Rmax) * (20.0 * pow(tanh(20.0 * (r[i]/Rmax) - 14.0), 2.0) - 20.0) * (0.8192 * M_PI * pow(((r[i]/Rmax) - 1.0), 6.0) * pow(((r[i]/Rmax) + 1.0), 6.0) * sin(theta) * cos(2.0 * M_PI * (r[i]/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r[i]/Rmax) * cos(theta)) - 0.8192 * M_PI * pow(((r[i]/Rmax) - 1.0), 6.0) * pow(((r[i]/Rmax) + 1.0), 6.0) * sin(2.0 * M_PI * (r[i]/Rmax) * sin(theta)) * sin(2.0 * M_PI * (r[i]/Rmax) * cos(theta)) * cos(theta) + 2.4576 * pow(((r[i]/Rmax) - 1.0), 6.0) * pow(((r[i]/Rmax) + 1.0), 5.0) * sin(2.0 * M_PI * (r[i]/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r[i]/Rmax) * cos(theta)) + 2.4576 * pow(((r[i]/Rmax) - 1.0), 5.0) * pow(((r[i]/Rmax) + 1.0), 6.0) * sin(2.0 * M_PI * (r[i]/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r[i]/Rmax) * cos(theta))) * exp(-tanh(20.0 * (r[i]/Rmax) - 14.0)) + (r[i]/Rmax) * ((-1.6384) * (M_PI * M_PI) * pow(((r[i]/Rmax) - 1.0), 6.0) * pow(((r[i]/Rmax) + 1.0), 6.0) * pow(sin(theta), 2.0) * sin(2.0 * M_PI * (r[i]/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r[i]/Rmax) * cos(theta)) - 3.2768 * (M_PI * M_PI) * pow(((r[i]/Rmax) - 1.0), 6.0) * pow(((r[i]/Rmax) + 1.0), 6.0) * sin(theta) * sin(2.0 * M_PI * (r[i]/Rmax) * cos(theta)) * cos(theta) * cos(2.0 * M_PI * (r[i]/Rmax) * sin(theta)) - 1.6384 * (M_PI * M_PI) * pow(((r[i]/Rmax) - 1.0), 6.0) * pow(((r[i]/Rmax) + 1.0), 6.0) * sin(2.0 * M_PI * (r[i]/Rmax) * sin(theta)) * pow(cos(theta), 2.0) * cos(2.0 * M_PI * (r[i]/Rmax) * cos(theta)) + 9.8304 * M_PI * pow(((r[i]/Rmax) - 1.0), 6.0) * pow(((r[i]/Rmax) + 1.0), 5.0) * sin(theta) * cos(2.0 * M_PI * (r[i]/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r[i]/Rmax) * cos(theta)) - 9.8304 * M_PI * pow(((r[i]/Rmax) - 1.0), 6.0) * pow(((r[i]/Rmax) + 1.0), 5.0) * sin(2.0 * M_PI * (r[i]/Rmax) * sin(theta)) * sin(2.0 * M_PI * (r[i]/Rmax) * cos(theta)) * cos(theta) + 12.288 * pow(((r[i]/Rmax) - 1.0), 6.0) * pow(((r[i]/Rmax) + 1.0), 4.0) * sin(2.0 * M_PI * (r[i]/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r[i]/Rmax) * cos(theta)) + 9.8304 * M_PI * pow(((r[i]/Rmax) - 1.0), 5.0) * pow(((r[i]/Rmax) + 1.0), 6.0) * sin(theta) * cos(2.0 * M_PI * (r[i]/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r[i]/Rmax) * cos(theta)) - 9.8304 * M_PI * pow(((r[i]/Rmax) - 1.0), 5.0) * pow(((r[i]/Rmax) + 1.0), 6.0) * sin(2.0 * M_PI * (r[i]/Rmax) * sin(theta)) * sin(2.0 * M_PI * (r[i]/Rmax) * cos(theta)) * cos(theta) + 29.4912 * pow(((r[i]/Rmax) - 1.0), 5.0) * pow(((r[i]/Rmax) + 1.0), 5.0) * sin(2.0 * M_PI * (r[i]/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r[i]/Rmax) * cos(theta)) + 12.288 * pow(((r[i]/Rmax) - 1.0), 4.0) * pow(((r[i]/Rmax) + 1.0), 6.0) * sin(2.0 * M_PI * (r[i]/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r[i]/Rmax) * cos(theta))) * exp(-tanh(20.0 * (r[i]/Rmax) - 14.0)) + (0.8192 * M_PI * pow(((r[i]/Rmax) - 1.0), 6.0) * pow(((r[i]/Rmax) + 1.0), 6.0) * sin(theta) * cos(2.0 * M_PI * (r[i]/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r[i]/Rmax) * cos(theta)) - 0.8192 * M_PI * pow(((r[i]/Rmax) - 1.0), 6.0) * pow(((r[i]/Rmax) + 1.0), 6.0) * sin(2.0 * M_PI * (r[i]/Rmax) * sin(theta)) * sin(2.0 * M_PI * (r[i]/Rmax) * cos(theta)) * cos(theta) + 2.4576 * pow(((r[i]/Rmax) - 1.0), 6.0) * pow(((r[i]/Rmax) + 1.0), 5.0) * sin(2.0 * M_PI * (r[i]/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r[i]/Rmax) * cos(theta)) + 2.4576 * pow(((r[i]/Rmax) - 1.0), 5.0) * pow(((r[i]/Rmax) + 1.0), 6.0) * sin(2.0 * M_PI * (r[i]/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r[i]/Rmax) * cos(theta))) * exp(-tanh(20.0 * (r[i]/Rmax) - 14.0)) + ((-1.6384) * (M_PI * M_PI) * (r[i]/Rmax) * pow(((r[i]/Rmax) - 1.0), 6.0) * pow(((r[i]/Rmax) + 1.0), 6.0) * pow(sin(theta), 2.0) * sin(2.0 * M_PI * (r[i]/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r[i]/Rmax) * cos(theta)) + 3.2768 * (M_PI * M_PI) * (r[i]/Rmax) * pow(((r[i]/Rmax) - 1.0), 6.0) * pow(((r[i]/Rmax) + 1.0), 6.0) * sin(theta) * sin(2.0 * M_PI * (r[i]/Rmax) * cos(theta)) * cos(theta) * cos(2.0 * M_PI * (r[i]/Rmax) * sin(theta)) - 1.6384 * (M_PI * M_PI) * (r[i]/Rmax) * pow(((r[i]/Rmax) - 1.0), 6.0) * pow(((r[i]/Rmax) + 1.0), 6.0) * sin(2.0 * M_PI * (r[i]/Rmax) * sin(theta)) * pow(cos(theta), 2.0) * cos(2.0 * M_PI * (r[i]/Rmax) * cos(theta)) - 0.8192 * M_PI * pow(((r[i]/Rmax) - 1.0), 6.0) * pow(((r[i]/Rmax) + 1.0), 6.0) * sin(theta) * cos(2.0 * M_PI * (r[i]/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r[i]/Rmax) * cos(theta)) + 0.8192 * M_PI * pow(((r[i]/Rmax) - 1.0), 6.0) * pow(((r[i]/Rmax) + 1.0), 6.0) * sin(2.0 * M_PI * (r[i]/Rmax) * sin(theta)) * sin(2.0 * M_PI * (r[i]/Rmax) * cos(theta)) * cos(theta)) * exp(-tanh(20.0 * (r[i]/Rmax) - 14.0)))) / (r[i]/Rmax);
    }
}
/*........................................*/
void CartesianR6ZoniShiftedCircular::rho_glob(double r, std::vector<double> const& theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (-((r/Rmax) * (20.0 * pow(tanh(20.0 * (r/Rmax) - 14.0), 2.0) - 20.0) * (0.8192 * M_PI * pow(((r/Rmax) - 1.0), 6.0) * pow(((r/Rmax) + 1.0), 6.0) * sin_theta[i] * cos(2.0 * M_PI * (r/Rmax) * sin_theta[i]) * cos(2.0 * M_PI * (r/Rmax) * cos_theta[i]) - 0.8192 * M_PI * pow(((r/Rmax) - 1.0), 6.0) * pow(((r/Rmax) + 1.0), 6.0) * sin(2.0 * M_PI * (r/Rmax) * sin_theta[i]) * sin(2.0 * M_PI * (r/Rmax) * cos_theta[i]) * cos_theta[i] + 2.4576 * pow(((r/Rmax) - 1.0), 6.0) * pow(((r/Rmax) + 1.0), 5.0) * sin(2.0 * M_PI * (r/Rmax) * sin_theta[i]) * cos(2.0 * M_PI * (r/Rmax) * cos_theta[i]) + 2.4576 * pow(((r/Rmax) - 1.0), 5.0) * pow(((r/Rmax) + 1.0), 6.0) * sin(2.0 * M_PI * (r/Rmax) * sin_theta[i]) * cos(2.0 * M_PI * (r/Rmax) * cos_theta[i])) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) + (r/Rmax) * ((-1.6384) * (M_PI * M_PI) * pow(((r/Rmax) - 1.0), 6.0) * pow(((r/Rmax) + 1.0), 6.0) * pow(sin_theta[i], 2.0) * sin(2.0 * M_PI * (r/Rmax) * sin_theta[i]) * cos(2.0 * M_PI * (r/Rmax) * cos_theta[i]) - 3.2768 * (M_PI * M_PI) * pow(((r/Rmax) - 1.0), 6.0) * pow(((r/Rmax) + 1.0), 6.0) * sin_theta[i] * sin(2.0 * M_PI * (r/Rmax) * cos_theta[i]) * cos_theta[i] * cos(2.0 * M_PI * (r/Rmax) * sin_theta[i]) - 1.6384 * (M_PI * M_PI) * pow(((r/Rmax) - 1.0), 6.0) * pow(((r/Rmax) + 1.0), 6.0) * sin(2.0 * M_PI * (r/Rmax) * sin_theta[i]) * pow(cos_theta[i], 2.0) * cos(2.0 * M_PI * (r/Rmax) * cos_theta[i]) + 9.8304 * M_PI * pow(((r/Rmax) - 1.0), 6.0) * pow(((r/Rmax) + 1.0), 5.0) * sin_theta[i] * cos(2.0 * M_PI * (r/Rmax) * sin_theta[i]) * cos(2.0 * M_PI * (r/Rmax) * cos_theta[i]) - 9.8304 * M_PI * pow(((r/Rmax) - 1.0), 6.0) * pow(((r/Rmax) + 1.0), 5.0) * sin(2.0 * M_PI * (r/Rmax) * sin_theta[i]) * sin(2.0 * M_PI * (r/Rmax) * cos_theta[i]) * cos_theta[i] + 12.288 * pow(((r/Rmax) - 1.0), 6.0) * pow(((r/Rmax) + 1.0), 4.0) * sin(2.0 * M_PI * (r/Rmax) * sin_theta[i]) * cos(2.0 * M_PI * (r/Rmax) * cos_theta[i]) + 9.8304 * M_PI * pow(((r/Rmax) - 1.0), 5.0) * pow(((r/Rmax) + 1.0), 6.0) * sin_theta[i] * cos(2.0 * M_PI * (r/Rmax) * sin_theta[i]) * cos(2.0 * M_PI * (r/Rmax) * cos_theta[i]) - 9.8304 * M_PI * pow(((r/Rmax) - 1.0), 5.0) * pow(((r/Rmax) + 1.0), 6.0) * sin(2.0 * M_PI * (r/Rmax) * sin_theta[i]) * sin(2.0 * M_PI * (r/Rmax) * cos_theta[i]) * cos_theta[i] + 29.4912 * pow(((r/Rmax) - 1.0), 5.0) * pow(((r/Rmax) + 1.0), 5.0) * sin(2.0 * M_PI * (r/Rmax) * sin_theta[i]) * cos(2.0 * M_PI * (r/Rmax) * cos_theta[i]) + 12.288 * pow(((r/Rmax) - 1.0), 4.0) * pow(((r/Rmax) + 1.0), 6.0) * sin(2.0 * M_PI * (r/Rmax) * sin_theta[i]) * cos(2.0 * M_PI * (r/Rmax) * cos_theta[i])) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) + (0.8192 * M_PI * pow(((r/Rmax) - 1.0), 6.0) * pow(((r/Rmax) + 1.0), 6.0) * sin_theta[i] * cos(2.0 * M_PI * (r/Rmax) * sin_theta[i]) * cos(2.0 * M_PI * (r/Rmax) * cos_theta[i]) - 0.8192 * M_PI * pow(((r/Rmax) - 1.0), 6.0) * pow(((r/Rmax) + 1.0), 6.0) * sin(2.0 * M_PI * (r/Rmax) * sin_theta[i]) * sin(2.0 * M_PI * (r/Rmax) * cos_theta[i]) * cos_theta[i] + 2.4576 * pow(((r/Rmax) - 1.0), 6.0) * pow(((r/Rmax) + 1.0), 5.0) * sin(2.0 * M_PI * (r/Rmax) * sin_theta[i]) * cos(2.0 * M_PI * (r/Rmax) * cos_theta[i]) + 2.4576 * pow(((r/Rmax) - 1.0), 5.0) * pow(((r/Rmax) + 1.0), 6.0) * sin(2.0 * M_PI * (r/Rmax) * sin_theta[i]) * cos(2.0 * M_PI * (r/Rmax) * cos_theta[i])) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) + ((-1.6384) * (M_PI * M_PI) * (r/Rmax) * pow(((r/Rmax) - 1.0), 6.0) * pow(((r/Rmax) + 1.0), 6.0) * pow(sin_theta[i], 2.0) * sin(2.0 * M_PI * (r/Rmax) * sin_theta[i]) * cos(2.0 * M_PI * (r/Rmax) * cos_theta[i]) + 3.2768 * (M_PI * M_PI) * (r/Rmax) * pow(((r/Rmax) - 1.0), 6.0) * pow(((r/Rmax) + 1.0), 6.0) * sin_theta[i] * sin(2.0 * M_PI * (r/Rmax) * cos_theta[i]) * cos_theta[i] * cos(2.0 * M_PI * (r/Rmax) * sin_theta[i]) - 1.6384 * (M_PI * M_PI) * (r/Rmax) * pow(((r/Rmax) - 1.0), 6.0) * pow(((r/Rmax) + 1.0), 6.0) * sin(2.0 * M_PI * (r/Rmax) * sin_theta[i]) * pow(cos_theta[i], 2.0) * cos(2.0 * M_PI * (r/Rmax) * cos_theta[i]) - 0.8192 * M_PI * pow(((r/Rmax) - 1.0), 6.0) * pow(((r/Rmax) + 1.0), 6.0) * sin_theta[i] * cos(2.0 * M_PI * (r/Rmax) * sin_theta[i]) * cos(2.0 * M_PI * (r/Rmax) * cos_theta[i]) + 0.8192 * M_PI * pow(((r/Rmax) - 1.0), 6.0) * pow(((r/Rmax) + 1.0), 6.0) * sin(2.0 * M_PI * (r/Rmax) * sin_theta[i]) * sin(2.0 * M_PI * (r/Rmax) * cos_theta[i]) * cos_theta[i]) * exp(-tanh(20.0 * (r/Rmax) - 14.0)))) / (r/Rmax);
    }
}
/*........................................*/
double CartesianR6ZoniShiftedCircular::rho_pole(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return 0.0;
}
/*........................................*/
void CartesianR6ZoniShiftedCircular::rho_pole(std::vector<double> const& r, double theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = 0.0;
    }
}
/*........................................*/
void CartesianR6ZoniShiftedCircular::rho_pole(double r, std::vector<double> const& theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = 0.0;
    }
}
/*........................................*/
double CartesianR6ZoniShiftedCircular::coeffs1(double r, double Rmax) const
{
    return exp(-tanh(20.0 * (r/Rmax) - 14.0));
}
/*........................................*/
void CartesianR6ZoniShiftedCircular::coeffs1(std::vector<double> const& r, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = exp(-tanh(20.0 * (r[i]/Rmax) - 14.0));
    }
}
/*........................................*/
double CartesianR6ZoniShiftedCircular::coeffs2(double r, double Rmax) const
{
    return 0.0;
}
/*........................................*/
void CartesianR6ZoniShiftedCircular::coeffs2(std::vector<double> const& r, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = 0.0;
    }
}
/*........................................*/
double CartesianR6ZoniShiftedCircular::phi_exact(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return 0.4096 * pow(((r/Rmax) - 1.0), 6.0) * pow(((r/Rmax) + 1.0), 6.0) * sin(2.0 * M_PI * (r/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r/Rmax) * cos(theta));
}
/*........................................*/
void CartesianR6ZoniShiftedCircular::phi_exact(std::vector<double> const& r, double theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = 0.4096 * pow(((r[i]/Rmax) - 1.0), 6.0) * pow(((r[i]/Rmax) + 1.0), 6.0) * sin(2.0 * M_PI * (r[i]/Rmax) * sin(theta)) * cos(2.0 * M_PI * (r[i]/Rmax) * cos(theta));
    }
}
/*........................................*/
void CartesianR6ZoniShiftedCircular::phi_exact(double r, std::vector<double> const& theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = 0.4096 * pow(((r/Rmax) - 1.0), 6.0) * pow(((r/Rmax) + 1.0), 6.0) * sin(2.0 * M_PI * (r/Rmax) * sin_theta[i]) * cos(2.0 * M_PI * (r/Rmax) * cos_theta[i]);
    }
}
/*........................................*/
