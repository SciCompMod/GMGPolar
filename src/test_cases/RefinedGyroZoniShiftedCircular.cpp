#include "RefinedGyroZoniShiftedCircular.h"
#include <stdlib.h>
#include <math.h>


/*........................................*/
double RefinedGyroZoniShiftedCircular::x(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (r/Rmax) * cos(theta);
}
/*........................................*/
void RefinedGyroZoniShiftedCircular::x(std::vector<double> const& r, double theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (r[i]/Rmax) * cos(theta);
    }
}
/*........................................*/
void RefinedGyroZoniShiftedCircular::x(double r, std::vector<double> const& theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (r/Rmax) * cos_theta[i];
    }
}
/*........................................*/
double RefinedGyroZoniShiftedCircular::y(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (r/Rmax) * sin(theta);
}
/*........................................*/
void RefinedGyroZoniShiftedCircular::y(std::vector<double> const& r, double theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (r[i]/Rmax) * sin(theta);
    }
}
/*........................................*/
void RefinedGyroZoniShiftedCircular::y(double r, std::vector<double> const& theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (r/Rmax) * sin_theta[i];
    }
}
/*........................................*/
double RefinedGyroZoniShiftedCircular::J_rr(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (cos(theta))/Rmax;
}
/*........................................*/
void RefinedGyroZoniShiftedCircular::J_rr(std::vector<double> const& r, double theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (cos(theta))/Rmax;
    }
}
/*........................................*/
void RefinedGyroZoniShiftedCircular::J_rr(double r, std::vector<double> const& theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (cos_theta[i])/Rmax;
    }
}
/*........................................*/
double RefinedGyroZoniShiftedCircular::J_rt(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (-(r/Rmax)) * sin(theta);
}
/*........................................*/
void RefinedGyroZoniShiftedCircular::J_rt(std::vector<double> const& r, double theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (-(r[i]/Rmax)) * sin(theta);
    }
}
/*........................................*/
void RefinedGyroZoniShiftedCircular::J_rt(double r, std::vector<double> const& theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (-(r/Rmax)) * sin_theta[i];
    }
}
/*........................................*/
double RefinedGyroZoniShiftedCircular::J_tr(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (sin(theta))/Rmax;
}
/*........................................*/
void RefinedGyroZoniShiftedCircular::J_tr(std::vector<double> const& r, double theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (sin(theta))/Rmax;
    }
}
/*........................................*/
void RefinedGyroZoniShiftedCircular::J_tr(double r, std::vector<double> const& theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (sin_theta[i])/Rmax;
    }
}
/*........................................*/
double RefinedGyroZoniShiftedCircular::J_tt(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (r/Rmax) * cos(theta);
}
/*........................................*/
void RefinedGyroZoniShiftedCircular::J_tt(std::vector<double> const& r, double theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (r[i]/Rmax) * cos(theta);
    }
}
/*........................................*/
void RefinedGyroZoniShiftedCircular::J_tt(double r, std::vector<double> const& theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (r/Rmax) * cos_theta[i];
    }
}
/*........................................*/
double RefinedGyroZoniShiftedCircular::J_xs(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (-pow(sin(theta), 2.0)) / cos(theta) + pow(cos(theta), (double)((-1)));
}
/*........................................*/
void RefinedGyroZoniShiftedCircular::J_xs(std::vector<double> const& r, double theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (-pow(sin(theta), 2.0)) / cos(theta) + pow(cos(theta), (double)((-1)));
    }
}
/*........................................*/
void RefinedGyroZoniShiftedCircular::J_xs(double r, std::vector<double> const& theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (-pow(sin_theta[i], 2.0)) / cos_theta[i] + pow(cos_theta[i], (double)((-1)));
    }
}
/*........................................*/
double RefinedGyroZoniShiftedCircular::J_xt(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return sin(theta);
}
/*........................................*/
void RefinedGyroZoniShiftedCircular::J_xt(std::vector<double> const& r, double theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = sin(theta);
    }
}
/*........................................*/
void RefinedGyroZoniShiftedCircular::J_xt(double r, std::vector<double> const& theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = sin_theta[i];
    }
}
/*........................................*/
double RefinedGyroZoniShiftedCircular::J_ys(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return (-sin(theta)) / (r/Rmax);
}
/*........................................*/
void RefinedGyroZoniShiftedCircular::J_ys(std::vector<double> const& r, double theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (-sin(theta)) / (r[i]/Rmax);
    }
}
/*........................................*/
void RefinedGyroZoniShiftedCircular::J_ys(double r, std::vector<double> const& theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (-sin_theta[i]) / (r/Rmax);
    }
}
/*........................................*/
double RefinedGyroZoniShiftedCircular::J_yt(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return cos(theta) / (r/Rmax);
}
/*........................................*/
void RefinedGyroZoniShiftedCircular::J_yt(std::vector<double> const& r, double theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = cos(theta) / (r[i]/Rmax);
    }
}
/*........................................*/
void RefinedGyroZoniShiftedCircular::J_yt(double r, std::vector<double> const& theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = cos_theta[i] / (r/Rmax);
    }
}
/*........................................*/
double RefinedGyroZoniShiftedCircular::rho_glob(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return 1.0 * (((-3.33823779536505e-15) * ((r/Rmax) * (r/Rmax)) - 0.0 * (r/Rmax) - 0.0 + exp((-3333.33333333333) * pow(((r/Rmax) - 0.9), 2.0))) * cos(21.0 * theta) + (0.00184273372222541 * ((r/Rmax) * (r/Rmax)) - 0.0018029383826828 * (r/Rmax) - 4.00652973929511e-05 + exp((-50.0) * pow(((r/Rmax) - 0.45), 2.0))) * cos(9.0 * theta)) * exp(tanh(20.0 * (r/Rmax) - 14.0)) - ((r/Rmax) * (((-6.67647559073009e-15) * (r/Rmax) + (6000.0 - 6666.66666666667 * (r/Rmax)) * exp((-3333.33333333333) * pow(((r/Rmax) - 0.9), 2.0))) * cos(21.0 * theta) + (0.00368546744445083 * (r/Rmax) + (45.0 - 100.0 * (r/Rmax)) * exp((-50.0) * pow(((r/Rmax) - 0.45), 2.0)) - 0.0018029383826828) * cos(9.0 * theta)) * (20.0 * pow(tanh(20.0 * (r/Rmax) - 14.0), 2.0) - 20.0) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) + (r/Rmax) * ((10000.0 * pow((0.45 - (r/Rmax)), 2.0) * exp((-50.0) * pow(((r/Rmax) - 0.45), 2.0)) + 0.00368546744445083 - 100.0 * exp((-50.0) * pow(((r/Rmax) - 0.45), 2.0))) * cos(9.0 * theta) + (44444444.4444444 * pow((0.9 - (r/Rmax)), 2.0) * exp((-3333.33333333333) * pow(((r/Rmax) - 0.9), 2.0)) - 6.67647559073009e-15 - 6666.66666666667 * exp((-3333.33333333333) * pow(((r/Rmax) - 0.9), 2.0))) * cos(21.0 * theta)) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) + (((-6.67647559073009e-15) * (r/Rmax) + (6000.0 - 6666.66666666667 * (r/Rmax)) * exp((-3333.33333333333) * pow(((r/Rmax) - 0.9), 2.0))) * cos(21.0 * theta) + (0.00368546744445083 * (r/Rmax) + (45.0 - 100.0 * (r/Rmax)) * exp((-50.0) * pow(((r/Rmax) - 0.45), 2.0)) - 0.0018029383826828) * cos(9.0 * theta)) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) + (21.0 * (7.0102993702666e-14 * ((r/Rmax) * (r/Rmax)) - 21.0 * exp((-3333.33333333333) * pow(((r/Rmax) - 0.9), 2.0))) * cos(21.0 * theta) + 9.0 * ((-0.0165846035000287) * ((r/Rmax) * (r/Rmax)) + 0.0162264454441452 * (r/Rmax) + 0.00036058767653656 - 9.0 * exp((-50.0) * pow(((r/Rmax) - 0.45), 2.0))) * cos(9.0 * theta)) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) / (r/Rmax)) / (r/Rmax);
}
/*........................................*/
void RefinedGyroZoniShiftedCircular::rho_glob(std::vector<double> const& r, double theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = 1.0 * (((-3.33823779536505e-15) * ((r[i]/Rmax) * (r[i]/Rmax)) - 0.0 * (r[i]/Rmax) - 0.0 + exp((-3333.33333333333) * pow(((r[i]/Rmax) - 0.9), 2.0))) * cos(21.0 * theta) + (0.00184273372222541 * ((r[i]/Rmax) * (r[i]/Rmax)) - 0.0018029383826828 * (r[i]/Rmax) - 4.00652973929511e-05 + exp((-50.0) * pow(((r[i]/Rmax) - 0.45), 2.0))) * cos(9.0 * theta)) * exp(tanh(20.0 * (r[i]/Rmax) - 14.0)) - ((r[i]/Rmax) * (((-6.67647559073009e-15) * (r[i]/Rmax) + (6000.0 - 6666.66666666667 * (r[i]/Rmax)) * exp((-3333.33333333333) * pow(((r[i]/Rmax) - 0.9), 2.0))) * cos(21.0 * theta) + (0.00368546744445083 * (r[i]/Rmax) + (45.0 - 100.0 * (r[i]/Rmax)) * exp((-50.0) * pow(((r[i]/Rmax) - 0.45), 2.0)) - 0.0018029383826828) * cos(9.0 * theta)) * (20.0 * pow(tanh(20.0 * (r[i]/Rmax) - 14.0), 2.0) - 20.0) * exp(-tanh(20.0 * (r[i]/Rmax) - 14.0)) + (r[i]/Rmax) * ((10000.0 * pow((0.45 - (r[i]/Rmax)), 2.0) * exp((-50.0) * pow(((r[i]/Rmax) - 0.45), 2.0)) + 0.00368546744445083 - 100.0 * exp((-50.0) * pow(((r[i]/Rmax) - 0.45), 2.0))) * cos(9.0 * theta) + (44444444.4444444 * pow((0.9 - (r[i]/Rmax)), 2.0) * exp((-3333.33333333333) * pow(((r[i]/Rmax) - 0.9), 2.0)) - 6.67647559073009e-15 - 6666.66666666667 * exp((-3333.33333333333) * pow(((r[i]/Rmax) - 0.9), 2.0))) * cos(21.0 * theta)) * exp(-tanh(20.0 * (r[i]/Rmax) - 14.0)) + (((-6.67647559073009e-15) * (r[i]/Rmax) + (6000.0 - 6666.66666666667 * (r[i]/Rmax)) * exp((-3333.33333333333) * pow(((r[i]/Rmax) - 0.9), 2.0))) * cos(21.0 * theta) + (0.00368546744445083 * (r[i]/Rmax) + (45.0 - 100.0 * (r[i]/Rmax)) * exp((-50.0) * pow(((r[i]/Rmax) - 0.45), 2.0)) - 0.0018029383826828) * cos(9.0 * theta)) * exp(-tanh(20.0 * (r[i]/Rmax) - 14.0)) + (21.0 * (7.0102993702666e-14 * ((r[i]/Rmax) * (r[i]/Rmax)) - 21.0 * exp((-3333.33333333333) * pow(((r[i]/Rmax) - 0.9), 2.0))) * cos(21.0 * theta) + 9.0 * ((-0.0165846035000287) * ((r[i]/Rmax) * (r[i]/Rmax)) + 0.0162264454441452 * (r[i]/Rmax) + 0.00036058767653656 - 9.0 * exp((-50.0) * pow(((r[i]/Rmax) - 0.45), 2.0))) * cos(9.0 * theta)) * exp(-tanh(20.0 * (r[i]/Rmax) - 14.0)) / (r[i]/Rmax)) / (r[i]/Rmax);
    }
}
/*........................................*/
void RefinedGyroZoniShiftedCircular::rho_glob(double r, std::vector<double> const& theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = 1.0 * (((-3.33823779536505e-15) * ((r/Rmax) * (r/Rmax)) - 0.0 * (r/Rmax) - 0.0 + exp((-3333.33333333333) * pow(((r/Rmax) - 0.9), 2.0))) * cos(21.0 * theta[i]) + (0.00184273372222541 * ((r/Rmax) * (r/Rmax)) - 0.0018029383826828 * (r/Rmax) - 4.00652973929511e-05 + exp((-50.0) * pow(((r/Rmax) - 0.45), 2.0))) * cos(9.0 * theta[i])) * exp(tanh(20.0 * (r/Rmax) - 14.0)) - ((r/Rmax) * (((-6.67647559073009e-15) * (r/Rmax) + (6000.0 - 6666.66666666667 * (r/Rmax)) * exp((-3333.33333333333) * pow(((r/Rmax) - 0.9), 2.0))) * cos(21.0 * theta[i]) + (0.00368546744445083 * (r/Rmax) + (45.0 - 100.0 * (r/Rmax)) * exp((-50.0) * pow(((r/Rmax) - 0.45), 2.0)) - 0.0018029383826828) * cos(9.0 * theta[i])) * (20.0 * pow(tanh(20.0 * (r/Rmax) - 14.0), 2.0) - 20.0) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) + (r/Rmax) * ((10000.0 * pow((0.45 - (r/Rmax)), 2.0) * exp((-50.0) * pow(((r/Rmax) - 0.45), 2.0)) + 0.00368546744445083 - 100.0 * exp((-50.0) * pow(((r/Rmax) - 0.45), 2.0))) * cos(9.0 * theta[i]) + (44444444.4444444 * pow((0.9 - (r/Rmax)), 2.0) * exp((-3333.33333333333) * pow(((r/Rmax) - 0.9), 2.0)) - 6.67647559073009e-15 - 6666.66666666667 * exp((-3333.33333333333) * pow(((r/Rmax) - 0.9), 2.0))) * cos(21.0 * theta[i])) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) + (((-6.67647559073009e-15) * (r/Rmax) + (6000.0 - 6666.66666666667 * (r/Rmax)) * exp((-3333.33333333333) * pow(((r/Rmax) - 0.9), 2.0))) * cos(21.0 * theta[i]) + (0.00368546744445083 * (r/Rmax) + (45.0 - 100.0 * (r/Rmax)) * exp((-50.0) * pow(((r/Rmax) - 0.45), 2.0)) - 0.0018029383826828) * cos(9.0 * theta[i])) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) + (21.0 * (7.0102993702666e-14 * ((r/Rmax) * (r/Rmax)) - 21.0 * exp((-3333.33333333333) * pow(((r/Rmax) - 0.9), 2.0))) * cos(21.0 * theta[i]) + 9.0 * ((-0.0165846035000287) * ((r/Rmax) * (r/Rmax)) + 0.0162264454441452 * (r/Rmax) + 0.00036058767653656 - 9.0 * exp((-50.0) * pow(((r/Rmax) - 0.45), 2.0))) * cos(9.0 * theta[i])) * exp(-tanh(20.0 * (r/Rmax) - 14.0)) / (r/Rmax)) / (r/Rmax);
    }
}
/*........................................*/
double RefinedGyroZoniShiftedCircular::rho_pole(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return 0.0;
}
/*........................................*/
void RefinedGyroZoniShiftedCircular::rho_pole(std::vector<double> const& r, double theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = 0.0;
    }
}
/*........................................*/
void RefinedGyroZoniShiftedCircular::rho_pole(double r, std::vector<double> const& theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = 0.0;
    }
}
/*........................................*/
double RefinedGyroZoniShiftedCircular::coeffs1(double r, double Rmax) const
{
    return exp(-tanh(20.0 * (r/Rmax) - 14.0));
}
/*........................................*/
void RefinedGyroZoniShiftedCircular::coeffs1(std::vector<double> const& r, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = exp(-tanh(20.0 * (r[i]/Rmax) - 14.0));
    }
}
/*........................................*/
double RefinedGyroZoniShiftedCircular::coeffs2(double r, double Rmax) const
{
    return 1.0 * exp(tanh(20.0 * (r/Rmax) - 14.0));
}
/*........................................*/
void RefinedGyroZoniShiftedCircular::coeffs2(std::vector<double> const& r, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = 1.0 * exp(tanh(20.0 * (r[i]/Rmax) - 14.0));
    }
}
/*........................................*/
double RefinedGyroZoniShiftedCircular::phi_exact(double r, double theta, double unused_1, double unused_2, double Rmax) const
{
    return ((-3.33823779536505e-15) * ((r/Rmax) * (r/Rmax)) - 0.0 * (r/Rmax) - 0.0 + exp((-3333.33333333333) * pow(((r/Rmax) - 0.9), 2.0))) * cos(21.0 * theta) + (0.00184273372222541 * ((r/Rmax) * (r/Rmax)) - 0.0018029383826828 * (r/Rmax) - 4.00652973929511e-05 + exp((-50.0) * pow(((r/Rmax) - 0.45), 2.0))) * cos(9.0 * theta);
}
/*........................................*/
void RefinedGyroZoniShiftedCircular::phi_exact(std::vector<double> const& r, double theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = ((-3.33823779536505e-15) * ((r[i]/Rmax) * (r[i]/Rmax)) - 0.0 * (r[i]/Rmax) - 0.0 + exp((-3333.33333333333) * pow(((r[i]/Rmax) - 0.9), 2.0))) * cos(21.0 * theta) + (0.00184273372222541 * ((r[i]/Rmax) * (r[i]/Rmax)) - 0.0018029383826828 * (r[i]/Rmax) - 4.00652973929511e-05 + exp((-50.0) * pow(((r[i]/Rmax) - 0.45), 2.0))) * cos(9.0 * theta);
    }
}
/*........................................*/
void RefinedGyroZoniShiftedCircular::phi_exact(double r, std::vector<double> const& theta, double unused_1, double unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = ((-3.33823779536505e-15) * ((r/Rmax) * (r/Rmax)) - 0.0 * (r/Rmax) - 0.0 + exp((-3333.33333333333) * pow(((r/Rmax) - 0.9), 2.0))) * cos(21.0 * theta[i]) + (0.00184273372222541 * ((r/Rmax) * (r/Rmax)) - 0.0018029383826828 * (r/Rmax) - 4.00652973929511e-05 + exp((-50.0) * pow(((r/Rmax) - 0.45), 2.0))) * cos(9.0 * theta[i]);
    }
}
/*........................................*/
