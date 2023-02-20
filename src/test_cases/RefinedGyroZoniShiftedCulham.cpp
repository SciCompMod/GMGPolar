#include "RefinedGyroZoniShiftedCulham.h"
#include <stdlib.h>
#include <array>
#include <stdint.h>
#include <math.h>
#include <stdbool.h>


/*........................................*/
double RefinedGyroZoniShiftedCulham::my_sum(std::array<double, 1001>& f, int64_t start_idx, int64_t end_idx) const
{
    int64_t i;
    double result;
    result = 0.0;
    #pragma omp parallel for reduction(+: result)
    for (i = start_idx; i < end_idx; i += 1)
    {
        result += f[i];
    }
    return result;
}
/*........................................*/
double RefinedGyroZoniShiftedCulham::q(double rr) const
{
    return 0.8 - 0.1 * (rr * rr);
}
/*........................................*/
double RefinedGyroZoniShiftedCulham::dq(double rr) const
{
    return (-0.2) * rr;
}
/*........................................*/
double RefinedGyroZoniShiftedCulham::p(double rr) const
{
    return 100000.0 - 90000.0 * (rr * rr);
}
/*........................................*/
double RefinedGyroZoniShiftedCulham::dp(double rr) const
{
    return (-180000.0) * rr;
}
/*........................................*/
double RefinedGyroZoniShiftedCulham::dg(double rr, double g) const
{
    return ((-g) * (0.0625000000000001 * (rr * rr) / pow((1.0 - 0.125 * (rr * rr)), 2.0) + 2.0 / (4.0 - 0.5 * (rr * rr))) + 2.261946711816e-06 * (4.0 - 0.5 * (rr * rr)) / g) / (rr / (4.0 - 0.5 * (rr * rr)) + (4.0 - 0.5 * (rr * rr)) / (g * rr));
}
/*........................................*/
double RefinedGyroZoniShiftedCulham::double_deriv(double rr, double c, double g, double dg, double val, double d_val) const
{
    return c * val / (rr * rr) - d_val * (pow(rr, (double)((-1))) + (4.0 - 0.5 * (rr * rr)) * (2.0 * dg * rr / (4.0 - 0.5 * (rr * rr)) + 0.125 * g * (rr * rr) / pow((1.0 - 0.125 * (rr * rr)), 2.0) + 2.0 * g / (4.0 - 0.5 * (rr * rr))) / (g * rr));
}
/*........................................*/
double RefinedGyroZoniShiftedCulham::g(double rr) const
{
    int64_t ri;
    double dr;
    double m;
    double c;
    ri = (int64_t)(rr * 1000 / 1.0);
    dr = 1.0 / 1000.0;
    if (ri == 1000)
    {
        return g_array[ri];
    }
    else
    {
        m = (g_array[ri + 1] - g_array[ri]) / dr;
        c = g_array[ri] - m * ri * dr;
        return m * rr + c;
    }
}
/*........................................*/
double RefinedGyroZoniShiftedCulham::Delta(double rr) const
{
    int64_t ri;
    double dr;
    double m;
    double c;
    ri = (int64_t)(rr * 1000 / 1.0);
    dr = 1.0 / 1000.0;
    if (ri == 1000)
    {
        return Delta_array[ri];
    }
    else
    {
        m = (Delta_array[ri + 1] - Delta_array[ri]) / dr;
        c = Delta_array[ri] - m * ri * dr;
        return m * rr + c;
    }
}
/*........................................*/
double RefinedGyroZoniShiftedCulham::Delta_prime(double rr) const
{
    int64_t ri;
    double dr;
    double m;
    double c;
    ri = (int64_t)(rr * 1000 / 1.0);
    dr = 1.0 / 1000.0;
    if (ri == 1000)
    {
        return Delta_prime_array[ri];
    }
    else
    {
        m = (Delta_prime_array[ri + 1] - Delta_prime_array[ri]) / dr;
        c = Delta_prime_array[ri] - m * ri * dr;
        return m * rr + c;
    }
}
/*........................................*/
double RefinedGyroZoniShiftedCulham::E(double rr) const
{
    int64_t ri;
    double dr;
    double m;
    double c;
    ri = (int64_t)(rr * 1000 / 1.0);
    dr = 1.0 / 1000.0;
    if (ri == 1000)
    {
        return E_array[ri];
    }
    else
    {
        m = (E_array[ri + 1] - E_array[ri]) / dr;
        c = E_array[ri] - m * ri * dr;
        return m * rr + c;
    }
}
/*........................................*/
double RefinedGyroZoniShiftedCulham::T(double rr) const
{
    int64_t ri;
    double dr;
    double m;
    double c;
    ri = (int64_t)(rr * 1000 / 1.0);
    dr = 1.0 / 1000.0;
    if (ri == 1000)
    {
        return T_array[ri];
    }
    else
    {
        m = (T_array[ri + 1] - T_array[ri]) / dr;
        c = T_array[ri] - m * ri * dr;
        return m * rr + c;
    }
}
/*........................................*/
double RefinedGyroZoniShiftedCulham::E_prime(double rr) const
{
    int64_t ri;
    double dr;
    double m;
    double c;
    ri = (int64_t)(rr * 1000 / 1.0);
    dr = 1.0 / 1000.0;
    if (ri == 1000)
    {
        return E_prime_array[ri];
    }
    else
    {
        m = (E_prime_array[ri + 1] - E_prime_array[ri]) / dr;
        c = E_prime_array[ri] - m * ri * dr;
        return m * rr + c;
    }
}
/*........................................*/
double RefinedGyroZoniShiftedCulham::T_prime(double rr) const
{
    int64_t ri;
    double dr;
    double m;
    double c;
    ri = (int64_t)(rr * 1000 / 1.0);
    dr = 1.0 / 1000.0;
    if (ri == 1000)
    {
        return T_prime_array[ri];
    }
    else
    {
        m = (T_prime_array[ri + 1] - T_prime_array[ri]) / dr;
        c = T_prime_array[ri] - m * ri * dr;
        return m * rr + c;
    }
}
/*........................................*/
double RefinedGyroZoniShiftedCulham::P(double rr) const
{
    if (rr == 0)
    {
        return 0.0;
    }
    else
    {
        return 0.005 * pow(rr, 3.0) + 0.1 * rr * Delta(rr) - 0.1 * pow(E(rr), 2.0) - pow(T(rr), 2.0) / rr;
    }
}
/*........................................*/
double RefinedGyroZoniShiftedCulham::dP(double rr) const
{
    if (rr == 0)
    {
        return 0.0;
    }
    else
    {
        return 0.015 * (rr * rr) + 0.1 * rr * Delta_prime(rr) + 0.1 * Delta(rr) - 0.2 * E(rr) * E_prime(rr) - 2.0 * T(rr) * T_prime(rr) / rr + pow(T(rr), 2.0) / (rr * rr);
    }
}
/*........................................*/
RefinedGyroZoniShiftedCulham::RefinedGyroZoniShiftedCulham()
{
    double rr;
    double dr;
    double dr_h;
    std::array<double, 1001> r;
    int64_t i;
    double dg_1;
    double dE_1;
    double dT_1;
    double ddE_1;
    double ddT_1;
    double r2;
    double g_2;
    double dg_2;
    double E_2;
    double T_2;
    double dE_2;
    double dT_2;
    double ddE_2;
    double ddT_2;
    double g_3;
    double dg_3;
    double E_3;
    double T_3;
    double dE_3;
    double dT_3;
    double ddE_3;
    double ddT_3;
    double g_4;
    double dg_4;
    double E_4;
    double T_4;
    double dE_4;
    double dT_4;
    double ddE_4;
    double ddT_4;
    double current_Ea;
    double current_Ta;
    std::array<double, 1001> f;
    std::array<double, 1001> integ_contents;
    double integral;
    double current_Delta_a;
    size_t i_0001;
    rr = 0.0;
    dr = 1.0 / 1000;
    dr_h = dr * 0.5;
    g_array[0] = 1.0;
    E_array[0] = 0.0;
    E_prime_array[0] = 1.0;
    T_array[0] = 0.0;
    T_prime_array[0] = 0.0;
    r[0] = rr;
    rr += dr;
    r[1] = rr;
    g_array[1] = 1.0;
    E_array[1] = rr;
    E_prime_array[1] = 1.0;
    T_array[1] = rr * rr;
    T_prime_array[1] = 2 * rr;
    for (i = 1; i < 1000; i += 1)
    {
        /*Step 1*/
        dg_1 = dg(rr, g_array[i]);
        dE_1 = E_prime_array[i];
        dT_1 = T_prime_array[i];
        ddE_1 = double_deriv(rr, 3.0, g_array[i], dg_1, E_array[i], E_prime_array[i]);
        ddT_1 = double_deriv(rr, 8.0, g_array[i], dg_1, T_array[i], T_prime_array[i]);
        /*Step 2*/
        r2 = rr + dr_h;
        g_2 = g_array[i] + dr_h * dg_1;
        dg_2 = dg(r2, g_2);
        E_2 = E_array[i] + dE_1 * dr_h;
        T_2 = T_array[i] + dT_1 * dr_h;
        dE_2 = E_prime_array[i] + ddE_1 * dr_h;
        dT_2 = T_prime_array[i] + ddT_1 * dr_h;
        ddE_2 = double_deriv(r2, 3.0, g_2, dg_2, E_2, dE_2);
        ddT_2 = double_deriv(r2, 8.0, g_2, dg_2, T_2, dT_2);
        /*Step 3*/
        g_3 = g_array[i] + dr_h * dg_2;
        dg_3 = dg(r2, g_3);
        E_3 = E_array[i] + dE_2 * dr_h;
        T_3 = T_array[i] + dT_2 * dr_h;
        dE_3 = E_prime_array[i] + ddE_2 * dr_h;
        dT_3 = T_prime_array[i] + ddT_2 * dr_h;
        ddE_3 = double_deriv(r2, 3.0, g_3, dg_3, E_3, dE_3);
        ddT_3 = double_deriv(r2, 8.0, g_3, dg_3, T_3, dT_3);
        /*Step 4*/
        rr = rr + dr;
        g_4 = g_array[i] + dr * dg_3;
        dg_4 = dg(rr, g_4);
        E_4 = E_array[i] + dE_3 * dr;
        T_4 = T_array[i] + dT_3 * dr;
        dE_4 = E_prime_array[i] + ddE_3 * dr;
        dT_4 = T_prime_array[i] + ddT_3 * dr;
        ddE_4 = double_deriv(rr, 3.0, g_4, dg_4, E_4, dE_4);
        ddT_4 = double_deriv(rr, 8.0, g_4, dg_4, T_4, dT_4);
        g_array[i + 1] = g_array[i] + dr * (dg_1 + 2 * dg_2 + 2 * dg_3 + dg_4) / 6.0;
        E_array[i + 1] = E_array[i] + dr * (dE_1 + 2 * dE_2 + 2 * dE_3 + dE_4) / 6.0;
        T_array[i + 1] = T_array[i] + dr * (dT_1 + 2 * dT_2 + 2 * dT_3 + dT_4) / 6.0;
        E_prime_array[i + 1] = E_prime_array[i] + dr * (ddE_1 + 2 * ddE_2 + 2 * ddE_3 + ddE_4) / 6.0;
        T_prime_array[i + 1] = T_prime_array[i] + dr * (ddT_1 + 2 * ddT_2 + 2 * ddT_3 + ddT_4) / 6.0;
        r[i + 1] = rr;
    }
    current_Ea = E(1.0);
    current_Ta = T(1.0);
    for (i_0001 = 0; i_0001 < E_array.size(); i_0001 += 1)
    {
        E_array[i_0001] = 0.25 * E_array[i_0001] / current_Ea;
    }
    for (i_0001 = 0; i_0001 < T_array.size(); i_0001 += 1)
    {
        T_array[i_0001] = 0.1 * T_array[i_0001] / current_Ta;
    }
    for (i_0001 = 0; i_0001 < E_prime_array.size(); i_0001 += 1)
    {
        E_prime_array[i_0001] = 0.25 * E_prime_array[i_0001] / current_Ea;
    }
    for (i_0001 = 0; i_0001 < T_prime_array.size(); i_0001 += 1)
    {
        T_prime_array[i_0001] = 0.1 * T_prime_array[i_0001] / current_Ta;
    }
    for (i_0001 = 0; i_0001 < f.size(); i_0001 += 1)
    {
        f[i_0001] = r[i_0001] * g_array[i_0001] / 5.0 / q(r[i_0001]);
    }
    for (i_0001 = 0; i_0001 < integ_contents.size(); i_0001 += 1)
    {
        integ_contents[i_0001] = r[i_0001] * (f[i_0001] * f[i_0001]) - 2 * (r[i_0001] * r[i_0001]) * 1.25663706212e-06 * dp(r[i_0001]) / (1.0 * 1.0);
    }
    Delta_prime_array[0] = 0;
    Delta_array[0] = 0;
    integral = 0.0;
    for (i = 1; i < 1000 + 1; i += 1)
    {
        integral += dr * (integ_contents[i - 1] + integ_contents[i]) * 0.5;
        Delta_prime_array[i] = (-integral) / (5.0 * r[i] * pow(f[i], 2.0));
        Delta_array[i] = Delta_array[i - 1] + dr * 0.5 * (Delta_prime_array[i - 1] + Delta_prime_array[i]);
    }
    current_Delta_a = Delta(1.0);
    for (i_0001 = 0; i_0001 < Delta_array.size(); i_0001 += 1)
    {
        Delta_array[i_0001] -= current_Delta_a;
    }
}
/*........................................*/
double RefinedGyroZoniShiftedCulham::x(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const
{
    return (r/Rmax) * cos(theta) + Delta((r/Rmax)) - E((r/Rmax)) * cos(theta) - P((r/Rmax)) * cos(theta) + T((r/Rmax)) * cos(2.0 * theta) + 5.0;
}
/*........................................*/
void RefinedGyroZoniShiftedCulham::x(std::vector<double> const& r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (r[i]/Rmax) * cos(theta) + Delta((r[i]/Rmax)) - E((r[i]/Rmax)) * cos(theta) - P((r[i]/Rmax)) * cos(theta) + T((r[i]/Rmax)) * cos(2.0 * theta) + 5.0;
    }
}
/*........................................*/
void RefinedGyroZoniShiftedCulham::x(double r, std::vector<double> const& theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (r/Rmax) * cos_theta[i] + Delta((r/Rmax)) - E((r/Rmax)) * cos_theta[i] - P((r/Rmax)) * cos_theta[i] + T((r/Rmax)) * cos(2.0 * theta[i]) + 5.0;
    }
}
/*........................................*/
double RefinedGyroZoniShiftedCulham::y(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const
{
    return (r/Rmax) * sin(theta) - E((r/Rmax)) * sin(theta) - P((r/Rmax)) * sin(theta) - T((r/Rmax)) * sin(2.0 * theta);
}
/*........................................*/
void RefinedGyroZoniShiftedCulham::y(std::vector<double> const& r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (r[i]/Rmax) * sin(theta) - E((r[i]/Rmax)) * sin(theta) - P((r[i]/Rmax)) * sin(theta) - T((r[i]/Rmax)) * sin(2.0 * theta);
    }
}
/*........................................*/
void RefinedGyroZoniShiftedCulham::y(double r, std::vector<double> const& theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (r/Rmax) * sin_theta[i] - E((r/Rmax)) * sin_theta[i] - P((r/Rmax)) * sin_theta[i] - T((r/Rmax)) * sin(2.0 * theta[i]);
    }
}
/*........................................*/
double RefinedGyroZoniShiftedCulham::J_rr(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const
{
    return (Delta_prime((r/Rmax)) - E_prime((r/Rmax)) * cos(theta) + T_prime((r/Rmax)) * cos(2.0 * theta) - dP((r/Rmax)) * cos(theta) + cos(theta))/Rmax;
}
/*........................................*/
void RefinedGyroZoniShiftedCulham::J_rr(std::vector<double> const& r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (Delta_prime((r[i]/Rmax)) - E_prime((r[i]/Rmax)) * cos(theta) + T_prime((r[i]/Rmax)) * cos(2.0 * theta) - dP((r[i]/Rmax)) * cos(theta) + cos(theta))/Rmax;
    }
}
/*........................................*/
void RefinedGyroZoniShiftedCulham::J_rr(double r, std::vector<double> const& theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (Delta_prime((r/Rmax)) - E_prime((r/Rmax)) * cos_theta[i] + T_prime((r/Rmax)) * cos(2.0 * theta[i]) - dP((r/Rmax)) * cos_theta[i] + cos_theta[i])/Rmax;
    }
}
/*........................................*/
double RefinedGyroZoniShiftedCulham::J_rt(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const
{
    return (-(r/Rmax)) * sin(theta) + E((r/Rmax)) * sin(theta) + P((r/Rmax)) * sin(theta) - 2.0 * T((r/Rmax)) * sin(2.0 * theta);
}
/*........................................*/
void RefinedGyroZoniShiftedCulham::J_rt(std::vector<double> const& r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (-(r[i]/Rmax)) * sin(theta) + E((r[i]/Rmax)) * sin(theta) + P((r[i]/Rmax)) * sin(theta) - 2.0 * T((r[i]/Rmax)) * sin(2.0 * theta);
    }
}
/*........................................*/
void RefinedGyroZoniShiftedCulham::J_rt(double r, std::vector<double> const& theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (-(r/Rmax)) * sin_theta[i] + E((r/Rmax)) * sin_theta[i] + P((r/Rmax)) * sin_theta[i] - 2.0 * T((r/Rmax)) * sin(2.0 * theta[i]);
    }
}
/*........................................*/
double RefinedGyroZoniShiftedCulham::J_tr(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const
{
    return ((-E_prime((r/Rmax))) * sin(theta) - T_prime((r/Rmax)) * sin(2.0 * theta) - dP((r/Rmax)) * sin(theta) + sin(theta))/Rmax;
}
/*........................................*/
void RefinedGyroZoniShiftedCulham::J_tr(std::vector<double> const& r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = ((-E_prime((r[i]/Rmax))) * sin(theta) - T_prime((r[i]/Rmax)) * sin(2.0 * theta) - dP((r[i]/Rmax)) * sin(theta) + sin(theta))/Rmax;
    }
}
/*........................................*/
void RefinedGyroZoniShiftedCulham::J_tr(double r, std::vector<double> const& theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = ((-E_prime((r/Rmax))) * sin_theta[i] - T_prime((r/Rmax)) * sin(2.0 * theta[i]) - dP((r/Rmax)) * sin_theta[i] + sin_theta[i])/Rmax;
    }
}
/*........................................*/
double RefinedGyroZoniShiftedCulham::J_tt(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const
{
    return (r/Rmax) * cos(theta) - E((r/Rmax)) * cos(theta) - P((r/Rmax)) * cos(theta) - 2.0 * T((r/Rmax)) * cos(2.0 * theta);
}
/*........................................*/
void RefinedGyroZoniShiftedCulham::J_tt(std::vector<double> const& r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (r[i]/Rmax) * cos(theta) - E((r[i]/Rmax)) * cos(theta) - P((r[i]/Rmax)) * cos(theta) - 2.0 * T((r[i]/Rmax)) * cos(2.0 * theta);
    }
}
/*........................................*/
void RefinedGyroZoniShiftedCulham::J_tt(double r, std::vector<double> const& theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = (r/Rmax) * cos_theta[i] - E((r/Rmax)) * cos_theta[i] - P((r/Rmax)) * cos_theta[i] - 2.0 * T((r/Rmax)) * cos(2.0 * theta[i]);
    }
}
/*........................................*/
double RefinedGyroZoniShiftedCulham::J_xr(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const
{
    return ((-(r/Rmax)) * cos(theta) + E((r/Rmax)) * cos(theta) + P((r/Rmax)) * cos(theta) + 2.0 * T((r/Rmax)) * cos(2.0 * theta)) / ((-(r/Rmax)) * Delta_prime((r/Rmax)) * cos(theta) + (r/Rmax) * E_prime((r/Rmax)) * pow(sin(theta), 2.0) + (r/Rmax) * E_prime((r/Rmax)) * pow(cos(theta), 2.0) + (r/Rmax) * T_prime((r/Rmax)) * sin(theta) * sin(2.0 * theta) - (r/Rmax) * T_prime((r/Rmax)) * cos(theta) * cos(2.0 * theta) + (r/Rmax) * dP((r/Rmax)) * pow(sin(theta), 2.0) + (r/Rmax) * dP((r/Rmax)) * pow(cos(theta), 2.0) - (r/Rmax) * pow(sin(theta), 2.0) - (r/Rmax) * pow(cos(theta), 2.0) + Delta_prime((r/Rmax)) * E((r/Rmax)) * cos(theta) + Delta_prime((r/Rmax)) * P((r/Rmax)) * cos(theta) + 2.0 * Delta_prime((r/Rmax)) * T((r/Rmax)) * cos(2.0 * theta) - E((r/Rmax)) * E_prime((r/Rmax)) * pow(sin(theta), 2.0) - E((r/Rmax)) * E_prime((r/Rmax)) * pow(cos(theta), 2.0) - E((r/Rmax)) * T_prime((r/Rmax)) * sin(theta) * sin(2.0 * theta) + E((r/Rmax)) * T_prime((r/Rmax)) * cos(theta) * cos(2.0 * theta) - E((r/Rmax)) * dP((r/Rmax)) * pow(sin(theta), 2.0) - E((r/Rmax)) * dP((r/Rmax)) * pow(cos(theta), 2.0) + E((r/Rmax)) * pow(sin(theta), 2.0) + E((r/Rmax)) * pow(cos(theta), 2.0) - E_prime((r/Rmax)) * P((r/Rmax)) * pow(sin(theta), 2.0) - E_prime((r/Rmax)) * P((r/Rmax)) * pow(cos(theta), 2.0) + 2.0 * E_prime((r/Rmax)) * T((r/Rmax)) * sin(theta) * sin(2.0 * theta) - 2.0 * E_prime((r/Rmax)) * T((r/Rmax)) * cos(theta) * cos(2.0 * theta) - P((r/Rmax)) * T_prime((r/Rmax)) * sin(theta) * sin(2.0 * theta) + P((r/Rmax)) * T_prime((r/Rmax)) * cos(theta) * cos(2.0 * theta) - P((r/Rmax)) * dP((r/Rmax)) * pow(sin(theta), 2.0) - P((r/Rmax)) * dP((r/Rmax)) * pow(cos(theta), 2.0) + P((r/Rmax)) * pow(sin(theta), 2.0) + P((r/Rmax)) * pow(cos(theta), 2.0) + 2.0 * T((r/Rmax)) * T_prime((r/Rmax)) * pow(sin(2.0 * theta), 2.0) + 2.0 * T((r/Rmax)) * T_prime((r/Rmax)) * pow(cos(2.0 * theta), 2.0) + 2.0 * T((r/Rmax)) * dP((r/Rmax)) * sin(theta) * sin(2.0 * theta) - 2.0 * T((r/Rmax)) * dP((r/Rmax)) * cos(theta) * cos(2.0 * theta) - 2.0 * T((r/Rmax)) * sin(theta) * sin(2.0 * theta) + 2.0 * T((r/Rmax)) * cos(theta) * cos(2.0 * theta));
}
/*........................................*/
void RefinedGyroZoniShiftedCulham::J_xr(std::vector<double> const& r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = ((-(r[i]/Rmax)) * cos(theta) + E((r[i]/Rmax)) * cos(theta) + P((r[i]/Rmax)) * cos(theta) + 2.0 * T((r[i]/Rmax)) * cos(2.0 * theta)) / ((-(r[i]/Rmax)) * Delta_prime((r[i]/Rmax)) * cos(theta) + (r[i]/Rmax) * E_prime((r[i]/Rmax)) * pow(sin(theta), 2.0) + (r[i]/Rmax) * E_prime((r[i]/Rmax)) * pow(cos(theta), 2.0) + (r[i]/Rmax) * T_prime((r[i]/Rmax)) * sin(theta) * sin(2.0 * theta) - (r[i]/Rmax) * T_prime((r[i]/Rmax)) * cos(theta) * cos(2.0 * theta) + (r[i]/Rmax) * dP((r[i]/Rmax)) * pow(sin(theta), 2.0) + (r[i]/Rmax) * dP((r[i]/Rmax)) * pow(cos(theta), 2.0) - (r[i]/Rmax) * pow(sin(theta), 2.0) - (r[i]/Rmax) * pow(cos(theta), 2.0) + Delta_prime((r[i]/Rmax)) * E((r[i]/Rmax)) * cos(theta) + Delta_prime((r[i]/Rmax)) * P((r[i]/Rmax)) * cos(theta) + 2.0 * Delta_prime((r[i]/Rmax)) * T((r[i]/Rmax)) * cos(2.0 * theta) - E((r[i]/Rmax)) * E_prime((r[i]/Rmax)) * pow(sin(theta), 2.0) - E((r[i]/Rmax)) * E_prime((r[i]/Rmax)) * pow(cos(theta), 2.0) - E((r[i]/Rmax)) * T_prime((r[i]/Rmax)) * sin(theta) * sin(2.0 * theta) + E((r[i]/Rmax)) * T_prime((r[i]/Rmax)) * cos(theta) * cos(2.0 * theta) - E((r[i]/Rmax)) * dP((r[i]/Rmax)) * pow(sin(theta), 2.0) - E((r[i]/Rmax)) * dP((r[i]/Rmax)) * pow(cos(theta), 2.0) + E((r[i]/Rmax)) * pow(sin(theta), 2.0) + E((r[i]/Rmax)) * pow(cos(theta), 2.0) - E_prime((r[i]/Rmax)) * P((r[i]/Rmax)) * pow(sin(theta), 2.0) - E_prime((r[i]/Rmax)) * P((r[i]/Rmax)) * pow(cos(theta), 2.0) + 2.0 * E_prime((r[i]/Rmax)) * T((r[i]/Rmax)) * sin(theta) * sin(2.0 * theta) - 2.0 * E_prime((r[i]/Rmax)) * T((r[i]/Rmax)) * cos(theta) * cos(2.0 * theta) - P((r[i]/Rmax)) * T_prime((r[i]/Rmax)) * sin(theta) * sin(2.0 * theta) + P((r[i]/Rmax)) * T_prime((r[i]/Rmax)) * cos(theta) * cos(2.0 * theta) - P((r[i]/Rmax)) * dP((r[i]/Rmax)) * pow(sin(theta), 2.0) - P((r[i]/Rmax)) * dP((r[i]/Rmax)) * pow(cos(theta), 2.0) + P((r[i]/Rmax)) * pow(sin(theta), 2.0) + P((r[i]/Rmax)) * pow(cos(theta), 2.0) + 2.0 * T((r[i]/Rmax)) * T_prime((r[i]/Rmax)) * pow(sin(2.0 * theta), 2.0) + 2.0 * T((r[i]/Rmax)) * T_prime((r[i]/Rmax)) * pow(cos(2.0 * theta), 2.0) + 2.0 * T((r[i]/Rmax)) * dP((r[i]/Rmax)) * sin(theta) * sin(2.0 * theta) - 2.0 * T((r[i]/Rmax)) * dP((r[i]/Rmax)) * cos(theta) * cos(2.0 * theta) - 2.0 * T((r[i]/Rmax)) * sin(theta) * sin(2.0 * theta) + 2.0 * T((r[i]/Rmax)) * cos(theta) * cos(2.0 * theta));
    }
}
/*........................................*/
void RefinedGyroZoniShiftedCulham::J_xr(double r, std::vector<double> const& theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = ((-(r/Rmax)) * cos_theta[i] + E((r/Rmax)) * cos_theta[i] + P((r/Rmax)) * cos_theta[i] + 2.0 * T((r/Rmax)) * cos(2.0 * theta[i])) / ((-(r/Rmax)) * Delta_prime((r/Rmax)) * cos_theta[i] + (r/Rmax) * E_prime((r/Rmax)) * pow(sin_theta[i], 2.0) + (r/Rmax) * E_prime((r/Rmax)) * pow(cos_theta[i], 2.0) + (r/Rmax) * T_prime((r/Rmax)) * sin_theta[i] * sin(2.0 * theta[i]) - (r/Rmax) * T_prime((r/Rmax)) * cos_theta[i] * cos(2.0 * theta[i]) + (r/Rmax) * dP((r/Rmax)) * pow(sin_theta[i], 2.0) + (r/Rmax) * dP((r/Rmax)) * pow(cos_theta[i], 2.0) - (r/Rmax) * pow(sin_theta[i], 2.0) - (r/Rmax) * pow(cos_theta[i], 2.0) + Delta_prime((r/Rmax)) * E((r/Rmax)) * cos_theta[i] + Delta_prime((r/Rmax)) * P((r/Rmax)) * cos_theta[i] + 2.0 * Delta_prime((r/Rmax)) * T((r/Rmax)) * cos(2.0 * theta[i]) - E((r/Rmax)) * E_prime((r/Rmax)) * pow(sin_theta[i], 2.0) - E((r/Rmax)) * E_prime((r/Rmax)) * pow(cos_theta[i], 2.0) - E((r/Rmax)) * T_prime((r/Rmax)) * sin_theta[i] * sin(2.0 * theta[i]) + E((r/Rmax)) * T_prime((r/Rmax)) * cos_theta[i] * cos(2.0 * theta[i]) - E((r/Rmax)) * dP((r/Rmax)) * pow(sin_theta[i], 2.0) - E((r/Rmax)) * dP((r/Rmax)) * pow(cos_theta[i], 2.0) + E((r/Rmax)) * pow(sin_theta[i], 2.0) + E((r/Rmax)) * pow(cos_theta[i], 2.0) - E_prime((r/Rmax)) * P((r/Rmax)) * pow(sin_theta[i], 2.0) - E_prime((r/Rmax)) * P((r/Rmax)) * pow(cos_theta[i], 2.0) + 2.0 * E_prime((r/Rmax)) * T((r/Rmax)) * sin_theta[i] * sin(2.0 * theta[i]) - 2.0 * E_prime((r/Rmax)) * T((r/Rmax)) * cos_theta[i] * cos(2.0 * theta[i]) - P((r/Rmax)) * T_prime((r/Rmax)) * sin_theta[i] * sin(2.0 * theta[i]) + P((r/Rmax)) * T_prime((r/Rmax)) * cos_theta[i] * cos(2.0 * theta[i]) - P((r/Rmax)) * dP((r/Rmax)) * pow(sin_theta[i], 2.0) - P((r/Rmax)) * dP((r/Rmax)) * pow(cos_theta[i], 2.0) + P((r/Rmax)) * pow(sin_theta[i], 2.0) + P((r/Rmax)) * pow(cos_theta[i], 2.0) + 2.0 * T((r/Rmax)) * T_prime((r/Rmax)) * pow(sin(2.0 * theta[i]), 2.0) + 2.0 * T((r/Rmax)) * T_prime((r/Rmax)) * pow(cos(2.0 * theta[i]), 2.0) + 2.0 * T((r/Rmax)) * dP((r/Rmax)) * sin_theta[i] * sin(2.0 * theta[i]) - 2.0 * T((r/Rmax)) * dP((r/Rmax)) * cos_theta[i] * cos(2.0 * theta[i]) - 2.0 * T((r/Rmax)) * sin_theta[i] * sin(2.0 * theta[i]) + 2.0 * T((r/Rmax)) * cos_theta[i] * cos(2.0 * theta[i]));
    }
}
/*........................................*/
double RefinedGyroZoniShiftedCulham::J_xq(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const
{
    return ((-(r/Rmax)) * sin(theta) + E((r/Rmax)) * sin(theta) + P((r/Rmax)) * sin(theta) - 2.0 * T((r/Rmax)) * sin(2.0 * theta)) / ((-(r/Rmax)) * Delta_prime((r/Rmax)) * cos(theta) + (r/Rmax) * E_prime((r/Rmax)) * pow(sin(theta), 2.0) + (r/Rmax) * E_prime((r/Rmax)) * pow(cos(theta), 2.0) + (r/Rmax) * T_prime((r/Rmax)) * sin(theta) * sin(2.0 * theta) - (r/Rmax) * T_prime((r/Rmax)) * cos(theta) * cos(2.0 * theta) + (r/Rmax) * dP((r/Rmax)) * pow(sin(theta), 2.0) + (r/Rmax) * dP((r/Rmax)) * pow(cos(theta), 2.0) - (r/Rmax) * pow(sin(theta), 2.0) - (r/Rmax) * pow(cos(theta), 2.0) + Delta_prime((r/Rmax)) * E((r/Rmax)) * cos(theta) + Delta_prime((r/Rmax)) * P((r/Rmax)) * cos(theta) + 2.0 * Delta_prime((r/Rmax)) * T((r/Rmax)) * cos(2.0 * theta) - E((r/Rmax)) * E_prime((r/Rmax)) * pow(sin(theta), 2.0) - E((r/Rmax)) * E_prime((r/Rmax)) * pow(cos(theta), 2.0) - E((r/Rmax)) * T_prime((r/Rmax)) * sin(theta) * sin(2.0 * theta) + E((r/Rmax)) * T_prime((r/Rmax)) * cos(theta) * cos(2.0 * theta) - E((r/Rmax)) * dP((r/Rmax)) * pow(sin(theta), 2.0) - E((r/Rmax)) * dP((r/Rmax)) * pow(cos(theta), 2.0) + E((r/Rmax)) * pow(sin(theta), 2.0) + E((r/Rmax)) * pow(cos(theta), 2.0) - E_prime((r/Rmax)) * P((r/Rmax)) * pow(sin(theta), 2.0) - E_prime((r/Rmax)) * P((r/Rmax)) * pow(cos(theta), 2.0) + 2.0 * E_prime((r/Rmax)) * T((r/Rmax)) * sin(theta) * sin(2.0 * theta) - 2.0 * E_prime((r/Rmax)) * T((r/Rmax)) * cos(theta) * cos(2.0 * theta) - P((r/Rmax)) * T_prime((r/Rmax)) * sin(theta) * sin(2.0 * theta) + P((r/Rmax)) * T_prime((r/Rmax)) * cos(theta) * cos(2.0 * theta) - P((r/Rmax)) * dP((r/Rmax)) * pow(sin(theta), 2.0) - P((r/Rmax)) * dP((r/Rmax)) * pow(cos(theta), 2.0) + P((r/Rmax)) * pow(sin(theta), 2.0) + P((r/Rmax)) * pow(cos(theta), 2.0) + 2.0 * T((r/Rmax)) * T_prime((r/Rmax)) * pow(sin(2.0 * theta), 2.0) + 2.0 * T((r/Rmax)) * T_prime((r/Rmax)) * pow(cos(2.0 * theta), 2.0) + 2.0 * T((r/Rmax)) * dP((r/Rmax)) * sin(theta) * sin(2.0 * theta) - 2.0 * T((r/Rmax)) * dP((r/Rmax)) * cos(theta) * cos(2.0 * theta) - 2.0 * T((r/Rmax)) * sin(theta) * sin(2.0 * theta) + 2.0 * T((r/Rmax)) * cos(theta) * cos(2.0 * theta));
}
/*........................................*/
void RefinedGyroZoniShiftedCulham::J_xq(std::vector<double> const& r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = ((-(r[i]/Rmax)) * sin(theta) + E((r[i]/Rmax)) * sin(theta) + P((r[i]/Rmax)) * sin(theta) - 2.0 * T((r[i]/Rmax)) * sin(2.0 * theta)) / ((-(r[i]/Rmax)) * Delta_prime((r[i]/Rmax)) * cos(theta) + (r[i]/Rmax) * E_prime((r[i]/Rmax)) * pow(sin(theta), 2.0) + (r[i]/Rmax) * E_prime((r[i]/Rmax)) * pow(cos(theta), 2.0) + (r[i]/Rmax) * T_prime((r[i]/Rmax)) * sin(theta) * sin(2.0 * theta) - (r[i]/Rmax) * T_prime((r[i]/Rmax)) * cos(theta) * cos(2.0 * theta) + (r[i]/Rmax) * dP((r[i]/Rmax)) * pow(sin(theta), 2.0) + (r[i]/Rmax) * dP((r[i]/Rmax)) * pow(cos(theta), 2.0) - (r[i]/Rmax) * pow(sin(theta), 2.0) - (r[i]/Rmax) * pow(cos(theta), 2.0) + Delta_prime((r[i]/Rmax)) * E((r[i]/Rmax)) * cos(theta) + Delta_prime((r[i]/Rmax)) * P((r[i]/Rmax)) * cos(theta) + 2.0 * Delta_prime((r[i]/Rmax)) * T((r[i]/Rmax)) * cos(2.0 * theta) - E((r[i]/Rmax)) * E_prime((r[i]/Rmax)) * pow(sin(theta), 2.0) - E((r[i]/Rmax)) * E_prime((r[i]/Rmax)) * pow(cos(theta), 2.0) - E((r[i]/Rmax)) * T_prime((r[i]/Rmax)) * sin(theta) * sin(2.0 * theta) + E((r[i]/Rmax)) * T_prime((r[i]/Rmax)) * cos(theta) * cos(2.0 * theta) - E((r[i]/Rmax)) * dP((r[i]/Rmax)) * pow(sin(theta), 2.0) - E((r[i]/Rmax)) * dP((r[i]/Rmax)) * pow(cos(theta), 2.0) + E((r[i]/Rmax)) * pow(sin(theta), 2.0) + E((r[i]/Rmax)) * pow(cos(theta), 2.0) - E_prime((r[i]/Rmax)) * P((r[i]/Rmax)) * pow(sin(theta), 2.0) - E_prime((r[i]/Rmax)) * P((r[i]/Rmax)) * pow(cos(theta), 2.0) + 2.0 * E_prime((r[i]/Rmax)) * T((r[i]/Rmax)) * sin(theta) * sin(2.0 * theta) - 2.0 * E_prime((r[i]/Rmax)) * T((r[i]/Rmax)) * cos(theta) * cos(2.0 * theta) - P((r[i]/Rmax)) * T_prime((r[i]/Rmax)) * sin(theta) * sin(2.0 * theta) + P((r[i]/Rmax)) * T_prime((r[i]/Rmax)) * cos(theta) * cos(2.0 * theta) - P((r[i]/Rmax)) * dP((r[i]/Rmax)) * pow(sin(theta), 2.0) - P((r[i]/Rmax)) * dP((r[i]/Rmax)) * pow(cos(theta), 2.0) + P((r[i]/Rmax)) * pow(sin(theta), 2.0) + P((r[i]/Rmax)) * pow(cos(theta), 2.0) + 2.0 * T((r[i]/Rmax)) * T_prime((r[i]/Rmax)) * pow(sin(2.0 * theta), 2.0) + 2.0 * T((r[i]/Rmax)) * T_prime((r[i]/Rmax)) * pow(cos(2.0 * theta), 2.0) + 2.0 * T((r[i]/Rmax)) * dP((r[i]/Rmax)) * sin(theta) * sin(2.0 * theta) - 2.0 * T((r[i]/Rmax)) * dP((r[i]/Rmax)) * cos(theta) * cos(2.0 * theta) - 2.0 * T((r[i]/Rmax)) * sin(theta) * sin(2.0 * theta) + 2.0 * T((r[i]/Rmax)) * cos(theta) * cos(2.0 * theta));
    }
}
/*........................................*/
void RefinedGyroZoniShiftedCulham::J_xq(double r, std::vector<double> const& theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = ((-(r/Rmax)) * sin_theta[i] + E((r/Rmax)) * sin_theta[i] + P((r/Rmax)) * sin_theta[i] - 2.0 * T((r/Rmax)) * sin(2.0 * theta[i])) / ((-(r/Rmax)) * Delta_prime((r/Rmax)) * cos_theta[i] + (r/Rmax) * E_prime((r/Rmax)) * pow(sin_theta[i], 2.0) + (r/Rmax) * E_prime((r/Rmax)) * pow(cos_theta[i], 2.0) + (r/Rmax) * T_prime((r/Rmax)) * sin_theta[i] * sin(2.0 * theta[i]) - (r/Rmax) * T_prime((r/Rmax)) * cos_theta[i] * cos(2.0 * theta[i]) + (r/Rmax) * dP((r/Rmax)) * pow(sin_theta[i], 2.0) + (r/Rmax) * dP((r/Rmax)) * pow(cos_theta[i], 2.0) - (r/Rmax) * pow(sin_theta[i], 2.0) - (r/Rmax) * pow(cos_theta[i], 2.0) + Delta_prime((r/Rmax)) * E((r/Rmax)) * cos_theta[i] + Delta_prime((r/Rmax)) * P((r/Rmax)) * cos_theta[i] + 2.0 * Delta_prime((r/Rmax)) * T((r/Rmax)) * cos(2.0 * theta[i]) - E((r/Rmax)) * E_prime((r/Rmax)) * pow(sin_theta[i], 2.0) - E((r/Rmax)) * E_prime((r/Rmax)) * pow(cos_theta[i], 2.0) - E((r/Rmax)) * T_prime((r/Rmax)) * sin_theta[i] * sin(2.0 * theta[i]) + E((r/Rmax)) * T_prime((r/Rmax)) * cos_theta[i] * cos(2.0 * theta[i]) - E((r/Rmax)) * dP((r/Rmax)) * pow(sin_theta[i], 2.0) - E((r/Rmax)) * dP((r/Rmax)) * pow(cos_theta[i], 2.0) + E((r/Rmax)) * pow(sin_theta[i], 2.0) + E((r/Rmax)) * pow(cos_theta[i], 2.0) - E_prime((r/Rmax)) * P((r/Rmax)) * pow(sin_theta[i], 2.0) - E_prime((r/Rmax)) * P((r/Rmax)) * pow(cos_theta[i], 2.0) + 2.0 * E_prime((r/Rmax)) * T((r/Rmax)) * sin_theta[i] * sin(2.0 * theta[i]) - 2.0 * E_prime((r/Rmax)) * T((r/Rmax)) * cos_theta[i] * cos(2.0 * theta[i]) - P((r/Rmax)) * T_prime((r/Rmax)) * sin_theta[i] * sin(2.0 * theta[i]) + P((r/Rmax)) * T_prime((r/Rmax)) * cos_theta[i] * cos(2.0 * theta[i]) - P((r/Rmax)) * dP((r/Rmax)) * pow(sin_theta[i], 2.0) - P((r/Rmax)) * dP((r/Rmax)) * pow(cos_theta[i], 2.0) + P((r/Rmax)) * pow(sin_theta[i], 2.0) + P((r/Rmax)) * pow(cos_theta[i], 2.0) + 2.0 * T((r/Rmax)) * T_prime((r/Rmax)) * pow(sin(2.0 * theta[i]), 2.0) + 2.0 * T((r/Rmax)) * T_prime((r/Rmax)) * pow(cos(2.0 * theta[i]), 2.0) + 2.0 * T((r/Rmax)) * dP((r/Rmax)) * sin_theta[i] * sin(2.0 * theta[i]) - 2.0 * T((r/Rmax)) * dP((r/Rmax)) * cos_theta[i] * cos(2.0 * theta[i]) - 2.0 * T((r/Rmax)) * sin_theta[i] * sin(2.0 * theta[i]) + 2.0 * T((r/Rmax)) * cos_theta[i] * cos(2.0 * theta[i]));
    }
}
/*........................................*/
double RefinedGyroZoniShiftedCulham::J_yr(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const
{
    return ((-E_prime((r/Rmax))) * sin(theta) - T_prime((r/Rmax)) * sin(2.0 * theta) - dP((r/Rmax)) * sin(theta) + sin(theta)) / ((-(r/Rmax)) * Delta_prime((r/Rmax)) * cos(theta) + (r/Rmax) * E_prime((r/Rmax)) * pow(sin(theta), 2.0) + (r/Rmax) * E_prime((r/Rmax)) * pow(cos(theta), 2.0) + (r/Rmax) * T_prime((r/Rmax)) * sin(theta) * sin(2.0 * theta) - (r/Rmax) * T_prime((r/Rmax)) * cos(theta) * cos(2.0 * theta) + (r/Rmax) * dP((r/Rmax)) * pow(sin(theta), 2.0) + (r/Rmax) * dP((r/Rmax)) * pow(cos(theta), 2.0) - (r/Rmax) * pow(sin(theta), 2.0) - (r/Rmax) * pow(cos(theta), 2.0) + Delta_prime((r/Rmax)) * E((r/Rmax)) * cos(theta) + Delta_prime((r/Rmax)) * P((r/Rmax)) * cos(theta) + 2.0 * Delta_prime((r/Rmax)) * T((r/Rmax)) * cos(2.0 * theta) - E((r/Rmax)) * E_prime((r/Rmax)) * pow(sin(theta), 2.0) - E((r/Rmax)) * E_prime((r/Rmax)) * pow(cos(theta), 2.0) - E((r/Rmax)) * T_prime((r/Rmax)) * sin(theta) * sin(2.0 * theta) + E((r/Rmax)) * T_prime((r/Rmax)) * cos(theta) * cos(2.0 * theta) - E((r/Rmax)) * dP((r/Rmax)) * pow(sin(theta), 2.0) - E((r/Rmax)) * dP((r/Rmax)) * pow(cos(theta), 2.0) + E((r/Rmax)) * pow(sin(theta), 2.0) + E((r/Rmax)) * pow(cos(theta), 2.0) - E_prime((r/Rmax)) * P((r/Rmax)) * pow(sin(theta), 2.0) - E_prime((r/Rmax)) * P((r/Rmax)) * pow(cos(theta), 2.0) + 2.0 * E_prime((r/Rmax)) * T((r/Rmax)) * sin(theta) * sin(2.0 * theta) - 2.0 * E_prime((r/Rmax)) * T((r/Rmax)) * cos(theta) * cos(2.0 * theta) - P((r/Rmax)) * T_prime((r/Rmax)) * sin(theta) * sin(2.0 * theta) + P((r/Rmax)) * T_prime((r/Rmax)) * cos(theta) * cos(2.0 * theta) - P((r/Rmax)) * dP((r/Rmax)) * pow(sin(theta), 2.0) - P((r/Rmax)) * dP((r/Rmax)) * pow(cos(theta), 2.0) + P((r/Rmax)) * pow(sin(theta), 2.0) + P((r/Rmax)) * pow(cos(theta), 2.0) + 2.0 * T((r/Rmax)) * T_prime((r/Rmax)) * pow(sin(2.0 * theta), 2.0) + 2.0 * T((r/Rmax)) * T_prime((r/Rmax)) * pow(cos(2.0 * theta), 2.0) + 2.0 * T((r/Rmax)) * dP((r/Rmax)) * sin(theta) * sin(2.0 * theta) - 2.0 * T((r/Rmax)) * dP((r/Rmax)) * cos(theta) * cos(2.0 * theta) - 2.0 * T((r/Rmax)) * sin(theta) * sin(2.0 * theta) + 2.0 * T((r/Rmax)) * cos(theta) * cos(2.0 * theta));
}
/*........................................*/
void RefinedGyroZoniShiftedCulham::J_yr(std::vector<double> const& r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = ((-E_prime((r[i]/Rmax))) * sin(theta) - T_prime((r[i]/Rmax)) * sin(2.0 * theta) - dP((r[i]/Rmax)) * sin(theta) + sin(theta)) / ((-(r[i]/Rmax)) * Delta_prime((r[i]/Rmax)) * cos(theta) + (r[i]/Rmax) * E_prime((r[i]/Rmax)) * pow(sin(theta), 2.0) + (r[i]/Rmax) * E_prime((r[i]/Rmax)) * pow(cos(theta), 2.0) + (r[i]/Rmax) * T_prime((r[i]/Rmax)) * sin(theta) * sin(2.0 * theta) - (r[i]/Rmax) * T_prime((r[i]/Rmax)) * cos(theta) * cos(2.0 * theta) + (r[i]/Rmax) * dP((r[i]/Rmax)) * pow(sin(theta), 2.0) + (r[i]/Rmax) * dP((r[i]/Rmax)) * pow(cos(theta), 2.0) - (r[i]/Rmax) * pow(sin(theta), 2.0) - (r[i]/Rmax) * pow(cos(theta), 2.0) + Delta_prime((r[i]/Rmax)) * E((r[i]/Rmax)) * cos(theta) + Delta_prime((r[i]/Rmax)) * P((r[i]/Rmax)) * cos(theta) + 2.0 * Delta_prime((r[i]/Rmax)) * T((r[i]/Rmax)) * cos(2.0 * theta) - E((r[i]/Rmax)) * E_prime((r[i]/Rmax)) * pow(sin(theta), 2.0) - E((r[i]/Rmax)) * E_prime((r[i]/Rmax)) * pow(cos(theta), 2.0) - E((r[i]/Rmax)) * T_prime((r[i]/Rmax)) * sin(theta) * sin(2.0 * theta) + E((r[i]/Rmax)) * T_prime((r[i]/Rmax)) * cos(theta) * cos(2.0 * theta) - E((r[i]/Rmax)) * dP((r[i]/Rmax)) * pow(sin(theta), 2.0) - E((r[i]/Rmax)) * dP((r[i]/Rmax)) * pow(cos(theta), 2.0) + E((r[i]/Rmax)) * pow(sin(theta), 2.0) + E((r[i]/Rmax)) * pow(cos(theta), 2.0) - E_prime((r[i]/Rmax)) * P((r[i]/Rmax)) * pow(sin(theta), 2.0) - E_prime((r[i]/Rmax)) * P((r[i]/Rmax)) * pow(cos(theta), 2.0) + 2.0 * E_prime((r[i]/Rmax)) * T((r[i]/Rmax)) * sin(theta) * sin(2.0 * theta) - 2.0 * E_prime((r[i]/Rmax)) * T((r[i]/Rmax)) * cos(theta) * cos(2.0 * theta) - P((r[i]/Rmax)) * T_prime((r[i]/Rmax)) * sin(theta) * sin(2.0 * theta) + P((r[i]/Rmax)) * T_prime((r[i]/Rmax)) * cos(theta) * cos(2.0 * theta) - P((r[i]/Rmax)) * dP((r[i]/Rmax)) * pow(sin(theta), 2.0) - P((r[i]/Rmax)) * dP((r[i]/Rmax)) * pow(cos(theta), 2.0) + P((r[i]/Rmax)) * pow(sin(theta), 2.0) + P((r[i]/Rmax)) * pow(cos(theta), 2.0) + 2.0 * T((r[i]/Rmax)) * T_prime((r[i]/Rmax)) * pow(sin(2.0 * theta), 2.0) + 2.0 * T((r[i]/Rmax)) * T_prime((r[i]/Rmax)) * pow(cos(2.0 * theta), 2.0) + 2.0 * T((r[i]/Rmax)) * dP((r[i]/Rmax)) * sin(theta) * sin(2.0 * theta) - 2.0 * T((r[i]/Rmax)) * dP((r[i]/Rmax)) * cos(theta) * cos(2.0 * theta) - 2.0 * T((r[i]/Rmax)) * sin(theta) * sin(2.0 * theta) + 2.0 * T((r[i]/Rmax)) * cos(theta) * cos(2.0 * theta));
    }
}
/*........................................*/
void RefinedGyroZoniShiftedCulham::J_yr(double r, std::vector<double> const& theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = ((-E_prime((r/Rmax))) * sin_theta[i] - T_prime((r/Rmax)) * sin(2.0 * theta[i]) - dP((r/Rmax)) * sin_theta[i] + sin_theta[i]) / ((-(r/Rmax)) * Delta_prime((r/Rmax)) * cos_theta[i] + (r/Rmax) * E_prime((r/Rmax)) * pow(sin_theta[i], 2.0) + (r/Rmax) * E_prime((r/Rmax)) * pow(cos_theta[i], 2.0) + (r/Rmax) * T_prime((r/Rmax)) * sin_theta[i] * sin(2.0 * theta[i]) - (r/Rmax) * T_prime((r/Rmax)) * cos_theta[i] * cos(2.0 * theta[i]) + (r/Rmax) * dP((r/Rmax)) * pow(sin_theta[i], 2.0) + (r/Rmax) * dP((r/Rmax)) * pow(cos_theta[i], 2.0) - (r/Rmax) * pow(sin_theta[i], 2.0) - (r/Rmax) * pow(cos_theta[i], 2.0) + Delta_prime((r/Rmax)) * E((r/Rmax)) * cos_theta[i] + Delta_prime((r/Rmax)) * P((r/Rmax)) * cos_theta[i] + 2.0 * Delta_prime((r/Rmax)) * T((r/Rmax)) * cos(2.0 * theta[i]) - E((r/Rmax)) * E_prime((r/Rmax)) * pow(sin_theta[i], 2.0) - E((r/Rmax)) * E_prime((r/Rmax)) * pow(cos_theta[i], 2.0) - E((r/Rmax)) * T_prime((r/Rmax)) * sin_theta[i] * sin(2.0 * theta[i]) + E((r/Rmax)) * T_prime((r/Rmax)) * cos_theta[i] * cos(2.0 * theta[i]) - E((r/Rmax)) * dP((r/Rmax)) * pow(sin_theta[i], 2.0) - E((r/Rmax)) * dP((r/Rmax)) * pow(cos_theta[i], 2.0) + E((r/Rmax)) * pow(sin_theta[i], 2.0) + E((r/Rmax)) * pow(cos_theta[i], 2.0) - E_prime((r/Rmax)) * P((r/Rmax)) * pow(sin_theta[i], 2.0) - E_prime((r/Rmax)) * P((r/Rmax)) * pow(cos_theta[i], 2.0) + 2.0 * E_prime((r/Rmax)) * T((r/Rmax)) * sin_theta[i] * sin(2.0 * theta[i]) - 2.0 * E_prime((r/Rmax)) * T((r/Rmax)) * cos_theta[i] * cos(2.0 * theta[i]) - P((r/Rmax)) * T_prime((r/Rmax)) * sin_theta[i] * sin(2.0 * theta[i]) + P((r/Rmax)) * T_prime((r/Rmax)) * cos_theta[i] * cos(2.0 * theta[i]) - P((r/Rmax)) * dP((r/Rmax)) * pow(sin_theta[i], 2.0) - P((r/Rmax)) * dP((r/Rmax)) * pow(cos_theta[i], 2.0) + P((r/Rmax)) * pow(sin_theta[i], 2.0) + P((r/Rmax)) * pow(cos_theta[i], 2.0) + 2.0 * T((r/Rmax)) * T_prime((r/Rmax)) * pow(sin(2.0 * theta[i]), 2.0) + 2.0 * T((r/Rmax)) * T_prime((r/Rmax)) * pow(cos(2.0 * theta[i]), 2.0) + 2.0 * T((r/Rmax)) * dP((r/Rmax)) * sin_theta[i] * sin(2.0 * theta[i]) - 2.0 * T((r/Rmax)) * dP((r/Rmax)) * cos_theta[i] * cos(2.0 * theta[i]) - 2.0 * T((r/Rmax)) * sin_theta[i] * sin(2.0 * theta[i]) + 2.0 * T((r/Rmax)) * cos_theta[i] * cos(2.0 * theta[i]));
    }
}
/*........................................*/
double RefinedGyroZoniShiftedCulham::J_yq(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const
{
    return ((-Delta_prime((r/Rmax))) + E_prime((r/Rmax)) * cos(theta) - T_prime((r/Rmax)) * cos(2.0 * theta) + dP((r/Rmax)) * cos(theta) - cos(theta)) / ((-(r/Rmax)) * Delta_prime((r/Rmax)) * cos(theta) + (r/Rmax) * E_prime((r/Rmax)) * pow(sin(theta), 2.0) + (r/Rmax) * E_prime((r/Rmax)) * pow(cos(theta), 2.0) + (r/Rmax) * T_prime((r/Rmax)) * sin(theta) * sin(2.0 * theta) - (r/Rmax) * T_prime((r/Rmax)) * cos(theta) * cos(2.0 * theta) + (r/Rmax) * dP((r/Rmax)) * pow(sin(theta), 2.0) + (r/Rmax) * dP((r/Rmax)) * pow(cos(theta), 2.0) - (r/Rmax) * pow(sin(theta), 2.0) - (r/Rmax) * pow(cos(theta), 2.0) + Delta_prime((r/Rmax)) * E((r/Rmax)) * cos(theta) + Delta_prime((r/Rmax)) * P((r/Rmax)) * cos(theta) + 2.0 * Delta_prime((r/Rmax)) * T((r/Rmax)) * cos(2.0 * theta) - E((r/Rmax)) * E_prime((r/Rmax)) * pow(sin(theta), 2.0) - E((r/Rmax)) * E_prime((r/Rmax)) * pow(cos(theta), 2.0) - E((r/Rmax)) * T_prime((r/Rmax)) * sin(theta) * sin(2.0 * theta) + E((r/Rmax)) * T_prime((r/Rmax)) * cos(theta) * cos(2.0 * theta) - E((r/Rmax)) * dP((r/Rmax)) * pow(sin(theta), 2.0) - E((r/Rmax)) * dP((r/Rmax)) * pow(cos(theta), 2.0) + E((r/Rmax)) * pow(sin(theta), 2.0) + E((r/Rmax)) * pow(cos(theta), 2.0) - E_prime((r/Rmax)) * P((r/Rmax)) * pow(sin(theta), 2.0) - E_prime((r/Rmax)) * P((r/Rmax)) * pow(cos(theta), 2.0) + 2.0 * E_prime((r/Rmax)) * T((r/Rmax)) * sin(theta) * sin(2.0 * theta) - 2.0 * E_prime((r/Rmax)) * T((r/Rmax)) * cos(theta) * cos(2.0 * theta) - P((r/Rmax)) * T_prime((r/Rmax)) * sin(theta) * sin(2.0 * theta) + P((r/Rmax)) * T_prime((r/Rmax)) * cos(theta) * cos(2.0 * theta) - P((r/Rmax)) * dP((r/Rmax)) * pow(sin(theta), 2.0) - P((r/Rmax)) * dP((r/Rmax)) * pow(cos(theta), 2.0) + P((r/Rmax)) * pow(sin(theta), 2.0) + P((r/Rmax)) * pow(cos(theta), 2.0) + 2.0 * T((r/Rmax)) * T_prime((r/Rmax)) * pow(sin(2.0 * theta), 2.0) + 2.0 * T((r/Rmax)) * T_prime((r/Rmax)) * pow(cos(2.0 * theta), 2.0) + 2.0 * T((r/Rmax)) * dP((r/Rmax)) * sin(theta) * sin(2.0 * theta) - 2.0 * T((r/Rmax)) * dP((r/Rmax)) * cos(theta) * cos(2.0 * theta) - 2.0 * T((r/Rmax)) * sin(theta) * sin(2.0 * theta) + 2.0 * T((r/Rmax)) * cos(theta) * cos(2.0 * theta));
}
/*........................................*/
void RefinedGyroZoniShiftedCulham::J_yq(std::vector<double> const& r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = ((-Delta_prime((r[i]/Rmax))) + E_prime((r[i]/Rmax)) * cos(theta) - T_prime((r[i]/Rmax)) * cos(2.0 * theta) + dP((r[i]/Rmax)) * cos(theta) - cos(theta)) / ((-(r[i]/Rmax)) * Delta_prime((r[i]/Rmax)) * cos(theta) + (r[i]/Rmax) * E_prime((r[i]/Rmax)) * pow(sin(theta), 2.0) + (r[i]/Rmax) * E_prime((r[i]/Rmax)) * pow(cos(theta), 2.0) + (r[i]/Rmax) * T_prime((r[i]/Rmax)) * sin(theta) * sin(2.0 * theta) - (r[i]/Rmax) * T_prime((r[i]/Rmax)) * cos(theta) * cos(2.0 * theta) + (r[i]/Rmax) * dP((r[i]/Rmax)) * pow(sin(theta), 2.0) + (r[i]/Rmax) * dP((r[i]/Rmax)) * pow(cos(theta), 2.0) - (r[i]/Rmax) * pow(sin(theta), 2.0) - (r[i]/Rmax) * pow(cos(theta), 2.0) + Delta_prime((r[i]/Rmax)) * E((r[i]/Rmax)) * cos(theta) + Delta_prime((r[i]/Rmax)) * P((r[i]/Rmax)) * cos(theta) + 2.0 * Delta_prime((r[i]/Rmax)) * T((r[i]/Rmax)) * cos(2.0 * theta) - E((r[i]/Rmax)) * E_prime((r[i]/Rmax)) * pow(sin(theta), 2.0) - E((r[i]/Rmax)) * E_prime((r[i]/Rmax)) * pow(cos(theta), 2.0) - E((r[i]/Rmax)) * T_prime((r[i]/Rmax)) * sin(theta) * sin(2.0 * theta) + E((r[i]/Rmax)) * T_prime((r[i]/Rmax)) * cos(theta) * cos(2.0 * theta) - E((r[i]/Rmax)) * dP((r[i]/Rmax)) * pow(sin(theta), 2.0) - E((r[i]/Rmax)) * dP((r[i]/Rmax)) * pow(cos(theta), 2.0) + E((r[i]/Rmax)) * pow(sin(theta), 2.0) + E((r[i]/Rmax)) * pow(cos(theta), 2.0) - E_prime((r[i]/Rmax)) * P((r[i]/Rmax)) * pow(sin(theta), 2.0) - E_prime((r[i]/Rmax)) * P((r[i]/Rmax)) * pow(cos(theta), 2.0) + 2.0 * E_prime((r[i]/Rmax)) * T((r[i]/Rmax)) * sin(theta) * sin(2.0 * theta) - 2.0 * E_prime((r[i]/Rmax)) * T((r[i]/Rmax)) * cos(theta) * cos(2.0 * theta) - P((r[i]/Rmax)) * T_prime((r[i]/Rmax)) * sin(theta) * sin(2.0 * theta) + P((r[i]/Rmax)) * T_prime((r[i]/Rmax)) * cos(theta) * cos(2.0 * theta) - P((r[i]/Rmax)) * dP((r[i]/Rmax)) * pow(sin(theta), 2.0) - P((r[i]/Rmax)) * dP((r[i]/Rmax)) * pow(cos(theta), 2.0) + P((r[i]/Rmax)) * pow(sin(theta), 2.0) + P((r[i]/Rmax)) * pow(cos(theta), 2.0) + 2.0 * T((r[i]/Rmax)) * T_prime((r[i]/Rmax)) * pow(sin(2.0 * theta), 2.0) + 2.0 * T((r[i]/Rmax)) * T_prime((r[i]/Rmax)) * pow(cos(2.0 * theta), 2.0) + 2.0 * T((r[i]/Rmax)) * dP((r[i]/Rmax)) * sin(theta) * sin(2.0 * theta) - 2.0 * T((r[i]/Rmax)) * dP((r[i]/Rmax)) * cos(theta) * cos(2.0 * theta) - 2.0 * T((r[i]/Rmax)) * sin(theta) * sin(2.0 * theta) + 2.0 * T((r[i]/Rmax)) * cos(theta) * cos(2.0 * theta));
    }
}
/*........................................*/
void RefinedGyroZoniShiftedCulham::J_yq(double r, std::vector<double> const& theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = ((-Delta_prime((r/Rmax))) + E_prime((r/Rmax)) * cos_theta[i] - T_prime((r/Rmax)) * cos(2.0 * theta[i]) + dP((r/Rmax)) * cos_theta[i] - cos_theta[i]) / ((-(r/Rmax)) * Delta_prime((r/Rmax)) * cos_theta[i] + (r/Rmax) * E_prime((r/Rmax)) * pow(sin_theta[i], 2.0) + (r/Rmax) * E_prime((r/Rmax)) * pow(cos_theta[i], 2.0) + (r/Rmax) * T_prime((r/Rmax)) * sin_theta[i] * sin(2.0 * theta[i]) - (r/Rmax) * T_prime((r/Rmax)) * cos_theta[i] * cos(2.0 * theta[i]) + (r/Rmax) * dP((r/Rmax)) * pow(sin_theta[i], 2.0) + (r/Rmax) * dP((r/Rmax)) * pow(cos_theta[i], 2.0) - (r/Rmax) * pow(sin_theta[i], 2.0) - (r/Rmax) * pow(cos_theta[i], 2.0) + Delta_prime((r/Rmax)) * E((r/Rmax)) * cos_theta[i] + Delta_prime((r/Rmax)) * P((r/Rmax)) * cos_theta[i] + 2.0 * Delta_prime((r/Rmax)) * T((r/Rmax)) * cos(2.0 * theta[i]) - E((r/Rmax)) * E_prime((r/Rmax)) * pow(sin_theta[i], 2.0) - E((r/Rmax)) * E_prime((r/Rmax)) * pow(cos_theta[i], 2.0) - E((r/Rmax)) * T_prime((r/Rmax)) * sin_theta[i] * sin(2.0 * theta[i]) + E((r/Rmax)) * T_prime((r/Rmax)) * cos_theta[i] * cos(2.0 * theta[i]) - E((r/Rmax)) * dP((r/Rmax)) * pow(sin_theta[i], 2.0) - E((r/Rmax)) * dP((r/Rmax)) * pow(cos_theta[i], 2.0) + E((r/Rmax)) * pow(sin_theta[i], 2.0) + E((r/Rmax)) * pow(cos_theta[i], 2.0) - E_prime((r/Rmax)) * P((r/Rmax)) * pow(sin_theta[i], 2.0) - E_prime((r/Rmax)) * P((r/Rmax)) * pow(cos_theta[i], 2.0) + 2.0 * E_prime((r/Rmax)) * T((r/Rmax)) * sin_theta[i] * sin(2.0 * theta[i]) - 2.0 * E_prime((r/Rmax)) * T((r/Rmax)) * cos_theta[i] * cos(2.0 * theta[i]) - P((r/Rmax)) * T_prime((r/Rmax)) * sin_theta[i] * sin(2.0 * theta[i]) + P((r/Rmax)) * T_prime((r/Rmax)) * cos_theta[i] * cos(2.0 * theta[i]) - P((r/Rmax)) * dP((r/Rmax)) * pow(sin_theta[i], 2.0) - P((r/Rmax)) * dP((r/Rmax)) * pow(cos_theta[i], 2.0) + P((r/Rmax)) * pow(sin_theta[i], 2.0) + P((r/Rmax)) * pow(cos_theta[i], 2.0) + 2.0 * T((r/Rmax)) * T_prime((r/Rmax)) * pow(sin(2.0 * theta[i]), 2.0) + 2.0 * T((r/Rmax)) * T_prime((r/Rmax)) * pow(cos(2.0 * theta[i]), 2.0) + 2.0 * T((r/Rmax)) * dP((r/Rmax)) * sin_theta[i] * sin(2.0 * theta[i]) - 2.0 * T((r/Rmax)) * dP((r/Rmax)) * cos_theta[i] * cos(2.0 * theta[i]) - 2.0 * T((r/Rmax)) * sin_theta[i] * sin(2.0 * theta[i]) + 2.0 * T((r/Rmax)) * cos_theta[i] * cos(2.0 * theta[i]));
    }
}
/*........................................*/
double RefinedGyroZoniShiftedCulham::coeffs1(double r, double Rmax) const
{
    return exp(-tanh(20.0 * (r/Rmax) - 14.0));
}
/*........................................*/
void RefinedGyroZoniShiftedCulham::coeffs1(std::vector<double> const& r, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = exp(-tanh(20.0 * (r[i]/Rmax) - 14.0));
    }
}
/*........................................*/
double RefinedGyroZoniShiftedCulham::coeffs2(double r, double Rmax) const
{
    return 1.0 * exp(tanh(20.0 * (r/Rmax) - 14.0));
}
/*........................................*/
void RefinedGyroZoniShiftedCulham::coeffs2(std::vector<double> const& r, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = 1.0 * exp(tanh(20.0 * (r[i]/Rmax) - 14.0));
    }
}
/*........................................*/
double RefinedGyroZoniShiftedCulham::rho_glob(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const
{
    return ((-3.33823779536505e-15) * ((r/Rmax) * (r/Rmax)) - 0.0 * (r/Rmax) - 0.0 + exp((-3333.33333333333) * pow(((r/Rmax) - 0.9), 2.0))) * cos(21.0 * theta) + (0.00184273372222541 * ((r/Rmax) * (r/Rmax)) - 0.0018029383826828 * (r/Rmax) - 4.00652973929511e-05 + exp((-50.0) * pow(((r/Rmax) - 0.45), 2.0))) * cos(9.0 * theta);
}
/*........................................*/
void RefinedGyroZoniShiftedCulham::rho_glob(std::vector<double> const& r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = ((-3.33823779536505e-15) * ((r[i]/Rmax) * (r[i]/Rmax)) - 0.0 * (r[i]/Rmax) - 0.0 + exp((-3333.33333333333) * pow(((r[i]/Rmax) - 0.9), 2.0))) * cos(21.0 * theta) + (0.00184273372222541 * ((r[i]/Rmax) * (r[i]/Rmax)) - 0.0018029383826828 * (r[i]/Rmax) - 4.00652973929511e-05 + exp((-50.0) * pow(((r[i]/Rmax) - 0.45), 2.0))) * cos(9.0 * theta);
    }
}
/*........................................*/
void RefinedGyroZoniShiftedCulham::rho_glob(double r, std::vector<double> const& theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = ((-3.33823779536505e-15) * ((r/Rmax) * (r/Rmax)) - 0.0 * (r/Rmax) - 0.0 + exp((-3333.33333333333) * pow(((r/Rmax) - 0.9), 2.0))) * cos(21.0 * theta[i]) + (0.00184273372222541 * ((r/Rmax) * (r/Rmax)) - 0.0018029383826828 * (r/Rmax) - 4.00652973929511e-05 + exp((-50.0) * pow(((r/Rmax) - 0.45), 2.0))) * cos(9.0 * theta[i]);
    }
}
/*........................................*/
double RefinedGyroZoniShiftedCulham::rho_pole(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const
{
    return 0.0;
}
/*........................................*/
void RefinedGyroZoniShiftedCulham::rho_pole(std::vector<double> const& r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = 0.0;
    }
}
/*........................................*/
void RefinedGyroZoniShiftedCulham::rho_pole(double r, std::vector<double> const& theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = 0.0;
    }
}
/*........................................*/
double RefinedGyroZoniShiftedCulham::phi_exact(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const
{
    return 0.0;
}
/*........................................*/
void RefinedGyroZoniShiftedCulham::phi_exact(std::vector<double> const& r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = 0.0;
    }
}
/*........................................*/
void RefinedGyroZoniShiftedCulham::phi_exact(double r, std::vector<double> const& theta, double map3_unused_1, double map3_unused_2, double Rmax, std::vector<double>& sol, std::vector<double>& sin_theta, std::vector<double>& cos_theta) const
{
    for (std::size_t i=0; i < sol.size(); ++i)
    {
        sol[i] = 0.0;
    }
}
/*........................................*/
