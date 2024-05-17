// PolarR6 simulates solution (22) of Bourne et al. https://doi.org/10.1016/j.jcp.2023.112249
#include "../../include/test_cases/PolarR6GyroZoniShiftedCulham.h"

inline double PolarR6GyroZoniShiftedCulham::q(double rr) const
{
    return 0.8 - 0.1 * (rr * rr);
}

inline double PolarR6GyroZoniShiftedCulham::dq(double rr) const
{
    return (-0.2) * rr;
}

inline double PolarR6GyroZoniShiftedCulham::p(double rr) const
{
    return 100000.0 - 90000.0 * (rr * rr);
}

inline double PolarR6GyroZoniShiftedCulham::dp(double rr) const
{
    return (-180000.0) * rr;
}

inline double PolarR6GyroZoniShiftedCulham::dg(double rr, double g) const
{
    return ((-g) * (0.0625000000000001 * (rr * rr) / pow((1.0 - 0.125 * (rr * rr)), 2.0) + 2.0 / (4.0 - 0.5 * (rr * rr))) + 2.261946711816e-06 * (4.0 - 0.5 * (rr * rr)) / g) / (rr / (4.0 - 0.5 * (rr * rr)) + (4.0 - 0.5 * (rr * rr)) / (g * rr));
}

inline double PolarR6GyroZoniShiftedCulham::double_deriv(double rr, double c, double g, double dg, double val, double d_val) const
{
    return c * val / (rr * rr) - d_val * (pow(rr, (double)((-1))) + (4.0 - 0.5 * (rr * rr)) * (2.0 * dg * rr / (4.0 - 0.5 * (rr * rr)) + 0.125 * g * (rr * rr) / pow((1.0 - 0.125 * (rr * rr)), 2.0) + 2.0 * g / (4.0 - 0.5 * (rr * rr))) / (g * rr));
}

inline double PolarR6GyroZoniShiftedCulham::g(double rr) const
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

inline double PolarR6GyroZoniShiftedCulham::Delta(double rr) const
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

inline double PolarR6GyroZoniShiftedCulham::Delta_prime(double rr) const
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

inline double PolarR6GyroZoniShiftedCulham::E(double rr) const
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

inline double PolarR6GyroZoniShiftedCulham::T(double rr) const
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

inline double PolarR6GyroZoniShiftedCulham::E_prime(double rr) const
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

inline double PolarR6GyroZoniShiftedCulham::T_prime(double rr) const
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

inline double PolarR6GyroZoniShiftedCulham::P(double rr) const
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

inline double PolarR6GyroZoniShiftedCulham::dP(double rr) const
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

/* ---------------------------------------------------------------------------------------------------------------------------------*/

inline double PolarR6GyroZoniShiftedCulham::x(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const
{
    return (r/Rmax) * cos(theta) + Delta((r/Rmax)) - E((r/Rmax)) * cos(theta) - P((r/Rmax)) * cos(theta) + T((r/Rmax)) * cos(2.0 * theta) + 5.0;
}

inline double PolarR6GyroZoniShiftedCulham::y(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const
{
    return (r/Rmax) * sin(theta) - E((r/Rmax)) * sin(theta) - P((r/Rmax)) * sin(theta) - T((r/Rmax)) * sin(2.0 * theta);
}

inline double PolarR6GyroZoniShiftedCulham::J_rr(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const
{
    return (Delta_prime((r/Rmax)) - E_prime((r/Rmax)) * cos(theta) + T_prime((r/Rmax)) * cos(2.0 * theta) - dP((r/Rmax)) * cos(theta) + cos(theta))/Rmax;
}

inline double PolarR6GyroZoniShiftedCulham::J_rt(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const
{
    return (-(r/Rmax)) * sin(theta) + E((r/Rmax)) * sin(theta) + P((r/Rmax)) * sin(theta) - 2.0 * T((r/Rmax)) * sin(2.0 * theta);
}

inline double PolarR6GyroZoniShiftedCulham::J_tr(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const
{
    return ((-E_prime((r/Rmax))) * sin(theta) - T_prime((r/Rmax)) * sin(2.0 * theta) - dP((r/Rmax)) * sin(theta) + sin(theta))/Rmax;
}

inline double PolarR6GyroZoniShiftedCulham::J_tt(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const
{
    return (r/Rmax) * cos(theta) - E((r/Rmax)) * cos(theta) - P((r/Rmax)) * cos(theta) - 2.0 * T((r/Rmax)) * cos(2.0 * theta);
}

// inline double PolarR6GyroZoniShiftedCulham::J_xr(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const
// {
//     return ((-(r/Rmax)) * cos(theta) + E((r/Rmax)) * cos(theta) + P((r/Rmax)) * cos(theta) + 2.0 * T((r/Rmax)) * cos(2.0 * theta)) / ((-(r/Rmax)) * Delta_prime((r/Rmax)) * cos(theta) + (r/Rmax) * E_prime((r/Rmax)) * pow(sin(theta), 2.0) + (r/Rmax) * E_prime((r/Rmax)) * pow(cos(theta), 2.0) + (r/Rmax) * T_prime((r/Rmax)) * sin(theta) * sin(2.0 * theta) - (r/Rmax) * T_prime((r/Rmax)) * cos(theta) * cos(2.0 * theta) + (r/Rmax) * dP((r/Rmax)) * pow(sin(theta), 2.0) + (r/Rmax) * dP((r/Rmax)) * pow(cos(theta), 2.0) - (r/Rmax) * pow(sin(theta), 2.0) - (r/Rmax) * pow(cos(theta), 2.0) + Delta_prime((r/Rmax)) * E((r/Rmax)) * cos(theta) + Delta_prime((r/Rmax)) * P((r/Rmax)) * cos(theta) + 2.0 * Delta_prime((r/Rmax)) * T((r/Rmax)) * cos(2.0 * theta) - E((r/Rmax)) * E_prime((r/Rmax)) * pow(sin(theta), 2.0) - E((r/Rmax)) * E_prime((r/Rmax)) * pow(cos(theta), 2.0) - E((r/Rmax)) * T_prime((r/Rmax)) * sin(theta) * sin(2.0 * theta) + E((r/Rmax)) * T_prime((r/Rmax)) * cos(theta) * cos(2.0 * theta) - E((r/Rmax)) * dP((r/Rmax)) * pow(sin(theta), 2.0) - E((r/Rmax)) * dP((r/Rmax)) * pow(cos(theta), 2.0) + E((r/Rmax)) * pow(sin(theta), 2.0) + E((r/Rmax)) * pow(cos(theta), 2.0) - E_prime((r/Rmax)) * P((r/Rmax)) * pow(sin(theta), 2.0) - E_prime((r/Rmax)) * P((r/Rmax)) * pow(cos(theta), 2.0) + 2.0 * E_prime((r/Rmax)) * T((r/Rmax)) * sin(theta) * sin(2.0 * theta) - 2.0 * E_prime((r/Rmax)) * T((r/Rmax)) * cos(theta) * cos(2.0 * theta) - P((r/Rmax)) * T_prime((r/Rmax)) * sin(theta) * sin(2.0 * theta) + P((r/Rmax)) * T_prime((r/Rmax)) * cos(theta) * cos(2.0 * theta) - P((r/Rmax)) * dP((r/Rmax)) * pow(sin(theta), 2.0) - P((r/Rmax)) * dP((r/Rmax)) * pow(cos(theta), 2.0) + P((r/Rmax)) * pow(sin(theta), 2.0) + P((r/Rmax)) * pow(cos(theta), 2.0) + 2.0 * T((r/Rmax)) * T_prime((r/Rmax)) * pow(sin(2.0 * theta), 2.0) + 2.0 * T((r/Rmax)) * T_prime((r/Rmax)) * pow(cos(2.0 * theta), 2.0) + 2.0 * T((r/Rmax)) * dP((r/Rmax)) * sin(theta) * sin(2.0 * theta) - 2.0 * T((r/Rmax)) * dP((r/Rmax)) * cos(theta) * cos(2.0 * theta) - 2.0 * T((r/Rmax)) * sin(theta) * sin(2.0 * theta) + 2.0 * T((r/Rmax)) * cos(theta) * cos(2.0 * theta));
// }

// inline double PolarR6GyroZoniShiftedCulham::J_xq(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const
// {
//     return ((-(r/Rmax)) * sin(theta) + E((r/Rmax)) * sin(theta) + P((r/Rmax)) * sin(theta) - 2.0 * T((r/Rmax)) * sin(2.0 * theta)) / ((-(r/Rmax)) * Delta_prime((r/Rmax)) * cos(theta) + (r/Rmax) * E_prime((r/Rmax)) * pow(sin(theta), 2.0) + (r/Rmax) * E_prime((r/Rmax)) * pow(cos(theta), 2.0) + (r/Rmax) * T_prime((r/Rmax)) * sin(theta) * sin(2.0 * theta) - (r/Rmax) * T_prime((r/Rmax)) * cos(theta) * cos(2.0 * theta) + (r/Rmax) * dP((r/Rmax)) * pow(sin(theta), 2.0) + (r/Rmax) * dP((r/Rmax)) * pow(cos(theta), 2.0) - (r/Rmax) * pow(sin(theta), 2.0) - (r/Rmax) * pow(cos(theta), 2.0) + Delta_prime((r/Rmax)) * E((r/Rmax)) * cos(theta) + Delta_prime((r/Rmax)) * P((r/Rmax)) * cos(theta) + 2.0 * Delta_prime((r/Rmax)) * T((r/Rmax)) * cos(2.0 * theta) - E((r/Rmax)) * E_prime((r/Rmax)) * pow(sin(theta), 2.0) - E((r/Rmax)) * E_prime((r/Rmax)) * pow(cos(theta), 2.0) - E((r/Rmax)) * T_prime((r/Rmax)) * sin(theta) * sin(2.0 * theta) + E((r/Rmax)) * T_prime((r/Rmax)) * cos(theta) * cos(2.0 * theta) - E((r/Rmax)) * dP((r/Rmax)) * pow(sin(theta), 2.0) - E((r/Rmax)) * dP((r/Rmax)) * pow(cos(theta), 2.0) + E((r/Rmax)) * pow(sin(theta), 2.0) + E((r/Rmax)) * pow(cos(theta), 2.0) - E_prime((r/Rmax)) * P((r/Rmax)) * pow(sin(theta), 2.0) - E_prime((r/Rmax)) * P((r/Rmax)) * pow(cos(theta), 2.0) + 2.0 * E_prime((r/Rmax)) * T((r/Rmax)) * sin(theta) * sin(2.0 * theta) - 2.0 * E_prime((r/Rmax)) * T((r/Rmax)) * cos(theta) * cos(2.0 * theta) - P((r/Rmax)) * T_prime((r/Rmax)) * sin(theta) * sin(2.0 * theta) + P((r/Rmax)) * T_prime((r/Rmax)) * cos(theta) * cos(2.0 * theta) - P((r/Rmax)) * dP((r/Rmax)) * pow(sin(theta), 2.0) - P((r/Rmax)) * dP((r/Rmax)) * pow(cos(theta), 2.0) + P((r/Rmax)) * pow(sin(theta), 2.0) + P((r/Rmax)) * pow(cos(theta), 2.0) + 2.0 * T((r/Rmax)) * T_prime((r/Rmax)) * pow(sin(2.0 * theta), 2.0) + 2.0 * T((r/Rmax)) * T_prime((r/Rmax)) * pow(cos(2.0 * theta), 2.0) + 2.0 * T((r/Rmax)) * dP((r/Rmax)) * sin(theta) * sin(2.0 * theta) - 2.0 * T((r/Rmax)) * dP((r/Rmax)) * cos(theta) * cos(2.0 * theta) - 2.0 * T((r/Rmax)) * sin(theta) * sin(2.0 * theta) + 2.0 * T((r/Rmax)) * cos(theta) * cos(2.0 * theta));
// }

// inline double PolarR6GyroZoniShiftedCulham::J_yr(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const
// {
//     return ((-E_prime((r/Rmax))) * sin(theta) - T_prime((r/Rmax)) * sin(2.0 * theta) - dP((r/Rmax)) * sin(theta) + sin(theta)) / ((-(r/Rmax)) * Delta_prime((r/Rmax)) * cos(theta) + (r/Rmax) * E_prime((r/Rmax)) * pow(sin(theta), 2.0) + (r/Rmax) * E_prime((r/Rmax)) * pow(cos(theta), 2.0) + (r/Rmax) * T_prime((r/Rmax)) * sin(theta) * sin(2.0 * theta) - (r/Rmax) * T_prime((r/Rmax)) * cos(theta) * cos(2.0 * theta) + (r/Rmax) * dP((r/Rmax)) * pow(sin(theta), 2.0) + (r/Rmax) * dP((r/Rmax)) * pow(cos(theta), 2.0) - (r/Rmax) * pow(sin(theta), 2.0) - (r/Rmax) * pow(cos(theta), 2.0) + Delta_prime((r/Rmax)) * E((r/Rmax)) * cos(theta) + Delta_prime((r/Rmax)) * P((r/Rmax)) * cos(theta) + 2.0 * Delta_prime((r/Rmax)) * T((r/Rmax)) * cos(2.0 * theta) - E((r/Rmax)) * E_prime((r/Rmax)) * pow(sin(theta), 2.0) - E((r/Rmax)) * E_prime((r/Rmax)) * pow(cos(theta), 2.0) - E((r/Rmax)) * T_prime((r/Rmax)) * sin(theta) * sin(2.0 * theta) + E((r/Rmax)) * T_prime((r/Rmax)) * cos(theta) * cos(2.0 * theta) - E((r/Rmax)) * dP((r/Rmax)) * pow(sin(theta), 2.0) - E((r/Rmax)) * dP((r/Rmax)) * pow(cos(theta), 2.0) + E((r/Rmax)) * pow(sin(theta), 2.0) + E((r/Rmax)) * pow(cos(theta), 2.0) - E_prime((r/Rmax)) * P((r/Rmax)) * pow(sin(theta), 2.0) - E_prime((r/Rmax)) * P((r/Rmax)) * pow(cos(theta), 2.0) + 2.0 * E_prime((r/Rmax)) * T((r/Rmax)) * sin(theta) * sin(2.0 * theta) - 2.0 * E_prime((r/Rmax)) * T((r/Rmax)) * cos(theta) * cos(2.0 * theta) - P((r/Rmax)) * T_prime((r/Rmax)) * sin(theta) * sin(2.0 * theta) + P((r/Rmax)) * T_prime((r/Rmax)) * cos(theta) * cos(2.0 * theta) - P((r/Rmax)) * dP((r/Rmax)) * pow(sin(theta), 2.0) - P((r/Rmax)) * dP((r/Rmax)) * pow(cos(theta), 2.0) + P((r/Rmax)) * pow(sin(theta), 2.0) + P((r/Rmax)) * pow(cos(theta), 2.0) + 2.0 * T((r/Rmax)) * T_prime((r/Rmax)) * pow(sin(2.0 * theta), 2.0) + 2.0 * T((r/Rmax)) * T_prime((r/Rmax)) * pow(cos(2.0 * theta), 2.0) + 2.0 * T((r/Rmax)) * dP((r/Rmax)) * sin(theta) * sin(2.0 * theta) - 2.0 * T((r/Rmax)) * dP((r/Rmax)) * cos(theta) * cos(2.0 * theta) - 2.0 * T((r/Rmax)) * sin(theta) * sin(2.0 * theta) + 2.0 * T((r/Rmax)) * cos(theta) * cos(2.0 * theta));
// }

// inline double PolarR6GyroZoniShiftedCulham::J_yq(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const
// {
//     return ((-Delta_prime((r/Rmax))) + E_prime((r/Rmax)) * cos(theta) - T_prime((r/Rmax)) * cos(2.0 * theta) + dP((r/Rmax)) * cos(theta) - cos(theta)) / ((-(r/Rmax)) * Delta_prime((r/Rmax)) * cos(theta) + (r/Rmax) * E_prime((r/Rmax)) * pow(sin(theta), 2.0) + (r/Rmax) * E_prime((r/Rmax)) * pow(cos(theta), 2.0) + (r/Rmax) * T_prime((r/Rmax)) * sin(theta) * sin(2.0 * theta) - (r/Rmax) * T_prime((r/Rmax)) * cos(theta) * cos(2.0 * theta) + (r/Rmax) * dP((r/Rmax)) * pow(sin(theta), 2.0) + (r/Rmax) * dP((r/Rmax)) * pow(cos(theta), 2.0) - (r/Rmax) * pow(sin(theta), 2.0) - (r/Rmax) * pow(cos(theta), 2.0) + Delta_prime((r/Rmax)) * E((r/Rmax)) * cos(theta) + Delta_prime((r/Rmax)) * P((r/Rmax)) * cos(theta) + 2.0 * Delta_prime((r/Rmax)) * T((r/Rmax)) * cos(2.0 * theta) - E((r/Rmax)) * E_prime((r/Rmax)) * pow(sin(theta), 2.0) - E((r/Rmax)) * E_prime((r/Rmax)) * pow(cos(theta), 2.0) - E((r/Rmax)) * T_prime((r/Rmax)) * sin(theta) * sin(2.0 * theta) + E((r/Rmax)) * T_prime((r/Rmax)) * cos(theta) * cos(2.0 * theta) - E((r/Rmax)) * dP((r/Rmax)) * pow(sin(theta), 2.0) - E((r/Rmax)) * dP((r/Rmax)) * pow(cos(theta), 2.0) + E((r/Rmax)) * pow(sin(theta), 2.0) + E((r/Rmax)) * pow(cos(theta), 2.0) - E_prime((r/Rmax)) * P((r/Rmax)) * pow(sin(theta), 2.0) - E_prime((r/Rmax)) * P((r/Rmax)) * pow(cos(theta), 2.0) + 2.0 * E_prime((r/Rmax)) * T((r/Rmax)) * sin(theta) * sin(2.0 * theta) - 2.0 * E_prime((r/Rmax)) * T((r/Rmax)) * cos(theta) * cos(2.0 * theta) - P((r/Rmax)) * T_prime((r/Rmax)) * sin(theta) * sin(2.0 * theta) + P((r/Rmax)) * T_prime((r/Rmax)) * cos(theta) * cos(2.0 * theta) - P((r/Rmax)) * dP((r/Rmax)) * pow(sin(theta), 2.0) - P((r/Rmax)) * dP((r/Rmax)) * pow(cos(theta), 2.0) + P((r/Rmax)) * pow(sin(theta), 2.0) + P((r/Rmax)) * pow(cos(theta), 2.0) + 2.0 * T((r/Rmax)) * T_prime((r/Rmax)) * pow(sin(2.0 * theta), 2.0) + 2.0 * T((r/Rmax)) * T_prime((r/Rmax)) * pow(cos(2.0 * theta), 2.0) + 2.0 * T((r/Rmax)) * dP((r/Rmax)) * sin(theta) * sin(2.0 * theta) - 2.0 * T((r/Rmax)) * dP((r/Rmax)) * cos(theta) * cos(2.0 * theta) - 2.0 * T((r/Rmax)) * sin(theta) * sin(2.0 * theta) + 2.0 * T((r/Rmax)) * cos(theta) * cos(2.0 * theta));
// }

inline double PolarR6GyroZoniShiftedCulham::coeffs1(double r, double Rmax) const
{ // With Rmax=1, equals alpha(r) from equation (18) of Bourne et al. https://doi.org/10.1016/j.jcp.2023.112249
    return exp(-tanh(20.0 * (r/Rmax) - 14.0));
}

inline double PolarR6GyroZoniShiftedCulham::coeffs2(double r, double Rmax) const
{
    return 1.0 * exp(tanh(20.0 * (r/Rmax) - 14.0));
}

inline double PolarR6GyroZoniShiftedCulham::rho_glob(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const
{
    return 0.4096 * pow((r/Rmax), 6.0) * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta);
}

inline double PolarR6GyroZoniShiftedCulham::rho_pole(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const
{
    return 0.0;
}

inline double PolarR6GyroZoniShiftedCulham::phi_exact(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax) const
{
    return 0.0;
}




inline double PolarR6GyroZoniShiftedCulham::x(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return (r/Rmax) * cos_theta + Delta((r/Rmax)) - E((r/Rmax)) * cos_theta - P((r/Rmax)) * cos_theta + T((r/Rmax)) * cos(2.0 * theta) + 5.0;
}

inline double PolarR6GyroZoniShiftedCulham::y(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return (r/Rmax) * sin_theta - E((r/Rmax)) * sin_theta - P((r/Rmax)) * sin_theta - T((r/Rmax)) * sin(2.0 * theta);
}

inline double PolarR6GyroZoniShiftedCulham::J_rr(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return (Delta_prime((r/Rmax)) - E_prime((r/Rmax)) * cos_theta + T_prime((r/Rmax)) * cos(2.0 * theta) - dP((r/Rmax)) * cos_theta + cos_theta)/Rmax;
}

inline double PolarR6GyroZoniShiftedCulham::J_rt(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return (-(r/Rmax)) * sin_theta + E((r/Rmax)) * sin_theta + P((r/Rmax)) * sin_theta - 2.0 * T((r/Rmax)) * sin(2.0 * theta);
}

inline double PolarR6GyroZoniShiftedCulham::J_tr(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return ((-E_prime((r/Rmax))) * sin_theta - T_prime((r/Rmax)) * sin(2.0 * theta) - dP((r/Rmax)) * sin_theta + sin_theta)/Rmax;
}

inline double PolarR6GyroZoniShiftedCulham::J_tt(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return (r/Rmax) * cos_theta - E((r/Rmax)) * cos_theta - P((r/Rmax)) * cos_theta - 2.0 * T((r/Rmax)) * cos(2.0 * theta);
}

// inline double PolarR6GyroZoniShiftedCulham::J_xr(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, double sin_theta, double cos_theta) const
// {
//     return ((-(r/Rmax)) * cos_theta + E((r/Rmax)) * cos_theta + P((r/Rmax)) * cos_theta + 2.0 * T((r/Rmax)) * cos(2.0 * theta)) / ((-(r/Rmax)) * Delta_prime((r/Rmax)) * cos_theta + (r/Rmax) * E_prime((r/Rmax)) * pow(sin_theta, 2.0) + (r/Rmax) * E_prime((r/Rmax)) * pow(cos_theta, 2.0) + (r/Rmax) * T_prime((r/Rmax)) * sin_theta * sin(2.0 * theta) - (r/Rmax) * T_prime((r/Rmax)) * cos_theta * cos(2.0 * theta) + (r/Rmax) * dP((r/Rmax)) * pow(sin_theta, 2.0) + (r/Rmax) * dP((r/Rmax)) * pow(cos_theta, 2.0) - (r/Rmax) * pow(sin_theta, 2.0) - (r/Rmax) * pow(cos_theta, 2.0) + Delta_prime((r/Rmax)) * E((r/Rmax)) * cos_theta + Delta_prime((r/Rmax)) * P((r/Rmax)) * cos_theta + 2.0 * Delta_prime((r/Rmax)) * T((r/Rmax)) * cos(2.0 * theta) - E((r/Rmax)) * E_prime((r/Rmax)) * pow(sin_theta, 2.0) - E((r/Rmax)) * E_prime((r/Rmax)) * pow(cos_theta, 2.0) - E((r/Rmax)) * T_prime((r/Rmax)) * sin_theta * sin(2.0 * theta) + E((r/Rmax)) * T_prime((r/Rmax)) * cos_theta * cos(2.0 * theta) - E((r/Rmax)) * dP((r/Rmax)) * pow(sin_theta, 2.0) - E((r/Rmax)) * dP((r/Rmax)) * pow(cos_theta, 2.0) + E((r/Rmax)) * pow(sin_theta, 2.0) + E((r/Rmax)) * pow(cos_theta, 2.0) - E_prime((r/Rmax)) * P((r/Rmax)) * pow(sin_theta, 2.0) - E_prime((r/Rmax)) * P((r/Rmax)) * pow(cos_theta, 2.0) + 2.0 * E_prime((r/Rmax)) * T((r/Rmax)) * sin_theta * sin(2.0 * theta) - 2.0 * E_prime((r/Rmax)) * T((r/Rmax)) * cos_theta * cos(2.0 * theta) - P((r/Rmax)) * T_prime((r/Rmax)) * sin_theta * sin(2.0 * theta) + P((r/Rmax)) * T_prime((r/Rmax)) * cos_theta * cos(2.0 * theta) - P((r/Rmax)) * dP((r/Rmax)) * pow(sin_theta, 2.0) - P((r/Rmax)) * dP((r/Rmax)) * pow(cos_theta, 2.0) + P((r/Rmax)) * pow(sin_theta, 2.0) + P((r/Rmax)) * pow(cos_theta, 2.0) + 2.0 * T((r/Rmax)) * T_prime((r/Rmax)) * pow(sin(2.0 * theta), 2.0) + 2.0 * T((r/Rmax)) * T_prime((r/Rmax)) * pow(cos(2.0 * theta), 2.0) + 2.0 * T((r/Rmax)) * dP((r/Rmax)) * sin_theta * sin(2.0 * theta) - 2.0 * T((r/Rmax)) * dP((r/Rmax)) * cos_theta * cos(2.0 * theta) - 2.0 * T((r/Rmax)) * sin_theta * sin(2.0 * theta) + 2.0 * T((r/Rmax)) * cos_theta * cos(2.0 * theta));
// }

// inline double PolarR6GyroZoniShiftedCulham::J_xq(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, double sin_theta, double cos_theta) const
// {
//     return ((-(r/Rmax)) * sin_theta + E((r/Rmax)) * sin_theta + P((r/Rmax)) * sin_theta - 2.0 * T((r/Rmax)) * sin(2.0 * theta)) / ((-(r/Rmax)) * Delta_prime((r/Rmax)) * cos_theta + (r/Rmax) * E_prime((r/Rmax)) * pow(sin_theta, 2.0) + (r/Rmax) * E_prime((r/Rmax)) * pow(cos_theta, 2.0) + (r/Rmax) * T_prime((r/Rmax)) * sin_theta * sin(2.0 * theta) - (r/Rmax) * T_prime((r/Rmax)) * cos_theta * cos(2.0 * theta) + (r/Rmax) * dP((r/Rmax)) * pow(sin_theta, 2.0) + (r/Rmax) * dP((r/Rmax)) * pow(cos_theta, 2.0) - (r/Rmax) * pow(sin_theta, 2.0) - (r/Rmax) * pow(cos_theta, 2.0) + Delta_prime((r/Rmax)) * E((r/Rmax)) * cos_theta + Delta_prime((r/Rmax)) * P((r/Rmax)) * cos_theta + 2.0 * Delta_prime((r/Rmax)) * T((r/Rmax)) * cos(2.0 * theta) - E((r/Rmax)) * E_prime((r/Rmax)) * pow(sin_theta, 2.0) - E((r/Rmax)) * E_prime((r/Rmax)) * pow(cos_theta, 2.0) - E((r/Rmax)) * T_prime((r/Rmax)) * sin_theta * sin(2.0 * theta) + E((r/Rmax)) * T_prime((r/Rmax)) * cos_theta * cos(2.0 * theta) - E((r/Rmax)) * dP((r/Rmax)) * pow(sin_theta, 2.0) - E((r/Rmax)) * dP((r/Rmax)) * pow(cos_theta, 2.0) + E((r/Rmax)) * pow(sin_theta, 2.0) + E((r/Rmax)) * pow(cos_theta, 2.0) - E_prime((r/Rmax)) * P((r/Rmax)) * pow(sin_theta, 2.0) - E_prime((r/Rmax)) * P((r/Rmax)) * pow(cos_theta, 2.0) + 2.0 * E_prime((r/Rmax)) * T((r/Rmax)) * sin_theta * sin(2.0 * theta) - 2.0 * E_prime((r/Rmax)) * T((r/Rmax)) * cos_theta * cos(2.0 * theta) - P((r/Rmax)) * T_prime((r/Rmax)) * sin_theta * sin(2.0 * theta) + P((r/Rmax)) * T_prime((r/Rmax)) * cos_theta * cos(2.0 * theta) - P((r/Rmax)) * dP((r/Rmax)) * pow(sin_theta, 2.0) - P((r/Rmax)) * dP((r/Rmax)) * pow(cos_theta, 2.0) + P((r/Rmax)) * pow(sin_theta, 2.0) + P((r/Rmax)) * pow(cos_theta, 2.0) + 2.0 * T((r/Rmax)) * T_prime((r/Rmax)) * pow(sin(2.0 * theta), 2.0) + 2.0 * T((r/Rmax)) * T_prime((r/Rmax)) * pow(cos(2.0 * theta), 2.0) + 2.0 * T((r/Rmax)) * dP((r/Rmax)) * sin_theta * sin(2.0 * theta) - 2.0 * T((r/Rmax)) * dP((r/Rmax)) * cos_theta * cos(2.0 * theta) - 2.0 * T((r/Rmax)) * sin_theta * sin(2.0 * theta) + 2.0 * T((r/Rmax)) * cos_theta * cos(2.0 * theta));
// }

// inline double PolarR6GyroZoniShiftedCulham::J_yr(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, double sin_theta, double cos_theta) const
// {
//     return ((-E_prime((r/Rmax))) * sin_theta - T_prime((r/Rmax)) * sin(2.0 * theta) - dP((r/Rmax)) * sin_theta + sin_theta) / ((-(r/Rmax)) * Delta_prime((r/Rmax)) * cos_theta + (r/Rmax) * E_prime((r/Rmax)) * pow(sin_theta, 2.0) + (r/Rmax) * E_prime((r/Rmax)) * pow(cos_theta, 2.0) + (r/Rmax) * T_prime((r/Rmax)) * sin_theta * sin(2.0 * theta) - (r/Rmax) * T_prime((r/Rmax)) * cos_theta * cos(2.0 * theta) + (r/Rmax) * dP((r/Rmax)) * pow(sin_theta, 2.0) + (r/Rmax) * dP((r/Rmax)) * pow(cos_theta, 2.0) - (r/Rmax) * pow(sin_theta, 2.0) - (r/Rmax) * pow(cos_theta, 2.0) + Delta_prime((r/Rmax)) * E((r/Rmax)) * cos_theta + Delta_prime((r/Rmax)) * P((r/Rmax)) * cos_theta + 2.0 * Delta_prime((r/Rmax)) * T((r/Rmax)) * cos(2.0 * theta) - E((r/Rmax)) * E_prime((r/Rmax)) * pow(sin_theta, 2.0) - E((r/Rmax)) * E_prime((r/Rmax)) * pow(cos_theta, 2.0) - E((r/Rmax)) * T_prime((r/Rmax)) * sin_theta * sin(2.0 * theta) + E((r/Rmax)) * T_prime((r/Rmax)) * cos_theta * cos(2.0 * theta) - E((r/Rmax)) * dP((r/Rmax)) * pow(sin_theta, 2.0) - E((r/Rmax)) * dP((r/Rmax)) * pow(cos_theta, 2.0) + E((r/Rmax)) * pow(sin_theta, 2.0) + E((r/Rmax)) * pow(cos_theta, 2.0) - E_prime((r/Rmax)) * P((r/Rmax)) * pow(sin_theta, 2.0) - E_prime((r/Rmax)) * P((r/Rmax)) * pow(cos_theta, 2.0) + 2.0 * E_prime((r/Rmax)) * T((r/Rmax)) * sin_theta * sin(2.0 * theta) - 2.0 * E_prime((r/Rmax)) * T((r/Rmax)) * cos_theta * cos(2.0 * theta) - P((r/Rmax)) * T_prime((r/Rmax)) * sin_theta * sin(2.0 * theta) + P((r/Rmax)) * T_prime((r/Rmax)) * cos_theta * cos(2.0 * theta) - P((r/Rmax)) * dP((r/Rmax)) * pow(sin_theta, 2.0) - P((r/Rmax)) * dP((r/Rmax)) * pow(cos_theta, 2.0) + P((r/Rmax)) * pow(sin_theta, 2.0) + P((r/Rmax)) * pow(cos_theta, 2.0) + 2.0 * T((r/Rmax)) * T_prime((r/Rmax)) * pow(sin(2.0 * theta), 2.0) + 2.0 * T((r/Rmax)) * T_prime((r/Rmax)) * pow(cos(2.0 * theta), 2.0) + 2.0 * T((r/Rmax)) * dP((r/Rmax)) * sin_theta * sin(2.0 * theta) - 2.0 * T((r/Rmax)) * dP((r/Rmax)) * cos_theta * cos(2.0 * theta) - 2.0 * T((r/Rmax)) * sin_theta * sin(2.0 * theta) + 2.0 * T((r/Rmax)) * cos_theta * cos(2.0 * theta));
// }

// inline double PolarR6GyroZoniShiftedCulham::J_yq(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, double sin_theta, double cos_theta) const
// {
//     return ((-Delta_prime((r/Rmax))) + E_prime((r/Rmax)) * cos_theta - T_prime((r/Rmax)) * cos(2.0 * theta) + dP((r/Rmax)) * cos_theta - cos_theta) / ((-(r/Rmax)) * Delta_prime((r/Rmax)) * cos_theta + (r/Rmax) * E_prime((r/Rmax)) * pow(sin_theta, 2.0) + (r/Rmax) * E_prime((r/Rmax)) * pow(cos_theta, 2.0) + (r/Rmax) * T_prime((r/Rmax)) * sin_theta * sin(2.0 * theta) - (r/Rmax) * T_prime((r/Rmax)) * cos_theta * cos(2.0 * theta) + (r/Rmax) * dP((r/Rmax)) * pow(sin_theta, 2.0) + (r/Rmax) * dP((r/Rmax)) * pow(cos_theta, 2.0) - (r/Rmax) * pow(sin_theta, 2.0) - (r/Rmax) * pow(cos_theta, 2.0) + Delta_prime((r/Rmax)) * E((r/Rmax)) * cos_theta + Delta_prime((r/Rmax)) * P((r/Rmax)) * cos_theta + 2.0 * Delta_prime((r/Rmax)) * T((r/Rmax)) * cos(2.0 * theta) - E((r/Rmax)) * E_prime((r/Rmax)) * pow(sin_theta, 2.0) - E((r/Rmax)) * E_prime((r/Rmax)) * pow(cos_theta, 2.0) - E((r/Rmax)) * T_prime((r/Rmax)) * sin_theta * sin(2.0 * theta) + E((r/Rmax)) * T_prime((r/Rmax)) * cos_theta * cos(2.0 * theta) - E((r/Rmax)) * dP((r/Rmax)) * pow(sin_theta, 2.0) - E((r/Rmax)) * dP((r/Rmax)) * pow(cos_theta, 2.0) + E((r/Rmax)) * pow(sin_theta, 2.0) + E((r/Rmax)) * pow(cos_theta, 2.0) - E_prime((r/Rmax)) * P((r/Rmax)) * pow(sin_theta, 2.0) - E_prime((r/Rmax)) * P((r/Rmax)) * pow(cos_theta, 2.0) + 2.0 * E_prime((r/Rmax)) * T((r/Rmax)) * sin_theta * sin(2.0 * theta) - 2.0 * E_prime((r/Rmax)) * T((r/Rmax)) * cos_theta * cos(2.0 * theta) - P((r/Rmax)) * T_prime((r/Rmax)) * sin_theta * sin(2.0 * theta) + P((r/Rmax)) * T_prime((r/Rmax)) * cos_theta * cos(2.0 * theta) - P((r/Rmax)) * dP((r/Rmax)) * pow(sin_theta, 2.0) - P((r/Rmax)) * dP((r/Rmax)) * pow(cos_theta, 2.0) + P((r/Rmax)) * pow(sin_theta, 2.0) + P((r/Rmax)) * pow(cos_theta, 2.0) + 2.0 * T((r/Rmax)) * T_prime((r/Rmax)) * pow(sin(2.0 * theta), 2.0) + 2.0 * T((r/Rmax)) * T_prime((r/Rmax)) * pow(cos(2.0 * theta), 2.0) + 2.0 * T((r/Rmax)) * dP((r/Rmax)) * sin_theta * sin(2.0 * theta) - 2.0 * T((r/Rmax)) * dP((r/Rmax)) * cos_theta * cos(2.0 * theta) - 2.0 * T((r/Rmax)) * sin_theta * sin(2.0 * theta) + 2.0 * T((r/Rmax)) * cos_theta * cos(2.0 * theta));
// }

inline double PolarR6GyroZoniShiftedCulham::rho_glob(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return 0.4096 * pow((r/Rmax), 6.0) * pow(((r/Rmax) - 1.0), 6.0) * cos(11.0 * theta);
}

inline double PolarR6GyroZoniShiftedCulham::rho_pole(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return 0.0;
}

inline double PolarR6GyroZoniShiftedCulham::phi_exact(double r, double theta, double map3_unused_1, double map3_unused_2, double Rmax, double sin_theta, double cos_theta) const
{
    return 0.0;
}
