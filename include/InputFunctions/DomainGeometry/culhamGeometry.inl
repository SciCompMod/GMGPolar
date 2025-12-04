#pragma once

#include "culhamGeometry.h"

// In earlier versions denoted by 'x'
inline double CulhamGeometry::Fx(double r, double theta) const {
    double sin_theta = std::sin(theta);
    double cos_theta = std::cos(theta);
    const double cos_two_theta = 1.0 - sin_theta * sin_theta;
    return (r/Rmax) * cos_theta + Delta((r/Rmax)) - E((r/Rmax)) * cos_theta - P((r/Rmax)) * cos_theta + T((r/Rmax)) * cos_two_theta + 5.0;
}

// In earlier versions denoted by 'y'
inline double CulhamGeometry::Fy(double r, double theta) const {
    double sin_theta = std::sin(theta);
    double cos_theta = std::cos(theta);
    const double sin_two_theta = 2.0 * sin_theta * cos_theta;
    return (r/Rmax) * sin_theta - E((r/Rmax)) * sin_theta - P((r/Rmax)) * sin_theta - T((r/Rmax)) * sin_two_theta;
}


// In earlier versions denoted by 'Jrr'
inline double CulhamGeometry::dFx_dr(double r, double theta) const {
    double sin_theta = std::sin(theta);
    double cos_theta = std::cos(theta);
    const double cos_two_theta = 1.0 - sin_theta * sin_theta;
    return (Delta_prime((r/Rmax)) - E_prime((r/Rmax)) * cos_theta + T_prime((r/Rmax)) * cos_two_theta - dP((r/Rmax)) * cos_theta + cos_theta)/Rmax;
}

// In earlier versions denoted by 'Jtr'
inline double CulhamGeometry::dFy_dr(double r, double theta) const {
    double sin_theta = std::sin(theta);
    double cos_theta = std::cos(theta);
    const double sin_two_theta = 2.0 * sin_theta * cos_theta;
    return ((-E_prime((r/Rmax))) * sin_theta - T_prime((r/Rmax)) * sin_two_theta - dP((r/Rmax)) * sin_theta + sin_theta)/Rmax;
}

// In earlier versions denoted by 'Jrt'
inline double CulhamGeometry::dFx_dt(double r, double theta) const {
    double sin_theta = std::sin(theta);
    double cos_theta = std::cos(theta);
    const double sin_two_theta = 2.0 * sin_theta * cos_theta;
    return (-(r/Rmax)) * sin_theta + E((r/Rmax)) * sin_theta + P((r/Rmax)) * sin_theta - 2.0 * T((r/Rmax)) * sin_two_theta;
}

// In earlier versions denoted by 'Jtt'
inline double CulhamGeometry::dFy_dt(double r, double theta) const {
    double sin_theta = std::sin(theta);
    double cos_theta = std::cos(theta);
    const double cos_two_theta = 1.0 - sin_theta * sin_theta;
    return (r/Rmax) * cos_theta - E((r/Rmax)) * cos_theta - P((r/Rmax)) * cos_theta - 2.0 * T((r/Rmax)) * cos_two_theta;
}



inline double CulhamGeometry::my_sum(std::array<double, 1001>& f, int64_t start_idx, int64_t end_idx) const
{
    int64_t i;
    double result;
    result = 0.0;
    for (i = start_idx; i < end_idx; i += 1)
    {
        result += f[i];
    }
    return result;
}

inline double CulhamGeometry::q(double rr) const
{
    return 0.8 - 0.1 * (rr * rr);
}

inline double CulhamGeometry::dq(double rr) const
{
    return (-0.2) * rr;
}

inline double CulhamGeometry::p(double rr) const
{
    return 100000.0 - 90000.0 * (rr * rr);
}

inline double CulhamGeometry::dp(double rr) const
{
    return (-180000.0) * rr;
}

inline double CulhamGeometry::dg(double rr, double g) const
{
    return ((-g) * (0.0625000000000001 * (rr * rr) / pow((1.0 - 0.125 * (rr * rr)), 2.0) + 2.0 / (4.0 - 0.5 * (rr * rr))) + 2.261946711816e-06 * (4.0 - 0.5 * (rr * rr)) / g) / (rr / (4.0 - 0.5 * (rr * rr)) + (4.0 - 0.5 * (rr * rr)) / (g * rr));
}

inline double CulhamGeometry::double_deriv(double rr, double c, double g, double dg, double val, double d_val) const
{
    return c * val / (rr * rr) - d_val * (pow(rr, (double)((-1))) + (4.0 - 0.5 * (rr * rr)) * (2.0 * dg * rr / (4.0 - 0.5 * (rr * rr)) + 0.125 * g * (rr * rr) / pow((1.0 - 0.125 * (rr * rr)), 2.0) + 2.0 * g / (4.0 - 0.5 * (rr * rr))) / (g * rr));
}

inline double CulhamGeometry::g(double rr) const
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

inline double CulhamGeometry::Delta(double rr) const
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

inline double CulhamGeometry::Delta_prime(double rr) const
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

inline double CulhamGeometry::E(double rr) const
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

inline double CulhamGeometry::T(double rr) const
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

inline double CulhamGeometry::E_prime(double rr) const
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

inline double CulhamGeometry::T_prime(double rr) const
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

inline double CulhamGeometry::P(double rr) const
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

inline double CulhamGeometry::dP(double rr) const
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
