#include "../include/InputFunctions/DomainGeometry/culhamGeometry.h"

CulhamGeometry::CulhamGeometry(const double& Rmax):
    Rmax(Rmax)
{
    initializeGeometry();
}


CulhamGeometry::CulhamGeometry(){
    initializeGeometry();
}

void CulhamGeometry::initializeGeometry() {
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
    for (i = 1; i < 1000; i += 1){
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




// // In earlier versions denoted by 'x'
// double CulhamGeometry::Fx(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
//     return (r/Rmax) * cos_theta + Delta((r/Rmax)) - E((r/Rmax)) * cos_theta - P((r/Rmax)) * cos_theta + T((r/Rmax)) * cos(2.0 * theta) + 5.0;
// }

// // In earlier versions denoted by 'y'
// double CulhamGeometry::Fy(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
//     return (r/Rmax) * sin_theta - E((r/Rmax)) * sin_theta - P((r/Rmax)) * sin_theta - T((r/Rmax)) * sin(2.0 * theta);
// }


// // In earlier versions denoted by 'Jrr'
// double CulhamGeometry::dFx_dr(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
//     return (Delta_prime((r/Rmax)) - E_prime((r/Rmax)) * cos_theta + T_prime((r/Rmax)) * cos(2.0 * theta) - dP((r/Rmax)) * cos_theta + cos_theta)/Rmax;
// }

// // In earlier versions denoted by 'Jtr'
// double CulhamGeometry::dFy_dr(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
//     return ((-E_prime((r/Rmax))) * sin_theta - T_prime((r/Rmax)) * sin(2.0 * theta) - dP((r/Rmax)) * sin_theta + sin_theta)/Rmax;
// }

// // In earlier versions denoted by 'Jrt'
// double CulhamGeometry::dFx_dt(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
//     return (-(r/Rmax)) * sin_theta + E((r/Rmax)) * sin_theta + P((r/Rmax)) * sin_theta - 2.0 * T((r/Rmax)) * sin(2.0 * theta);
// }

// // In earlier versions denoted by 'Jtt'
// double CulhamGeometry::dFy_dt(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
//     return (r/Rmax) * cos_theta - E((r/Rmax)) * cos_theta - P((r/Rmax)) * cos_theta - 2.0 * T((r/Rmax)) * cos(2.0 * theta);
// }


// double CulhamGeometry::my_sum(std::array<double, 1001>& f, int64_t start_idx, int64_t end_idx) const
// {
//     int64_t i;
//     double result;
//     result = 0.0;
//     #pragma omp parallel for reduction(+: result)
//     for (i = start_idx; i < end_idx; i += 1)
//     {
//         result += f[i];
//     }
//     return result;
// }

// double CulhamGeometry::q(double rr) const
// {
//     return 0.8 - 0.1 * (rr * rr);
// }

// double CulhamGeometry::dq(double rr) const
// {
//     return (-0.2) * rr;
// }

// double CulhamGeometry::p(double rr) const
// {
//     return 100000.0 - 90000.0 * (rr * rr);
// }

// double CulhamGeometry::dp(double rr) const
// {
//     return (-180000.0) * rr;
// }

// double CulhamGeometry::dg(double rr, double g) const
// {
//     return ((-g) * (0.0625000000000001 * (rr * rr) / pow((1.0 - 0.125 * (rr * rr)), 2.0) + 2.0 / (4.0 - 0.5 * (rr * rr))) + 2.261946711816e-06 * (4.0 - 0.5 * (rr * rr)) / g) / (rr / (4.0 - 0.5 * (rr * rr)) + (4.0 - 0.5 * (rr * rr)) / (g * rr));
// }

// double CulhamGeometry::double_deriv(double rr, double c, double g, double dg, double val, double d_val) const
// {
//     return c * val / (rr * rr) - d_val * (pow(rr, (double)((-1))) + (4.0 - 0.5 * (rr * rr)) * (2.0 * dg * rr / (4.0 - 0.5 * (rr * rr)) + 0.125 * g * (rr * rr) / pow((1.0 - 0.125 * (rr * rr)), 2.0) + 2.0 * g / (4.0 - 0.5 * (rr * rr))) / (g * rr));
// }

// double CulhamGeometry::g(double rr) const
// {
//     int64_t ri;
//     double dr;
//     double m;
//     double c;
//     ri = (int64_t)(rr * 1000 / 1.0);
//     dr = 1.0 / 1000.0;
//     if (ri == 1000)
//     {
//         return g_array[ri];
//     }
//     else
//     {
//         m = (g_array[ri + 1] - g_array[ri]) / dr;
//         c = g_array[ri] - m * ri * dr;
//         return m * rr + c;
//     }
// }

// double CulhamGeometry::Delta(double rr) const
// {
//     int64_t ri;
//     double dr;
//     double m;
//     double c;
//     ri = (int64_t)(rr * 1000 / 1.0);
//     dr = 1.0 / 1000.0;
//     if (ri == 1000)
//     {
//         return Delta_array[ri];
//     }
//     else
//     {
//         m = (Delta_array[ri + 1] - Delta_array[ri]) / dr;
//         c = Delta_array[ri] - m * ri * dr;
//         return m * rr + c;
//     }
// }

// double CulhamGeometry::Delta_prime(double rr) const
// {
//     int64_t ri;
//     double dr;
//     double m;
//     double c;
//     ri = (int64_t)(rr * 1000 / 1.0);
//     dr = 1.0 / 1000.0;
//     if (ri == 1000)
//     {
//         return Delta_prime_array[ri];
//     }
//     else
//     {
//         m = (Delta_prime_array[ri + 1] - Delta_prime_array[ri]) / dr;
//         c = Delta_prime_array[ri] - m * ri * dr;
//         return m * rr + c;
//     }
// }

// double CulhamGeometry::E(double rr) const
// {
//     int64_t ri;
//     double dr;
//     double m;
//     double c;
//     ri = (int64_t)(rr * 1000 / 1.0);
//     dr = 1.0 / 1000.0;
//     if (ri == 1000)
//     {
//         return E_array[ri];
//     }
//     else
//     {
//         m = (E_array[ri + 1] - E_array[ri]) / dr;
//         c = E_array[ri] - m * ri * dr;
//         return m * rr + c;
//     }
// }

// double CulhamGeometry::T(double rr) const
// {
//     int64_t ri;
//     double dr;
//     double m;
//     double c;
//     ri = (int64_t)(rr * 1000 / 1.0);
//     dr = 1.0 / 1000.0;
//     if (ri == 1000)
//     {
//         return T_array[ri];
//     }
//     else
//     {
//         m = (T_array[ri + 1] - T_array[ri]) / dr;
//         c = T_array[ri] - m * ri * dr;
//         return m * rr + c;
//     }
// }

// double CulhamGeometry::E_prime(double rr) const
// {
//     int64_t ri;
//     double dr;
//     double m;
//     double c;
//     ri = (int64_t)(rr * 1000 / 1.0);
//     dr = 1.0 / 1000.0;
//     if (ri == 1000)
//     {
//         return E_prime_array[ri];
//     }
//     else
//     {
//         m = (E_prime_array[ri + 1] - E_prime_array[ri]) / dr;
//         c = E_prime_array[ri] - m * ri * dr;
//         return m * rr + c;
//     }
// }

// double CulhamGeometry::T_prime(double rr) const
// {
//     int64_t ri;
//     double dr;
//     double m;
//     double c;
//     ri = (int64_t)(rr * 1000 / 1.0);
//     dr = 1.0 / 1000.0;
//     if (ri == 1000)
//     {
//         return T_prime_array[ri];
//     }
//     else
//     {
//         m = (T_prime_array[ri + 1] - T_prime_array[ri]) / dr;
//         c = T_prime_array[ri] - m * ri * dr;
//         return m * rr + c;
//     }
// }

// double CulhamGeometry::P(double rr) const
// {
//     if (rr == 0)
//     {
//         return 0.0;
//     }
//     else
//     {
//         return 0.005 * pow(rr, 3.0) + 0.1 * rr * Delta(rr) - 0.1 * pow(E(rr), 2.0) - pow(T(rr), 2.0) / rr;
//     }
// }

// double CulhamGeometry::dP(double rr) const
// {
//     if (rr == 0)
//     {
//         return 0.0;
//     }
//     else
//     {
//         return 0.015 * (rr * rr) + 0.1 * rr * Delta_prime(rr) + 0.1 * Delta(rr) - 0.2 * E(rr) * E_prime(rr) - 2.0 * T(rr) * T_prime(rr) / rr + pow(T(rr), 2.0) / (rr * rr);
//     }
// }
