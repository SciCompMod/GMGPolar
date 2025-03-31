#include "../include/InputFunctions/DomainGeometry/culhamGeometry.h"

CulhamGeometry::CulhamGeometry(const double& Rmax)
    : Rmax(Rmax)
{
    initializeGeometry();
}

CulhamGeometry::CulhamGeometry()
{
    initializeGeometry();
}

void CulhamGeometry::initializeGeometry()
{
    rr               = 0.0;
    dr               = 1.0 / 1000;
    dr_h             = dr * 0.5;
    g_array[0]       = 1.0;
    E_array[0]       = 0.0;
    E_prime_array[0] = 1.0;
    T_array[0]       = 0.0;
    T_prime_array[0] = 0.0;
    r[0]             = rr;
    rr += dr;
    r[1]             = rr;
    g_array[1]       = 1.0;
    E_array[1]       = rr;
    E_prime_array[1] = 1.0;
    T_array[1]       = rr * rr;
    T_prime_array[1] = 2 * rr;
    for (i = 1; i < 1000; i += 1) {
        /*Step 1*/
        dg_1  = dg(rr, g_array[i]);
        dE_1  = E_prime_array[i];
        dT_1  = T_prime_array[i];
        ddE_1 = double_deriv(rr, 3.0, g_array[i], dg_1, E_array[i], E_prime_array[i]);
        ddT_1 = double_deriv(rr, 8.0, g_array[i], dg_1, T_array[i], T_prime_array[i]);
        /*Step 2*/
        r2    = rr + dr_h;
        g_2   = g_array[i] + dr_h * dg_1;
        dg_2  = dg(r2, g_2);
        E_2   = E_array[i] + dE_1 * dr_h;
        T_2   = T_array[i] + dT_1 * dr_h;
        dE_2  = E_prime_array[i] + ddE_1 * dr_h;
        dT_2  = T_prime_array[i] + ddT_1 * dr_h;
        ddE_2 = double_deriv(r2, 3.0, g_2, dg_2, E_2, dE_2);
        ddT_2 = double_deriv(r2, 8.0, g_2, dg_2, T_2, dT_2);
        /*Step 3*/
        g_3   = g_array[i] + dr_h * dg_2;
        dg_3  = dg(r2, g_3);
        E_3   = E_array[i] + dE_2 * dr_h;
        T_3   = T_array[i] + dT_2 * dr_h;
        dE_3  = E_prime_array[i] + ddE_2 * dr_h;
        dT_3  = T_prime_array[i] + ddT_2 * dr_h;
        ddE_3 = double_deriv(r2, 3.0, g_3, dg_3, E_3, dE_3);
        ddT_3 = double_deriv(r2, 8.0, g_3, dg_3, T_3, dT_3);
        /*Step 4*/
        rr                   = rr + dr;
        g_4                  = g_array[i] + dr * dg_3;
        dg_4                 = dg(rr, g_4);
        E_4                  = E_array[i] + dE_3 * dr;
        T_4                  = T_array[i] + dT_3 * dr;
        dE_4                 = E_prime_array[i] + ddE_3 * dr;
        dT_4                 = T_prime_array[i] + ddT_3 * dr;
        ddE_4                = double_deriv(rr, 3.0, g_4, dg_4, E_4, dE_4);
        ddT_4                = double_deriv(rr, 8.0, g_4, dg_4, T_4, dT_4);
        g_array[i + 1]       = g_array[i] + dr * (dg_1 + 2 * dg_2 + 2 * dg_3 + dg_4) / 6.0;
        E_array[i + 1]       = E_array[i] + dr * (dE_1 + 2 * dE_2 + 2 * dE_3 + dE_4) / 6.0;
        T_array[i + 1]       = T_array[i] + dr * (dT_1 + 2 * dT_2 + 2 * dT_3 + dT_4) / 6.0;
        E_prime_array[i + 1] = E_prime_array[i] + dr * (ddE_1 + 2 * ddE_2 + 2 * ddE_3 + ddE_4) / 6.0;
        T_prime_array[i + 1] = T_prime_array[i] + dr * (ddT_1 + 2 * ddT_2 + 2 * ddT_3 + ddT_4) / 6.0;
        r[i + 1]             = rr;
    }
    current_Ea = E(1.0);
    current_Ta = T(1.0);
    for (i_0001 = 0; i_0001 < E_array.size(); i_0001 += 1) {
        E_array[i_0001] = 0.25 * E_array[i_0001] / current_Ea;
    }
    for (i_0001 = 0; i_0001 < T_array.size(); i_0001 += 1) {
        T_array[i_0001] = 0.1 * T_array[i_0001] / current_Ta;
    }
    for (i_0001 = 0; i_0001 < E_prime_array.size(); i_0001 += 1) {
        E_prime_array[i_0001] = 0.25 * E_prime_array[i_0001] / current_Ea;
    }
    for (i_0001 = 0; i_0001 < T_prime_array.size(); i_0001 += 1) {
        T_prime_array[i_0001] = 0.1 * T_prime_array[i_0001] / current_Ta;
    }
    for (i_0001 = 0; i_0001 < f.size(); i_0001 += 1) {
        f[i_0001] = r[i_0001] * g_array[i_0001] / 5.0 / q(r[i_0001]);
    }
    for (i_0001 = 0; i_0001 < integ_contents.size(); i_0001 += 1) {
        integ_contents[i_0001] = r[i_0001] * (f[i_0001] * f[i_0001]) -
                                 2 * (r[i_0001] * r[i_0001]) * 1.25663706212e-06 * dp(r[i_0001]) / (1.0 * 1.0);
    }
    Delta_prime_array[0] = 0;
    Delta_array[0]       = 0;
    integral             = 0.0;
    for (i = 1; i < 1000 + 1; i += 1) {
        integral += dr * (integ_contents[i - 1] + integ_contents[i]) * 0.5;
        Delta_prime_array[i] = (-integral) / (5.0 * r[i] * pow(f[i], 2.0));
        Delta_array[i]       = Delta_array[i - 1] + dr * 0.5 * (Delta_prime_array[i - 1] + Delta_prime_array[i]);
    }
    current_Delta_a = Delta(1.0);
    for (i_0001 = 0; i_0001 < Delta_array.size(); i_0001 += 1) {
        Delta_array[i_0001] -= current_Delta_a;
    }
}
