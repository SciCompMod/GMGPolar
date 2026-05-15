#pragma once

#include <cmath>
#include <memory>
#include <array>
#include <cstdint>
#include <Kokkos_Core.hpp>

#include "../domainGeometry.h"

namespace gmgpolar
{

class CulhamGeometry
{
public:
    CulhamGeometry();
    explicit CulhamGeometry(double Rmax);

    KOKKOS_INLINE_FUNCTION double Fx(double r, double theta) const;
    KOKKOS_INLINE_FUNCTION double Fy(double r, double theta) const;
    KOKKOS_INLINE_FUNCTION double dFx_dr(double r, double theta) const;
    KOKKOS_INLINE_FUNCTION double dFy_dr(double r, double theta) const;
    KOKKOS_INLINE_FUNCTION double dFx_dt(double r, double theta) const;
    KOKKOS_INLINE_FUNCTION double dFy_dt(double r, double theta) const;

private:
    const double Rmax = 1.3;

    void initializeGeometry();

    KOKKOS_INLINE_FUNCTION double my_sum(std::array<double, 1001>& f, int64_t start_idx, int64_t end_idx) const;
    KOKKOS_INLINE_FUNCTION double q(double rr) const;
    KOKKOS_INLINE_FUNCTION double dq(double rr) const;
    KOKKOS_INLINE_FUNCTION double p(double rr) const;
    KOKKOS_INLINE_FUNCTION double dp(double rr) const;
    KOKKOS_INLINE_FUNCTION double dg(double rr, double g) const;
    KOKKOS_INLINE_FUNCTION double double_deriv(double rr, double c, double g, double dg, double val,
                                               double d_val) const;
    KOKKOS_INLINE_FUNCTION double g(double rr) const;
    KOKKOS_INLINE_FUNCTION double Delta(double rr) const;
    KOKKOS_INLINE_FUNCTION double Delta_prime(double rr) const;
    KOKKOS_INLINE_FUNCTION double E(double rr) const;
    KOKKOS_INLINE_FUNCTION double T(double rr) const;
    KOKKOS_INLINE_FUNCTION double E_prime(double rr) const;
    KOKKOS_INLINE_FUNCTION double T_prime(double rr) const;
    KOKKOS_INLINE_FUNCTION double P(double rr) const;
    KOKKOS_INLINE_FUNCTION double dP(double rr) const;

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
    std::array<double, 1001> g_array;
    std::array<double, 1001> Delta_array;
    std::array<double, 1001> Delta_prime_array;
    std::array<double, 1001> E_array;
    std::array<double, 1001> T_array;
    std::array<double, 1001> E_prime_array;
    std::array<double, 1001> T_prime_array;
};

#include "culhamGeometry.inl"
} // namespace gmgpolar
