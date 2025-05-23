#include "../include/InputFunctions/SourceTerms/refined_ZoniShiftedGyro_CzarnyGeometry.h"

void Refined_ZoniShiftedGyro_CzarnyGeometry::initializeGeometry()
{
    factor_xi = 1.0 / sqrt(1.0 - inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon / 4.0);
}

Refined_ZoniShiftedGyro_CzarnyGeometry::Refined_ZoniShiftedGyro_CzarnyGeometry(
    const double& Rmax, const double& inverse_aspect_ratio_epsilon, const double& ellipticity_e)
    : Rmax(Rmax)
    , inverse_aspect_ratio_epsilon(inverse_aspect_ratio_epsilon)
    , ellipticity_e(ellipticity_e)
{
    initializeGeometry();
}

double Refined_ZoniShiftedGyro_CzarnyGeometry::rhs_f(const double& r, const double& theta, const double& sin_theta,
                                                     const double& cos_theta) const
{
    double temp =
        sqrt(inverse_aspect_ratio_epsilon * (inverse_aspect_ratio_epsilon + 2.0 * (r / Rmax) * cos_theta) + 1.0);
    double sin_theta_pow2 = pow(sin_theta, 2.0);
    double cos_theta_pow2 = pow(cos_theta, 2.0);
    double temp1          = pow((2.0 - temp), 2.0);
    double temp2          = pow((2.0 - temp), 3.0);

    return 1.0 *
               (((-3.33823779536505e-15) * ((r / Rmax) * (r / Rmax)) - 0.0 * (r / Rmax) - 0.0 +
                 exp((-3333.33333333333) * pow(((r / Rmax) - 0.9), 2.0))) *
                    cos(21.0 * theta) +
                (0.00184273372222541 * ((r / Rmax) * (r / Rmax)) - 0.0018029383826828 * (r / Rmax) -
                 4.00652973929511e-05 + exp((-50.0) * pow(((r / Rmax) - 0.45), 2.0))) *
                    cos(9.0 * theta)) *
               exp(tanh(20.0 * (r / Rmax) - 14.0)) -
           ((r / Rmax) *
                (((-6.67647559073009e-15) * (r / Rmax) +
                  (6000.0 - 6666.66666666667 * (r / Rmax)) * exp((-3333.33333333333) * pow(((r / Rmax) - 0.9), 2.0))) *
                     cos(21.0 * theta) +
                 (0.00368546744445083 * (r / Rmax) +
                  (45.0 - 100.0 * (r / Rmax)) * exp((-50.0) * pow(((r / Rmax) - 0.45), 2.0)) - 0.0018029383826828) *
                     cos(9.0 * theta)) *
                ((-2.0) * inverse_aspect_ratio_epsilon * sin_theta_pow2 * cos_theta / pow((temp * temp), 2.0) +
                 ((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                      (temp1 * temp) +
                  ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                     (2.0 * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) * (r / Rmax) *
                          sin_theta_pow2 * cos_theta * factor_xi / (temp1 * pow((temp * temp), (3.0 / 2.0))) -
                      4.0 * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) * (r / Rmax) *
                          sin_theta_pow2 * cos_theta * factor_xi / (temp2 * (temp * temp)) -
                      2.0 * ellipticity_e * inverse_aspect_ratio_epsilon * sin_theta_pow2 * factor_xi / (temp1 * temp) +
                      2.0 * ellipticity_e * inverse_aspect_ratio_epsilon * cos_theta_pow2 * factor_xi /
                          (temp1 * temp))) *
                exp(-tanh(20.0 * (r / Rmax) - 14.0)) /
                sqrt((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                                 (temp1 * temp) +
                             ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                                (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                     factor_xi / (temp1 * temp) +
                                 ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                            sin_theta * cos_theta / (temp * temp)),
                           2.0)) +
                     (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                               (temp1 * temp) +
                           ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                          2.0) +
                      sin_theta_pow2 / (temp * temp)) *
                         (pow((ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                   factor_xi / (temp1 * temp) +
                               ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))),
                              2.0) +
                          cos_theta_pow2 / (temp * temp))) +
            (r / Rmax) *
                (((-6.67647559073009e-15) * (r / Rmax) +
                  (6000.0 - 6666.66666666667 * (r / Rmax)) * exp((-3333.33333333333) * pow(((r / Rmax) - 0.9), 2.0))) *
                     cos(21.0 * theta) +
                 (0.00368546744445083 * (r / Rmax) +
                  (45.0 - 100.0 * (r / Rmax)) * exp((-50.0) * pow(((r / Rmax) - 0.45), 2.0)) - 0.0018029383826828) *
                     cos(9.0 * theta)) *
                (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                          (temp1 * temp) +
                      ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                     2.0) +
                 sin_theta_pow2 / (temp * temp)) *
                (20.0 * pow(tanh(20.0 * (r / Rmax) - 14.0), 2.0) - 20.0) * exp(-tanh(20.0 * (r / Rmax) - 14.0)) /
                sqrt((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                                 (temp1 * temp) +
                             ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                                (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                     factor_xi / (temp1 * temp) +
                                 ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                            sin_theta * cos_theta / (temp * temp)),
                           2.0)) +
                     (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                               (temp1 * temp) +
                           ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                          2.0) +
                      sin_theta_pow2 / (temp * temp)) *
                         (pow((ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                   factor_xi / (temp1 * temp) +
                               ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))),
                              2.0) +
                          cos_theta_pow2 / (temp * temp))) +
            (r / Rmax) *
                (((-6.67647559073009e-15) * (r / Rmax) +
                  (6000.0 - 6666.66666666667 * (r / Rmax)) * exp((-3333.33333333333) * pow(((r / Rmax) - 0.9), 2.0))) *
                     cos(21.0 * theta) +
                 (0.00368546744445083 * (r / Rmax) +
                  (45.0 - 100.0 * (r / Rmax)) * exp((-50.0) * pow(((r / Rmax) - 0.45), 2.0)) - 0.0018029383826828) *
                     cos(9.0 * theta)) *
                (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                          (temp1 * temp) +
                      ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                     2.0) +
                 sin_theta_pow2 / (temp * temp)) *
                (1.0 / 2.0 *
                     (((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                           (temp1 * temp) +
                       ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                          (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                               factor_xi / (temp1 * temp) +
                           ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                      sin_theta * cos_theta / (temp * temp)) *
                     (4.0 * inverse_aspect_ratio_epsilon * sin_theta * cos_theta_pow2 / pow((temp * temp), 2.0) +
                      2.0 *
                          ((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                               (temp1 * temp) +
                           ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                          ((-ellipticity_e) * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                               (r / Rmax) * sin_theta * cos_theta_pow2 * factor_xi /
                               (temp1 * pow((temp * temp), (3.0 / 2.0))) +
                           2.0 * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                               (r / Rmax) * sin_theta * cos_theta_pow2 * factor_xi / (temp2 * (temp * temp)) +
                           2.0 * ellipticity_e * inverse_aspect_ratio_epsilon * sin_theta * cos_theta * factor_xi /
                               (temp1 * temp)) +
                      2.0 *
                          (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                               factor_xi / (temp1 * temp) +
                           ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) *
                          (ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) * (r / Rmax) *
                               sin_theta_pow2 * cos_theta * factor_xi / (temp1 * pow((temp * temp), (3.0 / 2.0))) -
                           2.0 * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                               (r / Rmax) * sin_theta_pow2 * cos_theta * factor_xi / (temp2 * (temp * temp)) -
                           ellipticity_e * inverse_aspect_ratio_epsilon * sin_theta_pow2 * factor_xi / (temp1 * temp) +
                           ellipticity_e * inverse_aspect_ratio_epsilon * cos_theta_pow2 * factor_xi /
                               (temp1 * temp))) -
                 1.0 / 2.0 *
                     ((-2.0) * inverse_aspect_ratio_epsilon * pow(cos_theta, 3.0) / pow((temp * temp), 2.0) +
                      (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta * factor_xi /
                           (temp1 * temp) +
                       ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) *
                          ((-2.0) * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                               (r / Rmax) * sin_theta * cos_theta_pow2 * factor_xi /
                               (temp1 * pow((temp * temp), (3.0 / 2.0))) +
                           4.0 * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                               (r / Rmax) * sin_theta * cos_theta_pow2 * factor_xi / (temp2 * (temp * temp)) +
                           4.0 * ellipticity_e * inverse_aspect_ratio_epsilon * sin_theta * cos_theta * factor_xi /
                               (temp1 * temp))) *
                     (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                               (temp1 * temp) +
                           ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                          2.0) +
                      sin_theta_pow2 / (temp * temp)) -
                 1.0 / 2.0 *
                     ((-2.0) * inverse_aspect_ratio_epsilon * sin_theta_pow2 * cos_theta / pow((temp * temp), 2.0) +
                      ((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                           (temp1 * temp) +
                       ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                          (2.0 * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                               (r / Rmax) * sin_theta_pow2 * cos_theta * factor_xi /
                               (temp1 * pow((temp * temp), (3.0 / 2.0))) -
                           4.0 * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                               (r / Rmax) * sin_theta_pow2 * cos_theta * factor_xi / (temp2 * (temp * temp)) -
                           2.0 * ellipticity_e * inverse_aspect_ratio_epsilon * sin_theta_pow2 * factor_xi /
                               (temp1 * temp) +
                           2.0 * ellipticity_e * inverse_aspect_ratio_epsilon * cos_theta_pow2 * factor_xi /
                               (temp1 * temp))) *
                     (pow((ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                               factor_xi / (temp1 * temp) +
                           ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))),
                          2.0) +
                      cos_theta_pow2 / (temp * temp))) *
                exp(-tanh(20.0 * (r / Rmax) - 14.0)) /
                pow(((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                                 (temp1 * temp) +
                             ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                                (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                     factor_xi / (temp1 * temp) +
                                 ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                            sin_theta * cos_theta / (temp * temp)),
                           2.0)) +
                     (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                               (temp1 * temp) +
                           ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                          2.0) +
                      sin_theta_pow2 / (temp * temp)) *
                         (pow((ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                   factor_xi / (temp1 * temp) +
                               ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))),
                              2.0) +
                          cos_theta_pow2 / (temp * temp))),
                    (3.0 / 2.0)) +
            (r / Rmax) *
                ((10000.0 * pow((0.45 - (r / Rmax)), 2.0) * exp((-50.0) * pow(((r / Rmax) - 0.45), 2.0)) +
                  0.00368546744445083 - 100.0 * exp((-50.0) * pow(((r / Rmax) - 0.45), 2.0))) *
                     cos(9.0 * theta) +
                 (44444444.4444444 * pow((0.9 - (r / Rmax)), 2.0) *
                      exp((-3333.33333333333) * pow(((r / Rmax) - 0.9), 2.0)) -
                  6.67647559073009e-15 - 6666.66666666667 * exp((-3333.33333333333) * pow(((r / Rmax) - 0.9), 2.0))) *
                     cos(21.0 * theta)) *
                (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                          (temp1 * temp) +
                      ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                     2.0) +
                 sin_theta_pow2 / (temp * temp)) *
                exp(-tanh(20.0 * (r / Rmax) - 14.0)) /
                sqrt((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                                 (temp1 * temp) +
                             ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                                (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                     factor_xi / (temp1 * temp) +
                                 ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                            sin_theta * cos_theta / (temp * temp)),
                           2.0)) +
                     (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                               (temp1 * temp) +
                           ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                          2.0) +
                      sin_theta_pow2 / (temp * temp)) *
                         (pow((ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                   factor_xi / (temp1 * temp) +
                               ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))),
                              2.0) +
                          cos_theta_pow2 / (temp * temp))) -
            ((-21.0) *
                 ((-6.67647559073009e-15) * (r / Rmax) +
                  (6000.0 - 6666.66666666667 * (r / Rmax)) * exp((-3333.33333333333) * pow(((r / Rmax) - 0.9), 2.0))) *
                 sin(21.0 * theta) -
             9.0 *
                 (0.00368546744445083 * (r / Rmax) +
                  (45.0 - 100.0 * (r / Rmax)) * exp((-50.0) * pow(((r / Rmax) - 0.45), 2.0)) - 0.0018029383826828) *
                 sin(9.0 * theta)) *
                (((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                      (temp1 * temp) +
                  ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                     (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta * factor_xi /
                          (temp1 * temp) +
                      ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                 sin_theta * cos_theta / (temp * temp)) *
                exp(-tanh(20.0 * (r / Rmax) - 14.0)) /
                sqrt((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                                 (temp1 * temp) +
                             ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                                (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                     factor_xi / (temp1 * temp) +
                                 ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                            sin_theta * cos_theta / (temp * temp)),
                           2.0)) +
                     (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                               (temp1 * temp) +
                           ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                          2.0) +
                      sin_theta_pow2 / (temp * temp)) *
                         (pow((ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                   factor_xi / (temp1 * temp) +
                               ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))),
                              2.0) +
                          cos_theta_pow2 / (temp * temp))) -
            (((-6.67647559073009e-15) * (r / Rmax) +
              (6000.0 - 6666.66666666667 * (r / Rmax)) * exp((-3333.33333333333) * pow(((r / Rmax) - 0.9), 2.0))) *
                 cos(21.0 * theta) +
             (0.00368546744445083 * (r / Rmax) +
              (45.0 - 100.0 * (r / Rmax)) * exp((-50.0) * pow(((r / Rmax) - 0.45), 2.0)) - 0.0018029383826828) *
                 cos(9.0 * theta)) *
                (((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                      (temp1 * temp) +
                  ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                     (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta * factor_xi /
                          (temp1 * temp) +
                      ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                 sin_theta * cos_theta / (temp * temp)) *
                (1.0 / 2.0 *
                     (((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                           (temp1 * temp) +
                       ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                          (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                               factor_xi / (temp1 * temp) +
                           ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                      sin_theta * cos_theta / (temp * temp)) *
                     ((-4.0) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * cos_theta /
                          pow((temp * temp), 2.0) +
                      2.0 *
                          ((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                               (temp1 * temp) +
                           ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                          (ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                               ((r / Rmax) * (r / Rmax)) * sin_theta_pow2 * cos_theta * factor_xi /
                               (temp1 * pow((temp * temp), (3.0 / 2.0))) -
                           2.0 * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                               ((r / Rmax) * (r / Rmax)) * sin_theta_pow2 * cos_theta * factor_xi /
                               (temp2 * (temp * temp)) -
                           2.0 * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                               factor_xi / (temp1 * temp) +
                           ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * cos_theta_pow2 * factor_xi /
                               (temp1 * temp) +
                           ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) +
                      2.0 *
                          (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                               factor_xi / (temp1 * temp) +
                           ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) *
                          ((-ellipticity_e) * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                               ((r / Rmax) * (r / Rmax)) * pow(sin_theta, 3.0) * factor_xi /
                               (temp1 * pow((temp * temp), (3.0 / 2.0))) +
                           2.0 * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                               ((r / Rmax) * (r / Rmax)) * pow(sin_theta, 3.0) * factor_xi / (temp2 * (temp * temp)) -
                           3.0 * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                               factor_xi / (temp1 * temp) -
                           ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) +
                      2.0 * sin_theta_pow2 / (temp * temp) - 2.0 * cos_theta_pow2 / (temp * temp)) -
                 1.0 / 2.0 *
                     (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                               (temp1 * temp) +
                           ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                          2.0) +
                      sin_theta_pow2 / (temp * temp)) *
                     (2.0 * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta_pow2 /
                          pow((temp * temp), 2.0) +
                      (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta * factor_xi /
                           (temp1 * temp) +
                       ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) *
                          (2.0 * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                               ((r / Rmax) * (r / Rmax)) * sin_theta_pow2 * cos_theta * factor_xi /
                               (temp1 * pow((temp * temp), (3.0 / 2.0))) -
                           4.0 * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                               ((r / Rmax) * (r / Rmax)) * sin_theta_pow2 * cos_theta * factor_xi /
                               (temp2 * (temp * temp)) -
                           4.0 * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                               factor_xi / (temp1 * temp) +
                           2.0 * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * cos_theta_pow2 *
                               factor_xi / (temp1 * temp) +
                           2.0 * ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) -
                      2.0 * sin_theta * cos_theta / (temp * temp)) -
                 1.0 / 2.0 *
                     (pow((ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                               factor_xi / (temp1 * temp) +
                           ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))),
                          2.0) +
                      cos_theta_pow2 / (temp * temp)) *
                     (2.0 * inverse_aspect_ratio_epsilon * (r / Rmax) * pow(sin_theta, 3.0) / pow((temp * temp), 2.0) +
                      ((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                           (temp1 * temp) +
                       ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                          ((-2.0) * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                               ((r / Rmax) * (r / Rmax)) * pow(sin_theta, 3.0) * factor_xi /
                               (temp1 * pow((temp * temp), (3.0 / 2.0))) +
                           4.0 * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                               ((r / Rmax) * (r / Rmax)) * pow(sin_theta, 3.0) * factor_xi / (temp2 * (temp * temp)) -
                           6.0 * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                               factor_xi / (temp1 * temp) -
                           2.0 * ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) +
                      2.0 * sin_theta * cos_theta / (temp * temp))) *
                exp(-tanh(20.0 * (r / Rmax) - 14.0)) /
                pow(((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                                 (temp1 * temp) +
                             ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                                (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                     factor_xi / (temp1 * temp) +
                                 ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                            sin_theta * cos_theta / (temp * temp)),
                           2.0)) +
                     (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                               (temp1 * temp) +
                           ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                          2.0) +
                      sin_theta_pow2 / (temp * temp)) *
                         (pow((ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                   factor_xi / (temp1 * temp) +
                               ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))),
                              2.0) +
                          cos_theta_pow2 / (temp * temp))),
                    (3.0 / 2.0)) +
            (((-6.67647559073009e-15) * (r / Rmax) +
              (6000.0 - 6666.66666666667 * (r / Rmax)) * exp((-3333.33333333333) * pow(((r / Rmax) - 0.9), 2.0))) *
                 cos(21.0 * theta) +
             (0.00368546744445083 * (r / Rmax) +
              (45.0 - 100.0 * (r / Rmax)) * exp((-50.0) * pow(((r / Rmax) - 0.45), 2.0)) - 0.0018029383826828) *
                 cos(9.0 * theta)) *
                (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                          (temp1 * temp) +
                      ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                     2.0) +
                 sin_theta_pow2 / (temp * temp)) *
                exp(-tanh(20.0 * (r / Rmax) - 14.0)) /
                sqrt((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                                 (temp1 * temp) +
                             ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                                (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                     factor_xi / (temp1 * temp) +
                                 ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                            sin_theta * cos_theta / (temp * temp)),
                           2.0)) +
                     (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                               (temp1 * temp) +
                           ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                          2.0) +
                      sin_theta_pow2 / (temp * temp)) *
                         (pow((ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                   factor_xi / (temp1 * temp) +
                               ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))),
                              2.0) +
                          cos_theta_pow2 / (temp * temp))) -
            (((-6.67647559073009e-15) * (r / Rmax) +
              (6000.0 - 6666.66666666667 * (r / Rmax)) * exp((-3333.33333333333) * pow(((r / Rmax) - 0.9), 2.0))) *
                 cos(21.0 * theta) +
             (0.00368546744445083 * (r / Rmax) +
              (45.0 - 100.0 * (r / Rmax)) * exp((-50.0) * pow(((r / Rmax) - 0.45), 2.0)) - 0.0018029383826828) *
                 cos(9.0 * theta)) *
                ((-2.0) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * cos_theta /
                     pow((temp * temp), 2.0) +
                 ((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                      (temp1 * temp) +
                  ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                     (ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                          ((r / Rmax) * (r / Rmax)) * sin_theta_pow2 * cos_theta * factor_xi /
                          (temp1 * pow((temp * temp), (3.0 / 2.0))) -
                      2.0 * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                          ((r / Rmax) * (r / Rmax)) * sin_theta_pow2 * cos_theta * factor_xi / (temp2 * (temp * temp)) -
                      2.0 * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                          (temp1 * temp) +
                      ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * cos_theta_pow2 * factor_xi /
                          (temp1 * temp) +
                      ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) +
                 (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta * factor_xi /
                      (temp1 * temp) +
                  ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) *
                     ((-ellipticity_e) * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                          ((r / Rmax) * (r / Rmax)) * pow(sin_theta, 3.0) * factor_xi /
                          (temp1 * pow((temp * temp), (3.0 / 2.0))) +
                      2.0 * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                          ((r / Rmax) * (r / Rmax)) * pow(sin_theta, 3.0) * factor_xi / (temp2 * (temp * temp)) -
                      3.0 * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                          factor_xi / (temp1 * temp) -
                      ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) +
                 sin_theta_pow2 / (temp * temp) - cos_theta_pow2 / (temp * temp)) *
                exp(-tanh(20.0 * (r / Rmax) - 14.0)) /
                sqrt((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                                 (temp1 * temp) +
                             ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                                (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                     factor_xi / (temp1 * temp) +
                                 ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                            sin_theta * cos_theta / (temp * temp)),
                           2.0)) +
                     (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                               (temp1 * temp) +
                           ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                          2.0) +
                      sin_theta_pow2 / (temp * temp)) *
                         (pow((ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                   factor_xi / (temp1 * temp) +
                               ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))),
                              2.0) +
                          cos_theta_pow2 / (temp * temp))) -
            ((1.40205987405332e-13 * (r / Rmax) - 21.0 * (6000.0 - 6666.66666666667 * (r / Rmax)) *
                                                      exp((-3333.33333333333) * pow(((r / Rmax) - 0.9), 2.0))) *
                 sin(21.0 * theta) +
             ((-0.0331692070000574) * (r / Rmax) -
              9.0 * (45.0 - 100.0 * (r / Rmax)) * exp((-50.0) * pow(((r / Rmax) - 0.45), 2.0)) + 0.0162264454441452) *
                 sin(9.0 * theta)) *
                (((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                      (temp1 * temp) +
                  ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                     (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta * factor_xi /
                          (temp1 * temp) +
                      ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                 sin_theta * cos_theta / (temp * temp)) *
                exp(-tanh(20.0 * (r / Rmax) - 14.0)) /
                sqrt((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                                 (temp1 * temp) +
                             ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                                (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                     factor_xi / (temp1 * temp) +
                                 ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                            sin_theta * cos_theta / (temp * temp)),
                           2.0)) +
                     (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                               (temp1 * temp) +
                           ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                          2.0) +
                      sin_theta_pow2 / (temp * temp)) *
                         (pow((ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                   factor_xi / (temp1 * temp) +
                               ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))),
                              2.0) +
                          cos_theta_pow2 / (temp * temp))) -
            ((-21.0) *
                 ((-3.33823779536505e-15) * ((r / Rmax) * (r / Rmax)) +
                  exp((-3333.33333333333) * pow(((r / Rmax) - 0.9), 2.0))) *
                 sin(21.0 * theta) -
             9.0 *
                 (0.00184273372222541 * ((r / Rmax) * (r / Rmax)) - 0.0018029383826828 * (r / Rmax) -
                  4.00652973929511e-05 + exp((-50.0) * pow(((r / Rmax) - 0.45), 2.0))) *
                 sin(9.0 * theta)) *
                (((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                      (temp1 * temp) +
                  ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                     (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta * factor_xi /
                          (temp1 * temp) +
                      ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                 sin_theta * cos_theta / (temp * temp)) *
                (20.0 * pow(tanh(20.0 * (r / Rmax) - 14.0), 2.0) - 20.0) * exp(-tanh(20.0 * (r / Rmax) - 14.0)) /
                sqrt((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                                 (temp1 * temp) +
                             ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                                (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                     factor_xi / (temp1 * temp) +
                                 ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                            sin_theta * cos_theta / (temp * temp)),
                           2.0)) +
                     (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                               (temp1 * temp) +
                           ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                          2.0) +
                      sin_theta_pow2 / (temp * temp)) *
                         (pow((ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                   factor_xi / (temp1 * temp) +
                               ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))),
                              2.0) +
                          cos_theta_pow2 / (temp * temp))) -
            ((-21.0) *
                 ((-3.33823779536505e-15) * ((r / Rmax) * (r / Rmax)) +
                  exp((-3333.33333333333) * pow(((r / Rmax) - 0.9), 2.0))) *
                 sin(21.0 * theta) -
             9.0 *
                 (0.00184273372222541 * ((r / Rmax) * (r / Rmax)) - 0.0018029383826828 * (r / Rmax) -
                  4.00652973929511e-05 + exp((-50.0) * pow(((r / Rmax) - 0.45), 2.0))) *
                 sin(9.0 * theta)) *
                (((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                      (temp1 * temp) +
                  ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                     (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta * factor_xi /
                          (temp1 * temp) +
                      ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                 sin_theta * cos_theta / (temp * temp)) *
                (1.0 / 2.0 *
                     (((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                           (temp1 * temp) +
                       ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                          (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                               factor_xi / (temp1 * temp) +
                           ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                      sin_theta * cos_theta / (temp * temp)) *
                     (4.0 * inverse_aspect_ratio_epsilon * sin_theta * cos_theta_pow2 / pow((temp * temp), 2.0) +
                      2.0 *
                          ((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                               (temp1 * temp) +
                           ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                          ((-ellipticity_e) * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                               (r / Rmax) * sin_theta * cos_theta_pow2 * factor_xi /
                               (temp1 * pow((temp * temp), (3.0 / 2.0))) +
                           2.0 * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                               (r / Rmax) * sin_theta * cos_theta_pow2 * factor_xi / (temp2 * (temp * temp)) +
                           2.0 * ellipticity_e * inverse_aspect_ratio_epsilon * sin_theta * cos_theta * factor_xi /
                               (temp1 * temp)) +
                      2.0 *
                          (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                               factor_xi / (temp1 * temp) +
                           ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) *
                          (ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) * (r / Rmax) *
                               sin_theta_pow2 * cos_theta * factor_xi / (temp1 * pow((temp * temp), (3.0 / 2.0))) -
                           2.0 * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                               (r / Rmax) * sin_theta_pow2 * cos_theta * factor_xi / (temp2 * (temp * temp)) -
                           ellipticity_e * inverse_aspect_ratio_epsilon * sin_theta_pow2 * factor_xi / (temp1 * temp) +
                           ellipticity_e * inverse_aspect_ratio_epsilon * cos_theta_pow2 * factor_xi /
                               (temp1 * temp))) -
                 1.0 / 2.0 *
                     ((-2.0) * inverse_aspect_ratio_epsilon * pow(cos_theta, 3.0) / pow((temp * temp), 2.0) +
                      (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta * factor_xi /
                           (temp1 * temp) +
                       ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) *
                          ((-2.0) * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                               (r / Rmax) * sin_theta * cos_theta_pow2 * factor_xi /
                               (temp1 * pow((temp * temp), (3.0 / 2.0))) +
                           4.0 * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                               (r / Rmax) * sin_theta * cos_theta_pow2 * factor_xi / (temp2 * (temp * temp)) +
                           4.0 * ellipticity_e * inverse_aspect_ratio_epsilon * sin_theta * cos_theta * factor_xi /
                               (temp1 * temp))) *
                     (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                               (temp1 * temp) +
                           ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                          2.0) +
                      sin_theta_pow2 / (temp * temp)) -
                 1.0 / 2.0 *
                     ((-2.0) * inverse_aspect_ratio_epsilon * sin_theta_pow2 * cos_theta / pow((temp * temp), 2.0) +
                      ((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                           (temp1 * temp) +
                       ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                          (2.0 * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                               (r / Rmax) * sin_theta_pow2 * cos_theta * factor_xi /
                               (temp1 * pow((temp * temp), (3.0 / 2.0))) -
                           4.0 * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                               (r / Rmax) * sin_theta_pow2 * cos_theta * factor_xi / (temp2 * (temp * temp)) -
                           2.0 * ellipticity_e * inverse_aspect_ratio_epsilon * sin_theta_pow2 * factor_xi /
                               (temp1 * temp) +
                           2.0 * ellipticity_e * inverse_aspect_ratio_epsilon * cos_theta_pow2 * factor_xi /
                               (temp1 * temp))) *
                     (pow((ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                               factor_xi / (temp1 * temp) +
                           ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))),
                          2.0) +
                      cos_theta_pow2 / (temp * temp))) *
                exp(-tanh(20.0 * (r / Rmax) - 14.0)) /
                pow(((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                                 (temp1 * temp) +
                             ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                                (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                     factor_xi / (temp1 * temp) +
                                 ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                            sin_theta * cos_theta / (temp * temp)),
                           2.0)) +
                     (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                               (temp1 * temp) +
                           ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                          2.0) +
                      sin_theta_pow2 / (temp * temp)) *
                         (pow((ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                   factor_xi / (temp1 * temp) +
                               ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))),
                              2.0) +
                          cos_theta_pow2 / (temp * temp))),
                    (3.0 / 2.0)) -
            ((-21.0) *
                 ((-3.33823779536505e-15) * ((r / Rmax) * (r / Rmax)) +
                  exp((-3333.33333333333) * pow(((r / Rmax) - 0.9), 2.0))) *
                 sin(21.0 * theta) -
             9.0 *
                 (0.00184273372222541 * ((r / Rmax) * (r / Rmax)) - 0.0018029383826828 * (r / Rmax) -
                  4.00652973929511e-05 + exp((-50.0) * pow(((r / Rmax) - 0.45), 2.0))) *
                 sin(9.0 * theta)) *
                (2.0 * inverse_aspect_ratio_epsilon * sin_theta * cos_theta_pow2 / pow((temp * temp), 2.0) +
                 ((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                      (temp1 * temp) +
                  ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                     ((-ellipticity_e) * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) * (r / Rmax) *
                          sin_theta * cos_theta_pow2 * factor_xi / (temp1 * pow((temp * temp), (3.0 / 2.0))) +
                      2.0 * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) * (r / Rmax) *
                          sin_theta * cos_theta_pow2 * factor_xi / (temp2 * (temp * temp)) +
                      2.0 * ellipticity_e * inverse_aspect_ratio_epsilon * sin_theta * cos_theta * factor_xi /
                          (temp1 * temp)) +
                 (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta * factor_xi /
                      (temp1 * temp) +
                  ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) *
                     (ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) * (r / Rmax) *
                          sin_theta_pow2 * cos_theta * factor_xi / (temp1 * pow((temp * temp), (3.0 / 2.0))) -
                      2.0 * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) * (r / Rmax) *
                          sin_theta_pow2 * cos_theta * factor_xi / (temp2 * (temp * temp)) -
                      ellipticity_e * inverse_aspect_ratio_epsilon * sin_theta_pow2 * factor_xi / (temp1 * temp) +
                      ellipticity_e * inverse_aspect_ratio_epsilon * cos_theta_pow2 * factor_xi / (temp1 * temp))) *
                exp(-tanh(20.0 * (r / Rmax) - 14.0)) /
                sqrt((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                                 (temp1 * temp) +
                             ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                                (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                     factor_xi / (temp1 * temp) +
                                 ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                            sin_theta * cos_theta / (temp * temp)),
                           2.0)) +
                     (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                               (temp1 * temp) +
                           ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                          2.0) +
                      sin_theta_pow2 / (temp * temp)) *
                         (pow((ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                   factor_xi / (temp1 * temp) +
                               ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))),
                              2.0) +
                          cos_theta_pow2 / (temp * temp))) +
            ((-21.0) *
                 ((-3.33823779536505e-15) * ((r / Rmax) * (r / Rmax)) +
                  exp((-3333.33333333333) * pow(((r / Rmax) - 0.9), 2.0))) *
                 sin(21.0 * theta) -
             9.0 *
                 (0.00184273372222541 * ((r / Rmax) * (r / Rmax)) - 0.0018029383826828 * (r / Rmax) -
                  4.00652973929511e-05 + exp((-50.0) * pow(((r / Rmax) - 0.45), 2.0))) *
                 sin(9.0 * theta)) *
                (2.0 * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta_pow2 /
                     pow((temp * temp), 2.0) +
                 (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta * factor_xi /
                      (temp1 * temp) +
                  ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) *
                     (2.0 * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                          ((r / Rmax) * (r / Rmax)) * sin_theta_pow2 * cos_theta * factor_xi /
                          (temp1 * pow((temp * temp), (3.0 / 2.0))) -
                      4.0 * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                          ((r / Rmax) * (r / Rmax)) * sin_theta_pow2 * cos_theta * factor_xi / (temp2 * (temp * temp)) -
                      4.0 * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                          (temp1 * temp) +
                      2.0 * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * cos_theta_pow2 * factor_xi /
                          (temp1 * temp) +
                      2.0 * ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) -
                 2.0 * sin_theta * cos_theta / (temp * temp)) *
                exp(-tanh(20.0 * (r / Rmax) - 14.0)) /
                ((r / Rmax) * sqrt((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) *
                                               sin_theta_pow2 * factor_xi / (temp1 * temp) +
                                           ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                                              (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta *
                                                   cos_theta * factor_xi / (temp1 * temp) +
                                               ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                                          sin_theta * cos_theta / (temp * temp)),
                                         2.0)) +
                                   (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                             factor_xi / (temp1 * temp) +
                                         ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                                        2.0) +
                                    sin_theta_pow2 / (temp * temp)) *
                                       (pow((ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta *
                                                 cos_theta * factor_xi / (temp1 * temp) +
                                             ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))),
                                            2.0) +
                                        cos_theta_pow2 / (temp * temp)))) +
            ((-21.0) *
                 ((-3.33823779536505e-15) * ((r / Rmax) * (r / Rmax)) +
                  exp((-3333.33333333333) * pow(((r / Rmax) - 0.9), 2.0))) *
                 sin(21.0 * theta) -
             9.0 *
                 (0.00184273372222541 * ((r / Rmax) * (r / Rmax)) - 0.0018029383826828 * (r / Rmax) -
                  4.00652973929511e-05 + exp((-50.0) * pow(((r / Rmax) - 0.45), 2.0))) *
                 sin(9.0 * theta)) *
                (pow((ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta * factor_xi /
                          (temp1 * temp) +
                      ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))),
                     2.0) +
                 cos_theta_pow2 / (temp * temp)) *
                (1.0 / 2.0 *
                     (((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                           (temp1 * temp) +
                       ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                          (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                               factor_xi / (temp1 * temp) +
                           ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                      sin_theta * cos_theta / (temp * temp)) *
                     ((-4.0) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * cos_theta /
                          pow((temp * temp), 2.0) +
                      2.0 *
                          ((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                               (temp1 * temp) +
                           ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                          (ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                               ((r / Rmax) * (r / Rmax)) * sin_theta_pow2 * cos_theta * factor_xi /
                               (temp1 * pow((temp * temp), (3.0 / 2.0))) -
                           2.0 * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                               ((r / Rmax) * (r / Rmax)) * sin_theta_pow2 * cos_theta * factor_xi /
                               (temp2 * (temp * temp)) -
                           2.0 * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                               factor_xi / (temp1 * temp) +
                           ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * cos_theta_pow2 * factor_xi /
                               (temp1 * temp) +
                           ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) +
                      2.0 *
                          (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                               factor_xi / (temp1 * temp) +
                           ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) *
                          ((-ellipticity_e) * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                               ((r / Rmax) * (r / Rmax)) * pow(sin_theta, 3.0) * factor_xi /
                               (temp1 * pow((temp * temp), (3.0 / 2.0))) +
                           2.0 * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                               ((r / Rmax) * (r / Rmax)) * pow(sin_theta, 3.0) * factor_xi / (temp2 * (temp * temp)) -
                           3.0 * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                               factor_xi / (temp1 * temp) -
                           ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) +
                      2.0 * sin_theta_pow2 / (temp * temp) - 2.0 * cos_theta_pow2 / (temp * temp)) -
                 1.0 / 2.0 *
                     (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                               (temp1 * temp) +
                           ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                          2.0) +
                      sin_theta_pow2 / (temp * temp)) *
                     (2.0 * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta_pow2 /
                          pow((temp * temp), 2.0) +
                      (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta * factor_xi /
                           (temp1 * temp) +
                       ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) *
                          (2.0 * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                               ((r / Rmax) * (r / Rmax)) * sin_theta_pow2 * cos_theta * factor_xi /
                               (temp1 * pow((temp * temp), (3.0 / 2.0))) -
                           4.0 * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                               ((r / Rmax) * (r / Rmax)) * sin_theta_pow2 * cos_theta * factor_xi /
                               (temp2 * (temp * temp)) -
                           4.0 * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                               factor_xi / (temp1 * temp) +
                           2.0 * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * cos_theta_pow2 *
                               factor_xi / (temp1 * temp) +
                           2.0 * ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) -
                      2.0 * sin_theta * cos_theta / (temp * temp)) -
                 1.0 / 2.0 *
                     (pow((ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                               factor_xi / (temp1 * temp) +
                           ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))),
                          2.0) +
                      cos_theta_pow2 / (temp * temp)) *
                     (2.0 * inverse_aspect_ratio_epsilon * (r / Rmax) * pow(sin_theta, 3.0) / pow((temp * temp), 2.0) +
                      ((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                           (temp1 * temp) +
                       ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                          ((-2.0) * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                               ((r / Rmax) * (r / Rmax)) * pow(sin_theta, 3.0) * factor_xi /
                               (temp1 * pow((temp * temp), (3.0 / 2.0))) +
                           4.0 * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                               ((r / Rmax) * (r / Rmax)) * pow(sin_theta, 3.0) * factor_xi / (temp2 * (temp * temp)) -
                           6.0 * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                               factor_xi / (temp1 * temp) -
                           2.0 * ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) +
                      2.0 * sin_theta * cos_theta / (temp * temp))) *
                exp(-tanh(20.0 * (r / Rmax) - 14.0)) /
                ((r / Rmax) * pow(((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) *
                                               sin_theta_pow2 * factor_xi / (temp1 * temp) +
                                           ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                                              (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta *
                                                   cos_theta * factor_xi / (temp1 * temp) +
                                               ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                                          sin_theta * cos_theta / (temp * temp)),
                                         2.0)) +
                                   (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                             factor_xi / (temp1 * temp) +
                                         ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                                        2.0) +
                                    sin_theta_pow2 / (temp * temp)) *
                                       (pow((ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta *
                                                 cos_theta * factor_xi / (temp1 * temp) +
                                             ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))),
                                            2.0) +
                                        cos_theta_pow2 / (temp * temp))),
                                  (3.0 / 2.0))) +
            (21.0 *
                 (7.0102993702666e-14 * ((r / Rmax) * (r / Rmax)) -
                  21.0 * exp((-3333.33333333333) * pow(((r / Rmax) - 0.9), 2.0))) *
                 cos(21.0 * theta) +
             9.0 *
                 ((-0.0165846035000287) * ((r / Rmax) * (r / Rmax)) + 0.0162264454441452 * (r / Rmax) +
                  0.00036058767653656 - 9.0 * exp((-50.0) * pow(((r / Rmax) - 0.45), 2.0))) *
                 cos(9.0 * theta)) *
                (pow((ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta * factor_xi /
                          (temp1 * temp) +
                      ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))),
                     2.0) +
                 cos_theta_pow2 / (temp * temp)) *
                exp(-tanh(20.0 * (r / Rmax) - 14.0)) /
                ((r / Rmax) * sqrt((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) *
                                               sin_theta_pow2 * factor_xi / (temp1 * temp) +
                                           ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                                              (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta *
                                                   cos_theta * factor_xi / (temp1 * temp) +
                                               ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                                          sin_theta * cos_theta / (temp * temp)),
                                         2.0)) +
                                   (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                             factor_xi / (temp1 * temp) +
                                         ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                                        2.0) +
                                    sin_theta_pow2 / (temp * temp)) *
                                       (pow((ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta *
                                                 cos_theta * factor_xi / (temp1 * temp) +
                                             ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))),
                                            2.0) +
                                        cos_theta_pow2 / (temp * temp))))) /
               ((r / Rmax) * sqrt((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) *
                                              sin_theta_pow2 * factor_xi / (temp1 * temp) +
                                          ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                                             (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta *
                                                  cos_theta * factor_xi / (temp1 * temp) +
                                              ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                                         sin_theta * cos_theta / (temp * temp)),
                                        2.0)) +
                                  (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                            factor_xi / (temp1 * temp) +
                                        ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                                       2.0) +
                                   sin_theta_pow2 / (temp * temp)) *
                                      (pow((ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta *
                                                cos_theta * factor_xi / (temp1 * temp) +
                                            ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))),
                                           2.0) +
                                       cos_theta_pow2 / (temp * temp))));
}