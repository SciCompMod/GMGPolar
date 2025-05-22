#include "../include/InputFunctions/SourceTerms/polarR6_ZoniShiftedGyro_CzarnyGeometry.h"

void PolarR6_ZoniShiftedGyro_CzarnyGeometry::initializeGeometry()
{
    factor_xi = 1.0 / sqrt(1.0 - inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon / 4.0);
}

PolarR6_ZoniShiftedGyro_CzarnyGeometry::PolarR6_ZoniShiftedGyro_CzarnyGeometry(
    const double& Rmax, const double& inverse_aspect_ratio_epsilon, const double& ellipticity_e)
    : Rmax(Rmax)
    , inverse_aspect_ratio_epsilon(inverse_aspect_ratio_epsilon)
    , ellipticity_e(ellipticity_e)
{
    initializeGeometry();
}

double PolarR6_ZoniShiftedGyro_CzarnyGeometry::rhs_f(const double& r, const double& theta, const double& sin_theta,
                                                     const double& cos_theta) const
{
    // return 1.0; // use for memcheck
    double temp =
        sqrt(inverse_aspect_ratio_epsilon * (inverse_aspect_ratio_epsilon + 2.0 * (r / Rmax) * cos_theta) + 1.0);
    double sin_theta_pow2 = pow(sin_theta, 2.0);
    double cos_theta_pow2 = pow(cos_theta, 2.0);
    double temp1          = pow((2.0 - temp), 2.0);
    double temp2          = pow((2.0 - temp), 3.0);

    return 0.4096 * pow((r / Rmax), 6.0) * pow(((r / Rmax) - 1.0), 6.0) * exp(tanh(20.0 * (r / Rmax) - 14.0)) *
               cos(11.0 * theta) -
           pow((r / Rmax), 4.0) *
               (4.5056 * (r / Rmax) * pow(((r / Rmax) - 1.0), 6.0) *
                    (((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                          (temp1 * temp) +
                      ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                         (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                              factor_xi / (temp1 * temp) +
                          ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                     sin_theta * cos_theta / (temp * temp)) *
                    (20.0 * pow(tanh(20.0 * (r / Rmax) - 14.0), 2.0) - 20.0) * exp(-tanh(20.0 * (r / Rmax) - 14.0)) *
                    sin(11.0 * theta) /
                    sqrt((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                     factor_xi / (temp1 * temp) +
                                 ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                                    (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                         factor_xi / (temp1 * temp) +
                                     ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                                sin_theta * cos_theta / (temp * temp)),
                               2.0)) +
                         (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                   factor_xi / (temp1 * temp) +
                               ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                              2.0) +
                          sin_theta_pow2 / (temp * temp)) *
                             (pow((ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                       factor_xi / (temp1 * temp) +
                                   ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))),
                                  2.0) +
                              cos_theta_pow2 / (temp * temp))) +
                4.5056 * (r / Rmax) * pow(((r / Rmax) - 1.0), 6.0) *
                    (((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                          (temp1 * temp) +
                      ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                         (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                              factor_xi / (temp1 * temp) +
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
                              ((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                   factor_xi / (temp1 * temp) +
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
                              (ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                                   (r / Rmax) * sin_theta_pow2 * cos_theta * factor_xi /
                                   (temp1 * pow((temp * temp), (3.0 / 2.0))) -
                               2.0 * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                                   (r / Rmax) * sin_theta_pow2 * cos_theta * factor_xi / (temp2 * (temp * temp)) -
                               ellipticity_e * inverse_aspect_ratio_epsilon * sin_theta_pow2 * factor_xi /
                                   (temp1 * temp) +
                               ellipticity_e * inverse_aspect_ratio_epsilon * cos_theta_pow2 * factor_xi /
                                   (temp1 * temp))) -
                     1.0 / 2.0 *
                         ((-2.0) * inverse_aspect_ratio_epsilon * pow(cos_theta, 3.0) / pow((temp * temp), 2.0) +
                          (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                               factor_xi / (temp1 * temp) +
                           ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) *
                              ((-2.0) * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                                   (r / Rmax) * sin_theta * cos_theta_pow2 * factor_xi /
                                   (temp1 * pow((temp * temp), (3.0 / 2.0))) +
                               4.0 * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                                   (r / Rmax) * sin_theta * cos_theta_pow2 * factor_xi / (temp2 * (temp * temp)) +
                               4.0 * ellipticity_e * inverse_aspect_ratio_epsilon * sin_theta * cos_theta * factor_xi /
                                   (temp1 * temp))) *
                         (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                   factor_xi / (temp1 * temp) +
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
                    exp(-tanh(20.0 * (r / Rmax) - 14.0)) * sin(11.0 * theta) /
                    pow(((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                     factor_xi / (temp1 * temp) +
                                 ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                                    (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                         factor_xi / (temp1 * temp) +
                                     ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                                sin_theta * cos_theta / (temp * temp)),
                               2.0)) +
                         (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                   factor_xi / (temp1 * temp) +
                               ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                              2.0) +
                          sin_theta_pow2 / (temp * temp)) *
                             (pow((ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                       factor_xi / (temp1 * temp) +
                                   ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))),
                                  2.0) +
                              cos_theta_pow2 / (temp * temp))),
                        (3.0 / 2.0)) +
                4.5056 * (r / Rmax) * pow(((r / Rmax) - 1.0), 6.0) *
                    (2.0 * inverse_aspect_ratio_epsilon * sin_theta * cos_theta_pow2 / pow((temp * temp), 2.0) +
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
                     (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta * factor_xi /
                          (temp1 * temp) +
                      ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) *
                         (ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) * (r / Rmax) *
                              sin_theta_pow2 * cos_theta * factor_xi / (temp1 * pow((temp * temp), (3.0 / 2.0))) -
                          2.0 * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                              (r / Rmax) * sin_theta_pow2 * cos_theta * factor_xi / (temp2 * (temp * temp)) -
                          ellipticity_e * inverse_aspect_ratio_epsilon * sin_theta_pow2 * factor_xi / (temp1 * temp) +
                          ellipticity_e * inverse_aspect_ratio_epsilon * cos_theta_pow2 * factor_xi / (temp1 * temp))) *
                    exp(-tanh(20.0 * (r / Rmax) - 14.0)) * sin(11.0 * theta) /
                    sqrt((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                     factor_xi / (temp1 * temp) +
                                 ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                                    (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                         factor_xi / (temp1 * temp) +
                                     ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                                sin_theta * cos_theta / (temp * temp)),
                               2.0)) +
                         (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                   factor_xi / (temp1 * temp) +
                               ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                              2.0) +
                          sin_theta_pow2 / (temp * temp)) *
                             (pow((ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                       factor_xi / (temp1 * temp) +
                                   ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))),
                                  2.0) +
                              cos_theta_pow2 / (temp * temp))) +
                27.0336 * (r / Rmax) * pow(((r / Rmax) - 1.0), 5.0) *
                    (((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                          (temp1 * temp) +
                      ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                         (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                              factor_xi / (temp1 * temp) +
                          ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                     sin_theta * cos_theta / (temp * temp)) *
                    exp(-tanh(20.0 * (r / Rmax) - 14.0)) * sin(11.0 * theta) /
                    sqrt((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                     factor_xi / (temp1 * temp) +
                                 ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                                    (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                         factor_xi / (temp1 * temp) +
                                     ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                                sin_theta * cos_theta / (temp * temp)),
                               2.0)) +
                         (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                   factor_xi / (temp1 * temp) +
                               ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                              2.0) +
                          sin_theta_pow2 / (temp * temp)) *
                             (pow((ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                       factor_xi / (temp1 * temp) +
                                   ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))),
                                  2.0) +
                              cos_theta_pow2 / (temp * temp))) +
                (r / Rmax) *
                    (12.288 * (r / Rmax) * pow(((r / Rmax) - 1.0), 4.0) * cos(11.0 * theta) +
                     17.2032 * pow(((r / Rmax) - 1.0), 5.0) * cos(11.0 * theta)) *
                    (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                              (temp1 * temp) +
                          ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                         2.0) +
                     sin_theta_pow2 / (temp * temp)) *
                    exp(-tanh(20.0 * (r / Rmax) - 14.0)) /
                    sqrt((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                     factor_xi / (temp1 * temp) +
                                 ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                                    (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                         factor_xi / (temp1 * temp) +
                                     ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                                sin_theta * cos_theta / (temp * temp)),
                               2.0)) +
                         (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                   factor_xi / (temp1 * temp) +
                               ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                              2.0) +
                          sin_theta_pow2 / (temp * temp)) *
                             (pow((ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                       factor_xi / (temp1 * temp) +
                                   ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))),
                                  2.0) +
                              cos_theta_pow2 / (temp * temp))) +
                (r / Rmax) *
                    (2.4576 * (r / Rmax) * pow(((r / Rmax) - 1.0), 5.0) * cos(11.0 * theta) +
                     2.4576 * pow(((r / Rmax) - 1.0), 6.0) * cos(11.0 * theta)) *
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
                    exp(-tanh(20.0 * (r / Rmax) - 14.0)) /
                    sqrt((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                     factor_xi / (temp1 * temp) +
                                 ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                                    (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                         factor_xi / (temp1 * temp) +
                                     ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                                sin_theta * cos_theta / (temp * temp)),
                               2.0)) +
                         (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                   factor_xi / (temp1 * temp) +
                               ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                              2.0) +
                          sin_theta_pow2 / (temp * temp)) *
                             (pow((ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                       factor_xi / (temp1 * temp) +
                                   ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))),
                                  2.0) +
                              cos_theta_pow2 / (temp * temp))) +
                (r / Rmax) *
                    (2.4576 * (r / Rmax) * pow(((r / Rmax) - 1.0), 5.0) * cos(11.0 * theta) +
                     2.4576 * pow(((r / Rmax) - 1.0), 6.0) * cos(11.0 * theta)) *
                    (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                              (temp1 * temp) +
                          ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                         2.0) +
                     sin_theta_pow2 / (temp * temp)) *
                    (20.0 * pow(tanh(20.0 * (r / Rmax) - 14.0), 2.0) - 20.0) * exp(-tanh(20.0 * (r / Rmax) - 14.0)) /
                    sqrt((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                     factor_xi / (temp1 * temp) +
                                 ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                                    (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                         factor_xi / (temp1 * temp) +
                                     ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                                sin_theta * cos_theta / (temp * temp)),
                               2.0)) +
                         (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                   factor_xi / (temp1 * temp) +
                               ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                              2.0) +
                          sin_theta_pow2 / (temp * temp)) *
                             (pow((ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                       factor_xi / (temp1 * temp) +
                                   ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))),
                                  2.0) +
                              cos_theta_pow2 / (temp * temp))) +
                (r / Rmax) *
                    (2.4576 * (r / Rmax) * pow(((r / Rmax) - 1.0), 5.0) * cos(11.0 * theta) +
                     2.4576 * pow(((r / Rmax) - 1.0), 6.0) * cos(11.0 * theta)) *
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
                              ((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                   factor_xi / (temp1 * temp) +
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
                              (ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                                   (r / Rmax) * sin_theta_pow2 * cos_theta * factor_xi /
                                   (temp1 * pow((temp * temp), (3.0 / 2.0))) -
                               2.0 * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                                   (r / Rmax) * sin_theta_pow2 * cos_theta * factor_xi / (temp2 * (temp * temp)) -
                               ellipticity_e * inverse_aspect_ratio_epsilon * sin_theta_pow2 * factor_xi /
                                   (temp1 * temp) +
                               ellipticity_e * inverse_aspect_ratio_epsilon * cos_theta_pow2 * factor_xi /
                                   (temp1 * temp))) -
                     1.0 / 2.0 *
                         ((-2.0) * inverse_aspect_ratio_epsilon * pow(cos_theta, 3.0) / pow((temp * temp), 2.0) +
                          (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                               factor_xi / (temp1 * temp) +
                           ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) *
                              ((-2.0) * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                                   (r / Rmax) * sin_theta * cos_theta_pow2 * factor_xi /
                                   (temp1 * pow((temp * temp), (3.0 / 2.0))) +
                               4.0 * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                                   (r / Rmax) * sin_theta * cos_theta_pow2 * factor_xi / (temp2 * (temp * temp)) +
                               4.0 * ellipticity_e * inverse_aspect_ratio_epsilon * sin_theta * cos_theta * factor_xi /
                                   (temp1 * temp))) *
                         (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                   factor_xi / (temp1 * temp) +
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
                    pow(((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                     factor_xi / (temp1 * temp) +
                                 ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                                    (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                         factor_xi / (temp1 * temp) +
                                     ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                                sin_theta * cos_theta / (temp * temp)),
                               2.0)) +
                         (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                   factor_xi / (temp1 * temp) +
                               ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                              2.0) +
                          sin_theta_pow2 / (temp * temp)) *
                             (pow((ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                       factor_xi / (temp1 * temp) +
                                   ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))),
                                  2.0) +
                              cos_theta_pow2 / (temp * temp))),
                        (3.0 / 2.0)) +
                27.0336 * pow(((r / Rmax) - 1.0), 6.0) *
                    (((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                          (temp1 * temp) +
                      ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                         (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                              factor_xi / (temp1 * temp) +
                          ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                     sin_theta * cos_theta / (temp * temp)) *
                    exp(-tanh(20.0 * (r / Rmax) - 14.0)) * sin(11.0 * theta) /
                    sqrt((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                     factor_xi / (temp1 * temp) +
                                 ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                                    (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                         factor_xi / (temp1 * temp) +
                                     ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                                sin_theta * cos_theta / (temp * temp)),
                               2.0)) +
                         (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                   factor_xi / (temp1 * temp) +
                               ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                              2.0) +
                          sin_theta_pow2 / (temp * temp)) *
                             (pow((ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                       factor_xi / (temp1 * temp) +
                                   ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))),
                                  2.0) +
                              cos_theta_pow2 / (temp * temp))) -
                49.5616 * pow(((r / Rmax) - 1.0), 6.0) *
                    (pow((ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                              factor_xi / (temp1 * temp) +
                          ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))),
                         2.0) +
                     cos_theta_pow2 / (temp * temp)) *
                    exp(-tanh(20.0 * (r / Rmax) - 14.0)) * cos(11.0 * theta) /
                    sqrt((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                     factor_xi / (temp1 * temp) +
                                 ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                                    (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                         factor_xi / (temp1 * temp) +
                                     ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                                sin_theta * cos_theta / (temp * temp)),
                               2.0)) +
                         (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                   factor_xi / (temp1 * temp) +
                               ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                              2.0) +
                          sin_theta_pow2 / (temp * temp)) *
                             (pow((ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                       factor_xi / (temp1 * temp) +
                                   ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))),
                                  2.0) +
                              cos_theta_pow2 / (temp * temp))) -
                4.5056 * pow(((r / Rmax) - 1.0), 6.0) *
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
                          4.0 * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                              (temp1 * temp) +
                          2.0 * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * cos_theta_pow2 * factor_xi /
                              (temp1 * temp) +
                          2.0 * ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) -
                     2.0 * sin_theta * cos_theta / (temp * temp)) *
                    exp(-tanh(20.0 * (r / Rmax) - 14.0)) * sin(11.0 * theta) /
                    sqrt((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                     factor_xi / (temp1 * temp) +
                                 ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                                    (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                         factor_xi / (temp1 * temp) +
                                     ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                                sin_theta * cos_theta / (temp * temp)),
                               2.0)) +
                         (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                   factor_xi / (temp1 * temp) +
                               ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                              2.0) +
                          sin_theta_pow2 / (temp * temp)) *
                             (pow((ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                       factor_xi / (temp1 * temp) +
                                   ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))),
                                  2.0) +
                              cos_theta_pow2 / (temp * temp))) -
                4.5056 * pow(((r / Rmax) - 1.0), 6.0) *
                    (pow((ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                              factor_xi / (temp1 * temp) +
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
                              ((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                   factor_xi / (temp1 * temp) +
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
                                   ((r / Rmax) * (r / Rmax)) * pow(sin_theta, 3.0) * factor_xi /
                                   (temp2 * (temp * temp)) -
                               3.0 * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                   factor_xi / (temp1 * temp) -
                               ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) +
                          2.0 * sin_theta_pow2 / (temp * temp) - 2.0 * cos_theta_pow2 / (temp * temp)) -
                     1.0 / 2.0 *
                         (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                   factor_xi / (temp1 * temp) +
                               ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                              2.0) +
                          sin_theta_pow2 / (temp * temp)) *
                         (2.0 * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta_pow2 /
                              pow((temp * temp), 2.0) +
                          (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                               factor_xi / (temp1 * temp) +
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
                         (2.0 * inverse_aspect_ratio_epsilon * (r / Rmax) * pow(sin_theta, 3.0) /
                              pow((temp * temp), 2.0) +
                          ((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                               (temp1 * temp) +
                           ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                              ((-2.0) * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                                   ((r / Rmax) * (r / Rmax)) * pow(sin_theta, 3.0) * factor_xi /
                                   (temp1 * pow((temp * temp), (3.0 / 2.0))) +
                               4.0 * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                                   ((r / Rmax) * (r / Rmax)) * pow(sin_theta, 3.0) * factor_xi /
                                   (temp2 * (temp * temp)) -
                               6.0 * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                   factor_xi / (temp1 * temp) -
                               2.0 * ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) +
                          2.0 * sin_theta * cos_theta / (temp * temp))) *
                    exp(-tanh(20.0 * (r / Rmax) - 14.0)) * sin(11.0 * theta) /
                    pow(((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                     factor_xi / (temp1 * temp) +
                                 ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                                    (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                         factor_xi / (temp1 * temp) +
                                     ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                                sin_theta * cos_theta / (temp * temp)),
                               2.0)) +
                         (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                   factor_xi / (temp1 * temp) +
                               ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                              2.0) +
                          sin_theta_pow2 / (temp * temp)) *
                             (pow((ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                       factor_xi / (temp1 * temp) +
                                   ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))),
                                  2.0) +
                              cos_theta_pow2 / (temp * temp))),
                        (3.0 / 2.0)) -
                (((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                      (temp1 * temp) +
                  ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                     (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta * factor_xi /
                          (temp1 * temp) +
                      ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                 sin_theta * cos_theta / (temp * temp)) *
                    ((-27.0336) * (r / Rmax) * pow(((r / Rmax) - 1.0), 5.0) * sin(11.0 * theta) -
                     27.0336 * pow(((r / Rmax) - 1.0), 6.0) * sin(11.0 * theta)) *
                    exp(-tanh(20.0 * (r / Rmax) - 14.0)) /
                    sqrt((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                     factor_xi / (temp1 * temp) +
                                 ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                                    (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                         factor_xi / (temp1 * temp) +
                                     ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                                sin_theta * cos_theta / (temp * temp)),
                               2.0)) +
                         (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                   factor_xi / (temp1 * temp) +
                               ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                              2.0) +
                          sin_theta_pow2 / (temp * temp)) *
                             (pow((ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                       factor_xi / (temp1 * temp) +
                                   ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))),
                                  2.0) +
                              cos_theta_pow2 / (temp * temp))) -
                (((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                      (temp1 * temp) +
                  ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                     (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta * factor_xi /
                          (temp1 * temp) +
                      ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                 sin_theta * cos_theta / (temp * temp)) *
                    (2.4576 * (r / Rmax) * pow(((r / Rmax) - 1.0), 5.0) * cos(11.0 * theta) +
                     2.4576 * pow(((r / Rmax) - 1.0), 6.0) * cos(11.0 * theta)) *
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
                              ((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                   factor_xi / (temp1 * temp) +
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
                                   ((r / Rmax) * (r / Rmax)) * pow(sin_theta, 3.0) * factor_xi /
                                   (temp2 * (temp * temp)) -
                               3.0 * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                   factor_xi / (temp1 * temp) -
                               ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) +
                          2.0 * sin_theta_pow2 / (temp * temp) - 2.0 * cos_theta_pow2 / (temp * temp)) -
                     1.0 / 2.0 *
                         (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                   factor_xi / (temp1 * temp) +
                               ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                              2.0) +
                          sin_theta_pow2 / (temp * temp)) *
                         (2.0 * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta_pow2 /
                              pow((temp * temp), 2.0) +
                          (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                               factor_xi / (temp1 * temp) +
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
                         (2.0 * inverse_aspect_ratio_epsilon * (r / Rmax) * pow(sin_theta, 3.0) /
                              pow((temp * temp), 2.0) +
                          ((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                               (temp1 * temp) +
                           ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                              ((-2.0) * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                                   ((r / Rmax) * (r / Rmax)) * pow(sin_theta, 3.0) * factor_xi /
                                   (temp1 * pow((temp * temp), (3.0 / 2.0))) +
                               4.0 * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                                   ((r / Rmax) * (r / Rmax)) * pow(sin_theta, 3.0) * factor_xi /
                                   (temp2 * (temp * temp)) -
                               6.0 * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                   factor_xi / (temp1 * temp) -
                               2.0 * ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) +
                          2.0 * sin_theta * cos_theta / (temp * temp))) *
                    exp(-tanh(20.0 * (r / Rmax) - 14.0)) /
                    pow(((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                     factor_xi / (temp1 * temp) +
                                 ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                                    (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                         factor_xi / (temp1 * temp) +
                                     ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                                sin_theta * cos_theta / (temp * temp)),
                               2.0)) +
                         (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                   factor_xi / (temp1 * temp) +
                               ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                              2.0) +
                          sin_theta_pow2 / (temp * temp)) *
                             (pow((ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                       factor_xi / (temp1 * temp) +
                                   ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))),
                                  2.0) +
                              cos_theta_pow2 / (temp * temp))),
                        (3.0 / 2.0)) +
                6.0 *
                    (2.4576 * (r / Rmax) * pow(((r / Rmax) - 1.0), 5.0) * cos(11.0 * theta) +
                     2.4576 * pow(((r / Rmax) - 1.0), 6.0) * cos(11.0 * theta)) *
                    (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                              (temp1 * temp) +
                          ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                         2.0) +
                     sin_theta_pow2 / (temp * temp)) *
                    exp(-tanh(20.0 * (r / Rmax) - 14.0)) /
                    sqrt((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                     factor_xi / (temp1 * temp) +
                                 ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                                    (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                         factor_xi / (temp1 * temp) +
                                     ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                                sin_theta * cos_theta / (temp * temp)),
                               2.0)) +
                         (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                   factor_xi / (temp1 * temp) +
                               ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                              2.0) +
                          sin_theta_pow2 / (temp * temp)) *
                             (pow((ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                       factor_xi / (temp1 * temp) +
                                   ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))),
                                  2.0) +
                              cos_theta_pow2 / (temp * temp))) -
                (2.4576 * (r / Rmax) * pow(((r / Rmax) - 1.0), 5.0) * cos(11.0 * theta) +
                 2.4576 * pow(((r / Rmax) - 1.0), 6.0) * cos(11.0 * theta)) *
                    ((-2.0) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * cos_theta /
                         pow((temp * temp), 2.0) +
                     ((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                          (temp1 * temp) +
                      ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                         (ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                              ((r / Rmax) * (r / Rmax)) * sin_theta_pow2 * cos_theta * factor_xi /
                              (temp1 * pow((temp * temp), (3.0 / 2.0))) -
                          2.0 * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                              ((r / Rmax) * (r / Rmax)) * sin_theta_pow2 * cos_theta * factor_xi /
                              (temp2 * (temp * temp)) -
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
                    sqrt((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                     factor_xi / (temp1 * temp) +
                                 ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                                    (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                         factor_xi / (temp1 * temp) +
                                     ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                                sin_theta * cos_theta / (temp * temp)),
                               2.0)) +
                         (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                   factor_xi / (temp1 * temp) +
                               ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                              2.0) +
                          sin_theta_pow2 / (temp * temp)) *
                             (pow((ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                       factor_xi / (temp1 * temp) +
                                   ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))),
                                  2.0) +
                              cos_theta_pow2 / (temp * temp)))) /
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
                         cos_theta_pow2 / (temp * temp)));
}