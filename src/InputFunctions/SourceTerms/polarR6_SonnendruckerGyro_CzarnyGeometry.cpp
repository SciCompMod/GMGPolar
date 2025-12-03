#include "../include/InputFunctions/SourceTerms/polarR6_SonnendruckerGyro_CzarnyGeometry.h"

void PolarR6_SonnendruckerGyro_CzarnyGeometry::initializeGeometry()
{
    factor_xi = 1.0 / sqrt(1.0 - inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon / 4.0);
}

PolarR6_SonnendruckerGyro_CzarnyGeometry::PolarR6_SonnendruckerGyro_CzarnyGeometry(
    const double& Rmax, const double& inverse_aspect_ratio_epsilon, const double& ellipticity_e)
    : Rmax(Rmax)
    , inverse_aspect_ratio_epsilon(inverse_aspect_ratio_epsilon)
    , ellipticity_e(ellipticity_e)
{
    initializeGeometry();
}

double PolarR6_SonnendruckerGyro_CzarnyGeometry::rhs_f(const double& r, const double& theta)const
{
    double sin_theta = std::sin(theta);
    double cos_theta = std::cos(theta);
    double temp =
        sqrt(inverse_aspect_ratio_epsilon * (inverse_aspect_ratio_epsilon + 2.0 * (r / Rmax) * cos_theta) + 1.0);
    double sin_theta_pow2 = pow(sin_theta, 2.0);
    double cos_theta_pow2 = pow(cos_theta, 2.0);
    double temp1          = pow((2.0 - temp), 2.0);
    double temp2          = pow((2.0 - temp), 3.0);

    return 0.4096 * pow((r / Rmax), 6.0) * pow(((r / Rmax) - 1.0), 6.0) * cos(11.0 * theta) /
               (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r / Rmax) - 11.1111111111111)) -
           pow((r / Rmax), 4.0) *
               (4.5056 * (r / Rmax) *
                    (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r / Rmax) - 11.1111111111111)) *
                    pow(((r / Rmax) - 1.0), 6.0) *
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
                    sin(11.0 * theta) /
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
                4.5056 * (r / Rmax) *
                    (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r / Rmax) - 11.1111111111111)) *
                    pow(((r / Rmax) - 1.0), 6.0) *
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
                27.0336 * (r / Rmax) *
                    (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r / Rmax) - 11.1111111111111)) *
                    pow(((r / Rmax) - 1.0), 5.0) *
                    (((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                          (temp1 * temp) +
                      ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                         (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                              factor_xi / (temp1 * temp) +
                          ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                     sin_theta * cos_theta / (temp * temp)) *
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
                (r / Rmax) *
                    (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r / Rmax) - 11.1111111111111)) *
                    (12.288 * (r / Rmax) * pow(((r / Rmax) - 1.0), 4.0) * cos(11.0 * theta) +
                     17.2032 * pow(((r / Rmax) - 1.0), 5.0) * cos(11.0 * theta)) *
                    (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                              (temp1 * temp) +
                          ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                         2.0) +
                     sin_theta_pow2 / (temp * temp)) /
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
                    (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r / Rmax) - 11.1111111111111)) *
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
                              (temp1 * temp))) /
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
                    (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r / Rmax) - 11.1111111111111)) *
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
                          cos_theta_pow2 / (temp * temp))) /
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
                22.6762679055362 * (r / Rmax) * pow(((r / Rmax) - 1.0), 6.0) *
                    (((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                          (temp1 * temp) +
                      ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                         (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                              factor_xi / (temp1 * temp) +
                          ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                     sin_theta * cos_theta / (temp * temp)) *
                    sin(11.0 * theta) /
                    ((208.641975308642 * pow(((r / Rmax) - 0.769230769230769), 2.0) + 1.0) *
                     sqrt((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                      factor_xi / (temp1 * temp) +
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
                              (pow((ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                        factor_xi / (temp1 * temp) +
                                    ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))),
                                   2.0) +
                               cos_theta_pow2 / (temp * temp)))) -
                5.03290747193186 * (r / Rmax) *
                    (2.4576 * (r / Rmax) * pow(((r / Rmax) - 1.0), 5.0) * cos(11.0 * theta) +
                     2.4576 * pow(((r / Rmax) - 1.0), 6.0) * cos(11.0 * theta)) *
                    (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                              (temp1 * temp) +
                          ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                         2.0) +
                     sin_theta_pow2 / (temp * temp)) /
                    ((208.641975308642 * pow(((r / Rmax) - 0.769230769230769), 2.0) + 1.0) *
                     sqrt((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                      factor_xi / (temp1 * temp) +
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
                              (pow((ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                                        factor_xi / (temp1 * temp) +
                                    ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))),
                                   2.0) +
                               cos_theta_pow2 / (temp * temp)))) +
                27.0336 *
                    (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r / Rmax) - 11.1111111111111)) *
                    pow(((r / Rmax) - 1.0), 6.0) *
                    (((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                          (temp1 * temp) +
                      ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                         (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                              factor_xi / (temp1 * temp) +
                          ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                     sin_theta * cos_theta / (temp * temp)) *
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
                              cos_theta_pow2 / (temp * temp))) -
                49.5616 *
                    (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r / Rmax) - 11.1111111111111)) *
                    pow(((r / Rmax) - 1.0), 6.0) *
                    (pow((ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                              factor_xi / (temp1 * temp) +
                          ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))),
                         2.0) +
                     cos_theta_pow2 / (temp * temp)) *
                    cos(11.0 * theta) /
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
                4.5056 *
                    (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r / Rmax) - 11.1111111111111)) *
                    pow(((r / Rmax) - 1.0), 6.0) *
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
                              cos_theta_pow2 / (temp * temp))) -
                4.5056 *
                    (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r / Rmax) - 11.1111111111111)) *
                    pow(((r / Rmax) - 1.0), 6.0) *
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
                    sin(11.0 * theta) /
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
                (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r / Rmax) - 11.1111111111111)) *
                    (((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                          (temp1 * temp) +
                      ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                         (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                              factor_xi / (temp1 * temp) +
                          ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                     sin_theta * cos_theta / (temp * temp)) *
                    ((-27.0336) * (r / Rmax) * pow(((r / Rmax) - 1.0), 5.0) * sin(11.0 * theta) -
                     27.0336 * pow(((r / Rmax) - 1.0), 6.0) * sin(11.0 * theta)) /
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
                (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r / Rmax) - 11.1111111111111)) *
                    (((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                          (temp1 * temp) +
                      ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                         (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                              factor_xi / (temp1 * temp) +
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
                          2.0 * sin_theta * cos_theta / (temp * temp))) /
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
                6.0 * (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r / Rmax) - 11.1111111111111)) *
                    (2.4576 * (r / Rmax) * pow(((r / Rmax) - 1.0), 5.0) * cos(11.0 * theta) +
                     2.4576 * pow(((r / Rmax) - 1.0), 6.0) * cos(11.0 * theta)) *
                    (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                              (temp1 * temp) +
                          ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                         2.0) +
                     sin_theta_pow2 / (temp * temp)) /
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
                (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r / Rmax) - 11.1111111111111)) *
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
                     sin_theta_pow2 / (temp * temp) - cos_theta_pow2 / (temp * temp)) /
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