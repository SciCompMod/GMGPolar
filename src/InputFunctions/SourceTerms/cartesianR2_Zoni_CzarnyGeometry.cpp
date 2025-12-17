#include "../include/InputFunctions/SourceTerms/cartesianR2_Zoni_CzarnyGeometry.h"

void CartesianR2_Zoni_CzarnyGeometry::initializeGeometry()
{
    factor_xi = 1.0 / sqrt(1.0 - inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon / 4.0);
}

CartesianR2_Zoni_CzarnyGeometry::CartesianR2_Zoni_CzarnyGeometry(PolarGrid const& grid, double Rmax,
                                                                 double inverse_aspect_ratio_epsilon,
                                                                 double ellipticity_e)
    : grid_(grid)
    , Rmax(Rmax)
    , inverse_aspect_ratio_epsilon(inverse_aspect_ratio_epsilon)
    , ellipticity_e(ellipticity_e)
{
    initializeGeometry();
}

double CartesianR2_Zoni_CzarnyGeometry::operator()(int i_r, int i_theta) const
{
    double r         = grid_.radius(i_r);
    double theta     = grid_.theta(i_theta);
    double sin_theta = std::sin(theta);
    double cos_theta = std::cos(theta);
    double temp =
        sqrt(inverse_aspect_ratio_epsilon * (inverse_aspect_ratio_epsilon + 2.0 * (r / Rmax) * cos_theta) + 1.0);
    double sin_theta_pow2 = pow(sin_theta, 2.0);
    double cos_theta_pow2 = pow(cos_theta, 2.0);
    double temp1          = pow((2.0 - temp), 2.0);
    double temp2          = pow((2.0 - temp), 3.0);

    return (-((-(r / Rmax)) *
                  (((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                        (temp1 * temp) +
                    ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                       (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta * factor_xi /
                            (temp1 * temp) +
                        ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                   sin_theta * cos_theta / (temp * temp)) *
                  ((1.0 - (r / Rmax) * (r / Rmax)) *
                       ((-2.0) * M_PI * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                            factor_xi / (temp1 * temp) +
                        2.0 * M_PI * ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                       cos(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) *
                       cos(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) -
                   2.0 * M_PI * (1.0 - (r / Rmax) * (r / Rmax)) * sin_theta *
                       sin(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) *
                       sin(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) / temp) *
                  (10.0 * pow(tanh(10.0 * (r / Rmax) - 5.0), 2.0) - 10.0) * exp(-tanh(10.0 * (r / Rmax) - 5.0)) /
                  sqrt((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                   factor_xi / (temp1 * temp) +
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
              (r / Rmax) *
                  (((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                        (temp1 * temp) +
                    ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                       (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta * factor_xi /
                            (temp1 * temp) +
                        ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                   sin_theta * cos_theta / (temp * temp)) *
                  ((1.0 - (r / Rmax) * (r / Rmax)) *
                       ((-2.0) * M_PI * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                            factor_xi / (temp1 * temp) +
                        2.0 * M_PI * ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                       cos(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) *
                       cos(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) -
                   2.0 * M_PI * (1.0 - (r / Rmax) * (r / Rmax)) * sin_theta *
                       sin(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) *
                       sin(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) / temp) *
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
                  exp(-tanh(10.0 * (r / Rmax) - 5.0)) /
                  pow(((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                   factor_xi / (temp1 * temp) +
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
              (r / Rmax) *
                  (((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                        (temp1 * temp) +
                    ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                       (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta * factor_xi /
                            (temp1 * temp) +
                        ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
                   sin_theta * cos_theta / (temp * temp)) *
                  (2.0 * M_PI * inverse_aspect_ratio_epsilon * (1.0 - (r / Rmax) * (r / Rmax)) * sin_theta *
                       sin(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) *
                       sin(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) *
                       cos_theta / pow((temp * temp), (3.0 / 2.0)) -
                   2.0 * (r / Rmax) *
                       ((-2.0) * M_PI * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                            factor_xi / (temp1 * temp) +
                        2.0 * M_PI * ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                       cos(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) *
                       cos(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) +
                   4.0 * M_PI * (r / Rmax) * sin_theta * sin(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) *
                       sin(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) / temp -
                   (1.0 - (r / Rmax) * (r / Rmax)) *
                       ((-2.0) * M_PI * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                            factor_xi / (temp1 * temp) +
                        2.0 * M_PI * ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                       (2.0 * M_PI * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                            factor_xi / (temp1 * temp) +
                        2.0 * M_PI * ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) *
                       sin(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) *
                       cos(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) +
                   (1.0 - (r / Rmax) * (r / Rmax)) *
                       (2.0 * M_PI * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                            (r / Rmax) * sin_theta_pow2 * cos_theta * factor_xi /
                            (temp1 * pow((temp * temp), (3.0 / 2.0))) -
                        4.0 * M_PI * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                            (r / Rmax) * sin_theta_pow2 * cos_theta * factor_xi / (temp2 * (temp * temp)) -
                        2.0 * M_PI * ellipticity_e * inverse_aspect_ratio_epsilon * sin_theta_pow2 * factor_xi /
                            (temp1 * temp) +
                        2.0 * M_PI * ellipticity_e * inverse_aspect_ratio_epsilon * cos_theta_pow2 * factor_xi /
                            (temp1 * temp)) *
                       cos(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) *
                       cos(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) +
                   4.0 * (M_PI * M_PI) * (1.0 - (r / Rmax) * (r / Rmax)) * sin_theta *
                       sin(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) *
                       cos_theta * cos(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) / (temp * temp) +
                   2.0 * M_PI * (1.0 - (r / Rmax) * (r / Rmax)) *
                       ((-2.0) * M_PI * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                            factor_xi / (temp1 * temp) +
                        2.0 * M_PI * ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                       sin(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) * cos_theta *
                       cos(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) / temp -
                   2.0 * M_PI * (1.0 - (r / Rmax) * (r / Rmax)) *
                       (2.0 * M_PI * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                            factor_xi / (temp1 * temp) +
                        2.0 * M_PI * ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) *
                       sin_theta * sin(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) *
                       cos(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) / temp) *
                  exp(-tanh(10.0 * (r / Rmax) - 5.0)) /
                  sqrt((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                   factor_xi / (temp1 * temp) +
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
                  ((-2.0) * (r / Rmax) *
                       sin(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) *
                       cos(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) +
                   (1.0 - (r / Rmax) * (r / Rmax)) *
                       (2.0 * M_PI * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                            factor_xi / (temp1 * temp) +
                        2.0 * M_PI * ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) *
                       cos(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) *
                       cos(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) +
                   2.0 * M_PI * (1.0 - (r / Rmax) * (r / Rmax)) *
                       sin(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) *
                       sin(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) *
                       cos_theta / temp) *
                  exp(-tanh(10.0 * (r / Rmax) - 5.0)) /
                  sqrt((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                   factor_xi / (temp1 * temp) +
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
              (r / Rmax) *
                  ((1.0 - (r / Rmax) * (r / Rmax)) *
                       ((-2.0) * M_PI * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                            factor_xi / (temp1 * temp) +
                        2.0 * M_PI * ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                       cos(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) *
                       cos(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) -
                   2.0 * M_PI * (1.0 - (r / Rmax) * (r / Rmax)) * sin_theta *
                       sin(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) *
                       sin(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) / temp) *
                  (2.0 * inverse_aspect_ratio_epsilon * sin_theta * cos_theta_pow2 / pow((temp * temp), 2.0) +
                   ((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                        (temp1 * temp) +
                    ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                       ((-ellipticity_e) * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) * (r / Rmax) *
                            sin_theta * cos_theta_pow2 * factor_xi / (temp1 * pow((temp * temp), (3.0 / 2.0))) +
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
                  exp(-tanh(10.0 * (r / Rmax) - 5.0)) /
                  sqrt((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                   factor_xi / (temp1 * temp) +
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
                  (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                            (temp1 * temp) +
                        ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                       2.0) +
                   sin_theta_pow2 / (temp * temp)) *
                  (10.0 * pow(tanh(10.0 * (r / Rmax) - 5.0), 2.0) - 10.0) *
                  ((-2.0) * (r / Rmax) *
                       sin(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) *
                       cos(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) +
                   (1.0 - (r / Rmax) * (r / Rmax)) *
                       (2.0 * M_PI * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                            factor_xi / (temp1 * temp) +
                        2.0 * M_PI * ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) *
                       cos(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) *
                       cos(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) +
                   2.0 * M_PI * (1.0 - (r / Rmax) * (r / Rmax)) *
                       sin(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) *
                       sin(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) *
                       cos_theta / temp) *
                  exp(-tanh(10.0 * (r / Rmax) - 5.0)) /
                  sqrt((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                   factor_xi / (temp1 * temp) +
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
                  (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                            (temp1 * temp) +
                        ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                       2.0) +
                   sin_theta_pow2 / (temp * temp)) *
                  ((-2.0) * M_PI * inverse_aspect_ratio_epsilon * (1.0 - (r / Rmax) * (r / Rmax)) *
                       sin(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) *
                       sin(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) *
                       cos_theta_pow2 / pow((temp * temp), (3.0 / 2.0)) -
                   4.0 * (r / Rmax) *
                       (2.0 * M_PI * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                            factor_xi / (temp1 * temp) +
                        2.0 * M_PI * ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) *
                       cos(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) *
                       cos(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) -
                   8.0 * M_PI * (r / Rmax) * sin(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) *
                       sin(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) *
                       cos_theta / temp -
                   (1.0 - (r / Rmax) * (r / Rmax)) *
                       pow((2.0 * M_PI * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta *
                                cos_theta * factor_xi / (temp1 * temp) +
                            2.0 * M_PI * ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))),
                           2.0) *
                       sin(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) *
                       cos(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) +
                   (1.0 - (r / Rmax) * (r / Rmax)) *
                       ((-2.0) * M_PI * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                            (r / Rmax) * sin_theta * cos_theta_pow2 * factor_xi /
                            (temp1 * pow((temp * temp), (3.0 / 2.0))) +
                        4.0 * M_PI * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                            (r / Rmax) * sin_theta * cos_theta_pow2 * factor_xi / (temp2 * (temp * temp)) +
                        4.0 * M_PI * ellipticity_e * inverse_aspect_ratio_epsilon * sin_theta * cos_theta * factor_xi /
                            (temp1 * temp)) *
                       cos(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) *
                       cos(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) -
                   4.0 * (M_PI * M_PI) * (1.0 - (r / Rmax) * (r / Rmax)) *
                       sin(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) *
                       cos_theta_pow2 * cos(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) / (temp * temp) +
                   4.0 * M_PI * (1.0 - (r / Rmax) * (r / Rmax)) *
                       (2.0 * M_PI * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                            factor_xi / (temp1 * temp) +
                        2.0 * M_PI * ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) *
                       sin(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) * cos_theta *
                       cos(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) / temp -
                   2.0 * sin(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) *
                       cos(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon)) *
                  exp(-tanh(10.0 * (r / Rmax) - 5.0)) /
                  sqrt((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                   factor_xi / (temp1 * temp) +
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
                  ((-2.0) * (r / Rmax) *
                       sin(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) *
                       cos(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) +
                   (1.0 - (r / Rmax) * (r / Rmax)) *
                       (2.0 * M_PI * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                            factor_xi / (temp1 * temp) +
                        2.0 * M_PI * ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) *
                       cos(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) *
                       cos(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) +
                   2.0 * M_PI * (1.0 - (r / Rmax) * (r / Rmax)) *
                       sin(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) *
                       sin(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) *
                       cos_theta / temp) *
                  exp(-tanh(10.0 * (r / Rmax) - 5.0)) /
                  pow(((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                   factor_xi / (temp1 * temp) +
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
              (((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                    (temp1 * temp) +
                ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                   (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta * factor_xi /
                        (temp1 * temp) +
                    ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
               sin_theta * cos_theta / (temp * temp)) *
                  ((1.0 - (r / Rmax) * (r / Rmax)) *
                       ((-2.0) * M_PI * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                            factor_xi / (temp1 * temp) +
                        2.0 * M_PI * ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                       cos(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) *
                       cos(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) -
                   2.0 * M_PI * (1.0 - (r / Rmax) * (r / Rmax)) * sin_theta *
                       sin(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) *
                       sin(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) / temp) *
                  exp(-tanh(10.0 * (r / Rmax) - 5.0)) /
                  sqrt((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                   factor_xi / (temp1 * temp) +
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
              (((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                    (temp1 * temp) +
                ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                   (ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta * factor_xi /
                        (temp1 * temp) +
                    ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) -
               sin_theta * cos_theta / (temp * temp)) *
                  (2.0 * M_PI * inverse_aspect_ratio_epsilon * (r / Rmax) * (1.0 - (r / Rmax) * (r / Rmax)) *
                       sin_theta * sin(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) *
                       sin(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) *
                       cos_theta / pow((temp * temp), (3.0 / 2.0)) -
                   2.0 * ((r / Rmax) * (r / Rmax)) *
                       ((-2.0) * M_PI * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                            factor_xi / (temp1 * temp) +
                        2.0 * M_PI * ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                       cos(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) *
                       cos(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) +
                   4.0 * M_PI * ((r / Rmax) * (r / Rmax)) * sin_theta *
                       sin(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) *
                       sin(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) / temp -
                   (r / Rmax) * (1.0 - (r / Rmax) * (r / Rmax)) *
                       ((-2.0) * M_PI * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                            factor_xi / (temp1 * temp) +
                        2.0 * M_PI * ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                       (2.0 * M_PI * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                            factor_xi / (temp1 * temp) +
                        2.0 * M_PI * ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) *
                       sin(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) *
                       cos(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) +
                   4.0 * (M_PI * M_PI) * (r / Rmax) * (1.0 - (r / Rmax) * (r / Rmax)) * sin_theta *
                       sin(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) *
                       cos_theta * cos(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) / (temp * temp) +
                   2.0 * M_PI * (r / Rmax) * (1.0 - (r / Rmax) * (r / Rmax)) *
                       ((-2.0) * M_PI * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                            factor_xi / (temp1 * temp) +
                        2.0 * M_PI * ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                       sin(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) * cos_theta *
                       cos(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) / temp -
                   2.0 * M_PI * (r / Rmax) * (1.0 - (r / Rmax) * (r / Rmax)) *
                       (2.0 * M_PI * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                            factor_xi / (temp1 * temp) +
                        2.0 * M_PI * ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) *
                       sin_theta * sin(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) *
                       cos(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) / temp +
                   (1.0 - (r / Rmax) * (r / Rmax)) *
                       (2.0 * M_PI * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                            ((r / Rmax) * (r / Rmax)) * sin_theta_pow2 * cos_theta * factor_xi /
                            (temp1 * pow((temp * temp), (3.0 / 2.0))) -
                        4.0 * M_PI * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                            ((r / Rmax) * (r / Rmax)) * sin_theta_pow2 * cos_theta * factor_xi /
                            (temp2 * (temp * temp)) -
                        4.0 * M_PI * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                            factor_xi / (temp1 * temp) +
                        2.0 * M_PI * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * cos_theta_pow2 *
                            factor_xi / (temp1 * temp) +
                        2.0 * M_PI * ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                       cos(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) *
                       cos(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) -
                   2.0 * M_PI * (1.0 - (r / Rmax) * (r / Rmax)) * sin_theta *
                       sin(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) *
                       sin(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) / temp) *
                  exp(-tanh(10.0 * (r / Rmax) - 5.0)) /
                  sqrt((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                   factor_xi / (temp1 * temp) +
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
                       (2.0 * inverse_aspect_ratio_epsilon * (r / Rmax) * pow(sin_theta, 3.0) /
                            pow((temp * temp), 2.0) +
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
                  ((-2.0) * (r / Rmax) *
                       sin(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) *
                       cos(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) +
                   (1.0 - (r / Rmax) * (r / Rmax)) *
                       (2.0 * M_PI * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                            factor_xi / (temp1 * temp) +
                        2.0 * M_PI * ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) *
                       cos(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) *
                       cos(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) +
                   2.0 * M_PI * (1.0 - (r / Rmax) * (r / Rmax)) *
                       sin(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) *
                       sin(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) *
                       cos_theta / temp) *
                  exp(-tanh(10.0 * (r / Rmax) - 5.0)) /
                  pow(((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                   factor_xi / (temp1 * temp) +
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
              ((1.0 - (r / Rmax) * (r / Rmax)) *
                   ((-2.0) * M_PI * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                        factor_xi / (temp1 * temp) +
                    2.0 * M_PI * ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                   cos(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) *
                   cos(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) -
               2.0 * M_PI * (1.0 - (r / Rmax) * (r / Rmax)) * sin_theta *
                   sin(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) *
                   sin(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) / temp) *
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
                  exp(-tanh(10.0 * (r / Rmax) - 5.0)) /
                  sqrt((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                   factor_xi / (temp1 * temp) +
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
              ((1.0 - (r / Rmax) * (r / Rmax)) *
                   ((-2.0) * M_PI * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                        factor_xi / (temp1 * temp) +
                    2.0 * M_PI * ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                   cos(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) *
                   cos(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) -
               2.0 * M_PI * (1.0 - (r / Rmax) * (r / Rmax)) * sin_theta *
                   sin(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) *
                   sin(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) / temp) *
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
                       (2.0 * inverse_aspect_ratio_epsilon * (r / Rmax) * pow(sin_theta, 3.0) /
                            pow((temp * temp), 2.0) +
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
                  exp(-tanh(10.0 * (r / Rmax) - 5.0)) /
                  pow(((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                   factor_xi / (temp1 * temp) +
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
              (pow(((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 * factor_xi /
                        (temp1 * temp) +
                    ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                   2.0) +
               sin_theta_pow2 / (temp * temp)) *
                  ((-2.0) * (r / Rmax) *
                       sin(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) *
                       cos(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) +
                   (1.0 - (r / Rmax) * (r / Rmax)) *
                       (2.0 * M_PI * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                            factor_xi / (temp1 * temp) +
                        2.0 * M_PI * ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) *
                       cos(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) *
                       cos(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) +
                   2.0 * M_PI * (1.0 - (r / Rmax) * (r / Rmax)) *
                       sin(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) *
                       sin(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) *
                       cos_theta / temp) *
                  exp(-tanh(10.0 * (r / Rmax) - 5.0)) /
                  sqrt((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                   factor_xi / (temp1 * temp) +
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
              (pow((ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta * factor_xi /
                        (temp1 * temp) +
                    ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))),
                   2.0) +
               cos_theta_pow2 / (temp * temp)) *
                  ((-2.0) * M_PI * inverse_aspect_ratio_epsilon * (r / Rmax) * (1.0 - (r / Rmax) * (r / Rmax)) *
                       sin_theta_pow2 * sin(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) *
                       sin(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) /
                       pow((temp * temp), (3.0 / 2.0)) -
                   (r / Rmax) * (1.0 - (r / Rmax) * (r / Rmax)) *
                       pow(((-2.0) * M_PI * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                factor_xi / (temp1 * temp) +
                            2.0 * M_PI * ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))),
                           2.0) *
                       sin(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) *
                       cos(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) -
                   4.0 * (M_PI * M_PI) * (r / Rmax) * (1.0 - (r / Rmax) * (r / Rmax)) * sin_theta_pow2 *
                       sin(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) *
                       cos(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) / (temp * temp) -
                   4.0 * M_PI * (r / Rmax) * (1.0 - (r / Rmax) * (r / Rmax)) *
                       ((-2.0) * M_PI * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                            factor_xi / (temp1 * temp) +
                        2.0 * M_PI * ellipticity_e * cos_theta * factor_xi / ((2.0 - temp))) *
                       sin_theta * sin(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) *
                       cos(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) / temp +
                   (1.0 - (r / Rmax) * (r / Rmax)) *
                       ((-2.0) * M_PI * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                            ((r / Rmax) * (r / Rmax)) * pow(sin_theta, 3.0) * factor_xi /
                            (temp1 * pow((temp * temp), (3.0 / 2.0))) +
                        4.0 * M_PI * ellipticity_e * (inverse_aspect_ratio_epsilon * inverse_aspect_ratio_epsilon) *
                            ((r / Rmax) * (r / Rmax)) * pow(sin_theta, 3.0) * factor_xi / (temp2 * (temp * temp)) -
                        6.0 * M_PI * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                            factor_xi / (temp1 * temp) -
                        2.0 * M_PI * ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) *
                       cos(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) *
                       cos(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) -
                   2.0 * M_PI * (1.0 - (r / Rmax) * (r / Rmax)) *
                       sin(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) *
                       sin(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) *
                       cos_theta / temp) *
                  exp(-tanh(10.0 * (r / Rmax) - 5.0)) /
                  sqrt((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                   factor_xi / (temp1 * temp) +
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
              ((-2.0) * (r / Rmax) *
                   sin(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) *
                   cos(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) +
               (1.0 - (r / Rmax) * (r / Rmax)) *
                   (2.0 * M_PI * ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta * cos_theta *
                        factor_xi / (temp1 * temp) +
                    2.0 * M_PI * ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))) *
                   cos(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) *
                   cos(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) +
               2.0 * M_PI * (1.0 - (r / Rmax) * (r / Rmax)) *
                   sin(2.0 * M_PI * (1.0 - temp) / inverse_aspect_ratio_epsilon) *
                   sin(2.0 * M_PI * ellipticity_e * (r / Rmax) * sin_theta * factor_xi / ((2.0 - temp))) * cos_theta /
                   temp) *
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
                  exp(-tanh(10.0 * (r / Rmax) - 5.0)) /
                  sqrt((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
                                   factor_xi / (temp1 * temp) +
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
                            cos_theta_pow2 / (temp * temp))))) /
           ((r / Rmax) * sqrt((-pow((((-ellipticity_e) * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta_pow2 *
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
                                  (pow((ellipticity_e * inverse_aspect_ratio_epsilon * (r / Rmax) * sin_theta *
                                            cos_theta * factor_xi / (temp1 * temp) +
                                        ellipticity_e * sin_theta * factor_xi / ((2.0 - temp))),
                                       2.0) +
                                   cos_theta_pow2 / (temp * temp))));
}