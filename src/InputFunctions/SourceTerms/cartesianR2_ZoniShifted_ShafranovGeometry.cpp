#include "../include/InputFunctions/SourceTerms/cartesianR2_ZoniShifted_ShafranovGeometry.h"

CartesianR2_ZoniShifted_ShafranovGeometry::CartesianR2_ZoniShifted_ShafranovGeometry(const double& Rmax,
                                                                                     const double& elongation_kappa,
                                                                                     const double& shift_delta)
    : Rmax(Rmax)
    , elongation_kappa(elongation_kappa)
    , shift_delta(shift_delta)
{
}

double CartesianR2_ZoniShifted_ShafranovGeometry::rhs_f(const double& r, const double& theta, const double& sin_theta,
                                                        const double& cos_theta) const
{
    return (-(2.0 * shift_delta * (r / Rmax) * (elongation_kappa - 1.0) *
                  ((-M_PI) * (1.0 - (r / Rmax) * (r / Rmax)) * (2.0 * elongation_kappa - 2.0) * sin_theta *
                       sin(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
                       sin(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta) +
                   M_PI * (1.0 - (r / Rmax) * (r / Rmax)) * (2.0 * elongation_kappa + 2.0) * cos_theta *
                       cos(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
                       cos(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta)) *
                  exp(-tanh(20.0 * (r / Rmax) - 14.0)) * sin_theta /
                  sqrt((pow((elongation_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) +
                        pow(((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta), 2.0)) *
                           (elongation_kappa * elongation_kappa - 4.0 * elongation_kappa * pow(sin_theta, 2.0) +
                            2.0 * elongation_kappa + 1.0) -
                       pow(((elongation_kappa - 1.0) *
                                ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) *
                                sin_theta +
                            1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)),
                           2.0)) -
              (r / Rmax) *
                  ((elongation_kappa - 1.0) *
                       ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) * sin_theta +
                   1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)) *
                  ((-M_PI) * (1.0 - (r / Rmax) * (r / Rmax)) * (2.0 * elongation_kappa - 2.0) * sin_theta *
                       sin(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
                       sin(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta) +
                   M_PI * (1.0 - (r / Rmax) * (r / Rmax)) * (2.0 * elongation_kappa + 2.0) * cos_theta *
                       cos(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
                       cos(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta)) *
                  (20.0 * pow(tanh(20.0 * (r / Rmax) - 14.0), 2.0) - 20.0) * exp(-tanh(20.0 * (r / Rmax) - 14.0)) /
                  sqrt((pow((elongation_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) +
                        pow(((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta), 2.0)) *
                           (elongation_kappa * elongation_kappa - 4.0 * elongation_kappa * pow(sin_theta, 2.0) +
                            2.0 * elongation_kappa + 1.0) -
                       pow(((elongation_kappa - 1.0) *
                                ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) *
                                sin_theta +
                            1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)),
                           2.0)) -
              (r / Rmax) *
                  ((elongation_kappa - 1.0) *
                       ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) * sin_theta +
                   1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)) *
                  (2.0 * M_PI * (r / Rmax) * (2.0 * elongation_kappa - 2.0) * sin_theta *
                       sin(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
                       sin(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta) -
                   2.0 * M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * cos_theta *
                       cos(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
                       cos(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta) -
                   4.0 * (M_PI * M_PI) * (1.0 - (r / Rmax) * (r / Rmax)) * pow((elongation_kappa + 1.0), 2.0) *
                       sin_theta * sin(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta) * cos_theta *
                       cos(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) -
                   M_PI * M_PI * (1.0 - (r / Rmax) * (r / Rmax)) * (2.0 * elongation_kappa - 2.0) *
                       (2.0 * elongation_kappa + 2.0) * pow(sin_theta, 2.0) *
                       sin(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
                       cos(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta) -
                   M_PI * (1.0 - (r / Rmax) * (r / Rmax)) * (2.0 * elongation_kappa - 2.0) *
                       ((-2.0) * M_PI * shift_delta * (r / Rmax) +
                        M_PI * ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta +
                                2.0 * cos_theta)) *
                       sin_theta * sin(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta) *
                       cos(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) -
                   M_PI * (1.0 - (r / Rmax) * (r / Rmax)) * (2.0 * elongation_kappa + 2.0) *
                       ((-2.0) * M_PI * shift_delta * (r / Rmax) +
                        M_PI * ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta +
                                2.0 * cos_theta)) *
                       sin(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
                       cos_theta * cos(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta)) *
                  exp(-tanh(20.0 * (r / Rmax) - 14.0)) /
                  sqrt((pow((elongation_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) +
                        pow(((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta), 2.0)) *
                           (elongation_kappa * elongation_kappa - 4.0 * elongation_kappa * pow(sin_theta, 2.0) +
                            2.0 * elongation_kappa + 1.0) -
                       pow(((elongation_kappa - 1.0) *
                                ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) *
                                sin_theta +
                            1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)),
                           2.0)) +
              (r / Rmax) * (20.0 * pow(tanh(20.0 * (r / Rmax) - 14.0), 2.0) - 20.0) *
                  ((-2.0) * (r / Rmax) * sin(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta) *
                       cos(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) +
                   M_PI * (1.0 - (r / Rmax) * (r / Rmax)) * (2.0 * elongation_kappa + 2.0) * sin_theta *
                       cos(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
                       cos(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta) -
                   M_PI * (1.0 - (r / Rmax) * (r / Rmax)) *
                       ((-4.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta) *
                       sin(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
                       sin(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta)) *
                  (elongation_kappa * elongation_kappa - 4.0 * elongation_kappa * pow(sin_theta, 2.0) +
                   2.0 * elongation_kappa + 1.0) *
                  exp(-tanh(20.0 * (r / Rmax) - 14.0)) /
                  sqrt((pow((elongation_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) +
                        pow(((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta), 2.0)) *
                           (elongation_kappa * elongation_kappa - 4.0 * elongation_kappa * pow(sin_theta, 2.0) +
                            2.0 * elongation_kappa + 1.0) -
                       pow(((elongation_kappa - 1.0) *
                                ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) *
                                sin_theta +
                            1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)),
                           2.0)) +
              (r / Rmax) *
                  (elongation_kappa * elongation_kappa - 4.0 * elongation_kappa * pow(sin_theta, 2.0) +
                   2.0 * elongation_kappa + 1.0) *
                  (4.0 * M_PI * shift_delta * (1.0 - (r / Rmax) * (r / Rmax)) *
                       sin(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
                       sin(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta) -
                   4.0 * M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta *
                       cos(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
                       cos(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta) +
                   2.0 * (r / Rmax) *
                       ((-2.0) * M_PI * shift_delta * (r / Rmax) +
                        M_PI * ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta +
                                2.0 * cos_theta)) *
                       sin(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
                       sin(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta) +
                   2.0 * M_PI * (r / Rmax) *
                       ((-4.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta) *
                       sin(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
                       sin(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta) -
                   4.0 * (M_PI * M_PI) * (1.0 - (r / Rmax) * (r / Rmax)) * pow((elongation_kappa + 1.0), 2.0) *
                       pow(sin_theta, 2.0) * sin(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta) *
                       cos(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) -
                   M_PI * (1.0 - (r / Rmax) * (r / Rmax)) * (2.0 * elongation_kappa + 2.0) *
                       ((-2.0) * M_PI * shift_delta * (r / Rmax) +
                        M_PI * ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta +
                                2.0 * cos_theta)) *
                       sin_theta *
                       sin(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
                       cos(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta) -
                   M_PI * M_PI * (1.0 - (r / Rmax) * (r / Rmax)) * (2.0 * elongation_kappa + 2.0) *
                       ((-4.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta) *
                       sin_theta *
                       sin(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
                       cos(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta) -
                   M_PI * (1.0 - (r / Rmax) * (r / Rmax)) *
                       ((-2.0) * M_PI * shift_delta * (r / Rmax) +
                        M_PI * ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta +
                                2.0 * cos_theta)) *
                       ((-4.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta) *
                       sin(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta) *
                       cos(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) -
                   2.0 * sin(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta) *
                       cos(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta +
                            2.0 * cos_theta))) *
                  exp(-tanh(20.0 * (r / Rmax) - 14.0)) /
                  sqrt((pow((elongation_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) +
                        pow(((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta), 2.0)) *
                           (elongation_kappa * elongation_kappa - 4.0 * elongation_kappa * pow(sin_theta, 2.0) +
                            2.0 * elongation_kappa + 1.0) -
                       pow(((elongation_kappa - 1.0) *
                                ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) *
                                sin_theta +
                            1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)),
                           2.0)) -
              (r / Rmax) *
                  ((elongation_kappa - 1.0) *
                       ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) * sin_theta +
                   1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)) *
                  ((-2.0) * shift_delta * (elongation_kappa - 1.0) *
                       ((elongation_kappa - 1.0) *
                            ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) * sin_theta +
                        1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)) *
                       sin_theta +
                   2.0 * shift_delta * ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) *
                       (elongation_kappa * elongation_kappa - 4.0 * elongation_kappa * pow(sin_theta, 2.0) +
                        2.0 * elongation_kappa + 1.0)) *
                  ((-M_PI) * (1.0 - (r / Rmax) * (r / Rmax)) * (2.0 * elongation_kappa - 2.0) * sin_theta *
                       sin(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
                       sin(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta) +
                   M_PI * (1.0 - (r / Rmax) * (r / Rmax)) * (2.0 * elongation_kappa + 2.0) * cos_theta *
                       cos(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
                       cos(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta)) *
                  exp(-tanh(20.0 * (r / Rmax) - 14.0)) /
                  pow(((pow((elongation_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) +
                        pow(((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta), 2.0)) *
                           (elongation_kappa * elongation_kappa - 4.0 * elongation_kappa * pow(sin_theta, 2.0) +
                            2.0 * elongation_kappa + 1.0) -
                       pow(((elongation_kappa - 1.0) *
                                ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) *
                                sin_theta +
                            1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)),
                           2.0)),
                      (3.0 / 2.0)) +
              (r / Rmax) *
                  ((-2.0) * shift_delta * (elongation_kappa - 1.0) *
                       ((elongation_kappa - 1.0) *
                            ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) * sin_theta +
                        1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)) *
                       sin_theta +
                   2.0 * shift_delta * ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) *
                       (elongation_kappa * elongation_kappa - 4.0 * elongation_kappa * pow(sin_theta, 2.0) +
                        2.0 * elongation_kappa + 1.0)) *
                  ((-2.0) * (r / Rmax) * sin(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta) *
                       cos(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) +
                   M_PI * (1.0 - (r / Rmax) * (r / Rmax)) * (2.0 * elongation_kappa + 2.0) * sin_theta *
                       cos(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
                       cos(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta) -
                   M_PI * (1.0 - (r / Rmax) * (r / Rmax)) *
                       ((-4.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta) *
                       sin(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
                       sin(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta)) *
                  (elongation_kappa * elongation_kappa - 4.0 * elongation_kappa * pow(sin_theta, 2.0) +
                   2.0 * elongation_kappa + 1.0) *
                  exp(-tanh(20.0 * (r / Rmax) - 14.0)) /
                  pow(((pow((elongation_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) +
                        pow(((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta), 2.0)) *
                           (elongation_kappa * elongation_kappa - 4.0 * elongation_kappa * pow(sin_theta, 2.0) +
                            2.0 * elongation_kappa + 1.0) -
                       pow(((elongation_kappa - 1.0) *
                                ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) *
                                sin_theta +
                            1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)),
                           2.0)),
                      (3.0 / 2.0)) +
              (pow((elongation_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) +
               pow(((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta), 2.0)) *
                  ((-4.0) * (M_PI * M_PI) * (r / Rmax) * (1.0 - (r / Rmax) * (r / Rmax)) *
                       pow((elongation_kappa + 1.0), 2.0) *
                       sin(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta) * pow(cos_theta, 2.0) *
                       cos(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) -
                   M_PI * M_PI * (r / Rmax) * (1.0 - (r / Rmax) * (r / Rmax)) *
                       pow((2.0 * elongation_kappa - 2.0), 2.0) * pow(sin_theta, 2.0) *
                       sin(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta) *
                       cos(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) -
                   2.0 * (M_PI * M_PI) * (r / Rmax) * (1.0 - (r / Rmax) * (r / Rmax)) * (2.0 * elongation_kappa - 2.0) *
                       (2.0 * elongation_kappa + 2.0) * sin_theta *
                       sin(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
                       cos_theta * cos(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta) -
                   M_PI * (1.0 - (r / Rmax) * (r / Rmax)) * (2.0 * elongation_kappa - 2.0) *
                       sin(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
                       sin(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta) * cos_theta -
                   M_PI * (1.0 - (r / Rmax) * (r / Rmax)) * (2.0 * elongation_kappa + 2.0) * sin_theta *
                       cos(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
                       cos(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta)) *
                  exp(-tanh(20.0 * (r / Rmax) - 14.0)) /
                  sqrt((pow((elongation_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) +
                        pow(((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta), 2.0)) *
                           (elongation_kappa * elongation_kappa - 4.0 * elongation_kappa * pow(sin_theta, 2.0) +
                            2.0 * elongation_kappa + 1.0) -
                       pow(((elongation_kappa - 1.0) *
                                ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) *
                                sin_theta +
                            1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)),
                           2.0)) +
              (pow((elongation_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) +
               pow(((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta), 2.0)) *
                  ((-M_PI) * (1.0 - (r / Rmax) * (r / Rmax)) * (2.0 * elongation_kappa - 2.0) * sin_theta *
                       sin(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
                       sin(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta) +
                   M_PI * (1.0 - (r / Rmax) * (r / Rmax)) * (2.0 * elongation_kappa + 2.0) * cos_theta *
                       cos(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
                       cos(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta)) *
                  (4.0 * elongation_kappa *
                       (pow((elongation_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) +
                        pow(((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta), 2.0)) *
                       sin_theta * cos_theta -
                   1.0 / 2.0 *
                       (pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta) +
                        (2.0 * elongation_kappa - 2.0) *
                            ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) *
                            sin_theta) *
                       (elongation_kappa * elongation_kappa - 4.0 * elongation_kappa * pow(sin_theta, 2.0) +
                        2.0 * elongation_kappa + 1.0) +
                   1.0 / 2.0 *
                       ((elongation_kappa - 1.0) *
                            ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) * sin_theta +
                        1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)) *
                       (2.0 * pow((elongation_kappa - 1.0), 2.0) * pow(sin_theta, 2.0) +
                        2.0 * (elongation_kappa - 1.0) *
                            ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) * cos_theta +
                        2.0 * pow((elongation_kappa + 1.0), 2.0) * cos(2.0 * theta))) *
                  exp(-tanh(20.0 * (r / Rmax) - 14.0)) /
                  pow(((pow((elongation_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) +
                        pow(((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta), 2.0)) *
                           (elongation_kappa * elongation_kappa - 4.0 * elongation_kappa * pow(sin_theta, 2.0) +
                            2.0 * elongation_kappa + 1.0) -
                       pow(((elongation_kappa - 1.0) *
                                ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) *
                                sin_theta +
                            1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)),
                           2.0)),
                      (3.0 / 2.0)) +
              (pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta) +
               (2.0 * elongation_kappa - 2.0) *
                   ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) * sin_theta) *
                  ((-M_PI) * (1.0 - (r / Rmax) * (r / Rmax)) * (2.0 * elongation_kappa - 2.0) * sin_theta *
                       sin(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
                       sin(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta) +
                   M_PI * (1.0 - (r / Rmax) * (r / Rmax)) * (2.0 * elongation_kappa + 2.0) * cos_theta *
                       cos(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
                       cos(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta)) *
                  exp(-tanh(20.0 * (r / Rmax) - 14.0)) /
                  sqrt((pow((elongation_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) +
                        pow(((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta), 2.0)) *
                           (elongation_kappa * elongation_kappa - 4.0 * elongation_kappa * pow(sin_theta, 2.0) +
                            2.0 * elongation_kappa + 1.0) -
                       pow(((elongation_kappa - 1.0) *
                                ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) *
                                sin_theta +
                            1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)),
                           2.0)) -
              ((elongation_kappa - 1.0) *
                   ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) * sin_theta +
               1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)) *
                  ((-M_PI) * (1.0 - (r / Rmax) * (r / Rmax)) * (2.0 * elongation_kappa - 2.0) * sin_theta *
                       sin(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
                       sin(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta) +
                   M_PI * (1.0 - (r / Rmax) * (r / Rmax)) * (2.0 * elongation_kappa + 2.0) * cos_theta *
                       cos(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
                       cos(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta)) *
                  exp(-tanh(20.0 * (r / Rmax) - 14.0)) /
                  sqrt((pow((elongation_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) +
                        pow(((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta), 2.0)) *
                           (elongation_kappa * elongation_kappa - 4.0 * elongation_kappa * pow(sin_theta, 2.0) +
                            2.0 * elongation_kappa + 1.0) -
                       pow(((elongation_kappa - 1.0) *
                                ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) *
                                sin_theta +
                            1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)),
                           2.0)) -
              ((elongation_kappa - 1.0) *
                   ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) * sin_theta +
               1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)) *
                  (2.0 * M_PI * ((r / Rmax) * (r / Rmax)) * (2.0 * elongation_kappa - 2.0) * sin_theta *
                       sin(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
                       sin(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta) -
                   2.0 * M_PI * ((r / Rmax) * (r / Rmax)) * (2.0 * elongation_kappa + 2.0) * cos_theta *
                       cos(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
                       cos(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta) -
                   4.0 * (M_PI * M_PI) * (r / Rmax) * (1.0 - (r / Rmax) * (r / Rmax)) *
                       pow((elongation_kappa + 1.0), 2.0) * sin_theta *
                       sin(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta) * cos_theta *
                       cos(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) -
                   M_PI * M_PI * (r / Rmax) * (1.0 - (r / Rmax) * (r / Rmax)) * (2.0 * elongation_kappa - 2.0) *
                       (2.0 * elongation_kappa + 2.0) * pow(sin_theta, 2.0) *
                       sin(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
                       cos(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta) -
                   M_PI * M_PI * (r / Rmax) * (1.0 - (r / Rmax) * (r / Rmax)) * (2.0 * elongation_kappa - 2.0) *
                       ((-4.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta) *
                       sin_theta * sin(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta) *
                       cos(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) -
                   M_PI * M_PI * (r / Rmax) * (1.0 - (r / Rmax) * (r / Rmax)) * (2.0 * elongation_kappa + 2.0) *
                       ((-4.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta) *
                       sin(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
                       cos_theta * cos(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta) -
                   M_PI * (1.0 - (r / Rmax) * (r / Rmax)) * (2.0 * elongation_kappa - 2.0) * sin_theta *
                       sin(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
                       sin(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta) +
                   M_PI * (1.0 - (r / Rmax) * (r / Rmax)) * (2.0 * elongation_kappa + 2.0) * cos_theta *
                       cos(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
                       cos(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta)) *
                  exp(-tanh(20.0 * (r / Rmax) - 14.0)) /
                  sqrt((pow((elongation_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) +
                        pow(((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta), 2.0)) *
                           (elongation_kappa * elongation_kappa - 4.0 * elongation_kappa * pow(sin_theta, 2.0) +
                            2.0 * elongation_kappa + 1.0) -
                       pow(((elongation_kappa - 1.0) *
                                ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) *
                                sin_theta +
                            1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)),
                           2.0)) -
              (pow((elongation_kappa - 1.0), 2.0) * pow(sin_theta, 2.0) +
               (elongation_kappa - 1.0) *
                   ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) * cos_theta +
               pow((elongation_kappa + 1.0), 2.0) * cos(2.0 * theta)) *
                  ((-2.0) * (r / Rmax) * sin(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta) *
                       cos(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) +
                   M_PI * (1.0 - (r / Rmax) * (r / Rmax)) * (2.0 * elongation_kappa + 2.0) * sin_theta *
                       cos(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
                       cos(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta) -
                   M_PI * (1.0 - (r / Rmax) * (r / Rmax)) *
                       ((-4.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta) *
                       sin(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
                       sin(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta)) *
                  exp(-tanh(20.0 * (r / Rmax) - 14.0)) /
                  sqrt((pow((elongation_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) +
                        pow(((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta), 2.0)) *
                           (elongation_kappa * elongation_kappa - 4.0 * elongation_kappa * pow(sin_theta, 2.0) +
                            2.0 * elongation_kappa + 1.0) -
                       pow(((elongation_kappa - 1.0) *
                                ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) *
                                sin_theta +
                            1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)),
                           2.0)) +
              ((-2.0) * (r / Rmax) * sin(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta) *
                   cos(M_PI * (r / Rmax) *
                       ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) +
               M_PI * (1.0 - (r / Rmax) * (r / Rmax)) * (2.0 * elongation_kappa + 2.0) * sin_theta *
                   cos(M_PI * (r / Rmax) *
                       ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
                   cos(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta) -
               M_PI * (1.0 - (r / Rmax) * (r / Rmax)) *
                   ((-4.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta) *
                   sin(M_PI * (r / Rmax) *
                       ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
                   sin(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta)) *
                  (elongation_kappa * elongation_kappa - 4.0 * elongation_kappa * pow(sin_theta, 2.0) +
                   2.0 * elongation_kappa + 1.0) *
                  exp(-tanh(20.0 * (r / Rmax) - 14.0)) /
                  sqrt((pow((elongation_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) +
                        pow(((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta), 2.0)) *
                           (elongation_kappa * elongation_kappa - 4.0 * elongation_kappa * pow(sin_theta, 2.0) +
                            2.0 * elongation_kappa + 1.0) -
                       pow(((elongation_kappa - 1.0) *
                                ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) *
                                sin_theta +
                            1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)),
                           2.0)) -
              ((elongation_kappa - 1.0) *
                   ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) * sin_theta +
               1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)) *
                  ((-2.0) * (r / Rmax) * sin(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta) *
                       cos(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) +
                   M_PI * (1.0 - (r / Rmax) * (r / Rmax)) * (2.0 * elongation_kappa + 2.0) * sin_theta *
                       cos(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
                       cos(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta) -
                   M_PI * (1.0 - (r / Rmax) * (r / Rmax)) *
                       ((-4.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta) *
                       sin(M_PI * (r / Rmax) *
                           ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
                       sin(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta)) *
                  (4.0 * elongation_kappa *
                       (pow((elongation_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) +
                        pow(((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta), 2.0)) *
                       sin_theta * cos_theta -
                   1.0 / 2.0 *
                       (pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta) +
                        (2.0 * elongation_kappa - 2.0) *
                            ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) *
                            sin_theta) *
                       (elongation_kappa * elongation_kappa - 4.0 * elongation_kappa * pow(sin_theta, 2.0) +
                        2.0 * elongation_kappa + 1.0) +
                   1.0 / 2.0 *
                       ((elongation_kappa - 1.0) *
                            ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) * sin_theta +
                        1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)) *
                       (2.0 * pow((elongation_kappa - 1.0), 2.0) * pow(sin_theta, 2.0) +
                        2.0 * (elongation_kappa - 1.0) *
                            ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) * cos_theta +
                        2.0 * pow((elongation_kappa + 1.0), 2.0) * cos(2.0 * theta))) *
                  exp(-tanh(20.0 * (r / Rmax) - 14.0)) /
                  pow(((pow((elongation_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) +
                        pow(((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta), 2.0)) *
                           (elongation_kappa * elongation_kappa - 4.0 * elongation_kappa * pow(sin_theta, 2.0) +
                            2.0 * elongation_kappa + 1.0) -
                       pow(((elongation_kappa - 1.0) *
                                ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) *
                                sin_theta +
                            1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)),
                           2.0)),
                      (3.0 / 2.0)))) /
           ((r / Rmax) *
            sqrt((pow((elongation_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) +
                  pow(((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta), 2.0)) *
                     (elongation_kappa * elongation_kappa - 4.0 * elongation_kappa * pow(sin_theta, 2.0) +
                      2.0 * elongation_kappa + 1.0) -
                 pow(((elongation_kappa - 1.0) *
                          ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) * sin_theta +
                      1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)),
                     2.0)));
}