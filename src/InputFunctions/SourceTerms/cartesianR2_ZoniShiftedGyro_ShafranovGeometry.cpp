#include "../include/InputFunctions/SourceTerms/cartesianR2_ZoniShiftedGyro_ShafranovGeometry.h"

CartesianR2_ZoniShiftedGyro_ShafranovGeometry::CartesianR2_ZoniShiftedGyro_ShafranovGeometry(PolarGrid const& grid,
                                                                                             double Rmax,
                                                                                             double elongation_kappa,
                                                                                             double shift_delta)
    : grid_(grid)
    , Rmax(Rmax)
    , elongation_kappa(elongation_kappa)
    , shift_delta(shift_delta)
{
}

double CartesianR2_ZoniShiftedGyro_ShafranovGeometry::operator()(int i_r, int i_theta) const
{
    double r         = grid_.radius(i_r);
    double theta     = grid_.theta(i_theta);
    double sin_theta = std::sin(theta);
    double cos_theta = std::cos(theta);
    return (1.0 - (r / Rmax) * (r / Rmax)) * exp(tanh(20.0 * (r / Rmax) - 14.0)) *
               sin(M_PI * (2.0 * elongation_kappa * (r / Rmax) * sin_theta + 2.0 * (r / Rmax) * sin_theta)) *
               cos(M_PI * ((-2.0) * shift_delta * ((r / Rmax) * (r / Rmax)) -
                           2.0 * elongation_kappa * (r / Rmax) * cos_theta + 2.0 * (r / Rmax) * cos_theta)) -
           (2.0 * shift_delta * (r / Rmax) * (elongation_kappa - 1.0) *
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
                      M_PI *
                          ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
                     sin_theta * sin(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta) *
                     cos(M_PI * (r / Rmax) *
                         ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) -
                 M_PI * (1.0 - (r / Rmax) * (r / Rmax)) * (2.0 * elongation_kappa + 2.0) *
                     ((-2.0) * M_PI * shift_delta * (r / Rmax) +
                      M_PI *
                          ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
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
                      M_PI *
                          ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
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
                      M_PI *
                          ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
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
                      M_PI *
                          ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) *
                     ((-4.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta) *
                     sin(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta) *
                     cos(M_PI * (r / Rmax) *
                         ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta)) -
                 2.0 * sin(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta) *
                     cos(M_PI * (r / Rmax) *
                         ((-2.0) * shift_delta * (r / Rmax) - 2.0 * elongation_kappa * cos_theta + 2.0 * cos_theta))) *
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
                 M_PI * M_PI * (r / Rmax) * (1.0 - (r / Rmax) * (r / Rmax)) * pow((2.0 * elongation_kappa - 2.0), 2.0) *
                     pow(sin_theta, 2.0) * sin(M_PI * (r / Rmax) * (2.0 * elongation_kappa + 2.0) * sin_theta) *
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
                          ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) * sin_theta) *
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
            ((elongation_kappa - 1.0) * ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) *
                 sin_theta +
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
            ((elongation_kappa - 1.0) * ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) *
                 sin_theta +
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
             (elongation_kappa - 1.0) * ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) *
                 cos_theta +
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
            ((elongation_kappa - 1.0) * ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) *
                 sin_theta +
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
                          ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) * sin_theta) *
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
                    (3.0 / 2.0))) /
               ((r / Rmax) *
                sqrt((pow((elongation_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) +
                      pow(((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta), 2.0)) *
                         (elongation_kappa * elongation_kappa - 4.0 * elongation_kappa * pow(sin_theta, 2.0) +
                          2.0 * elongation_kappa + 1.0) -
                     pow(((elongation_kappa - 1.0) *
                              ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) *
                              sin_theta +
                          1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)),
                         2.0)));
}