#include "../include/InputFunctions/SourceTerms/refined_ZoniShiftedGyro_ShafranovGeometry.h"

Refined_ZoniShiftedGyro_ShafranovGeometry::Refined_ZoniShiftedGyro_ShafranovGeometry(PolarGrid const& grid, double Rmax,
                                                                                     double elongation_kappa,
                                                                                     double shift_delta)
    : grid_(grid)
    , Rmax(Rmax)
    , elongation_kappa(elongation_kappa)
    , shift_delta(shift_delta)
{
}

double Refined_ZoniShiftedGyro_ShafranovGeometry::operator()(std::size_t i_r, std::size_t i_theta) const
{
    double r         = grid_.radius(i_r);
    double theta     = grid_.theta(i_theta);
    double sin_theta = std::sin(theta);
    double cos_theta = std::cos(theta);
    return 1.0 *
               (((-3.33823779536505e-15) * ((r / Rmax) * (r / Rmax)) - 0.0 * (r / Rmax) - 0.0 +
                 exp((-3333.33333333333) * pow(((r / Rmax) - 0.9), 2.0))) *
                    cos(21.0 * theta) +
                (0.00184273372222541 * ((r / Rmax) * (r / Rmax)) - 0.0018029383826828 * (r / Rmax) -
                 4.00652973929511e-05 + exp((-50.0) * pow(((r / Rmax) - 0.45), 2.0))) *
                    cos(9.0 * theta)) *
               exp(tanh(20.0 * (r / Rmax) - 14.0)) -
           (2.0 * shift_delta * (elongation_kappa - 1.0) *
                ((-21.0) *
                     ((-3.33823779536505e-15) * ((r / Rmax) * (r / Rmax)) +
                      exp((-3333.33333333333) * pow(((r / Rmax) - 0.9), 2.0))) *
                     sin(21.0 * theta) -
                 9.0 *
                     (0.00184273372222541 * ((r / Rmax) * (r / Rmax)) - 0.0018029383826828 * (r / Rmax) -
                      4.00652973929511e-05 + exp((-50.0) * pow(((r / Rmax) - 0.45), 2.0))) *
                     sin(9.0 * theta)) *
                exp(-tanh(20.0 * (r / Rmax) - 14.0)) * sin_theta /
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
                (((-6.67647559073009e-15) * (r / Rmax) +
                  (6000.0 - 6666.66666666667 * (r / Rmax)) * exp((-3333.33333333333) * pow(((r / Rmax) - 0.9), 2.0))) *
                     cos(21.0 * theta) +
                 (0.00368546744445083 * (r / Rmax) +
                  (45.0 - 100.0 * (r / Rmax)) * exp((-50.0) * pow(((r / Rmax) - 0.45), 2.0)) - 0.0018029383826828) *
                     cos(9.0 * theta)) *
                (20.0 * pow(tanh(20.0 * (r / Rmax) - 14.0), 2.0) - 20.0) *
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
                (((-6.67647559073009e-15) * (r / Rmax) +
                  (6000.0 - 6666.66666666667 * (r / Rmax)) * exp((-3333.33333333333) * pow(((r / Rmax) - 0.9), 2.0))) *
                     cos(21.0 * theta) +
                 (0.00368546744445083 * (r / Rmax) +
                  (45.0 - 100.0 * (r / Rmax)) * exp((-50.0) * pow(((r / Rmax) - 0.45), 2.0)) - 0.0018029383826828) *
                     cos(9.0 * theta)) *
                ((-2.0) * shift_delta * (elongation_kappa - 1.0) *
                     ((elongation_kappa - 1.0) *
                          ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) * sin_theta +
                      1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)) *
                     sin_theta +
                 2.0 * shift_delta * ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) *
                     (elongation_kappa * elongation_kappa - 4.0 * elongation_kappa * pow(sin_theta, 2.0) +
                      2.0 * elongation_kappa + 1.0)) *
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
            (r / Rmax) *
                ((10000.0 * pow((0.45 - (r / Rmax)), 2.0) * exp((-50.0) * pow(((r / Rmax) - 0.45), 2.0)) +
                  0.00368546744445083 - 100.0 * exp((-50.0) * pow(((r / Rmax) - 0.45), 2.0))) *
                     cos(9.0 * theta) +
                 (44444444.4444444 * pow((0.9 - (r / Rmax)), 2.0) *
                      exp((-3333.33333333333) * pow(((r / Rmax) - 0.9), 2.0)) -
                  6.67647559073009e-15 - 6666.66666666667 * exp((-3333.33333333333) * pow(((r / Rmax) - 0.9), 2.0))) *
                     cos(21.0 * theta)) *
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
            ((-21.0) *
                 ((-6.67647559073009e-15) * (r / Rmax) +
                  (6000.0 - 6666.66666666667 * (r / Rmax)) * exp((-3333.33333333333) * pow(((r / Rmax) - 0.9), 2.0))) *
                 sin(21.0 * theta) -
             9.0 *
                 (0.00368546744445083 * (r / Rmax) +
                  (45.0 - 100.0 * (r / Rmax)) * exp((-50.0) * pow(((r / Rmax) - 0.45), 2.0)) - 0.0018029383826828) *
                 sin(9.0 * theta)) *
                ((elongation_kappa - 1.0) *
                     ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) * sin_theta +
                 1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)) *
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
            (((-6.67647559073009e-15) * (r / Rmax) +
              (6000.0 - 6666.66666666667 * (r / Rmax)) * exp((-3333.33333333333) * pow(((r / Rmax) - 0.9), 2.0))) *
                 cos(21.0 * theta) +
             (0.00368546744445083 * (r / Rmax) +
              (45.0 - 100.0 * (r / Rmax)) * exp((-50.0) * pow(((r / Rmax) - 0.45), 2.0)) - 0.0018029383826828) *
                 cos(9.0 * theta)) *
                (pow((elongation_kappa - 1.0), 2.0) * pow(sin_theta, 2.0) +
                 (elongation_kappa - 1.0) *
                     ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) * cos_theta +
                 pow((elongation_kappa + 1.0), 2.0) * cos(2.0 * theta)) *
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
            (((-6.67647559073009e-15) * (r / Rmax) +
              (6000.0 - 6666.66666666667 * (r / Rmax)) * exp((-3333.33333333333) * pow(((r / Rmax) - 0.9), 2.0))) *
                 cos(21.0 * theta) +
             (0.00368546744445083 * (r / Rmax) +
              (45.0 - 100.0 * (r / Rmax)) * exp((-50.0) * pow(((r / Rmax) - 0.45), 2.0)) - 0.0018029383826828) *
                 cos(9.0 * theta)) *
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
            (((-6.67647559073009e-15) * (r / Rmax) +
              (6000.0 - 6666.66666666667 * (r / Rmax)) * exp((-3333.33333333333) * pow(((r / Rmax) - 0.9), 2.0))) *
                 cos(21.0 * theta) +
             (0.00368546744445083 * (r / Rmax) +
              (45.0 - 100.0 * (r / Rmax)) * exp((-50.0) * pow(((r / Rmax) - 0.45), 2.0)) - 0.0018029383826828) *
                 cos(9.0 * theta)) *
                ((elongation_kappa - 1.0) *
                     ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) * sin_theta +
                 1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)) *
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
                    (3.0 / 2.0)) -
            ((1.40205987405332e-13 * (r / Rmax) - 21.0 * (6000.0 - 6666.66666666667 * (r / Rmax)) *
                                                      exp((-3333.33333333333) * pow(((r / Rmax) - 0.9), 2.0))) *
                 sin(21.0 * theta) +
             ((-0.0331692070000574) * (r / Rmax) -
              9.0 * (45.0 - 100.0 * (r / Rmax)) * exp((-50.0) * pow(((r / Rmax) - 0.45), 2.0)) + 0.0162264454441452) *
                 sin(9.0 * theta)) *
                ((elongation_kappa - 1.0) *
                     ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) * sin_theta +
                 1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)) *
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
            ((-21.0) *
                 ((-3.33823779536505e-15) * ((r / Rmax) * (r / Rmax)) +
                  exp((-3333.33333333333) * pow(((r / Rmax) - 0.9), 2.0))) *
                 sin(21.0 * theta) -
             9.0 *
                 (0.00184273372222541 * ((r / Rmax) * (r / Rmax)) - 0.0018029383826828 * (r / Rmax) -
                  4.00652973929511e-05 + exp((-50.0) * pow(((r / Rmax) - 0.45), 2.0))) *
                 sin(9.0 * theta)) *
                ((elongation_kappa - 1.0) *
                     ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) * sin_theta +
                 1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)) *
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
            ((-21.0) *
                 ((-3.33823779536505e-15) * ((r / Rmax) * (r / Rmax)) +
                  exp((-3333.33333333333) * pow(((r / Rmax) - 0.9), 2.0))) *
                 sin(21.0 * theta) -
             9.0 *
                 (0.00184273372222541 * ((r / Rmax) * (r / Rmax)) - 0.0018029383826828 * (r / Rmax) -
                  4.00652973929511e-05 + exp((-50.0) * pow(((r / Rmax) - 0.45), 2.0))) *
                 sin(9.0 * theta)) *
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
                ((-21.0) *
                     ((-3.33823779536505e-15) * ((r / Rmax) * (r / Rmax)) +
                      exp((-3333.33333333333) * pow(((r / Rmax) - 0.9), 2.0))) *
                     sin(21.0 * theta) -
                 9.0 *
                     (0.00184273372222541 * ((r / Rmax) * (r / Rmax)) - 0.0018029383826828 * (r / Rmax) -
                      4.00652973929511e-05 + exp((-50.0) * pow(((r / Rmax) - 0.45), 2.0))) *
                     sin(9.0 * theta)) *
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
                ((r / Rmax) *
                 pow(((pow((elongation_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) +
                       pow(((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta), 2.0)) *
                          (elongation_kappa * elongation_kappa - 4.0 * elongation_kappa * pow(sin_theta, 2.0) +
                           2.0 * elongation_kappa + 1.0) -
                      pow(((elongation_kappa - 1.0) *
                               ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) *
                               sin_theta +
                           1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)),
                          2.0)),
                     (3.0 / 2.0))) +
            (pow((elongation_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) +
             pow(((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta), 2.0)) *
                (21.0 *
                     (7.0102993702666e-14 * ((r / Rmax) * (r / Rmax)) -
                      21.0 * exp((-3333.33333333333) * pow(((r / Rmax) - 0.9), 2.0))) *
                     cos(21.0 * theta) +
                 9.0 *
                     ((-0.0165846035000287) * ((r / Rmax) * (r / Rmax)) + 0.0162264454441452 * (r / Rmax) +
                      0.00036058767653656 - 9.0 * exp((-50.0) * pow(((r / Rmax) - 0.45), 2.0))) *
                     cos(9.0 * theta)) *
                exp(-tanh(20.0 * (r / Rmax) - 14.0)) /
                ((r / Rmax) *
                 sqrt((pow((elongation_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) +
                       pow(((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta), 2.0)) *
                          (elongation_kappa * elongation_kappa - 4.0 * elongation_kappa * pow(sin_theta, 2.0) +
                           2.0 * elongation_kappa + 1.0) -
                      pow(((elongation_kappa - 1.0) *
                               ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) *
                               sin_theta +
                           1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)),
                          2.0))) +
            (pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta) +
             (2.0 * elongation_kappa - 2.0) *
                 ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) * sin_theta) *
                ((-21.0) *
                     ((-3.33823779536505e-15) * ((r / Rmax) * (r / Rmax)) +
                      exp((-3333.33333333333) * pow(((r / Rmax) - 0.9), 2.0))) *
                     sin(21.0 * theta) -
                 9.0 *
                     (0.00184273372222541 * ((r / Rmax) * (r / Rmax)) - 0.0018029383826828 * (r / Rmax) -
                      4.00652973929511e-05 + exp((-50.0) * pow(((r / Rmax) - 0.45), 2.0))) *
                     sin(9.0 * theta)) *
                exp(-tanh(20.0 * (r / Rmax) - 14.0)) /
                ((r / Rmax) *
                 sqrt((pow((elongation_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) +
                       pow(((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta), 2.0)) *
                          (elongation_kappa * elongation_kappa - 4.0 * elongation_kappa * pow(sin_theta, 2.0) +
                           2.0 * elongation_kappa + 1.0) -
                      pow(((elongation_kappa - 1.0) *
                               ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) *
                               sin_theta +
                           1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)),
                          2.0)))) /
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