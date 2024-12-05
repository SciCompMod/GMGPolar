#include "../include/InputFunctions/SourceTerms/polarR6_Poisson_ShafranovGeometry.h"

PolarR6_Poisson_ShafranovGeometry::PolarR6_Poisson_ShafranovGeometry(const double& Rmax, const double& elongation_kappa,
                                                                     const double& shift_delta)
    : Rmax(Rmax)
    , elongation_kappa(elongation_kappa)
    , shift_delta(shift_delta)
{
}

double PolarR6_Poisson_ShafranovGeometry::rhs_f(const double& r, const double& theta, const double& sin_theta,
                                                const double& cos_theta) const
{
    return (-pow((r / Rmax), 4.0)) *
           ((-9.0112) * shift_delta * (r / Rmax) * (elongation_kappa - 1.0) * pow(((r / Rmax) - 1.0), 6.0) * sin_theta *
                sin(11.0 * theta) /
                sqrt((pow((elongation_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) +
                      pow(((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta), 2.0)) *
                         (elongation_kappa * elongation_kappa - 4.0 * elongation_kappa * pow(sin_theta, 2.0) +
                          2.0 * elongation_kappa + 1.0) -
                     pow(((elongation_kappa - 1.0) *
                              ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) *
                              sin_theta +
                          1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)),
                         2.0)) +
            4.5056 * (r / Rmax) * pow(((r / Rmax) - 1.0), 6.0) *
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
                sin(11.0 * theta) /
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
            27.0336 * (r / Rmax) * pow(((r / Rmax) - 1.0), 5.0) *
                ((elongation_kappa - 1.0) *
                     ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) * sin_theta +
                 1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)) *
                sin(11.0 * theta) /
                sqrt((pow((elongation_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) +
                      pow(((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta), 2.0)) *
                         (elongation_kappa * elongation_kappa - 4.0 * elongation_kappa * pow(sin_theta, 2.0) +
                          2.0 * elongation_kappa + 1.0) -
                     pow(((elongation_kappa - 1.0) *
                              ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) *
                              sin_theta +
                          1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)),
                         2.0)) +
            1.0 * (r / Rmax) *
                (12.288 * (r / Rmax) * pow(((r / Rmax) - 1.0), 4.0) * cos(11.0 * theta) +
                 17.2032 * pow(((r / Rmax) - 1.0), 5.0) * cos(11.0 * theta)) *
                (elongation_kappa * elongation_kappa - 4.0 * elongation_kappa * pow(sin_theta, 2.0) +
                 2.0 * elongation_kappa + 1.0) /
                sqrt((pow((elongation_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) +
                      pow(((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta), 2.0)) *
                         (elongation_kappa * elongation_kappa - 4.0 * elongation_kappa * pow(sin_theta, 2.0) +
                          2.0 * elongation_kappa + 1.0) -
                     pow(((elongation_kappa - 1.0) *
                              ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) *
                              sin_theta +
                          1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)),
                         2.0)) +
            1.0 * (r / Rmax) *
                (2.4576 * (r / Rmax) * pow(((r / Rmax) - 1.0), 5.0) * cos(11.0 * theta) +
                 2.4576 * pow(((r / Rmax) - 1.0), 6.0) * cos(11.0 * theta)) *
                ((-2.0) * shift_delta * (elongation_kappa - 1.0) *
                     ((elongation_kappa - 1.0) *
                          ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) * sin_theta +
                      1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)) *
                     sin_theta +
                 2.0 * shift_delta * ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) *
                     (elongation_kappa * elongation_kappa - 4.0 * elongation_kappa * pow(sin_theta, 2.0) +
                      2.0 * elongation_kappa + 1.0)) *
                (elongation_kappa * elongation_kappa - 4.0 * elongation_kappa * pow(sin_theta, 2.0) +
                 2.0 * elongation_kappa + 1.0) /
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
            49.5616 * pow(((r / Rmax) - 1.0), 6.0) *
                (pow((elongation_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) +
                 pow(((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta), 2.0)) *
                cos(11.0 * theta) /
                sqrt((pow((elongation_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) +
                      pow(((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta), 2.0)) *
                         (elongation_kappa * elongation_kappa - 4.0 * elongation_kappa * pow(sin_theta, 2.0) +
                          2.0 * elongation_kappa + 1.0) -
                     pow(((elongation_kappa - 1.0) *
                              ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) *
                              sin_theta +
                          1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)),
                         2.0)) -
            4.5056 * pow(((r / Rmax) - 1.0), 6.0) *
                (pow((elongation_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) +
                 pow(((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta), 2.0)) *
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
                sin(11.0 * theta) /
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
            4.5056 * pow(((r / Rmax) - 1.0), 6.0) *
                (pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta) +
                 (2.0 * elongation_kappa - 2.0) *
                     ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) * sin_theta) *
                sin(11.0 * theta) /
                sqrt((pow((elongation_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) +
                      pow(((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta), 2.0)) *
                         (elongation_kappa * elongation_kappa - 4.0 * elongation_kappa * pow(sin_theta, 2.0) +
                          2.0 * elongation_kappa + 1.0) -
                     pow(((elongation_kappa - 1.0) *
                              ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) *
                              sin_theta +
                          1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)),
                         2.0)) +
            27.0336 * pow(((r / Rmax) - 1.0), 6.0) *
                ((elongation_kappa - 1.0) *
                     ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) * sin_theta +
                 1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)) *
                sin(11.0 * theta) /
                sqrt((pow((elongation_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) +
                      pow(((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta), 2.0)) *
                         (elongation_kappa * elongation_kappa - 4.0 * elongation_kappa * pow(sin_theta, 2.0) +
                          2.0 * elongation_kappa + 1.0) -
                     pow(((elongation_kappa - 1.0) *
                              ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) *
                              sin_theta +
                          1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)),
                         2.0)) -
            1.0 *
                ((-27.0336) * (r / Rmax) * pow(((r / Rmax) - 1.0), 5.0) * sin(11.0 * theta) -
                 27.0336 * pow(((r / Rmax) - 1.0), 6.0) * sin(11.0 * theta)) *
                ((elongation_kappa - 1.0) *
                     ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) * sin_theta +
                 1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)) /
                sqrt((pow((elongation_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) +
                      pow(((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta), 2.0)) *
                         (elongation_kappa * elongation_kappa - 4.0 * elongation_kappa * pow(sin_theta, 2.0) +
                          2.0 * elongation_kappa + 1.0) -
                     pow(((elongation_kappa - 1.0) *
                              ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) *
                              sin_theta +
                          1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)),
                         2.0)) -
            1.0 *
                (2.4576 * (r / Rmax) * pow(((r / Rmax) - 1.0), 5.0) * cos(11.0 * theta) +
                 2.4576 * pow(((r / Rmax) - 1.0), 6.0) * cos(11.0 * theta)) *
                (pow((elongation_kappa - 1.0), 2.0) * pow(sin_theta, 2.0) +
                 (elongation_kappa - 1.0) *
                     ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) * cos_theta +
                 pow((elongation_kappa + 1.0), 2.0) * cos(2.0 * theta)) /
                sqrt((pow((elongation_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) +
                      pow(((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta), 2.0)) *
                         (elongation_kappa * elongation_kappa - 4.0 * elongation_kappa * pow(sin_theta, 2.0) +
                          2.0 * elongation_kappa + 1.0) -
                     pow(((elongation_kappa - 1.0) *
                              ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) *
                              sin_theta +
                          1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)),
                         2.0)) +
            6.0 *
                (2.4576 * (r / Rmax) * pow(((r / Rmax) - 1.0), 5.0) * cos(11.0 * theta) +
                 2.4576 * pow(((r / Rmax) - 1.0), 6.0) * cos(11.0 * theta)) *
                (elongation_kappa * elongation_kappa - 4.0 * elongation_kappa * pow(sin_theta, 2.0) +
                 2.0 * elongation_kappa + 1.0) /
                sqrt((pow((elongation_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) +
                      pow(((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta), 2.0)) *
                         (elongation_kappa * elongation_kappa - 4.0 * elongation_kappa * pow(sin_theta, 2.0) +
                          2.0 * elongation_kappa + 1.0) -
                     pow(((elongation_kappa - 1.0) *
                              ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) *
                              sin_theta +
                          1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)),
                         2.0)) -
            1.0 *
                (2.4576 * (r / Rmax) * pow(((r / Rmax) - 1.0), 5.0) * cos(11.0 * theta) +
                 2.4576 * pow(((r / Rmax) - 1.0), 6.0) * cos(11.0 * theta)) *
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
                      2.0 * pow((elongation_kappa + 1.0), 2.0) * cos(2.0 * theta))) /
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
           sqrt((pow((elongation_kappa + 1.0), 2.0) * pow(sin_theta, 2.0) +
                 pow(((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta), 2.0)) *
                    (elongation_kappa * elongation_kappa - 4.0 * elongation_kappa * pow(sin_theta, 2.0) +
                     2.0 * elongation_kappa + 1.0) -
                pow(((elongation_kappa - 1.0) *
                         ((-2.0) * shift_delta * (r / Rmax) - elongation_kappa * cos_theta + cos_theta) * sin_theta +
                     1.0 / 2.0 * pow((elongation_kappa + 1.0), 2.0) * sin(2.0 * theta)),
                    2.0));
}