#include "../include/InputFunctions/SourceTerms/refined_ZoniShiftedGyro_CircularGeometry.h"

Refined_ZoniShiftedGyro_CircularGeometry::Refined_ZoniShiftedGyro_CircularGeometry(PolarGrid const& grid, double Rmax)
    : grid_(grid) , Rmax(Rmax)
{
}

double Refined_ZoniShiftedGyro_CircularGeometry::operator()(double r, double theta) const
{
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
                (20.0 * pow(tanh(20.0 * (r / Rmax) - 14.0), 2.0) - 20.0) * exp(-tanh(20.0 * (r / Rmax) - 14.0)) +
            (r / Rmax) *
                ((10000.0 * pow((0.45 - (r / Rmax)), 2.0) * exp((-50.0) * pow(((r / Rmax) - 0.45), 2.0)) +
                  0.00368546744445083 - 100.0 * exp((-50.0) * pow(((r / Rmax) - 0.45), 2.0))) *
                     cos(9.0 * theta) +
                 (44444444.4444444 * pow((0.9 - (r / Rmax)), 2.0) *
                      exp((-3333.33333333333) * pow(((r / Rmax) - 0.9), 2.0)) -
                  6.67647559073009e-15 - 6666.66666666667 * exp((-3333.33333333333) * pow(((r / Rmax) - 0.9), 2.0))) *
                     cos(21.0 * theta)) *
                exp(-tanh(20.0 * (r / Rmax) - 14.0)) +
            (((-6.67647559073009e-15) * (r / Rmax) +
              (6000.0 - 6666.66666666667 * (r / Rmax)) * exp((-3333.33333333333) * pow(((r / Rmax) - 0.9), 2.0))) *
                 cos(21.0 * theta) +
             (0.00368546744445083 * (r / Rmax) +
              (45.0 - 100.0 * (r / Rmax)) * exp((-50.0) * pow(((r / Rmax) - 0.45), 2.0)) - 0.0018029383826828) *
                 cos(9.0 * theta)) *
                exp(-tanh(20.0 * (r / Rmax) - 14.0)) +
            (21.0 *
                 (7.0102993702666e-14 * ((r / Rmax) * (r / Rmax)) -
                  21.0 * exp((-3333.33333333333) * pow(((r / Rmax) - 0.9), 2.0))) *
                 cos(21.0 * theta) +
             9.0 *
                 ((-0.0165846035000287) * ((r / Rmax) * (r / Rmax)) + 0.0162264454441452 * (r / Rmax) +
                  0.00036058767653656 - 9.0 * exp((-50.0) * pow(((r / Rmax) - 0.45), 2.0))) *
                 cos(9.0 * theta)) *
                exp(-tanh(20.0 * (r / Rmax) - 14.0)) / (r / Rmax)) /
               (r / Rmax);
}
