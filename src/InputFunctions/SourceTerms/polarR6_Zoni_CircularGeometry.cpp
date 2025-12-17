#include "../include/InputFunctions/SourceTerms/polarR6_Zoni_CircularGeometry.h"

PolarR6_Zoni_CircularGeometry::PolarR6_Zoni_CircularGeometry(double Rmax)
    : Rmax(Rmax)
{
}

double PolarR6_Zoni_CircularGeometry::operator()(double r, double theta) const
{
    return (-pow((r / Rmax), 4.0)) *
           ((r / Rmax) *
                (12.288 * (r / Rmax) * pow(((r / Rmax) - 1.0), 4.0) * cos(11.0 * theta) +
                 17.2032 * pow(((r / Rmax) - 1.0), 5.0) * cos(11.0 * theta)) *
                exp(-tanh(10.0 * (r / Rmax) - 5.0)) +
            (r / Rmax) *
                (2.4576 * (r / Rmax) * pow(((r / Rmax) - 1.0), 5.0) * cos(11.0 * theta) +
                 2.4576 * pow(((r / Rmax) - 1.0), 6.0) * cos(11.0 * theta)) *
                (10.0 * pow(tanh(10.0 * (r / Rmax) - 5.0), 2.0) - 10.0) * exp(-tanh(10.0 * (r / Rmax) - 5.0)) -
            49.5616 * pow(((r / Rmax) - 1.0), 6.0) * exp(-tanh(10.0 * (r / Rmax) - 5.0)) * cos(11.0 * theta) +
            6.0 *
                (2.4576 * (r / Rmax) * pow(((r / Rmax) - 1.0), 5.0) * cos(11.0 * theta) +
                 2.4576 * pow(((r / Rmax) - 1.0), 6.0) * cos(11.0 * theta)) *
                exp(-tanh(10.0 * (r / Rmax) - 5.0)));
}
