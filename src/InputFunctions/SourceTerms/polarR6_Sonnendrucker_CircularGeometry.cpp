#include "../include/InputFunctions/SourceTerms/polarR6_Sonnendrucker_CircularGeometry.h"

PolarR6_Sonnendrucker_CircularGeometry::PolarR6_Sonnendrucker_CircularGeometry(PolarGrid const& grid, double Rmax)
    : grid_(grid) , Rmax(Rmax)
{
}

double PolarR6_Sonnendrucker_CircularGeometry::operator()(double r, double theta) const
{
    return (-pow((r / Rmax), 4.0)) *
           ((r / Rmax) *
                (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r / Rmax) - 11.1111111111111)) *
                (12.288 * (r / Rmax) * pow(((r / Rmax) - 1.0), 4.0) * cos(11.0 * theta) +
                 17.2032 * pow(((r / Rmax) - 1.0), 5.0) * cos(11.0 * theta)) -
            5.03290747193186 * (r / Rmax) *
                (2.4576 * (r / Rmax) * pow(((r / Rmax) - 1.0), 5.0) * cos(11.0 * theta) +
                 2.4576 * pow(((r / Rmax) - 1.0), 6.0) * cos(11.0 * theta)) /
                (208.641975308642 * pow(((r / Rmax) - 0.769230769230769), 2.0) + 1.0) -
            49.5616 * (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r / Rmax) - 11.1111111111111)) *
                pow(((r / Rmax) - 1.0), 6.0) * cos(11.0 * theta) +
            6.0 * (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r / Rmax) - 11.1111111111111)) *
                (2.4576 * (r / Rmax) * pow(((r / Rmax) - 1.0), 5.0) * cos(11.0 * theta) +
                 2.4576 * pow(((r / Rmax) - 1.0), 6.0) * cos(11.0 * theta)));
}
