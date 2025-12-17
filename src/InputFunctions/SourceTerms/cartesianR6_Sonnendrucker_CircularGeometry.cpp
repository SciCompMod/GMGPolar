#include "../include/InputFunctions/SourceTerms/cartesianR6_Sonnendrucker_CircularGeometry.h"

CartesianR6_Sonnendrucker_CircularGeometry::CartesianR6_Sonnendrucker_CircularGeometry(PolarGrid const& grid, double Rmax)
    : grid_(grid) , Rmax(Rmax)
{
}

double CartesianR6_Sonnendrucker_CircularGeometry::operator()(double r, double theta) const
{
    double sin_theta = std::sin(theta);
    double cos_theta = std::cos(theta);
    return (-((r / Rmax) *
                  (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r / Rmax) - 11.1111111111111)) *
                  ((-1.6384) * (M_PI * M_PI) * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 6.0) *
                       pow(sin_theta, 2.0) * sin(2.0 * M_PI * (r / Rmax) * sin_theta) *
                       cos(2.0 * M_PI * (r / Rmax) * cos_theta) -
                   3.2768 * (M_PI * M_PI) * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 6.0) * sin_theta *
                       sin(2.0 * M_PI * (r / Rmax) * cos_theta) * cos_theta * cos(2.0 * M_PI * (r / Rmax) * sin_theta) -
                   1.6384 * (M_PI * M_PI) * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 6.0) *
                       sin(2.0 * M_PI * (r / Rmax) * sin_theta) * pow(cos_theta, 2.0) *
                       cos(2.0 * M_PI * (r / Rmax) * cos_theta) +
                   9.8304 * M_PI * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 5.0) * sin_theta *
                       cos(2.0 * M_PI * (r / Rmax) * sin_theta) * cos(2.0 * M_PI * (r / Rmax) * cos_theta) -
                   9.8304 * M_PI * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 5.0) *
                       sin(2.0 * M_PI * (r / Rmax) * sin_theta) * sin(2.0 * M_PI * (r / Rmax) * cos_theta) * cos_theta +
                   12.288 * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 4.0) *
                       sin(2.0 * M_PI * (r / Rmax) * sin_theta) * cos(2.0 * M_PI * (r / Rmax) * cos_theta) +
                   9.8304 * M_PI * pow(((r / Rmax) - 1.0), 5.0) * pow(((r / Rmax) + 1.0), 6.0) * sin_theta *
                       cos(2.0 * M_PI * (r / Rmax) * sin_theta) * cos(2.0 * M_PI * (r / Rmax) * cos_theta) -
                   9.8304 * M_PI * pow(((r / Rmax) - 1.0), 5.0) * pow(((r / Rmax) + 1.0), 6.0) *
                       sin(2.0 * M_PI * (r / Rmax) * sin_theta) * sin(2.0 * M_PI * (r / Rmax) * cos_theta) * cos_theta +
                   29.4912 * pow(((r / Rmax) - 1.0), 5.0) * pow(((r / Rmax) + 1.0), 5.0) *
                       sin(2.0 * M_PI * (r / Rmax) * sin_theta) * cos(2.0 * M_PI * (r / Rmax) * cos_theta) +
                   12.288 * pow(((r / Rmax) - 1.0), 4.0) * pow(((r / Rmax) + 1.0), 6.0) *
                       sin(2.0 * M_PI * (r / Rmax) * sin_theta) * cos(2.0 * M_PI * (r / Rmax) * cos_theta)) -
              5.03290747193186 * (r / Rmax) *
                  (0.8192 * M_PI * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 6.0) * sin_theta *
                       cos(2.0 * M_PI * (r / Rmax) * sin_theta) * cos(2.0 * M_PI * (r / Rmax) * cos_theta) -
                   0.8192 * M_PI * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 6.0) *
                       sin(2.0 * M_PI * (r / Rmax) * sin_theta) * sin(2.0 * M_PI * (r / Rmax) * cos_theta) * cos_theta +
                   2.4576 * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 5.0) *
                       sin(2.0 * M_PI * (r / Rmax) * sin_theta) * cos(2.0 * M_PI * (r / Rmax) * cos_theta) +
                   2.4576 * pow(((r / Rmax) - 1.0), 5.0) * pow(((r / Rmax) + 1.0), 6.0) *
                       sin(2.0 * M_PI * (r / Rmax) * sin_theta) * cos(2.0 * M_PI * (r / Rmax) * cos_theta)) /
                  (208.641975308642 * pow(((r / Rmax) - 0.769230769230769), 2.0) + 1.0) +
              (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r / Rmax) - 11.1111111111111)) *
                  (0.8192 * M_PI * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 6.0) * sin_theta *
                       cos(2.0 * M_PI * (r / Rmax) * sin_theta) * cos(2.0 * M_PI * (r / Rmax) * cos_theta) -
                   0.8192 * M_PI * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 6.0) *
                       sin(2.0 * M_PI * (r / Rmax) * sin_theta) * sin(2.0 * M_PI * (r / Rmax) * cos_theta) * cos_theta +
                   2.4576 * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 5.0) *
                       sin(2.0 * M_PI * (r / Rmax) * sin_theta) * cos(2.0 * M_PI * (r / Rmax) * cos_theta) +
                   2.4576 * pow(((r / Rmax) - 1.0), 5.0) * pow(((r / Rmax) + 1.0), 6.0) *
                       sin(2.0 * M_PI * (r / Rmax) * sin_theta) * cos(2.0 * M_PI * (r / Rmax) * cos_theta)) +
              (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r / Rmax) - 11.1111111111111)) *
                  ((-1.6384) * (M_PI * M_PI) * (r / Rmax) * pow(((r / Rmax) - 1.0), 6.0) *
                       pow(((r / Rmax) + 1.0), 6.0) * pow(sin_theta, 2.0) * sin(2.0 * M_PI * (r / Rmax) * sin_theta) *
                       cos(2.0 * M_PI * (r / Rmax) * cos_theta) +
                   3.2768 * (M_PI * M_PI) * (r / Rmax) * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 6.0) *
                       sin_theta * sin(2.0 * M_PI * (r / Rmax) * cos_theta) * cos_theta *
                       cos(2.0 * M_PI * (r / Rmax) * sin_theta) -
                   1.6384 * (M_PI * M_PI) * (r / Rmax) * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 6.0) *
                       sin(2.0 * M_PI * (r / Rmax) * sin_theta) * pow(cos_theta, 2.0) *
                       cos(2.0 * M_PI * (r / Rmax) * cos_theta) -
                   0.8192 * M_PI * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 6.0) * sin_theta *
                       cos(2.0 * M_PI * (r / Rmax) * sin_theta) * cos(2.0 * M_PI * (r / Rmax) * cos_theta) +
                   0.8192 * M_PI * pow(((r / Rmax) - 1.0), 6.0) * pow(((r / Rmax) + 1.0), 6.0) *
                       sin(2.0 * M_PI * (r / Rmax) * sin_theta) * sin(2.0 * M_PI * (r / Rmax) * cos_theta) *
                       cos_theta))) /
           (r / Rmax);
}