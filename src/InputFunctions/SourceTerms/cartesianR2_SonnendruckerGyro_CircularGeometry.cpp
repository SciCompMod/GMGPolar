#include "../include/InputFunctions/SourceTerms/cartesianR2_SonnendruckerGyro_CircularGeometry.h"

CartesianR2_SonnendruckerGyro_CircularGeometry::CartesianR2_SonnendruckerGyro_CircularGeometry(const double& Rmax)
    : Rmax(Rmax)
{
}

double CartesianR2_SonnendruckerGyro_CircularGeometry::rhs_f(const double& r, const double& theta)const
{
    double sin_theta = std::sin(theta);
    double cos_theta = std::cos(theta);
    return (1.0 - (r / Rmax) * (r / Rmax)) * sin(2.0 * M_PI * (r / Rmax) * sin_theta) *
               cos(2.0 * M_PI * (r / Rmax) * cos_theta) /
               (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r / Rmax) - 11.1111111111111)) -
           ((r / Rmax) *
                (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r / Rmax) - 11.1111111111111)) *
                ((-8.0) * M_PI * (r / Rmax) * sin(theta) * cos(2.0 * M_PI * (r / Rmax) * sin(theta)) *
                     cos(2.0 * M_PI * (r / Rmax) * cos(theta)) +
                 8.0 * M_PI * (r / Rmax) * sin(2.0 * M_PI * (r / Rmax) * sin(theta)) *
                     sin(2.0 * M_PI * (r / Rmax) * cos(theta)) * cos(theta) -
                 4.0 * (M_PI * M_PI) * (1.0 - (r / Rmax) * (r / Rmax)) * pow(sin(theta), 2.0) *
                     sin(2.0 * M_PI * (r / Rmax) * sin(theta)) * cos(2.0 * M_PI * (r / Rmax) * cos(theta)) -
                 8.0 * (M_PI * M_PI) * (1.0 - (r / Rmax) * (r / Rmax)) * sin(theta) *
                     sin(2.0 * M_PI * (r / Rmax) * cos(theta)) * cos(theta) *
                     cos(2.0 * M_PI * (r / Rmax) * sin(theta)) -
                 4.0 * (M_PI * M_PI) * (1.0 - (r / Rmax) * (r / Rmax)) * sin(2.0 * M_PI * (r / Rmax) * sin(theta)) *
                     pow(cos(theta), 2.0) * cos(2.0 * M_PI * (r / Rmax) * cos(theta)) -
                 2.0 * sin(2.0 * M_PI * (r / Rmax) * sin(theta)) * cos(2.0 * M_PI * (r / Rmax) * cos(theta))) -
            5.03290747193186 * (r / Rmax) *
                ((-2.0) * (r / Rmax) * sin(2.0 * M_PI * (r / Rmax) * sin(theta)) *
                     cos(2.0 * M_PI * (r / Rmax) * cos(theta)) +
                 2.0 * M_PI * (1.0 - (r / Rmax) * (r / Rmax)) * sin(theta) * cos(2.0 * M_PI * (r / Rmax) * sin(theta)) *
                     cos(2.0 * M_PI * (r / Rmax) * cos(theta)) -
                 2.0 * M_PI * (1.0 - (r / Rmax) * (r / Rmax)) * sin(2.0 * M_PI * (r / Rmax) * sin(theta)) *
                     sin(2.0 * M_PI * (r / Rmax) * cos(theta)) * cos(theta)) /
                (208.641975308642 * pow(((r / Rmax) - 0.769230769230769), 2.0) + 1.0) +
            (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r / Rmax) - 11.1111111111111)) *
                ((-2.0) * (r / Rmax) * sin(2.0 * M_PI * (r / Rmax) * sin(theta)) *
                     cos(2.0 * M_PI * (r / Rmax) * cos(theta)) +
                 2.0 * M_PI * (1.0 - (r / Rmax) * (r / Rmax)) * sin(theta) * cos(2.0 * M_PI * (r / Rmax) * sin(theta)) *
                     cos(2.0 * M_PI * (r / Rmax) * cos(theta)) -
                 2.0 * M_PI * (1.0 - (r / Rmax) * (r / Rmax)) * sin(2.0 * M_PI * (r / Rmax) * sin(theta)) *
                     sin(2.0 * M_PI * (r / Rmax) * cos(theta)) * cos(theta)) +
            (0.452961672473868 - 0.348432055749129 * atan(14.4444444444444 * (r / Rmax) - 11.1111111111111)) *
                ((-4.0) * (M_PI * M_PI) * (r / Rmax) * (1.0 - (r / Rmax) * (r / Rmax)) * pow(sin(theta), 2.0) *
                     sin(2.0 * M_PI * (r / Rmax) * sin(theta)) * cos(2.0 * M_PI * (r / Rmax) * cos(theta)) +
                 8.0 * (M_PI * M_PI) * (r / Rmax) * (1.0 - (r / Rmax) * (r / Rmax)) * sin(theta) *
                     sin(2.0 * M_PI * (r / Rmax) * cos(theta)) * cos(theta) *
                     cos(2.0 * M_PI * (r / Rmax) * sin(theta)) -
                 4.0 * (M_PI * M_PI) * (r / Rmax) * (1.0 - (r / Rmax) * (r / Rmax)) *
                     sin(2.0 * M_PI * (r / Rmax) * sin(theta)) * pow(cos(theta), 2.0) *
                     cos(2.0 * M_PI * (r / Rmax) * cos(theta)) -
                 2.0 * M_PI * (1.0 - (r / Rmax) * (r / Rmax)) * sin(theta) * cos(2.0 * M_PI * (r / Rmax) * sin(theta)) *
                     cos(2.0 * M_PI * (r / Rmax) * cos(theta)) +
                 2.0 * M_PI * (1.0 - (r / Rmax) * (r / Rmax)) * sin(2.0 * M_PI * (r / Rmax) * sin(theta)) *
                     sin(2.0 * M_PI * (r / Rmax) * cos(theta)) * cos(theta))) /
               (r / Rmax);
}
