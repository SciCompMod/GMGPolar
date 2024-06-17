// CartesianR6 simulates solution (23) of Bourne et al. https://doi.org/10.1016/j.jcp.2023.112249
#include "../../include/test_cases/CartesianR6PoissonTriangular.h"

double CartesianR6PoissonTriangular::x(double r, double theta, double map2_epsilon, double map2_e, double Rmax) const
{
    return (1.0 - sqrt(map2_epsilon * (map2_epsilon + 2.0 * (r/Rmax) * cos(theta)) + 1.0)) / map2_epsilon;
}

double CartesianR6PoissonTriangular::y(double r, double theta, double map2_epsilon, double map2_e, double Rmax) const
{
    return map2_e * (r/Rmax) * sin(theta) / (sqrt(1.0 - 1.0 / 4.0 * (map2_epsilon * map2_epsilon)) * (2.0 - sqrt(map2_epsilon * (map2_epsilon + 2.0 * (r/Rmax) * cos(theta)) + 1.0)));
}

double CartesianR6PoissonTriangular::J_rr(double r, double theta, double map2_epsilon, double map2_e, double Rmax) const
{
    return ((-cos(theta)) / sqrt(map2_epsilon * (map2_epsilon + 2.0 * (r/Rmax) * cos(theta)) + 1.0))/Rmax;
}

double CartesianR6PoissonTriangular::J_rt(double r, double theta, double map2_epsilon, double map2_e, double Rmax) const
{
    return (r/Rmax) * sin(theta) / sqrt(map2_epsilon * (map2_epsilon + 2.0 * (r/Rmax) * cos(theta)) + 1.0);
}

double CartesianR6PoissonTriangular::J_tr(double r, double theta, double map2_epsilon, double map2_e, double Rmax) const
{
    return (map2_e * map2_epsilon * (r/Rmax) * sin(theta) * cos(theta) / (sqrt(1.0 - 1.0 / 4.0 * (map2_epsilon * map2_epsilon)) * pow((2.0 - sqrt(map2_epsilon * (map2_epsilon + 2.0 * (r/Rmax) * cos(theta)) + 1.0)), 2.0) * sqrt(map2_epsilon * (map2_epsilon + 2.0 * (r/Rmax) * cos(theta)) + 1.0)) + map2_e * sin(theta) / (sqrt(1.0 - 1.0 / 4.0 * (map2_epsilon * map2_epsilon)) * (2.0 - sqrt(map2_epsilon * (map2_epsilon + 2.0 * (r/Rmax) * cos(theta)) + 1.0))))/Rmax;
}

double CartesianR6PoissonTriangular::J_tt(double r, double theta, double map2_epsilon, double map2_e, double Rmax) const
{
    return (r/Rmax) * ((-map2_e) * map2_epsilon * (r/Rmax) * pow(sin(theta), 2.0) / (sqrt(1.0 - 1.0 / 4.0 * (map2_epsilon * map2_epsilon)) * pow((2.0 - sqrt(map2_epsilon * (map2_epsilon + 2.0 * (r/Rmax) * cos(theta)) + 1.0)), 2.0) * sqrt(map2_epsilon * (map2_epsilon + 2.0 * (r/Rmax) * cos(theta)) + 1.0)) + map2_e * cos(theta) / (sqrt(1.0 - 1.0 / 4.0 * (map2_epsilon * map2_epsilon)) * (2.0 - sqrt(map2_epsilon * (map2_epsilon + 2.0 * (r/Rmax) * cos(theta)) + 1.0))));
}

double CartesianR6PoissonTriangular::J_xs(double r, double theta, double map2_epsilon, double map2_e, double Rmax) const
{
    return (pow(map2_epsilon, 6.0) * pow(sin(theta), 2.0) - pow(map2_epsilon, 6.0) + 5.0 * pow(map2_epsilon, 5.0) * (r/Rmax) * pow(sin(theta), 2.0) * cos(theta) - 6.0 * pow(map2_epsilon, 5.0) * (r/Rmax) * cos(theta) + 8.0 * pow(map2_epsilon, 4.0) * ((r/Rmax) * (r/Rmax)) * pow(sin(theta), 2.0) * pow(cos(theta), 2.0) - 12.0 * pow(map2_epsilon, 4.0) * ((r/Rmax) * (r/Rmax)) * pow(cos(theta), 2.0) - 8.0 * pow(map2_epsilon, 4.0) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos(theta) + 1.0) * pow(sin(theta), 2.0) + 8.0 * pow(map2_epsilon, 4.0) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos(theta) + 1.0) + 27.0 * pow(map2_epsilon, 4.0) * pow(sin(theta), 2.0) - 27.0 * pow(map2_epsilon, 4.0) + 4.0 * pow(map2_epsilon, 3.0) * pow((r/Rmax), 3.0) * pow(sin(theta), 2.0) * pow(cos(theta), 3.0) - 8.0 * pow(map2_epsilon, 3.0) * pow((r/Rmax), 3.0) * pow(cos(theta), 3.0) - 26.0 * pow(map2_epsilon, 3.0) * (r/Rmax) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos(theta) + 1.0) * pow(sin(theta), 2.0) * cos(theta) + 32.0 * pow(map2_epsilon, 3.0) * (r/Rmax) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos(theta) + 1.0) * cos(theta) + 94.0 * pow(map2_epsilon, 3.0) * (r/Rmax) * pow(sin(theta), 2.0) * cos(theta) - 108.0 * pow(map2_epsilon, 3.0) * (r/Rmax) * cos(theta) - 20.0 * (map2_epsilon * map2_epsilon) * ((r/Rmax) * (r/Rmax)) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos(theta) + 1.0) * pow(sin(theta), 2.0) * pow(cos(theta), 2.0) + 32.0 * (map2_epsilon * map2_epsilon) * ((r/Rmax) * (r/Rmax)) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos(theta) + 1.0) * pow(cos(theta), 2.0) + 80.0 * (map2_epsilon * map2_epsilon) * ((r/Rmax) * (r/Rmax)) * pow(sin(theta), 2.0) * pow(cos(theta), 2.0) - 108.0 * (map2_epsilon * map2_epsilon) * ((r/Rmax) * (r/Rmax)) * pow(cos(theta), 2.0) - 48.0 * (map2_epsilon * map2_epsilon) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos(theta) + 1.0) * pow(sin(theta), 2.0) + 48.0 * (map2_epsilon * map2_epsilon) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos(theta) + 1.0) + 67.0 * (map2_epsilon * map2_epsilon) * pow(sin(theta), 2.0) - 67.0 * (map2_epsilon * map2_epsilon) - 82.0 * map2_epsilon * (r/Rmax) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos(theta) + 1.0) * pow(sin(theta), 2.0) * cos(theta) + 96.0 * map2_epsilon * (r/Rmax) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos(theta) + 1.0) * cos(theta) + 121.0 * map2_epsilon * (r/Rmax) * pow(sin(theta), 2.0) * cos(theta) - 134.0 * map2_epsilon * (r/Rmax) * cos(theta) - 40.0 * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos(theta) + 1.0) * pow(sin(theta), 2.0) + 40.0 * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos(theta) + 1.0) + 41.0 * pow(sin(theta), 2.0) - 41.0) / (pow(map2_epsilon, 4.0) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos(theta) + 1.0) * cos(theta) - 8.0 * pow(map2_epsilon, 4.0) * cos(theta) + 4.0 * pow(map2_epsilon, 3.0) * (r/Rmax) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos(theta) + 1.0) * pow(cos(theta), 2.0) - 32.0 * pow(map2_epsilon, 3.0) * (r/Rmax) * pow(cos(theta), 2.0) + 4.0 * (map2_epsilon * map2_epsilon) * ((r/Rmax) * (r/Rmax)) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos(theta) + 1.0) * pow(cos(theta), 3.0) - 32.0 * (map2_epsilon * map2_epsilon) * ((r/Rmax) * (r/Rmax)) * pow(cos(theta), 3.0) + 26.0 * (map2_epsilon * map2_epsilon) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos(theta) + 1.0) * cos(theta) - 48.0 * (map2_epsilon * map2_epsilon) * cos(theta) + 52.0 * map2_epsilon * (r/Rmax) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos(theta) + 1.0) * pow(cos(theta), 2.0) - 96.0 * map2_epsilon * (r/Rmax) * pow(cos(theta), 2.0) + 41.0 * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos(theta) + 1.0) * cos(theta) - 40.0 * cos(theta));
}

double CartesianR6PoissonTriangular::J_xt(double r, double theta, double map2_epsilon, double map2_e, double Rmax) const
{
    return 1.0 / 2.0 * ((-(map2_epsilon * map2_epsilon)) * sqrt(4.0 - map2_epsilon * map2_epsilon) * sin(theta) - 2.0 * map2_epsilon * (r/Rmax) * sqrt(4.0 - map2_epsilon * map2_epsilon) * sin(theta) * cos(theta) + 2.0 * sqrt(4.0 - map2_epsilon * map2_epsilon) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos(theta) + 1.0) * sin(theta) - sqrt(4.0 - map2_epsilon * map2_epsilon) * sin(theta)) / (map2_e * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos(theta) + 1.0));
}

double CartesianR6PoissonTriangular::J_ys(double r, double theta, double map2_epsilon, double map2_e, double Rmax) const
{
    return (pow(map2_epsilon, 4.0) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos(theta) + 1.0) * sin(theta) - 6.0 * pow(map2_epsilon, 4.0) * sin(theta) + 3.0 * pow(map2_epsilon, 3.0) * (r/Rmax) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos(theta) + 1.0) * sin(theta) * cos(theta) - 20.0 * pow(map2_epsilon, 3.0) * (r/Rmax) * sin(theta) * cos(theta) + 2.0 * (map2_epsilon * map2_epsilon) * ((r/Rmax) * (r/Rmax)) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos(theta) + 1.0) * sin(theta) * pow(cos(theta), 2.0) - 16.0 * (map2_epsilon * map2_epsilon) * ((r/Rmax) * (r/Rmax)) * sin(theta) * pow(cos(theta), 2.0) + 14.0 * (map2_epsilon * map2_epsilon) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos(theta) + 1.0) * sin(theta) - 20.0 * (map2_epsilon * map2_epsilon) * sin(theta) + 23.0 * map2_epsilon * (r/Rmax) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos(theta) + 1.0) * sin(theta) * cos(theta) - 36.0 * map2_epsilon * (r/Rmax) * sin(theta) * cos(theta) + 13.0 * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos(theta) + 1.0) * sin(theta) - 14.0 * sin(theta)) / (pow(map2_epsilon, 4.0) * (r/Rmax) + 4.0 * pow(map2_epsilon, 3.0) * ((r/Rmax) * (r/Rmax)) * cos(theta) + 4.0 * (map2_epsilon * map2_epsilon) * pow((r/Rmax), 3.0) * pow(cos(theta), 2.0) - 6.0 * (map2_epsilon * map2_epsilon) * (r/Rmax) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos(theta) + 1.0) + 14.0 * (map2_epsilon * map2_epsilon) * (r/Rmax) - 12.0 * map2_epsilon * ((r/Rmax) * (r/Rmax)) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos(theta) + 1.0) * cos(theta) + 28.0 * map2_epsilon * ((r/Rmax) * (r/Rmax)) * cos(theta) - 14.0 * (r/Rmax) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos(theta) + 1.0) + 13.0 * (r/Rmax));
}

double CartesianR6PoissonTriangular::J_yt(double r, double theta, double map2_epsilon, double map2_e, double Rmax) const
{
    return 1.0 / 2.0 * ((-sqrt(4.0 - map2_epsilon * map2_epsilon)) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos(theta) + 1.0) * cos(theta) + 2.0 * sqrt(4.0 - map2_epsilon * map2_epsilon) * cos(theta)) / (map2_e * (r/Rmax));
}

double CartesianR6PoissonTriangular::rho_glob(double r, double theta, double map2_epsilon, double map2_e, double Rmax) const
{
    return 0.0;
}

double CartesianR6PoissonTriangular::rho_pole(double r, double theta, double map2_epsilon, double map2_e, double Rmax) const
{
    return 0.0;
}

double CartesianR6PoissonTriangular::coeffs1(double r, double Rmax) const
{
    return 0.0;
}

double CartesianR6PoissonTriangular::coeffs2(double r, double Rmax) const
{
    return 0.0;
}

double CartesianR6PoissonTriangular::phi_exact(double r, double theta, double map2_epsilon, double map2_e, double Rmax) const
{
    return 0.4096 * pow(((r/Rmax) - 1.0), 6.0) * pow(((r/Rmax) + 1.0), 6.0) * sin(2.0 * M_PI * map2_e * (r/Rmax) * sin(theta) / (sqrt(1.0 - 1.0 / 4.0 * (map2_epsilon * map2_epsilon)) * (2.0 - sqrt(map2_epsilon * (map2_epsilon + 2.0 * (r/Rmax) * cos(theta)) + 1.0)))) * cos(2.0 * M_PI * (1.0 - sqrt(map2_epsilon * (map2_epsilon + 2.0 * (r/Rmax) * cos(theta)) + 1.0)) / map2_epsilon);
}



double CartesianR6PoissonTriangular::x(double r, double theta, double map2_epsilon, double map2_e, double Rmax, double sin_theta, double cos_theta) const
{
    return (1.0 - sqrt(map2_epsilon * (map2_epsilon + 2.0 * (r/Rmax) * cos_theta) + 1.0)) / map2_epsilon;
}

double CartesianR6PoissonTriangular::y(double r, double theta, double map2_epsilon, double map2_e, double Rmax, double sin_theta, double cos_theta) const
{
    return map2_e * (r/Rmax) * sin_theta / (sqrt(1.0 - 1.0 / 4.0 * (map2_epsilon * map2_epsilon)) * (2.0 - sqrt(map2_epsilon * (map2_epsilon + 2.0 * (r/Rmax) * cos_theta) + 1.0)));
}

double CartesianR6PoissonTriangular::J_rr(double r, double theta, double map2_epsilon, double map2_e, double Rmax, double sin_theta, double cos_theta) const
{
    return ((-cos_theta) / sqrt(map2_epsilon * (map2_epsilon + 2.0 * (r/Rmax) * cos_theta) + 1.0))/Rmax;
}

double CartesianR6PoissonTriangular::J_rt(double r, double theta, double map2_epsilon, double map2_e, double Rmax, double sin_theta, double cos_theta) const
{
    return (r/Rmax) * sin_theta / sqrt(map2_epsilon * (map2_epsilon + 2.0 * (r/Rmax) * cos_theta) + 1.0);
}

double CartesianR6PoissonTriangular::J_tr(double r, double theta, double map2_epsilon, double map2_e, double Rmax, double sin_theta, double cos_theta) const
{
    return (map2_e * map2_epsilon * (r/Rmax) * sin_theta * cos_theta / (sqrt(1.0 - 1.0 / 4.0 * (map2_epsilon * map2_epsilon)) * pow((2.0 - sqrt(map2_epsilon * (map2_epsilon + 2.0 * (r/Rmax) * cos_theta) + 1.0)), 2.0) * sqrt(map2_epsilon * (map2_epsilon + 2.0 * (r/Rmax) * cos_theta) + 1.0)) + map2_e * sin_theta / (sqrt(1.0 - 1.0 / 4.0 * (map2_epsilon * map2_epsilon)) * (2.0 - sqrt(map2_epsilon * (map2_epsilon + 2.0 * (r/Rmax) * cos_theta) + 1.0))))/Rmax;
}

double CartesianR6PoissonTriangular::J_tt(double r, double theta, double map2_epsilon, double map2_e, double Rmax, double sin_theta, double cos_theta) const
{
    return (r/Rmax) * ((-map2_e) * map2_epsilon * (r/Rmax) * pow(sin_theta, 2.0) / (sqrt(1.0 - 1.0 / 4.0 * (map2_epsilon * map2_epsilon)) * pow((2.0 - sqrt(map2_epsilon * (map2_epsilon + 2.0 * (r/Rmax) * cos_theta) + 1.0)), 2.0) * sqrt(map2_epsilon * (map2_epsilon + 2.0 * (r/Rmax) * cos_theta) + 1.0)) + map2_e * cos_theta / (sqrt(1.0 - 1.0 / 4.0 * (map2_epsilon * map2_epsilon)) * (2.0 - sqrt(map2_epsilon * (map2_epsilon + 2.0 * (r/Rmax) * cos_theta) + 1.0))));
}

double CartesianR6PoissonTriangular::J_xs(double r, double theta, double map2_epsilon, double map2_e, double Rmax, double sin_theta, double cos_theta) const
{
    return (pow(map2_epsilon, 6.0) * pow(sin_theta, 2.0) - pow(map2_epsilon, 6.0) + 5.0 * pow(map2_epsilon, 5.0) * (r/Rmax) * pow(sin_theta, 2.0) * cos_theta - 6.0 * pow(map2_epsilon, 5.0) * (r/Rmax) * cos_theta + 8.0 * pow(map2_epsilon, 4.0) * ((r/Rmax) * (r/Rmax)) * pow(sin_theta, 2.0) * pow(cos_theta, 2.0) - 12.0 * pow(map2_epsilon, 4.0) * ((r/Rmax) * (r/Rmax)) * pow(cos_theta, 2.0) - 8.0 * pow(map2_epsilon, 4.0) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos_theta + 1.0) * pow(sin_theta, 2.0) + 8.0 * pow(map2_epsilon, 4.0) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos_theta + 1.0) + 27.0 * pow(map2_epsilon, 4.0) * pow(sin_theta, 2.0) - 27.0 * pow(map2_epsilon, 4.0) + 4.0 * pow(map2_epsilon, 3.0) * pow((r/Rmax), 3.0) * pow(sin_theta, 2.0) * pow(cos_theta, 3.0) - 8.0 * pow(map2_epsilon, 3.0) * pow((r/Rmax), 3.0) * pow(cos_theta, 3.0) - 26.0 * pow(map2_epsilon, 3.0) * (r/Rmax) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos_theta + 1.0) * pow(sin_theta, 2.0) * cos_theta + 32.0 * pow(map2_epsilon, 3.0) * (r/Rmax) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos_theta + 1.0) * cos_theta + 94.0 * pow(map2_epsilon, 3.0) * (r/Rmax) * pow(sin_theta, 2.0) * cos_theta - 108.0 * pow(map2_epsilon, 3.0) * (r/Rmax) * cos_theta - 20.0 * (map2_epsilon * map2_epsilon) * ((r/Rmax) * (r/Rmax)) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos_theta + 1.0) * pow(sin_theta, 2.0) * pow(cos_theta, 2.0) + 32.0 * (map2_epsilon * map2_epsilon) * ((r/Rmax) * (r/Rmax)) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos_theta + 1.0) * pow(cos_theta, 2.0) + 80.0 * (map2_epsilon * map2_epsilon) * ((r/Rmax) * (r/Rmax)) * pow(sin_theta, 2.0) * pow(cos_theta, 2.0) - 108.0 * (map2_epsilon * map2_epsilon) * ((r/Rmax) * (r/Rmax)) * pow(cos_theta, 2.0) - 48.0 * (map2_epsilon * map2_epsilon) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos_theta + 1.0) * pow(sin_theta, 2.0) + 48.0 * (map2_epsilon * map2_epsilon) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos_theta + 1.0) + 67.0 * (map2_epsilon * map2_epsilon) * pow(sin_theta, 2.0) - 67.0 * (map2_epsilon * map2_epsilon) - 82.0 * map2_epsilon * (r/Rmax) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos_theta + 1.0) * pow(sin_theta, 2.0) * cos_theta + 96.0 * map2_epsilon * (r/Rmax) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos_theta + 1.0) * cos_theta + 121.0 * map2_epsilon * (r/Rmax) * pow(sin_theta, 2.0) * cos_theta - 134.0 * map2_epsilon * (r/Rmax) * cos_theta - 40.0 * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos_theta + 1.0) * pow(sin_theta, 2.0) + 40.0 * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos_theta + 1.0) + 41.0 * pow(sin_theta, 2.0) - 41.0) / (pow(map2_epsilon, 4.0) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos_theta + 1.0) * cos_theta - 8.0 * pow(map2_epsilon, 4.0) * cos_theta + 4.0 * pow(map2_epsilon, 3.0) * (r/Rmax) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos_theta + 1.0) * pow(cos_theta, 2.0) - 32.0 * pow(map2_epsilon, 3.0) * (r/Rmax) * pow(cos_theta, 2.0) + 4.0 * (map2_epsilon * map2_epsilon) * ((r/Rmax) * (r/Rmax)) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos_theta + 1.0) * pow(cos_theta, 3.0) - 32.0 * (map2_epsilon * map2_epsilon) * ((r/Rmax) * (r/Rmax)) * pow(cos_theta, 3.0) + 26.0 * (map2_epsilon * map2_epsilon) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos_theta + 1.0) * cos_theta - 48.0 * (map2_epsilon * map2_epsilon) * cos_theta + 52.0 * map2_epsilon * (r/Rmax) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos_theta + 1.0) * pow(cos_theta, 2.0) - 96.0 * map2_epsilon * (r/Rmax) * pow(cos_theta, 2.0) + 41.0 * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos_theta + 1.0) * cos_theta - 40.0 * cos_theta);
}

double CartesianR6PoissonTriangular::J_xt(double r, double theta, double map2_epsilon, double map2_e, double Rmax, double sin_theta, double cos_theta) const
{
    return 1.0 / 2.0 * ((-(map2_epsilon * map2_epsilon)) * sqrt(4.0 - map2_epsilon * map2_epsilon) * sin_theta - 2.0 * map2_epsilon * (r/Rmax) * sqrt(4.0 - map2_epsilon * map2_epsilon) * sin_theta * cos_theta + 2.0 * sqrt(4.0 - map2_epsilon * map2_epsilon) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos_theta + 1.0) * sin_theta - sqrt(4.0 - map2_epsilon * map2_epsilon) * sin_theta) / (map2_e * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos_theta + 1.0));
}

double CartesianR6PoissonTriangular::J_ys(double r, double theta, double map2_epsilon, double map2_e, double Rmax, double sin_theta, double cos_theta) const
{
    return (pow(map2_epsilon, 4.0) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos_theta + 1.0) * sin_theta - 6.0 * pow(map2_epsilon, 4.0) * sin_theta + 3.0 * pow(map2_epsilon, 3.0) * (r/Rmax) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos_theta + 1.0) * sin_theta * cos_theta - 20.0 * pow(map2_epsilon, 3.0) * (r/Rmax) * sin_theta * cos_theta + 2.0 * (map2_epsilon * map2_epsilon) * ((r/Rmax) * (r/Rmax)) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos_theta + 1.0) * sin_theta * pow(cos_theta, 2.0) - 16.0 * (map2_epsilon * map2_epsilon) * ((r/Rmax) * (r/Rmax)) * sin_theta * pow(cos_theta, 2.0) + 14.0 * (map2_epsilon * map2_epsilon) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos_theta + 1.0) * sin_theta - 20.0 * (map2_epsilon * map2_epsilon) * sin_theta + 23.0 * map2_epsilon * (r/Rmax) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos_theta + 1.0) * sin_theta * cos_theta - 36.0 * map2_epsilon * (r/Rmax) * sin_theta * cos_theta + 13.0 * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos_theta + 1.0) * sin_theta - 14.0 * sin_theta) / (pow(map2_epsilon, 4.0) * (r/Rmax) + 4.0 * pow(map2_epsilon, 3.0) * ((r/Rmax) * (r/Rmax)) * cos_theta + 4.0 * (map2_epsilon * map2_epsilon) * pow((r/Rmax), 3.0) * pow(cos_theta, 2.0) - 6.0 * (map2_epsilon * map2_epsilon) * (r/Rmax) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos_theta + 1.0) + 14.0 * (map2_epsilon * map2_epsilon) * (r/Rmax) - 12.0 * map2_epsilon * ((r/Rmax) * (r/Rmax)) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos_theta + 1.0) * cos_theta + 28.0 * map2_epsilon * ((r/Rmax) * (r/Rmax)) * cos_theta - 14.0 * (r/Rmax) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos_theta + 1.0) + 13.0 * (r/Rmax));
}

double CartesianR6PoissonTriangular::J_yt(double r, double theta, double map2_epsilon, double map2_e, double Rmax, double sin_theta, double cos_theta) const
{
    return 1.0 / 2.0 * ((-sqrt(4.0 - map2_epsilon * map2_epsilon)) * sqrt(map2_epsilon * map2_epsilon + 2.0 * map2_epsilon * (r/Rmax) * cos_theta + 1.0) * cos_theta + 2.0 * sqrt(4.0 - map2_epsilon * map2_epsilon) * cos_theta) / (map2_e * (r/Rmax));
}

double CartesianR6PoissonTriangular::rho_glob(double r, double theta, double map2_epsilon, double map2_e, double Rmax, double sin_theta, double cos_theta) const
{
    return 0.0;
}

double CartesianR6PoissonTriangular::rho_pole(double r, double theta, double map2_epsilon, double map2_e, double Rmax, double sin_theta, double cos_theta) const
{
    return 0.0;
}

double CartesianR6PoissonTriangular::phi_exact(double r, double theta, double map2_epsilon, double map2_e, double Rmax, double sin_theta, double cos_theta) const
{
    return 0.4096 * pow(((r/Rmax) - 1.0), 6.0) * pow(((r/Rmax) + 1.0), 6.0) * sin(2.0 * M_PI * map2_e * (r/Rmax) * sin_theta / (sqrt(1.0 - 1.0 / 4.0 * (map2_epsilon * map2_epsilon)) * (2.0 - sqrt(map2_epsilon * (map2_epsilon + 2.0 * (r/Rmax) * cos_theta) + 1.0)))) * cos(2.0 * M_PI * (1.0 - sqrt(map2_epsilon * (map2_epsilon + 2.0 * (r/Rmax) * cos_theta) + 1.0)) / map2_epsilon);
}
