#include "../../include/InputFunctions/exact_solution.h"

// In earlier versions denoted by 'exact_phi'
inline double exact_solution_Functor::operator()(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
    static constexpr double Rmax = 1.3;
    static constexpr double map2_epsilon = 0.3;
    static constexpr double map2_e = 1.4;
    return (1.0 - (r/Rmax) * (r/Rmax)) * sin(2.0 * M_PI * map2_e * (r/Rmax) * sin_theta / (sqrt(1.0 - 1.0 / 4.0 * (map2_epsilon * map2_epsilon)) * (2.0 - sqrt(map2_epsilon * (map2_epsilon + 2.0 * (r/Rmax) * cos_theta) + 1.0)))) * cos(2.0 * M_PI * (1.0 - sqrt(map2_epsilon * (map2_epsilon + 2.0 * (r/Rmax) * cos_theta) + 1.0)) / map2_epsilon);
}
