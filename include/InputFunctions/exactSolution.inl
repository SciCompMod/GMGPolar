#include "../../include/InputFunctions/exactSolution.h"

// In earlier versions denoted by 'exact_phi'
inline double ExactSolution::exact_solution(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
    return (1.0 - (r/Rmax) * (r/Rmax)) * sin(2.0 * M_PI * map2_e * (r/Rmax) * sin_theta / (sqrt(1.0 - 1.0 / 4.0 * (map2_epsilon * map2_epsilon)) * (2.0 - sqrt(map2_epsilon * (map2_epsilon + 2.0 * (r/Rmax) * cos_theta) + 1.0)))) * cos(2.0 * M_PI * (1.0 - sqrt(map2_epsilon * (map2_epsilon + 2.0 * (r/Rmax) * cos_theta) + 1.0)) / map2_epsilon);
}
