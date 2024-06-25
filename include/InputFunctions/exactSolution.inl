#include "../../include/InputFunctions/exactSolution.h"

// In earlier versions denoted by 'exact_phi'
inline double ExactSolution::exact_solution(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
    return 0.4096 * pow(((r/Rmax) - 1.0), 6.0) * pow(((r/Rmax) + 1.0), 6.0) * sin(M_PI * (2.0 * map1_kappa * (r/Rmax) * sin_theta + 2.0 * (r/Rmax) * sin_theta)) * cos(M_PI * ((-2.0) * map1_delta * ((r/Rmax) * (r/Rmax)) - 2.0 * map1_kappa * (r/Rmax) * cos_theta + 2.0 * (r/Rmax) * cos_theta));
}
