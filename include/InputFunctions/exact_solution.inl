#include "../../include/InputFunctions/exact_solution.h"

// In earlier versions denoted by 'exact_phi'
inline double exact_solution_Functor::operator()(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
    static constexpr double Rmax = 1.3;
    return (1.0 - (r/Rmax) * (r/Rmax)) * sin(2.0 * M_PI * (r/Rmax) * sin_theta) * cos(2.0 * M_PI * (r/Rmax) * cos_theta);
}
