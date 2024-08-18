#include "../include/InputFunctions/ExactSolution/refined_ShafranovGeometry.h"

Refined_ShafranovGeometry::Refined_ShafranovGeometry(const double& Rmax, const double& elongation_kappa, const double& shift_delta) : 
    Rmax(Rmax),
    elongation_kappa(elongation_kappa),
    shift_delta(shift_delta)
{}

double Refined_ShafranovGeometry::exact_solution(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
    return ((-3.33823779536505e-15) * ((r/Rmax) * (r/Rmax)) - 0.0 * (r/Rmax) - 0.0 + exp((-3333.33333333333) * pow(((r/Rmax) - 0.9), 2.0))) * cos(21.0 * theta) + 
        (0.00184273372222541 * ((r/Rmax) * (r/Rmax)) - 0.0018029383826828 * (r/Rmax) - 4.00652973929511e-05 + exp((-50.0) * pow(((r/Rmax) - 0.45), 2.0))) * cos(9.0 * theta);
}
