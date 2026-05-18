#pragma once

#include <cmath>

#include "../exactSolution.h"

namespace gmgpolar
{

class Refined_ShafranovGeometry
{
public:
    Refined_ShafranovGeometry() = default;
    explicit Refined_ShafranovGeometry(double Rmax, double elongation_kappa, double shift_delta);
    virtual ~Refined_ShafranovGeometry() = default;

    double exact_solution(double r, double theta) const;

private:
    const double Rmax             = 1.3;
    const double elongation_kappa = 0.3;
    const double shift_delta      = 0.2;
};
} // namespace gmgpolar
