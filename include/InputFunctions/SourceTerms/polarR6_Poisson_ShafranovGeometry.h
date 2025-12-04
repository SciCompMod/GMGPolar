#pragma once

#include <cmath>

#include "../sourceTerm.h"

class PolarR6_Poisson_ShafranovGeometry : public SourceTerm
{
public:
    PolarR6_Poisson_ShafranovGeometry() = default;
    explicit PolarR6_Poisson_ShafranovGeometry(const double& Rmax, const double& elongation_kappa,
                                               const double& shift_delta);
    virtual ~PolarR6_Poisson_ShafranovGeometry() = default;

    double rhs_f(const double& r, const double& theta) const override;

private:
    const double Rmax             = 1.3;
    const double elongation_kappa = 0.3;
    const double shift_delta      = 0.2;
};
