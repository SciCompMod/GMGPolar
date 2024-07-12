#pragma once

#include <memory>

#include "densityProfileCoefficients.h"
#include "sourceTerm.h"
#include "boundaryConditions.h"

class SystemParameters {
public:
    SystemParameters(
        std::unique_ptr<DensityProfileCoefficients> coeffs, 
        std::unique_ptr<BoundaryConditions> u_D, 
        std::unique_ptr<SourceTerm> rhs_f
    ) : 
        densityProfileCoeffs(std::move(coeffs)),
        boundaryConditions(std::move(u_D)),
        sourceTerm(std::move(rhs_f)) 
    {}

    // DensityProfileCoeffs functions
    double alpha(const double& r) const {
        return densityProfileCoeffs->alpha(r);
    }

    double beta(const double& r) const {
        return densityProfileCoeffs->beta(r);
    }

    double getAlphaJump() const {
        return densityProfileCoeffs->getAlphaJump();
    }

    // BoundaryConditions functions
    double u_D(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
        return boundaryConditions->u_D(r, theta, sin_theta, cos_theta);
    }

    double u_D_Interior(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
        return boundaryConditions->u_D_Interior(r, theta, sin_theta, cos_theta);
    }

    // SourceTerm functions
    double rhs_f(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const {
        return sourceTerm->rhs_f(r, theta, sin_theta, cos_theta);
    }

private:
    std::unique_ptr<DensityProfileCoefficients> densityProfileCoeffs;
    std::unique_ptr<BoundaryConditions> boundaryConditions;
    std::unique_ptr<SourceTerm> sourceTerm;
};