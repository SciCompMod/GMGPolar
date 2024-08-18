#pragma once

class SourceTerm {
public:
    SourceTerm() = default;
    virtual ~SourceTerm() = default;

    virtual double rhs_f(const double& r, const double& theta, const double& sin_theta, const double& cos_theta) const = 0;
};