#pragma once

class SourceTerm
{
public:
    SourceTerm()          = default;
    virtual ~SourceTerm() = default;

    virtual double rhs_f(const double& r, const double& theta) const = 0;
};
