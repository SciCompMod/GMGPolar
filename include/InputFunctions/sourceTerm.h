#pragma once

class SourceTerm
{
public:
    SourceTerm()          = default;
    virtual ~SourceTerm() = default;

    virtual double rhs_f(double r, double theta) const = 0;
};
