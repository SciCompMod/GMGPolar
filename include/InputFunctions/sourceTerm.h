#pragma once

class SourceTerm
{
public:
    SourceTerm()          = default;
    virtual ~SourceTerm() = default;

    virtual double operator()(double r, double theta) const = 0;
};
