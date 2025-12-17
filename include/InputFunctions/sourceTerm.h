#pragma once

class SourceTerm
{
public:
    SourceTerm()          = default;
    virtual ~SourceTerm() = default;

    virtual double operator()(int i_r, int i_theta) const = 0;
};
