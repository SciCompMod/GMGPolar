#pragma once

class SourceTerm
{
public:
    SourceTerm()          = default;
    virtual ~SourceTerm() = default;

    virtual double operator()(std::size_t i_r, std::size_t i_theta) const = 0;
};
