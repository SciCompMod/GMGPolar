#pragma once

namespace gmgpolar
{

class ExactSolution
{
public:
    ExactSolution()          = default;
    virtual ~ExactSolution() = default;

    virtual double exact_solution(double r, double theta) const = 0;
};
} // namespace gmgpolar
