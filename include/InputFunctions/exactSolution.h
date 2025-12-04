#pragma once

class ExactSolution
{
public:
    ExactSolution()          = default;
    virtual ~ExactSolution() = default;

    virtual double exact_solution(const double& r, const double& theta) const = 0;
};
