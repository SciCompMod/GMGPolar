#pragma once

class ExactSolution
{
public:
    ExactSolution()          = default;
    virtual ~ExactSolution() = default;

    virtual double exact_solution(const double& r, const double& theta, const double& sin_theta,
                                  const double& cos_theta) const = 0;
};