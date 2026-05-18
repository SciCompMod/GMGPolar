#pragma once

namespace concepts
{

template <typename T>
concept ExactSolution = requires(const T solution, double r, double theta) {
    { solution.exact_solution(r, theta) } -> std::convertible_to<double>;
};

} // namespace concepts