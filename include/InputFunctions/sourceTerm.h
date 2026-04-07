#pragma once

namespace concepts
{

template <typename T>
concept SourceTerm = requires(const T source, std::size_t i_r, std::size_t i_theta) {
    { source(i_r, i_theta) } -> std::convertible_to<double>;
};

} // namespace concepts
