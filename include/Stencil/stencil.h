#pragma once

#include <array>
#include <initializer_list>

enum class StencilType
{
    TopLeft,
    Top,
    TopRight,
    Left,
    Center,
    Right,
    BottomLeft,
    Bottom,
    BottomRight,
};

class Stencil
{
public:
    Stencil(std::initializer_list<int> init);
    int operator[](StencilType type) const;

private:
    std::array<int, 9> values_;
    int stencil_size_;
};