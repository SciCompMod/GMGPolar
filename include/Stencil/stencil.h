#pragma once

#include <array>
#include <initializer_list>

enum class StencilPosition
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
    int operator[](StencilPosition type) const;

private:
    std::array<int, 9> values_;
    int stencil_size_;
};