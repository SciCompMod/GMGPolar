#pragma once

#include <initializer_list>

#include <Kokkos_Core.hpp>

#include "../include/LinearAlgebra/Vector/vector.h"

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

/**
 * @brief Represents a stencil pattern used in sparse matrix construction 
 * (see p. 44 table 4 of Julian Litz Master Thesis)
 * 
 * The Stencil class helps define neighborhood structure of a grid node.
 * It maps each `StencilPosition` to an integer index that represents its 
 * inclusion in the stencil pattern.
 *
 * Non-zero indices are obtained via `ptr + offset`, where:
 * - The Stencil class stores the offset for each position.
 * - A value of `-1` means the position is not included in the stencil pattern.
 * - Other values (0, 1, 2, ..., stencil_size - 1) correspond to valid stencil indices.
 */
class Stencil
{
public:
    Stencil(std::initializer_list<int> init)
        : values_("stencil_values", 9)
        , stencil_size_(0)
    {
        int i = 0;
        for (int v : init) {
            values_(i++) = v;
            if (v != -1)
                stencil_size_++;
        }
    }

    KOKKOS_INLINE_FUNCTION int operator[](StencilPosition type) const
    {
        return values_(static_cast<int>(type));
    }

    KOKKOS_INLINE_FUNCTION int stencil_size() const
    {
        return stencil_size_;
    }

private:
    AllocatableVector<int> values_;
    int stencil_size_;
};
