#pragma once

#include <array>
#include <initializer_list>

#include "../PolarGrid/polargrid.h"

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


struct Stencil {
    std::array<int, 9> values;

    Stencil(std::initializer_list<int> init);

    int operator[](StencilType type) const;
};

namespace matrixA_Stencil {
    const Stencil& get_stencil(const PolarGrid& grid, int i_r, int DirBC_Interior);

    int nnz_matrixA(const PolarGrid& grid, const int DirBC_Interior);
    int ptr_nz_index_matrixA(const PolarGrid& grid, const int i_r, const int i_theta, const int DirBC_Interior);
}

namespace symmetric_matrixA_Stencil {
    const Stencil& get_stencil(const PolarGrid& grid, int i_r, int DirBC_Interior);
    int nnz_matrixA(const PolarGrid& grid, const int DirBC_Interior);
    int ptr_nz_index_matrixA(const PolarGrid& grid, const int i_r, const int i_theta, const int DirBC_Interior);
}