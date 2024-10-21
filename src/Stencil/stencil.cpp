#include "../../include/Stencil/stencil.h"

Stencil::Stencil(std::initializer_list<int> init) : values_{} {
    std::copy(init.begin(), init.end(), values_.begin());
    stencil_size_ = 0;
    for (std::size_t i = 0; i < init.size(); i++){
        if(values_[i] != -1) stencil_size_ ++;
    }
}

int Stencil::operator[](StencilType type) const {
    return values_[static_cast<size_t>(type)];
}