#include "../../include/Interpolation/interpolation.h"

#include "../../include/LinearAlgebra/Vector/vector_operations.h"

Interpolation::Interpolation(const bool DirBC_Interior)
    : DirBC_Interior_(DirBC_Interior)
{}