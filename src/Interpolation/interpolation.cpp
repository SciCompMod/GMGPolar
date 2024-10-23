#include "../../include/Interpolation/interpolation.h"

Interpolation::Interpolation(const std::vector<int>& threads_per_level, const bool DirBC_Interior)
    : threads_per_level_(threads_per_level)
    , DirBC_Interior_(DirBC_Interior)
{
}