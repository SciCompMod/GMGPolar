#include "../../include/Interpolation/interpolation.h"

Interpolation::Interpolation(const std::vector<int>& threads_per_level) :
    threads_per_level_(threads_per_level)
{}