#include "../../include/Interpolation/interpolation.h"

Interpolation::Interpolation(const int maxOpenMPThreads, const std::vector<int>& taskingThreads) :
    maxOpenMPThreads_(maxOpenMPThreads),
    taskingThreads_(taskingThreads)
{}