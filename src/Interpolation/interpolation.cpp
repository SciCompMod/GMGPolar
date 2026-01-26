#include "../../include/Interpolation/interpolation.h"

Interpolation::Interpolation(int max_omp_threads, bool DirBC_Interior)
    : max_omp_threads_(max_omp_threads)
    , DirBC_Interior_(DirBC_Interior)
{
}
