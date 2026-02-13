#pragma once

#include <vector>
#include <cassert>
#include <iostream>
#include <omp.h>

#include "../PolarGrid/polargrid.h"
#include "../LinearAlgebra/Vector/vector.h"
#include "../LinearAlgebra/Vector/vector_operations.h"
#include "../Definitions/global_definitions.h"

class Interpolation
{
public:
    explicit Interpolation(int max_omp_threads, bool DirBC_Interior);

    /* Remark: This injection is not scaled. */
    void applyInjection(const PolarGrid& fine_grid, const PolarGrid& coarse_grid, Vector<double> coarse_result,
                        ConstVector<double> fine_values) const;

    /* Bilinear interpolation operator */
    void applyProlongation(const PolarGrid& coarse_grid, const PolarGrid& fine_grid, Vector<double> fine_result,
                           ConstVector<double> coarse_values) const;

    void applyExtrapolatedProlongation(const PolarGrid& coarse_grid, const PolarGrid& fine_grid,
                                       Vector<double> fine_result, ConstVector<double> coarse_values) const;

    /* Scaled full weighting (FW) restriction operator. */
    void applyRestriction(const PolarGrid& fine_grid, const PolarGrid& coarse_grid, Vector<double> coarse_result,
                          ConstVector<double> fine_values) const;

    void applyExtrapolatedRestriction(const PolarGrid& fine_grid, const PolarGrid& coarse_grid,
                                      Vector<double> coarse_result, ConstVector<double> fine_values) const;

    /* Bicubic FMG interpolator 1/16 * [-1, 9, 9, -1] */
    void applyFMGInterpolation(const PolarGrid& coarse_grid, const PolarGrid& fine_grid, Vector<double> fine_result,
                               ConstVector<double> coarse_values) const;

private:
    const int max_omp_threads_;
    const bool DirBC_Interior_;
};
