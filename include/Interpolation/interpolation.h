#pragma once

#include <vector>
#include <cassert>
#include <iostream>

#include <omp.h>

#include "../PolarGrid/polargrid.h"
#include "../Level/level.h"

#include "../LinearAlgebra/vector.h"
#include "../LinearAlgebra/vector_operations.h"

#include "../common/global_definitions.h"

class Interpolation
{
public:
    explicit Interpolation(const std::vector<int>& threads_per_level, const bool DirBC_Interior);

    /* Remark: This injection is not scaled. */
    void applyInjection(const Level& fromLevel, const Level& toLevel, Vector<double> result,
                        ConstVector<double> x) const;

    /* Bilinear interpolation operator */
    void applyProlongation0(const Level& fromLevel, const Level& toLevel, Vector<double> result,
                            ConstVector<double> x) const;
    void applyProlongation(const Level& fromLevel, const Level& toLevel, Vector<double> result,
                           ConstVector<double> x) const;
    void applyExtrapolatedProlongation0(const Level& fromLevel, const Level& toLevel, Vector<double> result,
                                        ConstVector<double> x) const;
    void applyExtrapolatedProlongation(const Level& fromLevel, const Level& toLevel, Vector<double> result,
                                       ConstVector<double> x) const;

    /* Scaled full weighting (FW) restriction operator. */
    void applyRestriction0(const Level& fromLevel, const Level& toLevel, Vector<double> result,
                           ConstVector<double> x) const;
    void applyRestriction(const Level& fromLevel, const Level& toLevel, Vector<double> result,
                          ConstVector<double> x) const;
    void applyExtrapolatedRestriction0(const Level& fromLevel, const Level& toLevel, Vector<double> result,
                                       ConstVector<double> x) const;
    void applyExtrapolatedRestriction(const Level& fromLevel, const Level& toLevel, Vector<double> result,
                                      ConstVector<double> x) const;

    /* Bicubic FMG interpolator 1/16 * [-1, 9, 9, -1] */
    void applyFMGInterpolation(const Level& fromLevel, const Level& toLevel, Vector<double> result,
                               ConstVector<double> x) const;

private:
    const std::vector<int>& threads_per_level_;
    const bool DirBC_Interior_;
};