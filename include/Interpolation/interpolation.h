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
                        ConstVector<double> x) const
    {
        assert(toLevel.level_depth() == fromLevel.level_depth() - 1);
        applyInjection_(fromLevel.grid(), toLevel.grid(), result, x, threads_per_level_[toLevel.level_depth()]);
    }

    /* Bilinear interpolation operator */
    void applyProlongation0(const Level& fromLevel, const Level& toLevel, Vector<double> result,
                            ConstVector<double> x) const
    {
        assert(toLevel.level_depth() == fromLevel.level_depth() - 1);
        applyProlongation0_(fromLevel.grid(), toLevel.grid(), result, x, threads_per_level_[toLevel.level_depth()]);
    }
    void applyProlongation(const Level& fromLevel, const Level& toLevel, Vector<double> result,
                           ConstVector<double> x) const
    {
        assert(toLevel.level_depth() == fromLevel.level_depth() - 1);
        applyProlongation_(fromLevel.grid(), toLevel.grid(), result, x, threads_per_level_[toLevel.level_depth()]);
    }
    void applyExtrapolatedProlongation0(const Level& fromLevel, const Level& toLevel, Vector<double> result,
                                        ConstVector<double> x) const
    {
        assert(toLevel.level_depth() == fromLevel.level_depth() - 1);
        applyExtrapolatedProlongation0_(fromLevel.grid(), toLevel.grid(), result, x,
                                        threads_per_level_[toLevel.level_depth()]);
    }
    void applyExtrapolatedProlongation(const Level& fromLevel, const Level& toLevel, Vector<double> result,
                                       ConstVector<double> x) const
    {
        assert(toLevel.level_depth() == fromLevel.level_depth() - 1);
        applyExtrapolatedProlongation_(fromLevel.grid(), toLevel.grid(), result, x,
                                       threads_per_level_[toLevel.level_depth()]);
    }

    /* Scaled full weighting (FW) restriction operator. */
    void applyRestriction0(const Level& fromLevel, const Level& toLevel, Vector<double> result,
                           ConstVector<double> x) const
    {
        assert(toLevel.level_depth() == fromLevel.level_depth() - 1);
        applyRestriction0_(fromLevel.grid(), toLevel.grid(), result, x, threads_per_level_[toLevel.level_depth()]);
    }
    void applyRestriction(const Level& fromLevel, const Level& toLevel, Vector<double> result,
                          ConstVector<double> x) const
    {
        assert(toLevel.level_depth() == fromLevel.level_depth() - 1);
        applyRestriction_(fromLevel.grid(), toLevel.grid(), result, x, threads_per_level_[toLevel.level_depth()]);
    }
    void applyExtrapolatedRestriction0(const Level& fromLevel, const Level& toLevel, Vector<double> result,
                                       ConstVector<double> x) const
    {
        assert(toLevel.level_depth() == fromLevel.level_depth() - 1);
        applyExtrapolatedRestriction0_(fromLevel.grid(), toLevel.grid(), result, x,
                                       threads_per_level_[toLevel.level_depth()]);
    }
    void applyExtrapolatedRestriction(const Level& fromLevel, const Level& toLevel, Vector<double> result,
                                      ConstVector<double> x) const
    {
        assert(toLevel.level_depth() == fromLevel.level_depth() - 1);
        applyExtrapolatedRestriction_(fromLevel.grid(), toLevel.grid(), result, x,
                                      threads_per_level_[toLevel.level_depth()]);
    }

    /* Bicubic FMG interpolator 1/16 * [-1, 9, 9, -1] */
    void applyFMGInterpolation(const Level& fromLevel, const Level& toLevel, Vector<double> result,
                               ConstVector<double> x) const
    {
        assert(toLevel.level_depth() == fromLevel.level_depth() - 1);
        applyFMGInterpolation_(fromLevel.grid(), toLevel.grid(), result, x, threads_per_level_[toLevel.level_depth()]);
    }

private:
    /* Remark: This injection is not scaled. */
    void applyInjection_(const PolarGrid& coarseGrid, const PolarGrid& fineGrid, Vector<double> result,
                         ConstVector<double> x, int nthreads) const;

    /* Bilinear interpolation operator */
    void applyProlongation0_(const PolarGrid& coarseGrid, const PolarGrid& fineGrid, Vector<double> result,
                             ConstVector<double> x, int nthreads) const;
    void applyProlongation_(const PolarGrid& coarseGrid, const PolarGrid& fineGrid, Vector<double> result,
                            ConstVector<double> x, int nthreads) const;
    void applyExtrapolatedProlongation0_(const PolarGrid& coarseGrid, const PolarGrid& fineGrid, Vector<double> result,
                                         ConstVector<double> x, int nthreads) const;
    void applyExtrapolatedProlongation_(const PolarGrid& coarseGrid, const PolarGrid& fineGrid, Vector<double> result,
                                        ConstVector<double> x, int nthreads) const;

    /* Scaled full weighting (FW) restriction operator. */
    void applyRestriction0_(const PolarGrid& coarseGrid, const PolarGrid& fineGrid, Vector<double> result,
                            ConstVector<double> x, int nthreads) const;
    void applyRestriction_(const PolarGrid& coarseGrid, const PolarGrid& fineGrid, Vector<double> result,
                           ConstVector<double> x, int nthreads) const;
    void applyExtrapolatedRestriction0_(const PolarGrid& coarseGrid, const PolarGrid& fineGrid, Vector<double> result,
                                        ConstVector<double> x, int nthreads) const;
    void applyExtrapolatedRestriction_(const PolarGrid& coarseGrid, const PolarGrid& fineGrid, Vector<double> result,
                                       ConstVector<double> x, int nthreads) const;

    /* Bicubic FMG interpolator 1/16 * [-1, 9, 9, -1] */
    void applyFMGInterpolation_(const PolarGrid& coarseGrid, const PolarGrid& fineGrid, Vector<double> result,
                                ConstVector<double> x, int nthreads) const;

private:
    const std::vector<int>& threads_per_level_;
    const bool DirBC_Interior_;
};
