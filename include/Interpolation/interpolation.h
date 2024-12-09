#pragma once

#include <vector>
#include <cassert>
#include <iostream>

#include <omp.h>

#include "../PolarGrid/polargrid.h"
#include "../Level/level.h"

#include "../LinearAlgebra/Vector/vector.h"

#include "../LinearAlgebra/Vector/gpu_vector.h"

#include "../common/constants.h"

class Interpolation {
public:
    explicit Interpolation(const bool DirBC_Interior);

    /* Remark: This injection is not scaled. */
    void applyInjection(const Level& fromLevel, const Level& toLevel, Vector<double>& result, const Vector<double>& x) const;
    void applyInjection(const Level& fromLevel, const Level& toLevel, GPU_Vector<double>& result, const GPU_Vector<double>& x) const;

    /* Bilinear interpolation operator */
    void applyProlongation(const Level& fromLevel, const Level& toLevel, Vector<double>& result, const Vector<double>& x) const;
    void applyProlongation(const Level& fromLevel, const Level& toLevel, GPU_Vector<double>& result, const GPU_Vector<double>& x) const;
    void applyExtrapolatedProlongation(const Level& fromLevel, const Level& toLevel, Vector<double>& result, const Vector<double>& x) const;
    void applyExtrapolatedProlongation(const Level& fromLevel, const Level& toLevel, GPU_Vector<double>& result, const GPU_Vector<double>& x) const;

    /* Scaled full weighting (FW) restriction operator. */
    void applyRestriction(const Level& fromLevel, const Level& toLevel, Vector<double>& result, const Vector<double>& x) const;
    void applyRestriction(const Level& fromLevel, const Level& toLevel, GPU_Vector<double>& result, const GPU_Vector<double>& x) const;
    void applyExtrapolatedRestriction(const Level& fromLevel, const Level& toLevel, Vector<double>& result, const Vector<double>& x) const;
    void applyExtrapolatedRestriction(const Level& fromLevel, const Level& toLevel, GPU_Vector<double>& result, const GPU_Vector<double>& x) const;

    /* Bicubic FMG interpolator 1/16 * [-1, 9, 9, -1] */
    void applyFMGInterpolation(const Level& fromLevel, const Level& toLevel, Vector<double>& result, const Vector<double>& x) const;
    void applyFMGInterpolation(const Level& fromLevel, const Level& toLevel, GPU_Vector<double>& result, const GPU_Vector<double>& x) const;

private:
    const bool DirBC_Interior_;
};