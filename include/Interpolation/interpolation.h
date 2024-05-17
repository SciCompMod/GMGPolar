#pragma once

#include <vector>
#include <cassert>

#include <omp.h>

#include "../include/PolarGrid/polargrid.h"
#include "../include/Level/level.h"
#include "../common/scalar.h"
#include "../include/TaskDistribution/TaskDistribution.h"
#include "../include/linear_algebra/vector.h"
#include "../include/linear_algebra/operations.h"

class Interpolation {
public:
    void applyProlongation0(const Level& fromLevel, const Level& toLevel, Vector<scalar_t>& result, const Vector<scalar_t>& x) const;
    void applyProlongation(const Level& fromLevel, const Level& toLevel, Vector<scalar_t>& result, const Vector<scalar_t>& x) const;

    void applyRestrictionTake0(const Level& fromLevel, const Level& toLevel, Vector<scalar_t>& result, const Vector<scalar_t>& x) const;
    void applyRestrictionTake(const Level& fromLevel, const Level& toLevel, Vector<scalar_t>& result, const Vector<scalar_t>& x) const;
    void applyRestrictionTakeTasks(const Level& fromLevel, const Level& toLevel, Vector<scalar_t>& result, const Vector<scalar_t>& x) const;

    void applyRestrictionGive0(const Level& fromLevel, const Level& toLevel, Vector<scalar_t>& result, const Vector<scalar_t>& x) const;
    void applyRestrictionGive(const Level& fromLevel, const Level& toLevel, Vector<scalar_t>& result, const Vector<scalar_t>& x) const;
    void applyRestrictionGiveTasks(const Level& fromLevel, const Level& toLevel, Vector<scalar_t>& result, const Vector<scalar_t>& x) const;
};
