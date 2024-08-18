#pragma once

#include <vector>
#include <cassert>
#include <iostream>

#include <omp.h>

#include "../PolarGrid/polargrid.h"
#include "../Level/level.h"

#include "../LinearAlgebra/vector.h"
#include "../LinearAlgebra/vector_operations.h"

#include "../common/constants.h"
#include "../TaskDistribution/taskDistribution.h"

class Interpolation {
public:
    explicit Interpolation(const std::vector<int>& threads_per_level);

    void applyInjection(const Level& fromLevel, const Level& toLevel, Vector<double>& result, const Vector<double>& x) const;

    void applyProlongation0(const Level& fromLevel, const Level& toLevel, Vector<double>& result, const Vector<double>& x) const;
    void applyProlongation(const Level& fromLevel, const Level& toLevel, Vector<double>& result, const Vector<double>& x) const;
    void applyExtrapolatedProlongation0(const Level& fromLevel, const Level& toLevel, Vector<double>& result, const Vector<double>& x) const;
    void applyExtrapolatedProlongation(const Level& fromLevel, const Level& toLevel, Vector<double>& result, const Vector<double>& x) const;

    void applyRestriction0(const Level& fromLevel, const Level& toLevel, Vector<double>& result, const Vector<double>& x) const;
    void applyRestriction(const Level& fromLevel, const Level& toLevel, Vector<double>& result, const Vector<double>& x) const;
    void applyExtrapolatedRestriction0(const Level& fromLevel, const Level& toLevel, Vector<double>& result, const Vector<double>& x) const;
    void applyExtrapolatedRestriction(const Level& fromLevel, const Level& toLevel, Vector<double>& result, const Vector<double>& x) const;

private:
    const std::vector<int>& threads_per_level_;
};