#pragma once

#include <vector>
#include <cassert>
#include <iostream>

#include <omp.h>

#include "../PolarGrid/polargrid.h"
#include "../Level/level.h"

#include "../LinearAlgebra/vector.h"
#include "../LinearAlgebra/operations.h"

#include "../common/constants.h"
#include "../TaskDistribution/taskDistribution.h"

class Interpolation {
public:
    explicit Interpolation(const int maxOpenMPThreads, const std::vector<int>& taskingThreads);

    void applyProlongation0(const Level& fromLevel, const Level& toLevel, Vector<double>& result, const Vector<double>& x) const;
    void applyProlongation(const Level& fromLevel, const Level& toLevel, Vector<double>& result, const Vector<double>& x) const;

    void applyRestrictionTake0(const Level& fromLevel, const Level& toLevel, Vector<double>& result, const Vector<double>& x) const;
    void applyRestrictionTake(const Level& fromLevel, const Level& toLevel, Vector<double>& result, const Vector<double>& x) const;
    void applyRestrictionTakeTasks(const Level& fromLevel, const Level& toLevel, Vector<double>& result, const Vector<double>& x) const;

    void applyRestrictionGive0(const Level& fromLevel, const Level& toLevel, Vector<double>& result, const Vector<double>& x) const;
    void applyRestrictionGive(const Level& fromLevel, const Level& toLevel, Vector<double>& result, const Vector<double>& x) const;
    void applyRestrictionGiveTasks(const Level& fromLevel, const Level& toLevel, Vector<double>& result, const Vector<double>& x) const;
    
private:
    const int maxOpenMPThreads_;
    const std::vector<int>& taskingThreads_;
};