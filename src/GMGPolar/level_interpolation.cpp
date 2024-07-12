#include "../../include/GMGPolar/gmgpolar.h"

void GMGPolar::prolongateToUpperLevel(const int current_level, Vector<double>& result, const Vector<double>& x) const {
    assert(current_level < numberOflevels_ && 1 <= current_level);
    if(!interpolation_) throw std::runtime_error("Interpolation not initialized.");

    interpolation_->applyProlongation(levels_[current_level], levels_[current_level-1], result, x);
}

void GMGPolar::restrictToLowerLevel(const int current_level, Vector<double>& result, const Vector<double>& x) const {
    assert(current_level < numberOflevels_ - 1 && 0 <= current_level);
    if(!interpolation_) throw std::runtime_error("Interpolation not initialized.");

    interpolation_->applyRestriction(levels_[current_level], levels_[current_level+1], result, x);
}

void GMGPolar::injectToLowerLevel(const int current_level, Vector<double>& result, const Vector<double>& x) const {
    assert(current_level < numberOflevels_ - 1 && 0 <= current_level);
    if(!interpolation_) throw std::runtime_error("Interpolation not initialized.");
    
    interpolation_->applyInjection(levels_[current_level], levels_[current_level+1], result, x);
}