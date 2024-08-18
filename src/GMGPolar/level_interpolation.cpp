#include "../../include/GMGPolar/gmgpolar.h"

void GMGPolar::prolongation(const int current_level, Vector<double>& result, const Vector<double>& x) const {
    assert(current_level < number_of_levels_ && 1 <= current_level);
    if(!interpolation_) throw std::runtime_error("Interpolation not initialized.");

    interpolation_->applyProlongation(levels_[current_level], levels_[current_level-1], result, x);
}

void GMGPolar::restriction(const int current_level, Vector<double>& result, const Vector<double>& x) const {
    assert(current_level < number_of_levels_ - 1 && 0 <= current_level);
    if(!interpolation_) throw std::runtime_error("Interpolation not initialized.");

    interpolation_->applyRestriction(levels_[current_level], levels_[current_level+1], result, x);
}

void GMGPolar::injection(const int current_level, Vector<double>& result, const Vector<double>& x) const {
    assert(current_level < number_of_levels_ - 1 && 0 <= current_level);
    if(!interpolation_) throw std::runtime_error("Interpolation not initialized.");
    
    interpolation_->applyInjection(levels_[current_level], levels_[current_level+1], result, x);
}

void GMGPolar::extrapolatedProlongation(const int current_level, Vector<double>& result, const Vector<double>& x) const {
    assert(current_level < number_of_levels_ && 1 <= current_level);
    if(!interpolation_) throw std::runtime_error("Interpolation not initialized.");

    interpolation_->applyExtrapolatedProlongation(levels_[current_level], levels_[current_level-1], result, x);
}

void GMGPolar::extrapolatedRestriction(const int current_level, Vector<double>& result, const Vector<double>& x) const {
    assert(current_level < number_of_levels_ - 1 && 0 <= current_level);
    if(!interpolation_) throw std::runtime_error("Interpolation not initialized.");

    interpolation_->applyExtrapolatedRestriction(levels_[current_level], levels_[current_level+1], result, x);
}