#include "../../include/GMGPolar/gmgpolar.h"

void IGMGPolar::prolongation(int current_level, Vector<double> result, ConstVector<double> x) const
{
    assert(current_level < number_of_levels_ && 1 <= current_level);
    if (!interpolation_)
        throw std::runtime_error("Interpolation not initialized.");

    interpolation_->applyProlongation(levels_[current_level].grid(), levels_[current_level - 1].grid(), result, x);
}

void IGMGPolar::restriction(int current_level, Vector<double> result, ConstVector<double> x) const
{
    assert(current_level < number_of_levels_ - 1 && 0 <= current_level);
    if (!interpolation_)
        throw std::runtime_error("Interpolation not initialized.");

    interpolation_->applyRestriction(levels_[current_level].grid(), levels_[current_level + 1].grid(), result, x);
}

void IGMGPolar::injection(int current_level, Vector<double> result, ConstVector<double> x) const
{
    assert(current_level < number_of_levels_ - 1 && 0 <= current_level);
    if (!interpolation_)
        throw std::runtime_error("Interpolation not initialized.");

    interpolation_->applyInjection(levels_[current_level].grid(), levels_[current_level + 1].grid(), result, x);
}

void IGMGPolar::extrapolatedProlongation(int current_level, Vector<double> result, ConstVector<double> x) const
{
    assert(current_level < number_of_levels_ && 1 <= current_level);
    if (!interpolation_)
        throw std::runtime_error("Interpolation not initialized.");

    interpolation_->applyExtrapolatedProlongation(levels_[current_level].grid(), levels_[current_level - 1].grid(),
                                                  result, x);
}

void IGMGPolar::extrapolatedRestriction(int current_level, Vector<double> result, ConstVector<double> x) const
{
    assert(current_level < number_of_levels_ - 1 && 0 <= current_level);
    if (!interpolation_)
        throw std::runtime_error("Interpolation not initialized.");

    interpolation_->applyExtrapolatedRestriction(levels_[current_level].grid(), levels_[current_level + 1].grid(),
                                                 result, x);
}

void IGMGPolar::FMGInterpolation(int current_level, Vector<double> result, ConstVector<double> x) const
{
    assert(current_level < number_of_levels_ && 1 <= current_level);
    if (!interpolation_)
        throw std::runtime_error("Interpolation not initialized.");

    interpolation_->applyFMGInterpolation(levels_[current_level].grid(), levels_[current_level - 1].grid(), result, x);
}
