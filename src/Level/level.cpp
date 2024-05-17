#include "../../include/Level/level.h"

Level::Level(const int level, std::unique_ptr<const PolarGrid> grid, std::shared_ptr<const ExactFunctions> exactFunctions)
    : level_(level), grid_(std::move(grid)), exactFunctions_(std::move(exactFunctions)) {}

const PolarGrid& Level::grid() const {
    return *grid_;
}

const ExactFunctions& Level::exactFunctions() const {
    return *exactFunctions_;
}

const Operator& Level::getOperator() const {
    return *operator_;
}

int Level::level() const {
    return level_;
}

void Level::setOperator(const GMGPolar& gmgpolar) {
    operator_ = std::make_unique<const Operator>(gmgpolar, *grid_);
}

void Level::applyA(Vector<scalar_t> &result, const Vector<scalar_t> &x) const{
    operator_->applyAGive(*this, result, x);
}

void Level::applyATasks(Vector<scalar_t> &result, const Vector<scalar_t> &x) const{
    operator_->applyAGiveTasks(*this, result, x);
}

void Level::applyATake0(Vector<scalar_t> &result, const Vector<scalar_t> &x) const{
    operator_->applyATake0(*this, result, x);
}