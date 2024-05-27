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
    operator_ = std::make_unique<Operator>(gmgpolar, *grid_);
}

void Level::applyA(Vector<scalar_t> &result, const Vector<scalar_t> &x, const scalar_t& scaleAx) const{
    operator_->applyAGive(*this, result, x, scaleAx);
}

void Level::applyATasks(Vector<scalar_t> &result, const Vector<scalar_t> &x, const scalar_t& scaleAx) const{
    operator_->applyAGiveTasks(*this, result, x, scaleAx);
}

void Level::applyATake0(Vector<scalar_t> &result, const Vector<scalar_t> &x, const scalar_t& scaleAx) const{
    operator_->applyATake0(*this, result, x, scaleAx);
}

void Level::applyAMutex(Vector<scalar_t> &result, const Vector<scalar_t> &x, const scalar_t& scaleAx){
    operator_->applyAGiveMutex(*this, result, x, scaleAx);
}


void Level::applySmoothing(Vector<scalar_t> &result, const Vector<scalar_t> &x) const{
    operator_->applySmoothing(*this, result, x);
}

std::pair<SparseMatrix<scalar_t>, Vector<scalar_t>> Level::build_system() const{
    return operator_->build_system(*this);
}