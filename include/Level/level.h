#pragma once

class Operator;
class GMGPolar;
class ExactFunctions;
class PolarGrid;

#include "../include/Operator/operator.h"

#include "../common/scalar.h"
#include "../linear_algebra/vector.h"

#include <memory>

class Level {
public:
    Level(const int level, std::unique_ptr<const PolarGrid> grid, std::shared_ptr<const ExactFunctions> exactFunctions);
    void setOperator(const GMGPolar& gmgpolar);

    const PolarGrid& grid() const;
    const ExactFunctions& exactFunctions() const;
    const Operator& getOperator() const;
    int level() const;

    void applyA(Vector<scalar_t> &result, const Vector<scalar_t> &x) const;
    void applyATasks(Vector<scalar_t> &result, const Vector<scalar_t> &x) const;
    void applyATake0(Vector<scalar_t> &result, const Vector<scalar_t> &x) const;
private:
    const int level_;
    std::unique_ptr<const PolarGrid> grid_;
    std::shared_ptr<const ExactFunctions> exactFunctions_;
    std::unique_ptr<const Operator> operator_;
};
