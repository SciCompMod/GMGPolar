#pragma once

class Operator;
class GMGPolar;
class ExactFunctions;
class PolarGrid;

#include "../include/Operator/operator.h"

#include "../common/scalar.h"
#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include <memory>

class Level {
public:
    Level(const int level, std::unique_ptr<const PolarGrid> grid, std::shared_ptr<const ExactFunctions> exactFunctions);
    void setOperator(const GMGPolar& gmgpolar);

    const PolarGrid& grid() const;
    const ExactFunctions& exactFunctions() const;
    const Operator& getOperator() const;
    int level() const;

    void applyA(Vector<scalar_t> &result, const Vector<scalar_t> &x, const scalar_t& scaleAx) const;
    void applyATasks(Vector<scalar_t> &result, const Vector<scalar_t> &x, const scalar_t& scaleAx) const;
    void applyAMutex(Vector<scalar_t> &result, const Vector<scalar_t> &x, const scalar_t& scaleAx);
    void applyATake0(Vector<scalar_t> &result, const Vector<scalar_t> &x, const scalar_t& scaleAx) const;

    void applySmoothing(Vector<scalar_t> &result, const Vector<scalar_t> &x) const;

    std::pair<SparseMatrix<scalar_t>, Vector<scalar_t>> build_system() const;

private:
    const int level_;
    std::unique_ptr<const PolarGrid> grid_;
    std::shared_ptr<const ExactFunctions> exactFunctions_;
    std::unique_ptr<Operator> operator_;
};
