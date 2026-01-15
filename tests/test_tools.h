#pragma once
#include <random>

#include "../include/PolarGrid/polargrid.h"
#include "../include/LinearAlgebra/vector.h"

inline Vector<double> generate_random_sample_data(const PolarGrid& grid, unsigned int seed, double min_val = -100.0,
                                                  double max_val = 100.0)
{
    Vector<double> x("x", grid.numberOfNodes());
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(min_val, max_val);
    for (uint i = 0; i < x.size(); ++i) {
        x(i) = dist(gen);
    }
    return x;
}
