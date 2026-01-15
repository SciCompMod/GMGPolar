#pragma once
#include <random>

#include "../include/PolarGrid/polargrid.h"
#include "../include/LinearAlgebra/vector.h"

Vector<double> generate_random_sample_data(const PolarGrid& grid, unsigned int seed)
{
    Vector<double> x("x", grid.numberOfNodes());
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(-100.0, 100.0);
    for (uint i = 0; i < x.size(); ++i) {
        x(i) = dist(gen);
    }
    return x;
}

