#pragma once
#include <random>

#include "../include/PolarGrid/polargrid.h"
#include "../include/LinearAlgebra/Vector/vector.h"

using namespace gmgpolar;

inline HostVector<double> generate_random_sample_data(const PolarGrid<Kokkos::HostSpace>& grid, unsigned int seed,
                                                      double min_val = -100.0, double max_val = 100.0)
{
    HostVector<double> x("x", grid.numberOfNodes());
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(min_val, max_val);
    for (std::size_t i = 0; i < x.size(); ++i) {
        x(i) = dist(gen);
    }
    return x;
}

inline Vector<double> generate_random_sample_data(const PolarGrid<DefaultMemorySpace>& grid, unsigned int seed,
                                                  double min_val = -100.0, double max_val = 100.0)
{
    HostVector<double> h_x = generate_random_sample_data(PolarGrid<Kokkos::HostSpace>(grid), seed, min_val, max_val);
    return Kokkos::create_mirror_view_and_copy(DefaultMemorySpace(), h_x);
}
