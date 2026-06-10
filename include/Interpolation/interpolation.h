#pragma once

#include <vector>
#include <cassert>
#include <iostream>
#include <omp.h>

#include "../PolarGrid/polargrid.h"
#include "../LinearAlgebra/Vector/vector.h"
#include "../LinearAlgebra/Vector/vector_operations.h"
#include "../Definitions/global_definitions.h"

namespace gmgpolar
{

class Interpolation
{
public:
    explicit Interpolation(bool DirBC_Interior);

    /* Remark: This injection is not scaled. */
    void applyInjection(const PolarGrid<DefaultMemorySpace>& fine_grid,
                        const PolarGrid<DefaultMemorySpace>& coarse_grid, Vector<double> coarse_result,
                        ConstVector<double> fine_values) const;

    /* Bilinear interpolation operator */
    void applyProlongation(const PolarGrid<DefaultMemorySpace>& coarse_grid,
                           const PolarGrid<DefaultMemorySpace>& fine_grid, Vector<double> fine_result,
                           ConstVector<double> coarse_values) const;

    void applyExtrapolatedProlongation(const PolarGrid<DefaultMemorySpace>& coarse_grid,
                                       const PolarGrid<DefaultMemorySpace>& fine_grid, Vector<double> fine_result,
                                       ConstVector<double> coarse_values) const;

    /* Scaled full weighting (FW) restriction operator. */
    void applyRestriction(const PolarGrid<DefaultMemorySpace>& fine_grid,
                          const PolarGrid<DefaultMemorySpace>& coarse_grid, Vector<double> coarse_result,
                          ConstVector<double> fine_values) const;

    void applyExtrapolatedRestriction(const PolarGrid<DefaultMemorySpace>& fine_grid,
                                      const PolarGrid<DefaultMemorySpace>& coarse_grid, Vector<double> coarse_result,
                                      ConstVector<double> fine_values) const;

    /* Bicubic FMG interpolator 1/16 * [-1, 9, 9, -1] */
    void applyFMGInterpolation(const PolarGrid<Kokkos::HostSpace>& coarse_grid,
                               const PolarGrid<Kokkos::HostSpace>& fine_grid, HostVector<double> fine_result,
                               HostConstVector<double> coarse_values) const;

private:
    const bool DirBC_Interior_;
};
} // namespace gmgpolar
