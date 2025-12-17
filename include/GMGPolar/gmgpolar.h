#pragma once

#include <filesystem>
#include <iostream>
#include <omp.h>
#include <utility>

#include "../InputFunctions/densityProfileCoefficients.h"
#include "../InputFunctions/domainGeometry.h"

#include "igmgpolar.h"

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
class GMGPolar : public IGMGPolar
{
public:
    /* ---------------------------------------------------------------------- */
    /* Constructor & Initialization                                           */
    /* ---------------------------------------------------------------------- */
    // Construct a polar PDE multigrid solver for the Poisson-like equation:
    // - \nabla \cdot (\alpha \nabla u) + \beta u = rhs_f  in \Omega,
    // with Dirichlet boundary condition        u = u_D    on \partial \Omega.
    // Parameters:
    // - grid: Cartesian mesh discretizing the computational domain.
    // - domain_geometry: Mapping from the reference domain to the physical domain \Omega.
    // - density_profile_coefficients: Coefficients \alpha and \beta defining the PDE.
    GMGPolar(const PolarGrid& grid, const DomainGeometry& domain_geometry,
             const DensityProfileCoefficients& density_profile_coefficients)
        : IGMGPolar(grid)
        , domain_geometry_(domain_geometry)
        , density_profile_coefficients_(density_profile_coefficients)
    {
    }

    /* ---------------------------------------------------------------------- */
    /* Setup & Solve                                                          */
    /* ---------------------------------------------------------------------- */
    // Finalize solver setup (allocate data, build operators, etc.).
    void setup() final;

private:
    /* ------------------------------------ */
    /* Grid Configuration & Input Functions */
    /* ------------------------------------ */
    const DomainGeometry& domain_geometry_;
    const DensityProfileCoefficients& density_profile_coefficients_;

    /* --------------- */
    /* Setup Functions */
    void discretize_rhs_f(const Level& level, Vector<double> rhs_f) final;

    /* ------------- */
    /* Visualization */
    void writeToVTK(const std::filesystem::path& file_path, const PolarGrid& grid) final;
    void writeToVTK(const std::filesystem::path& file_path, const Level& level,
                    ConstVector<double> grid_function) final;
};

#include "build_rhs_f.h"
#include "setup.h"
#include "writeToVTK.h"
