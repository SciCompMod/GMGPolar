#include "../../include/ConfigParser/config_parser.h"
#include "../../include/GMGPolar/gmgpolar.h"

std::unique_ptr<IGMGPolar> ConfigParser::solver() const
{
    // Create local aliases so the class doesn't need to be captured by the lamda
    // These are references, not copies.
    const PolarGrid& grid = grid_;

    // Create a solver specialized to the active domain geometry.
    return std::visit(
        [&grid](auto const& domain_geometry, auto const& density_profile_coefficients) {
            using DomainGeomType                 = std::decay_t<decltype(domain_geometry)>;
            using DensityProfileCoefficientsType = std::decay_t<decltype(density_profile_coefficients)>;

            // Construct the solver specialized for this geometry type.
            std::unique_ptr<IGMGPolar> solver =
                std::make_unique<GMGPolar<DomainGeomType, DensityProfileCoefficientsType>>(
                    grid, domain_geometry, density_profile_coefficients);

            // The lambdas must return objects of identical type
            return solver;
        },
        *domain_geometry_, *density_profile_coefficients_);
}

void ConfigParser::selectTestCase(GeometryType geometry_type, ProblemType problem_type, AlphaCoeff alpha_type,
                                  BetaCoeff beta_type, double Rmax, double kappa_eps, double delta_e, double alpha_jump)
{
    geometry_type_ = geometry_type;
    problem_type_  = problem_type;
    alpha_type_    = alpha_type;
    beta_type_     = beta_type;

    /* --------------- */
    /* Domain Geometry */
    switch (geometry_type) {
    case GeometryType::CIRCULAR:
        domain_geometry_ = std::make_unique<DomainGeometryVariant>(CircularGeometry(Rmax));
        break;

    case GeometryType::SHAFRANOV:
        domain_geometry_ = std::make_unique<DomainGeometryVariant>(ShafranovGeometry(Rmax, kappa_eps, delta_e));
        break;

    case GeometryType::CZARNY:
        domain_geometry_ = std::make_unique<DomainGeometryVariant>(CzarnyGeometry(Rmax, kappa_eps, delta_e));
        break;

    case GeometryType::CULHAM:
        domain_geometry_ = std::make_unique<DomainGeometryVariant>(CulhamGeometry(Rmax));
        break;

    default:
        throw std::runtime_error("Invalid geometry.\n");
    }

    /* ---------------------------- */
    /* Density Profile Coefficients */
    switch (alpha_type) {
    case AlphaCoeff::POISSON:
        density_profile_coefficients_ =
            std::make_unique<DensityProfileCoefficientsVariant>(PoissonCoefficients(Rmax, alpha_jump));
        break;

    case AlphaCoeff::SONNENDRUCKER:
        switch (beta_type) {
        case BetaCoeff::ZERO:
            density_profile_coefficients_ =
                std::make_unique<DensityProfileCoefficientsVariant>(SonnendruckerCoefficients(Rmax, alpha_jump));
            break;
        case BetaCoeff::ALPHA_INVERSE:
            density_profile_coefficients_ =
                std::make_unique<DensityProfileCoefficientsVariant>(SonnendruckerGyroCoefficients(Rmax, alpha_jump));
            break;
        default:
            throw std::runtime_error("Invalid beta.\n");
        }
        break;

    case AlphaCoeff::ZONI:
        switch (beta_type) {
        case BetaCoeff::ZERO:
            density_profile_coefficients_ =
                std::make_unique<DensityProfileCoefficientsVariant>(ZoniCoefficients(Rmax, alpha_jump));
            break;
        case BetaCoeff::ALPHA_INVERSE:
            density_profile_coefficients_ =
                std::make_unique<DensityProfileCoefficientsVariant>(ZoniGyroCoefficients(Rmax, alpha_jump));
            break;
        default:
            throw std::runtime_error("Invalid beta.\n");
        }
        break;

    case AlphaCoeff::ZONI_SHIFTED:
        switch (beta_type) {
        case BetaCoeff::ZERO:
            density_profile_coefficients_ =
                std::make_unique<DensityProfileCoefficientsVariant>(ZoniShiftedCoefficients(Rmax, alpha_jump));
            break;
        case BetaCoeff::ALPHA_INVERSE:
            density_profile_coefficients_ =
                std::make_unique<DensityProfileCoefficientsVariant>(ZoniShiftedGyroCoefficients(Rmax, alpha_jump));
            break;
        default:
            throw std::runtime_error("Invalid beta.\n");
        }
        break;

    default:
        throw std::runtime_error("Invalid alpha.\n");
    }

    /* -------------- */
    /* Exact Solution */
    switch (problem_type) {
    case ProblemType::CARTESIAN_R2:
        switch (geometry_type) {
        case GeometryType::CIRCULAR:
            exact_solution_ = std::make_unique<CartesianR2_CircularGeometry>(Rmax);
            break;
        case GeometryType::SHAFRANOV:
            exact_solution_ = std::make_unique<CartesianR2_ShafranovGeometry>(Rmax, kappa_eps, delta_e);
            break;
        case GeometryType::CZARNY:
            exact_solution_ = std::make_unique<CartesianR2_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
            break;
        default:
            throw std::runtime_error("Invalid geometry for configuration.\n");
        }
        break;

    case ProblemType::CARTESIAN_R6:
        switch (geometry_type) {
        case GeometryType::CIRCULAR:
            exact_solution_ = std::make_unique<CartesianR6_CircularGeometry>(Rmax);
            break;
        case GeometryType::SHAFRANOV:
            exact_solution_ = std::make_unique<CartesianR6_ShafranovGeometry>(Rmax, kappa_eps, delta_e);
            break;
        case GeometryType::CZARNY:
            exact_solution_ = std::make_unique<CartesianR6_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
            break;
        default:
            throw std::runtime_error("Invalid geometry for configuration.\n");
        }
        break;

    case ProblemType::POLAR_R6:
        switch (geometry_type) {
        case GeometryType::CIRCULAR:
            exact_solution_ = std::make_unique<PolarR6_CircularGeometry>(Rmax);
            break;
        case GeometryType::SHAFRANOV:
            exact_solution_ = std::make_unique<PolarR6_ShafranovGeometry>(Rmax, kappa_eps, delta_e);
            break;
        case GeometryType::CZARNY:
            exact_solution_ = std::make_unique<PolarR6_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
            break;
        case GeometryType::CULHAM:
            exact_solution_ = std::make_unique<PolarR6_CulhamGeometry>(Rmax);
            break;
        default:
            throw std::runtime_error("Invalid geometry for configuration.\n");
        }
        break;

    case ProblemType::REFINED_RADIUS:
        switch (geometry_type) {
        case GeometryType::CIRCULAR:
            exact_solution_ = std::make_unique<Refined_CircularGeometry>(Rmax);
            break;
        case GeometryType::SHAFRANOV:
            exact_solution_ = std::make_unique<Refined_ShafranovGeometry>(Rmax, kappa_eps, delta_e);
            break;
        case GeometryType::CZARNY:
            exact_solution_ = std::make_unique<Refined_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
            break;
        case GeometryType::CULHAM:
            exact_solution_ = std::make_unique<Refined_CulhamGeometry>(Rmax);
            break;
        default:
            throw std::runtime_error("Invalid geometry for configuration.\n");
        }
        break;

    default:
        throw std::runtime_error("Invalid problem.\n");
    }
}

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
void ConfigParser::solve(GMGPolar<DomainGeometry, DensityProfileCoefficients>& solver) const
{
    assert(domain_geometry_ != nullptr);
    /* ------------------------------------------ */
    /* Source Term (rhs_f) and BoundaryConditions */
    switch (problem_type_) {
    case ProblemType::CARTESIAN_R2:

        switch (geometry_type_) {
        case GeometryType::CIRCULAR: {
            CartesianR2_Boundary_CircularGeometry boundary_conditions(Rmax_);

            switch (alpha_type_) {
            case AlphaCoeff::POISSON: {
                CartesianR2_Poisson_CircularGeometry source_term(grid_, Rmax_);
                solver.solve(boundary_conditions, source_term);
                break;
            }
            case AlphaCoeff::SONNENDRUCKER:
                switch (beta_type_) {
                case BetaCoeff::ZERO: {
                    CartesianR2_Sonnendrucker_CircularGeometry source_term(grid_, Rmax_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                case BetaCoeff::ALPHA_INVERSE: {
                    CartesianR2_SonnendruckerGyro_CircularGeometry source_term(grid_, Rmax_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            case AlphaCoeff::ZONI:
                switch (beta_type_) {
                case BetaCoeff::ZERO: {
                    CartesianR2_Zoni_CircularGeometry source_term(grid_, Rmax_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                case BetaCoeff::ALPHA_INVERSE: {
                    CartesianR2_ZoniGyro_CircularGeometry source_term(grid_, Rmax_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            case AlphaCoeff::ZONI_SHIFTED:
                switch (beta_type_) {
                case BetaCoeff::ZERO: {
                    CartesianR2_ZoniShifted_CircularGeometry source_term(grid_, Rmax_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                case BetaCoeff::ALPHA_INVERSE: {
                    CartesianR2_ZoniShiftedGyro_CircularGeometry source_term(grid_, Rmax_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            default:
                throw std::runtime_error("Invalid alpha.\n");
            }
            break;
        }

        case GeometryType::SHAFRANOV: {
            CartesianR2_Boundary_ShafranovGeometry boundary_conditions(Rmax_, kappa_eps_, delta_e_);

            switch (alpha_type_) {
            case AlphaCoeff::POISSON: {
                CartesianR2_Poisson_ShafranovGeometry source_term(grid_, Rmax_, kappa_eps_, delta_e_);
                solver.solve(boundary_conditions, source_term);
                break;
            }
            case AlphaCoeff::SONNENDRUCKER:
                switch (beta_type_) {
                case BetaCoeff::ZERO: {
                    CartesianR2_Sonnendrucker_ShafranovGeometry source_term(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                case BetaCoeff::ALPHA_INVERSE: {
                    CartesianR2_SonnendruckerGyro_ShafranovGeometry source_term(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            case AlphaCoeff::ZONI:
                switch (beta_type_) {
                case BetaCoeff::ZERO: {
                    CartesianR2_Zoni_ShafranovGeometry source_term(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                case BetaCoeff::ALPHA_INVERSE: {
                    CartesianR2_ZoniGyro_ShafranovGeometry source_term(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            case AlphaCoeff::ZONI_SHIFTED:
                switch (beta_type_) {
                case BetaCoeff::ZERO: {
                    CartesianR2_ZoniShifted_ShafranovGeometry source_term(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                case BetaCoeff::ALPHA_INVERSE: {
                    CartesianR2_ZoniShiftedGyro_ShafranovGeometry source_term(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            default:
                throw std::runtime_error("Invalid alpha.\n");
            }
            break;
        }

        case GeometryType::CZARNY: {
            CartesianR2_Boundary_CzarnyGeometry boundary_conditions(Rmax_, kappa_eps_, delta_e_);

            switch (alpha_type_) {
            case AlphaCoeff::POISSON: {
                CartesianR2_Poisson_CzarnyGeometry source_term(grid_, Rmax_, kappa_eps_, delta_e_);
                solver.solve(boundary_conditions, source_term);
                break;
            }
            case AlphaCoeff::SONNENDRUCKER:
                switch (beta_type_) {
                case BetaCoeff::ZERO: {
                    CartesianR2_Sonnendrucker_CzarnyGeometry source_term(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                case BetaCoeff::ALPHA_INVERSE: {
                    CartesianR2_SonnendruckerGyro_CzarnyGeometry source_term(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            case AlphaCoeff::ZONI:
                switch (beta_type_) {
                case BetaCoeff::ZERO: {
                    CartesianR2_Zoni_CzarnyGeometry source_term(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                case BetaCoeff::ALPHA_INVERSE: {
                    CartesianR2_ZoniGyro_CzarnyGeometry source_term(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            case AlphaCoeff::ZONI_SHIFTED:
                switch (beta_type_) {
                case BetaCoeff::ZERO: {
                    CartesianR2_ZoniShifted_CzarnyGeometry source_term(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                case BetaCoeff::ALPHA_INVERSE: {
                    CartesianR2_ZoniShiftedGyro_CzarnyGeometry source_term(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            default:
                throw std::runtime_error("Invalid alpha.\n");
            }
            break;
        }

        default:
            throw std::runtime_error("Invalid geometry for configuration.\n");
        }
        break;

    case ProblemType::CARTESIAN_R6:

        switch (geometry_type_) {
        case GeometryType::CIRCULAR: {
            CartesianR6_Boundary_CircularGeometry boundary_conditions(Rmax_);

            switch (alpha_type_) {
            case AlphaCoeff::POISSON: {
                CartesianR6_Poisson_CircularGeometry source_term(grid_, Rmax_);
                solver.solve(boundary_conditions, source_term);
                break;
            }
            case AlphaCoeff::SONNENDRUCKER:
                switch (beta_type_) {
                case BetaCoeff::ZERO: {
                    CartesianR6_Sonnendrucker_CircularGeometry source_term(grid_, Rmax_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                case BetaCoeff::ALPHA_INVERSE: {
                    CartesianR6_SonnendruckerGyro_CircularGeometry source_term(grid_, Rmax_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            case AlphaCoeff::ZONI:
                switch (beta_type_) {
                case BetaCoeff::ZERO: {
                    CartesianR6_Zoni_CircularGeometry source_term(grid_, Rmax_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                case BetaCoeff::ALPHA_INVERSE: {
                    CartesianR6_ZoniGyro_CircularGeometry source_term(grid_, Rmax_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            case AlphaCoeff::ZONI_SHIFTED:
                switch (beta_type_) {
                case BetaCoeff::ZERO: {
                    CartesianR6_ZoniShifted_CircularGeometry source_term(grid_, Rmax_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                case BetaCoeff::ALPHA_INVERSE: {
                    CartesianR6_ZoniShiftedGyro_CircularGeometry source_term(grid_, Rmax_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            default:
                throw std::runtime_error("Invalid alpha.\n");
            }
            break;
        }

        case GeometryType::SHAFRANOV: {
            CartesianR6_Boundary_ShafranovGeometry boundary_conditions(Rmax_, kappa_eps_, delta_e_);

            switch (alpha_type_) {
            case AlphaCoeff::POISSON: {
                CartesianR6_Poisson_ShafranovGeometry source_term(grid_, Rmax_, kappa_eps_, delta_e_);
                solver.solve(boundary_conditions, source_term);
                break;
            }
            case AlphaCoeff::SONNENDRUCKER:
                switch (beta_type_) {
                case BetaCoeff::ZERO: {
                    CartesianR6_Sonnendrucker_ShafranovGeometry source_term(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                case BetaCoeff::ALPHA_INVERSE: {
                    CartesianR6_SonnendruckerGyro_ShafranovGeometry source_term(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            case AlphaCoeff::ZONI:
                switch (beta_type_) {
                case BetaCoeff::ZERO: {
                    CartesianR6_Zoni_ShafranovGeometry source_term(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                case BetaCoeff::ALPHA_INVERSE: {
                    CartesianR6_ZoniGyro_ShafranovGeometry source_term(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            case AlphaCoeff::ZONI_SHIFTED:
                switch (beta_type_) {
                case BetaCoeff::ZERO: {
                    CartesianR6_ZoniShifted_ShafranovGeometry source_term(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                case BetaCoeff::ALPHA_INVERSE: {
                    CartesianR6_ZoniShiftedGyro_ShafranovGeometry source_term(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            default:
                throw std::runtime_error("Invalid alpha.\n");
            }
            break;
        }

        case GeometryType::CZARNY: {
            CartesianR6_Boundary_CzarnyGeometry boundary_conditions(Rmax_, kappa_eps_, delta_e_);

            switch (alpha_type_) {
            case AlphaCoeff::POISSON: {
                CartesianR6_Poisson_CzarnyGeometry source_term(grid_, Rmax_, kappa_eps_, delta_e_);
                solver.solve(boundary_conditions, source_term);
                break;
            }
            case AlphaCoeff::SONNENDRUCKER:
                switch (beta_type_) {
                case BetaCoeff::ZERO: {
                    CartesianR6_Sonnendrucker_CzarnyGeometry source_term(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                case BetaCoeff::ALPHA_INVERSE: {
                    CartesianR6_SonnendruckerGyro_CzarnyGeometry source_term(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            case AlphaCoeff::ZONI:
                switch (beta_type_) {
                case BetaCoeff::ZERO: {
                    CartesianR6_Zoni_CzarnyGeometry source_term(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                case BetaCoeff::ALPHA_INVERSE: {
                    CartesianR6_ZoniGyro_CzarnyGeometry source_term(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            case AlphaCoeff::ZONI_SHIFTED:
                switch (beta_type_) {
                case BetaCoeff::ZERO: {
                    CartesianR6_ZoniShifted_CzarnyGeometry source_term(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                case BetaCoeff::ALPHA_INVERSE: {
                    CartesianR6_ZoniShiftedGyro_CzarnyGeometry source_term(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            default:
                throw std::runtime_error("Invalid alpha.\n");
            }
            break;
        }

        default:
            throw std::runtime_error("Invalid geometry for configuration.\n");
        }
        break;

    case ProblemType::POLAR_R6:

        switch (geometry_type_) {
        case GeometryType::CIRCULAR: {
            PolarR6_Boundary_CircularGeometry boundary_conditions(Rmax_);

            switch (alpha_type_) {
            case AlphaCoeff::POISSON: {
                PolarR6_Poisson_CircularGeometry source_term(grid_, Rmax_);
                solver.solve(boundary_conditions, source_term);
                break;
            }
            case AlphaCoeff::SONNENDRUCKER:
                switch (beta_type_) {
                case BetaCoeff::ZERO: {
                    PolarR6_Sonnendrucker_CircularGeometry source_term(grid_, Rmax_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                case BetaCoeff::ALPHA_INVERSE: {
                    PolarR6_SonnendruckerGyro_CircularGeometry source_term(grid_, Rmax_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            case AlphaCoeff::ZONI:
                switch (beta_type_) {
                case BetaCoeff::ZERO: {
                    PolarR6_Zoni_CircularGeometry source_term(grid_, Rmax_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                case BetaCoeff::ALPHA_INVERSE: {
                    PolarR6_ZoniGyro_CircularGeometry source_term(grid_, Rmax_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            case AlphaCoeff::ZONI_SHIFTED:
                switch (beta_type_) {
                case BetaCoeff::ZERO: {
                    PolarR6_ZoniShifted_CircularGeometry source_term(grid_, Rmax_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                case BetaCoeff::ALPHA_INVERSE: {
                    PolarR6_ZoniShiftedGyro_CircularGeometry source_term(grid_, Rmax_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            default:
                throw std::runtime_error("Invalid alpha.\n");
            }
            break;
        }

        case GeometryType::SHAFRANOV: {
            PolarR6_Boundary_ShafranovGeometry boundary_conditions(Rmax_, kappa_eps_, delta_e_);

            switch (alpha_type_) {
            case AlphaCoeff::POISSON: {
                PolarR6_Poisson_ShafranovGeometry source_term(grid_, Rmax_, kappa_eps_, delta_e_);
                solver.solve(boundary_conditions, source_term);
                break;
            }
            case AlphaCoeff::SONNENDRUCKER:
                switch (beta_type_) {
                case BetaCoeff::ZERO: {
                    PolarR6_Sonnendrucker_ShafranovGeometry source_term(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                case BetaCoeff::ALPHA_INVERSE: {
                    PolarR6_SonnendruckerGyro_ShafranovGeometry source_term(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            case AlphaCoeff::ZONI:
                switch (beta_type_) {
                case BetaCoeff::ZERO: {
                    PolarR6_Zoni_ShafranovGeometry source_term(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                case BetaCoeff::ALPHA_INVERSE: {
                    PolarR6_ZoniGyro_ShafranovGeometry source_term(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            case AlphaCoeff::ZONI_SHIFTED:
                switch (beta_type_) {
                case BetaCoeff::ZERO: {
                    PolarR6_ZoniShifted_ShafranovGeometry source_term(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                case BetaCoeff::ALPHA_INVERSE: {
                    PolarR6_ZoniShiftedGyro_ShafranovGeometry source_term(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            default:
                throw std::runtime_error("Invalid alpha.\n");
            }
            break;
        }

        case GeometryType::CZARNY: {
            PolarR6_Boundary_CzarnyGeometry boundary_conditions(Rmax_, kappa_eps_, delta_e_);

            switch (alpha_type_) {
            case AlphaCoeff::POISSON: {
                PolarR6_Poisson_CzarnyGeometry source_term(grid_, Rmax_, kappa_eps_, delta_e_);
                solver.solve(boundary_conditions, source_term);
                break;
            }
            case AlphaCoeff::SONNENDRUCKER:
                switch (beta_type_) {
                case BetaCoeff::ZERO: {
                    PolarR6_Sonnendrucker_CzarnyGeometry source_term(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                case BetaCoeff::ALPHA_INVERSE: {
                    PolarR6_SonnendruckerGyro_CzarnyGeometry source_term(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            case AlphaCoeff::ZONI:
                switch (beta_type_) {
                case BetaCoeff::ZERO: {
                    PolarR6_Zoni_CzarnyGeometry source_term(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                case BetaCoeff::ALPHA_INVERSE: {
                    PolarR6_ZoniGyro_CzarnyGeometry source_term(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            case AlphaCoeff::ZONI_SHIFTED:
                switch (beta_type_) {
                case BetaCoeff::ZERO: {
                    PolarR6_ZoniShifted_CzarnyGeometry source_term(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                case BetaCoeff::ALPHA_INVERSE: {
                    PolarR6_ZoniShiftedGyro_CzarnyGeometry source_term(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            default:
                throw std::runtime_error("Invalid alpha.\n");
            }
            break;
        }

        case GeometryType::CULHAM: {
            PolarR6_Boundary_CulhamGeometry boundary_conditions(Rmax_);
            switch (alpha_type_) {
            case AlphaCoeff::ZONI_SHIFTED:
                switch (beta_type_) {
                case BetaCoeff::ALPHA_INVERSE: {
                    PolarR6_ZoniShiftedGyro_CulhamGeometry source_term(grid_, Rmax_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                default:
                    throw std::runtime_error("Invalid beta for configuration.\n");
                }
                break;
            default:
                throw std::runtime_error("Invalid alpha for configuration.\n");
            }
            break;
        }

        default:
            throw std::runtime_error("Invalid geometry for configuration.\n");
        }
        break;

    case ProblemType::REFINED_RADIUS:

        switch (geometry_type_) {
        case GeometryType::CIRCULAR: {
            Refined_Boundary_CircularGeometry boundary_conditions(Rmax_);

            switch (alpha_type_) {
            case AlphaCoeff::ZONI_SHIFTED:
                switch (beta_type_) {
                case BetaCoeff::ALPHA_INVERSE: {
                    Refined_ZoniShiftedGyro_CircularGeometry source_term(grid_, Rmax_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                default:
                    throw std::runtime_error("Invalid beta for configuration.\n");
                }
                break;
            default:
                throw std::runtime_error("Invalid alpha for configuration.\n");
            }
            break;
        }

        case GeometryType::SHAFRANOV: {
            Refined_Boundary_ShafranovGeometry boundary_conditions(Rmax_, kappa_eps_, delta_e_);

            switch (alpha_type_) {
            case AlphaCoeff::ZONI_SHIFTED:
                switch (beta_type_) {
                case BetaCoeff::ALPHA_INVERSE: {
                    Refined_ZoniShiftedGyro_ShafranovGeometry source_term(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                default:
                    throw std::runtime_error("Invalid beta for configuration.\n");
                }
                break;
            default:
                throw std::runtime_error("Invalid alpha for configuration.\n");
            }
            break;
        }

        case GeometryType::CZARNY: {
            Refined_Boundary_CzarnyGeometry boundary_conditions(Rmax_, kappa_eps_, delta_e_);

            switch (alpha_type_) {
            case AlphaCoeff::ZONI_SHIFTED:
                switch (beta_type_) {
                case BetaCoeff::ALPHA_INVERSE: {
                    Refined_ZoniShiftedGyro_CzarnyGeometry source_term(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                default:
                    throw std::runtime_error("Invalid beta for configuration.\n");
                }
                break;
            default:
                throw std::runtime_error("Invalid alpha for configuration.\n");
            }
            break;
        }

        case GeometryType::CULHAM: {
            Refined_Boundary_CulhamGeometry boundary_conditions(Rmax_);

            switch (alpha_type_) {
            case AlphaCoeff::ZONI_SHIFTED:
                switch (beta_type_) {
                case BetaCoeff::ALPHA_INVERSE: {
                    Refined_ZoniShiftedGyro_CulhamGeometry source_term(grid_, Rmax_);
                    solver.solve(boundary_conditions, source_term);
                    break;
                }
                default:
                    throw std::runtime_error("Invalid beta for configuration.\n");
                }
                break;
            default:
                throw std::runtime_error("Invalid alpha for configuration.\n");
            }
            break;
        }

        default:
            throw std::runtime_error("Invalid geometry for configuration.\n");
        }
        break;

    default:
        throw std::runtime_error("Invalid problem.\n");
    }
}

template void ConfigParser::solve<CircularGeometry, PoissonCoefficients>(
    GMGPolar<CircularGeometry, PoissonCoefficients>& solver) const;
template void ConfigParser::solve<CircularGeometry, SonnendruckerCoefficients>(
    GMGPolar<CircularGeometry, SonnendruckerCoefficients>& solver) const;
template void ConfigParser::solve<CircularGeometry, SonnendruckerGyroCoefficients>(
    GMGPolar<CircularGeometry, SonnendruckerGyroCoefficients>& solver) const;
template void
ConfigParser::solve<CircularGeometry, ZoniCoefficients>(GMGPolar<CircularGeometry, ZoniCoefficients>& solver) const;
template void ConfigParser::solve<CircularGeometry, ZoniGyroCoefficients>(
    GMGPolar<CircularGeometry, ZoniGyroCoefficients>& solver) const;
template void ConfigParser::solve<CircularGeometry, ZoniShiftedCoefficients>(
    GMGPolar<CircularGeometry, ZoniShiftedCoefficients>& solver) const;
template void ConfigParser::solve<CircularGeometry, ZoniShiftedGyroCoefficients>(
    GMGPolar<CircularGeometry, ZoniShiftedGyroCoefficients>& solver) const;

template void ConfigParser::solve<ShafranovGeometry, PoissonCoefficients>(
    GMGPolar<ShafranovGeometry, PoissonCoefficients>& solver) const;
template void ConfigParser::solve<ShafranovGeometry, SonnendruckerCoefficients>(
    GMGPolar<ShafranovGeometry, SonnendruckerCoefficients>& solver) const;
template void ConfigParser::solve<ShafranovGeometry, SonnendruckerGyroCoefficients>(
    GMGPolar<ShafranovGeometry, SonnendruckerGyroCoefficients>& solver) const;
template void
ConfigParser::solve<ShafranovGeometry, ZoniCoefficients>(GMGPolar<ShafranovGeometry, ZoniCoefficients>& solver) const;
template void ConfigParser::solve<ShafranovGeometry, ZoniGyroCoefficients>(
    GMGPolar<ShafranovGeometry, ZoniGyroCoefficients>& solver) const;
template void ConfigParser::solve<ShafranovGeometry, ZoniShiftedCoefficients>(
    GMGPolar<ShafranovGeometry, ZoniShiftedCoefficients>& solver) const;
template void ConfigParser::solve<ShafranovGeometry, ZoniShiftedGyroCoefficients>(
    GMGPolar<ShafranovGeometry, ZoniShiftedGyroCoefficients>& solver) const;

template void
ConfigParser::solve<CzarnyGeometry, PoissonCoefficients>(GMGPolar<CzarnyGeometry, PoissonCoefficients>& solver) const;
template void ConfigParser::solve<CzarnyGeometry, SonnendruckerCoefficients>(
    GMGPolar<CzarnyGeometry, SonnendruckerCoefficients>& solver) const;
template void ConfigParser::solve<CzarnyGeometry, SonnendruckerGyroCoefficients>(
    GMGPolar<CzarnyGeometry, SonnendruckerGyroCoefficients>& solver) const;
template void
ConfigParser::solve<CzarnyGeometry, ZoniCoefficients>(GMGPolar<CzarnyGeometry, ZoniCoefficients>& solver) const;
template void
ConfigParser::solve<CzarnyGeometry, ZoniGyroCoefficients>(GMGPolar<CzarnyGeometry, ZoniGyroCoefficients>& solver) const;
template void ConfigParser::solve<CzarnyGeometry, ZoniShiftedCoefficients>(
    GMGPolar<CzarnyGeometry, ZoniShiftedCoefficients>& solver) const;
template void ConfigParser::solve<CzarnyGeometry, ZoniShiftedGyroCoefficients>(
    GMGPolar<CzarnyGeometry, ZoniShiftedGyroCoefficients>& solver) const;

template void
ConfigParser::solve<CulhamGeometry, PoissonCoefficients>(GMGPolar<CulhamGeometry, PoissonCoefficients>& solver) const;
template void ConfigParser::solve<CulhamGeometry, SonnendruckerCoefficients>(
    GMGPolar<CulhamGeometry, SonnendruckerCoefficients>& solver) const;
template void ConfigParser::solve<CulhamGeometry, SonnendruckerGyroCoefficients>(
    GMGPolar<CulhamGeometry, SonnendruckerGyroCoefficients>& solver) const;
template void
ConfigParser::solve<CulhamGeometry, ZoniCoefficients>(GMGPolar<CulhamGeometry, ZoniCoefficients>& solver) const;
template void
ConfigParser::solve<CulhamGeometry, ZoniGyroCoefficients>(GMGPolar<CulhamGeometry, ZoniGyroCoefficients>& solver) const;
template void ConfigParser::solve<CulhamGeometry, ZoniShiftedCoefficients>(
    GMGPolar<CulhamGeometry, ZoniShiftedCoefficients>& solver) const;
template void ConfigParser::solve<CulhamGeometry, ZoniShiftedGyroCoefficients>(
    GMGPolar<CulhamGeometry, ZoniShiftedGyroCoefficients>& solver) const;
