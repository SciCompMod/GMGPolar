#include "../../include/ConfigParser/config_parser.h"
#include "../../include/GMGPolar/gmgpolar.h"
using namespace gmgpolar;

std::unique_ptr<IGMGPolar> ConfigParser::solver() const
{
    // Create local aliases so the class doesn't need to be captured by the lamda
    // These are references, not copies.
    const PolarGrid<Kokkos::HostSpace>& grid = grid_;

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
    Rmax_          = Rmax;
    kappa_eps_     = kappa_eps;
    delta_e_       = delta_e;
    alpha_jump_    = alpha_jump;

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
