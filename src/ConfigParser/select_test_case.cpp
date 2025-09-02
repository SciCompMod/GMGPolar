#include "../include/ConfigParser/config_parser.h"

void ConfigParser::selectTestCase(GeometryType geometry_type, ProblemType problem_type, AlphaCoeff alpha_type,
                                  BetaCoeff beta_type, double Rmax, double kappa_eps, double delta_e, double alpha_jump)
{
    /* --------------- */
    /* Domain Geometry */
    switch (geometry_type) {
    case GeometryType::CIRCULAR:
        domain_geometry_ = std::make_unique<CircularGeometry>(Rmax);
        break;

    case GeometryType::SHAFRANOV:
        domain_geometry_ = std::make_unique<ShafranovGeometry>(Rmax, kappa_eps, delta_e);
        break;

    case GeometryType::CZARNY:
        domain_geometry_ = std::make_unique<CzarnyGeometry>(Rmax, kappa_eps, delta_e);
        break;

    case GeometryType::CULHAM:
        domain_geometry_ = std::make_unique<CulhamGeometry>(Rmax);
        break;

    default:
        throw std::runtime_error("Invalid geometry.\n");
    }

    /* ---------------------------- */
    /* Density Profile Coefficients */
    switch (alpha_type) {
    case AlphaCoeff::POISSON:
        density_profile_coefficients_ = std::make_unique<PoissonCoefficients>(Rmax, alpha_jump);
        break;

    case AlphaCoeff::SONNENDRUCKER:
        switch (beta_type) {
        case BetaCoeff::ZERO:
            density_profile_coefficients_ = std::make_unique<SonnendruckerCoefficients>(Rmax, alpha_jump);
            break;
        case BetaCoeff::ALPHA_INVERSE:
            density_profile_coefficients_ = std::make_unique<SonnendruckerGyroCoefficients>(Rmax, alpha_jump);
            break;
        default:
            throw std::runtime_error("Invalid beta.\n");
        }
        break;

    case AlphaCoeff::ZONI:
        switch (beta_type) {
        case BetaCoeff::ZERO:
            density_profile_coefficients_ = std::make_unique<ZoniCoefficients>(Rmax, alpha_jump);
            break;
        case BetaCoeff::ALPHA_INVERSE:
            density_profile_coefficients_ = std::make_unique<ZoniGyroCoefficients>(Rmax, alpha_jump);
            break;
        default:
            throw std::runtime_error("Invalid beta.\n");
        }
        break;

    case AlphaCoeff::ZONI_SHIFTED:
        switch (beta_type) {
        case BetaCoeff::ZERO:
            density_profile_coefficients_ = std::make_unique<ZoniShiftedCoefficients>(Rmax, alpha_jump);
            break;
        case BetaCoeff::ALPHA_INVERSE:
            density_profile_coefficients_ = std::make_unique<ZoniShiftedGyroCoefficients>(Rmax, alpha_jump);
            break;
        default:
            throw std::runtime_error("Invalid beta.\n");
        }
        break;

    default:
        throw std::runtime_error("Invalid alpha.\n");
    }

    /* ------------------------------------ */
    /* Exact Solution & Boundary Conditions */
    switch (problem_type) {
    case ProblemType::CARTESIAN_R2:
        switch (geometry_type) {
        case GeometryType::CIRCULAR:
            exact_solution_      = std::make_unique<CartesianR2_CircularGeometry>(Rmax);
            boundary_conditions_ = std::make_unique<CartesianR2_Boundary_CircularGeometry>(Rmax);
            break;
        case GeometryType::SHAFRANOV:
            exact_solution_      = std::make_unique<CartesianR2_ShafranovGeometry>(Rmax, kappa_eps, delta_e);
            boundary_conditions_ = std::make_unique<CartesianR2_Boundary_ShafranovGeometry>(Rmax, kappa_eps, delta_e);
            break;
        case GeometryType::CZARNY:
            exact_solution_      = std::make_unique<CartesianR2_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
            boundary_conditions_ = std::make_unique<CartesianR2_Boundary_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
            break;
        default:
            throw std::runtime_error("Invalid geometry for configuration.\n");
        }
        break;

    case ProblemType::CARTESIAN_R6:
        switch (geometry_type) {
        case GeometryType::CIRCULAR:
            exact_solution_      = std::make_unique<CartesianR6_CircularGeometry>(Rmax);
            boundary_conditions_ = std::make_unique<CartesianR6_Boundary_CircularGeometry>(Rmax);
            break;
        case GeometryType::SHAFRANOV:
            exact_solution_      = std::make_unique<CartesianR6_ShafranovGeometry>(Rmax, kappa_eps, delta_e);
            boundary_conditions_ = std::make_unique<CartesianR6_Boundary_ShafranovGeometry>(Rmax, kappa_eps, delta_e);
            break;
        case GeometryType::CZARNY:
            exact_solution_      = std::make_unique<CartesianR6_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
            boundary_conditions_ = std::make_unique<CartesianR6_Boundary_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
            break;
        default:
            throw std::runtime_error("Invalid geometry for configuration.\n");
        }
        break;

    case ProblemType::POLAR_R6:
        switch (geometry_type) {
        case GeometryType::CIRCULAR:
            exact_solution_      = std::make_unique<PolarR6_CircularGeometry>(Rmax);
            boundary_conditions_ = std::make_unique<PolarR6_Boundary_CircularGeometry>(Rmax);
            break;
        case GeometryType::SHAFRANOV:
            exact_solution_      = std::make_unique<PolarR6_ShafranovGeometry>(Rmax, kappa_eps, delta_e);
            boundary_conditions_ = std::make_unique<PolarR6_Boundary_ShafranovGeometry>(Rmax, kappa_eps, delta_e);
            break;
        case GeometryType::CZARNY:
            exact_solution_      = std::make_unique<PolarR6_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
            boundary_conditions_ = std::make_unique<PolarR6_Boundary_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
            break;
        case GeometryType::CULHAM:
            exact_solution_      = std::make_unique<PolarR6_CulhamGeometry>(Rmax);
            boundary_conditions_ = std::make_unique<PolarR6_Boundary_CulhamGeometry>(Rmax);
            break;
        default:
            throw std::runtime_error("Invalid geometry for configuration.\n");
        }
        break;

    case ProblemType::REFINED_RADIUS:
        switch (geometry_type) {
        case GeometryType::CIRCULAR:
            exact_solution_      = std::make_unique<Refined_CircularGeometry>(Rmax);
            boundary_conditions_ = std::make_unique<Refined_Boundary_CircularGeometry>(Rmax);
            break;
        case GeometryType::SHAFRANOV:
            exact_solution_      = std::make_unique<Refined_ShafranovGeometry>(Rmax, kappa_eps, delta_e);
            boundary_conditions_ = std::make_unique<Refined_Boundary_ShafranovGeometry>(Rmax, kappa_eps, delta_e);
            break;
        case GeometryType::CZARNY:
            exact_solution_      = std::make_unique<Refined_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
            boundary_conditions_ = std::make_unique<Refined_Boundary_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
            break;
        case GeometryType::CULHAM:
            exact_solution_      = std::make_unique<Refined_CulhamGeometry>(Rmax);
            boundary_conditions_ = std::make_unique<Refined_Boundary_CulhamGeometry>(Rmax);
            break;
        default:
            throw std::runtime_error("Invalid geometry for configuration.\n");
        }
        break;

    default:
        throw std::runtime_error("Invalid problem.\n");
    }

    /* ------------------- */
    /* Source Term (rhs_f) */
    switch (problem_type) {
    case ProblemType::CARTESIAN_R2:

        switch (geometry_type) {
        case GeometryType::CIRCULAR:

            switch (alpha_type) {
            case AlphaCoeff::POISSON:
                source_term_ = std::make_unique<CartesianR2_Poisson_CircularGeometry>(Rmax);
                break;
            case AlphaCoeff::SONNENDRUCKER:
                switch (beta_type) {
                case BetaCoeff::ZERO:
                    source_term_ = std::make_unique<CartesianR2_Sonnendrucker_CircularGeometry>(Rmax);
                    break;
                case BetaCoeff::ALPHA_INVERSE:
                    source_term_ = std::make_unique<CartesianR2_SonnendruckerGyro_CircularGeometry>(Rmax);
                    break;
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            case AlphaCoeff::ZONI:
                switch (beta_type) {
                case BetaCoeff::ZERO:
                    source_term_ = std::make_unique<CartesianR2_Zoni_CircularGeometry>(Rmax);
                    break;
                case BetaCoeff::ALPHA_INVERSE:
                    source_term_ = std::make_unique<CartesianR2_ZoniGyro_CircularGeometry>(Rmax);
                    break;
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            case AlphaCoeff::ZONI_SHIFTED:
                switch (beta_type) {
                case BetaCoeff::ZERO:
                    source_term_ = std::make_unique<CartesianR2_ZoniShifted_CircularGeometry>(Rmax);
                    break;
                case BetaCoeff::ALPHA_INVERSE:
                    source_term_ = std::make_unique<CartesianR2_ZoniShiftedGyro_CircularGeometry>(Rmax);
                    break;
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            default:
                throw std::runtime_error("Invalid alpha.\n");
            }
            break;

        case GeometryType::SHAFRANOV:

            switch (alpha_type) {
            case AlphaCoeff::POISSON:
                source_term_ = std::make_unique<CartesianR2_Poisson_ShafranovGeometry>(Rmax, kappa_eps, delta_e);
                break;
            case AlphaCoeff::SONNENDRUCKER:
                switch (beta_type) {
                case BetaCoeff::ZERO:
                    source_term_ =
                        std::make_unique<CartesianR2_Sonnendrucker_ShafranovGeometry>(Rmax, kappa_eps, delta_e);
                    break;
                case BetaCoeff::ALPHA_INVERSE:
                    source_term_ =
                        std::make_unique<CartesianR2_SonnendruckerGyro_ShafranovGeometry>(Rmax, kappa_eps, delta_e);
                    break;
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            case AlphaCoeff::ZONI:
                switch (beta_type) {
                case BetaCoeff::ZERO:
                    source_term_ = std::make_unique<CartesianR2_Zoni_ShafranovGeometry>(Rmax, kappa_eps, delta_e);
                    break;
                case BetaCoeff::ALPHA_INVERSE:
                    source_term_ = std::make_unique<CartesianR2_ZoniGyro_ShafranovGeometry>(Rmax, kappa_eps, delta_e);
                    break;
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            case AlphaCoeff::ZONI_SHIFTED:
                switch (beta_type) {
                case BetaCoeff::ZERO:
                    source_term_ =
                        std::make_unique<CartesianR2_ZoniShifted_ShafranovGeometry>(Rmax, kappa_eps, delta_e);
                    break;
                case BetaCoeff::ALPHA_INVERSE:
                    source_term_ =
                        std::make_unique<CartesianR2_ZoniShiftedGyro_ShafranovGeometry>(Rmax, kappa_eps, delta_e);
                    break;
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            default:
                throw std::runtime_error("Invalid alpha.\n");
            }
            break;

        case GeometryType::CZARNY:

            switch (alpha_type) {
            case AlphaCoeff::POISSON:
                source_term_ = std::make_unique<CartesianR2_Poisson_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
                break;
            case AlphaCoeff::SONNENDRUCKER:
                switch (beta_type) {
                case BetaCoeff::ZERO:
                    source_term_ = std::make_unique<CartesianR2_Sonnendrucker_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
                    break;
                case BetaCoeff::ALPHA_INVERSE:
                    source_term_ =
                        std::make_unique<CartesianR2_SonnendruckerGyro_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
                    break;
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            case AlphaCoeff::ZONI:
                switch (beta_type) {
                case BetaCoeff::ZERO:
                    source_term_ = std::make_unique<CartesianR2_Zoni_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
                    break;
                case BetaCoeff::ALPHA_INVERSE:
                    source_term_ = std::make_unique<CartesianR2_ZoniGyro_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
                    break;
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            case AlphaCoeff::ZONI_SHIFTED:
                switch (beta_type) {
                case BetaCoeff::ZERO:
                    source_term_ = std::make_unique<CartesianR2_ZoniShifted_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
                    break;
                case BetaCoeff::ALPHA_INVERSE:
                    source_term_ =
                        std::make_unique<CartesianR2_ZoniShiftedGyro_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
                    break;
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            default:
                throw std::runtime_error("Invalid alpha.\n");
            }
            break;

        default:
            throw std::runtime_error("Invalid geometry for configuration.\n");
        }
        break;

    case ProblemType::CARTESIAN_R6:

        switch (geometry_type) {
        case GeometryType::CIRCULAR:

            switch (alpha_type) {
            case AlphaCoeff::POISSON:
                source_term_ = std::make_unique<CartesianR6_Poisson_CircularGeometry>(Rmax);
                break;
            case AlphaCoeff::SONNENDRUCKER:
                switch (beta_type) {
                case BetaCoeff::ZERO:
                    source_term_ = std::make_unique<CartesianR6_Sonnendrucker_CircularGeometry>(Rmax);
                    break;
                case BetaCoeff::ALPHA_INVERSE:
                    source_term_ = std::make_unique<CartesianR6_SonnendruckerGyro_CircularGeometry>(Rmax);
                    break;
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            case AlphaCoeff::ZONI:
                switch (beta_type) {
                case BetaCoeff::ZERO:
                    source_term_ = std::make_unique<CartesianR6_Zoni_CircularGeometry>(Rmax);
                    break;
                case BetaCoeff::ALPHA_INVERSE:
                    source_term_ = std::make_unique<CartesianR6_ZoniGyro_CircularGeometry>(Rmax);
                    break;
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            case AlphaCoeff::ZONI_SHIFTED:
                switch (beta_type) {
                case BetaCoeff::ZERO:
                    source_term_ = std::make_unique<CartesianR6_ZoniShifted_CircularGeometry>(Rmax);
                    break;
                case BetaCoeff::ALPHA_INVERSE:
                    source_term_ = std::make_unique<CartesianR6_ZoniShiftedGyro_CircularGeometry>(Rmax);
                    break;
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            default:
                throw std::runtime_error("Invalid alpha.\n");
            }
            break;

        case GeometryType::SHAFRANOV:

            switch (alpha_type) {
            case AlphaCoeff::POISSON:
                source_term_ = std::make_unique<CartesianR6_Poisson_ShafranovGeometry>(Rmax, kappa_eps, delta_e);
                break;
            case AlphaCoeff::SONNENDRUCKER:
                switch (beta_type) {
                case BetaCoeff::ZERO:
                    source_term_ =
                        std::make_unique<CartesianR6_Sonnendrucker_ShafranovGeometry>(Rmax, kappa_eps, delta_e);
                    break;
                case BetaCoeff::ALPHA_INVERSE:
                    source_term_ =
                        std::make_unique<CartesianR6_SonnendruckerGyro_ShafranovGeometry>(Rmax, kappa_eps, delta_e);
                    break;
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            case AlphaCoeff::ZONI:
                switch (beta_type) {
                case BetaCoeff::ZERO:
                    source_term_ = std::make_unique<CartesianR6_Zoni_ShafranovGeometry>(Rmax, kappa_eps, delta_e);
                    break;
                case BetaCoeff::ALPHA_INVERSE:
                    source_term_ = std::make_unique<CartesianR6_ZoniGyro_ShafranovGeometry>(Rmax, kappa_eps, delta_e);
                    break;
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            case AlphaCoeff::ZONI_SHIFTED:
                switch (beta_type) {
                case BetaCoeff::ZERO:
                    source_term_ =
                        std::make_unique<CartesianR6_ZoniShifted_ShafranovGeometry>(Rmax, kappa_eps, delta_e);
                    break;
                case BetaCoeff::ALPHA_INVERSE:
                    source_term_ =
                        std::make_unique<CartesianR6_ZoniShiftedGyro_ShafranovGeometry>(Rmax, kappa_eps, delta_e);
                    break;
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            default:
                throw std::runtime_error("Invalid alpha.\n");
            }
            break;

        case GeometryType::CZARNY:

            switch (alpha_type) {
            case AlphaCoeff::POISSON:
                source_term_ = std::make_unique<CartesianR6_Poisson_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
                break;
            case AlphaCoeff::SONNENDRUCKER:
                switch (beta_type) {
                case BetaCoeff::ZERO:
                    source_term_ = std::make_unique<CartesianR6_Sonnendrucker_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
                    break;
                case BetaCoeff::ALPHA_INVERSE:
                    source_term_ =
                        std::make_unique<CartesianR6_SonnendruckerGyro_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
                    break;
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            case AlphaCoeff::ZONI:
                switch (beta_type) {
                case BetaCoeff::ZERO:
                    source_term_ = std::make_unique<CartesianR6_Zoni_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
                    break;
                case BetaCoeff::ALPHA_INVERSE:
                    source_term_ = std::make_unique<CartesianR6_ZoniGyro_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
                    break;
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            case AlphaCoeff::ZONI_SHIFTED:
                switch (beta_type) {
                case BetaCoeff::ZERO:
                    source_term_ = std::make_unique<CartesianR6_ZoniShifted_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
                    break;
                case BetaCoeff::ALPHA_INVERSE:
                    source_term_ =
                        std::make_unique<CartesianR6_ZoniShiftedGyro_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
                    break;
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            default:
                throw std::runtime_error("Invalid alpha.\n");
            }
            break;

        default:
            throw std::runtime_error("Invalid geometry for configuration.\n");
        }
        break;

    case ProblemType::POLAR_R6:

        switch (geometry_type) {
        case GeometryType::CIRCULAR:

            switch (alpha_type) {
            case AlphaCoeff::POISSON:
                source_term_ = std::make_unique<PolarR6_Poisson_CircularGeometry>(Rmax);
                break;
            case AlphaCoeff::SONNENDRUCKER:
                switch (beta_type) {
                case BetaCoeff::ZERO:
                    source_term_ = std::make_unique<PolarR6_Sonnendrucker_CircularGeometry>(Rmax);
                    break;
                case BetaCoeff::ALPHA_INVERSE:
                    source_term_ = std::make_unique<PolarR6_SonnendruckerGyro_CircularGeometry>(Rmax);
                    break;
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            case AlphaCoeff::ZONI:
                switch (beta_type) {
                case BetaCoeff::ZERO:
                    source_term_ = std::make_unique<PolarR6_Zoni_CircularGeometry>(Rmax);
                    break;
                case BetaCoeff::ALPHA_INVERSE:
                    source_term_ = std::make_unique<PolarR6_ZoniGyro_CircularGeometry>(Rmax);
                    break;
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            case AlphaCoeff::ZONI_SHIFTED:
                switch (beta_type) {
                case BetaCoeff::ZERO:
                    source_term_ = std::make_unique<PolarR6_ZoniShifted_CircularGeometry>(Rmax);
                    break;
                case BetaCoeff::ALPHA_INVERSE:
                    source_term_ = std::make_unique<PolarR6_ZoniShiftedGyro_CircularGeometry>(Rmax);
                    break;
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            default:
                throw std::runtime_error("Invalid alpha.\n");
            }
            break;

        case GeometryType::SHAFRANOV:

            switch (alpha_type) {
            case AlphaCoeff::POISSON:
                source_term_ = std::make_unique<PolarR6_Poisson_ShafranovGeometry>(Rmax, kappa_eps, delta_e);
                break;
            case AlphaCoeff::SONNENDRUCKER:
                switch (beta_type) {
                case BetaCoeff::ZERO:
                    source_term_ = std::make_unique<PolarR6_Sonnendrucker_ShafranovGeometry>(Rmax, kappa_eps, delta_e);
                    break;
                case BetaCoeff::ALPHA_INVERSE:
                    source_term_ =
                        std::make_unique<PolarR6_SonnendruckerGyro_ShafranovGeometry>(Rmax, kappa_eps, delta_e);
                    break;
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            case AlphaCoeff::ZONI:
                switch (beta_type) {
                case BetaCoeff::ZERO:
                    source_term_ = std::make_unique<PolarR6_Zoni_ShafranovGeometry>(Rmax, kappa_eps, delta_e);
                    break;
                case BetaCoeff::ALPHA_INVERSE:
                    source_term_ = std::make_unique<PolarR6_ZoniGyro_ShafranovGeometry>(Rmax, kappa_eps, delta_e);
                    break;
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            case AlphaCoeff::ZONI_SHIFTED:
                switch (beta_type) {
                case BetaCoeff::ZERO:
                    source_term_ = std::make_unique<PolarR6_ZoniShifted_ShafranovGeometry>(Rmax, kappa_eps, delta_e);
                    break;
                case BetaCoeff::ALPHA_INVERSE:
                    source_term_ =
                        std::make_unique<PolarR6_ZoniShiftedGyro_ShafranovGeometry>(Rmax, kappa_eps, delta_e);
                    break;
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            default:
                throw std::runtime_error("Invalid alpha.\n");
            }
            break;

        case GeometryType::CZARNY:

            switch (alpha_type) {
            case AlphaCoeff::POISSON:
                source_term_ = std::make_unique<PolarR6_Poisson_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
                break;
            case AlphaCoeff::SONNENDRUCKER:
                switch (beta_type) {
                case BetaCoeff::ZERO:
                    source_term_ = std::make_unique<PolarR6_Sonnendrucker_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
                    break;
                case BetaCoeff::ALPHA_INVERSE:
                    source_term_ = std::make_unique<PolarR6_SonnendruckerGyro_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
                    break;
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            case AlphaCoeff::ZONI:
                switch (beta_type) {
                case BetaCoeff::ZERO:
                    source_term_ = std::make_unique<PolarR6_Zoni_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
                    break;
                case BetaCoeff::ALPHA_INVERSE:
                    source_term_ = std::make_unique<PolarR6_ZoniGyro_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
                    break;
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            case AlphaCoeff::ZONI_SHIFTED:
                switch (beta_type) {
                case BetaCoeff::ZERO:
                    source_term_ = std::make_unique<PolarR6_ZoniShifted_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
                    break;
                case BetaCoeff::ALPHA_INVERSE:
                    source_term_ = std::make_unique<PolarR6_ZoniShiftedGyro_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
                    break;
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            default:
                throw std::runtime_error("Invalid alpha.\n");
            }
            break;

        case GeometryType::CULHAM:
            switch (alpha_type) {
            case AlphaCoeff::ZONI_SHIFTED:
                switch (beta_type) {
                case BetaCoeff::ALPHA_INVERSE:
                    source_term_ = std::make_unique<PolarR6_ZoniShiftedGyro_CulhamGeometry>(Rmax);
                    break;
                default:
                    throw std::runtime_error("Invalid beta for configuration.\n");
                }
                break;
            default:
                throw std::runtime_error("Invalid alpha for configuration.\n");
            }
            break;

        default:
            throw std::runtime_error("Invalid geometry for configuration.\n");
        }
        break;

    case ProblemType::REFINED_RADIUS:

        switch (geometry_type) {
        case GeometryType::CIRCULAR:

            switch (alpha_type) {
            case AlphaCoeff::ZONI_SHIFTED:
                switch (beta_type) {
                case BetaCoeff::ALPHA_INVERSE:
                    source_term_ = std::make_unique<Refined_ZoniShiftedGyro_CircularGeometry>(Rmax);
                    break;
                default:
                    throw std::runtime_error("Invalid beta for configuration.\n");
                }
                break;
            default:
                throw std::runtime_error("Invalid alpha for configuration.\n");
            }
            break;

        case GeometryType::SHAFRANOV:

            switch (alpha_type) {
            case AlphaCoeff::ZONI_SHIFTED:
                switch (beta_type) {
                case BetaCoeff::ALPHA_INVERSE:
                    source_term_ =
                        std::make_unique<Refined_ZoniShiftedGyro_ShafranovGeometry>(Rmax, kappa_eps, delta_e);
                    break;
                default:
                    throw std::runtime_error("Invalid beta for configuration.\n");
                }
                break;
            default:
                throw std::runtime_error("Invalid alpha for configuration.\n");
            }
            break;

        case GeometryType::CZARNY:

            switch (alpha_type) {
            case AlphaCoeff::ZONI_SHIFTED:
                switch (beta_type) {
                case BetaCoeff::ALPHA_INVERSE:
                    source_term_ = std::make_unique<Refined_ZoniShiftedGyro_CzarnyGeometry>(Rmax, kappa_eps, delta_e);
                    break;
                default:
                    throw std::runtime_error("Invalid beta for configuration.\n");
                }
                break;
            default:
                throw std::runtime_error("Invalid alpha for configuration.\n");
            }
            break;

        case GeometryType::CULHAM:

            switch (alpha_type) {
            case AlphaCoeff::ZONI_SHIFTED:
                switch (beta_type) {
                case BetaCoeff::ALPHA_INVERSE:
                    source_term_ = std::make_unique<Refined_ZoniShiftedGyro_CulhamGeometry>(Rmax);
                    break;
                default:
                    throw std::runtime_error("Invalid beta for configuration.\n");
                }
                break;
            default:
                throw std::runtime_error("Invalid alpha for configuration.\n");
            }
            break;

        default:
            throw std::runtime_error("Invalid geometry for configuration.\n");
        }
        break;

    default:
        throw std::runtime_error("Invalid problem.\n");
    }
}