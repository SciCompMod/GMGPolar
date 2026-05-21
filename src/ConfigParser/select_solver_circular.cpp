#include <ConfigParser/config_parser.h>
#include <GMGPolar/gmgpolar.h>

using namespace gmgpolar;

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
void ConfigParser::solve(GMGPolar<DomainGeometry, DensityProfileCoefficients>& solver) const
{
    if constexpr (std::is_same_v<DomainGeometry, CircularGeometry>) {

        switch (problem_type_) {

        /* ------------------------------------------------------------------ */
        case ProblemType::CARTESIAN_R2: {
            CartesianR2_Boundary_CircularGeometry bc(Rmax_);
            switch (alpha_type_) {
            case AlphaCoeff::POISSON: {
                CartesianR2_Poisson_CircularGeometry src(grid_, Rmax_);
                solver.solve(bc, src);
                break;
            }
            case AlphaCoeff::SONNENDRUCKER:
                switch (beta_type_) {
                case BetaCoeff::ZERO: {
                    CartesianR2_Sonnendrucker_CircularGeometry src(grid_, Rmax_);
                    solver.solve(bc, src);
                    break;
                }
                case BetaCoeff::ALPHA_INVERSE: {
                    CartesianR2_SonnendruckerGyro_CircularGeometry src(grid_, Rmax_);
                    solver.solve(bc, src);
                    break;
                }
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            case AlphaCoeff::ZONI:
                switch (beta_type_) {
                case BetaCoeff::ZERO: {
                    CartesianR2_Zoni_CircularGeometry src(grid_, Rmax_);
                    solver.solve(bc, src);
                    break;
                }
                case BetaCoeff::ALPHA_INVERSE: {
                    CartesianR2_ZoniGyro_CircularGeometry src(grid_, Rmax_);
                    solver.solve(bc, src);
                    break;
                }
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            case AlphaCoeff::ZONI_SHIFTED:
                switch (beta_type_) {
                case BetaCoeff::ZERO: {
                    CartesianR2_ZoniShifted_CircularGeometry src(grid_, Rmax_);
                    solver.solve(bc, src);
                    break;
                }
                case BetaCoeff::ALPHA_INVERSE: {
                    CartesianR2_ZoniShiftedGyro_CircularGeometry src(grid_, Rmax_);
                    solver.solve(bc, src);
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

        /* ------------------------------------------------------------------ */
        case ProblemType::CARTESIAN_R6: {
            CartesianR6_Boundary_CircularGeometry bc(Rmax_);
            switch (alpha_type_) {
            case AlphaCoeff::POISSON: {
                CartesianR6_Poisson_CircularGeometry src(grid_, Rmax_);
                solver.solve(bc, src);
                break;
            }
            case AlphaCoeff::SONNENDRUCKER:
                switch (beta_type_) {
                case BetaCoeff::ZERO: {
                    CartesianR6_Sonnendrucker_CircularGeometry src(grid_, Rmax_);
                    solver.solve(bc, src);
                    break;
                }
                case BetaCoeff::ALPHA_INVERSE: {
                    CartesianR6_SonnendruckerGyro_CircularGeometry src(grid_, Rmax_);
                    solver.solve(bc, src);
                    break;
                }
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            case AlphaCoeff::ZONI:
                switch (beta_type_) {
                case BetaCoeff::ZERO: {
                    CartesianR6_Zoni_CircularGeometry src(grid_, Rmax_);
                    solver.solve(bc, src);
                    break;
                }
                case BetaCoeff::ALPHA_INVERSE: {
                    CartesianR6_ZoniGyro_CircularGeometry src(grid_, Rmax_);
                    solver.solve(bc, src);
                    break;
                }
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            case AlphaCoeff::ZONI_SHIFTED:
                switch (beta_type_) {
                case BetaCoeff::ZERO: {
                    CartesianR6_ZoniShifted_CircularGeometry src(grid_, Rmax_);
                    solver.solve(bc, src);
                    break;
                }
                case BetaCoeff::ALPHA_INVERSE: {
                    CartesianR6_ZoniShiftedGyro_CircularGeometry src(grid_, Rmax_);
                    solver.solve(bc, src);
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

        /* ------------------------------------------------------------------ */
        case ProblemType::POLAR_R6: {
            PolarR6_Boundary_CircularGeometry bc(Rmax_);
            switch (alpha_type_) {
            case AlphaCoeff::POISSON: {
                PolarR6_Poisson_CircularGeometry src(grid_, Rmax_);
                solver.solve(bc, src);
                break;
            }
            case AlphaCoeff::SONNENDRUCKER:
                switch (beta_type_) {
                case BetaCoeff::ZERO: {
                    PolarR6_Sonnendrucker_CircularGeometry src(grid_, Rmax_);
                    solver.solve(bc, src);
                    break;
                }
                case BetaCoeff::ALPHA_INVERSE: {
                    PolarR6_SonnendruckerGyro_CircularGeometry src(grid_, Rmax_);
                    solver.solve(bc, src);
                    break;
                }
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            case AlphaCoeff::ZONI:
                switch (beta_type_) {
                case BetaCoeff::ZERO: {
                    PolarR6_Zoni_CircularGeometry src(grid_, Rmax_);
                    solver.solve(bc, src);
                    break;
                }
                case BetaCoeff::ALPHA_INVERSE: {
                    PolarR6_ZoniGyro_CircularGeometry src(grid_, Rmax_);
                    solver.solve(bc, src);
                    break;
                }
                default:
                    throw std::runtime_error("Invalid beta.\n");
                }
                break;
            case AlphaCoeff::ZONI_SHIFTED:
                switch (beta_type_) {
                case BetaCoeff::ZERO: {
                    PolarR6_ZoniShifted_CircularGeometry src(grid_, Rmax_);
                    solver.solve(bc, src);
                    break;
                }
                case BetaCoeff::ALPHA_INVERSE: {
                    PolarR6_ZoniShiftedGyro_CircularGeometry src(grid_, Rmax_);
                    solver.solve(bc, src);
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

        /* ------------------------------------------------------------------ */
        case ProblemType::REFINED_RADIUS: {
            Refined_Boundary_CircularGeometry bc(Rmax_);
            switch (alpha_type_) {
            case AlphaCoeff::ZONI_SHIFTED:
                switch (beta_type_) {
                case BetaCoeff::ALPHA_INVERSE: {
                    Refined_ZoniShiftedGyro_CircularGeometry src(grid_, Rmax_);
                    solver.solve(bc, src);
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
            throw std::runtime_error("Invalid problem.\n");
        }
    }
}

// Explicit instantiations — CircularGeometry for all coefficient types
template void
ConfigParser::solve<CircularGeometry, PoissonCoefficients>(GMGPolar<CircularGeometry, PoissonCoefficients>&) const;
template void ConfigParser::solve<CircularGeometry, SonnendruckerCoefficients>(
    GMGPolar<CircularGeometry, SonnendruckerCoefficients>&) const;
template void ConfigParser::solve<CircularGeometry, SonnendruckerGyroCoefficients>(
    GMGPolar<CircularGeometry, SonnendruckerGyroCoefficients>&) const;
template void
ConfigParser::solve<CircularGeometry, ZoniCoefficients>(GMGPolar<CircularGeometry, ZoniCoefficients>&) const;
template void
ConfigParser::solve<CircularGeometry, ZoniGyroCoefficients>(GMGPolar<CircularGeometry, ZoniGyroCoefficients>&) const;
template void ConfigParser::solve<CircularGeometry, ZoniShiftedCoefficients>(
    GMGPolar<CircularGeometry, ZoniShiftedCoefficients>&) const;
template void ConfigParser::solve<CircularGeometry, ZoniShiftedGyroCoefficients>(
    GMGPolar<CircularGeometry, ZoniShiftedGyroCoefficients>&) const;
