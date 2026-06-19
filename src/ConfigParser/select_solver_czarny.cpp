#include <ConfigParser/config_parser.h>
#include <GMGPolar/gmgpolar.h>

using namespace gmgpolar;

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
void ConfigParser::solve(GMGPolar<DomainGeometry, DensityProfileCoefficients>& solver) const
{
    if constexpr (std::is_same_v<DomainGeometry, CzarnyGeometry>) {

        switch (problem_type_) {

        /* ------------------------------------------------------------------ */
        case ProblemType::CARTESIAN_R2: {
            CartesianR2_Boundary_CzarnyGeometry bc(Rmax_, kappa_eps_, delta_e_);
            switch (alpha_type_) {
            case AlphaCoeff::POISSON: {
                CartesianR2_Poisson_CzarnyGeometry src(grid_, Rmax_, kappa_eps_, delta_e_);
                solver.solve(bc, src);
                break;
            }
            case AlphaCoeff::SONNENDRUCKER:
                switch (beta_type_) {
                case BetaCoeff::ZERO: {
                    CartesianR2_Sonnendrucker_CzarnyGeometry src(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(bc, src);
                    break;
                }
                case BetaCoeff::ALPHA_INVERSE: {
                    CartesianR2_SonnendruckerGyro_CzarnyGeometry src(grid_, Rmax_, kappa_eps_, delta_e_);
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
                    CartesianR2_Zoni_CzarnyGeometry src(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(bc, src);
                    break;
                }
                case BetaCoeff::ALPHA_INVERSE: {
                    CartesianR2_ZoniGyro_CzarnyGeometry src(grid_, Rmax_, kappa_eps_, delta_e_);
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
                    CartesianR2_ZoniShifted_CzarnyGeometry src(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(bc, src);
                    break;
                }
                case BetaCoeff::ALPHA_INVERSE: {
                    CartesianR2_ZoniShiftedGyro_CzarnyGeometry src(grid_, Rmax_, kappa_eps_, delta_e_);
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
            CartesianR6_Boundary_CzarnyGeometry bc(Rmax_, kappa_eps_, delta_e_);
            switch (alpha_type_) {
            case AlphaCoeff::POISSON: {
                CartesianR6_Poisson_CzarnyGeometry src(grid_, Rmax_, kappa_eps_, delta_e_);
                solver.solve(bc, src);
                break;
            }
            case AlphaCoeff::SONNENDRUCKER:
                switch (beta_type_) {
                case BetaCoeff::ZERO: {
                    CartesianR6_Sonnendrucker_CzarnyGeometry src(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(bc, src);
                    break;
                }
                case BetaCoeff::ALPHA_INVERSE: {
                    CartesianR6_SonnendruckerGyro_CzarnyGeometry src(grid_, Rmax_, kappa_eps_, delta_e_);
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
                    CartesianR6_Zoni_CzarnyGeometry src(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(bc, src);
                    break;
                }
                case BetaCoeff::ALPHA_INVERSE: {
                    CartesianR6_ZoniGyro_CzarnyGeometry src(grid_, Rmax_, kappa_eps_, delta_e_);
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
                    CartesianR6_ZoniShifted_CzarnyGeometry src(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(bc, src);
                    break;
                }
                case BetaCoeff::ALPHA_INVERSE: {
                    CartesianR6_ZoniShiftedGyro_CzarnyGeometry src(grid_, Rmax_, kappa_eps_, delta_e_);
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
            PolarR6_Boundary_CzarnyGeometry bc(Rmax_, kappa_eps_, delta_e_);
            switch (alpha_type_) {
            case AlphaCoeff::POISSON: {
                PolarR6_Poisson_CzarnyGeometry src(grid_, Rmax_, kappa_eps_, delta_e_);
                solver.solve(bc, src);
                break;
            }
            case AlphaCoeff::SONNENDRUCKER:
                switch (beta_type_) {
                case BetaCoeff::ZERO: {
                    PolarR6_Sonnendrucker_CzarnyGeometry src(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(bc, src);
                    break;
                }
                case BetaCoeff::ALPHA_INVERSE: {
                    PolarR6_SonnendruckerGyro_CzarnyGeometry src(grid_, Rmax_, kappa_eps_, delta_e_);
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
                    PolarR6_Zoni_CzarnyGeometry src(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(bc, src);
                    break;
                }
                case BetaCoeff::ALPHA_INVERSE: {
                    PolarR6_ZoniGyro_CzarnyGeometry src(grid_, Rmax_, kappa_eps_, delta_e_);
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
                    PolarR6_ZoniShifted_CzarnyGeometry src(grid_, Rmax_, kappa_eps_, delta_e_);
                    solver.solve(bc, src);
                    break;
                }
                case BetaCoeff::ALPHA_INVERSE: {
                    PolarR6_ZoniShiftedGyro_CzarnyGeometry src(grid_, Rmax_, kappa_eps_, delta_e_);
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
            Refined_Boundary_CzarnyGeometry bc(Rmax_, kappa_eps_, delta_e_);
            switch (alpha_type_) {
            case AlphaCoeff::ZONI_SHIFTED:
                switch (beta_type_) {
                case BetaCoeff::ALPHA_INVERSE: {
                    Refined_ZoniShiftedGyro_CzarnyGeometry src(grid_, Rmax_, kappa_eps_, delta_e_);
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

// Explicit instantiations — CzarnyGeometry for all coefficient types
template void
ConfigParser::solve<CzarnyGeometry, PoissonCoefficients>(GMGPolar<CzarnyGeometry, PoissonCoefficients>&) const;
template void ConfigParser::solve<CzarnyGeometry, SonnendruckerCoefficients>(
    GMGPolar<CzarnyGeometry, SonnendruckerCoefficients>&) const;
template void ConfigParser::solve<CzarnyGeometry, SonnendruckerGyroCoefficients>(
    GMGPolar<CzarnyGeometry, SonnendruckerGyroCoefficients>&) const;
template void ConfigParser::solve<CzarnyGeometry, ZoniCoefficients>(GMGPolar<CzarnyGeometry, ZoniCoefficients>&) const;
template void
ConfigParser::solve<CzarnyGeometry, ZoniGyroCoefficients>(GMGPolar<CzarnyGeometry, ZoniGyroCoefficients>&) const;
template void
ConfigParser::solve<CzarnyGeometry, ZoniShiftedCoefficients>(GMGPolar<CzarnyGeometry, ZoniShiftedCoefficients>&) const;
template void ConfigParser::solve<CzarnyGeometry, ZoniShiftedGyroCoefficients>(
    GMGPolar<CzarnyGeometry, ZoniShiftedGyroCoefficients>&) const;
