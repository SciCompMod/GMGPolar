#include <ConfigParser/config_parser.h>
#include <GMGPolar/gmgpolar.h>

using namespace gmgpolar;

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
void ConfigParser::solve(GMGPolar<DomainGeometry, DensityProfileCoefficients>& solver) const
{
    if constexpr (std::is_same_v<DomainGeometry, CulhamGeometry>) {

        switch (problem_type_) {

            /* ------------------------------------------------------------------ */
            /* Culham only supports POLAR_R6 and REFINED_RADIUS.                  */
            /* CARTESIAN_R2 and CARTESIAN_R6 are not defined for this geometry.   */
            /* ------------------------------------------------------------------ */

        case ProblemType::POLAR_R6: {
            PolarR6_Boundary_CulhamGeometry bc(Rmax_);
            switch (alpha_type_) {
            case AlphaCoeff::ZONI_SHIFTED:
                switch (beta_type_) {
                case BetaCoeff::ALPHA_INVERSE: {
                    PolarR6_ZoniShiftedGyro_CulhamGeometry src(grid_, Rmax_);
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

        /* ------------------------------------------------------------------ */
        case ProblemType::REFINED_RADIUS: {
            Refined_Boundary_CulhamGeometry bc(Rmax_);
            switch (alpha_type_) {
            case AlphaCoeff::ZONI_SHIFTED:
                switch (beta_type_) {
                case BetaCoeff::ALPHA_INVERSE: {
                    Refined_ZoniShiftedGyro_CulhamGeometry src(grid_, Rmax_);
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
            throw std::runtime_error("Invalid geometry for configuration.\n");
        }
    }
}

// Explicit instantiations — CulhamGeometry for all coefficient types
template void
ConfigParser::solve<CulhamGeometry, PoissonCoefficients>(GMGPolar<CulhamGeometry, PoissonCoefficients>&) const;
template void ConfigParser::solve<CulhamGeometry, SonnendruckerCoefficients>(
    GMGPolar<CulhamGeometry, SonnendruckerCoefficients>&) const;
template void ConfigParser::solve<CulhamGeometry, SonnendruckerGyroCoefficients>(
    GMGPolar<CulhamGeometry, SonnendruckerGyroCoefficients>&) const;
template void ConfigParser::solve<CulhamGeometry, ZoniCoefficients>(GMGPolar<CulhamGeometry, ZoniCoefficients>&) const;
template void
ConfigParser::solve<CulhamGeometry, ZoniGyroCoefficients>(GMGPolar<CulhamGeometry, ZoniGyroCoefficients>&) const;
template void
ConfigParser::solve<CulhamGeometry, ZoniShiftedCoefficients>(GMGPolar<CulhamGeometry, ZoniShiftedCoefficients>&) const;
template void ConfigParser::solve<CulhamGeometry, ZoniShiftedGyroCoefficients>(
    GMGPolar<CulhamGeometry, ZoniShiftedGyroCoefficients>&) const;
