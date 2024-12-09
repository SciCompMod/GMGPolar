#include "../../include/Residual/ResidualTakeGPU/residual.h"

ResidualTakeGPU::ResidualTakeGPU(
    const Level& level,
    const DomainGeometry& domain_geometry, const DensityProfileCoefficients& density_profile_coefficients,
    bool DirBC_Interior
) :
    level_(level),
    domain_geometry_(domain_geometry),
    density_profile_coefficients_(density_profile_coefficients),
    DirBC_Interior_(DirBC_Interior)
{}