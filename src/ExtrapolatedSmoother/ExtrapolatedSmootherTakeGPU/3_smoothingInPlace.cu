#include "../../../include/ExtrapolatedSmoother/ExtrapolatedSmootherTakeGPU/extrapolatedSmoother.h"

void ExtrapolatedSmootherTakeGPU::extrapolatedSmoothingInPlace(GPU_Vector<double>& x, const GPU_Vector<double>& rhs, GPU_Vector<double>& temp)
{
    const PolarGrid& grid = level_.grid();

    assert(x.size() == grid.numberOfNodes());
    assert(rhs.size() == grid.numberOfNodes());
    assert(temp.size() == grid.numberOfNodes());

    DomainGeometry* device_domain_geometry;
    cudaMalloc(&device_domain_geometry, sizeof(DomainGeometry));
    cudaMemcpy(device_domain_geometry, &domain_geometry_, sizeof(DomainGeometry), cudaMemcpyHostToDevice);

    /* We use precomputed DensityProfileCoefficients values. */
    // DensityProfileCoefficients* device_density_profile;
    // cudaMalloc(&device_density_profile, sizeof(DensityProfileCoefficients));
    // cudaMemcpy(device_density_profile, &density_profile_coefficients_, sizeof(DensityProfileCoefficients), cudaMemcpyHostToDevice);
    
    applyAscOrtho_BlackCircle(x, rhs, temp, device_domain_geometry);
    solveAsc_BlackCircle(x, rhs, temp); 

    applyAscOrtho_WhiteCircle(x, rhs, temp, device_domain_geometry);
    solveAsc_WhiteCircle(x, rhs, temp);

    applyAscOrtho_BlackRadial(x, rhs, device_domain_geometry);
    solveAsc_BlackRadial(x, rhs);

    applyAscOrtho_WhiteRadial(x, rhs, device_domain_geometry);
    solveAsc_WhiteRadial(x, rhs);

    /* We use precomputed DensityProfileCoefficients values. */
    cudaFree(device_domain_geometry);
    // cudaFree(device_density_profile);
}