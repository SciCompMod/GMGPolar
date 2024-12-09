#pragma once

#include "polargrid.h"

__host__ __device__ __forceinline__ int PolarGrid::wrapThetaIndex(const int unwrapped_theta_index) const
{
    // unwrapped_theta_index may be negative or larger than ntheta() to allow for periodicity.
    // If ntheta = 2^k, we can use the optimization idx & (ntheta - 1).
    if (is_ntheta_PowerOfTwo_) {
        return unwrapped_theta_index & (ntheta() - 1);
    }
    else {
        int wrapped_index = unwrapped_theta_index % ntheta();
        // Adjust for negative indices
        if (wrapped_index < 0) {
            wrapped_index += ntheta();
        }
        return wrapped_index;
    }
}

__host__ __device__ __forceinline__ int PolarGrid::index(const int r_index, const int unwrapped_theta_index) const
{
    assert(0 <= r_index && r_index < nr());
    const int theta_index = wrapThetaIndex(unwrapped_theta_index);
    assert(0 <= theta_index && theta_index < ntheta());
    return r_index < numberSmootherCircles() ? theta_index + ntheta() * r_index
                                             : numberCircularSmootherNodes() + r_index - numberSmootherCircles() +
                                                   lengthSmootherRadial() * theta_index;
}
__host__ __device__ __forceinline__ void PolarGrid::multiIndex(const int node_index, int& r_index,
                                                               int& theta_index) const
{
    assert(0 <= node_index && node_index < numberOfNodes());
    if (node_index < numberCircularSmootherNodes()) {
        r_index     = node_index / ntheta();
        theta_index = is_ntheta_PowerOfTwo_ ? node_index & (ntheta() - 1) : node_index % ntheta();
    }
    else {
        theta_index = (node_index - numberCircularSmootherNodes()) / lengthSmootherRadial();
        r_index     = numberSmootherCircles() + (node_index - numberCircularSmootherNodes()) % lengthSmootherRadial();
    }
}

__host__ __device__ __forceinline__ int PolarGrid::numberOfNodes() const
{
    return nr() * ntheta();
}
__host__ __device__ __forceinline__ int PolarGrid::nr() const
{
    return nr_;
}
__host__ __device__ __forceinline__ int PolarGrid::ntheta() const
{
    return ntheta_;
}

__host__ __device__ __forceinline__ double PolarGrid::radius(const int r_index) const
{
    assert(r_index >= 0 && r_index < nr());
#ifdef __CUDA_ARCH__
    return d_radii_[r_index];
#else
    return radii_[r_index];
#endif
}
__host__ __device__ __forceinline__ double PolarGrid::theta(const int theta_index) const
{
    assert(theta_index >= 0 && theta_index < ntheta() + 1);
#ifdef __CUDA_ARCH__
    return d_angles_[theta_index];
#else
    return angles_[theta_index];
#endif
}

__host__ __device__ __forceinline__ double PolarGrid::radialSpacing(const int r_index) const
{
    assert(r_index >= 0 && r_index < nr() - 1);
#ifdef __CUDA_ARCH__
    return d_radial_spacings_[r_index];
#else
    return radial_spacings_[r_index];
#endif
}
__host__ __device__ __forceinline__ double PolarGrid::angularSpacing(const int unwrapped_theta_index) const
{
    // unwrapped_theta_index may be negative or larger than ntheta() to allow for periodicity.
    const int theta_index = wrapThetaIndex(unwrapped_theta_index);
    assert(theta_index >= 0 && theta_index < ntheta());
#ifdef __CUDA_ARCH__
    return d_angular_spacings_[theta_index];
#else
    return angular_spacings_[theta_index];
#endif
}

__host__ __device__ __forceinline__ double PolarGrid::smootherSplittingRadius() const
{
    return smoother_splitting_radius_;
}

__host__ __device__ __forceinline__ int PolarGrid::numberSmootherCircles() const
{
    return number_smoother_circles_;
}
__host__ __device__ __forceinline__ int PolarGrid::lengthSmootherRadial() const
{
    return length_smoother_radial_;
}
__host__ __device__ __forceinline__ int PolarGrid::numberCircularSmootherNodes() const
{
    return number_circular_smoother_nodes_;
}
__host__ __device__ __forceinline__ int PolarGrid::numberRadialSmootherNodes() const
{
    return number_radial_smoother_nodes_;
}
