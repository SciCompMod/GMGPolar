#pragma once

#include <iterator>

#include "polargrid.h"

// Constructor to initialize grid using vectors of radii and angles.
template<class MemorySpace>
template<class MemorySpace2>
PolarGrid<MemorySpace>::PolarGrid(Vector<double, MemorySpace2> radii, Vector<double, MemorySpace2> angles, std::optional<double> splitting_radius)
    : nr_(radii.size())
    , ntheta_(angles.size() - 1)
    , is_ntheta_PowerOfTwo_((ntheta_ & (ntheta_ - 1)) == 0)
    , radii_(radii)
    , angles_(angles)
{
    // Check parameter validity
    checkParameters(radii, angles);
    // Store distances to adjacent neighboring nodes.
    // Initializes radial_spacings_, angular_spacings_
    initializeDistances();
    // Initializes smoothers splitting radius for circle/radial indexing.
    initializeLineSplitting(splitting_radius);
}

template<class MemorySpace>
KOKKOS_INLINE_FUNCTION int PolarGrid<MemorySpace>::nr() const
{
    return nr_;
}
template<class MemorySpace>
KOKKOS_INLINE_FUNCTION int PolarGrid<MemorySpace>::ntheta() const
{
    return ntheta_;
}

template<class MemorySpace>
KOKKOS_INLINE_FUNCTION int PolarGrid<MemorySpace>::numberOfNodes() const
{
    return nr() * ntheta();
}

template<class MemorySpace>
KOKKOS_INLINE_FUNCTION double PolarGrid<MemorySpace>::radius(const int r_index) const
{
    assert(r_index >= 0 && r_index < std::ssize(radii_));
    return radii_[r_index];
}

template<class MemorySpace>
KOKKOS_INLINE_FUNCTION double PolarGrid<MemorySpace>::theta(const int theta_index) const
{
    assert(theta_index >= 0 && theta_index < std::ssize(angles_));
    return angles_[theta_index];
}

template<class MemorySpace>
KOKKOS_INLINE_FUNCTION double PolarGrid<MemorySpace>::smootherSplittingRadius() const
{
    return smoother_splitting_radius_;
}

// Get the number of circles in the circular smoother.
template<class MemorySpace>
KOKKOS_INLINE_FUNCTION int PolarGrid<MemorySpace>::numberSmootherCircles() const
{
    return number_smoother_circles_;
}
// Get the length of the radial smoother lines.
template<class MemorySpace>
KOKKOS_INLINE_FUNCTION int PolarGrid<MemorySpace>::lengthRadialSmoother() const
{
    return length_smoother_radial_;
}

// Get the number of nodes in circular smoother.
template<class MemorySpace>
KOKKOS_INLINE_FUNCTION int PolarGrid<MemorySpace>::numberCircularSmootherNodes() const
{
    return number_circular_smoother_nodes_;
}
// Get the number of nodes in radial smoother.
template<class MemorySpace>
KOKKOS_INLINE_FUNCTION int PolarGrid<MemorySpace>::numberRadialSmootherNodes() const
{
    return number_radial_smoother_nodes_;
}

template<class MemorySpace>
KOKKOS_INLINE_FUNCTION double PolarGrid<MemorySpace>::radialSpacing(const int r_index) const
{
    assert(r_index >= 0 && r_index < std::ssize(radial_spacings_));
    return radial_spacings_[r_index];
}

template<class MemorySpace>
KOKKOS_INLINE_FUNCTION double PolarGrid<MemorySpace>::angularSpacing(const int unwrapped_theta_index) const
{
    // unwrapped_theta_index may be negative or larger than ntheta() to allow for periodicity.
    const int theta_index = wrapThetaIndex(unwrapped_theta_index);
    return angular_spacings_[theta_index];
}

template<class MemorySpace>
KOKKOS_INLINE_FUNCTION int PolarGrid<MemorySpace>::wrapThetaIndex(const int unwrapped_theta_index) const
{
    // The unwrapped_theta_index may be negative or exceed the number of theta steps (ntheta()),
    // so we need to wrap it into the valid range [0, ntheta() - 1] to maintain periodicity.
    //
    // When ntheta is a power of 2 (i.e., ntheta = 2^k), we can optimize the modulo operation:
    // - ntheta() - 1 yields a bitmask with the lower k bits set to 1.
    //   For example, if ntheta() is 8 (2^3), then ntheta() - 1 equals 7, which in binary is 0b111.
    // - Applying the bitwise AND operator (&) with this mask extracts the lower k bits of unwrapped_theta_index.
    //   This effectively computes unwrapped_theta_index % ntheta(), because it discards all higher bits.
    //
    // If ntheta is not a power of two, we use the standard modulo approach to handle wrapping.
    int theta_index = is_ntheta_PowerOfTwo_ ? unwrapped_theta_index & (ntheta() - 1)
                                            : (unwrapped_theta_index % ntheta() + ntheta()) % ntheta();
    assert(0 <= theta_index && theta_index < ntheta());
    return theta_index;
}

template<class MemorySpace>
KOKKOS_INLINE_FUNCTION int PolarGrid<MemorySpace>::index(const int r_index, const int unwrapped_theta_index) const
{
    // unwrapped_theta_index may be negative or larger than ntheta() to allow for periodicity.
    assert(0 <= r_index && r_index < nr());
    int theta_index = wrapThetaIndex(unwrapped_theta_index);
    int global_index =
        r_index < numberSmootherCircles()
            ? theta_index + ntheta() * r_index
            : numberCircularSmootherNodes() + r_index - numberSmootherCircles() + lengthRadialSmoother() * theta_index;
    assert(0 <= global_index && global_index < numberOfNodes());
    return global_index;
}

template<class MemorySpace>
KOKKOS_INLINE_FUNCTION void PolarGrid<MemorySpace>::multiIndex(const int node_index, int& r_index, int& theta_index) const
{
    assert(0 <= node_index && node_index < numberOfNodes());
    if (node_index < numberCircularSmootherNodes()) {
        r_index     = node_index / ntheta();
        theta_index = wrapThetaIndex(node_index);
    }
    else {
        theta_index = (node_index - numberCircularSmootherNodes()) / lengthRadialSmoother();
        r_index     = numberSmootherCircles() + (node_index - numberCircularSmootherNodes()) % lengthRadialSmoother();
    }
}
