#pragma once

#include "polargrid.h"

inline int PolarGrid::nr() const
{
    return nr_;
}
inline int PolarGrid::ntheta() const
{
    return ntheta_;
}

inline int PolarGrid::numberOfNodes() const
{
    return nr() * ntheta();
}

inline double PolarGrid::radius(const int r_index) const
{
    assert(r_index >= 0 && static_cast<size_t>(r_index) < radii_.size());
    return radii_[r_index];
}

inline double PolarGrid::theta(const int theta_index) const
{
    assert(theta_index >= 0 && static_cast<size_t>(theta_index) < angles_.size());
    return angles_[theta_index];
}

inline double PolarGrid::smootherSplittingRadius() const
{
    return smoother_splitting_radius_;
}

// Get the number of circles in the circular smoother.
inline int PolarGrid::numberSmootherCircles() const
{
    return number_smoother_circles_;
}
// Get the length of the radial smoother lines.
inline int PolarGrid::lengthSmootherRadial() const
{
    return length_smoother_radial_;
}

// Get the number of nodes in circular smoother.
inline int PolarGrid::numberCircularSmootherNodes() const
{
    return number_circular_smoother_nodes_;
}
// Get the number of nodes in radial smoother.
inline int PolarGrid::numberRadialSmootherNodes() const
{
    return number_radial_smoother_nodes_;
}

inline double PolarGrid::radialSpacing(const int r_index) const
{
    assert(r_index >= 0 && static_cast<size_t>(r_index) < radial_spacings_.size());
    return radial_spacings_[r_index];
}

inline double PolarGrid::angularSpacing(const int unwrapped_theta_index) const
{
    // unwrapped_theta_index may be negative or larger than ntheta() to allow for periodicity.
    const int theta_index = wrapThetaIndex(unwrapped_theta_index);
    return angular_spacings_[theta_index];
}

inline int PolarGrid::wrapThetaIndex(const int unwrapped_theta_index) const
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

inline int PolarGrid::index(const int r_index, const int unwrapped_theta_index) const
{
    // unwrapped_theta_index may be negative or larger than ntheta() to allow for periodicity.
    assert(0 <= r_index && r_index < nr());
    int theta_index = wrapThetaIndex(unwrapped_theta_index);
    int global_index =
        r_index < numberSmootherCircles()
            ? theta_index + ntheta() * r_index
            : numberCircularSmootherNodes() + r_index - numberSmootherCircles() + lengthSmootherRadial() * theta_index;
    assert(0 <= global_index && global_index < numberOfNodes());
    return global_index;
}

inline void PolarGrid::multiIndex(const int node_index, int& r_index, int& theta_index) const
{
    assert(0 <= node_index && node_index < numberOfNodes());
    if (node_index < numberCircularSmootherNodes()) {
        r_index     = node_index / ntheta();
        theta_index = wrapThetaIndex(node_index);
    }
    else {
        theta_index = (node_index - numberCircularSmootherNodes()) / lengthSmootherRadial();
        r_index     = numberSmootherCircles() + (node_index - numberCircularSmootherNodes()) % lengthSmootherRadial();
    }
}
