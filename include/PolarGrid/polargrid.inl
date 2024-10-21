#pragma once

#include "polargrid.h"

inline int PolarGrid::nr() const{
    return nr_;
}
inline int PolarGrid::ntheta() const{
    return ntheta_;
}

inline int PolarGrid::numberOfNodes() const{
    return nr() * ntheta();
}

inline const double& PolarGrid::radius(const int r_index) const { 
    assert(r_index >= 0 && static_cast<size_t>(r_index) < radii_.size());
    return radii_[r_index];
}

inline const double& PolarGrid::theta(const int theta_index) const { 
    assert(theta_index >= 0 && static_cast<size_t>(theta_index) < angles_.size());
    return angles_[theta_index]; 
}

// Get the number of circles in the circular smoother.
inline int PolarGrid::numberSmootherCircles() const{
    return number_smoother_circles_;
}
// Get the length of the radial smoother lines. 
inline int PolarGrid::lengthSmootherRadial() const{
    return length_smoother_radial_;
}

// Get the number of nodes in circular smoother.
inline int PolarGrid::numberCircularSmootherNodes() const{
    return number_circular_smoother_nodes_;
}
// Get the number of nodes in radial smoother.
inline int PolarGrid::numberRadialSmootherNodes() const{
    return number_radial_smoother_nodes_;
}

inline const double& PolarGrid::radialSpacing(const int r_index) const{ 
    assert(r_index >= 0 && static_cast<size_t>(r_index) < radial_spacings_.size());
    return radial_spacings_[r_index];
}

inline const double& PolarGrid::angularSpacing(const int unwrapped_theta_index) const{ 
    // unwrapped_theta_index may be negative or larger than ntheta() to allow for periodicity.
    const int theta_index = wrapThetaIndex(unwrapped_theta_index);
    assert(theta_index >= 0 && theta_index < ntheta());
    return angular_spacings_[theta_index];
}

// ------------------ //
// Optimized indexing //
// ------------------ //

inline int PolarGrid::wrapThetaIndex(const int unwrapped_theta_index) const
{   
    // unwrapped_theta_index may be negative or larger than ntheta() to allow for periodicity.
    // If ntheta = 2^k we can use the optimization idx & (ntheta-1)
    return is_ntheta_PowerOfTwo_ ? 
        unwrapped_theta_index & (ntheta() - 1) : 
        (unwrapped_theta_index % ntheta() + ntheta()) % ntheta();
}

inline int PolarGrid::index(const int r_index, const int unwrapped_theta_index) const
{
    // unwrapped_theta_index may be negative or larger than ntheta() to allow for periodicity.
    assert(0 <= r_index && r_index < nr());
    const int theta_index = wrapThetaIndex(unwrapped_theta_index);
    assert(0 <= theta_index && theta_index < ntheta());
    return r_index < numberSmootherCircles() ? 
        theta_index + ntheta() * r_index : 
        numberCircularSmootherNodes() + 
        r_index - numberSmootherCircles() + lengthSmootherRadial() * theta_index;
}

inline int PolarGrid::fastIndex(const int r_index, const int theta_index) const
{
    // unwrapped_theta_index may be negative or larger than ntheta() to allow for periodicity.
    assert(0 <= r_index && r_index < nr());
    assert(0 <= theta_index && theta_index < ntheta());
    return r_index < numberSmootherCircles() ? 
        theta_index + ntheta() * r_index : 
        numberCircularSmootherNodes() + 
        r_index - numberSmootherCircles() + lengthSmootherRadial() * theta_index;
}

inline void PolarGrid::multiIndex(const int node_index, int& r_index, int& theta_index) const
{
    assert(0 <= node_index && node_index < numberOfNodes());
    if(node_index < numberCircularSmootherNodes()){
        r_index = node_index / ntheta();
        theta_index = is_ntheta_PowerOfTwo_ ? node_index & (ntheta() - 1) : node_index % ntheta();
    }
    else{
        theta_index = (node_index-numberCircularSmootherNodes()) / lengthSmootherRadial();
        r_index = numberSmootherCircles() + (node_index-numberCircularSmootherNodes()) % lengthSmootherRadial();
    }
}
