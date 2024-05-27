inline int PolarGrid::nr() const{
    return nr_;
}
inline int PolarGrid::ntheta() const{
    return ntheta_;
}

inline int PolarGrid::number_of_nodes() const{
    return nr() * ntheta();
}

inline int PolarGrid::number_of_inner_nodes() const{
    if (nr() < 2 || ntheta() < 1 ) return 0;
    return (nr() - 2) * ntheta();
}

inline int PolarGrid::number_of_inner_boundary_nodes() const{
    return ntheta();
}

inline int PolarGrid::number_of_outer_boundary_nodes() const{
    return ntheta();
}

inline const scalar_t& PolarGrid::radius(const int r_index) const { 
    assert(r_index >= 0 && static_cast<size_t>(r_index) < radii_.size());
    return radii_[r_index];
}

inline const scalar_t& PolarGrid::theta(const int theta_index) const { 
    assert(theta_index >= 0 && static_cast<size_t>(theta_index) < angles_.size());
    return angles_[theta_index]; 
}

// Get the number of circles in the circular smoother.
inline int PolarGrid::numberSmootherCircles() const{
    return numberSmootherCircles_;
}
// Get the length of the radial smoother lines. 
inline int PolarGrid::lengthSmootherRadial() const{
    return lengthSmootherRadial_;
}

// Get the number of nodes in circular smoother.
inline int PolarGrid::numberCircularSmootherNodes() const{
    return numberCircularSmootherNodes_;
}
// Get the number of nodes in radial smoother.
inline int PolarGrid::numberRadialSmootherNodes() const{
    return numberRadialSmootherNodes_;
}


inline int PolarGrid::numberBlackSmootherCircles() const{
    return numberBlackSmootherCircles_;
}
inline int PolarGrid::numberWhiteSmootherCircles() const{
    return numberWhiteSmootherCircles_;
}
inline int PolarGrid::numberBlackSmootherRadials() const{
    return numberBlackSmootherRadials_;
}
inline int PolarGrid::numberWhiteSmootherRadials() const{
    return numberWhiteSmootherRadials_;
}

inline const scalar_t& PolarGrid::r_dist(const int r_index) const{ 
    assert(r_index >= 0 && static_cast<size_t>(r_index) < r_dist_.size());
    return r_dist_[r_index];
}

inline const scalar_t& PolarGrid::theta_dist(const int unwrapped_theta_index) const{ 
    // unwrapped_theta_index may be negative or larger than ntheta() to allow for periodicity.
    const int theta_index = is_ntheta_PowerOfTwo_ ? unwrapped_theta_index & (ntheta() - 1) : (unwrapped_theta_index % ntheta() + ntheta()) % ntheta();
    return theta_dist_[theta_index];
}

// ------------------ //
// Optimized indexing //
// ------------------ //

inline int PolarGrid::index(const int r_index, const int unwrapped_theta_index) const
{
    // unwrapped_theta_index may be negative or larger than ntheta() to allow for periodicity.
    // If ntheta = 2^k we can use the optimization idx & (ntheta-1)
    assert(r_index >= 0 && r_index < nr());
    const int theta_index = is_ntheta_PowerOfTwo_ ? unwrapped_theta_index & (ntheta() - 1) : (unwrapped_theta_index % ntheta() + ntheta()) % ntheta();
    return r_index < numberSmootherCircles() ? 
        theta_index + ntheta() * r_index : 
        numberCircularSmootherNodes() + 
        r_index - numberSmootherCircles() + lengthSmootherRadial() * theta_index;
}


// inline int PolarGrid::index(const int r_index, const int unwrapped_theta_index) const
// {
//     // unwrapped_theta_index may be negative or larger than ntheta() to allow for periodicity.
//     // If ntheta = 2^k we can use the optimization idx & (ntheta-1)
//     assert(r_index >= 0 && r_index < nr());
//     const int theta_index = is_ntheta_PowerOfTwo_ ? unwrapped_theta_index & (ntheta() - 1) : (unwrapped_theta_index % ntheta() + ntheta()) % ntheta();

//     return 
//     r_index < numberSmootherCircles() ?
//         r_index & 1 ?
//             /* White Circle */
//             (r_index >> 1) * ntheta() + theta_index + num_BC_
//             :
//             /* Black Circle */
//             (r_index >> 1) * ntheta() + theta_index
//     :
//         theta_index & 1 ?
//             /* White Radial */
//             numberCircularSmootherNodes() + (r_index-numberSmootherCircles() + lengthSmootherRadial() * (theta_index>>1)) + num_BR_
//             :
//             /* Black Radial */
//             numberCircularSmootherNodes() + (r_index-numberSmootherCircles() + lengthSmootherRadial() * (theta_index>>1));
// }



// inline int PolarGrid::index(const int r_index, const int unwrapped_theta_index) const
// {
//     // unwrapped_theta_index may be negative or larger than ntheta() to allow for periodicity.
//     // If ntheta = 2^k we can use the optimization idx & (ntheta-1)
//     assert(r_index >= 0 && r_index < nr());

//     int theta_index;
//     if(is_ntheta_PowerOfTwo_){
//         theta_index = unwrapped_theta_index & (ntheta() - 1);
//     } else{
//         theta_index = (unwrapped_theta_index % ntheta() + ntheta()) % ntheta();;
//     }

//     if (r_index < numberSmootherCircles()) {
//         if (r_index & 1) {
//             // White Circle
//             return (r_index >> 1) * ntheta() + theta_index + num_BC_;
//         } else {
//             // Black Circle
//             return (r_index >> 1) * ntheta() + theta_index;
//         }
//     } else {
//         if (theta_index & 1) {
//             // White Radial
//             return numberCircularSmootherNodes() + (r_index - numberSmootherCircles() + lengthSmootherRadial() * (theta_index >> 1)) + num_BR_;
//         } else {
//             // Black Radial
//             return numberCircularSmootherNodes() + (r_index - numberSmootherCircles() + lengthSmootherRadial() * (theta_index >> 1));
//         }
//     }
// }




// inline int PolarGrid::index(const int r_index, const int unwrapped_theta_index) const
// {
//     // unwrapped_theta_index may be negative or larger than ntheta() to allow for periodicity.
//     // If ntheta = 2^k we can use the optimization idx & (ntheta-1)
//     assert(r_index >= 0 && r_index < nr());

//     int theta_index;
//     if(is_ntheta_PowerOfTwo_){
//         theta_index = unwrapped_theta_index & (ntheta() - 1);
//     } else{
//         theta_index = (unwrapped_theta_index % ntheta() + ntheta()) % ntheta();;
//     }

//     if (r_index < numberSmootherCircles()) {

//         return (r_index & 1) ? 
//             // White Circle
//             (r_index >> 1) * ntheta() + theta_index + num_BC_
//             :
//             // Black Circle
//             (r_index >> 1) * ntheta() + theta_index;

//     } else {
//         return (theta_index & 1) ? 
//             // White Radial
//             numberCircularSmootherNodes() + (r_index - numberSmootherCircles() + lengthSmootherRadial() * (theta_index >> 1)) + num_BR_
//         :
//             // Black Radial
//             numberCircularSmootherNodes() + (r_index - numberSmootherCircles() + lengthSmootherRadial() * (theta_index >> 1));
//     }
// }
