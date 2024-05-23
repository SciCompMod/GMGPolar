#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <cassert>
#include <set>
#include <functional>
#include <iostream>
#include <iomanip>
#include <memory>
#include <optional>
#include <algorithm>

#include <omp.h>

#include "../common/equals.h"
#include "../common/scalar.h"
#include "../common/space_dimension.h"
#include "../common/constants.h"

#include "../PolarGrid/multiindex.h"
#include "../PolarGrid/point.h"

// The PolarGrid class implements a donut-shaped 2D grid.
// It is designed to handle polar coordinates, providing functionalities
// for storing data points and performing operations on them.

class PolarGrid {
public:
    // Default constructor.
    PolarGrid() = default;

    // Constructor to initialize grid using vectors of radii and angles.
    PolarGrid(const std::vector<scalar_t>& radii, const std::vector<scalar_t>& angles, 
        std::optional<scalar_t> splitting_radius = std::nullopt
    );

    // Constructor to initialize grid using data from text files containing radii and angles.
    PolarGrid(const std::string& file_grid_radii, const std::string& file_grid_angles, 
        std::optional<scalar_t> splitting_radius = std::nullopt
    );

    // Constructor to initialize grid using parameters from GMGPolar.
    // Parameters:
    //   R0: Inner radius of the grid.
    //   R: Outer radius of the grid.
    //   nr_exp: Exponential scaling factor for radial divisions.
    //   nr_exp: Exponential scaling factor for angular divisions.
    //   anisotropic_factor: Factor for adjusting anisotropy in radial divisions.
    //   alpha: Allow refinement of the grid at r_jump,
    //          the center point of the drop of the diffusion coefficient alpha.
    //   divideBy2: Number of times to divide both radial and angular divisions by 2.
    PolarGrid(
        const scalar_t R0, const scalar_t R, const int nr_exp, const int ntheta_exp, const double r_jump,
        const int anisotropic_factor, const alpha_coeff alpha, const int divideBy2, 
        std::optional<scalar_t> splitting_radius = std::nullopt
    );

    // Node Indexing (based on the circle/radial smoother)
    // Get the index of a node based on its position.
    int index(const MultiIndex& position) const;
    // Get the position of a node based on its index.
    MultiIndex multiindex(const int node_index) const;
    // Get the polar coordinates of a node based on its position.
    Point polar_coordinates(const MultiIndex& position) const;

    // Optimized, inlined version of index(). Allows for unwrapped theta indices.
    int index(const int r_index, const int unwrapped_theta_index) const;

    // Neighboring nodes
    // Get adjacent neighbors of a node.
    // If the neighbor index is -1, then there is no neighboring node in that direction.
    void adjacent_neighbors_of(const MultiIndex& position, 
        std::array<std::pair<int,int>, space_dimension>& neighbors
    ) const;
    // Get diagonal neighbors of a node.
    // If the neighbor index is -1, then there is no neighboring node in that direction.
    void diagonal_neighbors_of(
        const MultiIndex& position, 
        std::array<std::pair<int,int>, space_dimension>& neighbors
    ) const;
    // Neighbor distances
    // Get distances to adjacent neighbors of a node.
    // If there is no neighboring node in that direction, then the neighbor distance is 0.
    void adjacent_neighbor_distances(const MultiIndex& position, 
        std::array<std::pair<scalar_t,scalar_t>, space_dimension>& neighbor_distances
    ) const;

    // Total number of grid nodes
    int number_of_nodes() const;
    // Number of inner grid nodes
    int number_of_inner_nodes() const;
    // Number of inner boundary grid nodes
    int number_of_inner_boundary_nodes() const;
    // Number of outer boundary grid nodes
    int number_of_outer_boundary_nodes() const;

    // Grid Parameters
    // Get the number of grid points in radial direction
    int nr() const;
    // Get the number of angular divisions
    int ntheta() const;
    // Get the radius at a specific radial index
    const scalar_t& radius(const int r_index) const;
    // Get the angle at a specific angular index
    const scalar_t& theta(const int theta_index) const;
    // Get all radii and angles available which define the grid
    const std::vector<scalar_t>& radii() const;
    const std::vector<scalar_t>& angles() const;

    // Grid distances
    // Get the radial distance to the next consecutive radial node at a specified radial index.
    const scalar_t& r_dist(const int r_index) const;
    // Get the angular distance to the next consecutive angular node at a specified unwrapped angular index.
    const scalar_t& theta_dist(const int unwrapped_theta_index) const;
    // Get the vector of distances between consecutive radial nodes.
    const std::vector<scalar_t>& r_distances() const;
    // Get the vector of distances of the angular divisions.
    const std::vector<scalar_t>& theta_distances() const;

    // Circle/radial smoother division
    // Get the radius which splits the grid into circular and radial smoothing
    scalar_t smoother_splitting_radius() const;
    // Get the number of circles in the circular smoother.
    int numberSmootherCircles() const;
    // Get the length of the radial smoother lines. 
    int lengthSmootherRadial() const;
    // Get the number of nodes in circular smoother.
    int numberCircularSmootherNodes() const;
    // Get the number of nodes in radial smoother.
    int numberRadialSmootherNodes() const;

    int numberBlackSmootherCircles() const;
    int numberWhiteSmootherCircles() const;
    int numberBlackSmootherRadials() const;
    int numberWhiteSmootherRadials() const;
    
    // Implementation in src/PolarGrid/load_write_grid.cpp
    // Write the grid data to files specified for radii and angles with given precision.
    void writeToFile(
        const std::string& file_r, const std::string& file_theta, const int precision
    ) const;

private:
    // --------------- //
    // Private members //
    // --------------- //

    // We will use the convention.
    // radii = [R0, ..., R], angles = [0, ..., 2*pi]
    // Note that ntheta will be one less than the size of angles since 0 and 2pi or the same point.
    int nr_; // number of nodes in radial direction
    int ntheta_; // number of (unique) nodes in angular direction 
    bool is_ntheta_PowerOfTwo_;
    std::vector<scalar_t> radii_; // Vector of radial coordiantes
    std::vector<scalar_t> angles_; // Vector of angular coordinates

    // r_dist contains the distances between each consecutive radii division.
    // r_dist = [r_{1}-r_{0}, ..., r_{N}-r_{N-1}].
    std::vector<scalar_t> r_dist_; // size(r_dist) = nr() - 1

    // theta_dist contains the angles between each consecutive theta division.
    // Since we have a periodic boundary in theta direction, 
    // we have to make sure the index wraps around correctly when accessing it.
    // Here theta_0 = 0.0 and theta_N = 2*pi refer to the same point. 
    // theta_dist = [theta_{1}-theta_{0}, ..., theta_{N}-theta_{N-1}].
    std::vector<scalar_t> theta_dist_; // size(theta_dist) = ntheta()

    // Circle/radial smoother division
    scalar_t smoother_splitting_radius_; // Radius at which the grid is split into circular and radial smoothing
    int numberSmootherCircles_; // Number of smoother circles in the grid
    int lengthSmootherRadial_; // Length of the radial smoother lines. 
    int numberCircularSmootherNodes_; // Number of nodes in the circular smoother
    int numberRadialSmootherNodes_; // Number of nodes in the radial smoother

    int numberBlackSmootherCircles_;
    int numberWhiteSmootherCircles_;
    int numberBlackSmootherRadials_;
    int numberWhiteSmootherRadials_;

    // Relationship constraints:
    // radius(numberSmootherCircles) <= smoother_splitting_radius < radius(numberSmootherCircles + 1)
    // numberSmootherCircles + lengthSmootherRadial = nr()
    // numberCircularSmootherNodes + numberRadialSmootherNodes = number_of_nodes()

    // ------------------------ //
    // Private Helper Functions //
    // ------------------------ //

    // Check parameter validity
    void checkParameters(const std::vector<scalar_t>& radii, const std::vector<scalar_t>& angles) const;
    // Initialize r_dist_, theta_dist_
    void initializeDistances();

    // Initializes line splitting parameters for Circle/radial indexing.
    // splitting_radius: The radius value used for dividing the smoother into a circular and radial section.
    //      If std::nullopt, automatic line-splitting is enabled.
    //      If the splitting radius is less than R0, only Radial indexing is used.
    //      If the splitting radius is greater than or equal to R, only Circular indexing is used.
    void initializeLineSplitting(std::optional<scalar_t> splitting_radius);

    // Construct radial divisions for grid generation.
    // Parameters:
    //   R0: Inner radius of the grid.
    //   R: Outer radius of the grid.
    //   nr_exp: Exponential scaling factor for radial divisions.
    //   anisotropic_factor: Factor for adjusting anisotropy in radial divisions.
    //   alpha: Allow refining of the grid at r_jump, 
    //          the center point of the drop of the diffusion coefficient alpha.
    void constructRadialDivisions(
        const scalar_t R0, const scalar_t R, const int nr_exp, const double r_jump,
        const int anisotropic_factor, const alpha_coeff alpha);

    // Construct angular divisions for grid generation.
    // Parameters:
    //   ntheta_exp: Exponential scaling factor for angular divisions.
    //               Set to -1 to determine ntheta similar to nr.
    //   nr: Number of radial nodes (used for calculating ntheta if ntheta_exp is -1).
    void constructAngularDivisions(const int ntheta_exp, const int nr);

    // Refine the grid by dividing radial and angular divisions by 2.
    // Parameters:
    //   divideBy2: Number of times to divide both radial and angular divisions by 2.
    void refineGrid(const int divideBy2);
    std::vector<scalar_t> divideVector(const std::vector<scalar_t>& vec, const int divideBy2) const;

    // Help constrcut radii_ when an anisotropic radial division is requested
    // Implementation in src/PolarGrid/anisotropic_division.cpp
    void RadialAnisotropicDivision(std::vector<scalar_t>& r_temp, 
        const scalar_t R0, const scalar_t R, const int nr_exp, const double r_jump,
        const int anisotropic_factor, const alpha_coeff alpha
    ) const;

    // Implementation in src/PolarGrid/load_write_grid.cpp
    void writeVectorToFile(
        const std::string& filename, const std::vector<scalar_t>& vector, const int precision
    ) const;
    void loadVectorFromFile(
        const std::string& filename, std::vector<scalar_t>& vector
    ) const;
};

// ---------------------------------------------------- //
// Generates a coarser PolarGrid from a finer PolarGrid //
// ---------------------------------------------------- //
std::unique_ptr<PolarGrid> coarseningGrid(const PolarGrid& grid);

#include "polargrid.inl"  // Include the inline function definitions