#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <memory>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include <omp.h>

#include "../Definitions/equals.h"

class PolarGrid
{
public:
    // Default constructor.
    explicit PolarGrid() = default;

    // Constructor to initialize grid using vectors of radii and angles.
    PolarGrid(const std::vector<double>& radii, const std::vector<double>& angles,
              std::optional<double> splitting_radius = std::nullopt);

    // Constructor to initialize grid using parameters from GMGPolar.
    explicit PolarGrid(double R0, double Rmax, const int nr_exp, const int ntheta_exp, double refinement_radius,
                       const int anisotropic_factor, const int divideBy2,
                       std::optional<double> splitting_radius = std::nullopt);

    // Optimized, inlined indexing.
    int wrapThetaIndex(const int unwrapped_theta_index) const;
    int index(const int r_index, const int unwrapped_theta_index) const;
    void multiIndex(const int node_index, int& r_index, int& theta_index) const;

    // Grid Parameters
    // Number of grid nodes
    int numberOfNodes() const;
    // Get the number of grid points in radial direction
    int nr() const;
    // Get the number of angular divisions
    int ntheta() const;
    // Get the radius at a specific radial index
    double radius(const int r_index) const;
    // Get the angle at a specific angular index
    double theta(const int theta_index) const;

    // Grid distances
    // Get the radial distance to the next consecutive radial node at a specified radial index.
    double radialSpacing(const int r_index) const;
    // Get the angular distance to the next consecutive angular node at a specified unwrapped angular index.
    double angularSpacing(const int unwrapped_theta_index) const;

    // Circle/radial smoother division
    // Get the radius which splits the grid into circular and radial smoothing
    double smootherSplittingRadius() const;
    // Get the number of circles in the circular smoother.
    int numberSmootherCircles() const;
    // Get the length of the radial smoother lines.
    int lengthSmootherRadial() const;
    // Get the number of nodes in circular smoother.
    int numberCircularSmootherNodes() const;
    // Get the number of nodes in radial smoother.
    int numberRadialSmootherNodes() const;

private:
    // We use the convention:
    // radii = [R0, ..., R], angles = [0, ..., 2*pi]
    // Note that ntheta will be one less than the size of angles since 0 and 2pi are the same point.
    int nr_; // number of nodes in radial direction
    int ntheta_; // number of (unique) nodes in angular direction
    bool is_ntheta_PowerOfTwo_;
    std::vector<double> radii_; // Vector of radial coordiantes
    std::vector<double> angles_; // Vector of angular coordinates

    // radial_spacings_ contains the distances between each consecutive radii division.
    // radial_spacings_ = [r_{1}-r_{0}, ..., r_{N}-r_{N-1}].
    std::vector<double> radial_spacings_; // size(radial_spacings_) = nr() - 1

    // angular_spacings_ contains the angles between each consecutive theta division.
    // Since we have a periodic boundary in theta direction,
    // we have to make sure the index wraps around correctly when accessing it.
    // Here theta_0 = 0.0 and theta_N = 2*pi refer to the same point.
    // angular_spacings_ = [theta_{1}-theta_{0}, ..., theta_{N}-theta_{N-1}].
    std::vector<double> angular_spacings_; // size(angular_spacings_) = ntheta()

    // Circle/radial smoother division
    double smoother_splitting_radius_; // Radius at which the grid is split into circular and radial smoothing
    int number_smoother_circles_; // Number of smoother circles in the grid
    int length_smoother_radial_; // Length of the radial smoother lines.
    int number_circular_smoother_nodes_; // Number of nodes in the circular smoother
    int number_radial_smoother_nodes_; // Number of nodes in the radial smoother

    /*
     * Relationship constraints:
     * - radius(numberSmootherCircles) <= smoother_splitting_radius < radius(numberSmootherCircles + 1)
     * - numberSmootherCircles + lengthSmootherRadial = nr()
     * - numberCircularSmootherNodes + numberRadialSmootherNodes = number_of_nodes()
     */

    // Check parameter validity
    void checkParameters(const std::vector<double>& radii, const std::vector<double>& angles) const;

    // Initialize radial_spacings_, angular_spacings_
    void initializeDistances();
    // Initializes line splitting parameters for Circle/radial indexing.
    // splitting_radius: The radius value used for dividing the smoother into a circular and radial section.
    //      If std::nullopt, automatic line-splitting is enabled.
    //      If the splitting radius is less than R0, only Radial indexing is used.
    //      If the splitting radius is greater than or equal to R, only Circular indexing is used.
    void initializeLineSplitting(std::optional<double> splitting_radius);

    // Construct radial divisions for grid generation.
    void constructRadialDivisions(double R0, double R, const int nr_exp, double refinement_radius,
                                  const int anisotropic_factor);
    // Construct angular divisions for grid generation.
    void constructAngularDivisions(const int ntheta_exp, const int nr);

    // Refine the grid by dividing radial and angular divisions by 2.
    void refineGrid(const int divideBy2);
    std::vector<double> divideVector(const std::vector<double>& vec, const int divideBy2) const;

    // Help constrcut radii_ when an anisotropic radial division is requested
    // Implementation in src/PolarGrid/anisotropic_division.cpp
    void RadialAnisotropicDivision(std::vector<double>& r_temp, double R0, double R, const int nr_exp,
                                   double refinement_radius, const int anisotropic_factor) const;
};

// ---------------------------------------------------- //
// Generates a coarser PolarGrid from a finer PolarGrid //
// ---------------------------------------------------- //
PolarGrid coarseningGrid(const PolarGrid& grid);

#include "polargrid.inl" // Include the inline function definitions
