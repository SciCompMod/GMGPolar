#pragma once

#include <vector>
#include <optional>
#include <cassert>
#include <cuda_runtime.h>

class PolarGrid
{
public:
    /* Constructor to initialize grid using vectors of radii and angles. */
    PolarGrid(const std::vector<double>& radii, const std::vector<double>& angles,
              std::optional<double> splitting_radius = std::nullopt);

    /* Constructor to initialize grid using parameters from GMGPolar. */
    explicit PolarGrid(const double& R0, const double& Rmax, const int nr_exp, const int ntheta_exp,
                       const double& refinement_radius, const int anisotropic_factor, const int divideBy2,
                       std::optional<double> splitting_radius = std::nullopt);

    PolarGrid(const PolarGrid& other);
    PolarGrid(PolarGrid&& other) noexcept;

    PolarGrid& operator=(const PolarGrid& other);
    PolarGrid& operator=(PolarGrid&& other) noexcept;

    ~PolarGrid();

    __host__ __device__ __forceinline__ int wrapThetaIndex(const int unwrapped_theta_index) const;
    __host__ __device__ __forceinline__ int index(const int r_index, const int unwrapped_theta_index) const;
    __host__ __device__ __forceinline__ void multiIndex(const int node_index, int& r_index, int& theta_index) const;

    // Get the number of grid nodes
    __host__ __device__ __forceinline__ int numberOfNodes() const;
    // Get the number of grid points in radial direction
    __host__ __device__ __forceinline__ int nr() const;
    // Get the number of angular divisions
    __host__ __device__ __forceinline__ int ntheta() const;

    // Get the radius at a specific radial index
    __host__ __device__ __forceinline__ double radius(const int r_index) const;
    // Get the angle at a specific angular index
    __host__ __device__ __forceinline__ double theta(const int theta_index) const;

    // Get the radial distance to the next consecutive radial node at a specified radial index.
    __host__ __device__ __forceinline__ double radialSpacing(const int r_index) const;
    // Get the angular distance to the next consecutive angular node at a specified unwrapped angular index.
    __host__ __device__ __forceinline__ double angularSpacing(const int unwrapped_theta_index) const;

    // Circle/radial smoother division
    // Get the radius which splits the grid into circular and radial smoothing
    __host__ __device__ __forceinline__ double smootherSplittingRadius() const;
    // Get the number of circles in the circular smoother.
    __host__ __device__ __forceinline__ int numberSmootherCircles() const;
    // Get the length of the radial smoother lines.
    __host__ __device__ __forceinline__ int lengthSmootherRadial() const;
    // Get the number of nodes in circular smoother.
    __host__ __device__ __forceinline__ int numberCircularSmootherNodes() const;
    // Get the number of nodes in radial smoother.
    __host__ __device__ __forceinline__ int numberRadialSmootherNodes() const;

private:
    // We will use the convention:
    // radii = [R0, ..., R], angles = [0, ..., 2*pi]
    // Note that ntheta will be one less than the size of angles since 0 and 2pi are the same point.
    int nr_; // number of nodes in radial direction
    int ntheta_; // number of (unique) nodes in angular direction
    bool is_ntheta_PowerOfTwo_;
    std::vector<double> radii_; // Vector of radial coordiantes
    std::vector<double> angles_; // Vector of angular coordinates
    // GPU memory
    double* d_radii_;
    double* d_angles_;

    // radial_spacings_ contains the distances between each consecutive radii division.
    // radial_spacings_ = [r_{1}-r_{0}, ..., r_{N}-r_{N-1}].
    std::vector<double> radial_spacings_; // size(radial_spacings_) = nr() - 1
    // angular_spacings_ contains the angles between each consecutive theta division.
    // Since we have a periodic boundary in theta direction,
    // we have to make sure the index wraps around correctly when accessing it.
    // Here theta_0 = 0.0 and theta_N = 2*pi refer to the same point.
    // angular_spacings_ = [theta_{1}-theta_{0}, ..., theta_{N}-theta_{N-1}].
    std::vector<double> angular_spacings_; // size(angular_spacings_) = ntheta()
    // GPU memory
    double* d_radial_spacings_;
    double* d_angular_spacings_;

    // Circle/radial smoother division
    double smoother_splitting_radius_; // Radius at which the grid is split into circular and radial smoothing
    int number_smoother_circles_; // Number of smoother circles in the grid
    int length_smoother_radial_; // Length of the radial smoother lines.
    int number_circular_smoother_nodes_; // Number of nodes in the circular smoother
    int number_radial_smoother_nodes_; // Number of nodes in the radial smoother

    // ------------------------ //
    // Private Helper Functions //
    // ------------------------ //

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
    void constructRadialDivisions(const double& R0, const double& R, const int nr_exp, const double& refinement_radius,
                                  const int anisotropic_factor);
    // Construct angular divisions for grid generation.
    void constructAngularDivisions(const int ntheta_exp, const int nr);

    // Refine the grid by dividing radial and angular divisions by 2.
    void refineGrid(const int divideBy2);
    std::vector<double> divideVector(const std::vector<double>& vec, const int divideBy2) const;

    // Help constrcut radii_ when an anisotropic radial division is requested
    // Implementation in src/PolarGrid/anisotropic_division.cpp
    void RadialAnisotropicDivision(std::vector<double>& r_temp, const double& R0, const double& R, const int nr_exp,
                                   const double& refinement_radius, const int anisotropic_factor) const;

    void allocateDeviceMemory();
    void copyDataToDevice();
    void freeDeviceMemory();
};

#include "polargrid.inl"

// ---------------------------------------------------- //
// Generates a coarser PolarGrid from a finer PolarGrid //
// ---------------------------------------------------- //
PolarGrid coarseningGrid(const PolarGrid& grid);
