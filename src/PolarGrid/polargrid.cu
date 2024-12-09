#include "../../include/PolarGrid/polargrid.h"

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

#include "../../include/common/equals.h"

// ------------ //
// Constructors //
// ------------ //

// Constructor to initialize grid using vectors of radii and angles.
PolarGrid::PolarGrid(const std::vector<double>& radii, const std::vector<double>& angles, std::optional<double> splitting_radius)
    : nr_(radii.size())
    , ntheta_(angles.size() - 1)
    , is_ntheta_PowerOfTwo_((ntheta_ & (ntheta_ - 1)) == 0)
    , radii_(radii)
    , angles_(angles)
    , d_radii_(nullptr)
    , d_angles_(nullptr)
    , d_radial_spacings_(nullptr)
    , d_angular_spacings_(nullptr)
{
    // Check parameter validity
    checkParameters(radii, angles);
    // Store distances to adjacent neighboring nodes.
    // Initializes radial_spacings_, angular_spacings_
    initializeDistances();
    // Initializes smoothers splitting radius for circle/radial indexing.
    initializeLineSplitting(splitting_radius);

    allocateDeviceMemory();
    copyDataToDevice();
}

// Constructor to initialize grid using parameters from GMGPolar.
PolarGrid::PolarGrid(
    const double& R0,
    const double& Rmax,
    const int nr_exp,
    const int ntheta_exp,
    const double& refinement_radius,
    const int anisotropic_factor,
    const int divideBy2,
    std::optional<double> splitting_radius
)
    : d_radii_(nullptr)
    , d_angles_(nullptr)
    , d_radial_spacings_(nullptr)
    , d_angular_spacings_(nullptr)
{
    assert(0.0 < R0 && R0 < Rmax);
    // Construct radii_ and angles_
    constructRadialDivisions(R0, Rmax, nr_exp, refinement_radius, anisotropic_factor);
    constructAngularDivisions(ntheta_exp, nr_);
    refineGrid(divideBy2);
    // Check parameter validity
    checkParameters(radii_, angles_);
    // Store distances to adjacent neighboring nodes.
    // Initializes radial_spacings_, angular_spacings_
    initializeDistances();
    // Initializes smoothers splitting radius for circle/radial indexing.
    initializeLineSplitting(splitting_radius);

    allocateDeviceMemory();
    copyDataToDevice();
}

// Copy Constructor
PolarGrid::PolarGrid(const PolarGrid& other)
    : nr_(other.nr_), ntheta_(other.ntheta_), is_ntheta_PowerOfTwo_(other.is_ntheta_PowerOfTwo_),
        radii_(other.radii_), angles_(other.angles_), radial_spacings_(other.radial_spacings_),
        angular_spacings_(other.angular_spacings_),
        smoother_splitting_radius_(other.smoother_splitting_radius_),
        number_smoother_circles_(other.number_smoother_circles_),
        length_smoother_radial_(other.length_smoother_radial_),
        number_circular_smoother_nodes_(other.number_circular_smoother_nodes_),
        number_radial_smoother_nodes_(other.number_radial_smoother_nodes_)
{
    allocateDeviceMemory();
    copyDataToDevice();
}

// Copy Assignment Operator
PolarGrid& PolarGrid::operator=(const PolarGrid& other)
{
    if (this != &other)
    {
        freeDeviceMemory(); // Release existing resources

        // Copy data from the source
        nr_ = other.nr_;
        ntheta_ = other.ntheta_;
        is_ntheta_PowerOfTwo_ = other.is_ntheta_PowerOfTwo_;
        radii_ = other.radii_;
        angles_ = other.angles_;
        radial_spacings_ = other.radial_spacings_;
        angular_spacings_ = other.angular_spacings_;
        smoother_splitting_radius_ = other.smoother_splitting_radius_;
        number_smoother_circles_ = other.number_smoother_circles_;
        length_smoother_radial_ = other.length_smoother_radial_;
        number_circular_smoother_nodes_ = other.number_circular_smoother_nodes_;
        number_radial_smoother_nodes_ = other.number_radial_smoother_nodes_;

        allocateDeviceMemory();
        copyDataToDevice();
    }
    return *this;
}

// Move Constructor
PolarGrid::PolarGrid(PolarGrid&& other) noexcept
    : nr_(other.nr_), ntheta_(other.ntheta_), is_ntheta_PowerOfTwo_(other.is_ntheta_PowerOfTwo_),
        radii_(std::move(other.radii_)), angles_(std::move(other.angles_)),
        d_radii_(other.d_radii_), d_angles_(other.d_angles_),
        radial_spacings_(std::move(other.radial_spacings_)), angular_spacings_(std::move(other.angular_spacings_)),
        d_radial_spacings_(other.d_radial_spacings_), d_angular_spacings_(other.d_angular_spacings_),
        smoother_splitting_radius_(other.smoother_splitting_radius_),
        number_smoother_circles_(other.number_smoother_circles_),
        length_smoother_radial_(other.length_smoother_radial_),
        number_circular_smoother_nodes_(other.number_circular_smoother_nodes_),
        number_radial_smoother_nodes_(other.number_radial_smoother_nodes_)
{
    other.d_radii_ = nullptr;
    other.d_angles_ = nullptr;
    other.d_radial_spacings_ = nullptr;
    other.d_angular_spacings_ = nullptr;
}

// Move Assignment Operator
PolarGrid& PolarGrid::operator=(PolarGrid&& other) noexcept
{
    if (this != &other)
    {
        freeDeviceMemory(); // Release existing resources

        // Transfer ownership
        nr_ = other.nr_;
        ntheta_ = other.ntheta_;
        is_ntheta_PowerOfTwo_ = other.is_ntheta_PowerOfTwo_;
        radii_ = std::move(other.radii_);
        angles_ = std::move(other.angles_);
        d_radii_ = other.d_radii_;
        d_angles_ = other.d_angles_;
        radial_spacings_ = std::move(other.radial_spacings_);
        angular_spacings_ = std::move(other.angular_spacings_);
        d_radial_spacings_ = other.d_radial_spacings_;
        d_angular_spacings_ = other.d_angular_spacings_;
        smoother_splitting_radius_ = other.smoother_splitting_radius_;
        number_smoother_circles_ = other.number_smoother_circles_;
        length_smoother_radial_ = other.length_smoother_radial_;
        number_circular_smoother_nodes_ = other.number_circular_smoother_nodes_;
        number_radial_smoother_nodes_ = other.number_radial_smoother_nodes_;

        // Nullify source object
        other.d_radii_ = nullptr;
        other.d_angles_ = nullptr;
        other.d_radial_spacings_ = nullptr;
        other.d_angular_spacings_ = nullptr;
    }
    return *this;
}

// Destructor
PolarGrid::~PolarGrid()
{
    freeDeviceMemory();
}

// ------------------- //
// Constructor Helpers //
// ------------------- //

// Construct radial divisions for grid generation.
void PolarGrid::constructRadialDivisions(
    const double& R0, const double& R, const int nr_exp, const double& refinement_radius, const int anisotropic_factor)
{
    // r_temp contains the values before we refine one last time for extrapolation.
    // Therefore we first consider 2^(nr_exp-1) points.
    std::vector<double> r_temp;
    if (anisotropic_factor == 0)
    {
        int nr = pow(2, nr_exp - 1) + 1;
        double uniform_distance = (R - R0) / (nr - 1);
        assert(uniform_distance > 0.0);
        r_temp.resize(nr);
        for (int i = 0; i < nr - 1; i++)
        {
            r_temp[i] = R0 + i * uniform_distance;
        }
        r_temp[nr - 1] = R;
    }
    else
    {
        // Implementation in src/PolarGrid/anisotropic_division.cpp
        RadialAnisotropicDivision(r_temp, R0, R, nr_exp, refinement_radius, anisotropic_factor);
    }
    // Refine division in the middle for extrapolation
    nr_ = 2 * r_temp.size() - 1;
    radii_.resize(nr_);
    for (int i = 0; i < nr_; i++)
    {
        if (!(i % 2))
            radii_[i] = r_temp[i / 2];
        else
            radii_[i] = 0.5 * (r_temp[(i - 1) / 2] + r_temp[(i + 1) / 2]);
    }
}

// Construct angular divisions for grid generation.
// Currently we dont allow anisotropic refinement in angular direction
void PolarGrid::constructAngularDivisions(const int ntheta_exp, const int nr)
{
    if (ntheta_exp < 0)
    {
        // Choose number of theta divisions similar to radial divisions.
        ntheta_ = pow(2, ceil(log2(nr)));
        // ntheta_ = pow(2, ceil(log2(nr-1)));
    }
    else
    {
        ntheta_ = pow(2, ntheta_exp);
    }
    is_ntheta_PowerOfTwo_ = (ntheta_ & (ntheta_ - 1)) == 0;
    // Note that currently ntheta_ = 2^k which allows us to do some optimizations when indexing.
    double uniform_distance = 2 * M_PI / ntheta_;
    angles_.resize(ntheta_ + 1);
    for (int i = 0; i < ntheta_; i++)
    {
        angles_[i] = i * uniform_distance;
    }
    angles_[ntheta_] = 2 * M_PI;
}

// divideBy2: Number of times to divide both radial and angular divisions by 2.
void PolarGrid::refineGrid(const int divideBy2)
{
    radii_ = divideVector(radii_, divideBy2);
    angles_ = divideVector(angles_, divideBy2);
    nr_ = radii_.size();
    ntheta_ = angles_.size() - 1;
    is_ntheta_PowerOfTwo_ = (ntheta_ & (ntheta_ - 1)) == 0;
}

std::vector<double> PolarGrid::divideVector(const std::vector<double>& vec, const int divideBy2) const
{
    const double powerOfTwo = 1 << divideBy2;
    size_t vecSize = vec.size();
    size_t resultSize = vecSize + (vecSize - 1) * (powerOfTwo - 1);
    std::vector<double> result(resultSize);

    for (size_t i = 0; i < vecSize - 1; ++i)
    {
        size_t baseIndex = i * powerOfTwo;
        result[baseIndex] = vec[i]; // Add the original value
        for (int j = 1; j < powerOfTwo; ++j)
        {
            double interpolated_value = vec[i] + j * (vec[i + 1] - vec[i]) / powerOfTwo;
            result[baseIndex + j] = interpolated_value;
        }
    }
    result[resultSize - 1] = vec.back(); // Add the last value of the original vector
    return result;
}

void PolarGrid::initializeDistances()
{
    // radial_spacings contains the distances between each consecutive radii division.
    // radial_spacings = [R_1-R0, ..., R_{N} - R_{N-1}].
    radial_spacings_.resize(nr() - 1);
    for (int i = 0; i < nr() - 1; i++)
    {
        radial_spacings_[i] = radius(i + 1) - radius(i);
    }
    // angular_spacings contains the angles between each consecutive theta division.
    // Since we have a periodic boundary in theta direction,
    // we have to make sure the index wraps around correctly when accessing it.
    // Here theta_0 = 0.0 and theta_N = 2*pi refer to the same point.
    // angular_spacings = [theta_{1}-theta_{0}, ..., theta_{N}-theta_{N-1}].
    angular_spacings_.resize(ntheta());
    for (int i = 0; i < ntheta(); i++)
    {
        angular_spacings_[i] = theta(i + 1) - theta(i);
    }
}

// Initializes line splitting parameters for Circle/radial indexing.
// splitting_radius: The radius value used for dividing the smoother into a circular and radial section.
//      If std::nullopt, automatic line-splitting is enabled.
//      If the splitting radius is less than R0, only Radial indexing is used.
//      If the splitting radius is greater than or equal to R, only Circular indexing is used.
void PolarGrid::initializeLineSplitting(std::optional<double> splitting_radius)
{
    if (splitting_radius.has_value())
    {


        if (splitting_radius.value() < radii_.front())
        {
            number_smoother_circles_ = 0;
            length_smoother_radial_ = nr();
            smoother_splitting_radius_ = -1.0;
        }
        else
        {
            auto it = std::lower_bound(radii_.begin(), radii_.end(), splitting_radius.value());
            if (it != radii_.end())
            {
                number_smoother_circles_ = std::distance(radii_.begin(), it);
                length_smoother_radial_ = nr() - number_smoother_circles_;
                smoother_splitting_radius_ = splitting_radius.value();
            }
            else
            {
                number_smoother_circles_ = nr();
                length_smoother_radial_ = 0;
                smoother_splitting_radius_ = radii_.back() + 1.0;
            }
        }
    }
    else
    {
        number_smoother_circles_ = 2; /* We assume numberSmootherCircles_ >= 2 in the further implementation */
        for (int i_r = 2; i_r < nr() - 2; i_r++)
        { /* We assume lengthSmootherRadial_ >= 3 in the further implementation */
            double uniform_theta_k = (2 * M_PI) / ntheta();
            double radius_r = radius(i_r);
            double radial_dist_h = radius(i_r) - radius(i_r - 1);
            
            double q = uniform_theta_k / radial_dist_h;
            if (q * radius_r > 1.0)
            {
                number_smoother_circles_ = i_r;
                break;
            }
        }
        /* The ExtrapolatedSmoother requires numberSmootherCircles_ >= 3 and lengthSmootherRadial_ >= 3. */
        if (number_smoother_circles_ < 3 && nr() > 5) number_smoother_circles_ = 3;

        length_smoother_radial_ = nr() - number_smoother_circles_;
        smoother_splitting_radius_ = radius(number_smoother_circles_);
    }

    number_circular_smoother_nodes_ = number_smoother_circles_ * ntheta();
    number_radial_smoother_nodes_ = length_smoother_radial_ * ntheta();

    assert(numberSmootherCircles() + lengthSmootherRadial() == nr());
    assert(numberCircularSmootherNodes() + numberRadialSmootherNodes() == numberOfNodes());
}

// ------------------------ //
// Check parameter validity //
// ---------------------..- //

void PolarGrid::checkParameters(const std::vector<double>& radii, const std::vector<double>& angles) const
{
    if (radii.size() < 2)
    {
        throw std::invalid_argument("At least two radii are required.");
    }

    if (!std::all_of(radii.begin(), radii.end(), [](double r) { return r > 0.0; }))
    {
        throw std::invalid_argument("All radii must be greater than zero.");
    }

    if (std::adjacent_find(radii.begin(), radii.end(), std::greater_equal<double>()) != radii.end())
    {
        throw std::invalid_argument("Radii must be strictly increasing.");
    }

    if (angles.size() < 3)
    {
        throw std::invalid_argument("At least two angles are required.");
    }

    if (!std::all_of(angles.begin(), angles.end(), [](double theta) { return theta >= 0.0; }))
    {
        throw std::invalid_argument("All angles must be non-negative.");
    }

    if (std::adjacent_find(angles.begin(), angles.end(), std::greater_equal<double>()) != angles.end())
    {
        throw std::invalid_argument("Angles must be strictly increasing.");
    }

    if (!equals(angles.front(), 0.0))
    {
        throw std::invalid_argument("First angle must be 0.");
    }

    if (!equals(angles.back(), 2 * M_PI))
    {
        throw std::invalid_argument("Last angle must be 2*pi.");
    }

    // Additional constraint for our stencil. Not needed in general.
    if (!std::all_of(angles.begin(), angles.end(),
                     [&angles](double theta)
                     {
                         double opposite = theta + M_PI >= 2 * M_PI ? theta - M_PI : theta + M_PI;
                         return std::find_if(angles.begin(), angles.end(), [&opposite](double angle) { return equals(opposite, angle); }) !=
                                angles.end();
                     }))
    {
        throw std::invalid_argument("Each angle must have its opposite in the set:\n"
                                    "Every node in the interior ring needs to have an opposite neighboring node.");
    }
}

void PolarGrid::RadialAnisotropicDivision(
    std::vector<double>& r_temp,
    const double& R0,
    const double& R,
    const int nr_exp,
    const double& refinement_radius,
    const int anisotropic_factor) const
{
    // Calculate the percentage of refinement_radius.
    const double percentage = (refinement_radius - R0) / (R - R0);
    assert(percentage >= 0.0 && percentage <= 1.0);

    // 1) uniform division with nr=2^dummy_lognr - 2^aniso
    // 2) remaining nodes are added by refining the part centered around 2/3 of r
    std::set<double, std::greater<double>>::iterator itr, itr_p1;
    // very ugly anisotropy hack.... dividing recursively smaller and smaller number of cells

    /* uniform division of r in 2^nr_exp - 2^aniso */
    int dummy_lognr = nr_exp;
    int n_elems_equi = pow(2, dummy_lognr) - pow(2, anisotropic_factor);
    if (anisotropic_factor < 0 || n_elems_equi <= 0)
    {
        throw std::runtime_error("Please choose anisotropy factor a such that 2^fac_ani < 2^nr_exp.\n");
    }

    if ((anisotropic_factor % 2) == 1) // odd number of elements on an open circular disk is desired because of coarsening
        n_elems_equi++;
    double uniform_distance = (R - R0) / n_elems_equi;
    int nr = n_elems_equi + 1;
    std::vector<double> r_temp2 = std::vector<double>(nr);
    for (int i = 0; i < nr - 1; i++)
        r_temp2[i] = R0 + i * uniform_distance;
    r_temp2[nr - 1] = R;

    /* refine around 2/3 of r */
    int n_elems_refined = pow(2, anisotropic_factor);

    // edge
    int se;

    // Added by Allan Kuhn to fix a memory error
    if (floor(nr * percentage) > nr - (n_elems_refined / 2))
    {
        int new_aniso = log2(nr - floor(nr * percentage)) + 1;
        n_elems_refined = pow(2, new_aniso);
    }

    se = floor(nr * percentage) - n_elems_refined / 2;
    int ee = se + n_elems_refined;
    // takeout
    int st = ceil((double)n_elems_refined / 4.0 + 1) - 1;
    int et = floor(3 * ((double)n_elems_refined / 4.0));

    std::set<double> r_set;
    std::set<double> r_set_p1;
    int count = 0;
    for (int i = 0; i < n_elems_refined; i++)
    {
        r_set_p1.insert(r_temp2[se + i]);
        count++;
    }
    double half = uniform_distance / 2.0;
    for (int k = 0; k < anisotropic_factor; k++)
    {
        std::set<double> r_set_p1_tmp;
        itr_p1 = r_set_p1.begin();
        int r_size = count;
        count = 0;
        for (int i = 0; i < r_size - 1; i++)
        {
            r_set.insert((*itr_p1) + half);
            if (k < anisotropic_factor - 1 && i >= st && i < et)
            {
                r_set_p1_tmp.insert(*(itr_p1));
                r_set_p1_tmp.insert(*(itr_p1) + half);
                count += 2;
            }
            itr_p1++;
        }
        r_set_p1 = r_set_p1_tmp;
        half *= 0.5;
    }

    // such that the total size is 8*x+1 (or we do not refine)
    nr = nr + r_set.size();
    int shift = 0;
    shift = std::min(nr % 8 - 1, (int)r_set.size());
    itr = r_set.begin();
    std::advance(itr, shift);
    r_set.erase(r_set.begin(), itr);
    for (int i = 0; i < n_elems_refined; i++)
        r_set.insert(r_temp2[se + i]);

    // group all in r_tmp
    nr = n_elems_equi - n_elems_refined + r_set.size() + 1;

    r_temp.resize(nr);

    for (int i = 0; i < se; i++)
        r_temp[i] = r_temp2[i];
    itr = r_set.begin();
    for (int i = 0; i < (int)r_set.size(); i++)
    {
        r_temp[se + i] = *itr;
        itr++;
    }
    for (int i = 0; i < n_elems_equi - ee + 1; i++)
        r_temp[se + r_set.size() + i] = r_temp2[ee + i];
}

void PolarGrid::allocateDeviceMemory()
{
    if (cudaMalloc(&d_radii_, radii_.size() * sizeof(double)) != cudaSuccess) {
        throw std::runtime_error("Failed to allocate device memory for radii.");
    }
    if (cudaMalloc(&d_angles_, angles_.size() * sizeof(double)) != cudaSuccess) {
        throw std::runtime_error("Failed to allocate device memory for angles.");
    }

    if (cudaMalloc(&d_radial_spacings_, radial_spacings_.size() * sizeof(double)) != cudaSuccess) {
        throw std::runtime_error("Failed to allocate device memory for radial spacings.");
    }
    if (cudaMalloc(&d_angular_spacings_, angular_spacings_.size() * sizeof(double)) != cudaSuccess) {
        throw std::runtime_error("Failed to allocate device memory for angular spacings.");
    }
}

void PolarGrid::copyDataToDevice()
{
    if (cudaMemcpy(d_radii_, radii_.data(), radii_.size() * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) {
        throw std::runtime_error("Failed to copy radii to device.");
    }
    if (cudaMemcpy(d_angles_, angles_.data(), angles_.size() * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) {
        throw std::runtime_error("Failed to copy angles to device.");
    }
    if (cudaMemcpy(d_radial_spacings_, radial_spacings_.data(), radial_spacings_.size() * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) {
        throw std::runtime_error("Failed to copy radial spacings to device.");
    }
    if (cudaMemcpy(d_angular_spacings_, angular_spacings_.data(), angular_spacings_.size() * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) {
        throw std::runtime_error("Failed to copy angular spacings to device.");
    }
}
void PolarGrid::freeDeviceMemory()
{
    if (d_radii_) {
        cudaError_t err = cudaFree(d_radii_);
        if (err != cudaSuccess) {
            std::cerr << "Failed to free device memory for d_radii_: " << cudaGetErrorString(err) << std::endl;
        }
        d_radii_ = nullptr;
    }

    if (d_angles_) {
        cudaError_t err = cudaFree(d_angles_);
        if (err != cudaSuccess) {
            std::cerr << "Failed to free device memory for d_angles_: " << cudaGetErrorString(err) << std::endl;
        }
        d_angles_ = nullptr;
    }

    if (d_radial_spacings_) {
        cudaError_t err = cudaFree(d_radial_spacings_);
        if (err != cudaSuccess) {
            std::cerr << "Failed to free device memory for d_radial_spacings_: " << cudaGetErrorString(err) << std::endl;
        }
        d_radial_spacings_ = nullptr;
    }

    if (d_angular_spacings_) {
        cudaError_t err = cudaFree(d_angular_spacings_);
        if (err != cudaSuccess) {
            std::cerr << "Failed to free device memory for d_angular_spacings_: " << cudaGetErrorString(err) << std::endl;
        }
        d_angular_spacings_ = nullptr;
    }
}

// ---------------------------------------------------- //
// Generates a coarser PolarGrid from a finer PolarGrid //
// ---------------------------------------------------- //

PolarGrid coarseningGrid(const PolarGrid& fineGrid)
{
    assert((fineGrid.nr() - 1) % 2 == 0 && (fineGrid.ntheta()) % 2 == 0);
    const int coarse_nr = (fineGrid.nr() + 1) / 2;
    const int coarse_ntheta = fineGrid.ntheta() / 2;

    std::vector<double> coarse_r(coarse_nr);
    std::vector<double> coarse_theta(coarse_ntheta + 1);

    for (int i = 0; i < coarse_nr; i++)
    {
        coarse_r[i] = fineGrid.radius(2 * i);
    }
    for (int j = 0; j < coarse_ntheta + 1; j++)
    {
        coarse_theta[j] = fineGrid.theta(2 * j);
    }

    const bool use_same_splitting_radius = false;

    if (use_same_splitting_radius)
    {
        return PolarGrid(coarse_r, coarse_theta, fineGrid.smootherSplittingRadius());
    }
    else
    {        
        return PolarGrid(coarse_r, coarse_theta);
    }
}
