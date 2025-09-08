#include "../../include/PolarGrid/polargrid.h"

// ------------ //
// Constructors //
// ------------ //

// Constructor to initialize grid using vectors of radii and angles.
PolarGrid::PolarGrid(const std::vector<double>& radii, const std::vector<double>& angles,
                     std::optional<double> splitting_radius)
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

// Constructor to initialize grid using parameters from GMGPolar.
PolarGrid::PolarGrid(const double& R0, const double& Rmax, const int nr_exp, const int ntheta_exp,
                     const double& refinement_radius, const int anisotropic_factor, const int divideBy2,
                     std::optional<double> splitting_radius)
{
    assert(R0 > 0.0 && Rmax > R0 && !equals(R0, Rmax));
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
}

// ------------------- //
// Constructor Helpers //
// ------------------- //

// Construct radial divisions for grid generation.
void PolarGrid::constructRadialDivisions(const double& R0, const double& R, const int nr_exp,
                                         const double& refinement_radius, const int anisotropic_factor)
{
    // r_temp contains the values before we refine one last time for extrapolation.
    // Therefore we first consider 2^(nr_exp-1) points.
    std::vector<double> r_temp;
    if (anisotropic_factor == 0) {
        int nr                  = pow(2, nr_exp - 1) + 1;
        double uniform_distance = (R - R0) / (nr - 1);
        assert(uniform_distance > 0.0);
        r_temp.resize(nr);
        for (int i = 0; i < nr - 1; i++) {
            r_temp[i] = R0 + i * uniform_distance;
        }
        r_temp[nr - 1] = R;
    }
    else {
        // Implementation in src/PolarGrid/anisotropic_division.cpp
        RadialAnisotropicDivision(r_temp, R0, R, nr_exp, refinement_radius, anisotropic_factor);
    }
    // Refine division in the middle for extrapolation
    nr_ = 2 * r_temp.size() - 1;
    radii_.resize(nr_);
    for (int i = 0; i < nr_; i++) {
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
    if (ntheta_exp < 0) {
        // Choose number of theta divisions similar to radial divisions.
        ntheta_ = pow(2, ceil(log2(nr)));
        // ntheta_ = pow(2, ceil(log2(nr-1)));
    }
    else {
        ntheta_ = pow(2, ntheta_exp);
    }
    is_ntheta_PowerOfTwo_ = (ntheta_ & (ntheta_ - 1)) == 0;
    // Note that currently ntheta_ = 2^k which allows us to do some optimizations when indexing.
    double uniform_distance = 2 * M_PI / ntheta_;
    angles_.resize(ntheta_ + 1);
    for (int i = 0; i < ntheta_; i++) {
        angles_[i] = i * uniform_distance;
    }
    angles_[ntheta_] = 2 * M_PI;
}

// divideBy2: Number of times to divide both radial and angular divisions by 2.
void PolarGrid::refineGrid(const int divideBy2)
{
    radii_                = divideVector(radii_, divideBy2);
    angles_               = divideVector(angles_, divideBy2);
    nr_                   = radii_.size();
    ntheta_               = angles_.size() - 1;
    is_ntheta_PowerOfTwo_ = (ntheta_ & (ntheta_ - 1)) == 0;
}

std::vector<double> PolarGrid::divideVector(const std::vector<double>& vec, const int divideBy2) const
{
    const double powerOfTwo = 1 << divideBy2;
    size_t vecSize          = vec.size();
    size_t resultSize       = vecSize + (vecSize - 1) * (powerOfTwo - 1);
    std::vector<double> result(resultSize);

    for (size_t i = 0; i < vecSize - 1; ++i) {
        size_t baseIndex  = i * powerOfTwo;
        result[baseIndex] = vec[i]; // Add the original value
        for (int j = 1; j < powerOfTwo; ++j) {
            double interpolated_value = vec[i] + j * (vec[i + 1] - vec[i]) / powerOfTwo;
            result[baseIndex + j]     = interpolated_value;
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
    for (int i = 0; i < nr() - 1; i++) {
        radial_spacings_[i] = radius(i + 1) - radius(i);
    }
    // angular_spacings contains the angles between each consecutive theta division.
    // Since we have a periodic boundary in theta direction,
    // we have to make sure the index wraps around correctly when accessing it.
    // Here theta_0 = 0.0 and theta_N = 2*pi refer to the same point.
    // angular_spacings = [theta_{1}-theta_{0}, ..., theta_{N}-theta_{N-1}].
    angular_spacings_.resize(ntheta());
    for (int i = 0; i < ntheta(); i++) {
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
    if (splitting_radius.has_value()) {
        if (splitting_radius.value() < radii_.front()) {
            number_smoother_circles_   = 0;
            length_smoother_radial_    = nr();
            smoother_splitting_radius_ = -1.0;
        }
        else {
            auto it = std::lower_bound(radii_.begin(), radii_.end(), splitting_radius.value());
            if (it != radii_.end()) {
                number_smoother_circles_   = std::distance(radii_.begin(), it);
                length_smoother_radial_    = nr() - number_smoother_circles_;
                smoother_splitting_radius_ = splitting_radius.value();
            }
            else {
                number_smoother_circles_   = nr();
                length_smoother_radial_    = 0;
                smoother_splitting_radius_ = radii_.back() + 1.0;
            }
        }
    }
    else {
        number_smoother_circles_ = 2; /* We assume numberSmootherCircles_ >= 2 in the further implementation */
        for (int i_r = 2; i_r < nr() - 2;
             i_r++) { /* We assume lengthSmootherRadial_ >= 3 in the further implementation */
            double uniform_theta_k = (2 * M_PI) / ntheta();
            double radius_r        = radius(i_r);
            double radial_dist_h   = radius(i_r) - radius(i_r - 1);
            ;
            double q = uniform_theta_k / radial_dist_h;
            if (q * radius_r > 1.0) {
                number_smoother_circles_ = i_r;
                break;
            }
        }
        /* The ExtrapolatedSmoother requires numberSmootherCircles_ >= 3 and lengthSmootherRadial_ >= 3. */
        if (number_smoother_circles_ < 3 && nr() > 5)
            number_smoother_circles_ = 3;

        length_smoother_radial_    = nr() - number_smoother_circles_;
        smoother_splitting_radius_ = radius(number_smoother_circles_);
    }

    number_circular_smoother_nodes_ = number_smoother_circles_ * ntheta();
    number_radial_smoother_nodes_   = length_smoother_radial_ * ntheta();

    assert(numberSmootherCircles() + lengthSmootherRadial() == nr());
    assert(numberCircularSmootherNodes() + numberRadialSmootherNodes() == numberOfNodes());
}

// ---------------------------------------------------- //
// Generates a coarser PolarGrid from a finer PolarGrid //
// ---------------------------------------------------- //

PolarGrid coarseningGrid(const PolarGrid& fineGrid)
{
    assert((fineGrid.nr() - 1) % 2 == 0 && (fineGrid.ntheta()) % 2 == 0);
    const int coarse_nr     = (fineGrid.nr() + 1) / 2;
    const int coarse_ntheta = fineGrid.ntheta() / 2;

    std::vector<double> coarse_r(coarse_nr);
    std::vector<double> coarse_theta(coarse_ntheta + 1);

    for (int i = 0; i < coarse_nr; i++) {
        coarse_r[i] = fineGrid.radius(2 * i);
    }
    for (int j = 0; j < coarse_ntheta + 1; j++) {
        coarse_theta[j] = fineGrid.theta(2 * j);
    }

    const bool use_same_splitting_radius = false;

    if (use_same_splitting_radius) {
        return PolarGrid(coarse_r, coarse_theta, fineGrid.smootherSplittingRadius());
    }
    else {
        return PolarGrid(coarse_r, coarse_theta);
    }
}

// ---------------- //
// Getter Functions //
// ---------------- //

const std::vector<double>& PolarGrid::radii() const
{
    return radii_;
}
const std::vector<double>& PolarGrid::angles() const
{
    return angles_;
}

// Get the radius at which the grid is split into circular and radial smoothing
double PolarGrid::smootherSplittingRadius() const
{
    return smoother_splitting_radius_;
}

// ------------------------------------- //
// Definition of node indexing.          //
// Based on the circular-radial smoother //
// ------------------------------------- //

/* OPTIMIZED INDEXING IS DEFINED IN include/PolarGrid/polargrid.inl */

int PolarGrid::index(const MultiIndex& position) const
{
    assert(position[0] >= 0 && position[0] < nr());
    assert(position[1] >= 0 && position[1] < ntheta());
    if (position[0] < numberSmootherCircles()) {
        return position[1] + ntheta() * position[0];
    }
    else {
        return numberCircularSmootherNodes() +
               (position[0] - numberSmootherCircles() + lengthSmootherRadial() * position[1]);
    }
}

MultiIndex PolarGrid::multiIndex(const int node_index) const
{
    assert(0 <= node_index && node_index < numberOfNodes());
    if (node_index < numberCircularSmootherNodes()) {
        auto result = std::div(node_index, ntheta());
        return MultiIndex(result.quot, result.rem);
    }
    else {
        auto result = std::div(node_index - numberCircularSmootherNodes(), lengthSmootherRadial());
        return MultiIndex(numberSmootherCircles() + result.rem, result.quot);
    }
}

Point PolarGrid::polarCoordinates(const MultiIndex& position) const
{
    assert(position[0] >= 0 && position[0] < nr());
    assert(position[1] >= 0 && position[1] < ntheta());
    return Point(radii_[position[0]], angles_[(position[1])]);
}

void PolarGrid::adjacentNeighborDistances(
    const MultiIndex& position, std::array<std::pair<double, double>, space_dimension>& neighbor_distance) const
{
    assert(position[0] >= 0 && position[0] < nr());
    assert(position[1] >= 0 && position[1] < ntheta());

    neighbor_distance[0].first  = (position[0] <= 0) ? 0.0 : radialSpacing(position[0] - 1);
    neighbor_distance[0].second = (position[0] >= nr() - 1) ? 0.0 : radialSpacing(position[0]);

    neighbor_distance[1].first  = angularSpacing(position[1] - 1);
    neighbor_distance[1].second = angularSpacing(position[1]);
}

void PolarGrid::adjacentNeighborsOf(const MultiIndex& position,
                                    std::array<std::pair<int, int>, space_dimension>& neighbors) const
{
    assert(position[0] >= 0 && position[0] < nr());
    assert(position[1] >= 0 && position[1] < ntheta());

    MultiIndex neigbor_position = position;
    neigbor_position[0] -= 1;
    neighbors[0].first = (neigbor_position[0] < 0) ? -1 : index(neigbor_position);

    neigbor_position = position;
    neigbor_position[0] += 1;
    neighbors[0].second = (neigbor_position[0] >= nr()) ? -1 : index(neigbor_position);

    neigbor_position = position;
    neigbor_position[1] -= 1;
    if (neigbor_position[1] < 0)
        neigbor_position[1] += ntheta();
    neighbors[1].first = index(neigbor_position);

    neigbor_position = position;
    neigbor_position[1] += 1;
    if (neigbor_position[1] >= ntheta())
        neigbor_position[1] -= ntheta();
    neighbors[1].second = index(neigbor_position);
}

void PolarGrid::diagonalNeighborsOf(const MultiIndex& position,
                                    std::array<std::pair<int, int>, space_dimension>& neighbors) const
{
    assert(position[0] >= 0 && position[0] < nr());
    assert(position[1] >= 0 && position[1] < ntheta());

    MultiIndex neigbor_position = position;
    neigbor_position[0] -= 1;
    neigbor_position[1] -= 1;
    if (neigbor_position[1] < 0)
        neigbor_position[1] += ntheta();
    neighbors[0].first = (neigbor_position[0] < 0) ? -1 : index(neigbor_position);

    neigbor_position = position;
    neigbor_position[0] += 1;
    neigbor_position[1] -= 1;
    if (neigbor_position[1] < 0)
        neigbor_position[1] += ntheta();
    neighbors[0].second = (neigbor_position[0] >= nr()) ? -1 : index(neigbor_position);

    neigbor_position = position;
    neigbor_position[0] -= 1;
    neigbor_position[1] += 1;
    if (neigbor_position[1] >= ntheta())
        neigbor_position[1] -= ntheta();
    neighbors[1].first = (neigbor_position[0] < 0) ? -1 : index(neigbor_position);

    neigbor_position = position;
    neigbor_position[0] += 1;
    neigbor_position[1] += 1;
    if (neigbor_position[1] >= ntheta())
        neigbor_position[1] -= ntheta();
    neighbors[1].second = (neigbor_position[0] >= nr()) ? -1 : index(neigbor_position);
}

// ------------------------ //
// Check parameter validity //
// ---------------------..- //

void PolarGrid::checkParameters(const std::vector<double>& radii, const std::vector<double>& angles) const
{
    if (radii.size() < 2) {
        throw std::invalid_argument("At least two radii are required.");
    }

    if (!std::all_of(radii.begin(), radii.end(), [](double r) {
            return r > 0.0;
        })) {
        throw std::invalid_argument("All radii must be greater than zero.");
    }

    if (std::adjacent_find(radii.begin(), radii.end(), std::greater_equal<double>()) != radii.end()) {
        throw std::invalid_argument("Radii must be strictly increasing.");
    }

    if (angles.size() < 3) {
        throw std::invalid_argument("At least two angles are required.");
    }

    if (!std::all_of(angles.begin(), angles.end(), [](double theta) {
            return theta >= 0.0;
        })) {
        throw std::invalid_argument("All angles must be non-negative.");
    }

    if (std::adjacent_find(angles.begin(), angles.end(), std::greater_equal<double>()) != angles.end()) {
        throw std::invalid_argument("Angles must be strictly increasing.");
    }

    if (!equals(angles.front(), 0.0)) {
        throw std::invalid_argument("First angle must be 0.");
    }

    if (!equals(angles.back(), 2 * M_PI)) {
        throw std::invalid_argument("Last angle must be 2*pi.");
    }

    // Additional constraint for our stencil. Not needed in general.
    if (ntheta() % 2 != 0) {
        throw std::invalid_argument(
            "ntheta must be divisible by 2.");
    }
}
