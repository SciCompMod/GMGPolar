#include "../../include/PolarGrid/polargrid.h"
#include <Kokkos_StdAlgorithms.hpp>
using namespace gmgpolar;
// ------------ //
// Constructors //
// ------------ //

// Constructor to initialize grid using vectors of radii and angles.
template <class MemorySpace>
PolarGrid<MemorySpace>::PolarGrid(std::vector<double> radii, std::vector<double> angles,
                                  std::optional<double> splitting_radius)
    : nr_(radii.size())
    , ntheta_(angles.size() - 1)
    , is_ntheta_PowerOfTwo_((ntheta_ & (ntheta_ - 1)) == 0)
    , radii_("radii", nr_)
    , angles_("angles", angles.size())

{
    auto radii_host  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), radii_);
    auto angles_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), angles_);
    // Copy from std vector to Kokkos view
    for (std::size_t i(0); i < radii.size(); ++i) {
        radii_host(i) = radii[i];
    }
    for (std::size_t i(0); i < angles.size(); ++i) {
        angles_host(i) = angles[i];
    }
    Kokkos::deep_copy(radii_, radii_host);
    Kokkos::deep_copy(angles_, angles_host);

    // Check parameter validity
    checkParameters(radii_host, angles_host);
    // Store distances to adjacent neighboring nodes.
    // Initializes radial_spacings_, angular_spacings_
    initializeDistances();
    // Initializes smoothers splitting radius for circle/radial indexing.
    initializeLineSplitting(splitting_radius);
}

// Constructor to initialize grid using parameters from GMGPolar.
template <class MemorySpace>
PolarGrid<MemorySpace>::PolarGrid(double R0, double Rmax, const int nr_exp, const int ntheta_exp,
                                  double refinement_radius, const int anisotropic_factor, const int divideBy2,
                                  std::optional<double> splitting_radius)
{
    assert(R0 > 0.0 && Rmax > R0 && !equals(R0, Rmax));
    assert((std::is_same_v<MemorySpace, Kokkos::HostSpace>));
    // Construct radii_ and angles_
    constructRadialDivisions(R0, Rmax, nr_exp, refinement_radius, anisotropic_factor);
    constructAngularDivisions(ntheta_exp, nr_);
    refineGrid(divideBy2);
    // Check parameter validity
    checkParameters(Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), radii_),
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), angles_));
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
template <class MemorySpace>
void PolarGrid<MemorySpace>::constructRadialDivisions(double R0, double R, const int nr_exp, double refinement_radius,
                                                      const int anisotropic_factor)
{
    // r_temp contains the values before we refine one last time for extrapolation.
    // Therefore we first consider 2^(nr_exp-1) points.
    AllocatableVector<double, MemorySpace> r_temp;
    if (anisotropic_factor == 0) {
        // nr = 2**(nr_exp-1) + 1
        int nr                  = (1 << (nr_exp - 1)) + 1;
        double uniform_distance = (R - R0) / (nr - 1);
        assert(uniform_distance > 0.0);
        r_temp = Vector<double, MemorySpace>("r_temp", nr);
        for (int i = 0; i < nr - 1; i++) {
            r_temp[i] = R0 + i * uniform_distance;
        }
        r_temp[nr - 1] = R;
    }
    else {
        // Implementation in src/PolarGrid/anisotropic_division.cpp
        r_temp = RadialAnisotropicDivision(R0, R, nr_exp, refinement_radius, anisotropic_factor);
    }
    // Refine division in the middle for extrapolation
    nr_    = 2 * r_temp.size() - 1;
    radii_ = AllocatableVector<double, MemorySpace>("radii_", nr_);
    for (int i = 0; i < nr_; i++) {
        if (!(i % 2))
            radii_[i] = r_temp[i / 2];
        else
            radii_[i] = 0.5 * (r_temp[(i - 1) / 2] + r_temp[(i + 1) / 2]);
    }
}

// Construct angular divisions for grid generation.
// Currently we dont allow anisotropic refinement in angular direction
template <class MemorySpace>
void PolarGrid<MemorySpace>::constructAngularDivisions(const int ntheta_exp, const int nr)
{
    if (ntheta_exp < 0) {
        // Choose number of theta divisions similar to radial divisions.
        // ntheta_ = 2**ceil(log2(nr))
        ntheta_ = 1 << static_cast<int>(ceil(log2(nr)));
        // ntheta_ = 1 << static_cast<int>(ceil(log2(nr-1)));
    }
    else {
        // ntheta_ = 2**ntheta_exp
        ntheta_ = 1 << ntheta_exp;
    }
    is_ntheta_PowerOfTwo_ = (ntheta_ & (ntheta_ - 1)) == 0;
    // Note that currently ntheta_ = 2^k which allows us to do some optimizations when indexing.
    double uniform_distance = 2 * M_PI / ntheta_;
    Kokkos::resize(angles_, ntheta_ + 1);
    for (int i = 0; i < ntheta_; i++) {
        angles_[i] = i * uniform_distance;
    }
    angles_[ntheta_] = 2 * M_PI;
}

// divideBy2: Number of times to divide both radial and angular divisions by 2.
template <class MemorySpace>
void PolarGrid<MemorySpace>::refineGrid(const int divideBy2)
{
    radii_                = divideVector(radii_, divideBy2);
    angles_               = divideVector(angles_, divideBy2);
    nr_                   = radii_.size();
    ntheta_               = angles_.size() - 1;
    is_ntheta_PowerOfTwo_ = (ntheta_ & (ntheta_ - 1)) == 0;
}

template <class MemorySpace>
Vector<double, MemorySpace> PolarGrid<MemorySpace>::divideVector(Vector<double, MemorySpace> vec,
                                                                 const int divideBy2) const
{
    const int powerOfTwo = 1 << divideBy2;
    size_t vecSize       = vec.size();
    size_t resultSize    = vecSize + (vecSize - 1) * (powerOfTwo - 1);
    Vector<double, MemorySpace> result("result", resultSize);

    for (size_t i = 0; i < vecSize - 1; ++i) {
        size_t baseIndex  = i * powerOfTwo;
        result[baseIndex] = vec[i]; // Add the original value
        for (int j = 1; j < powerOfTwo; ++j) {
            double interpolated_value = vec[i] + j * (vec[i + 1] - vec[i]) / static_cast<double>(powerOfTwo);
            result[baseIndex + j]     = interpolated_value;
        }
    }
    result[resultSize - 1] = vec(vec.size() - 1); // Add the last value of the original vector
    return result;
}

template <class MemorySpace>
void PolarGrid<MemorySpace>::initializeDistances()
{
    using ExecSpace = std::conditional_t<std::is_same_v<MemorySpace, Kokkos::HostSpace>,
                                         Kokkos::DefaultHostExecutionSpace, Kokkos::DefaultExecutionSpace>;
    // radial_spacings contains the distances between each consecutive radii division.
    // radial_spacings = [R_1-R0, ..., R_{N} - R_{N-1}].
    Kokkos::resize(radial_spacings_, nr() - 1);
    const Vector<double, MemorySpace>& radii           = radii_;
    const Vector<double, MemorySpace>& radial_spacings = radial_spacings_;
    Kokkos::parallel_for(
        "init radial_spacings_", Kokkos::RangePolicy<ExecSpace>(0, nr() - 1),
        KOKKOS_LAMBDA(int i) { radial_spacings[i] = radii[i + 1] - radii[i]; });
    // angular_spacings contains the angles between each consecutive theta division.
    // Since we have a periodic boundary in theta direction,
    // we have to make sure the index wraps around correctly when accessing it.
    // Here theta_0 = 0.0 and theta_N = 2*pi refer to the same point.
    // angular_spacings = [theta_{1}-theta_{0}, ..., theta_{N}-theta_{N-1}].
    Kokkos::resize(angular_spacings_, ntheta());
    const Vector<double, MemorySpace>& angles           = angles_;
    const Vector<double, MemorySpace>& angular_spacings = angular_spacings_;
    Kokkos::parallel_for(
        "init angular_spacings_", Kokkos::RangePolicy<ExecSpace>(0, ntheta()),
        KOKKOS_LAMBDA(int i) { angular_spacings[i] = angles[i + 1] - angles[i]; });
    Kokkos::fence();
}

// Initializes line splitting parameters for Circle/radial indexing.
// splitting_radius: The radius value used for dividing the smoother into a circular and radial section.
//      If std::nullopt, automatic line-splitting is enabled.
//      If the splitting radius is less than R0, only Radial indexing is used.
//      If the splitting radius is greater than or equal to R, only Circular indexing is used.
template <class MemorySpace>
void PolarGrid<MemorySpace>::initializeLineSplitting(std::optional<double> splitting_radius)
{
    auto radii_host(Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), radii_));
    if (splitting_radius.has_value()) {
        if (splitting_radius.value() < radii_(0)) {
            number_smoother_circles_   = 0;
            length_smoother_radial_    = nr();
            smoother_splitting_radius_ = -1.0;
        }
        else {
            auto start = Kokkos::Experimental::begin(radii_host);
            auto end   = Kokkos::Experimental::end(radii_host);
            auto it    = std::lower_bound(start, end, splitting_radius.value());

            if (it != end) {
                number_smoother_circles_   = std::distance(start, it);
                length_smoother_radial_    = nr() - number_smoother_circles_;
                smoother_splitting_radius_ = splitting_radius.value();
            }
            else {
                number_smoother_circles_   = nr();
                length_smoother_radial_    = 0;
                smoother_splitting_radius_ = radii_host(radii_host.size() - 1) + 1.0;
            }
        }
    }
    else {
		auto h_radius = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, radii_);
        number_smoother_circles_ = 2; /* We assume numberSmootherCircles_ >= 2 in the further implementation */
        for (int i_r = 2; i_r < nr() - 2;
             i_r++) { /* We assume lengthRadialSmoother_ >= 3 in the further implementation */
            const double uniform_theta_k = (2 * M_PI) / ntheta();
            const double radial_dist_h   = h_radius[i_r + 1] - h_radius[i_r];

            const double q = uniform_theta_k / radial_dist_h;
            if (q * h_radius[i_r] > 1.0) {
                number_smoother_circles_ = i_r;
                break;
            }
        }
        /* The ExtrapolatedSmoother requires numberSmootherCircles_ >= 3 and lengthRadialSmoother_ >= 3. */
        if (number_smoother_circles_ < 3 && nr() > 5)
            number_smoother_circles_ = 3;

        length_smoother_radial_    = nr() - number_smoother_circles_;
        smoother_splitting_radius_ = radii_host(number_smoother_circles_);
    }

    number_circular_smoother_nodes_ = number_smoother_circles_ * ntheta();
    number_radial_smoother_nodes_   = length_smoother_radial_ * ntheta();

    assert(numberSmootherCircles() + lengthRadialSmoother() == nr());
    assert(numberCircularSmootherNodes() + numberRadialSmootherNodes() == numberOfNodes());
}

// ---------------------------------------------------- //
// Generates a coarser PolarGrid from a finer PolarGrid //
// ---------------------------------------------------- //

namespace gmgpolar
{

template <class MemorySpace>
PolarGrid<MemorySpace> coarseningGrid(const PolarGrid<MemorySpace>& fineGrid)
{
    using ExecSpace = std::conditional_t<std::is_same_v<MemorySpace, Kokkos::HostSpace>,
                                         Kokkos::DefaultHostExecutionSpace, Kokkos::DefaultExecutionSpace>;

    assert((fineGrid.nr() - 1) % 2 == 0 && (fineGrid.ntheta()) % 2 == 0);
    const int coarse_nr     = (fineGrid.nr() + 1) / 2;
    const int coarse_ntheta = fineGrid.ntheta() / 2;

    Vector<double, MemorySpace> coarse_r("coarse_r", coarse_nr);
    Vector<double, MemorySpace> coarse_theta("coarse_theta", coarse_ntheta + 1);

    Kokkos::parallel_for(
        "coarse r", Kokkos::RangePolicy<ExecSpace>(0, coarse_nr),
        KOKKOS_LAMBDA(int i) { coarse_r[i] = fineGrid.radius(2 * i); });
    Kokkos::parallel_for(
        "coarse theta", Kokkos::RangePolicy<ExecSpace>(0, coarse_ntheta + 1),
        KOKKOS_LAMBDA(int j) { coarse_theta[j] = fineGrid.theta(2 * j); });

    const bool use_same_splitting_radius = false;

    if (use_same_splitting_radius) {
        return PolarGrid<MemorySpace>(coarse_r, coarse_theta, fineGrid.smootherSplittingRadius());
    }
    else {
        return PolarGrid<MemorySpace>(coarse_r, coarse_theta);
    }
}

template PolarGrid<Kokkos::HostSpace> coarseningGrid<Kokkos::HostSpace>(const PolarGrid<Kokkos::HostSpace>& grid);
#ifdef KOKKOS_ENABLE_CUDA
template PolarGrid<DefaultMemorySpace> coarseningGrid<DefaultMemorySpace>(const PolarGrid<DefaultMemorySpace>& grid);
#endif

} // namespace gmgpolar

// ------------------------ //
// Check parameter validity //
// ------------------------ //

template <class MemorySpace>
void PolarGrid<MemorySpace>::checkParameters(Vector<double, Kokkos::HostSpace> radii,
                                             Vector<double, Kokkos::HostSpace> angles) const
{
    auto radii_start = Kokkos::Experimental::begin(radii);
    auto radii_end   = Kokkos::Experimental::end(radii);
    if (radii.size() < 2) {
        throw std::invalid_argument("At least two radii are required.");
    }

    if (!std::all_of(radii_start, radii_end, [](double r) {
            return r > 0.0;
        })) {
        throw std::invalid_argument("All radii must be greater than zero.");
    }

    if (std::adjacent_find(radii_start, radii_end, std::greater_equal<double>()) != radii_end) {
        throw std::invalid_argument("Radii must be strictly increasing.");
    }

    if (angles.size() < 3) {
        throw std::invalid_argument("At least two angles are required.");
    }
    auto angles_start = Kokkos::Experimental::begin(angles);
    auto angles_end   = Kokkos::Experimental::end(angles);
    if (!std::all_of(angles_start, angles_end, [](double theta) {
            return theta >= 0.0;
        })) {
        throw std::invalid_argument("All angles must be non-negative.");
    }

    if (std::adjacent_find(angles_start, angles_end, std::greater_equal<double>()) != angles_end) {
        throw std::invalid_argument("Angles must be strictly increasing.");
    }

    if (!equals(angles(0), 0.0)) {
        throw std::invalid_argument("First angle must be 0.");
    }

    if (!equals(angles(angles.size() - 1), 2 * M_PI)) {
        throw std::invalid_argument("Last angle must be 2*pi.");
    }

    // Additional constraint for our stencil. Not needed in general.
    if (ntheta() % 2 != 0) {
        throw std::invalid_argument("ntheta must be divisible by 2.");
    }
}
