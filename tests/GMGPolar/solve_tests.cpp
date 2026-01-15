#include <gtest/gtest.h>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <random>

// Including the necessary header from the project
#include "../../include/GMGPolar/gmgpolar.h"

template <class T>
class GMGPolarTestCase;

// clang-format off
template <
    // Parameters
    class DomainGeometryType,
    class DensityProfileCoefficientsType,
    class BoundaryConditionsType,
    class SourceTermType,
    class ExactSolutionType,
    double R0_,
    double Rmax_,
    int nrExp_,
    int nthetaExp_,
    double refinementRadius_,
    int anisotropicFactor_,
    int divideBy2_,
    int verbose_,
    int maxOpenMPThreads_,
    bool DirBC_Interior_,
    StencilDistributionMethod stencilDistributionMethod_,
    bool cacheDensityProfileCoefficients_,
    bool cacheDomainGeometry_,
    ExtrapolationType extrapolation_,
    int maxLevels_,
    MultigridCycleType multigridCycle_,
    bool FMG_,
    int FMG_iterations_,
    MultigridCycleType FMG_cycle_,
    int maxIterations_,
    ResidualNormType residualNormType_,
    double absoluteTolerance_,
    double relativeTolerance_,
    // Results
    int expected_iterations_,
    double expected_l2_error_,
    double expected_inf_error_,
    double expected_residual_reduction_
>
class GMGPolarTestCase<
    std::tuple<
        DomainGeometryType,
        DensityProfileCoefficientsType,
        BoundaryConditionsType,
        SourceTermType,
        ExactSolutionType,
        std::integral_constant<double, R0_>,
        std::integral_constant<double, Rmax_>,
        std::integral_constant<int, nrExp_>,
        std::integral_constant<int, nthetaExp_>,
        std::integral_constant<double, refinementRadius_>,
        std::integral_constant<int, anisotropicFactor_>,
        std::integral_constant<int, divideBy2_>,
        std::integral_constant<int, verbose_>,
        std::integral_constant<int, maxOpenMPThreads_>,
        std::integral_constant<bool, DirBC_Interior_>,
        std::integral_constant<StencilDistributionMethod, stencilDistributionMethod_>,
        std::integral_constant<bool, cacheDensityProfileCoefficients_>,
        std::integral_constant<bool, cacheDomainGeometry_>,
        std::integral_constant<ExtrapolationType, extrapolation_>,
        std::integral_constant<int, maxLevels_>,
        std::integral_constant<MultigridCycleType, multigridCycle_>,
        std::integral_constant<bool, FMG_>,
        std::integral_constant<int, FMG_iterations_>,
        std::integral_constant<MultigridCycleType, FMG_cycle_>,
        std::integral_constant<int, maxIterations_>,
        std::integral_constant<ResidualNormType, residualNormType_>,
        std::integral_constant<double, absoluteTolerance_>,
        std::integral_constant<double, relativeTolerance_>,
        std::integral_constant<int, expected_iterations_>,
        std::integral_constant<double, expected_l2_error_>,
        std::integral_constant<double, expected_inf_error_>,
        std::integral_constant<double, expected_residual_reduction_>
    >
> : public testing::Test { // clang-format on
public:
    // Renaming static constexpr variables to avoid shadowing the template parameters
    using DomainGeometry                                                 = DomainGeometryType;
    using DensityProfileCoefficients                                     = DensityProfileCoefficientsType;
    using BoundaryConditions                                             = BoundaryConditionsType;
    using SourceTerm                                                     = SourceTermType;
    using ExactSolution                                                  = ExactSolutionType;
    static constexpr double R0                                           = R0_;
    static constexpr double Rmax                                         = Rmax_;
    static constexpr int nrExp                                           = nrExp_;
    static constexpr int nthetaExp                                       = nthetaExp_;
    static constexpr double refinementRadius                             = refinementRadius_;
    static constexpr int anisotropicFactor                               = anisotropicFactor_;
    static constexpr int divideBy2                                       = divideBy2_;
    static constexpr int verbose                                         = verbose_;
    static constexpr int maxOpenMPThreads                                = maxOpenMPThreads_;
    static constexpr bool DirBC_Interior                                 = DirBC_Interior_;
    static constexpr StencilDistributionMethod stencilDistributionMethod = stencilDistributionMethod_;
    static constexpr bool cacheDensityProfileCoefficients                = cacheDensityProfileCoefficients_;
    static constexpr bool cacheDomainGeometry                            = cacheDomainGeometry_;
    static constexpr ExtrapolationType extrapolation                     = extrapolation_;
    static constexpr int maxLevels                                       = maxLevels_;
    static constexpr MultigridCycleType multigridCycle                   = multigridCycle_;
    static constexpr bool FMG                                            = FMG_;
    static constexpr int FMG_iterations                                  = FMG_iterations_;
    static constexpr MultigridCycleType FMG_cycle                        = FMG_cycle_;
    static constexpr int maxIterations                                   = maxIterations_;
    static constexpr ResidualNormType residualNormType                   = residualNormType_;
    static constexpr double absoluteTolerance                            = absoluteTolerance_;
    static constexpr double relativeTolerance                            = relativeTolerance_;
    static constexpr int expected_iterations                             = expected_iterations_;
    static constexpr double expected_l2_error                            = expected_l2_error_;
    static constexpr double expected_inf_error                           = expected_inf_error_;
    static constexpr double expected_residual_reduction                  = expected_residual_reduction_;
};

// clang-format off
using gmgpolar_test_suite = testing::Types<
    /* --------------- */
    /* 1. Verbose Test */
    /* --------------- */
    std::tuple<
        CzarnyGeometry,
        ZoniShiftedGyroCoefficients,
        PolarR6_Boundary_CzarnyGeometry,
        PolarR6_ZoniShiftedGyro_CzarnyGeometry,
        PolarR6_CzarnyGeometry,
        std::integral_constant<double, 1e-7>, // R0
        std::integral_constant<double, 1.3>, // Rmax
        std::integral_constant<int, 3>, // nrExp
        std::integral_constant<int, -1>, // nthetaExp
        std::integral_constant<double, 0.66>, // refinementRadius
        std::integral_constant<int, 2>, // anisotropicFactor
        std::integral_constant<int, 1>, // divideBy2
        std::integral_constant<int, 1>, // verbose
        std::integral_constant<int, 1>, // maxOpenMPThreads
        std::integral_constant<bool, true>, // DirBC_Interior
        std::integral_constant<StencilDistributionMethod, StencilDistributionMethod::CPU_GIVE>, // StencilDistributionMethod
        std::integral_constant<bool, false>, // cacheDensityProfileCoefficient
        std::integral_constant<bool, false>, // cacheDomainGeometry
        std::integral_constant<ExtrapolationType, ExtrapolationType::IMPLICIT_EXTRAPOLATION>, // extrapolation
        std::integral_constant<int, 5>, // maxLevels
        std::integral_constant<MultigridCycleType, MultigridCycleType::V_CYCLE>, // multigridCycle
        std::integral_constant<bool, false>, // FMG
        std::integral_constant<int, 2>, // FMG_iterations
        std::integral_constant<MultigridCycleType, MultigridCycleType::F_CYCLE>, // FMG_cycle
        std::integral_constant<int, 50>, // maxIterations
        std::integral_constant<ResidualNormType, ResidualNormType::WEIGHTED_EUCLIDEAN>, // residualNormType
        std::integral_constant<double, 1e-6>, // absoluteTolerance
        std::integral_constant<double, 1e-6>, // relativeTolerance
        std::integral_constant<int, 4>, // expected_iterations
        std::integral_constant<double, 3e-6>, // expected_l2_error
        std::integral_constant<double, 8e-6>, // expected_inf_error
        std::integral_constant<double, 0.33> // expected_residual_reduction
    >,
    /* ---------------------------------------- */
    /* 2.1 Stencil Test: Give (Without caching) */
    /* ---------------------------------------- */
    std::tuple<
        CzarnyGeometry,
        ZoniShiftedGyroCoefficients,
        PolarR6_Boundary_CzarnyGeometry,
        PolarR6_ZoniShiftedGyro_CzarnyGeometry,
        PolarR6_CzarnyGeometry,
        std::integral_constant<double, 1e-7>, // R0
        std::integral_constant<double, 1.3>, // Rmax
        std::integral_constant<int, 4>, // nrExp
        std::integral_constant<int, -1>, // nthetaExp
        std::integral_constant<double, 0.5>, // refinementRadius
        std::integral_constant<int, 1>, // anisotropicFactor
        std::integral_constant<int, 0>, // divideBy2
        std::integral_constant<int, 0>, // verbose
        std::integral_constant<int, 8>, // maxOpenMPThreads
        std::integral_constant<bool, false>, // DirBC_Interior
        std::integral_constant<StencilDistributionMethod, StencilDistributionMethod::CPU_GIVE>, // StencilDistributionMethod
        std::integral_constant<bool, false>, // cacheDensityProfileCoefficient
        std::integral_constant<bool, false>, // cacheDomainGeometry
        std::integral_constant<ExtrapolationType, ExtrapolationType::IMPLICIT_EXTRAPOLATION>, // extrapolation
        std::integral_constant<int, 2>, // maxLevels
        std::integral_constant<MultigridCycleType, MultigridCycleType::V_CYCLE>, // multigridCycle
        std::integral_constant<bool, true>, // FMG
        std::integral_constant<int, 1>, // FMG_iterations
        std::integral_constant<MultigridCycleType, MultigridCycleType::V_CYCLE>, // FMG_cycle
        std::integral_constant<int, 50>, // maxIterations
        std::integral_constant<ResidualNormType, ResidualNormType::EUCLIDEAN>, // residualNormType
        std::integral_constant<double, 1e-8>, // absoluteTolerance
        std::integral_constant<double, 1e-8>, // relativeTolerance
        std::integral_constant<int, 26>, // expected_iterations
        std::integral_constant<double, 2e-6>, // expected_l2_error
        std::integral_constant<double, 9e-6>, // expected_inf_error
        std::integral_constant<double, 0.7> // expected_residual_reduction
    >,
    /* ------------------------------------- */
    /* 2.2 Stencil Test: Give (With caching) */
    /* ------------------------------------- */
    std::tuple<
        CzarnyGeometry,
        ZoniShiftedGyroCoefficients,
        PolarR6_Boundary_CzarnyGeometry,
        PolarR6_ZoniShiftedGyro_CzarnyGeometry,
        PolarR6_CzarnyGeometry,
        std::integral_constant<double, 1e-7>, // R0
        std::integral_constant<double, 1.3>, // Rmax
        std::integral_constant<int, 4>, // nrExp
        std::integral_constant<int, -1>, // nthetaExp
        std::integral_constant<double, 0.5>, // refinementRadius
        std::integral_constant<int, 1>, // anisotropicFactor
        std::integral_constant<int, 0>, // divideBy2
        std::integral_constant<int, 0>, // verbose
        std::integral_constant<int, 8>, // maxOpenMPThreads
        std::integral_constant<bool, false>, // DirBC_Interior
        std::integral_constant<StencilDistributionMethod, StencilDistributionMethod::CPU_GIVE>, // StencilDistributionMethod
        std::integral_constant<bool, true>, // cacheDensityProfileCoefficient
        std::integral_constant<bool, true>, // cacheDomainGeometry
        std::integral_constant<ExtrapolationType, ExtrapolationType::IMPLICIT_EXTRAPOLATION>, // extrapolation
        std::integral_constant<int, 2>, // maxLevels
        std::integral_constant<MultigridCycleType, MultigridCycleType::V_CYCLE>, // multigridCycle
        std::integral_constant<bool, true>, // FMG
        std::integral_constant<int, 1>, // FMG_iterations
        std::integral_constant<MultigridCycleType, MultigridCycleType::V_CYCLE>, // FMG_cycle
        std::integral_constant<int, 50>, // maxIterations
        std::integral_constant<ResidualNormType, ResidualNormType::EUCLIDEAN>, // residualNormType
        std::integral_constant<double, 1e-8>, // absoluteTolerance
        std::integral_constant<double, 1e-8>, // relativeTolerance
        std::integral_constant<int, 26>, // expected_iterations
        std::integral_constant<double, 2e-6>, // expected_l2_error
        std::integral_constant<double, 9e-6>, // expected_inf_error
        std::integral_constant<double, 0.7> // expected_residual_reduction
    >,
    /* ------------------------------------- */
    /* 2.3 Stencil Test: Take (With caching) */
    /* ------------------------------------- */
    std::tuple<
        CzarnyGeometry,
        ZoniShiftedGyroCoefficients,
        PolarR6_Boundary_CzarnyGeometry,
        PolarR6_ZoniShiftedGyro_CzarnyGeometry,
        PolarR6_CzarnyGeometry,
        std::integral_constant<double, 1e-7>, // R0
        std::integral_constant<double, 1.3>, // Rmax
        std::integral_constant<int, 4>, // nrExp
        std::integral_constant<int, -1>, // nthetaExp
        std::integral_constant<double, 0.5>, // refinementRadius
        std::integral_constant<int, 1>, // anisotropicFactor
        std::integral_constant<int, 0>, // divideBy2
        std::integral_constant<int, 0>, // verbose
        std::integral_constant<int, 8>, // maxOpenMPThreads
        std::integral_constant<bool, false>, // DirBC_Interior
        std::integral_constant<StencilDistributionMethod, StencilDistributionMethod::CPU_TAKE>, // StencilDistributionMethod
        std::integral_constant<bool, true>, // cacheDensityProfileCoefficient
        std::integral_constant<bool, true>, // cacheDomainGeometry
        std::integral_constant<ExtrapolationType, ExtrapolationType::IMPLICIT_EXTRAPOLATION>, // extrapolation
        std::integral_constant<int, 2>, // maxLevels
        std::integral_constant<MultigridCycleType, MultigridCycleType::V_CYCLE>, // multigridCycle
        std::integral_constant<bool, true>, // FMG
        std::integral_constant<int, 1>, // FMG_iterations
        std::integral_constant<MultigridCycleType, MultigridCycleType::V_CYCLE>, // FMG_cycle
        std::integral_constant<int, 50>, // maxIterations
        std::integral_constant<ResidualNormType, ResidualNormType::EUCLIDEAN>, // residualNormType
        std::integral_constant<double, 1e-8>, // absoluteTolerance
        std::integral_constant<double, 1e-8>, // relativeTolerance
        std::integral_constant<int, 26>, // expected_iterations
        std::integral_constant<double, 2e-6>, // expected_l2_error
        std::integral_constant<double, 9e-6>, // expected_inf_error
        std::integral_constant<double, 0.7> // expected_residual_reduction
    >, 
    /* No Extrapolation */
    std::tuple<
        CzarnyGeometry,
        ZoniShiftedGyroCoefficients,
        PolarR6_Boundary_CzarnyGeometry,
        PolarR6_ZoniShiftedGyro_CzarnyGeometry,
        PolarR6_CzarnyGeometry,
        std::integral_constant<double, 1e-7>, // R0
        std::integral_constant<double, 1.3>, // Rmax
        std::integral_constant<int, 4>, // nrExp
        std::integral_constant<int, -1>, // nthetaExp
        std::integral_constant<double, 0.6>, // refinementRadius
        std::integral_constant<int, 3>, // anisotropicFactor
        std::integral_constant<int, 0>, // divideBy2
        std::integral_constant<int, 0>, // verbose
        std::integral_constant<int, 8>, // maxOpenMPThreads
        std::integral_constant<bool, false>, // DirBC_Interior
        std::integral_constant<StencilDistributionMethod, StencilDistributionMethod::CPU_TAKE>, // StencilDistributionMethod
        std::integral_constant<bool, true>, // cacheDensityProfileCoefficient
        std::integral_constant<bool, true>, // cacheDomainGeometry
        std::integral_constant<ExtrapolationType, ExtrapolationType::NONE>, // extrapolation
        std::integral_constant<int, 2>, // maxLevels
        std::integral_constant<MultigridCycleType, MultigridCycleType::V_CYCLE>, // multigridCycle
        std::integral_constant<bool, true>, // FMG
        std::integral_constant<int, 1>, // FMG_iterations
        std::integral_constant<MultigridCycleType, MultigridCycleType::W_CYCLE>, // FMG_cycle
        std::integral_constant<int, 50>, // maxIterations
        std::integral_constant<ResidualNormType, ResidualNormType::INFINITY_NORM>, // residualNormType
        std::integral_constant<double, 1e-12>, // absoluteTolerance
        std::integral_constant<double, 1e-10>, // relativeTolerance
        std::integral_constant<int, 12>, // expected_iterations
        std::integral_constant<double, 6e-6>, // expected_l2_error
        std::integral_constant<double, 2e-5>, // expected_inf_error
        std::integral_constant<double, 0.3> // expected_residual_reduction
    >,
    /* Implicit Extrapolation */
    std::tuple<
        CzarnyGeometry,
        ZoniShiftedGyroCoefficients,
        PolarR6_Boundary_CzarnyGeometry,
        PolarR6_ZoniShiftedGyro_CzarnyGeometry,
        PolarR6_CzarnyGeometry,
        std::integral_constant<double, 1e-6>, // R0
        std::integral_constant<double, 1.5>, // Rmax
        std::integral_constant<int, 3>, // nrExp
        std::integral_constant<int, -1>, // nthetaExp
        std::integral_constant<double, 0.66>, // refinementRadius
        std::integral_constant<int, 0>, // anisotropicFactor
        std::integral_constant<int, 0>, // divideBy2
        std::integral_constant<int, 0>, // verbose
        std::integral_constant<int, 1>, // maxOpenMPThreads
        std::integral_constant<bool, false>, // DirBC_Interior
        std::integral_constant<StencilDistributionMethod, StencilDistributionMethod::CPU_GIVE>, // StencilDistributionMethod
        std::integral_constant<bool, true>, // cacheDensityProfileCoefficient
        std::integral_constant<bool, false>, // cacheDomainGeometry
        std::integral_constant<ExtrapolationType, ExtrapolationType::IMPLICIT_EXTRAPOLATION>, // extrapolation
        std::integral_constant<int, 2>, // maxLevels
        std::integral_constant<MultigridCycleType, MultigridCycleType::W_CYCLE>, // multigridCycle
        std::integral_constant<bool, true>, // FMG
        std::integral_constant<int, 2>, // FMG_iterations
        std::integral_constant<MultigridCycleType, MultigridCycleType::F_CYCLE>, // FMG_cycle
        std::integral_constant<int, 50>, // maxIterations
        std::integral_constant<ResidualNormType, ResidualNormType::EUCLIDEAN>, // residualNormType
        std::integral_constant<double, 1e-8>, // absoluteTolerance
        std::integral_constant<double, 1e-6>, // relativeTolerance
        std::integral_constant<int, 32>, // expected_iterations
        std::integral_constant<double, 5e-4>, // expected_l2_error
        std::integral_constant<double, 2e-3>, // expected_inf_error
        std::integral_constant<double, 0.7> // expected_residual_reduction
    >,
    /* Combined Extrapolation */
    std::tuple<
        CzarnyGeometry,
        ZoniShiftedGyroCoefficients,
        PolarR6_Boundary_CzarnyGeometry,
        PolarR6_ZoniShiftedGyro_CzarnyGeometry,
        PolarR6_CzarnyGeometry,
        std::integral_constant<double, 1e-6>, // R0
        std::integral_constant<double, 1.5>, // Rmax
        std::integral_constant<int, 3>, // nrExp
        std::integral_constant<int, 3>, // nthetaExp
        std::integral_constant<double, 0.66>, // refinementRadius
        std::integral_constant<int, 1>, // anisotropicFactor
        std::integral_constant<int, 0>, // divideBy2
        std::integral_constant<int, 0>, // verbose
        std::integral_constant<int, 6>, // maxOpenMPThreads
        std::integral_constant<bool, true>, // DirBC_Interior
        std::integral_constant<StencilDistributionMethod, StencilDistributionMethod::CPU_GIVE>, // StencilDistributionMethod
        std::integral_constant<bool, false>, // cacheDensityProfileCoefficient
        std::integral_constant<bool, true>, // cacheDomainGeometry
        std::integral_constant<ExtrapolationType, ExtrapolationType::COMBINED>, // extrapolation
        std::integral_constant<int, 3>, // maxLevels
        std::integral_constant<MultigridCycleType, MultigridCycleType::F_CYCLE>, // multigridCycle
        std::integral_constant<bool, false>, // FMG
        std::integral_constant<int, 0>, // FMG_iterations
        std::integral_constant<MultigridCycleType, MultigridCycleType::F_CYCLE>, // FMG_cycle
        std::integral_constant<int, 50>, // maxIterations
        std::integral_constant<ResidualNormType, ResidualNormType::EUCLIDEAN>, // residualNormType
        std::integral_constant<double, 1e-8>, // absoluteTolerance
        std::integral_constant<double, 1e-6>, // relativeTolerance
        std::integral_constant<int, 38>, // expected_iterations
        std::integral_constant<double, 2e-3>, // expected_l2_error
        std::integral_constant<double, 3e-3>, // expected_inf_error
        std::integral_constant<double, 0.7> // expected_residual_reduction
    >,
    /* Implicit Extrapolation with full grid smoothing */
    std::tuple<
        CzarnyGeometry,
        ZoniShiftedGyroCoefficients,
        PolarR6_Boundary_CzarnyGeometry,
        PolarR6_ZoniShiftedGyro_CzarnyGeometry,
        PolarR6_CzarnyGeometry,
        std::integral_constant<double, 1e-5>, // R0
        std::integral_constant<double, 1.0>, // Rmax
        std::integral_constant<int, 4>, // nrExp
        std::integral_constant<int, 4>, // nthetaExp
        std::integral_constant<double, 0.66>, // refinementRadius
        std::integral_constant<int, 1>, // anisotropicFactor
        std::integral_constant<int, 0>, // divideBy2
        std::integral_constant<int, 0>, // verbose
        std::integral_constant<int, 6>, // maxOpenMPThreads
        std::integral_constant<bool, true>, // DirBC_Interior
        std::integral_constant<StencilDistributionMethod, StencilDistributionMethod::CPU_GIVE>, // StencilDistributionMethod
        std::integral_constant<bool, false>, // cacheDensityProfileCoefficient
        std::integral_constant<bool, true>, // cacheDomainGeometry
        std::integral_constant<ExtrapolationType, ExtrapolationType::IMPLICIT_FULL_GRID_SMOOTHING>, // extrapolation
        std::integral_constant<int, 3>, // maxLevels
        std::integral_constant<MultigridCycleType, MultigridCycleType::F_CYCLE>, // multigridCycle
        std::integral_constant<bool, false>, // FMG
        std::integral_constant<int, 0>, // FMG_iterations
        std::integral_constant<MultigridCycleType, MultigridCycleType::F_CYCLE>, // FMG_cycle
        std::integral_constant<int, 25>, // maxIterations
        std::integral_constant<ResidualNormType, ResidualNormType::WEIGHTED_EUCLIDEAN>, // residualNormType
        std::integral_constant<double, -1.0>, // absoluteTolerance
        std::integral_constant<double, -1.0>, // relativeTolerance
        std::integral_constant<int, 25>, // expected_iterations
        std::integral_constant<double, 2e-3>, // expected_l2_error
        std::integral_constant<double, 3e-3>, // expected_inf_error
        std::integral_constant<double, 0.7> // expected_residual_reduction
    >,
    std::tuple<
        CzarnyGeometry,
        ZoniShiftedGyroCoefficients,
        PolarR6_Boundary_CzarnyGeometry,
        PolarR6_ZoniShiftedGyro_CzarnyGeometry,
        PolarR6_CzarnyGeometry,
        std::integral_constant<double, 1e-5>, // R0
        std::integral_constant<double, 1.0>, // Rmax
        std::integral_constant<int, 4>, // nrExp
        std::integral_constant<int, 4>, // nthetaExp
        std::integral_constant<double, 0.66>, // refinementRadius
        std::integral_constant<int, 1>, // anisotropicFactor
        std::integral_constant<int, 0>, // divideBy2
        std::integral_constant<int, 0>, // verbose
        std::integral_constant<int, 12>, // maxOpenMPThreads
        std::integral_constant<bool, true>, // DirBC_Interior
        std::integral_constant<StencilDistributionMethod, StencilDistributionMethod::CPU_GIVE>, // StencilDistributionMethod
        std::integral_constant<bool, true>, // cacheDensityProfileCoefficient
        std::integral_constant<bool, true>, // cacheDomainGeometry
        std::integral_constant<ExtrapolationType, ExtrapolationType::NONE>, // extrapolation
        std::integral_constant<int, 2>, // maxLevels
        std::integral_constant<MultigridCycleType, MultigridCycleType::V_CYCLE>, // multigridCycle
        std::integral_constant<bool, true>, // FMG
        std::integral_constant<int, 0>, // FMG_iterations
        std::integral_constant<MultigridCycleType, MultigridCycleType::V_CYCLE>, // FMG_cycle
        std::integral_constant<int, 30>, // maxIterations
        std::integral_constant<ResidualNormType, ResidualNormType::INFINITY_NORM>, // residualNormType
        std::integral_constant<double, -1.0>, // absoluteTolerance
        std::integral_constant<double, 1e-7>, // relativeTolerance
        std::integral_constant<int, 10>, // expected_iterations
        std::integral_constant<double, 3e-4>, // expected_l2_error
        std::integral_constant<double, 9e-4>, // expected_inf_error
        std::integral_constant<double, 0.2> // expected_residual_reduction
    >,
    std::tuple<
        CzarnyGeometry,
        ZoniShiftedGyroCoefficients,
        PolarR6_Boundary_CzarnyGeometry,
        PolarR6_ZoniShiftedGyro_CzarnyGeometry,
        PolarR6_CzarnyGeometry,
        std::integral_constant<double, 1e-5>, // R0
        std::integral_constant<double, 1.3>, // Rmax
        std::integral_constant<int, 5>, // nrExp
        std::integral_constant<int, 5>, // nthetaExp
        std::integral_constant<double, 0.6>, // refinementRadius
        std::integral_constant<int, 2>, // anisotropicFactor
        std::integral_constant<int, 0>, // divideBy2
        std::integral_constant<int, 0>, // verbose
        std::integral_constant<int, 12>, // maxOpenMPThreads
        std::integral_constant<bool, false>, // DirBC_Interior
        std::integral_constant<StencilDistributionMethod, StencilDistributionMethod::CPU_TAKE>, // StencilDistributionMethod
        std::integral_constant<bool, true>, // cacheDensityProfileCoefficient
        std::integral_constant<bool, true>, // cacheDomainGeometry
        std::integral_constant<ExtrapolationType, ExtrapolationType::IMPLICIT_EXTRAPOLATION>, // extrapolation
        std::integral_constant<int, 2>, // maxLevels
        std::integral_constant<MultigridCycleType, MultigridCycleType::V_CYCLE>, // multigridCycle
        std::integral_constant<bool, true>, // FMG
        std::integral_constant<int, 2>, // FMG_iterations
        std::integral_constant<MultigridCycleType, MultigridCycleType::V_CYCLE>, // FMG_cycle
        std::integral_constant<int, 50>, // maxIterations
        std::integral_constant<ResidualNormType, ResidualNormType::WEIGHTED_EUCLIDEAN>, // residualNormType
        std::integral_constant<double, 1e-8>, // absoluteTolerance
        std::integral_constant<double, -1.0>, // relativeTolerance
        std::integral_constant<int, 14>, // expected_iterations
        std::integral_constant<double, 9e-5>, // expected_l2_error
        std::integral_constant<double, 3e-4>, // expected_inf_error
        std::integral_constant<double, 0.6> // expected_residual_reduction
    >,
    std::tuple<
        CzarnyGeometry,
        ZoniShiftedGyroCoefficients,
        PolarR6_Boundary_CzarnyGeometry,
        PolarR6_ZoniShiftedGyro_CzarnyGeometry,
        PolarR6_CzarnyGeometry,
        std::integral_constant<double, 1e-6>, // R0
        std::integral_constant<double, 1.5>, // Rmax
        std::integral_constant<int, 4>, // nrExp
        std::integral_constant<int, -1>, // nthetaExp
        std::integral_constant<double, 0.66>, // refinementRadius
        std::integral_constant<int, 2>, // anisotropicFactor
        std::integral_constant<int, 0>, // divideBy2
        std::integral_constant<int, 0>, // verbose
        std::integral_constant<int, 1>, // maxOpenMPThreads
        std::integral_constant<bool, false>, // DirBC_Interior
        std::integral_constant<StencilDistributionMethod, StencilDistributionMethod::CPU_GIVE>, // StencilDistributionMethod
        std::integral_constant<bool, true>, // cacheDensityProfileCoefficient
        std::integral_constant<bool, false>, // cacheDomainGeometry
        std::integral_constant<ExtrapolationType, ExtrapolationType::NONE>, // extrapolation
        std::integral_constant<int, 3>, // maxLevels
        std::integral_constant<MultigridCycleType, MultigridCycleType::W_CYCLE>, // multigridCycle
        std::integral_constant<bool, true>, // FMG
        std::integral_constant<int, 1>, // FMG_iterations
        std::integral_constant<MultigridCycleType, MultigridCycleType::W_CYCLE>, // FMG_cycle
        std::integral_constant<int, 50>, // maxIterations
        std::integral_constant<ResidualNormType, ResidualNormType::EUCLIDEAN>, // residualNormType
        std::integral_constant<double, 1e-9>, // absoluteTolerance
        std::integral_constant<double, 1e-8>, // relativeTolerance
        std::integral_constant<int, 7>, // expected_iterations
        std::integral_constant<double, 5e-6>, // expected_l2_error
        std::integral_constant<double, 2e-5>, // expected_inf_error
        std::integral_constant<double, 0.2> // expected_residual_reduction
    >
>;
// clang-format on

TYPED_TEST_SUITE(GMGPolarTestCase, gmgpolar_test_suite);

template <class TestFixture>
void run_gmgpolar()
{
    PolarGrid grid(TestFixture::R0, TestFixture::Rmax, TestFixture::nrExp, TestFixture::nthetaExp,
                   TestFixture::refinementRadius, TestFixture::anisotropicFactor, TestFixture::divideBy2);

    const double inverse_aspect_ratio_epsilon = 0.3;
    const double ellipticity_e                = 1.4;

    typename TestFixture::DomainGeometry domain(TestFixture::Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
    typename TestFixture::DensityProfileCoefficients profile_coefficients(TestFixture::Rmax, 0.0);
    typename TestFixture::BoundaryConditions boundary_conditions(TestFixture::Rmax, inverse_aspect_ratio_epsilon,
                                                                 ellipticity_e);
    typename TestFixture::SourceTerm source_term(TestFixture::Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);
    typename TestFixture::ExactSolution exact_solution(TestFixture::Rmax, inverse_aspect_ratio_epsilon, ellipticity_e);

    GMGPolar solver(grid, domain, profile_coefficients);

    bool paraview          = false;
    int preSmoothingSteps  = 1;
    int postSmoothingSteps = 1;

    // --- General solver output and visualization settings --- //
    solver.verbose(TestFixture::verbose);
    solver.paraview(paraview);
    // --- Parallelization and threading settings --- //
    solver.maxOpenMPThreads(TestFixture::maxOpenMPThreads);
    omp_set_num_threads(TestFixture::maxOpenMPThreads);
    // --- Numerical method setup --- //
    solver.DirBC_Interior(TestFixture::DirBC_Interior);
    solver.stencilDistributionMethod(TestFixture::stencilDistributionMethod);
    solver.cacheDensityProfileCoefficients(TestFixture::cacheDensityProfileCoefficients);
    solver.cacheDomainGeometry(TestFixture::cacheDomainGeometry);
    // --- Multigrid settings --- //
    solver.extrapolation(TestFixture::extrapolation);
    solver.maxLevels(TestFixture::maxLevels);
    solver.preSmoothingSteps(preSmoothingSteps);
    solver.postSmoothingSteps(postSmoothingSteps);
    solver.multigridCycle(TestFixture::multigridCycle);
    solver.FMG(TestFixture::FMG);
    solver.FMG_iterations(TestFixture::FMG_iterations);
    solver.FMG_cycle(TestFixture::FMG_cycle);
    // --- Iterative solver controls --- //
    solver.maxIterations(TestFixture::maxIterations);
    solver.residualNormType(TestFixture::residualNormType);
    solver.absoluteTolerance(TestFixture::absoluteTolerance);
    solver.relativeTolerance(TestFixture::relativeTolerance);

    // --- Finalize solver setup --- //
    solver.setup(); // (allocates internal data, prepares operators, etc.)

    // --- Provide optional exact solution --- //
    solver.setSolution(&exact_solution);
    // --- Solve Phase --- //
    solver.solve(boundary_conditions, source_term);

    // --- Retrieve solution and associated grid --- //
    Vector<double> solution        = solver.solution();
    const PolarGrid& solution_grid = solver.grid();

    if (TestFixture::verbose > 0) {
        solver.printTimings();
    }

    EXPECT_EQ(solver.verbose(), TestFixture::verbose);
    EXPECT_EQ(solver.paraview(), paraview);
    EXPECT_EQ(solver.maxOpenMPThreads(), TestFixture::maxOpenMPThreads);
    EXPECT_EQ(solver.DirBC_Interior(), TestFixture::DirBC_Interior);
    EXPECT_EQ(solver.stencilDistributionMethod(), TestFixture::stencilDistributionMethod);
    EXPECT_EQ(solver.cacheDensityProfileCoefficients(), TestFixture::cacheDensityProfileCoefficients);
    EXPECT_EQ(solver.cacheDomainGeometry(), TestFixture::cacheDomainGeometry);
    EXPECT_EQ(solver.extrapolation(), TestFixture::extrapolation);
    EXPECT_EQ(solver.maxLevels(), TestFixture::maxLevels);
    EXPECT_EQ(solver.preSmoothingSteps(), preSmoothingSteps);
    EXPECT_EQ(solver.postSmoothingSteps(), postSmoothingSteps);
    EXPECT_EQ(solver.multigridCycle(), TestFixture::multigridCycle);
    EXPECT_EQ(solver.FMG(), TestFixture::FMG);
    EXPECT_EQ(solver.FMG_iterations(), TestFixture::FMG_iterations);
    EXPECT_EQ(solver.FMG_cycle(), TestFixture::FMG_cycle);
    EXPECT_EQ(solver.maxIterations(), TestFixture::maxIterations);
    EXPECT_EQ(solver.residualNormType(), TestFixture::residualNormType);
    if (solver.absoluteTolerance().has_value()) {
        EXPECT_DOUBLE_EQ(solver.absoluteTolerance().value(), TestFixture::absoluteTolerance);
    }
    if (solver.relativeTolerance().has_value()) {
        EXPECT_DOUBLE_EQ(solver.relativeTolerance().value(), TestFixture::relativeTolerance);
    }

    // --- Verify results ---
    int number_of_iterations                             = solver.numberOfIterations();
    std::optional<double> exact_error_weighted_euclidean = solver.exactErrorWeightedEuclidean();
    std::optional<double> exact_infinity_error           = solver.exactErrorInfinity();
    double reductionFactor                               = solver.meanResidualReductionFactor();

    ASSERT_TRUE(exact_error_weighted_euclidean.has_value());
    ASSERT_TRUE(exact_infinity_error.has_value());

    EXPECT_EQ(number_of_iterations, TestFixture::expected_iterations);
    EXPECT_LE(exact_error_weighted_euclidean.value(), TestFixture::expected_l2_error);
    EXPECT_LE(exact_infinity_error.value(), TestFixture::expected_inf_error);
    if (solver.absoluteTolerance().has_value() || solver.relativeTolerance().has_value()) {
        EXPECT_LT(reductionFactor, TestFixture::expected_residual_reduction);
    }

    // --- Verify timings (all available must be non-negative) ---
    EXPECT_GE(solver.timeSetupTotal(), 0.0);
    EXPECT_GE(solver.timeSetupCreateLevels(), 0.0);
    EXPECT_GE(solver.timeSetupRHS(), 0.0);
    EXPECT_GE(solver.timeSetupSmoother(), 0.0);
    EXPECT_GE(solver.timeSetupDirectSolver(), 0.0);

    EXPECT_GE(solver.timeSolveTotal(), 0.0);
    EXPECT_GE(solver.timeSolveInitialApproximation(), 0.0);
    EXPECT_GE(solver.timeSolveMultigridIterations(), 0.0);
    EXPECT_GE(solver.timeCheckConvergence(), 0.0);
    EXPECT_GE(solver.timeCheckExactError(), 0.0);

    EXPECT_GE(solver.timeAvgMGCTotal(), 0.0);
    EXPECT_GE(solver.timeAvgMGCPreSmoothing(), 0.0);
    EXPECT_GE(solver.timeAvgMGCPostSmoothing(), 0.0);
    EXPECT_GE(solver.timeAvgMGCResidual(), 0.0);
    EXPECT_GE(solver.timeAvgMGCDirectSolver(), 0.0);
}

TYPED_TEST(GMGPolarTestCase, GeneralTest)
{
    run_gmgpolar<TestFixture>();
}
