#pragma once

class DirectSolver;
class ResidualTakeCPU;
class SmootherTakeCPU;
class ExtrapolatedSmootherTakeCPU;
class ResidualTakeGPU;
class SmootherTakeGPU;
class ExtrapolatedSmootherTakeGPU;


#include <memory>
#include <omp.h>
#include <vector>

#include "../common/constants.h"

#include "../PolarGrid/polargrid.h"

#include "../LinearAlgebra/Vector/vector.h"
#include "../LinearAlgebra/Vector/gpu_vector.h"

#include "../InputFunctions/boundaryConditions.h"
#include "../InputFunctions/densityProfileCoefficients.h"
#include "../InputFunctions/domainGeometry.h"
#include "../InputFunctions/sourceTerm.h"


#include "../DirectSolver/directSolver.h"

#include "../Residual/ResidualTakeCPU/residual.h"
#include "../Smoother/SmootherTakeCPU/smoother.h"
#include "../ExtrapolatedSmoother/ExtrapolatedSmootherTakeCPU/extrapolatedSmoother.h"

#include "../Residual/ResidualTakeGPU/residual.h"
#include "../Smoother/SmootherTakeGPU/smoother.h"
#include "../ExtrapolatedSmoother/ExtrapolatedSmootherTakeGPU/extrapolatedSmoother.h"

class LevelCache;

class Level
{
public:
    // ----------- //
    // Constructor //

    explicit Level(const int level,
        const ProcessingType processing_type,
        std::unique_ptr<const PolarGrid> grid,
        std::unique_ptr<const LevelCache> level_cache,
        const ExtrapolationType extrapolation,
        const bool FMG);

    ~Level();
    // Move Constructor / Assignment
    Level(Level&& other) noexcept;
    Level& operator=(Level&& other) noexcept;
    // Copy Constructor (Deleted to avoid accidental copies)
    Level(const Level& other) = delete;

    // ---------------- //
    // Getter Functions //
    int level() const;
    ProcessingType processingType() const;
    const PolarGrid& grid() const;
    PolarGrid* device_grid() const;
    const LevelCache& levelCache() const;

    Vector<double>& rhs();
    const Vector<double>& rhs() const;
    Vector<double>& solution();
    const Vector<double>& solution() const;
    Vector<double>& residual();
    const Vector<double>& residual() const;
    Vector<double>& error_correction();
    const Vector<double>& error_correction() const;

    GPU_Vector<double>& GPU_rhs();
    const GPU_Vector<double>& GPU_rhs() const;
    GPU_Vector<double>& GPU_solution();
    const GPU_Vector<double>& GPU_solution() const;
    GPU_Vector<double>& GPU_residual();
    const GPU_Vector<double>& GPU_residual() const;
    GPU_Vector<double>& GPU_error_correction();
    const GPU_Vector<double>& GPU_error_correction() const;


    // // -------------- //
    // // Apply Residual //
    // void initializeResidual(const DomainGeometry& domain_geometry,
    //                         const DensityProfileCoefficients& density_profile_coefficients,
    //                         const bool DirBC_Interior);
    // void computeResidual(GPU_Vector<double>& result, const GPU_Vector<double>& rhs, const GPU_Vector<double>& x) const;

    // // ------------------- //
    // // Solve coarse System //
    // void initializeDirectSolver(const DomainGeometry& domain_geometry,
    //                             const DensityProfileCoefficients& density_profile_coefficients,
    //                             const bool DirBC_Interior);
    // void directSolveInPlace(Vector<double>& x) const;

    // // --------------- //
    // // Apply Smoothing //
    // void initializeSmoothing(const DomainGeometry& domain_geometry,
    //                          const DensityProfileCoefficients& density_profile_coefficients,
    //                          const bool DirBC_Interior,
    //                          const int num_omp_threads,
    //                          const StencilDistributionMethod stencil_distribution_method);
    // void smoothingInPlace(Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp) const;

    // // ---------------------------- //
    // // Apply Extrapolated Smoothing //
    // void initializeExtrapolatedSmoothing(const DomainGeometry& domain_geometry,
    //                                      const DensityProfileCoefficients& density_profile_coefficients,
    //                                      const bool DirBC_Interior,
    //                                      const int num_omp_threads,
    //                                      const StencilDistributionMethod stencil_distribution_method);
    // void extrapolatedSmoothingInPlace(Vector<double>& x, const Vector<double>& rhs, Vector<double>& temp) const;

private:
    int level_;
    ProcessingType processing_type_;
    std::unique_ptr<const PolarGrid> grid_;
    PolarGrid* device_grid_;
    std::unique_ptr<const LevelCache> level_cache_;
    ExtrapolationType extrapolation_;
    int FMG_;

    std::unique_ptr<DirectSolver> op_directSolver_;

    std::unique_ptr<ResidualTakeCPU> op_residual_CPU_;
    std::unique_ptr<SmootherTakeCPU> op_smoother_CPU_;
    std::unique_ptr<ExtrapolatedSmootherTakeCPU> op_extrapolated_smoother_CPU_;

    std::unique_ptr<ResidualTakeGPU> op_residual_GPU_;
    std::unique_ptr<SmootherTakeGPU> op_smoother_GPU_;
    std::unique_ptr<ExtrapolatedSmootherTakeGPU> op_extrapolated_smoother_GPU_;

    Vector<double> rhs_;
    Vector<double> solution_;
    Vector<double> residual_;
    Vector<double> error_correction_;

    GPU_Vector<double> gpu_rhs_;
    GPU_Vector<double> gpu_solution_;
    GPU_Vector<double> gpu_residual_;
    GPU_Vector<double> gpu_error_correction_;
};


class LevelCache
{
public:
    explicit LevelCache(const ProcessingType processing_type,
                        const PolarGrid& grid,
                        const DensityProfileCoefficients& density_profile_coefficients,
                        const DomainGeometry& domain_geometry);

    const std::vector<double>& sin_theta() const;
    const std::vector<double>& cos_theta() const;

    const std::vector<double>& coeff_alpha() const;
    const std::vector<double>& coeff_beta() const;

    const Vector<double>& arr() const;
    const Vector<double>& att() const;
    const Vector<double>& art() const;
    const Vector<double>& detDF() const;

    const GPU_Vector<double>& GPU_sin_theta() const;
    const GPU_Vector<double>& GPU_cos_theta() const;

    const GPU_Vector<double>& GPU_coeff_alpha() const;
    const GPU_Vector<double>& GPU_coeff_beta() const;

private:
    const ProcessingType processing_type_;

    std::vector<double> sin_theta_;
    std::vector<double> cos_theta_;

    std::vector<double> coeff_alpha_;
    std::vector<double> coeff_beta_;

    Vector<double> arr_;
    Vector<double> att_;
    Vector<double> art_;
    Vector<double> detDF_;

    GPU_Vector<double> gpu_sin_theta_;
    GPU_Vector<double> gpu_cos_theta_;

    GPU_Vector<double> gpu_coeff_alpha_;
    GPU_Vector<double> gpu_coeff_beta_;
};
