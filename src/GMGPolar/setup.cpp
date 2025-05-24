#include "../../include/GMGPolar/gmgpolar.h"

void GMGPolar::setup()
{
    LIKWID_START("Setup");
    auto start_setup = std::chrono::high_resolution_clock::now();

    resetTimings();

    auto start_setup_createLevels = std::chrono::high_resolution_clock::now();

    if (stencil_distribution_method_ == StencilDistributionMethod::CPU_TAKE) {
        if(verbose_ > 0)
        {
            std::cout << "\nStencil version: TAKE" << std::endl;
        }
        if (!cache_density_profile_coefficients_ || !cache_domain_geometry_) {
            throw std::runtime_error("Error: Caching must be enabled for both density profile coefficients and domain "
                                     "geometry in 'Take' implementation strategy.");
        }
    }else{
        if(verbose_ > 0)
        {
            std::cout << "\nStencil version: GIVE" << std::endl;
        }
    }

    if(verbose_ > 0)
    {
        std::cout << "Geometry: ";
        if (typeid(*domain_geometry_) == typeid(CircularGeometry)) {
            std::cout << "Circular";
        }
        else if (typeid(*domain_geometry_) == typeid(ShafranovGeometry)) {
            std::cout << "Shafranov";
        }
        else if (typeid(*domain_geometry_) == typeid(CzarnyGeometry)) {
            std::cout << "Czarny";
        }
        else if (typeid(*domain_geometry_) == typeid(CulhamGeometry)) {
            std::cout << "Culham";
        }
        else {
            std::cout << "Unknown";
        }
        std::cout << std::endl;

        std::cout << "Coefficients: ";
        if (typeid(*density_profile_coefficients_) == typeid(PoissonCoefficients)) {
            std::cout << "PoissonCoefficients";
        }
        else if (typeid(*density_profile_coefficients_) == typeid(SonnendruckerCoefficients)) {
            std::cout << "SonnendruckerCoefficients, Alpha: Sonnendruecker, Beta: Zero";
        }
        else if (typeid(*density_profile_coefficients_) == typeid(SonnendruckerGyroCoefficients)) {
            std::cout << "SonnendruckerGyroCoefficients, Alpha: Sonnendruecker, Beta: Inverse";
        }
        else if (typeid(*density_profile_coefficients_) == typeid(ZoniCoefficients)) {
            std::cout << "ZoniCoefficients, Alpha: Zoni, Beta: Zero";
        }
        else if (typeid(*density_profile_coefficients_) == typeid(ZoniGyroCoefficients)) {
            std::cout << "ZoniGyroCoefficients, Alpha: Zoni, Beta: Inverse";
        }
        else if (typeid(*density_profile_coefficients_) == typeid(ZoniShiftedCoefficients)) {
            std::cout << "ZoniShiftedCoefficients, Alpha: ZoniShifted, Beta: Zero";
        }
        else if (typeid(*density_profile_coefficients_) == typeid(ZoniShiftedGyroCoefficients)) {
            std::cout << "ZoniShiftedGyroCoefficients, Alpha: ZoniShifted, Beta: Inverse";
        }        
        else {
            std::cout << "Unknown";
        }
        std::cout << std::endl;

        std::cout << "Problem: ";
        if (typeid(*exact_solution_) == typeid(CartesianR2_CircularGeometry) 
    || typeid(*exact_solution_) == typeid(CartesianR2_CzarnyGeometry)
    || typeid(*exact_solution_) == typeid(CartesianR2_ShafranovGeometry)) {
            std::cout << "CARTESIAN_R2";
        }
        else if (typeid(*exact_solution_) == typeid(CartesianR6_CircularGeometry) 
    || typeid(*exact_solution_) == typeid(CartesianR6_CzarnyGeometry)
    || typeid(*exact_solution_) == typeid(CartesianR6_ShafranovGeometry)) {
            std::cout << "CARTESIAN_R6";
        }
        else if (typeid(*exact_solution_) == typeid(PolarR6_CircularGeometry) 
    || typeid(*exact_solution_) == typeid(PolarR6_CulhamGeometry)
    || typeid(*exact_solution_) == typeid(PolarR6_CzarnyGeometry)
    || typeid(*exact_solution_) == typeid(PolarR6_ShafranovGeometry)) {
            std::cout << "POLAR_R6";
        }
        else if (typeid(*exact_solution_) == typeid(Refined_CircularGeometry) 
    || typeid(*exact_solution_) == typeid(Refined_CulhamGeometry)
    || typeid(*exact_solution_) == typeid(Refined_CzarnyGeometry)
    || typeid(*exact_solution_) == typeid(Refined_ShafranovGeometry)) {
            std::cout << "REFINED_RADIUS";
        }    
        else {
            std::cout << "Unknown";
        }
        std::cout << std::endl;  
    }


    // -------------------------------- //
    // Create the finest mesh (level 0) //
    // -------------------------------- //
    auto finest_grid = std::make_unique<PolarGrid>(createFinestGrid()); /* Implementation below */
    if (verbose_ > 0) {
        std::cout << "System of size (nr x ntheta) = (" << finest_grid->nr() << " x " << finest_grid->ntheta() << ")\n";
        std::cout << "on the coordinates (r x theta): (" << R0_ << ", " << Rmax_ << ") x (" << 0 << ", " << 2 * M_PI
                  << ")\n";

        std::cout << "Anisotropy factor: " << anisotropic_factor_ << std::endl;
        std::cout << "Dirichlet boundary (interior): " << DirBC_Interior_ << std::endl; 
    }
    if (paraview_)
        writeToVTK("output_finest_grid", *finest_grid);

    // ---------------------------------------------------------- //
    // Building PolarGrid and LevelCache for all multigrid levels //
    // ---------------------------------------------------------- //
    number_of_levels_ = chooseNumberOfLevels(*finest_grid); /* Implementation below */
    levels_.clear();
    levels_.reserve(number_of_levels_);

    if (verbose_ > 0)
        std::cout << "Number of levels: " << number_of_levels_ << "\n";

    int level_depth = 0;
    auto finest_levelCache =
        std::make_unique<LevelCache>(*finest_grid, *density_profile_coefficients_, *domain_geometry_,
                                     cache_density_profile_coefficients_, cache_domain_geometry_);
    levels_.emplace_back(level_depth, std::move(finest_grid), std::move(finest_levelCache), extrapolation_, FMG_);

    for (level_depth = 1; level_depth < number_of_levels_; level_depth++) {
        auto current_grid       = std::make_unique<PolarGrid>(coarseningGrid(levels_[level_depth - 1].grid()));
        auto current_levelCache = std::make_unique<LevelCache>(levels_[level_depth - 1], *current_grid);
        levels_.emplace_back(level_depth, std::move(current_grid), std::move(current_levelCache), extrapolation_, FMG_);
    }

    auto end_setup_createLevels = std::chrono::high_resolution_clock::now();
    t_setup_createLevels += std::chrono::duration<double>(end_setup_createLevels - start_setup_createLevels).count();

    if (paraview_)
        writeToVTK("output_coarsest_grid", levels_.back().grid());

    // ----------------------------------------------------------- //
    // Initializing the optimal number of threads for OpenMP tasks //
    // ----------------------------------------------------------- //
    threads_per_level_.resize(number_of_levels_, max_omp_threads_);
    for (int level_depth = 0; level_depth < number_of_levels_; level_depth++) {
        threads_per_level_[level_depth] = std::max(
            1,
            std::min(max_omp_threads_,
                     static_cast<int>(std::floor(max_omp_threads_ * std::pow(thread_reduction_factor_, level_depth)))));
    }

    if (verbose_ > 0)
        std::cout << "Maximum number of threads: " << max_omp_threads_ << "\n";

    interpolation_ = std::make_unique<Interpolation>(threads_per_level_, DirBC_Interior_);

    auto start_setup_rhs = std::chrono::high_resolution_clock::now();

    // ------------------------------------- //
    // Build rhs_f on Level 0 (finest Level) //
    // ------------------------------------- //
    LIKWID_STOP("Setup");
    build_rhs_f(levels_[0], levels_[0].rhs());
    LIKWID_START("Setup");

    /* ---------------- */
    /* Discretize rhs_f */
    /* ---------------- */
    int initial_rhs_f_levels = FMG_ ? number_of_levels_ : (extrapolation_ == ExtrapolationType::NONE ? 1 : 2);
    // Loop through the levels, injecting and discretizing rhs
    for (int level_depth = 0; level_depth < initial_rhs_f_levels; ++level_depth) {
        Level& current_level = levels_[level_depth];
        // Inject rhs if there is a next level
        if (level_depth + 1 < initial_rhs_f_levels) {
            Level& next_level = levels_[level_depth + 1];
            injection(level_depth, next_level.rhs(), current_level.rhs());
        }
        // Discretize the rhs for the current level
        discretize_rhs_f(current_level, current_level.rhs());
    }

    auto end_setup_rhs = std::chrono::high_resolution_clock::now();
    t_setup_rhs += std::chrono::duration<double>(end_setup_rhs - start_setup_rhs).count();

    // -------------------------------------------------------
    // Initializing various operators based on the level index
    for (int level_depth = 0; level_depth < number_of_levels_; level_depth++) {
        // ---------------------- //
        // Level 0 (finest Level) //
        // ---------------------- //
        if (level_depth == 0) {
            auto start_setup_smoother = std::chrono::high_resolution_clock::now();
            switch (extrapolation_) {
            case ExtrapolationType::NONE:
                full_grid_smoothing_ = true;
                levels_[level_depth].initializeSmoothing(*domain_geometry_, *density_profile_coefficients_,
                                                         DirBC_Interior_, threads_per_level_[level_depth],
                                                         stencil_distribution_method_);
                break;
            case ExtrapolationType::IMPLICIT_EXTRAPOLATION:
                full_grid_smoothing_ = false;
                levels_[level_depth].initializeExtrapolatedSmoothing(*domain_geometry_, *density_profile_coefficients_,
                                                                     DirBC_Interior_, threads_per_level_[level_depth],
                                                                     stencil_distribution_method_);
                break;
            case ExtrapolationType::IMPLICIT_FULL_GRID_SMOOTHING:
                full_grid_smoothing_ = true;
                levels_[level_depth].initializeSmoothing(*domain_geometry_, *density_profile_coefficients_,
                                                         DirBC_Interior_, threads_per_level_[level_depth],
                                                         stencil_distribution_method_);
                break;
            case ExtrapolationType::COMBINED:
                full_grid_smoothing_ = true;
                levels_[level_depth].initializeSmoothing(*domain_geometry_, *density_profile_coefficients_,
                                                         DirBC_Interior_, threads_per_level_[level_depth],
                                                         stencil_distribution_method_);
                levels_[level_depth].initializeExtrapolatedSmoothing(*domain_geometry_, *density_profile_coefficients_,
                                                                     DirBC_Interior_, threads_per_level_[level_depth],
                                                                     stencil_distribution_method_);
                break;
            default:
                full_grid_smoothing_ = false;
                levels_[level_depth].initializeSmoothing(*domain_geometry_, *density_profile_coefficients_,
                                                         DirBC_Interior_, threads_per_level_[level_depth],
                                                         stencil_distribution_method_);
                levels_[level_depth].initializeExtrapolatedSmoothing(*domain_geometry_, *density_profile_coefficients_,
                                                                     DirBC_Interior_, threads_per_level_[level_depth],
                                                                     stencil_distribution_method_);
                break;
            }
            auto end_setup_smoother = std::chrono::high_resolution_clock::now();
            t_setup_smoother += std::chrono::duration<double>(end_setup_smoother - start_setup_smoother).count();
            levels_[level_depth].initializeResidual(*domain_geometry_, *density_profile_coefficients_, DirBC_Interior_,
                                                    threads_per_level_[level_depth], stencil_distribution_method_);
        }
        // -------------------------- //
        // Level n-1 (coarsest Level) //
        // -------------------------- //
        else if (level_depth == number_of_levels_ - 1) {
            auto start_setup_directSolver = std::chrono::high_resolution_clock::now();
            levels_[level_depth].initializeDirectSolver(*domain_geometry_, *density_profile_coefficients_,
                                                        DirBC_Interior_, threads_per_level_[level_depth],
                                                        stencil_distribution_method_);
            auto end_setup_directSolver = std::chrono::high_resolution_clock::now();
            t_setup_directSolver +=
                std::chrono::duration<double>(end_setup_directSolver - start_setup_directSolver).count();
            levels_[level_depth].initializeResidual(*domain_geometry_, *density_profile_coefficients_, DirBC_Interior_,
                                                    threads_per_level_[level_depth], stencil_distribution_method_);
        }
        // ------------------- //
        // Intermediate levels //
        // ------------------- //
        else {
            auto start_setup_smoother = std::chrono::high_resolution_clock::now();
            levels_[level_depth].initializeSmoothing(*domain_geometry_, *density_profile_coefficients_, DirBC_Interior_,
                                                     threads_per_level_[level_depth], stencil_distribution_method_);
            auto end_setup_smoother = std::chrono::high_resolution_clock::now();
            t_setup_smoother += std::chrono::duration<double>(end_setup_smoother - start_setup_smoother).count();
            levels_[level_depth].initializeResidual(*domain_geometry_, *density_profile_coefficients_, DirBC_Interior_,
                                                    threads_per_level_[level_depth], stencil_distribution_method_);
        }
    }

    auto end_setup = std::chrono::high_resolution_clock::now();
    t_setup_total += std::chrono::duration<double>(end_setup - start_setup).count();
    LIKWID_STOP("Setup");
}

PolarGrid GMGPolar::createFinestGrid()
{
    PolarGrid finest_grid;

    if (load_grid_file_) {
        assert(!file_grid_radii_.empty() && !file_grid_angles_.empty());
        finest_grid = PolarGrid(file_grid_radii_, file_grid_angles_);
    }
    else {
        const double& refinement_radius =
            density_profile_coefficients_->getAlphaJump(); /* Radius of anisotropic grid refinement */
        std::optional<double> splitting_radius = std::nullopt; /* (Automatic) line splitting radius for the smoother */
        finest_grid = PolarGrid(R0_, Rmax_, nr_exp_, ntheta_exp_, refinement_radius, anisotropic_factor_, divideBy2_,
                                splitting_radius);
    }
    if (write_grid_file_) {
        const int precision = 18;
        finest_grid.writeToFile(file_grid_radii_, file_grid_angles_, precision);
    }
    return finest_grid;
}

int GMGPolar::chooseNumberOfLevels(const PolarGrid& finestGrid)
{
    const int minRadialNodes      = 5;
    const int minAngularDivisions = 4;

    // Minimum level for Multigrid
    const int multigridMinLevel = 2;

    // Calculate radial maximum level
    int radialNodes    = finestGrid.nr();
    int radialMaxLevel = 1;
    while ((radialNodes + 1) / 2 >= minRadialNodes && (radialNodes + 1) % 2 == 0) {
        radialNodes = (radialNodes + 1) / 2;
        radialMaxLevel++;
    }

    // Calculate angular maximum level
    int angularDivisions = finestGrid.ntheta();
    int angularMaxLevel  = 1;
    while (angularDivisions / 2 >= minAngularDivisions && angularDivisions % 2 == 0 &&
           (angularDivisions / 2) % 2 == 0) {
        angularDivisions = angularDivisions / 2;
        angularMaxLevel++;
    }

    /* Currently unused: Number of levels which guarantee linear scalability */
    const int linear_complexity_levels = std::ceil(
        (2.0 * std::log(static_cast<double>(finestGrid.numberOfNodes())) - std::log(3.0)) / (3.0 * std::log(4.0)));

    // Determine the number of levels as the minimum of radial maximum level, angular maximum level,
    // and the maximum levels specified.
    int levels = std::min(radialMaxLevel, angularMaxLevel);
    if (max_levels_ > 0)
        levels = std::min(max_levels_, levels);

    // Check if levels is less than Multigrid minimum level and throw an error
    if (levels < multigridMinLevel) {
        throw std::runtime_error("Number of possible levels is less than Multigrid minimum level");
    }

    return levels;
}

void GMGPolar::resetTimings()
{
    t_setup_total        = 0.0;
    t_setup_createLevels = 0.0;
    t_setup_rhs          = 0.0;
    t_setup_smoother     = 0.0;
    t_setup_directSolver = 0.0;

    t_solve_total                 = 0.0;
    t_solve_initial_approximation = 0.0;
    t_solve_multigrid_iterations  = 0.0;
    t_check_convergence           = 0.0;
    t_check_exact_error           = 0.0;

    t_avg_MGC_total         = 0.0;
    t_avg_MGC_preSmoothing  = 0.0;
    t_avg_MGC_postSmoothing = 0.0;
    t_avg_MGC_residual      = 0.0;
    t_avg_MGC_directSolver  = 0.0;
}