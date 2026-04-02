/* ---------------------------------------------------------------------- */
/* Interpolation                                                          */
/* ---------------------------------------------------------------------- */
template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
void GMGPolar<DomainGeometry, DensityProfileCoefficients>::prolongation(int current_level, Vector<double> result,
                                                                        ConstVector<double> x) const
{
    assert(current_level < number_of_levels_ && 1 <= current_level);
    if (!interpolation_)
        throw std::runtime_error("Interpolation not initialized.");

    interpolation_->applyProlongation(levels_[current_level].grid(), levels_[current_level - 1].grid(), result, x);
}

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
void GMGPolar<DomainGeometry, DensityProfileCoefficients>::restriction(int current_level, Vector<double> result,
                                                                       ConstVector<double> x) const
{
    assert(current_level < number_of_levels_ - 1 && 0 <= current_level);
    if (!interpolation_)
        throw std::runtime_error("Interpolation not initialized.");

    interpolation_->applyRestriction(levels_[current_level].grid(), levels_[current_level + 1].grid(), result, x);
}

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
void GMGPolar<DomainGeometry, DensityProfileCoefficients>::injection(int current_level, Vector<double> result,
                                                                     ConstVector<double> x) const
{
    assert(current_level < number_of_levels_ - 1 && 0 <= current_level);
    if (!interpolation_)
        throw std::runtime_error("Interpolation not initialized.");

    interpolation_->applyInjection(levels_[current_level].grid(), levels_[current_level + 1].grid(), result, x);
}

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
void GMGPolar<DomainGeometry, DensityProfileCoefficients>::extrapolatedProlongation(int current_level,
                                                                                    Vector<double> result,
                                                                                    ConstVector<double> x) const
{
    assert(current_level < number_of_levels_ && 1 <= current_level);
    if (!interpolation_)
        throw std::runtime_error("Interpolation not initialized.");

    interpolation_->applyExtrapolatedProlongation(levels_[current_level].grid(), levels_[current_level - 1].grid(),
                                                  result, x);
}

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
void GMGPolar<DomainGeometry, DensityProfileCoefficients>::extrapolatedRestriction(int current_level,
                                                                                   Vector<double> result,
                                                                                   ConstVector<double> x) const
{
    assert(current_level < number_of_levels_ - 1 && 0 <= current_level);
    if (!interpolation_)
        throw std::runtime_error("Interpolation not initialized.");

    interpolation_->applyExtrapolatedRestriction(levels_[current_level].grid(), levels_[current_level + 1].grid(),
                                                 result, x);
}

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
void GMGPolar<DomainGeometry, DensityProfileCoefficients>::FMGInterpolation(int current_level, Vector<double> result,
                                                                            ConstVector<double> x) const
{
    assert(current_level < number_of_levels_ && 1 <= current_level);
    if (!interpolation_)
        throw std::runtime_error("Interpolation not initialized.");

    interpolation_->applyFMGInterpolation(levels_[current_level].grid(), levels_[current_level - 1].grid(), result, x);
}

/* ---------------------------------------------------------------------- */
/* Solution & Grid Access                                                 */
/* ---------------------------------------------------------------------- */
template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
Vector<double> GMGPolar<DomainGeometry, DensityProfileCoefficients>::solution()
{
    int level_depth = 0;
    return levels_[level_depth].solution();
}

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
ConstVector<double> GMGPolar<DomainGeometry, DensityProfileCoefficients>::solution() const
{
    int level_depth = 0;
    return levels_[level_depth].solution();
}

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
const PolarGrid& GMGPolar<DomainGeometry, DensityProfileCoefficients>::grid() const
{
    return grid_;
}

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
void GMGPolar<DomainGeometry, DensityProfileCoefficients>::setSolution(const ExactSolution* exact_solution)
{
    exact_solution_ = exact_solution;
}

// Print timing breakdown for setup, smoothing, coarse solve, etc.
template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
void GMGPolar<DomainGeometry, DensityProfileCoefficients>::printTimings() const
{
    // t_setup_rhs_ is neither included in t_setup_total_ and t_solve_total_.
    std::cout << "\n------------------" << std::endl;
    std::cout << "Timing Information" << std::endl;
    std::cout << "------------------" << std::endl;
    std::cout << "Setup Time: " << t_setup_total_ << " seconds" << std::endl;
    std::cout << "    Create Levels: " << t_setup_createLevels_ << " seconds" << std::endl;
    std::cout << "    Smoother: " << t_setup_smoother_ << " seconds" << std::endl;
    std::cout << "    Direct Solver: " << t_setup_directSolver_ << " seconds" << std::endl;
    std::cout << "    (Build rhs: " << t_setup_rhs_ << " seconds)" << std::endl;
    std::cout << "\nSolve Time: " << t_solve_total_ << " seconds" << std::endl;
    std::cout << "    Initial Approximation: " << t_solve_initial_approximation_ << " seconds" << std::endl;
    if (!PCG_) {
        std::cout << "    Multigrid Iterations: " << t_solve_multigrid_iterations_ << " seconds" << std::endl;
    }
    else {
        std::cout << "    Preconditioned Conjugate Gradient: "
                  << std::max(t_conjugate_gradient_ - t_check_convergence_ - t_check_exact_error_, 0.0) << " seconds"
                  << std::endl;
    }
    std::cout << "    Check Convergence: " << t_check_convergence_ << " seconds" << std::endl;
    std::cout << "    (Check Exact Error: " << t_check_exact_error_ << " seconds)" << std::endl;

    if (!PCG_) {
        std::cout << "\nAverage Multigrid Iteration: " << t_avg_MGC_total_ << " seconds" << std::endl;
        std::cout << "    PreSmoothing: " << t_avg_MGC_preSmoothing_ << " seconds" << std::endl;
        std::cout << "    PostSmoothing: " << t_avg_MGC_postSmoothing_ << " seconds" << std::endl;
        std::cout << "    Residual: " << t_avg_MGC_residual_ << " seconds" << std::endl;
        std::cout << "    DirectSolve: " << t_avg_MGC_directSolver_ << " seconds" << std::endl;
        std::cout << "    Other Computations: "
                  << std::max(t_avg_MGC_total_ - t_avg_MGC_preSmoothing_ - t_avg_MGC_postSmoothing_ -
                                  t_avg_MGC_residual_ - t_avg_MGC_directSolver_,
                              0.0)
                  << " seconds" << std::endl;
    }
}

// Number of iterations taken by last solve.
template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
int GMGPolar<DomainGeometry, DensityProfileCoefficients>::numberOfIterations() const
{
    return number_of_iterations_;
}

// Mean residual reduction factor per iteration.
template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
double GMGPolar<DomainGeometry, DensityProfileCoefficients>::meanResidualReductionFactor() const
{
    return mean_residual_reduction_factor_;
}

// Error norms (only available if exact solution was set).
template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
std::optional<double> GMGPolar<DomainGeometry, DensityProfileCoefficients>::exactErrorWeightedEuclidean() const
{
    if (exact_solution_) {
        return exact_errors_.back().first;
    }
    return std::nullopt;
}
template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
std::optional<double> GMGPolar<DomainGeometry, DensityProfileCoefficients>::exactErrorInfinity() const
{
    if (exact_solution_) {
        return exact_errors_.back().second;
    }
    return std::nullopt;
}

/* ---------------------------------------------------------------------- */
/* Visualization                                                          */
/* ---------------------------------------------------------------------- */
template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
void GMGPolar<DomainGeometry, DensityProfileCoefficients>::writeToVTK(const std::filesystem::path& file_path,
                                                                      const PolarGrid& grid)
{
    const auto filename = file_path.stem().string() + ".vtu";

    std::ofstream file(file_path.parent_path() / filename);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file '" + (file_path.parent_path() / filename).string() + "'");
    }

    file << "<?xml version=\"1.0\"?>\n"
         << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
         << "<UnstructuredGrid>\n"
         << "<Piece NumberOfPoints=\"" << grid.numberOfNodes() << "\" NumberOfCells=\""
         << (grid.nr() - 1) * grid.ntheta() << "\">\n";

    // Write points
    file << "<Points>\n"
         << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    int i_r, i_theta;
    double r, theta;
    for (int index = 0; index < grid.numberOfNodes(); index++) {
        grid.multiIndex(index, i_r, i_theta);
        r     = grid.radius(i_r);
        theta = grid.theta(i_theta);
        file << domain_geometry_.Fx(r, theta) << " " << domain_geometry_.Fy(r, theta) << " " << 0 << "\n";
    }
    file << "</DataArray>\n"
         << "</Points>\n";

    // Write cells
    file << "<Cells>\n";
    file << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (int i_r = 0; i_r < grid.nr() - 1; i_r++) {
        for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
            file << grid.index(i_r, i_theta) << " " << grid.index(i_r + 1, i_theta) << " "
                 << grid.index(i_r + 1, i_theta + 1) << " " << grid.index(i_r, i_theta + 1) << "\n";
        }
    }
    file << "</DataArray>\n";

    file << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (int i = 0; i < (grid.nr() - 1) * grid.ntheta(); i++) {
        file << 4 * (i + 1) << " ";
    }
    file << "</DataArray>\n";

    file << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (int i = 0; i < (grid.nr() - 1) * grid.ntheta(); i++) {
        file << "9 "; // VTK_QUAD
    }
    file << "</DataArray>\n";

    file << "</Cells>\n"
         << "</Piece>\n"
         << "</UnstructuredGrid>\n"
         << "</VTKFile>\n";
}

template <concepts::DomainGeometry DomainGeometry, concepts::DensityProfileCoefficients DensityProfileCoefficients>
void GMGPolar<DomainGeometry, DensityProfileCoefficients>::writeToVTK(const std::filesystem::path& file_path,
                                                                      const Level<DomainGeometry>& level,
                                                                      ConstVector<double> grid_function)
{
    const PolarGrid& grid                         = level.grid();
    const LevelCache<DomainGeometry>& level_cache = level.levelCache();

    assert(std::ssize(grid_function) == grid.numberOfNodes());

    const auto filename = file_path.stem().string() + ".vtu";

    std::ofstream file(file_path.parent_path() / filename);
    if (!file.is_open())
        throw std::runtime_error("Failed to open file '" + (file_path.parent_path() / filename).string() + "'");

    file << "<?xml version=\"1.0\"?>\n"
         << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
         << "<UnstructuredGrid>\n"
         << "<Piece NumberOfPoints=\"" << grid.numberOfNodes() << "\" NumberOfCells=\""
         << (grid.nr() - 1) * grid.ntheta() << "\">\n";

    // Write points
    file << "<Points>\n"
         << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    int i_r, i_theta;
    double r, theta;
    for (int index = 0; index < grid.numberOfNodes(); index++) {
        grid.multiIndex(index, i_r, i_theta);
        r     = grid.radius(i_r);
        theta = grid.theta(i_theta);
        file << domain_geometry_.Fx(r, theta) << " " << domain_geometry_.Fy(r, theta) << " " << 0 << "\n";
    }

    file << "</DataArray>\n"
         << "</Points>\n";
    file << "<Cells>\n";

    // Connectivity
    file << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (int i_r = 0; i_r < grid.nr() - 1; i_r++) {
        for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
            file << grid.index(i_r, i_theta) << " " << grid.index(i_r + 1, i_theta) << " "
                 << grid.index(i_r + 1, i_theta + 1) << " " << grid.index(i_r, i_theta + 1) << "\n";
        }
    }

    file << "</DataArray>\n";

    file << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (int i = 0; i < (grid.nr() - 1) * grid.ntheta(); i++) {
        file << 4 * (i + 1) << " ";
    }

    file << "</DataArray>\n";
    // Cell types
    file << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (int i = 0; i < (grid.nr() - 1) * grid.ntheta(); ++i) {
        file << "9 "; // VTK_QUAD
    }

    file << "</DataArray>\n";
    file << "</Cells>\n";
    file << "<PointData Scalars=\"FunctionValues\">\n";

    file << "<DataArray type=\"Float32\" Name=\"FunctionValues\" format=\"ascii\">";

    for (int i = 0; i < grid.numberOfNodes(); i++) {
        file << grid_function[i] << "\n";
    }

    file << "</DataArray>\n"
         << "</PointData>\n"
         << "</Piece>\n"
         << "</UnstructuredGrid>\n"
         << "</VTKFile>\n";
}
