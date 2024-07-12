#include "../../include/GMGPolar/gmgpolar.h"

GMGPolar::GMGPolar() : 
    parser_()
{
    initializeGrid(); initializeGeometry();
    initializeMultigrid(); initializeGeneral();
    
    setParameters(0, nullptr);
}


void GMGPolar::setParameters(int argc, char* argv[]) {
    if(argc != 0){
        try {
            parser_.parse_check(argc, argv);
        } catch (const cmdline::cmdline_error& parse_error) {
            std::cerr << "Error: " << parse_error.what() << std::endl;
            std::cerr << "Usage: " << parser_.usage() << std::endl;
        }
    }
    
    parseGrid(); parseGeometry();
    parseMultigrid(); parseGeneral();
}


void GMGPolar::setSolution(std::unique_ptr<const ExactSolution> exact_solution) {
    exact_solution_ = std::move(exact_solution);
}


/* ----------------- */
/* GMGPolar Solution */
Vector<double>& GMGPolar::solution() {
    return levels_[0].solution();
}
const Vector<double>& GMGPolar::solution() const {
    return levels_[0].solution();
}

const PolarGrid& GMGPolar::grid() const {
    return levels_[0].grid();
}

void GMGPolar::printTimings() const {
    std::cout << "\nTiming Information:" << std::endl;
    std::cout << "Setup Time: " << t_setup_total << " seconds" << std::endl;
    std::cout << "    Create Levels: " << t_setup_createLevels << " seconds" << std::endl;
    std::cout << "    (Build rhs_f: " << t_setup_rhs << " seconds)" << std::endl;
    std::cout << "    Direct Solver: " << t_setup_directSolver << " seconds" << std::endl;
    std::cout << "    Smoother: " << t_setup_smoother << " seconds" << std::endl;
    std::cout << "\nSolve Time: " << t_solve_total << " seconds" << std::endl;
    std::cout << "    Multigrid Iteration: " << t_solve_multigrid_iterations << " seconds" << std::endl;
    std::cout << "    Check Convergence: " << t_check_convergence << " seconds" << std::endl;
    std::cout << "    Check Exact Error: " << t_check_exact_error << " seconds" << std::endl;
    std::cout << "\nAverage Multigrid Iteration: " << t_avg_MGC_total << " seconds" << std::endl;
    std::cout << "    PreSmoothing: " << t_avg_MGC_preSmoothing << " seconds" << std::endl;
    std::cout << "    PostSmoothing: " << t_avg_MGC_postSmoothing << " seconds" << std::endl;
    std::cout << "    Residual: " << t_avg_MGC_residual << " seconds" << std::endl;
    std::cout << "    DirectSolver: " << t_avg_MGC_directSolver << " seconds" << std::endl;
    std::cout<<"\n"<<std::endl;
}


/* --------------- */
/* Grid Parameters */
double GMGPolar::R0() const {
    return R0_;
}

void GMGPolar::R0(double R0) {
    R0_ = R0;
}

double GMGPolar::Rmax() const {
    return Rmax_;
}

void GMGPolar::Rmax(double Rmax) {
    Rmax_ = Rmax;
}

int GMGPolar::nr_exp() const {
    return nr_exp_;
}

void GMGPolar::nr_exp(int nr_exp) {
    nr_exp_ = nr_exp;
}

int GMGPolar::ntheta_exp() const {
    return ntheta_exp_;
}

void GMGPolar::ntheta_exp(int ntheta_exp) {
    ntheta_exp_ = ntheta_exp;
}

int GMGPolar::anisotropic_factor() const {
    return anisotropic_factor_;
}

void GMGPolar::anisotropic_factor(int anisotropic_factor) {
    anisotropic_factor_ = anisotropic_factor;
}

int GMGPolar::divideBy2() const {
    return divideBy2_;
}

void GMGPolar::divideBy2(int divideBy2) {
    divideBy2_ = divideBy2;
}

bool GMGPolar::write_grid_file() const {
    return write_grid_file_;
}

void GMGPolar::write_grid_file(bool write_grid_file) {
    write_grid_file_ = write_grid_file;
}

bool GMGPolar::load_grid_file() const {
    return load_grid_file_;
}

void GMGPolar::load_grid_file(bool load_grid_file) {
    load_grid_file_ = load_grid_file;
}

std::string GMGPolar::file_grid_radii() const {
    return file_grid_radii_;
}

void GMGPolar::file_grid_radii(const std::string& file_grid_radii) {
    file_grid_radii_ = file_grid_radii;
}

std::string GMGPolar::file_grid_angles() const {
    return file_grid_angles_;
}

void GMGPolar::file_grid_angles(const std::string& file_grid_angles) {
    file_grid_angles_ = file_grid_angles;
}

/* ------------------- */
/* Geometry Parameters */
bool GMGPolar::DirBC_Interior() const {
    return DirBC_Interior_;
}

void GMGPolar::DirBC_Interior(bool DirBC_Interior) {
    DirBC_Interior_ = DirBC_Interior;
}

/* -------------------- */
/* Multigrid Parameters */
int GMGPolar::extrapolation() const {
    return extrapolation_;
}

void GMGPolar::extrapolation(int extrapolation) {
    extrapolation_ = extrapolation;
}

int GMGPolar::maxLevels() const {
    return maxLevels_;
}

void GMGPolar::maxLevels(int maxLevels) {
    maxLevels_ = maxLevels;
}

MultigridCycleType GMGPolar::multigrid_cycle() const {
    return multigrid_cycle_;
}

void GMGPolar::multigrid_cycle(MultigridCycleType multigrid_cycle) {
    multigrid_cycle_ = multigrid_cycle;
}

int GMGPolar::preSmoothingSteps() const {
    return preSmoothingSteps_;
}

void GMGPolar::preSmoothingSteps(int preSmoothingSteps) {
    preSmoothingSteps_ = preSmoothingSteps;
}

int GMGPolar::postSmoothingSteps() const {
    return postSmoothingSteps_;
}

void GMGPolar::postSmoothingSteps(int postSmoothingSteps) {
    postSmoothingSteps_ = postSmoothingSteps;
}

int GMGPolar::maxIterations() const {
    return max_iterations_;
}

void GMGPolar::maxIterations(int maxIterations) {
    max_iterations_ = maxIterations;
}

double GMGPolar::absoluteTolerance() const {
    if (absolute_tolerance_.has_value()) {
        return absolute_tolerance_.value();
    } else {
        return -1.0;
    }
}

void GMGPolar::absoluteTolerance(double absoluteTolerance) {
    if(absoluteTolerance > 0 ){
        absolute_tolerance_ = absoluteTolerance;
    } else{
        absolute_tolerance_ = std::nullopt;
    }
}

double GMGPolar::relativeTolerance() const {
    if (relative_tolerance_.has_value()) {
        return relative_tolerance_.value();
    } else {
        return -1.0;
    }
}

void GMGPolar::relativeTolerance(double relativeTolerance) {
    if(relativeTolerance > 0 ){
        relative_tolerance_ = relativeTolerance;
    } else{
        relative_tolerance_ = std::nullopt;
    }
}

/* ------------------ */
/* Control Parameters */
int GMGPolar::maxOpenMPThreads() const {
    return maxOpenMPThreads_;
}

void GMGPolar::maxOpenMPThreads(int maxOpenMPThreads) {
    maxOpenMPThreads_ = maxOpenMPThreads;
}

int GMGPolar::finestLevelThreads() const {
    return finestLevelThreads_;
}

void GMGPolar::finestLevelThreads(int finestLevelThreads) {
    finestLevelThreads_ = finestLevelThreads;
}

double GMGPolar::threadReductionFactor() const {
    return threadReductionFactor_;
}

void GMGPolar::threadReductionFactor(double threadReductionFactor) {
    threadReductionFactor_ = threadReductionFactor;
}
