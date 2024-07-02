#include "../../include/GMGPolar/gmgpolar.h"

/* Dies ist ein Test */

GMGPolar::GMGPolar(const DomainGeometry& domain_geometry, const SystemParameters& system_parameters) : 
    domain_geometry_(domain_geometry),
    system_parameters_(system_parameters),
    parser_()
{
    initializeGrid(); initializeGeometry();
    initializeMultigrid(); initializeGeneral();
    setParameters(0, nullptr);
}


void GMGPolar::setSolution(const ExactSolution& exact_solution) {
    exact_solution_ = std::make_shared<ExactSolution>(exact_solution);
}

void GMGPolar::setParameters(int argc, char* argv[]) {
    if(argc != 0){
        try {
            parser_.parse_check(argc, argv);
        } catch (const cmdline::cmdline_error &parse_error) {
            std::cerr << "Error: " << parse_error.what() << std::endl;
            std::cerr << "Usage: " << parser_.usage() << std::endl;
        }
    }
    
    parseGrid(); parseGeometry();
    parseMultigrid(); parseGeneral();
}


void GMGPolar::prolongateToUpperLevel(const int current_level, Vector<double>& result, const Vector<double>& x) const {
    assert(static_cast<size_t>(current_level) < levels_.size() && 1 <= current_level);
    if(!interpolation_) throw std::runtime_error("Interpolation not initialized.");
    interpolation_->applyProlongation(levels_[current_level], levels_[current_level-1], result, x);
}

void GMGPolar::restrictToLowerLevel(const int current_level, Vector<double>& result, const Vector<double>& x) const {
    assert(static_cast<size_t>(current_level) < levels_.size() - 1 && 0 <= current_level);
    if(!interpolation_) throw std::runtime_error("Interpolation not initialized.");
    interpolation_->applyRestriction(levels_[current_level], levels_[current_level+1], result, x);
}