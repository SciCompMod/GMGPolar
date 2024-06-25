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
