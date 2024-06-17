#include "../../include/Smoother/smoother.h"

Smoother::Smoother(const PolarGrid& grid, const LevelCache& level_data, 
    const DomainGeometry& domain_geometry, const SystemParameters& system_parameters, const bool DirBC_Interior, 
    const int maxOpenMPThreads, const int openMPTaskThreads
) :
    grid_(grid), 
    sin_theta_(level_data.sin_theta()),
    cos_theta_(level_data.cos_theta()),
    domain_geometry_(domain_geometry),
    system_parameters_(system_parameters),
    DirBC_Interior_(DirBC_Interior),
    maxOpenMPThreads_(maxOpenMPThreads),
    openMPTaskThreads_(openMPTaskThreads),
    rhs_(grid.number_of_nodes())
{
    build_Asc_matrices(circle_Asc_matrix_, radial_Asc_matrix_);
    initializeMumps(circle_Asc_mumps_, circle_Asc_matrix_);
    initializeMumps(radial_Asc_mumps_, radial_Asc_matrix_);
}

Smoother::~Smoother(){
    deleteMumps(circle_Asc_mumps_);
    deleteMumps(radial_Asc_mumps_);
}

