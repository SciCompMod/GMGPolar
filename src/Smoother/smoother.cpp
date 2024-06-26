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
    openMPTaskThreads_(openMPTaskThreads)
{
    build_Asc_matrices();
    initializeMumps(inner_boundary_circle_Asc_mumps_, inner_boundary_circle_Asc_matrix_);    
}

Smoother::~Smoother(){
    deleteMumps(inner_boundary_circle_Asc_mumps_);
}

