#include "../../include/GMGPolar/gmgpolar.h"  

void GMGPolar::writeToVTK(const std::filesystem::path& file_path, const PolarGrid& grid) {
    const auto filename = file_path.stem().string() + ".vtu";

    std::ofstream file(file_path.parent_path() / filename);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file '" + (file_path.parent_path() / filename).string() + "'");
    }

    file << "<?xml version=\"1.0\"?>\n"
         << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
         << "<UnstructuredGrid>\n"
         << "<Piece NumberOfPoints=\"" << grid.numberOfNodes() << "\" NumberOfCells=\"" << (grid.nr() - 1) * grid.ntheta() << "\">\n";

    // Write points
    file << "<Points>\n"
         << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    int i_r, i_theta;
    double r, theta;
    for (int index = 0; index < grid.numberOfNodes(); index++) {
        grid.multiIndex(index, i_r, i_theta);
        r = grid.radius(i_r);
        theta = grid.theta(i_theta);
        double sin_theta = std::sin(theta);
        double cos_theta = std::cos(theta);
        file << domain_geometry_->Fx(r, theta, sin_theta, cos_theta) << " "
             << domain_geometry_->Fy(r, theta, sin_theta, cos_theta) << " "
             << 0 << "\n";
    }
    file << "</DataArray>\n"
         << "</Points>\n";

    // Write cells
    file << "<Cells>\n";
    file << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (int i_r = 0; i_r < grid.nr() - 1; i_r++) {
        for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
            file << grid.index(i_r, i_theta) << " "
                 << grid.index(i_r + 1, i_theta) << " "
                 << grid.index(i_r + 1, i_theta + 1) << " "
                 << grid.index(i_r, i_theta + 1) << "\n";
        }
    }
    file << "</DataArray>\n";

    file << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (size_t i = 0; i < (grid.nr() - 1) * grid.ntheta(); i++) {
        file << 4 * (i + 1) << " ";
    }
    file << "</DataArray>\n";

    file << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (size_t i = 0; i < (grid.nr() - 1) * grid.ntheta(); i++) {
        file << "9 "; // VTK_QUAD
    }
    file << "</DataArray>\n";

    file << "</Cells>\n"
         << "</Piece>\n"
         << "</UnstructuredGrid>\n"
         << "</VTKFile>\n";
}



void GMGPolar::writeToVTK(const std::filesystem::path& file_path, const Level& level, const Vector<double>& grid_function){
    const PolarGrid& grid = level.grid();
    const LevelCache& level_cache = level.levelCache();

    assert(grid.numberOfNodes() == grid_function.size());

    const Vector<double>& sin_theta_cache = level_cache.sin_theta();
    const Vector<double>& cos_theta_cache = level_cache.cos_theta();

    const auto filename = file_path.stem().string() + ".vtu";

    std::ofstream file(file_path.parent_path() / filename);
    if(!file.is_open()) throw std::runtime_error("Failed to open file '" + (file_path.parent_path() / filename).string() + "'");

   file << "<?xml version=\"1.0\"?>\n"
        << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
        << "<UnstructuredGrid>\n"
        << "<Piece NumberOfPoints=\"" << grid.numberOfNodes() << "\" NumberOfCells=\"" << (grid.nr()-1)*grid.ntheta() << "\">\n";

    // Write points
    file << "<Points>\n"
        << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    int i_r, i_theta;
    double r, theta;
    double sin_theta, cos_theta;
    for (int index = 0; index < grid.numberOfNodes(); index++) {
        grid.multiIndex(index, i_r, i_theta);
        r = grid.radius(i_r);
        theta = grid.theta(i_theta);
        sin_theta = sin_theta_cache[i_theta];
        cos_theta = cos_theta_cache[i_theta];
        file << domain_geometry_->Fx(r, theta, sin_theta, cos_theta) << " " << domain_geometry_->Fy(r, theta, sin_theta, cos_theta) << " " << 0 << "\n";
    }

    file << "</DataArray>\n"
        << "</Points>\n";
    file << "<Cells>\n";

    // Connectivity
    file << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (int i_r = 0; i_r < grid.nr() - 1; i_r++) {
        for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
            file << grid.index(i_r, i_theta) << " " << grid.index(i_r + 1, i_theta) << " " << grid.index(i_r + 1, i_theta + 1) << " " << grid.index(i_r, i_theta + 1) << "\n";
        }
    }

    file <<"</DataArray>\n";

    
    file << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (size_t i = 0; i < (grid.nr() - 1) * grid.ntheta(); i++) {
        file << 4 * (i + 1) << " ";
    }

    file << "</DataArray>\n";
    // Cell types
    file << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (size_t i = 0; i < (grid.nr() - 1) * grid.ntheta(); ++i) {
        file << "9 "; // VTK_QUAD
    }

    file << "</DataArray>\n";
    file << "</Cells>\n";
    file << "<PointData Scalars=\"FunctionValues\">\n";

    file << "<DataArray type=\"Float32\" Name=\"FunctionValues\" format=\"ascii\">";


    for (int i = 0; i < grid.numberOfNodes(); i++) {
        file << grid_function[i] << "\n";
    }

    file <<"</DataArray>\n"
        <<"</PointData>\n"
        << "</Piece>\n"
        << "</UnstructuredGrid>\n"
        << "</VTKFile>\n";
}

