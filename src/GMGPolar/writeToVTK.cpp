#include "../../include/GMGPolar/gmgpolar.h"

void GMGPolar::writeToVTK(const std::filesystem::path& file_path, const PolarGrid& grid)
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
        r                = grid.radius(i_r);
        theta            = grid.theta(i_theta);
        double sin_theta = std::sin(theta);
        double cos_theta = std::cos(theta);
        file << domain_geometry_->Fx(r, theta, sin_theta, cos_theta) << " "
             << domain_geometry_->Fy(r, theta, sin_theta, cos_theta) << " " << 0 << "\n";
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

void GMGPolar::writeToVTK(const std::filesystem::path& file_path, const Level& level,
                          const Vector<double>& grid_function)
{
    const PolarGrid& grid         = level.grid();
    const LevelCache& level_cache = level.levelCache();

    assert(grid.numberOfNodes() == grid_function.size());

    const Vector<double>& sin_theta_cache = level_cache.sin_theta();
    const Vector<double>& cos_theta_cache = level_cache.cos_theta();

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
    double sin_theta, cos_theta;
    for (int index = 0; index < grid.numberOfNodes(); index++) {
        grid.multiIndex(index, i_r, i_theta);
        r         = grid.radius(i_r);
        theta     = grid.theta(i_theta);
        sin_theta = sin_theta_cache[i_theta];
        cos_theta = cos_theta_cache[i_theta];
        file << domain_geometry_->Fx(r, theta, sin_theta, cos_theta) << " "
             << domain_geometry_->Fy(r, theta, sin_theta, cos_theta) << " " << 0 << "\n";
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


void GMGPolar::writeToVTKTorus(const std::filesystem::path& file_path, const Level& level,
    const Vector<double>& grid_function, int n_phi, double major_radius)
{
const PolarGrid& grid         = level.grid();
const LevelCache& level_cache = level.levelCache();

assert(grid.numberOfNodes() == grid_function.size());

const Vector<double>& sin_theta_cache = level_cache.sin_theta();
const Vector<double>& cos_theta_cache = level_cache.cos_theta();

const auto filename = file_path.stem().string() + "_torus.vtu";

std::ofstream file(file_path.parent_path() / filename);
if (!file.is_open())
throw std::runtime_error("Failed to open file '" + (file_path.parent_path() / filename).string() + "'");

const int total_points = grid.numberOfNodes() * n_phi;
const int total_cells = (grid.nr() - 1) * grid.ntheta() * n_phi;

file << "<?xml version=\"1.0\"?>\n"
<< "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
<< "<UnstructuredGrid>\n"
<< "<Piece NumberOfPoints=\"" << total_points << "\" NumberOfCells=\"" << total_cells << "\">\n";

// Write points
file << "<Points>\n"
<< "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";

for (int i_phi = 0; i_phi < n_phi; ++i_phi) {
double phi = 2.0 * M_PI * i_phi / n_phi;
double cos_phi = std::cos(phi);
double sin_phi = std::sin(phi);

for (int index = 0; index < grid.numberOfNodes(); index++) {
int i_r, i_theta;
grid.multiIndex(index, i_r, i_theta);
double r = grid.radius(i_r);
double theta = grid.theta(i_theta);
double sin_theta = sin_theta_cache[i_theta];
double cos_theta = cos_theta_cache[i_theta];

// Original 2D cross-section coordinates (X, Y)
double x2d = domain_geometry_->Fx(r, theta, sin_theta, cos_theta);
double y2d = domain_geometry_->Fy(r, theta, sin_theta, cos_theta);

// Map into torus (3D)
double X = (major_radius + x2d) * cos_phi;
double Y = (major_radius + x2d) * sin_phi;
double Z = y2d;

file << X << " " << Y << " " << Z << "\n";
}
}

file << "</DataArray>\n"
<< "</Points>\n";

// Cells
file << "<Cells>\n";

// Connectivity
file << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";

for (int i_phi = 0; i_phi < n_phi; ++i_phi) {
int next_phi = (i_phi + 1) % n_phi;

for (int i_r = 0; i_r < grid.nr() - 1; i_r++) {
for (int i_theta = 0; i_theta < grid.ntheta(); i_theta++) {
int i_theta_next = (i_theta + 1) % grid.ntheta();

// Indices in current phi-slice
int id0 = grid.index(i_r, i_theta) + i_phi * grid.numberOfNodes();
int id1 = grid.index(i_r + 1, i_theta) + i_phi * grid.numberOfNodes();
int id2 = grid.index(i_r + 1, i_theta_next) + i_phi * grid.numberOfNodes();
int id3 = grid.index(i_r, i_theta_next) + i_phi * grid.numberOfNodes();

// Same indices in next phi-slice for wrapping cells (optional 3D cells)

file << id0 << " " << id1 << " " << id2 << " " << id3 << "\n";
}
}
}

file << "</DataArray>\n";

// Offsets
file << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
for (int i = 0; i < total_cells; i++) {
file << 4 * (i + 1) << " ";
}
file << "\n</DataArray>\n";

// Cell types
file << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
for (int i = 0; i < total_cells; ++i) {
file << "9 "; // VTK_QUAD
}
file << "\n</DataArray>\n";

file << "</Cells>\n";

// Point Data
file << "<PointData Scalars=\"FunctionValues\">\n"
<< "<DataArray type=\"Float32\" Name=\"FunctionValues\" format=\"ascii\">\n";

for (int i_phi = 0; i_phi < n_phi; ++i_phi) {
for (int i = 0; i < grid.numberOfNodes(); i++) {
file << grid_function[i] << "\n";
}
}

file << "</DataArray>\n"
<< "</PointData>\n"
<< "</Piece>\n"
<< "</UnstructuredGrid>\n"
<< "</VTKFile>\n";
}
