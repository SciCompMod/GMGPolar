#include <gtest/gtest.h>
#include <random>

#include "../../include/GMGPolar/gmgpolar.h"

#include "../../include/Residual/ResidualTakeCPU/residual.h"
#include "../../include/Residual/ResidualTakeGPU/residual.h"

namespace ResidualTest
{

    Vector<double> generate_random_sample_data(const PolarGrid& grid, unsigned int seed)
    {
        Vector<double> x(grid.numberOfNodes());
        std::mt19937 gen(seed);
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        for (int i = 0; i < x.size(); ++i) {
            x[i] = dist(gen);
        }
        return x;
    }

    TEST(ResidualTest, applyResidual_DirBC_Interior)
    {
        std::vector<double> radii  = {0.1, 0.2, 0.25, 0.5, 0.8, 0.9, 1.3, 1.4, 2.0};
        std::vector<double> angles = {
            0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

        double Rmax = radii.back();

        // Parameters
        ExtrapolationType extrapolation = ExtrapolationType::NONE;
        int FMG = 0;
        bool DirBC_Interior = true;
        int num_omp_threads = 16;
        omp_set_num_threads(num_omp_threads);

        Vector<double> result;
        Vector<double> host_result;

        // CPU Processing
        {
            ProcessingType processing_type = ProcessingType::CPU;
            DomainGeometry domain_geometry;
            DensityProfileCoefficients density_profile_coefficients;

            auto grid = std::make_unique<PolarGrid>(radii, angles);
            auto CPU_levelCache = std::make_unique<LevelCache>(processing_type, *grid, density_profile_coefficients, domain_geometry);

            Level level(0, processing_type, std::move(grid), std::move(CPU_levelCache), extrapolation, FMG);

            // Generate random data using the level's grid
            Vector<double> x   = generate_random_sample_data(level.grid(), 42);
            Vector<double> rhs = generate_random_sample_data(level.grid(), 69);

            ResidualTakeCPU op_residualTakeCPU(level, domain_geometry, density_profile_coefficients, DirBC_Interior);

            result = Vector<double> (level.grid().numberOfNodes());
            op_residualTakeCPU.computeResidual(result, rhs, x);
        }

        // GPU Processing
        {
            ProcessingType processing_type = ProcessingType::GPU;
            DomainGeometry domain_geometry;
            DensityProfileCoefficients density_profile_coefficients;

            auto grid = std::make_unique<PolarGrid>(radii, angles);
            auto GPU_levelCache = std::make_unique<LevelCache>(processing_type, *grid, density_profile_coefficients, domain_geometry);

            Level level(0, processing_type, std::move(grid), std::move(GPU_levelCache), extrapolation, FMG);

            // Generate random data using the level's grid (GPU-compatible)
            Vector<double> x   = generate_random_sample_data(level.grid(), 42);
            Vector<double> rhs = generate_random_sample_data(level.grid(), 69);

            GPU_Vector<double> gpu_x(x.size());
            GPU_Vector<double> gpu_rhs(rhs.size());
            copyHostToDevice(x, gpu_x);
            copyHostToDevice(rhs, gpu_rhs);

            GPU_Vector<double> gpu_result(level.grid().numberOfNodes());
            ResidualTakeGPU op_residualTakeGPU(level, domain_geometry, density_profile_coefficients, DirBC_Interior);

            op_residualTakeGPU.computeResidual(gpu_result, gpu_rhs, gpu_x);

            host_result = Vector<double>(gpu_result.size());
            copyDeviceToHost(gpu_result, host_result);
        }

        ASSERT_EQ(result.size(), host_result.size());

        double tolerance = 1e-10;
        for (size_t i = 0; i < result.size(); ++i) {
            ASSERT_NEAR(result[i], host_result[i], tolerance);
        }
    }

    TEST(ResidualTest, applyResidualAcrossOrigin)
    {
        std::vector<double> radii  = {0.1, 0.2, 0.25, 0.5, 0.8, 0.9, 1.3, 1.4, 2.0};
        std::vector<double> angles = {
            0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

        double Rmax = radii.back();

        // Parameters
        ExtrapolationType extrapolation = ExtrapolationType::NONE;
        int FMG = 0;
        bool DirBC_Interior = false;
        int num_omp_threads = 16;
        omp_set_num_threads(num_omp_threads);

        Vector<double> result;
        Vector<double> host_result;

        // CPU Processing
        {
            ProcessingType processing_type = ProcessingType::CPU;
            DomainGeometry domain_geometry;
            DensityProfileCoefficients density_profile_coefficients;

            auto grid = std::make_unique<PolarGrid>(radii, angles);
            auto CPU_levelCache = std::make_unique<LevelCache>(processing_type, *grid, density_profile_coefficients, domain_geometry);

            Level level(0, processing_type, std::move(grid), std::move(CPU_levelCache), extrapolation, FMG);

            // Generate random data using the level's grid
            Vector<double> x   = generate_random_sample_data(level.grid(), 42);
            Vector<double> rhs = generate_random_sample_data(level.grid(), 69);

            ResidualTakeCPU op_residualTakeCPU(level, domain_geometry, density_profile_coefficients, DirBC_Interior);

            result = Vector<double> (level.grid().numberOfNodes());
            op_residualTakeCPU.computeResidual(result, rhs, x);
        }

        // GPU Processing
        {
            ProcessingType processing_type = ProcessingType::GPU;
            DomainGeometry domain_geometry;
            DensityProfileCoefficients density_profile_coefficients;

            auto grid = std::make_unique<PolarGrid>(radii, angles);
            auto GPU_levelCache = std::make_unique<LevelCache>(processing_type, *grid, density_profile_coefficients, domain_geometry);

            Level level(0, processing_type, std::move(grid), std::move(GPU_levelCache), extrapolation, FMG);

            // Generate random data using the level's grid (GPU-compatible)
            Vector<double> x   = generate_random_sample_data(level.grid(), 42);
            Vector<double> rhs = generate_random_sample_data(level.grid(), 69);

            GPU_Vector<double> gpu_x(x.size());
            GPU_Vector<double> gpu_rhs(rhs.size());
            copyHostToDevice(x, gpu_x);
            copyHostToDevice(rhs, gpu_rhs);

            GPU_Vector<double> gpu_result(level.grid().numberOfNodes());
            ResidualTakeGPU op_residualTakeGPU(level, domain_geometry, density_profile_coefficients, DirBC_Interior);

            op_residualTakeGPU.computeResidual(gpu_result, gpu_rhs, gpu_x);

            host_result = Vector<double>(gpu_result.size());
            copyDeviceToHost(gpu_result, host_result);
        }

        ASSERT_EQ(result.size(), host_result.size());

        double tolerance = 1e-10;
        for (size_t i = 0; i < result.size(); ++i) {
            ASSERT_NEAR(result[i], host_result[i], tolerance);
        }
    }

}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}