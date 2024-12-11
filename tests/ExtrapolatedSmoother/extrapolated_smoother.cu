#include <gtest/gtest.h>
#include <random>

#include "../../include/GMGPolar/gmgpolar.h"

#include "../../include/ExtrapolatedSmoother/ExtrapolatedSmootherTakeCPU/extrapolatedSmoother.h"
#include "../../include/ExtrapolatedSmoother/ExtrapolatedSmootherTakeGPU/extrapolatedSmoother.h"

#include <chrono>

namespace ExtrapolatedSmootherTest
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

    TEST(ExtrapolatedSmootherTest, applyExtrapolatedSmoother_DirBC_Interior)
    {
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

            const double R0 = 0.1;
            const double Rmax = 1.3;
            const int nr_exp = 4;
            const int ntheta_exp = -1;
            const double refinement_radius = 0.5;
            const int anisotropic_factor = 2;
            const int divideBy2 = 4;

            auto grid = std::make_unique<PolarGrid>(R0, Rmax, nr_exp, ntheta_exp, refinement_radius, anisotropic_factor, divideBy2, 0.5);

            auto CPU_levelCache = std::make_unique<LevelCache>(processing_type, *grid, density_profile_coefficients, domain_geometry);

            Level level(0, processing_type, std::move(grid), std::move(CPU_levelCache), extrapolation, FMG);

            // Generate random data using the level's grid
            Vector<double> x = generate_random_sample_data(level.grid(), 42);
            Vector<double> rhs = generate_random_sample_data(level.grid(), 69);
            Vector<double> temp = generate_random_sample_data(level.grid(), 187);

            ExtrapolatedSmootherTakeCPU op_extrapolatedSmootherTakeCPU(level, domain_geometry, density_profile_coefficients, DirBC_Interior);

            op_extrapolatedSmootherTakeCPU.extrapolatedSmoothingInPlace(x, rhs, temp);

            result = x;
        }

        // GPU Processing
        {
            ProcessingType processing_type = ProcessingType::GPU;
            DomainGeometry domain_geometry;
            DensityProfileCoefficients density_profile_coefficients;

            const double R0 = 0.1;
            const double Rmax = 1.3;
            const int nr_exp = 4;
            const int ntheta_exp = -1;
            const double refinement_radius = 0.5;
            const int anisotropic_factor = 2;
            const int divideBy2 = 4;

            auto grid = std::make_unique<PolarGrid>(R0, Rmax, nr_exp, ntheta_exp, refinement_radius, anisotropic_factor, divideBy2, 0.5);

            auto GPU_levelCache = std::make_unique<LevelCache>(processing_type, *grid, density_profile_coefficients, domain_geometry);

            Level level(0, processing_type, std::move(grid), std::move(GPU_levelCache), extrapolation, FMG);

            // Generate random data using the level's grid (GPU-compatible)
            Vector<double> x = generate_random_sample_data(level.grid(), 42);
            Vector<double> rhs = generate_random_sample_data(level.grid(), 69);
            Vector<double> temp = generate_random_sample_data(level.grid(), 187);

            GPU_Vector<double> gpu_x(x.size());
            GPU_Vector<double> gpu_rhs(rhs.size());
            GPU_Vector<double> gpu_temp(temp.size());
            copyHostToDevice(x, gpu_x);
            copyHostToDevice(rhs, gpu_rhs);
            copyHostToDevice(temp, gpu_temp);

            ExtrapolatedSmootherTakeGPU op_extrapolatedSmootherTakeGPU(level, domain_geometry, density_profile_coefficients, DirBC_Interior);

            op_extrapolatedSmootherTakeGPU.extrapolatedSmoothingInPlace(gpu_x, gpu_rhs, gpu_temp);

            host_result = Vector<double>(gpu_x.size());
            copyDeviceToHost(gpu_x, host_result);
        }

        ASSERT_EQ(result.size(), host_result.size());

        double tolerance = 1e-8;
        for (size_t i = 0; i < result.size(); ++i) {
            ASSERT_NEAR(result[i], host_result[i], tolerance);
        }
    }

    TEST(ExtrapolatedSmootherTest, applyExtrapolatedSmoother_AcrossOrigin)
    {
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

            const double R0 = 0.1;
            const double Rmax = 1.3;
            const int nr_exp = 4;
            const int ntheta_exp = -1;
            const double refinement_radius = 0.5;
            const int anisotropic_factor = 2;
            const int divideBy2 = 4;

            auto grid = std::make_unique<PolarGrid>(R0, Rmax, nr_exp, ntheta_exp, refinement_radius, anisotropic_factor, divideBy2, 0.5);

            auto CPU_levelCache = std::make_unique<LevelCache>(processing_type, *grid, density_profile_coefficients, domain_geometry);

            Level level(0, processing_type, std::move(grid), std::move(CPU_levelCache), extrapolation, FMG);

            // Generate random data using the level's grid
            Vector<double> x = generate_random_sample_data(level.grid(), 42);
            Vector<double> rhs = generate_random_sample_data(level.grid(), 69);
            Vector<double> temp = generate_random_sample_data(level.grid(), 187);

            ExtrapolatedSmootherTakeCPU op_extrapolatedSmootherTakeCPU(level, domain_geometry, density_profile_coefficients, DirBC_Interior);

            op_extrapolatedSmootherTakeCPU.extrapolatedSmoothingInPlace(x, rhs, temp);

            result = x;
        }

        // GPU Processing
        {
            ProcessingType processing_type = ProcessingType::GPU;
            DomainGeometry domain_geometry;
            DensityProfileCoefficients density_profile_coefficients;

            const double R0 = 0.1;
            const double Rmax = 1.3;
            const int nr_exp = 4;
            const int ntheta_exp = -1;
            const double refinement_radius = 0.5;
            const int anisotropic_factor = 2;
            const int divideBy2 = 4;

            auto grid = std::make_unique<PolarGrid>(R0, Rmax, nr_exp, ntheta_exp, refinement_radius, anisotropic_factor, divideBy2, 0.5);

            auto GPU_levelCache = std::make_unique<LevelCache>(processing_type, *grid, density_profile_coefficients, domain_geometry);

            Level level(0, processing_type, std::move(grid), std::move(GPU_levelCache), extrapolation, FMG);

            // Generate random data using the level's grid (GPU-compatible)
            Vector<double> x = generate_random_sample_data(level.grid(), 42);
            Vector<double> rhs = generate_random_sample_data(level.grid(), 69);
            Vector<double> temp = generate_random_sample_data(level.grid(), 187);

            GPU_Vector<double> gpu_x(x.size());
            GPU_Vector<double> gpu_rhs(rhs.size());
            GPU_Vector<double> gpu_temp(temp.size());
            copyHostToDevice(x, gpu_x);
            copyHostToDevice(rhs, gpu_rhs);
            copyHostToDevice(temp, gpu_temp);

            ExtrapolatedSmootherTakeGPU op_extrapolatedSmootherTakeGPU(level, domain_geometry, density_profile_coefficients, DirBC_Interior);

            op_extrapolatedSmootherTakeGPU.extrapolatedSmoothingInPlace(gpu_x, gpu_rhs, gpu_temp);

            host_result = Vector<double>(gpu_x.size());
            copyDeviceToHost(gpu_x, host_result);
        }

        ASSERT_EQ(result.size(), host_result.size());

        double tolerance = 1e-8;
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