#include <gtest/gtest.h>
#include <random>

#include "../../include/GMGPolar/gmgpolar.h"

#include "../../include/Smoother/SmootherTakeCPU/smoother.h"
#include "../../include/Smoother/SmootherTakeGPU/smoother.h"

#include <chrono>

namespace SmootherTest
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

    TEST(SmootherTest, applySmoother_DirBC_Interior)
    {
        std::vector<double> radii  = {0.1, 0.2, 0.25, 0.5, 0.8, 0.9, 1.3, 1.4, 2.0};
        std::vector<double> angles = {
            0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

        double Rmax = radii.back();

        // Parameters
        ExtrapolationType extrapolation = ExtrapolationType::NONE;
        int FMG = 0;
        bool DirBC_Interior = false;
        int num_omp_threads = 1;
        omp_set_num_threads(num_omp_threads);

        Vector<double> result;
        Vector<double> host_result;

        // CPU Processing
        {
            ProcessingType processing_type = ProcessingType::CPU;
            DomainGeometry domain_geometry;
            DensityProfileCoefficients density_profile_coefficients;

          auto grid = std::make_unique<PolarGrid>(radii, angles);

    //     const double R0 = 0.1;
    //     const double Rmax = 1.3;
    //     const int nr_exp = 3;
    //     const int ntheta_exp = 3;
    //     const double refinement_radius = 0.5;
    //     const int anisotropic_factor = 0;
    //     const int divideBy2 = 9;

    //    auto grid = std::make_unique<PolarGrid>(R0, Rmax, nr_exp, ntheta_exp, refinement_radius, anisotropic_factor, divideBy2, 0.6);

        std::cout<<grid->ntheta()<<std::endl;
         std::cout<<grid->numberSmootherCircles()<<std::endl;

            auto CPU_levelCache = std::make_unique<LevelCache>(processing_type, *grid, density_profile_coefficients, domain_geometry);

            Level level(0, processing_type, std::move(grid), std::move(CPU_levelCache), extrapolation, FMG);

            // Generate random data using the level's grid
            Vector<double> x = generate_random_sample_data(level.grid(), 42);
            Vector<double> rhs = generate_random_sample_data(level.grid(), 69);
            Vector<double> temp = generate_random_sample_data(level.grid(), 187);
            assign(x, 1.0);
            assign(rhs, 1.0);
            assign(temp, 1.0);

            SmootherTakeCPU op_smootherTakeCPU(level, domain_geometry, density_profile_coefficients, DirBC_Interior);

            auto start_time = std::chrono::high_resolution_clock::now();

            op_smootherTakeCPU.smoothingInPlace(x, rhs, temp);

            auto end_time = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end_time - start_time;
            
            std::cout << "CPU time: " << elapsed.count() << " seconds" << std::endl;


            result = x;
        }

        // GPU Processing
        {
            ProcessingType processing_type = ProcessingType::GPU;
            DomainGeometry domain_geometry;
            DensityProfileCoefficients density_profile_coefficients;

       auto grid = std::make_unique<PolarGrid>(radii, angles);

        // const double R0 = 0.1;
        // const double Rmax = 1.3;
        // const int nr_exp = 3;
        // const int ntheta_exp = 3;
        // const double refinement_radius = 0.5;
        // const int anisotropic_factor = 0;
        // const int divideBy2 = 9;

        // auto grid = std::make_unique<PolarGrid>(R0, Rmax, nr_exp, ntheta_exp, refinement_radius, anisotropic_factor, divideBy2, 0.6);



            auto GPU_levelCache = std::make_unique<LevelCache>(processing_type, *grid, density_profile_coefficients, domain_geometry);

            Level level(0, processing_type, std::move(grid), std::move(GPU_levelCache), extrapolation, FMG);

            // Generate random data using the level's grid (GPU-compatible)
            Vector<double> x = generate_random_sample_data(level.grid(), 42);
            Vector<double> rhs = generate_random_sample_data(level.grid(), 69);
            Vector<double> temp = generate_random_sample_data(level.grid(), 187);

            assign(x, 1.0);
            assign(rhs, 1.0);
            assign(temp, 1.0);

            GPU_Vector<double> gpu_x(x.size());
            GPU_Vector<double> gpu_rhs(rhs.size());
            GPU_Vector<double> gpu_temp(temp.size());
            copyHostToDevice(x, gpu_x);
            copyHostToDevice(rhs, gpu_rhs);
            copyHostToDevice(temp, gpu_temp);

            SmootherTakeGPU op_smootherTakeGPU(level, domain_geometry, density_profile_coefficients, DirBC_Interior);

            auto start_time = std::chrono::high_resolution_clock::now();

            op_smootherTakeGPU.smoothingInPlace(gpu_x, gpu_rhs, gpu_temp);

            auto end_time = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end_time - start_time;
            
            std::cout << "GPU time: " << elapsed.count() << " seconds" << std::endl;

            host_result = Vector<double>(gpu_x.size());
            copyDeviceToHost(gpu_x, host_result);
        }

        // std::cout<<result<<std::endl;

        // std::cout<<host_result<<std::endl;
        

        // ASSERT_EQ(result.size(), host_result.size());

        // double tolerance = 1e-10;
        // for (size_t i = 0; i < result.size(); ++i) {
        //     ASSERT_NEAR(result[i], host_result[i], tolerance);
        // }
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}