#include <gtest/gtest.h>
#include <random>

#include "../../include/GMGPolar/gmgpolar.h"

namespace ExtrapolatedProlongationTest
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

    TEST(ExtrapolatedProlongationTest, applyExtrapolatedProlongation)
    {
        std::vector<double> fine_radii  = {0.1, 0.2, 0.25, 0.5, 0.8, 0.9, 1.3, 1.4, 2.0};
        std::vector<double> fine_angles = {
            0, M_PI / 16, M_PI / 8, M_PI / 2, M_PI, M_PI + M_PI / 16, M_PI + M_PI / 8, M_PI + M_PI / 2, M_PI + M_PI};

        double Rmax = fine_radii.back();

        ProcessingType processing_type = ProcessingType::CPU_HYBRID;
        DomainGeometry domain_geometry;
        DensityProfileCoefficients density_profile_coefficients;

        ExtrapolationType extrapolation = ExtrapolationType::NONE;
        int FMG = 0;
        bool DirBC_Interior = true;

        auto finest_grid = std::make_unique<PolarGrid>(fine_radii, fine_angles);
        auto coarse_grid = std::make_unique<PolarGrid>(coarseningGrid(*finest_grid));

        auto finest_levelCache = std::make_unique<LevelCache>(processing_type, *finest_grid, density_profile_coefficients, domain_geometry);
        auto coarse_levelCache = std::make_unique<LevelCache>(processing_type, *coarse_grid, density_profile_coefficients, domain_geometry);

        Level finest_level(0, processing_type, std::move(finest_grid), std::move(finest_levelCache), extrapolation, FMG);
        Level coarse_level(1, processing_type, std::move(coarse_grid), std::move(coarse_levelCache), extrapolation, FMG);

        Interpolation interpolation_operator(DirBC_Interior);

        Vector<double> x = generate_random_sample_data(coarse_level.grid(), 42);
        GPU_Vector<double> gpu_x(x.size());
        copyHostToDevice(x, gpu_x);

        Vector<double> result(finest_level.grid().numberOfNodes());
        GPU_Vector<double> gpu_result(finest_level.grid().numberOfNodes());

        interpolation_operator.applyExtrapolatedProlongation(coarse_level, finest_level, result, x);
        interpolation_operator.applyExtrapolatedProlongation(coarse_level, finest_level, gpu_result, gpu_x);

        Vector<double> host_result(gpu_result.size());
        copyDeviceToHost(gpu_result, host_result);

        ASSERT_EQ(result.size(), host_result.size());
        for (int i = 0; i < result.size(); ++i) {
            ASSERT_DOUBLE_EQ(result[i], host_result[i]);
        }
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}