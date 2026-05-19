#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

// Define the main function only once
int main(int argc, char* argv[])
{
    Kokkos::InitializationSettings settings;
    // 8 (CPU) threads
    settings.set_num_threads(GMGPOLAR_TEST_THREADS);
    Kokkos::ScopeGuard kokkos_scope(settings);
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
