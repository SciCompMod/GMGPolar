#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

// Define the main function only once
int main(int argc, char* argv[])
{
    Kokkos::ScopeGuard kokkos_scope(argc, argv);
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
