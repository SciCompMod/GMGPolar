enable_testing()

set(TESTSOURCES comparetest.cpp test_prolongation.cpp test_restriction.cpp mockgrid.cpp)

add_executable(
  gmgpolar_tests
  ${TESTSOURCES}
)
target_link_libraries(
  gmgpolar_tests
  GTest::gtest_main
  GMGPolar
)


target_include_directories(gmgpolar_tests PRIVATE . ../include ../include/test_cases)
include(GoogleTest)
gtest_discover_tests(gmgpolar_tests)