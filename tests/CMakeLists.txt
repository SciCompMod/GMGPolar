enable_testing()

set(TESTSOURCES comparetest.cpp)

add_executable(
  test_parameters
  ${TESTSOURCES}
)
target_link_libraries(
  test_parameters
  GTest::gtest_main
  GMGPolarlib
)

target_include_directories(test_parameters PRIVATE ../include ../include/test_cases)
include(GoogleTest)
gtest_discover_tests(test_parameters)