cmake_minimum_required(VERSION 3.16.3)

project(gmgpolartest)


file(GLOB sources_src "./src/*.cpp" "./src/test_cases/*.cpp")

add_library(foruse ${sources_src}) 
add_executable(main_allantest ./src/main.cpp)
target_include_directories(main_allantest PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/include/test_cases ${CMAKE_SOURCE_DIR}/olderruns)
target_include_directories(foruse PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/include/test_cases ${CMAKE_SOURCE_DIR}/olderruns)

find_package(OpenMP)
if(OPENMP_FOUND)
    set_target_properties(main_allantest PROPERTIES COMPILE_FLAGS "-O2 -Wall -MMD -MP -Wwrite-strings -fopenmp")
    set_target_properties(main_allantest PROPERTIES LINK_FLAGS "-O2 -Wall -MMD -MP -Wwrite-strings -fopenmp")

    set_target_properties(foruse PROPERTIES COMPILE_FLAGS "-O2 -Wall -MMD -MP -Wwrite-strings -fopenmp")
    set_target_properties(foruse PROPERTIES LINK_FLAGS "-O2 -Wall -MMD -MP -Wwrite-strings -fopenmp")
else()
    message(FATAL_ERROR "OpenMP needed")
endif()

target_link_libraries(main_allantest PRIVATE foruse)

add_subdirectory(tests)

