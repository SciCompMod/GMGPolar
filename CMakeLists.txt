# cmake_minimum_required(VERSION 3.18)
# project(GMGPolarGPU LANGUAGES CXX CUDA)

# option(GMGPOLAR_BUILD_TESTS "Build GMGPolar unit tests." OFF)

# set(CMAKE_CXX_STANDARD 20)
# set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(CMAKE_CUDA_STANDARD 17)
# set(CMAKE_CUDA_STANDARD_REQUIRED ON)
# set(CMAKE_CUDA_ARCHITECTURES 70)
# enable_language(CUDA)

# set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -pedantic -Wno-unused")
# set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2 -mtune=generic")
# set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler -fopenmp -arch=sm_70")
# set(CMAKE_CUDA_FLAGS_RELEASE "${CMAKE_CUDA_FLAGS_RELEASE} -O3")
# set(CMAKE_CUDA_FLAGS_DEBUG "${CMAKE_CUDA_FLAGS_DEBUG} -g -O0")

# set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)
# set(CMAKE_CUDA_RESOLVE_DEVICE_SYMBOLS ON)

# # Find CUDA and OpenMP packages
# find_package(CUDA REQUIRED)
# find_package(OpenMP REQUIRED)

# # OpenMP for host code (C++):
# if(OpenMP_CXX_FOUND)
#     message(STATUS "OpenMP enabled for C++")
#     set(OpenMP_CXX_LIB OpenMP::OpenMP_CXX)
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
# else()
#     message(WARNING "OpenMP not found for C++")
# endif()

# # Mumps: Sparse Matrix Solver
# set(MUMPS_PREFIX_PATH "~/spack/opt/spack/linux-debian11-ivybridge/gcc-10.2.1/mumps-5.6.2-hm3vnybofcspkmjb6xdhttmdt3uetvo4")
# include_directories(${MUMPS_PREFIX_PATH}/include)
# link_directories(${MUMPS_PREFIX_PATH}/lib)

# # Metis: Matrix reordering for Mumps
# set(METIS_PREFIX_PATH "~/spack/opt/spack/linux-debian11-ivybridge/gcc-10.2.1/metis-5.1.0-wkmpso5x6ytk7u7t5z7a5h4r4xtm35lv")
# include_directories(${METIS_PREFIX_PATH}/include)
# link_directories(${METIS_PREFIX_PATH}/lib)

# set(MUMPS_LIBRARIES
#     mumps_common
#     smumps
#     dmumps
#     mpiseq
#     metis
# )

# # Include directories for the project
# include_directories(include)

# # 1. Create a library target for the PolarGrid module
# file(GLOB_RECURSE POLAR_GRID_SOURCES "src/PolarGrid/*.cpp" "src/PolarGrid/*.cu")
# add_library(PolarGrid STATIC ${POLAR_GRID_SOURCES})
# target_include_directories(PolarGrid PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/PolarGrid)

# # 2. Create a library target for the InputFunctions module
# file(GLOB_RECURSE INPUT_FUNCTIONS_SOURCES "src/InputFunctions/*.cpp" "src/InputFunctions/*.cu")
# add_library(InputFunctions STATIC ${INPUT_FUNCTIONS_SOURCES})
# target_include_directories(InputFunctions PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/InputFunctions)

# file(GLOB_RECURSE GMG_POLAR_SOURCES "src/GMGPolar/*.cpp" "src/GMGPolar/*.cu")
# file(GLOB_RECURSE LEVEL_SOURCES "src/Level/*.cpp" "src/Level/*.cu")
# file(GLOB_RECURSE INTERPOLATION_SOURCES "src/Interpolation/*.cpp" "src/Interpolation/*.cu")
# file(GLOB_RECURSE DIRECT_SOLVER_SOURCES "src/DirectSolver/*.cpp")
# file(GLOB_RECURSE STENCIL_SOURCES "src/Stencil/*.cpp")
# file(GLOB_RECURSE RESIDUAL_SOURCES "src/Residual/*.cpp" "src/Residual/*.cu")

# # 4. Create the GMGPolarLib library and link PolarGrid and InputFunctions
# add_library(GMGPolarLib STATIC
#     ${GMG_POLAR_SOURCES}
#     ${LEVEL_SOURCES}
#     ${INTERPOLATION_SOURCES}
#     ${DIRECT_SOLVER_SOURCES}
#     ${STENCIL_SOURCES}
#     ${RESIDUAL_SOURCES}
# )
# # Link PolarGrid and InputFunctions to GMGPolarLib
# target_link_libraries(GMGPolarLib PUBLIC PolarGrid InputFunctions)

# set(MAIN_SOURCE "src/main.cpp")
# add_executable(gmgpolar ${MAIN_SOURCE})

# target_link_libraries(gmgpolar PRIVATE GMGPolarLib ${MUMPS_LIBRARIES} ${OpenMP_CXX_LIB})







cmake_minimum_required(VERSION 3.18)
project(GMGPolarGPU LANGUAGES CXX CUDA)

option(GMGPOLAR_BUILD_TESTS "Build GMGPolar unit tests." ON)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_ARCHITECTURES 70)
enable_language(CUDA)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -pedantic -Wno-unused")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2 -mtune=generic")
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler -fopenmp -arch=sm_70 -Wno-deprecated-declarations")
set(CMAKE_CUDA_FLAGS_RELEASE "${CMAKE_CUDA_FLAGS_RELEASE} -O3")
set(CMAKE_CUDA_FLAGS_DEBUG "${CMAKE_CUDA_FLAGS_DEBUG} -g -O0")

set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)
set(CMAKE_CUDA_RESOLVE_DEVICE_SYMBOLS ON)

# Find CUDA and OpenMP packages
find_package(CUDA REQUIRED)
find_package(OpenMP REQUIRED)

# OpenMP for host code (C++):
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP enabled for C++")
    set(OpenMP_CXX_LIB OpenMP::OpenMP_CXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
else()
    message(WARNING "OpenMP not found for C++")
endif()

# Mumps: Sparse Matrix Solver
set(MUMPS_PREFIX_PATH "~/spack/opt/spack/linux-debian11-ivybridge/gcc-10.2.1/mumps-5.6.2-hm3vnybofcspkmjb6xdhttmdt3uetvo4")
include_directories(${MUMPS_PREFIX_PATH}/include)
link_directories(${MUMPS_PREFIX_PATH}/lib)

# Metis: Matrix reordering for Mumps
set(METIS_PREFIX_PATH "~/spack/opt/spack/linux-debian11-ivybridge/gcc-10.2.1/metis-5.1.0-wkmpso5x6ytk7u7t5z7a5h4r4xtm35lv")
include_directories(${METIS_PREFIX_PATH}/include)
link_directories(${METIS_PREFIX_PATH}/lib)

set(MUMPS_LIBRARIES
    mumps_common
    smumps
    dmumps
    mpiseq
    metis
)

# 1. Create a library target for the PolarGrid module
file(GLOB_RECURSE POLAR_GRID_SOURCES "src/PolarGrid/*.cpp" "src/PolarGrid/*.cu")
add_library(PolarGrid STATIC ${POLAR_GRID_SOURCES})
target_include_directories(PolarGrid PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/PolarGrid)

# 2. Create a library target for the InputFunctions module
file(GLOB_RECURSE INPUT_FUNCTIONS_SOURCES "src/InputFunctions/*.cpp" "src/InputFunctions/*.cu")
add_library(InputFunctions STATIC ${INPUT_FUNCTIONS_SOURCES})
target_include_directories(InputFunctions PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/InputFunctions)

file(GLOB_RECURSE GMG_POLAR_SOURCES "src/GMGPolar/*.cpp" "src/GMGPolar/*.cu" "src/GMGPolar/MultigridMethods/*.cpp" "src/GMGPolar/MultigridMethods/*.cu")
file(GLOB_RECURSE LEVEL_SOURCES "src/Level/*.cpp" "src/Level/*.cu")
file(GLOB_RECURSE INTERPOLATION_SOURCES "src/Interpolation/*.cpp" "src/Interpolation/*.cu")
file(GLOB_RECURSE DIRECT_SOLVER_SOURCES "src/DirectSolver/*.cpp")
file(GLOB_RECURSE STENCIL_SOURCES "src/Stencil/*.cpp")

file(GLOB_RECURSE RESIDUAL_TAKE_CPU_SOURCES "src/Residual/ResidualTakeCPU/*.cpp")
file(GLOB_RECURSE RESIDUAL_TAKE_GPU_SOURCES "src/Residual/ResidualTakeGPU/*.cpp" "src/Residual/ResidualTakeGPU/*.cu")

file(GLOB_RECURSE SMOOTHER_TAKE_CPU_SOURCES "src/Smoother/SmootherTakeCPU/*.cpp")
file(GLOB_RECURSE SMOOTHER_TAKE_GPU_SOURCES "src/Smoother/SmootherTakeGPU/*.cpp" "src/Smoother/SmootherTakeGPU/*.cu")

file(GLOB_RECURSE EXTRAPOLATED_SMOOTHER_TAKE_CPU_SOURCES "src/ExtrapolatedSmoother/ExtrapolatedSmootherTakeCPU/*.cpp")
file(GLOB_RECURSE EXTRAPOLATED_SMOOTHER_TAKE_GPU_SOURCES "src/ExtrapolatedSmoother/ExtrapolatedSmootherTakeGPU/*.cpp" "src/ExtrapolatedSmoother/ExtrapolatedSmootherTakeGPU/*.cu")

# 4. Create the GMGPolarLib library and link PolarGrid and InputFunctions
add_library(GMGPolarLib STATIC
    ${GMG_POLAR_SOURCES}
    ${LEVEL_SOURCES}
    ${INTERPOLATION_SOURCES}

    ${RESIDUAL_TAKE_CPU_SOURCES}
    ${RESIDUAL_TAKE_GPU_SOURCES}

    ${SMOOTHER_TAKE_CPU_SOURCES}
    ${SMOOTHER_TAKE_GPU_SOURCES}

    ${DIRECT_SOLVER_SOURCES}
    ${STENCIL_SOURCES}

    ${EXTRAPOLATED_SMOOTHER_TAKE_CPU_SOURCES}
    ${EXTRAPOLATED_SMOOTHER_TAKE_GPU_SOURCES}
)
# Link PolarGrid and InputFunctions to GMGPolarLib
target_link_libraries(GMGPolarLib PUBLIC PolarGrid InputFunctions ${MUMPS_LIBRARIES} ${OpenMP_CXX_LIB} cuda cudart cudadevrt cusparse cusolver)

set(MAIN_SOURCE "src/main.cpp")
add_executable(gmgpolar ${MAIN_SOURCE})

target_link_libraries(gmgpolar PRIVATE GMGPolarLib)

if(GMGPOLAR_BUILD_TESTS)
    enable_testing()
    add_subdirectory(third-party)
    add_subdirectory(tests)
endif()

