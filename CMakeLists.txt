cmake_minimum_required(VERSION 3.14)
project(GMGPolar LANGUAGES CXX CUDA)

# add_definitions(-DUSE_MUMPS)
add_definitions(-DCUDA)

find_package(OpenMP REQUIRED)

add_library(iexmg)

target_link_libraries(iexmg PUBLIC OpenMP::OpenMP_CXX)

target_include_directories(
  iexmg PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
  $<INSTALL_INTERFACE:include>)

if (CUDA)
    find_package(CUDAToolkit)
    target_link_libraries(iexmg PRIVATE CUDA::cudart)
endif()

find_library(CUSOLVER_LIBRARY NAMES cusolver HINTS $ENV{CUSPARSE_LIBRARY_PATH})
find_library(CUBLAS_LIBRARY   NAMES cublas   HINTS $ENV{CUSPARSE_LIBRARY_PATH})
find_library(CUBLASLT_LIBRARY NAMES cublasLt HINTS $ENV{CUSPARSE_LIBRARY_PATH})
find_library(CUSPARSE_LIBRARY NAMES cusparse HINTS $ENV{CUSPARSE_LIBRARY_PATH})
add_subdirectory(src)
target_link_libraries(iexmg PUBLIC ${CUSOLVER_LIBRARY})
target_link_libraries(iexmg PUBLIC ${CUBLAS_LIBRARY})
target_link_libraries(iexmg PUBLIC ${CUBLASLT_LIBRARY})
target_link_libraries(iexmg PUBLIC ${CUSPARSE_LIBRARY})


