cmake_minimum_required(VERSION 3.10)
project(AutocalibrationLibVM)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(ROOTBA_INSTANTIATIONS_DOUBLE "Instatiate templates for Scalar=double." ON)
option(ROOTBA_INSTANTIATIONS_FLOAT "Instatiate templates for Scalar=float." ON)

if(ROOTBA_INSTANTIATIONS_DOUBLE)
  add_compile_definitions(ROOTBA_INSTANTIATIONS_DOUBLE)
endif()
if(ROOTBA_INSTANTIATIONS_FLOAT)
  add_compile_definitions(ROOTBA_INSTANTIATIONS_FLOAT)
endif()

# Add Eigen include directory manually
include_directories(/Data/gautt/bundle_adjustment/robust_povar/external/eigen)

add_executable(autocalibration main.cpp libmv.cpp testing_functions.cpp utils_for_testing.cpp)

enable_testing()
add_subdirectory(tests)