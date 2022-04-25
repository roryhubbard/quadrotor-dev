cmake_minimum_required(VERSION 3.16)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(experiments
  DESCRIPTION "random experiments"
  LANGUAGES CXX
)

set(CMAKE_CXX_COMPILER /usr/bin/g++)

find_package(casadi REQUIRED)
find_package(Eigen3 REQUIRED)

include_directories("${EIGEN3_INCLUDE_DIR}")

add_executable(casadi_eigen_comparison casadi_eigen_comparison.cpp)

target_link_libraries(casadi_eigen_comparison casadi "${EIGEN_LIBRARIES}")

