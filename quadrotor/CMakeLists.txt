cmake_minimum_required(VERSION 3.1)

project(quadrotor-control
  DESCRIPTION "optimal control for quadrotor"
  LANGUAGES CXX
  VERSION 1.0
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(casadi REQUIRED)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../dynamics ${CMAKE_CURRENT_BINARY_DIR}/dynamics)
list(APPEND EXTRA_LIBS dynamics)

add_executable(quadrotor_control quadrotor_control.cpp quadrotor_plant.cpp quadrotor_plant.hpp)

target_link_libraries(quadrotor_control PRIVATE casadi ${EXTRA_LIBS})

