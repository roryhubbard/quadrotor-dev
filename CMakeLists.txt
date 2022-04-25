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

add_subdirectory("${PROJECT_SOURCE_DIR}/dynamics")
list(APPEND EXTRA_LIBS dynamics)

add_executable(quadrotor_control "${PROJECT_SOURCE_DIR}/control/quadrotor_control.cpp")

target_link_libraries(quadrotor_control PRIVATE casadi ${EXTRA_LIBS})

