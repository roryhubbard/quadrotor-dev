cmake_minimum_required(VERSION 3.16)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(quadrotor-control
  DESCRIPTION "optimal control for quadrotor"
  LANGUAGES CXX
)

find_package(casadi REQUIRED)

add_executable(quadrotor_control quadrotor_control.cpp)
target_link_libraries(quadrotor_control casadi)

