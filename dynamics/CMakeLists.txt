cmake_minimum_required(VERSION 3.1)

project(quadrotor-control
  DESCRIPTION "optimal control for quadrotor"
  LANGUAGES CXX
  VERSION 1.0
)

add_library(dynamics
            quadrotor_plant.cpp 
            quadrotor_plant.hpp)

add_executable(quadrotor_control quadrotor_control.cpp)

target_link_libraries(quadrotor_control PUBLIC casadi dynamics)
