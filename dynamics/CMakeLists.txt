cmake_minimum_required(VERSION 3.1)

project(quadrotor-control
  DESCRIPTION "optimal control for quadrotor"
  LANGUAGES CXX
  VERSION 1.0
)

add_library(dynamics
            quadrotor_plant.cpp 
            roll_pitch_yaw.cpp 
            quadrotor_plant.hpp
            roll_pitch_yaw.hpp)

target_link_libraries(dynamics PUBLIC casadi)
target_include_directories(dynamics INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

