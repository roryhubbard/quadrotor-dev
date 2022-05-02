cmake_minimum_required(VERSION 3.1)

add_library(dynamics
            roll_pitch_yaw.cpp 
            system.cpp 
            roll_pitch_yaw.hpp
            system.hpp)

target_link_libraries(dynamics PUBLIC casadi)
target_include_directories(dynamics INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

