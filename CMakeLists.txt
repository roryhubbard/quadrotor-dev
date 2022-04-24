cmake_minimum_required(VERSION 3.16)
set (CMAKE_CXX_STANDARD 20)

#export CMAKE_PREFIX_PATH=../build/
find_package(casadi)

#add_executable(quadrotor quadrotor.cpp)
#target_link_libraries(quadrotor casadi)

add_executable(tmp tmp.cpp)
target_link_libraries(tmp casadi)

