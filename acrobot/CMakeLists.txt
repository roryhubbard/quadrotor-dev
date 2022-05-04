cmake_minimum_required(VERSION 3.1)

project(acrobot-control
  DESCRIPTION "optimal control for acrobot"
  LANGUAGES CXX
  VERSION 1.0
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(casadi REQUIRED)
find_package(nlohmann_json 3.10.5 REQUIRED)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../dynamics ${CMAKE_CURRENT_BINARY_DIR}/dynamics)
list(APPEND EXTRA_LIBS dynamics)

add_executable(acrobot_control acrobot_control.cpp acrobot_plant.cpp acrobot_plant.hpp)
add_executable(acrobot_control2 acrobot_control2.cpp acrobot_plant.cpp acrobot_plant.hpp)
add_executable(acrobot_simulate acrobot_simulate.cpp acrobot_plant.cpp acrobot_plant.hpp)

target_link_libraries(acrobot_control PRIVATE casadi ${EXTRA_LIBS})
target_link_libraries(acrobot_control2 PRIVATE casadi nlohmann_json::nlohmann_json ${EXTRA_LIBS})
target_link_libraries(acrobot_simulate PRIVATE casadi nlohmann_json::nlohmann_json ${EXTRA_LIBS})

