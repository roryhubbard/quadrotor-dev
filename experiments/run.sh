#!/bin/bash

#export CMAKE_PREFIX_PATH=../casadi_install/casadi/cmake/
rm -rf build
mkdir build
pushd build
cmake ..
cmake --build .
popd
build/casadi_eigen_comparison

