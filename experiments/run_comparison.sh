#!/bin/bash

#export CMAKE_PREFIX_PATH=../casadi_install/casadi/cmake/
rm -rf build
mkdir build
pushd build
cmake ..
make -j8
popd
build/casadi_eigen_comparison

