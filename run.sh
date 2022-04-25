#!/bin/bash

rm -rf build
mkdir build
pushd build
cmake ..
cmake --build .
popd
build/quadrotor_control

