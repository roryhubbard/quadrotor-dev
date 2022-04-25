#!/bin/bash

rm -rf build
mkdir build
pushd build
cmake ..
cmake --build . -j8
popd
build/quadrotor_control

