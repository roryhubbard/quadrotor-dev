#!/bin/bash

rm -rf build
mkdir build
pushd build
cmake ..
make -j8
popd
build/acrobot_control

