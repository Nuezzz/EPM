#!/bin/bash
#
# Creates a configuration and compiles
# MC3D for the Release version in 
# "Source/Debug/bin/MC3D_debug"
vers="Debug"
if [ -d "./$vers" ]; then
    rm -r $vers
fi
mkdir $vers
cd $vers

cmake -DCMAKE_BUILD_TYPE=$vers ../
make
