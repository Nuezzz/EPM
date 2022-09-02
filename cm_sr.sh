#!/bin/bash
#
# Creates a configuration and compiles
# MC3D for the Release version in 
# "Source/Debug/bin/MC3D_debug"
vers="Release"
if [ -d "./$vers" ]; then
    rm -r $vers
fi
mkdir $vers
cd $vers

cmake3 -DCMAKE_BUILD_TYPE=$vers ../
make