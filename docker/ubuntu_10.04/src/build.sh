#!/bin/bash -fx

cd siconos
git checkout 62bdcb90c124eb41f74c6694f292d107f9ad73f7 -b noselecs_v0

# build numerics
cd /src
mkdir -p build/siconos/Numerics
cd build/siconos/Numerics
cmake /src/siconos/Numerics/
make -j8 install
# build kernel
cd /src
mkdir -p build/siconos/Kernel
cd build/siconos/Kernel
cmake /src/siconos/Kernel/
make -j8
cp Settings.cmake BuildSettings.cmake # it seems that a file is missing
make install

# build 
cd /src/noselecs
git checkout 7e9e845695097bb7a32c8a17e4aec2814a228fed -b noselecs_v0

cd /src/noselecs/trunk/src/Spiceparser/
ln -sf /src/siconos/Numerics/cmake .
mkdir build ; cd build
cmake ..
make -j8 install

cd /src/noselecs/trunk/src/ACE/
ln -sf /src/siconos/Numerics/cmake .
mkdir build ; cd build
cmake ..
make -j8 install

cd /src/noselecs/trunk/src/EXE/
ln -sf /src/siconos/Numerics/cmake .
mkdir build ; cd build
cmake ..
make -j8
