#!/bin/bash -fx

cd /src/siconos
# on host machine:
# git clone git@gricad-gitlab.univ-grenoble-alpes.fr:nonsmooth/siconos.git
git checkout --track origin/noselecs_v0

# build numerics
cd /src
mkdir -p build/siconos/Numerics
cd build/siconos/Numerics
cmake /src/siconos/Numerics/ -DWITH_CMAKE_BUILD_TYPE=Debug
make -j8 install
# build kernel
cd /src
mkdir -p build/siconos/Kernel
cd build/siconos/Kernel
cmake /src/siconos/Kernel/ -DWITH_CMAKE_BUILD_TYPE=Debug
make -j8
cp Settings.cmake BuildSettings.cmake # it seems that a file is missing
make install

# build 
cd /src/noselecs
# on host machine:
# git clone git@gricad-gitlab.univ-grenoble-alpes.fr:nonsmooth/siconos.git
git checkout --track origin/noselecs_v0

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
