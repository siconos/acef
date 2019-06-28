#!/bin/bash -fx

current_last_commit=`git log -1 --format="%H"`
echo $current_last_commit




cd /home/acary/build

\rm -rf bisect
mkdir bisect
cd bisect


\rm -rf siconos
mkdir siconos
cd siconos
\rm -rf /usr/local/include/Siconos
\rm -rf /usr/local/lib/libSiconos*
\rm -rf bisect
mkdir bisect
cd bisect

mkdir Numerics
cd Numerics
cmake /src/siconos/Numerics/
make -j8 install || exit 125

mkdir Kernel
cd Kernel
cmake /src/siconos/Kernel/
cp Settings.cmake BuildSettings.cmake # it seems that a file is missing
make -j8 install || exit 125


cd /home/acary/build/
cd bisect
cp -r /src/noselecs .

cd noselecs/trunk/src/Spiceparser/
ln -sf /src/siconos/Numerics/cmake .
mkdir build ; cd build
cmake ..
make -j8 install || exit 125

cd ../../ACE/
ln -sf /src/siconos/Numerics/cmake .
mkdir build ; cd build
cmake ..
make -j8 install || exit 125

cd ../../EXE/
ln -sf /src/siconos/Numerics/cmake .
mkdir build ; cd build
cmake ..
make -j8 install || exit 125


# test part

cd /home/acary/build/bisect/noselecs/Tests/ACE

/home/acary/build/bisect/noselecs/trunk/src/EXE/noselect LC.cir ENUM
status=$?

if [ $status -eq 0 ]
then
    echo $status
    exit $status
else
    git log -1
    echo $status
    exit 1
fi
