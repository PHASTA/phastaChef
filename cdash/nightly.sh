#!/bin/bash -x

# load environment variables
source /usr/local/etc/bash_profile
module load cmake/latest
module load mpich3/3.1.2-thread-multiple
module load parmetis/mpich3.1.2/4.0.3
module load zoltan/mpich3.1.2/3.8
module load simmetrix/simModSuite

#cdash output root
cd /lore/cwsmith/cdash

#download test data
set +e
rm -rf meshes phastaChefTests
rm *.tar.gz
set -e
wget www.scorec.rpi.edu/pumi/pumi_test_meshes.tar.gz
tar xf pumi_test_meshes.tar.gz
wget www.scorec.rpi.edu/~cwsmith/phastaChefTests.tar.gz
tar xf phastaChefTests.tar.gz
set +e

#remove compilation directories created by nightly.cmake
rm -rf build/

#run nightly.cmake script
ctest -VV -D Nightly -S /lore/cwsmith/cdash/phastaChefRef/cdash/nightly.cmake &> cmake_log.txt
