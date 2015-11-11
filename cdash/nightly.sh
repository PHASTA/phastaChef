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

#run nightly.cmake script
ctest -VV -D Nightly -S /lore/cwsmith/cdash/phastaChefRef/cdash/nightly.cmake &> cmake_log.txt

#remove compilation directories created by nightly.cmake
cd /lore/cwsmith/cdash
rm -rf build/
