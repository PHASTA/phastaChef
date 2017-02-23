#!/bin/bash -x

# load environment variables
source /usr/local/etc/bash_profile
module load cmake/latest
module load mpich3/3.1.2-thread-multiple
module load parmetis/mpich3.1.2/4.0.3
module load zoltan/mpich3.1.2/3.8
module load simmetrix/simModSuite
module load pumi

#cdash output root
cd /lore/cwsmith/cdash

#download test data
set +e
rm -rf phastaChefTests
rm *.tar.gz
set -e
wget www.scorec.rpi.edu/~cwsmith/phastaChefTests.tar.gz
tar xf phastaChefTests.tar.gz
set +e

#remove compilation directories created by nightly.cmake
rm -rf build/

#run nightly.cmake script
ctest -VV -D Nightly -S /lore/cwsmith/cdash/phastaChefRef/cdash/nightly.cmake &> cmake_log.txt

#core repository built by nightly.cmake
buildDir=/lore/cwsmith/cdash/build/phastaChef-sim
[ ! -e ${buildDir} ] && exit 0

cd $buildDir
#build the Doxygen html documentation
make doc
if [ -d "$PWD/doc/html" ]; then
  webdocs=$HOME/www/phastaChefDocs
  #remove the old web documentation
  rm -rf $webdocs
  #replace it with the generated one
  cp -r doc/html $webdocs
fi
