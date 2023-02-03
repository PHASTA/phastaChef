# Setup

Download required dependences and code:

```
root=$PWD
git clone https://github.com/SCOREC/core.git
cd core
git checkout shock_detection
cd $root
git clone https://github.com/SCOREC/phasta.git
cd phasta
git checkout shock_detection
cd $root
git clone https://github.com/PHASTA/phastaChef.git
cd phastaChef
git checkout shockDetection
cd $root
```

# Build and Test
Instructions for CCI erp

Start in directory where all repos are stored

## Environment
```
ln -snf /usr/lib64/libslurm.so.30 $HOME/libslurm.so.27
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME

export http_proxy=http://proxy:8888
export https_proxy=$http_proxy

module purge
module use /gpfs/u/software/erp-spack-install/lmod/linux-centos7-x86_64/Core/
module load \
  gcc/7.4.0-vuvgy2j \
  openmpi/4.0.1-jklawme \
  metis/5.1.0-int32-p5jssbp \
  parmetis/4.0.3-int32-z5vu5me \
  zoltan/3.83-int32-ruqfm3u \
  simmetrix-simmodsuite/14.0-190617dev-ufy24ll \
  cmake

root=$PWD
```
## SCOREC Core
```
cd $root
mkdir core/build && cd !$

flags=" -03 -DNDEBUG -Wall -Wextra"

cmake3 .. \
  -DCMAKE_C_COMPILER="mpicc" \
  -DCMAKE_CXX_COMPILER="mpicxx" \
  -DCMAKE_C_FLAGS="${flags}" \
  -DCMAKE_CXX_FLAGS="${flags}" \
  -DSIM_MPI=openmpi4.0.1 \
  -DENABLE_ZOLTAN=ON \
  -DCMAKE_INSTALL_PREFIX=$PWD/install \
  -DENABLE_SIMMETRIX=True \
  -DENABLE_FIELDSIM=True \
  -DSIM_PARASOLID=True \
  -DPCU_COMPRESS=OFF \
  -DPUMI_FORTRAN_INTERFACE=ON \
  -DCMAKE_Fortran_COMPILER="mpif90" \

make install -j8
```
## phastaChef
```
cd $root
mkdir phastaChef/build && cd !$

opt="-Wextra -pedantic -g -O3 "

cmake .. \
-DCMAKE_C_COMPILER=mpicc \
-DCMAKE_CXX_COMPILER=mpicxx \
-DCMAKE_Fortran_COMPILER=mpif90 \
-DCMAKE_BUILD_TYPE=Release \
-DSCOREC_PREFIX=/gpfs/u/home/MCFD/MCFDstsp/scratch/code/detection/core/build/install/ \
-DPHASTA_SRC_DIR=/gpfs/u/home/MCFD/MCFDstsp/scratch/code/detection/phasta \
-DPHASTA_INCOMPRESSIBLE=OFF \
-DPHASTA_COMPRESSIBLE=ON \

make -j8
```
