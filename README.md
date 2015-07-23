# setup
    git clone git@github.com:cwsmith/phastaChef.git
    cd phastaChef
    git clone git@github.com:cwsmith/phasta.git
    git clone git@github.com:SCOREC/core.git --branch chefLib core

    wget www.scorec.rpi.edu/~cwsmith/phastaChefTests.tar.gz .
    tar xzf phastaChefTests.tar.gz # use for CASES path

    wget https://www.scorec.rpi.edu/pumi/pumi_test_meshes.tar.gz
    tar xzf pumi_test_meshes.tar.gz # use for MESHES path

# build and test

    mpi=/path/to/mpi/install
    export PATH=$mpi/bin:$PATH
    opt="-Wextra -pedantic -g -O2 -isystem $mpi/include "

    cmake .. \
    -DCMAKE_C_COMPILER=mpicc \
    -DCMAKE_CXX_COMPILER=mpicxx \
    -DCMAKE_Fortran_COMPILER=gfortran \
    -DCMAKE_C_FLAGS="$opt" \
    -DCMAKE_CXX_FLAGS="$opt" \
    -DCMAKE_EXE_LINKER_FLAGS="-ldl $opt" \
    -DCMAKE_INSTALL_PREFIX=$PWD/install_nothread/ \
    \
    -DPCU_COMPRESS=ON \
    -DENABLE_THREADS=OFF \
    -DIS_TESTING=True \
    -DMESHES=/path/to/core/testMeshes \
    \
    -DPHASTA_INCOMPRESSIBLE=ON \
    -DPHASTA_COMPRESSIBLE=ON \
    -DACUSOLVE_LIB=/path/to/libles.a \
    -DCASES=/path/to/phastaChefTests/ \
    ..

    make
    make test


# development discussion

http://fluid.colorado.edu/wiki/index.php/Work_Plan_In_Memory_Adapt#Detailed_Description_of_Phases
