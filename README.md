# Code Status
Check the nightly status on the
[CDash page](http://my.cdash.org/index.php?project=phastaChef)

# Setup
    git clone git@github.com:cwsmith/phastaChef.git
    cd phastaChef
    git clone git@github.com:PHASTA/phasta.git
    git clone git@github.com:SCOREC/core.git

    wget www.scorec.rpi.edu/~cwsmith/phastaChefTests.tar.gz .
    tar xzf phastaChefTests.tar.gz # use for CASES path below

    wget https://www.scorec.rpi.edu/pumi/pumi_test_meshes.tar.gz
    tar xzf pumi_test_meshes.tar.gz # use for MESHES path below

# Build and Test

The following example compilation instructions assume that both the phasta and core repos were
downloaded following the Setup instructions above.  If you already have the
repos checked out the variables `-DCORE_SRC_DIR=/path/to/core` and
`-DPHASTA_SRC_DIR=/path/to/phasta` can be passed to
the cmake command to specify their locations.

Note, the following disables the SVLS and PETSC solvers and relies on LESLIB for the incompressible solver and the native compressible solver.

    opt="-Wextra -pedantic -g -O2 "
    
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
    -DENABLE_ZOLTAN=ON \
    -DIS_TESTING=True \
    -DMESHES=/path/to/core/testMeshes \
    \
    -DPHASTA_INCOMPRESSIBLE=ON \
    -DPHASTA_COMPRESSIBLE=ON \
    -DPHASTA_USE_SVLS=OFF \
    -DPHASTA_USE_PETSC=OFF \
    -DLESLIB=/path/to/libles.a \
    -DPHASTA_TESTING=ON \
    -DCASES=/path/to/phastaChefTests/ \
    ..

    make
    make test


# development discussion

http://fluid.colorado.edu/wiki/index.php/Work_Plan_In_Memory_Adapt#Detailed_Description_of_Phases
