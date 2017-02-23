# Code Status
Check the nightly status on the
[CDash page](http://my.cdash.org/index.php?project=phastaChef)

# Setup
    git clone git@github.com:cwsmith/phastaChef.git
    cd phastaChef
    git clone git@github.com:PHASTA/phasta.git

    wget www.scorec.rpi.edu/~cwsmith/phastaChefTests.tar.gz .
    tar xzf phastaChefTests.tar.gz # use for CASES path below

# Build and Test

Build SCOREC/core; see [SCOREC/core build
instructions](https://github.com/SCOREC/core/wiki/General-Build-instructions).  Note, you must specify `CMAKE_INSTALL_PREFIX` when running cmake and run `make install` to install the necessary CMake package files.

The following example compilation instructions assume that the phasta repo was
downloaded following the Setup instructions above.  If you already have the
phasta repo checked out the variable `-DPHASTA_SRC_DIR=/path/to/phasta` can 
be passed to the cmake command to specify its locations.

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
    -DSCOREC_PREFIX=/path/to/SCOREC/core/install/ \
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

# Documentation

After running CMake Doxygen documentation can be generated with

  make doc

The latest version of the docs are also available here

  https://www.scorec.rpi.edu/~cwsmith/phastaChefDocs

# development discussion

http://fluid.colorado.edu/wiki/index.php/Work_Plan_In_Memory_Adapt#Detailed_Description_of_Phases
