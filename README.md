# Tlab

Tools to simulate and analyze turbulent flows in 2D and 3D configurations. The numerical schemes are based on compact finite differences with structured meshes, where grid stretching allowed. Time advancement in based on Runge-Kutta schemes. There are two possible modes [temporal|spatial], which corresponds to temporally evolving flows and spatially evolving flows, respectively. This parameter refers to the statistical homogeneities of the configuration, and it partly defines the boundary conditions and the calculation of statistical properties. Hybrid parallelization using MPI in two directions and OpenMPI in big loops.

See [`doc/manual.pdf`](./doc/manual.pdf) for more information.

## Install

In order to compile the code, run the following commands:

```shell
cd ${PATH_TO_TLAB}
mkdir build
cd build
cmake ../src -DSYST={mpipc,juqueen,...} -DBUILD_TYPE={BIG,LITTLE,PARALLEL,NONBLOCKING}
make
```
Instead of mpipc or juqueen, you have to use the corresponding file from the directory ${PATH_TO_TLAB}/config

You can also run `./configure.sh`, which would create the different build_* directories automatically for your system if appropriately set up.

To clean the tree, simply delete the directories build*

**Prerequisites**
* cmake (at a computing centre, try `module load cmake`)
* fortran compiler
* [FFTW](http://www.fftw.org/)

## Check

In order to check the code, run the following commands:

```shell
cd ${PATH_TO_DNS}
cd examples
make check BINS_PATH=${bins_LITTLE,bins_BIG,...} PRJS=${Case01,Case02,...}
```

Instead of bins_LITTLE or bins_BIG, you have to use the corresponding directory whose executables you want to check.

If the variable PRJS is not passed, all projects in the directory examples will be checked.

To clean the examples tree, run `make clean`.

Use valgrind to check for memory leaks.

## Directory organization

* [`config`](./config): configuration files for different architectures and compilers
* [`doc`](./doc): documentation files
* [`examples`](./examples): sample cases for reference and for development testing
* [`scripts`](./scripts): shell and python scripts for job management in supercomputing centers and postprocessing
* [`src`](./src): source files  
  * [`src/tools`](./src/tools): of executables
  * [`src/mappings`](./src/mappings): of library with level1 mappings that depend on level0 operators but not on fdm routines
  * [`src/operators`](./src/operators): of library with level0 operators that depend only on fdm routines
  * [`src/fdm`](./src/fdm): of library with finite difference schemes
  * [`src/utils`](./src/utils): of library with basic operations with dependencies
  * [`src/external`](./src/external): of external libraries
  * ...

## Data structure

**Restart files**

* IEEE double precision floating point representation
* Sequential access mode
* Unformatted records
* Endianness determined in the compiler options
* No record length information: first 4-bytes are an integer pointing to the field array after the header

**Arrays**

* The inner-most index is the streamwise direction
* The outer-most index is the spanwise direction
* MPI domain decomposition in these two directions
