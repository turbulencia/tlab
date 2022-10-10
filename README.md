# Tlab

Tools to simulate and analyze turbulent flows in 2D and 3D configurations. The numerical schemes are based on compact finite differences and Runge-Kutta methods. Meshes are structured and allow grid stretching. There are two possible simulation modes [temporal|spatial], which correspond to temporally evolving flows and spatially evolving flows, respectively. This parameter refers to the statistical homogeneities of the configuration, and it partly defines the boundary conditions and the calculation of statistical properties. The code implements a hybrid parallelization: MPI for a domain decomposition in the first and third directions, and OpenMPI in big loops.

Some examples of applications can be found in this [website](https://jpmellado.github.io/gallery.html).

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
* cmake
* fortran compiler
* [FFTW](http://www.fftw.org/)
* Optionally, [NetCDF](https://docs.unidata.ucar.edu/netcdf-c/current/building_netcdf_fortran.html) for statistical data.

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
* [`examples`](./examples): sample cases for reference and for debugging
* [`scripts`](./scripts): shell and python scripts for job management and postprocessing
* [`src`](./src): source files  
  * [`src/tools`](./src/tools): executables
  * [`src/mappings`](./src/mappings): library with mappings from 3D fields to 3D, 2D and 1D data
  * [`src/operators`](./src/operators): library with operators that depend only on fdm routines
  * [`src/fdm`](./src/fdm): library with finite difference schemes
  * [`src/utils`](./src/utils): library with basic generic operators
  * [`src/external`](./src/external): external libraries, if any
  * [`src/modules`](./src/modules): constants, variables, arrays and basic procedures
  * ...

## Library dependencies

tools → mappings → operators → {fdm,filter,io,thermo,utils,external} → tlab

## Data structure

**Restart files**

* IEEE double precision floating point representation
* Sequential access mode
* Unformatted records
* Endianness determined in the compiler options
* No record length information: first 4-bytes are an integer pointing to the field array after the header

**Arrays**

* The inner-most index is the first (streamwise) direction
* The outer-most index is the third (spanwise) direction
* MPI domain decomposition in these two directions
