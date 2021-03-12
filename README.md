# Tlab

Tools to simulate and analyze turbulent flows in 2D and 3D configurations. The numerical schemes are based on compact finite differences with structured meshes, where grid stretching allowed. Time advancement in based on Runge-Kutta schemes. There are two possible modes [temporal|spatial], which corresponds to temporally evolving flows and spatially evolving flows, respectively. This parameter refers to the statistical homogeneities of the configuration, and it partly defines the boundary conditions and the calculation of statistical properties.

See doc/manual.pdf for more information.

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
