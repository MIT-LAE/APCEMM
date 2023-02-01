[![DOI](https://zenodo.org/badge/256520978.svg)](https://zenodo.org/badge/latestdoi/256520978)

# README for the APCEMM Source code repository

APCEMM stands for Aircraft Plume Chemistry Emission and Microphysics Model. The model aims to assess the chemical and microphysical perturbations introduced by a conventional aircraft, equipped with gas turbine engines. Global chemistry transport models commonly assume that emissions are released in grid-boxes that can be several orders of magnitude greater than an airplane's typical dimensions. APCEMM accounts for the fine-scale representation of an aircraft plume. To account for
plume-scale processes, APCEMM has the option of computing effective emissions that correspond to what should be released in a grid-box model to match the plume model's output.

## APCEMM Development

The development of APCEMM in C++ started in September 2018. 

This repository contains multiple branches. Each branch pertains to a specific function.

* The __master__ branch always contains the most up-to-date and stable version. New code should never be added to that branch directly. Instead, a new branch, forked from master, should be created.
* The __dev*__ and __feature/*__ branches contain in-development code for future versions.

## Dependencies 

APCEMM requires the following libraries to build:

netcdf-c (requires HDF5 and zlib)

netcdf-cxx4

Catch2

FFTW3

OpenMP

Boost libraries

yaml-cpp

## APCEMM: Installation instructions

APCEMM can be built using CMake. Previously, the dependency structure and compile instructions were specified using manually generated Makefiles. CMake generates these Makefiles automatically, and should lead to a more pleasant software build experience. 

CMake will generate a single executable `APCEMM` that can receive an input file `input.yaml`. To compile this executable, you can call cmake as follows

```
cd APCEMM/Code.v05-00
cmake .
cmake --build .
```
This will generate the `APCEMM` executable in the `Code.v05-00` directory. You can also perform an "out-of-source" build, by simply calling cmake from within a different directory. For example,
```
cd APCEMM/rundirs/SampleRunDir/
cmake ../../Code.v05-00
cmake --build .
```
will generate the executable in the `rundirs/SampleRunDir/` directory. 

APCEMM now requires a working installation of yaml-cpp. If you did not install yaml-cpp in a standard location, pass the argument `-DCMAKE_PREFIX_PATH=/path/to/yaml-cpp` as part of the original `cmake` command.

## Getting Started
To start a run from the aforementioned `rundirs/SampleRunDir`, simply call:
```
./../../Code.v05-00 input.yaml
```
Three examples and their accompanying jupyter notebooks for postprocessing tutorials are provided in the `examples` folder. The first example is one where the contrail doesn't persists, and only focuses on analyzing the output of the early plume model (EPM) module of APCEMM. The second example is a persistent contrail simulation where the ice supersaturated layer depth is specified. The third example features using a meteorological input file.

The input file options are explained via comments in the file `rundirs/SampleRunDir/input.yaml`

Advanced simulation parameters hidden in the input files (e.g. Domain size, grid spacing, aerosol bin size ratios, etc) can be modified in `Code.v05-00/src/include/Parameters.hpp`