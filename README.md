[![DOI](https://zenodo.org/badge/256520978.svg)](https://zenodo.org/badge/latestdoi/256520978)

# README for the APCEMM Source code repository

APCEMM stands for Aircraft Plume Chemistry Emission and Microphysics Model. The model aims to assess the chemical and microphysical perturbations introduced by a conventional aircraft, equipped with gas turbine engines. Global chemistry transport models commonly assume that emissions are released in grid-boxes that can be several orders of magnitude greater than an airplane's typical dimensions. APCEMM accounts for the fine-scale representation of an aircraft plume. To account for plume-scale processes, APCEMM has the option of computing effective emissions that correspond to what should be released in a grid-box model to match the plume model's output.

Currently, we are focusing on the development of the contrail modeling components of APCEMM. The chemistry modules are for the time being incompatible with the current version of the code and hence should be disabled until development resumes. 

## APCEMM Development

The development of APCEMM in C++ started in September 2018. 

This repository contains multiple branches. Each branch pertains to a specific function.

* The __main__ branch always contains the most up-to-date and stable version. New code should never be added to that branch directly. Instead, a new branch, forked from master, should be created.
* The __dev*__ branch contains in-development code for future versions.

For VSCode users, a Docker Dev Container is defined in `.devcontainer`. See [the tutorial](https://code.visualstudio.com/docs/devcontainers/tutorial) to develop inside a containerized environment.

## Dependencies 

These are all managed using the `vcpkg` tool (see below) so do not need to be installed explicitly.

- netcdf-c (requires HDF5 and zlib)
- netcdf-cxx4
- Catch2
- FFTW3
- OpenMP
- Boost libraries
- yaml-cpp
- Eigen3

See the [Dockerfile](.devcontainer/Dockerfile.apcemm) in the .devcontainer directory for specifics.

## APCEMM: Installation instructions

APCEMM can be built using CMake. Previously, the dependency structure and compile instructions were specified using manually generated Makefiles. CMake generates these Makefiles automatically, and should lead to a more pleasant software build experience. Dependencies on external libraries are managed using the [vcpkg](https://vcpkg.io/en/) tool, which is installed as a Git submodule. (This means that you just need to run the `git submodule update` command below to set it up.)

CMake will generate a single executable `APCEMM` that can receive an input file `input.yaml`. To compile this executable, you can call CMake as follows:

```
git submodule update --init --recursive
mkdir build
cd build
cmake ../Code.v05-00
cmake --build .
```

The `git submodule update` command installs the `vcpkg` dependency management tool, and the first time that you run CMake, all of the C++ dependencies will be installed. This will take some time, but subsequent runs of CMake will use cached binary builds of the dependencies, so will be much quicker.

The above commands will generate the `APCEMM` executable in the `build` directory (an "out-of-source" build). It is also possible to perform a build directly in the `Code.v05-00` directory, but this is not preferred. You can perform an "out-of-source" build anywhere that it's convenient, simply by calling CMake from within a different directory. For example,
```
cd APCEMM/rundirs/SampleRunDir/
cmake ../../Code.v05-00
cmake --build .
```
will generate the executable in the `rundirs/SampleRunDir/` directory. 

## Getting Started
To start a run from the aforementioned `rundirs/SampleRunDir`, simply call:
```
./../../Code.v05-00 input.yaml
```
Three examples and their accompanying jupyter notebooks for postprocessing tutorials are provided in the `examples` folder. The first example is one where the contrail doesn't persists, and only focuses on analyzing the output of the early plume model (EPM) module of APCEMM. The second example is a persistent contrail simulation where the ice supersaturated layer depth is specified. The third example features using a meteorological input file.

The input file options are explained via comments in the file `rundirs/SampleRunDir/input.yaml`

Advanced simulation parameters hidden in the input files (e.g. Aerosol bin size ratios, minimum/max bin aerosol sizes, etc) can be modified in `Code.v05-00/src/include/Parameters.hpp`. 
