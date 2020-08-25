[![DOI](https://zenodo.org/badge/256520978.svg)](https://zenodo.org/badge/latestdoi/256520978)

# README for the APCEMM Source code repository

APCEMM stands for Aircraft Plume Chemistry Emission and Microphysics Model. The model aims to assess the chemical and microphysical perturbations introduced by a conventional aircraft, equipped with gas turbine engines. Global chemistry transport models commonly assume that emissions are released in grid-boxes that can be several orders of magnitude greater than airplane's typical dimensions. APCEMM account for the fine scale representation of an aircraft plume. To account for
plume-scale processes, APCEMM has the option of computing effective emissions that correspond to what should be released in a grid-box model to match the plume model's output.

The source code for APCEMM can be found in this repository (https://github.com/fritzt/APCEMM).

## APCEMM Development


The development of APCEMM in C++ started in September 2018. 

This repository contains multiple branches. Each branch pertains to a specific function.

* The __master__ branch always contains the most up-to-date and stable version. New code should never be added to that branch directly. Instead, a new branch, forked from master, should be created.
* The __dev*__ and __feature/*__ branches contain in-development code for future versions.

