/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Diag_Mod Header File                                             */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 12/17/2018                                */
/* File                 : Diag_Mod.hpp                              */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef DIAG_MOD_H_INCLUDED
#define DIAG_MOD_H_INCLUDED

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include "Core/Interface.hpp"
#include "Core/FileHandler.hpp"
#include "Core/Structure.hpp"
#include "Core/Mesh.hpp"
#include "Core/Util.hpp"

/* ================================================================== */
/* ---- Timeseries Diagnostics -------------------------------------- */
/* ================================================================== */

/* Timeseries diagnostic files must be of the form:
 *      *hhmm.nc */

bool Diag_TS( const char* rootName, const int hh, const int mm, \ 
              const std::vector<int> speciesIndices, \
              const Solution& Data, const Mesh& m );

#endif /* DIAG_MOD_H_INCLUDED */
