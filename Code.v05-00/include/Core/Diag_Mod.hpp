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
#include "Core/Structure.hpp"
#include "Core/Mesh.hpp"
#include "Core/Meteorology.hpp"
#include "Core/Util.hpp"
#include "KPP/KPP_Global.h"
#include <netcdf>

using namespace netCDF;
using namespace netCDF::exceptions;

static const int SAVE_SUCCESS = 1;
static const int SAVE_FAILURE = 0;

/* ================================================================== */
/* ---- Timeseries Diagnostics -------------------------------------- */
/* ================================================================== */

/* Timeseries diagnostic files must be of the form:
 *      *hhmmss.nc or *hhmm.nc */

bool Diag_TS_Chem( const char* ROOTNAME,                     \
                   const std::vector<int> speciesIndices,    \
                   const int hh, const int mm, const int ss, \
                   const Solution& Data, const Mesh& m );

bool Diag_TS_Phys( const char* ROOTNAME,                     \
                   const std::vector<int> aerosolIndices,    \
                   const int hh, const int mm, const int ss, \
                   const Solution& Data, const Mesh& m,      \
                   const Meteorology &Met,                   \
                   const int outputPDF = 0,                  \
	           const float partNum_lost = 1.0,             \
	           const float iceMass_lost = 1.0 );

/* ================================================================== */
/* ---- Prod & Loss Rates Diagnostics ------------------------------- */
/* ================================================================== */

/* If chemistry is performed at the grid cell level, then the rates
 * are stored as:
 * NY x NX x NFAM 
 * into netCDF files at a frequency specified by the input file */

//bool Diag_PL( const char* ROOTNAME,                     \
//              const int hh, const int mm, const int ss, \
//              const Solution& Data,                     \
//              const Mesh& m );

#endif /* DIAG_MOD_H_INCLUDED */
