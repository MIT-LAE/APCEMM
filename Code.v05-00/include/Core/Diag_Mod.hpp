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
#include "Util/VectorUtils.hpp"
#include <netcdf>
#include <filesystem>

namespace Diag {
    using namespace netCDF;
    using namespace netCDF::exceptions;
    using std::string;
    using std::vector;


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

    void Diag_TS_Phys( const char* rootName,
                    const int hh, const int mm, const int ss,
                    const AIM::Grid_Aerosol& iceAer, const Vector_2D& H2O,
                    const Vector_1D& xCoord, const Vector_1D& yCoord,
                    const Vector_1D& xEdges, const Vector_1D& yEdges,
                    const Meteorology &met);
    
    void add0DVar(NcFile& currFile, const float toSave, const NcDim& dim, const string& name, const string& desc, const string& units);
    void add1DVar(NcFile& currFile, const Vector_1D& toSave, const NcDim& dim, const string& name, const string& desc, const string& units);
    void add2DVar(NcFile& currFile, const Vector_2D& toSave, const vector<NcDim> dims, const string& name, const string& desc, const string& units);
    void replace_hhmmss(string& fileName, int hh, int mm, int ss);

    /* ================================================================== */
    /* ---- Prod & Loss Rates Diagnostics ------------------------------- */
    /* ================================================================== */

    /* If chemistry is performed at the grid cell level, then the rates
    * are stored as:
    * NY x NX x NFAM 
    * into netCDF files at a frequency specified by the input file */

    /*
    bool Diag_PL( const char* ROOTNAME,                     \
                 const int hh, const int mm, const int ss, \
                 const Solution& Data,                     \
                 const Mesh& m );
    */

    #endif /* DIAG_MOD_H_INCLUDED */

}
