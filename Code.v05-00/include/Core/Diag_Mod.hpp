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

#include <string>
#include <vector>
#include <netcdf>
#include "AIM/Aerosol.hpp"
#include "Core/Meteorology.hpp"
#include "KPP/KPP_Global.h"

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

    void Diag_TS_Phys( const char* rootName,
                    const int hh, const int mm, const int ss,
                    const AIM::Grid_Aerosol& iceAer, const Vector_2D& H2O,
                    const Vector_1D& xCoord, const Vector_1D& yCoord,
                    const Vector_1D& xEdges, const Vector_1D& yEdges,
                    const Meteorology &met);
    
    void add0DVar(NcFile& currFile, const float toSave, const NcDim& dim, const string& name, const string& desc, const string& units);
    void add1DVar(NcFile& currFile, const Vector_1D& toSave, const NcDim& dim, const string& name, const string& desc, const string& units);
    void add2DVar(NcFile& currFile, const Vector_2D& toSave, const vector<NcDim> dims, const string& name, const string& desc, const string& units);
    void add3DVar(NcFile& currFile, const Vector_3D& toSave, const vector<NcDim> dims, const string& name, const string& desc, const string& units);
    void replace_hhmmss(string& fileName, int hh, int mm, int ss);
    void set_storePSD(bool val);

} // namespace Diag

#endif /* DIAG_MOD_H_INCLUDED */
