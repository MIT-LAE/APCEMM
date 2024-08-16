/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* ReadJRates Header File                                           */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 1/16/2019                                 */
/* File                 : ReadJRates.hpp                            */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef READJRATES_H_INCLUDED
#define READJRATES_H_INCLUDED

#include <cstring>
#include <netcdf>

using namespace netCDF;
using namespace netCDF::exceptions;

void ReadJRates( const char* ROOTDIR,                          \
                 const unsigned int MM, const unsigned int DD, \
                 const double LON, const double LAT,           \
                 const double P_hPa,                           \
                 double NOON_JRATES[] );

#endif /* READJREATES_H_INCLUDED */
