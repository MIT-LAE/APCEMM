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

#include <iostream>
#include <iomanip>
#include <cstring>
#include <algorithm>
#include "mat.h"

#include "Util/ForwardDecl.hpp"
#include "KPP/KPP_Parameters.h"

Vector_1D ReadJRates( const std::string ROOTDIR,                    \
                      const unsigned int MM, const unsigned int DD, \
                      const double LON, const double LAT,           \
                      const double P_hPa );

#endif /* READJREATES_H_INCLUDED */
