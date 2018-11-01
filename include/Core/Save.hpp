/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Save Header File                                                 */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : Save.hpp                                  */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef SAVE_H_INCLUDED
#define SAVE_H_INCLUDED

#include <vector>

#include "Core/Parameters.hpp"
#include "Core/Interface.hpp"
#include "Core/Output.hpp"
#include "Core/Species.hpp"
#include "Core/Ambient.hpp"
#include "Core/Cluster.hpp"
#include "Core/FileHandler.hpp"
#include "Core/Util.hpp"

namespace output
{
    
    static const int SAVE_SUCCESS = 1;
    static const int SAVE_FAILURE = 0;

    int Write( const SpeciesArray &ringSpecies, const Ambient ambientData, const Cluster &ringCluster, \
               const std::vector<double> &timeArray, const double &temperature_K, const double &pressure_Pa, \
               const double &airDens, const double &relHumidity_w, const double &relHumidity_i, \
               const double &longitude_deg, const double &latitude_deg, const double &sunRise, const double &sunSet );
    int Write_MicroPhys( const std::vector<std::vector<std::vector<std::vector<double>>>> &output_MicroPhys, \
                         const std::vector<double> &timeArray, const std::vector<double> &binCenters, \
                         const std::vector<double> &horizDim, const std::vector<double> &verticDim, \
                         const double temperature_K, const double pressure_Pa, const double lapseRate,
                         const double relHumidity_w, const double relHumidity_i );

}

#endif /* SAVE_H_INCLUDED */
