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

#include "Parameters.hpp"
#include "Interface.hpp"
#include "Output.hpp"
#include "Species.hpp"
#include "Ambient.hpp"
#include "Cluster.hpp"
#include "FileHandler.hpp"
#include "Util.hpp"

namespace output
{
    
    static const int SAVE_SUCCESS = 1;
    static const int SAVE_FAILURE = 0;

    int Write( SpeciesArray &ringSpecies, Ambient ambientData, Cluster &ringCluster, std::vector<double> &timeArray, \
               double &temperature_K, double &pressure_Pa, double &airDens, double &relHumidity_w, double &relHumidity_i, \
               double &longitude_deg, double &latitude_deg, double &sunRise, double &sunSet );

}
#endif /* SAVE_H_INCLUDED */
