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
#include "Core/Input_Mod.hpp"
#include "Core/Input.hpp"
#include "Core/Species.hpp"
#include "Core/Ambient.hpp"
#include "Core/Cluster.hpp"

namespace output
{
    
    static const int SAVE_SUCCESS = 1;
    static const int SAVE_FAILURE = 0;

    int Write( const char* outFile,                                              \
               const OptInput &Input_Opt,                                        \
               const std::vector<int> speciesIndices,                            \
               const SpeciesArray &ringSpecies, const Ambient &ambientData,      \
               const Cluster &ringCluster, const std::vector<double> &timeArray, \
               const Input &input,                                               \
               const double &airDens, const double &relHumidity_i,               \
               const double &sunRise, const double &sunSet,                      \
               const Vector_3D &plumeRates, const Vector_2D &ambientRates );
    int Write_MicroPhys( const char* outputFile, \
                         const std::vector<std::vector<std::vector<std::vector<double>>>> &output_MicroPhys, \
                         const std::vector<double> &timeArray, const std::vector<double> &binCenters,        \
                         const std::vector<double> &horizDim, const std::vector<double> &verticDim,          \
                         const double temperature_K, const double pressure_Pa, const double lapseRate,       \
                         const double relHumidity_w, const double relHumidity_i );
    int Write_Adjoint( const char* outputFile,                                       \
                       const std::vector<int> speciesIndices,                        \
                       const SpeciesArray &ringSpecies, const Ambient &ambientData,  \
                       const Ambient &adjointData,                                   \
                       const std::vector<double> &ringArea, const double totArea,    \
                       const std::vector<double> &timeArray,                         \
                       const Input &input,                                           \
                       const double &airDens, const double &relHumidity_i );
    int Write_Box( const char* outputFile,               \
                   const std::vector<int> speciesIndices,\
                   const Ambient &boxData,               \
                   const std::vector<double> &timeArray, \
                   const Input &input,                   \
                   const double &airDens,                \
                   const double &relHumidity_i );

}

#endif /* SAVE_H_INCLUDED */
