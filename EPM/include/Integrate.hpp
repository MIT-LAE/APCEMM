/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*                      Early Plume Microphysics                    */
/*                              (EPM)                               */
/*                                                                  */
/* Integrate Header File                                            */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 9/27/2018                                 */
/* File                 : Integrate.hpp                             */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef INTEGRATE_H_INCLUDED
#define INTEGRATE_H_INCLUDED

#include <iostream>
#include <cmath>
#include <vector>

#include "../../Headers/PhysConstant.hpp"
#include "../../Headers/PhysFunction.hpp"
#include "../../Core/include/Interface.hpp"
#include "../../Core/include/Parameters.hpp"
#include "../../Core/include/Monitor.hpp"
#include "../../Core/include/Aircraft.hpp"
#include "../../Core/include/Emission.hpp"
#include "../../AIM/include/Coagulation.hpp"
#include "../../AIM/include/Aerosol.hpp"

namespace EPM
{

    static const int EPM_SUCCESS = 1;
    static const int EPM_FAILURE = 0;

    /* Vortex sinking timescales, taken from Unterstrasser et al., 2008 */
    const double t_Vortex_0 = 8.0;
    const double t_Vortex_1 = 110.0;

    int Integrate( double temperature_K, double pressure_Pa, double relHumidity_w, double varArray[], double fixArray[], const Aircraft &ac, const Emission &EI );
    int RunMicrophysics( double temperature_K, double pressure_Pa, double relHumidity_w, double varArray[], double fixArray[], const Aircraft &ac, const Emission &EI, double delta_T_ad, double delta_T );
    double temp_Vortex( double time, double delta_T );



}


#endif /* INTEGRATE_H_INCLUDED */
