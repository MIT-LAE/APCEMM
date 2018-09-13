/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* UpdateTime Program File                                          */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : UpdateTime.cpp                            */
/* Working directory    : /home/fritzt/APCEMM-SourceCode            */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <iostream>
#include <cmath>

#include "Parameters.hpp"

double UpdateTime( double time, double tStart, double sunRise, double sunSet )
{

    double default_TimeStep = DT;
    double timeStep;

    if (( time - tStart ) < 3600.0 && ( time - tStart ) >= 0.0) {
        if (( time - tStart ) < 800.0)
            timeStep = (double) 50.0;
        else if (( time - tStart ) >= 800.0 && ( time - tStart) < 1200.0)
            timeStep = (double) 100.0;
        else if (( time - tStart ) >= 1200.0 && ( time - tStart) < 3600.0)
            timeStep = (double) 300.0;
        else
            timeStep = default_TimeStep;
    }
    else
        timeStep = default_TimeStep;

    if ( (std::fmod((time),(24.0*3600.0)) < (sunSet)) && (std::fmod((time + timeStep),(24.0*3600.0)) > sunSet) )
        timeStep = std::max( sunSet - std::fmod((time),(24.0*3600.0)), 1.0 );
    if ( (std::fmod((time),(24.0*3600.0)) < (sunRise)) && (std::fmod((time + timeStep),(24.0*3600.0)) > sunRise) )
        timeStep = std::max( sunRise - std::fmod((time),(24.0*3600.0)), 1.0 );

    return timeStep;
} /* End of Updatetime */
