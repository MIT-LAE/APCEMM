/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* BuildTime Program File                                           */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : BuildTime.cpp                             */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <iostream>
#include <vector>
#include <cmath>
#include "Core/Parameters.hpp"

double UpdateTime( double time, const double tStart, \
                   const double sunRise, const double sunSet );

std::vector<double> BuildTime( const double tStart, const double tEnd, \
                               const double sunRise, const double sunSet )
{

    unsigned int nT = 0;

    std::vector<double> timeArray;
    double time = tStart;
    double timeStep;

    while ( time < tEnd ) {
        timeArray.push_back((double) 0.0);
        timeArray[nT] = time;
        timeStep = UpdateTime( time, tStart, sunRise, sunSet );
        //std::cout << time/3600 << ", " << timeStep << std::endl;
        time += std::min( timeStep, std::abs( ( tEnd - time ) ) );
        nT++;
    }

    timeArray.push_back((double) 0.0);
    timeArray[nT] = time;

    return timeArray;

} /* End of BuildTime */

double UpdateTime( double time, const double tStart, \
                   const double sunRise, const double sunSet )
{

    const double default_TimeStep = DT;
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

    /* Add a small bit to ensure that we get passed the sunrise/set */
    if ( (std::fmod((time),(24.0*3600.0)) < (sunSet)) \
      && (std::fmod((time + timeStep),(24.0*3600.0)) > sunSet) )
        timeStep = sunSet - std::fmod((time),(24.0*3600.0)) + 1.00E-00;
    if ( (std::fmod((time),(24.0*3600.0)) < (sunRise)) \
      && (std::fmod((time + timeStep),(24.0*3600.0)) > sunRise) )
        timeStep = sunRise - std::fmod((time),(24.0*3600.0)) + 1.00E-00;

    return timeStep;

} /* End of Updatetime */

/* End of BuildTime.cpp */
