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

double UpdateTime( double time, const double tStart, \
                   const double sunRise, const double sunSet, \
                   const double DYN_DT, double& nextTimeStep );

std::vector<double> BuildTime( const double tStart, const double tEnd, \
                               const double sunRise, const double sunSet, \
                               const double DYN_DT )
{

    unsigned int nT = 0;

    std::vector<double> timeArray;
    double time = tStart;
    double timeStep, nextTimeStep;

    timeStep = 0.0E+00;
    nextTimeStep = 0.0E+00;

    while ( time < tEnd ) {

        /* For debug purposes, uncomment the next line */
        // std::cout << time/3600 << ", " << timeStep << "  ";

        timeArray.push_back( time );
        timeStep = UpdateTime( time, tStart, sunRise, sunSet, \
                               DYN_DT, nextTimeStep );
        time += std::min( timeStep, std::abs( ( tEnd - time ) ) );
        nT++;
    }

    timeArray.push_back( time );
    return timeArray;

} /* End of BuildTime */

double UpdateTime( double time, const double tStart, \
                   const double sunRise, const double sunSet, \
                   const double DYN_DT, double &nextTimeStep )
{

    const double default_TimeStep = DYN_DT;
    double timeStep;
    double currTimeStep;
    bool sunCorrection = 0;

    if ( nextTimeStep > 0.0E+00 ) {
        timeStep = nextTimeStep;
        nextTimeStep = 0.0E+00;
        return timeStep;
    }

    if (( time - tStart ) < 3600.0 && ( time - tStart ) >= 0.0) {
        if ( ( time - tStart ) < default_TimeStep )
            timeStep = default_TimeStep / (double) 10.0;
        else if ( ( time - tStart ) >= default_TimeStep && \
                  ( time - tStart ) < 2*default_TimeStep )
            timeStep = default_TimeStep / (double) 5.0;
        else if (( time - tStart ) >= 2*default_TimeStep && \
                 ( time - tStart) < 6*default_TimeStep )
            timeStep = default_TimeStep / (double) 2.0;
        else
            timeStep = default_TimeStep;
    }
    else
        timeStep = default_TimeStep;

    /* Add a small bit to ensure that we get passed the sunrise/set */
    currTimeStep = timeStep;
    if ( (std::fmod((time),(24.0*3600.0)) < (sunSet)) \
      && (std::fmod((time + timeStep),(24.0*3600.0)) > sunSet) ) {
        timeStep = sunSet - std::fmod((time),(24.0*3600.0)) + 1.00E-00;
        sunCorrection = 1;
    }
    if ( (std::fmod((time),(24.0*3600.0)) < (sunRise)) \
      && (std::fmod((time + timeStep),(24.0*3600.0)) > sunRise) ) {
        timeStep = sunRise - std::fmod((time),(24.0*3600.0)) + 1.00E-00;
        sunCorrection = 1;
    }

    if ( sunCorrection == 1 )
        nextTimeStep = currTimeStep - timeStep;

    return timeStep;

} /* End of Updatetime */

/* End of BuildTime.cpp */
