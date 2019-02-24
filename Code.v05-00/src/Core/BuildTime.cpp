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
    const bool INIT_TIMESTEP_REFINEMENT = 0;

    if ( nextTimeStep > 0.0E+00 ) {
        timeStep = nextTimeStep;
        nextTimeStep = 0.0E+00;
        return timeStep;
    }

    if (( time - tStart ) < 3600.0 && ( time - tStart ) >= 0.0) {
        if ( INIT_TIMESTEP_REFINEMENT ) {
            /* If default_TimeStep = 10 mins, then:
             * 0.0 - .5 - 1.0 - 1.5 - 2.0 - 2.5 - 3.0 - 3.5 - ... - 59.5 - 60.0 */
            timeStep = default_TimeStep / (double) 20.0; 
        } else {
            /* If default_TimeStep = 10 mins, then:
             * 0.0 - 5.0 - 10.0 - ... - 55.0 - 60.0 */
            timeStep = default_TimeStep / (double) 2.0;
        }
    } else
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
