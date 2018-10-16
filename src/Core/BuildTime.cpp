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

#define ABS(x)	 ( ((x) >= 0 )  ? (x):(-x) )
#define MIN(a,b) ( ((a) <= b )  ? (a):(b)  )

double UpdateTime( double time, double tStart, double sunRise, double sunSet );

std::vector<double> BuildTime( double tStart, double tEnd, \
                               double sunRise, double sunSet )
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
        time += MIN( timeStep, ABS( ( tEnd - time ) ) );
        nT++;
    }

    timeArray.push_back((double) 0.0);
    timeArray[nT] = time;

    return timeArray;

} /* End of BuildTime */
