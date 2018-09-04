/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* DiffParam Program File                                           */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : DiffParam.cpp                             */
/* Working directory    : /home/fritzt/APCEMM-SourceCode            */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <iostream>
#include <cmath>
#include "Parameters.h"

void DiffParam( double time, double &d_H, double &d_V )
{

    if ( DPROF == 0 ) {
        if ( time <= tH0 )
            d_H = DH0;
        else
            d_H = DH;
        if ( time <= tV0 )
            d_V = DV0;
        else
            d_V = DV0;
    }
    else if ( DPROF == 1 ) {
        d_H = std::max( DH, DH0 + (DH - DH0)/tH0 * time );
        d_V = std::max( DV, DV0 + (DV - DV0)/tV0 * time );
    }
    else if (DPROF == 2 ) {
        d_H = DH + (DH - DH0)*exp( -time / tH0 );
        d_V = DV + (DV - DV0)*exp( -time / tV0 );
    }
    else {
        std::string const currFile("DiffParam.cpp");
        std::cout << "ERROR: In " << currFile << ": DPROF set to " << DPROF << std::endl;
    }

} /* End of DiffParam */
