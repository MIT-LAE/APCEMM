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
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <iostream>
#include <cmath>
#include "Parameters.hpp"

void DiffParam( double time, double &d_x, double &d_y )
{

    /* DiffParam:
     * Computes diffusion parameters
     *
     * INPUTS:
     * (double) time: current time since simulation started in [s]
     *
     * OUTPUTS:
     * (double) d_x: current horizontal diffusion coefficient at time in [m^2/s]
     * (double) d_y: current vertical diffusion coefficient at time in [m^2/s]
     */
    
     if ( DPROF == 0 ) {
        if ( time <= tH0 )
            d_x = DH0;
        else
            d_x = DH;
        if ( time <= tV0 )
            d_y = DV0;
        else
            d_y = DV0;
    }
    else if ( DPROF == 1 ) {
        d_x = std::max( DH, DH0 + (DH - DH0) * time / tH0 );
        d_y = std::max( DV, DV0 + (DV - DV0) * time / tV0 );
    }
    else if (DPROF == 2 ) {
        d_x = DH + (DH0 - DH) * exp( -time / tH0 );
        d_y = DV + (DV0 - DV) * exp( -time / tV0 );
    }
    else {
        std::string const currFile("DiffParam.cpp");
        std::cout << "ERROR: In " << currFile << ": DPROF set to " << DPROF << "\n";
    }

    if ( d_x < 0.0 ) {
        std::cout << "d_x is negative: d_x = " << d_x << " [m^2/s]" << "\n";
        std::cout << "Setting d_x to 0.0" << "\n";
        d_x = 0.0;
    }

    if ( d_y < 0.0 ) {
        std::cout << "d_y is negative: d_y = " << d_y << " [m^2/s]" << "\n";
        std::cout << "Setting d_y to 0.0" << "\n";
        d_y = 0.0;
    }


} /* End of DiffParam */
