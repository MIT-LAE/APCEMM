/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Diffusion Program File                                           */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : Diffusion.cpp                             */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <iostream>
#include <cmath>
#include "Core/Parameters.hpp"

void DiffParam( const double time, double &d_x, double &d_y, \
                const double D_X, const double D_Y )
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
            d_x = 1.13 * D_X;
        else
            d_x = D_X;
        if ( time <= tV0 )
            d_y = 7.0 * D_Y;
        else
            d_y = D_Y;
    }
    else if ( DPROF == 1 ) {
        d_x = std::max( D_X, 1.13 * D_X + (D_X - 1.13 * D_X) * time / tH0 );
        d_y = std::max( D_Y, 7.00 * D_Y + (D_Y - 7.00 * D_Y) * time / tV0 );
    }
    else if ( DPROF == 2 ) {
        d_x = D_X + (1.13 * D_X - D_X) * exp( -time / tH0 );
        d_y = D_Y + (7.00 * D_Y - D_Y) * exp( -time / tV0 );
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

/* End of Diffusion.cpp */
