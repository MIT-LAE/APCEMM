/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Advection Program File                                           */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : Advection.cpp                             */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <iostream>
#include <cmath>
#include "Core/Parameters.hpp"

void AdvGlobal( const double time, const double T_UPDRAFT, \
                const double V_UPDRAFT,                    \
                double &v_x, double &v_y,                  \
                double &dTrav_x, double &dTrav_y )
{

    /* AdvGlobal:
     * Computes domain advection parameters
     *
     * INPUTS:
     * (double) time     : current time since simulation started in [s]
     * (double) T_UPDRAFT: updraft timescale in [s]
     * (double) V_UPDRAFT: initial upward velocity in [s]
     *
     * OUTPUTS:
     * (double) v_x: current horizontal velocity at time in [m/s]
     * (double) v_y: current vertical velocity at time in [m/s]
     * (double) dTrav_x: vertical distance traveled since simulation started [m]
     * (double) dTrav_y: vertical distance traveled since simulation started [m]
     */


    v_x = VX;
    dTrav_x = VX * time;

    if ( SYNPROF == 0 ) {

        /* Velocity step */
        if ( time < T_UPDRAFT ) {
            v_y = V_UPDRAFT;
            dTrav_y = v_y * time;
        }
        else
            v_y = 0.0;

    }
    else if ( SYNPROF == 1 ) {

        /* Exponentially decaying velocity */
        v_y = V_UPDRAFT * exp ( - time / T_UPDRAFT );

        dTrav_y = V_UPDRAFT * T_UPDRAFT * ( 1.0 - exp ( - time / T_UPDRAFT ) );

    }
    else {
        std::string const currFile("AdvGlobal.cpp");
        std::cout << "ERROR: In " << currFile << ":SYNPROF set to " << SYNPROF << std::endl;
    }

} /* End of AdvGlobal */

void AdvParam( double time, double &v_x, double &v_y )
{

    /* AdvParam:
     * Computes advection parameters
     *
     * INPUTS:
     * (double) time: current time since simulation started in [s]
     *
     * OUTPUTS:
     * (double) v_x: current horizontal velocity at time in [m/s]
     * (double) v_y: current vertical velocity at time in [m/s]
     */

    v_x = 0.0;

    v_y = 0.0;


} /* End of AdvParam */

/* End of Advection.cpp */
