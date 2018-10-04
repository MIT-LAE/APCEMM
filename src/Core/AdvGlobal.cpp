/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* AdvGlobal Program File                                           */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : AdvGlobal.cpp                             */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <iostream>
#include <cmath>
#include "Core/Parameters.hpp"

void AdvGlobal( double time, double &v_x, double &v_y, double &dTrav_x, double &dTrav_y )
{

    /* AdvGlobal:
     * Computes domain advection parameters
     *
     * INPUTS:
     * (double) time: current time since simulation started in [s]
     *
     * OUTPUTS:
     * (double) v_x: current horizontal velocity at time in [m/s]
     * (double) v_y: current vertical velocity at time in [m/s]
     * (double) dTrav_x: vertical distance traveled since simulation started [m]
     * (double) dTrav_y: vertical distance traveled since simulation started [m]
     */


    v_x = VX;
    dTrav_x = VX * time;

    if ( SYNLIFT ) {
       
        if ( SYNPROF == 0 ) {

            /* Velocity step */
            if ( time < T_SYN ) {
                v_y = V_SYN;
                dTrav_y = v_y * time;
            }
            else {
                v_y = 0.0;
            }

        }
        else if ( SYNPROF == 1 ) {

            /* Exponentially decaying velocity */
            v_y = V_SYN * exp ( - time / T_SYN );

            dTrav_y = V_SYN * T_SYN * ( 1.0 - exp ( - time / T_SYN ) );

        }
        else {
            std::string const currFile("AdvGlobal.cpp");
            std::cout << "ERROR: In " << currFile << ":SYNPROF set to " << SYNPROF << std::endl;
        }

    }
    else {

        v_y = VY;
        dTrav_y = VY * time;

    }

} /* End of AdvGlobal */
