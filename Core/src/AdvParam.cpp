/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* AdvParam Program File                                            */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : AdvParam.cpp                              */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <iostream>
#include <cmath>
#include "Parameters.hpp"

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
