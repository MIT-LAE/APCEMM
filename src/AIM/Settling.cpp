/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*                        AIrcraft Microphysics                     */
/*                              (AIM)                               */
/*                                                                  */
/* Settling Program File                                            */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 10/31/2018                                */
/* File                 : Settling.cpp                              */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "AIM/Settling.hpp"

namespace AIM
{

    Vector_1D SettlingVelocity( const Vector_1D binCenters, const RealDouble T, const RealDouble P )
    {
   
        /* DESCRIPTION: Computes the fall speed of a particle through air accounting for slip flow 
         * correction factor */

        /* INPUTS:
         * - Vector_1D binCenters : Centers of each bin in m
         * - RealDouble T         : Temperature in K
         * - RealDouble P         : Pressure in Pa
         *
         * OUTPUT:
         * - Vector_1D vFall      : Settling velocity in m/s
         *
         * The settling velocity is computed as follows: 
         * v = 2 * g * ( \rho_p - \rho_a ) * r^2 / ( 9 \mu ) * G( Kn )
         * where \mu is the dynamic viscosity of air, Kn the Knudsen number
         * and G the slip correction factor. */

        Vector_1D vFall( binCenters.size(), 0.0E+00 );

        for ( UInt iBin = 0; iBin < binCenters.size(); iBin++ ) {
            vFall[iBin] = 2.0E+00 * physConst::g * ( physConst::RHO_ICE - physFunc::rhoAir( T, P ) ) * \
                          binCenters[iBin] * binCenters[iBin] /\
                          9.0E+00 * physFunc::dynVisc( T ) * \
                          physFunc::slip_flowCorrection( physFunc::Kn( binCenters[iBin], T, P ) );
        }

        return vFall;

    }

}


/* End of Settling.cpp */
