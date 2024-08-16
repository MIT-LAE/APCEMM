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

#include <iostream>
#include "Util/PhysFunction.hpp"

#include "AIM/Settling.hpp"

namespace AIM
{

    Vector_1D SettlingVelocity( const Vector_1D binCenters, const double T, const double P )
    {
   
        /* DESCRIPTION: Computes the fall speed of a particle through air accounting for slip flow 
         * correction factor */

        /* INPUTS:
         * - Vector_1D binCenters : Centers of each bin in m
         * - double T         : Temperature in K
         * - double P         : Pressure in Pa
         *
         * OUTPUT:
         * - Vector_1D vFall      : Settling velocity in m/s
         *
         * The settling velocity is computed as follows: 
         * v = 2 * g * ( \rho_p - \rho_a ) * r^2 / ( 9 \mu ) * G( Kn )
         * where \mu is the dynamic viscosity of air, Kn the Knudsen number
         * and G the slip correction factor. */

        const bool DBG    = 0;
        const bool Stokes = 1;

        Vector_1D vFall( binCenters.size(), 0.0E+00 );

        if ( DBG ) {
            std::cout << " Acceleration due to gravity: " << physConst::g << " [m/s^2]\n";
            std::cout << " Density of air: " << physFunc::rhoAir( T, P ) << " [kg/m^3]\n";
            std::cout << " Dynamic viscosity: " << physFunc::dynVisc( T ) << " [kg/(m.s)]\n";
            std::cout << " Minimum bin radius: " << binCenters[0] << " [m]\n";
            std::cout << " Maximum bin radius: " << binCenters.back() << " [m]\n";
        }

        if ( Stokes ) {

            /* Terminal fall velocity computed according to Stokes' law:
             * v = 2 g / ( 9 * \mu ) * ( \rho_i - \rho_air ) * R^2 * C
             *
             * where C represents the Cunningham correction factor. */

            for ( UInt iBin = 0; iBin < binCenters.size(); iBin++ ) {
                vFall[iBin] = 2.0E+00 * physConst::g * ( physConst::RHO_ICE - physFunc::rhoAir( T, P ) ) * \
                              binCenters[iBin] * binCenters[iBin] /\
                              ( 9.0E+00 * physFunc::dynVisc( T ) ) * \
                              physFunc::slip_flowCorrection( physFunc::Kn( binCenters[iBin], T, P ) );
            }
        } else {

            /* Terminal fall velocity taken from Sölch and Kärcher, 2010
             * Sölch, Ingo, and Bernd Kärcher. "A large‐eddy model for
             * cirrus clouds with explicit aerosol and ice microphysics
             * and Lagrangian ice particle tracking." Quarterly Journal
             * of the Royal Meteorological Society 136.653 (2010): 2074-2093.*/

            const bool Heymsfield   = 0;

            const double a0     = 1.70E-03;
            const double b0     = 8.00E-01;
            const double C0     = 6.00E-01;
            const double delta0 = 5.83E+00;
            const double C1     = 4.00E+00 / ( delta0 * delta0 * sqrt( C0 ) );
            const double C2     = 2.50E-01 * delta0 * delta0;
            const double eta    = physFunc::dynVisc(T);
            const double eta2   = eta * eta;
            const double rhoA   = physFunc::rhoAir(T, P);

            double X     = 0.0E+00;
            double mi_Ai = 0.0E+00;
            double Re    = 0.0E+00;

            for ( UInt iBin = 0; iBin < binCenters.size(); iBin++ ) {

                if ( Heymsfield ) {
                    mi_Ai = 2.28E-02 * pow( 2.0E+00 * binCenters[iBin], 0.59 );
                } else {
                    mi_Ai = binCenters[iBin] * physConst::RHO_ICE / 3.0E+00;
                }

                X = 8.0E+00 * physConst::g * rhoA / eta2 * binCenters[iBin] * binCenters[iBin] * mi_Ai;

                Re = C2 * pow(sqrt(1.0E+00 + C1 * sqrt(X)) - 1.0E+00, 2.0) - a0 * pow( X, b0 );

                vFall[iBin] = Re * eta / ( rhoA * 2.0E+00 * binCenters[iBin] );

            }
        }

        if ( DBG ) {
            for ( UInt iBin = 0; iBin < binCenters.size(); iBin++ ) {
                std::cout << " Bin radius: " << binCenters[iBin] << " [m] => Corr. factor: " << physFunc::slip_flowCorrection( physFunc::Kn( binCenters[iBin], T, P ) ) << " [-]\n";
                std::cout << " Bin radius: " << binCenters[iBin] << " [m] => Fall speed  : " << vFall[iBin] << " [m/s]\n";
            }
        }

        return vFall;

    }

}


/* End of Settling.cpp */
