/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* PhysFunction Program File                                        */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : PhysFunction.cpp                          */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "PhysFunction.hpp"

namespace physFunc
{

    double pSat_H2Ol( double T )
    {
        
        /* Returns water liquid saturation pressure in Pascal. */
        
        return 100.0 * exp( - 6096.9385 / T \
                            + 16.635794 \
                            - 0.02711193 * T \
                            + 1.673952E-5 * T * T \
                            + 2.433502 * log( T ) );

    } /* End of pSat_H2Ol */

    double pSat_H2Os( double T )
    {
        
        /* Returns water solid saturation pressure in Pascal. */
        
        return 100.0 * exp( - 6024.5282 / T \
                            + 24.7219 \
                            + 0.010613868 * T \
                            - 1.3198825E-5 * T * T \
                            - 0.49382577 * log( T ) );

    } /* End of pSat_H2Os */

    double pSat_H2SO4( double T )
    {
        
        /* Returns sulfuric acid saturation pressure in Pascal. */
        
        return 100.0 * exp( + 23.1885 \
                            + 10156.0 * ( \
                                          - 1.0 / T\
                                          + 0.38/(545.0)*( \
                                                           + 1.0 \
                                                           + log( 360.0 * 1.0 / T ) \
                                                           - 360.0 * 1.0 / T ) ) );

    } /* End of pSat_H2SO4 */

    double pSat_HNO3( double T , double PPH2O ) 
    {
        
        /* Returns nitric acid saturation pressure in Pascal. */
        
        return ( physConst::ATM / 760.0 ) * pow( 10.0, ( ( ( - 2.7836 - 0.00088 * T ) * log10( PPH2O * ( 760.0 / physConst::ATM ) ) ) + ( 38.9855 - 11397.0 / T + 0.009179 * T ) ) );  

    } /* End of pSat_HNO3 */

    double rhoAir( double T, double P )
    {

        /* Returns the density of air in kg/m^3 */

        return P / ( physConst::R_Air * T);

    } /* End of airDens */ 

    double dynVisc( double T )
    {

        /* Returns the dynamic viscosity of air in kg/m/s */

        return 1.8325E-05 * ( 416.16 / ( T + 120.0) ) * pow( T / 296.16, 1.5 );

    } /* End of dynVisc */

    double kinVisc( double T, double P )
    {

        /* Returns the kinematic viscosity of air in m^2/s */

        return dynVisc( T ) / rhoAir( T, P );

    } /* End of kinVisc */

    double thermalSpeed( double T, double m )
    {

        /* Returns the thermal speed of an air molecule/particle in m/s */

        return sqrt( 8.0 * physConst::kB * T / ( physConst::PI * m ) );

    } /* End of thermalSpeed */

    double lambda( double T, double P )
    {

        /* Returns the mean free path of an air molecule in m */

        return 2.0 * kinVisc( T, P ) / thermalSpeed( T );

    } /* End of lambda */

    double mass_sphere( double r, double rho )
    {
    
        /* Returns the mass of a sphere in kg */

        return 4.0 / double(3.0) * physConst::PI * rho * r * r * r;

    } /* End of mass_sphere */

    double vFall( double r, double rho, double T, double P )
    {

        /* Returns the terminal fall speed of a spherical particle in m/s */

        return 2.0 * rho * physConst::g * r * r / ( 9.0 * dynVisc( T ) ) * slip_flowCorrection( Kn( r, T, P ) ); 

    } /* End of vFall */
    
    double Kn( double r, double T, double P )
    {

        /* Returns the Knudsen number of air in - */

        return lambda( T, P ) / r;

    } /* End of Kn */

    double partDiffCoef( double r, double T, double P )
    {

        /* Returns the particle diffusion coefficient in m^2/s */

        return physConst::kB * T / ( 6.0 * physConst::PI * dynVisc( T ) * r ) * slip_flowCorrection( Kn ( r, T, P ) );

    } /* End of partDiffCoef */

    double slip_flowCorrection( double Kn )
    {

        /* Returns the Cunningham slip-flow correction */

        static const double A = 1.249;
        static const double B = 0.42;
        static const double C = 0.87;

        return 1 + Kn * ( A + B * exp( - C / Kn ) );

    } /* End of slip_flowCorrection */
    
    double lambda_p( double r, double m, double T, double P )
    {

        /* Returns the particle mean path in air in m */

        return 8.0 * partDiffCoef( r, T, P ) / ( physConst::PI * thermalSpeed( T, m ) );

    } /* End of lambda_p */

    double delta_p( double r, double m, double T, double P )
    {

        /* Returns the mean distance in m from the center of a sphere reached by particles 
         * leaving the sphere's surface and traveling a distance of particle mean free 
         * path lambda_p */

        double l = lambda_p( r, m, T, P );

        return ( pow( 2.0 * r + l , 3.0 ) - pow( 4.0 * r * r + l * l , 1.5 ) ) / ( 6.0 * r * l ) - 2.0 * r;

    } /* End of delta_p */

}

/* End of PhysFunction.cpp */

