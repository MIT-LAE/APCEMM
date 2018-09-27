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

static const double PI = 3.141592653589793238460; /* \pi */

physFunc::physFunc()
{

    /* Default Constructor */

} /* End of physFunc::physFunc */

physFunc::~physFunc()
{

    /* Destructor */

} /* End of physFunc::~physFunc */

double physFunc::pSat_H2Ol( double T )
{
    
    /* Returns water liquid saturation pressure in Pascal. */
    
    return 100.0 * exp( - 6096.9385 / T \
                        + 16.635794 \
                        - 0.02711193 * T \
                        + 1.673952E-5 * T * T \
                        + 2.433502 * log( T ) );

} /* End of physFunc::pSat_H2Ol */

double physFunc::pSat_H2Os( double T )
{
    
    /* Returns water solid saturation pressure in Pascal. */
    
    return 100.0 * exp( - 6024.5282 / T \
                        + 24.7219 \
                        + 0.010613868 * T \
                        - 1.3198825E-5 * T * T \
                        - 0.49382577 * log( T ) );

} /* End of physFunc::pSat_H2Os */

double physFunc::pSat_H2SO4( double T )
{
    
    /* Returns sulfuric acid saturation pressure in Pascal. */
    
    return 100.0 * exp( + 23.1885 \
                        + 10156.0 * ( \
                                      - 1.0 / T\
                                      + 0.38/(545.0)*( \
                                                       + 1.0 \
                                                       + log( 360.0 * 1.0 / T ) \
                                                       - 360.0 * 1.0 / T ) ) );

} /* End of physFunc::pSat_H2SO4 */

double physFunc::pSat_HNO3( double T , double PPH2O ) 
{
    
    /* Returns nitric acid saturation pressure in Pascal. */
    
    return ( ATM / 760.0 ) * pow( 10.0, ( ( ( - 2.7836 - 0.00088 * T ) * log10( PPH2O * ( 760.0 / ATM ) ) ) + ( 38.9855 - 11397.0 / T + 0.009179 * T ) ) );  

} /* End of physFunc::pSat_HNO3 */

double physFunc::rhoAir( double T, double P )
{

    /* Returns the density of air in kg/m^3 */

    return P / ( R_Air * T);

} /* End of physFunc::airDens */ 

double physFunc::dynVisc( double T )
{

    /* Returns the dynamic viscosity of air in kg/m/s */

    return 1.8325E-05 * ( 416.16 / ( T + 120.0) ) * pow( T / 296.16, 1.5 );

} /* End of physFunc::dynVisc */

double physFunc::kinVisc( double T, double P )
{

    /* Returns the kinematic viscosity of air in m^2/s */

    return dynVisc( T ) / rhoAir( T, P );

}

double physFunc::thermalSpeed( double T )
{

    /* Returns the thermal speed of an air molecule in m/s */

    return sqrt( 8.0 * kB * T / ( PI * M_Air ) );

}

double physFunc::lambda( double T, double P )
{

    /* Returns the mean free path of an air molecule in m */

    return 2.0 * kinVisc( T, P ) / thermalSpeed( T );

} /* End of physFunc::lambda */

double physFunc::Kn( double r, double T, double P )
{

    /* Returns the Knudsen number of air in - */

    return lambda( T, P ) / r;

} /* End of physFunc::Kn */

double physFunc::partDiffCoef( double r, double T, double P )
{

    /* Returns the particle diffusion coefficient in m^2/s */

    return kB * T / ( 6.0 * PI * dynVisc( T ) * r ) * slip_flowCorrection( Kn ( r, T, P ) );

} /* End of physFunc::partDiffCoef */

double physFunc::slip_flowCorrection( double Kn )
{

    /* Returns the Cunningham slip-flow correction */

    static const double A = 1.249;
    static const double B = 0.42;
    static const double C = 0.87;

    return 1 + Kn * ( A + B * exp( - C / Kn ) );

} /* End of physFunc::slip_flowCorrection */

/* End of PhysFunction.cpp */

