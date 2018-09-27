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

/* End of PhysFunction.cpp */

