/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* SatPressures Program File                                        */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : SatPressures.cpp                          */
/* Working directory    : /home/fritzt/APCEMM-SourceCode            */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <math.h>

#define HUN	    (double)100.0
#define ONE	    (double)1.0
#define TWO	    (double)2.0
#define TEN     (double)10.0
#define ATM     (double)1.01325e+05

double pSat_H2Ol( double T )
{
    /* Returns water liquid saturation pressure in Pascal. */

    return HUN * exp( - 6096.9385 / T \
                      + 16.635794 \
                      - 0.02711193 * T \
                      + 1.673952E-5 * pow( T, TWO ) \
                      + 2.433502 * log( T ) );

} /* End of pSat_H2Ol */

double pSat_H2Os( double T )
{
    /* Returns water solid saturation pressure in Pascal. */

    return HUN * exp( - 6024.5282 / T \
                      + 24.7219 \
                      + 0.010613868 * T \
                      - 1.3198825E-5 * pow( T, TWO ) \
                      - 0.49382577 * log( T ) );

} /* End of pSat_H2Os */

double pSat_H2SO4( double T )
{
    /* Returns sulfuric acid saturation pressure in Pascal. */

    return HUN * exp( + 23.1885 \
                      + 10156.0 * ( \
                                 - ONE / T\
                                 + 0.38/(545.0)*( \
                                                 + ONE \
                                                 + log( 360.0 * ONE / T )\
                                                 - 360.0 * ONE / T ) ) );

} /* End of pSat_H2SO4 */

double pSat_HNO3( double T , double PPH2O )
{
    /* Returns nitric acid saturation pressure in Pascal. */

    return (ATM / 760.0) * pow( TEN, ( ( ( - 2.7836 - 0.00088 * T) * log10(PPH2O * (760.0 / ATM) ) ) + ( 38.9855 - 11397.0 / T + 0.009179 * T ) ) );

} /* End of pSat_HNO3 */

