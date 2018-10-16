/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*                        AIrcraft Microphysics                     */
/*                              (AIM)                               */
/*                                                                  */
/* Nucleation Header File                                           */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 10/1/2018                                 */
/* File                 : Nucleation.hpp                            */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef NUCLEATION_H_INCLUDED
#define NUCLEATION_H_INCLUDED

#include <iostream>
#include <cmath>

#include "Util/ForwardsDecl.hpp"

namespace AIM
{

    RealDouble sigma( RealDouble x_m, RealDouble T );
    RealDouble rho( RealDouble x_m, RealDouble T );
    RealDouble x_star( RealDouble T, RealDouble RH, RealDouble nSulf );
    RealDouble nuclRate( RealDouble T, RealDouble x_m, RealDouble RH, RealDouble nSulf );
    RealDouble nTot( RealDouble T, RealDouble x_m, RealDouble RH, RealDouble nSulf );
    RealDouble radCluster( RealDouble x_m, RealDouble n );
    RealDouble nThresh( RealDouble T, RealDouble RH );

}

#endif /* NUCLEATION_H_INCLUDED */
