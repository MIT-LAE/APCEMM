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

#include <cmath>

namespace AIM
{

    double sigma( double x_m, double T );
    double rho( double x_m, double T );
    double x_star( double T, double RH, double nSulf );
    double nuclRate( double T, double x_m, double RH, double nSulf );
    double nTot( double T, double x_m, double RH, double nSulf );
    double radCluster( double x_m, double n );
    double nThresh( double T, double RH );

}

#endif /* NUCLEATION_H_INCLUDED */
