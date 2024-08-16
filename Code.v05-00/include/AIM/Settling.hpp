/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*                        AIrcraft Microphysics                     */
/*                              (AIM)                               */
/*                                                                  */
/* Settling Header File                                             */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 10/31/2018                                */
/* File                 : Settling.hpp                              */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef SETTLING_H_INCLUDED
#define SETTLING_H_INCLUDED

#include "Util/ForwardDecl.hpp"

namespace AIM
{

    Vector_1D SettlingVelocity( const Vector_1D binCenters, const double T, const double P );

}

#endif /* SETTLING_H_INCLUDED */
