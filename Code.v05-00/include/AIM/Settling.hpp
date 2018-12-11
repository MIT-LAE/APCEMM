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

#include <iostream>

#include "Util/ForwardDecl.hpp"
#include "Util/PhysConstant.hpp"
#include "Util/PhysFunction.hpp"

namespace AIM
{

    Vector_1D SettlingVelocity( const Vector_1D binCenters, const RealDouble T, const RealDouble P );

}

#endif /* SETTLING_H_INCLUDED */
