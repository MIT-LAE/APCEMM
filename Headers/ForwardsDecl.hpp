/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* ForwardsDecl Header File                                         */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 10/1/2018                                 */
/* File                 : ForwardsDecl.hpp                          */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef FORWARDSDECL_H_INCLUDED
#define FORWARDSDECL_H_INCLUDED

#include <vector>

/* List typedefs */
typedef unsigned int UInt;
typedef unsigned long ULong;
typedef double RealDouble;

/* List vector typedefs */
typedef std::vector<RealDouble> Vector_1D;
typedef std::vector<Vector_1D> Vector_2D;
typedef std::vector<Vector_2D> Vector_3D;
typedef std::vector<Vector_3D> Vector_4D;


#endif /* FORWARDSDECL_H_INCLUDED */
