/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* ForwardDecl Header File                                          */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 10/1/2018                                 */
/* File                 : ForwardDecl.hpp                           */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef FORWARDDECL_H_INCLUDED
#define FORWARDDECL_H_INCLUDED

#include <vector>

/* List typedefs */
typedef unsigned int UInt;

/* List vector typedefs */
typedef std::vector<double> Vector_1D;
typedef std::vector<Vector_1D> Vector_2D;
typedef std::vector<Vector_2D> Vector_3D;
typedef std::vector<UInt> Vector_1Dui;
typedef std::vector<Vector_1Dui> Vector_2Dui;

#endif /* FORWARDDECL_H_INCLUDED */
