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
#include <complex>

/* List typedefs */
typedef unsigned int UInt;
typedef unsigned long ULong;
typedef double RealDouble;
typedef float RealFloat;
typedef long double LongDouble;
typedef std::complex<double> ComplexDouble;
typedef std::complex<float> ComplexFloat;
typedef std::complex<long double> ComplexlDouble;

/* List vector typedefs */
typedef std::vector<RealDouble> Vector_1D;
typedef std::vector<Vector_1D> Vector_2D;
typedef std::vector<Vector_2D> Vector_3D;
typedef std::vector<Vector_3D> Vector_4D;
typedef std::vector<RealFloat> Vector_1Df;
typedef std::vector<Vector_1Df> Vector_2Df;
typedef std::vector<Vector_2Df> Vector_3Df;
typedef std::vector<Vector_3Df> Vector_4Df;
typedef std::vector<LongDouble> Vector_1Dl;
typedef std::vector<Vector_1Dl> Vector_2Dl;
typedef std::vector<Vector_2Dl> Vector_3Dl;
typedef std::vector<Vector_3Dl> Vector_4Dl;
typedef std::vector<ComplexDouble> Vector_1Dc;
typedef std::vector<Vector_1Dc> Vector_2Dc;
typedef std::vector<Vector_2Dc> Vector_3Dc;
typedef std::vector<Vector_3Dc> Vector_4Dc;
typedef std::vector<ComplexFloat> Vector_1Dcf;
typedef std::vector<Vector_1Dcf> Vector_2Dcf;
typedef std::vector<Vector_2Dcf> Vector_3Dcf;
typedef std::vector<Vector_3Dcf> Vector_4Dcf;
typedef std::vector<ComplexlDouble> Vector_1Dcl;
typedef std::vector<Vector_1Dcl> Vector_2Dcl;
typedef std::vector<Vector_2Dcl> Vector_3Dcl;
typedef std::vector<Vector_3Dcl> Vector_4Dcl;

#endif /* FORWARDDECL_H_INCLUDED */
