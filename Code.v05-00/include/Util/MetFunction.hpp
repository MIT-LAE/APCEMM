/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* MetFunction Header File                                          */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 11/2/2018                                 */
/* File                 : MetFunction.hpp                           */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef METFUNCTION_H_INCLUDED
#define METFUNCTION_H_INCLUDED

#include <iostream>
#include <cmath>
#include "ForwardDecl.hpp"
#include "PhysConstant.hpp"
#include "PhysFunction.hpp"
#include "../Core/Parameters.hpp"

#define ABS(x)   ( ((x) >=  0 ) ?(x):(-x) ) 

namespace met
{

    /* Defined meteorological functions */

    /* International Standard Atmosphere (ISA) */

    const RealDouble ISA_LAPSERATE = 6.5E-03;
    const RealDouble ISA_HTS = 1.1E+04;
    const RealDouble ISA_HTP = 2.0E+04;
    const RealDouble ISA_RHO0 = 1.225E+00;
    const RealDouble ISA_P0 = 101325;
    const RealDouble ISA_T0 = 288.15;

    void ISA( const RealDouble &z, RealDouble &pressure, \
              RealDouble &temperature );
    void ISA( const Vector_1D &z, Vector_1D &pressure, \
              Vector_1D &temperature );
    void ISA( const RealDouble &z, RealDouble &pressure );
    void ISA( const Vector_1D &z, Vector_1D &pressure );
    void ISA_pAlt( RealDouble &z, const RealDouble &pressure );
    void ISA_pAlt( Vector_1D &z, const Vector_1D &pressure );


    RealDouble ComputeLapseRate( const RealDouble TEMP, const RealDouble RHi, \
                                 const RealDouble DEPTH );
    float linearInterp( float xq[], float yq[], const float &x, size_t xq_size);
    UInt nearestNeighbor( float xq[], const float &x, size_t xq_size);
    UInt nearestNeighbor( Vector_1D xq, const float &x );
    RealDouble satdepth_calc( float RHw[], float T[], float alt[], UInt iFlight, UInt var_length );

}

#endif /* METFUNCTION_H_INCLUDED */
