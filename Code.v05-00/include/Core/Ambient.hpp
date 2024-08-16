/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Ambient Header File                                              */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : Ambient.hpp                               */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef AMBIENT_H_INCLUDED
#define AMBIENT_H_INCLUDED

#include <iostream>
#include <vector>

#include "KPP/KPP_Parameters.h"
#include "Util/ForwardDecl.hpp"
#include "Core/Parameters.hpp"

class Ambient
{
    public:

        Ambient( );
        Ambient( UInt nTime, Vector_1D ambientVector, \
                 Vector_2D aerVector, Vector_1D liqVector );
        ~Ambient( );
        Ambient( const Ambient &a );
        Ambient& operator=( const Ambient &a );
        void getData(  RealDouble aerArray[][2], UInt iTime ) const;
        void FillIn( UInt iTime );
        UInt getnTime() const;

        Vector_2D Species;

        /* */
        Vector_1D NIT;

        /* Tracers */
        Vector_1D SO4T;

        /* Aerosols */
        Vector_1D sootDens, sootRadi, sootArea, \
                            iceDens , iceRadi , iceArea , \
                            sulfDens, sulfRadi, sulfArea;

        /* Cosigne solar zenight angle */
        Vector_1D cosSZA;

    protected:

        UInt nTime;

    private:

};


#endif /* AMBIENT_H_INCLUDED */


