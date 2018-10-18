/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* SZA Header File                                                  */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 10/18/2018                                */
/* File                 : SZA.hpp                                   */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef SZA_H_INCLUDED
#define SZA_H_INCLUDED

#include <iostream>
#include <cmath>

#include "Util/PhysConstant.hpp"

class SZA
{

    public:
        
        SZA( const double latitude, const unsigned int dayGMT );
        ~SZA( );

        void Update( const double solarTime );

        const double latitude;
        const unsigned int dayGMT;

        double sunRise;
        double sunSet;

        double CSZA;

    private:
        
        double const A0 = 0.006918;
        double const A1 = 0.399912;
        double const A2 = 0.006758;
        double const A3 = 0.002697;
        double const B1 = 0.070257;
        double const B2 = 0.000907;
        double const B3 = 0.000148;

        double DEC;
        double SINLAT;
        double COSLAT;
        double SINDEC;
        double COSDEC;

};

#endif /* SZA_H_INCLUDED */
