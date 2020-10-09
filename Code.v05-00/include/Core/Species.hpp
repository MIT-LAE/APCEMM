/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Species Header File                                              */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : Species.hpp                               */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef SPECIES_H_INCLUDED
#define SPECIES_H_INCLUDED

#include <iostream> 
#include <vector>

#include "KPP/KPP_Parameters.h"
#include "Core/Interface.hpp"
#include "Core/Structure.hpp"
#include "Util/ForwardDecl.hpp"

class SpeciesArray
{

    public:

        SpeciesArray( );

        SpeciesArray( const UInt nRing, const UInt nTime, const bool halfRing = 0 );

        SpeciesArray( const SpeciesArray &sp );

        SpeciesArray& operator=( const SpeciesArray &sp );

        SpeciesArray& operator+( const SpeciesArray &sp );

        SpeciesArray& operator-( const SpeciesArray &sp );

        ~SpeciesArray( );

        void FillIn( Solution &Data,             \
                     const Vector_3D &weights,   \
                     UInt nCounter );

        void FillIn( UInt iTime, UInt iRing );

        void getData( UInt iTime, UInt iRing );

        Vector_1D RingAverage( const Vector_1D ringArea, \
                               const RealDouble totArea, \
                               const UInt iNt ) const;

        Vector_2D RingAverage( const Vector_1D ringArea, \
                               const RealDouble totArea ) const;

        UInt getnRing() const { return nRing; }
        UInt getnTime() const { return nTime; }
        bool gethalfRing() const { return halfRing; }

        /* Reactive species */
        Vector_3D Species;

        /* Aerosols */
        Vector_2D sootDens, sootRadi, sootArea, \
                  iceDens , iceRadi , iceArea,  \
                  sulfDens, sulfRadi, sulfArea;

    protected:

        UInt nRing;
        UInt nTime;
        bool halfRing;

    private:

};

#endif /* SPECIES_H_INCLUDED */

