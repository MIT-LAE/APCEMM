/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Ambient Program File                                             */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : Ambient.cpp                               */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "Core/Ambient.hpp"

Ambient::Ambient( )
{

    /* Default Constructor */

} /* End of Ambient::Ambient */

Ambient::Ambient( UInt nTime_, Vector_1D ambientVector, \
                  Vector_2D aerVector, Vector_1D liqVector )
{

    /* Constructor */

    UInt N = 0;

    nTime = nTime_;


    for ( N = 0; N < NSPECREACT; N++ )
        Species.push_back( Vector_1D( nTime, ambientVector[N] ) );

    NIT.assign   ( nTime, 0.0 );

    SO4T.assign  ( nTime, ambientVector[ind_SO4] + liqVector[0] );
    
    sootDens.assign( nTime, aerVector[  0][0] );
    sootRadi.assign( nTime, aerVector[  0][1] );
    sootArea.assign( nTime, aerVector[  0][2] );
    iceDens.assign ( nTime, aerVector[  1][0] );
    iceRadi.assign ( nTime, aerVector[  1][1] );
    iceArea.assign ( nTime, aerVector[  1][2] );
    sulfDens.assign( nTime, aerVector[  2][0] );
    sulfRadi.assign( nTime, aerVector[  2][1] );
    sulfArea.assign( nTime, aerVector[  1][2] );

    cosSZA.assign( nTime - 1, 0.0 );

} /* End of Ambient::Ambient */

Ambient::~Ambient( )
{

    /* Destructor */

} /* End of Ambient::~Ambient */ 


Ambient::Ambient( const Ambient &a )
{

    /* Copy size */
    nTime = a.nTime;
    
    /* Copy arrays */
    Species  = a.Species;

    NIT      = a.NIT      ;
    SO4T     = a.SO4T     ;

    sootDens = a.sootDens ;
    sootRadi = a.sootRadi ;
    iceDens  = a.iceDens  ;
    sulfDens = a.sulfDens ;
    sulfRadi = a.sulfRadi ;

    cosSZA = a.cosSZA;

} /* End of Ambient::Ambient */

Ambient& Ambient::operator=( const Ambient &a )
{

    if ( &a == this )
        return *this;

    /* Assign size */
    nTime = a.nTime;
    
    /* Assign arrays */
    Species  = a.Species;

    NIT      = a.NIT      ;
    SO4T     = a.SO4T     ;

    sootDens = a.sootDens ;
    sootRadi = a.sootRadi ;
    iceDens  = a.iceDens  ;
    sulfDens = a.sulfDens ;
    sulfRadi = a.sulfRadi ;

    cosSZA = a.cosSZA;

    return *this;

} /* End of Ambient::operator= */

void Ambient::getData( RealDouble aerArray[][2], UInt iTime ) const
{

    UInt N = 0;

    for ( N = 0; N < NVAR; N++ )
        VAR[N] = Species[N][iTime];

    for ( N = 0; N < NFIX; N++ )
        FIX[N] = Species[NVAR+N][iTime];

    aerArray[  0][0] = sootDens[iTime];
    aerArray[  0][1] = sootRadi[iTime];
    aerArray[  1][0] = iceDens[iTime];
    aerArray[  1][1] = iceRadi[iTime];
    aerArray[  2][0] = sulfDens[iTime];
    aerArray[  2][1] = sulfRadi[iTime];

    /* Ensure positiveness */
    for ( N = 0; N < NVAR; N++ ) {
        if ( VAR[N] <= 0.0 ) {
            VAR[N] = 1.0E-50;
        }
    }

} /* End of Ambient::getData */

void Ambient::FillIn( UInt iTime )
{

    UInt N = 0;

    /* Ensure positiveness */
    for ( N = 0; N < NVAR; N++ ) {
        if ( VAR[N] <= 0.0 ) {
            VAR[N] = 1.0E-50;
        }
    }

    for ( N = 0; N < NVAR; N++ )
        Species[N][iTime] = VAR[N];

} /* End of Ambient::FillIn */

UInt Ambient::getnTime() const
{

    return nTime;

} /* End of Ambient::getnTime */


/* End of Ambient.cpp */
