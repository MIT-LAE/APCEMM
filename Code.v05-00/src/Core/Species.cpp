/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Species Program File                                             */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : Species.cpp                               */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "Core/Species.hpp"

static const double ZERO = 1.00E-50;

SpeciesArray::SpeciesArray( )
{
    /* Default Constructor */

} /* End of SpeciesArray::SpeciesArray */

SpeciesArray::SpeciesArray( const UInt nRing_, const UInt nTime_, const bool halfRing_ )
{

    /* Constructor */

    nRing = nRing_;
    nTime = nTime_;
    halfRing = halfRing_;

    Vector_2D v2d( nTime, Vector_1D( nRing, 0.0E+00 ) );
    Vector_1D v1d( nRing, 0.0E+00 );

    for ( UInt N = 0; N < NSPECREACT; N++ )
        Species.push_back( v2d );

    for ( UInt i = 0; i < nTime; i++ ) {
        sootDens.push_back( v1d );
        sootRadi.push_back( v1d );
        sootArea.push_back( v1d );
        iceDens.push_back( v1d );
        iceRadi.push_back( v1d );
        iceArea.push_back( v1d );
        sulfDens.push_back( v1d );
        sulfRadi.push_back( v1d );
        sulfArea.push_back( v1d );
    }

} /* End of SpeciesArray::SpeciesArray */


SpeciesArray::SpeciesArray( const SpeciesArray &sp )
{
    
    nRing    = sp.getnRing();
    nTime    = sp.getnTime();
    halfRing = sp.gethalfRing();

    Species  = sp.Species;

    sootDens = sp.sootDens;
    sootRadi = sp.sootRadi;
    sootArea = sp.sootArea;
    iceDens  = sp.iceDens;
    iceRadi  = sp.iceRadi;
    iceArea  = sp.iceArea;
    sulfDens = sp.sulfDens;
    sulfRadi = sp.sulfRadi;
    sulfArea = sp.sulfArea;


} /* End of SpeciesArray::SpeciesArray */
    
SpeciesArray& SpeciesArray::operator=( const SpeciesArray &sp )
{

    if ( &sp == this )
        return *this;

    nRing    = sp.getnRing();
    nTime    = sp.getnTime();
    halfRing = sp.gethalfRing();

    Species  = sp.Species;

    sootDens = sp.sootDens;
    sootRadi = sp.sootRadi;
    sootArea = sp.sootArea;
    iceDens  = sp.iceDens;
    iceRadi  = sp.iceRadi;
    iceArea  = sp.iceArea;
    sulfDens = sp.sulfDens;
    sulfRadi = sp.sulfRadi;
    sulfArea = sp.sulfArea;
    return *this;

} /* End of SpeciesArray::operator= */

SpeciesArray& SpeciesArray::operator+( const SpeciesArray &sp )
{

    UInt iRing, iTime, N;

    if ( nRing != sp.getnRing() ) {
        std::cout << "Can't perform + on SpeciesArray: nRing exception: " << nRing << " != " << sp.getnRing() << std::endl;
        return *this;
    }
    
    if ( nTime != sp.getnTime() ) {
        std::cout << "Can't perform + on SpeciesArray: nTime exception: " << nTime << " != " << sp.getnTime() << std::endl;
        return *this;
    }

    for ( iRing = 0; iRing < nRing; iRing++ ) {
        for ( iTime = 0; iTime < nTime; iTime++ ) {
            for ( N = 0; N < NSPECREACT; N++ )
                Species[N][iTime][iRing] += sp.Species[N][iTime][iRing];

            sootDens[iTime][iRing] += sp.sootDens[iTime][iRing];
            sootRadi[iTime][iRing] += sp.sootRadi[iTime][iRing];
            sootArea[iTime][iRing] += sp.sootArea[iTime][iRing];
            iceDens[iTime][iRing]  += sp.iceDens[iTime][iRing];
            iceRadi[iTime][iRing]  += sp.iceRadi[iTime][iRing];
            iceArea[iTime][iRing]  += sp.iceArea[iTime][iRing];
            sulfDens[iTime][iRing] += sp.sulfDens[iTime][iRing];
            sulfRadi[iTime][iRing] += sp.sulfRadi[iTime][iRing];
            sulfArea[iTime][iRing] += sp.sulfArea[iTime][iRing];

        }
    }
    return *this;

} /* End of SpeciesArray::operator+ */

SpeciesArray& SpeciesArray::operator-( const SpeciesArray &sp )
{

    if ( nRing != sp.getnRing() ) {
        std::cout << "Can't perform + on SpeciesArray: nRing exception: " << nRing << " != " << sp.nRing << std::endl;
        return *this;
    }
    
    if ( nTime != sp.getnTime() ) {
        std::cout << "Can't perform + on SpeciesArray: nTime exception: " << nTime << " != " << sp.nTime << std::endl;
        return *this;
    }

    for ( UInt iRing = 0; iRing < nRing; iRing++ ) {
        for ( UInt iTime = 0; iTime < nTime; iTime++ ) {
            for ( UInt N = 0; N < NSPECREACT; N++ )
                Species[N][iTime][iRing] -= sp.Species[N][iTime][iRing];

            sootDens[iTime][iRing] -= sp.sootDens[iTime][iRing];
            sootRadi[iTime][iRing] -= sp.sootRadi[iTime][iRing];
            sootArea[iTime][iRing] -= sp.sootArea[iTime][iRing];
            iceDens[iTime][iRing]  -= sp.iceDens[iTime][iRing];
            iceRadi[iTime][iRing]  -= sp.iceRadi[iTime][iRing];
            iceArea[iTime][iRing]  -= sp.iceArea[iTime][iRing];
            sulfDens[iTime][iRing] -= sp.sulfDens[iTime][iRing];
            sulfRadi[iTime][iRing] -= sp.sulfRadi[iTime][iRing];
            sulfArea[iTime][iRing] -= sp.sulfArea[iTime][iRing];

        }
    }
    return *this;

} /* End of SpeciesArray::operator- */

SpeciesArray::~SpeciesArray( )
{
    /* Destructor */

} /* End of SpeciesArray::~SpeciesArray */

void SpeciesArray::FillIn( Solution &Data,             \
                           const Vector_3D &weights,   \
                           UInt nCounter )
{

    Vector_1D properties_LA( 4, 0.0E+00 );
    Vector_1D properties_PA( 4, 0.0E+00 );

    UInt jNy, iNx, iRing, N;
    double w = 0.0E+00;
    double totW = 0.0E+00;

    for ( iRing = 0; iRing < nRing; iRing++ ) {

        /* Precompute total weights */
        totW = 0.0E+00;
        for ( jNy = 0; jNy < Data.Species[0].size(); jNy++ ) {
            for ( iNx = 0; iNx < Data.Species[0][0].size(); iNx++ )
                totW += weights[iRing][jNy][iNx];
        }
        for ( jNy = 0; jNy < Data.Species[0].size(); jNy++ ) {
            for ( iNx = 0; iNx < Data.Species[0][0].size(); iNx++ ) {

                w = weights[iRing][jNy][iNx] / totW;

                for ( N = 0; N < NSPECREACT; N++ )
                    Species[N][nCounter][iRing] += Data.Species[N][jNy][iNx] * w;

                sootDens[nCounter][iRing] += Data.sootDens[jNy][iNx] * w;
                sootRadi[nCounter][iRing] += Data.sootRadi[jNy][iNx] * w;
                sootArea[nCounter][iRing] += Data.sootArea[jNy][iNx] * w;

            }
        }

        properties_LA = Data.liquidAerosol.Average( weights[iRing], \
                                                    totW );
        sulfDens[nCounter][iRing] = properties_LA[0];
        sulfRadi[nCounter][iRing] = properties_LA[1];
        sulfArea[nCounter][iRing] = properties_LA[2];

        properties_PA = Data.solidAerosol.Average( weights[iRing], \
                                                   totW );
        iceDens[nCounter][iRing] = properties_PA[0];
        iceRadi[nCounter][iRing] = properties_PA[1];
        iceArea[nCounter][iRing] = properties_PA[2];

    } 

} /* End of SpeciesArray::FillIn */

void SpeciesArray::FillIn( UInt iTime, UInt iRing, double* varSpeciesArray )
{

    /* Ensure positiveness */
    for ( UInt N = 0; N < NVAR; N++ ) {
        if ( varSpeciesArray[N] <= 0.0 ) {
            varSpeciesArray[N] = ZERO;
        }
        Species[N][iTime][iRing] = varSpeciesArray[N];
    }

} /* End of SpeciesArray::FillIn */

void SpeciesArray::getData( UInt iTime, UInt iRing, double* varSpeciesArray, double* fixSpeciesArray )
{

    for ( UInt N = 0; N < NVAR; N++ )
        varSpeciesArray[N] = Species[N][iTime][iRing];

    for ( UInt N = 0; N < NFIX; N++ )
        fixSpeciesArray[N] = Species[NVAR+N][iTime][iRing];

    /* Ensure positiveness */
    for ( UInt N = 0; N < NVAR; N++ ) {
        if ( varSpeciesArray[N] <= 0.0 ) {
            varSpeciesArray[N] = ZERO;
        }
    }

} /* End of SpeciesArray::getData */

Vector_1D SpeciesArray::RingAverage( const Vector_1D ringArea, const double totArea, \
                                     const UInt iNt ) const
{

    Vector_1D ringAverage( NVAR, 0.0E+00 );
    UInt iRing, N;
    double area;

    for ( iRing = 0; iRing < nRing; iRing++ ) {
        area = ringArea[iRing] / totArea;
        for ( N = 0; N < NVAR; N++ )
            ringAverage[N] += Species[N][iNt][iRing] * area;
    }

    return ringAverage;

} /* End of SpeciesArray::RingAverage */

Vector_2D SpeciesArray::RingAverage( const Vector_1D ringArea, \
                                     const double totArea ) const
{

    Vector_2D ringAverage( nTime, Vector_1D( NVAR, 0.0E+00 ) );
    UInt iRing, iTime, N;
    double area;

    for ( iRing = 0; iRing < nRing; iRing++ ) {
        area = ringArea[iRing] / totArea;
        for ( iTime = 0; iTime < nTime; iTime++ ) {
            for ( N = 0; N < NVAR; N++ )
                ringAverage[iTime][N] += Species[N][iTime][iRing] * area;
        }
    }

    return ringAverage;

} /* End of SpeciesArray::RingAverage */

/* End of Species.cpp */
