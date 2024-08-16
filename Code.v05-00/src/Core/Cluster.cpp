/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Cluster Program File                                             */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 8/12/2018                                 */
/* File                 : Cluster.cpp                               */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <iostream>
#include <iomanip>
#include "Core/Cluster.hpp"

Cluster::Cluster( )
{

    /* Default Constructor */
    nR = 0;
    semiRing = 0;
    sigmaX = 0;
    sigmaY = 0;
    dH = 0;
    dV = 0;
    for ( UInt iR = 0; iR < nR; iR++ ) {
        Rings.push_back( Ring() );
        ringIndices.push_back( iR );
        ringAreas.push_back( 0.0E+00 );
    }

} /* End of Cluster::Cluster */

Cluster::Cluster( const UInt n, const bool sRing,                   \
                  const double sigma1, const double sigma2, \
                  const double d1, const double d2 )
{

    /* Constructor */
    nR = n;
    semiRing = sRing;
    sigmaX = sigma1;
    sigmaY = sigma2;
    dH = d1;
    dV = d2;

    if ( sigmaX <= 0 )
        sigmaX = 73.796;
    if ( sigmaY <= 0 )
        sigmaY = 27.361;
    
    if ( dH <= 0.0 )
        dH = 15.0;
    if ( dV <= 0.0 )
        dV = 0.15;

    Vector_1D ringSizes;
    if ( !semiRing )
        ringSizes = {0.0, 0.25, 1.0, 2.0, 4.0, 6.0, 8.0, 12.0, 16.0, 20.0, 25.0, 30.0, 35.0, 40.0, 48.0};
    else {
        nR *= 2;
        ringSizes = {0.0, 0.0, 0.25, 0.25, 1.0, 1.0, 2.0, 2.0, 4.0, 4.0, 6.0, 6.0, 8.0, 8.0, 12.0, 12.0, \
            16.0, 16.0, 20.0, 20.0, 25.0, 25.0, 30.0, 30.0, 35.0, 35.0, 40.0, 40.0, 48.0, 48.0};
    }

    /* Create rings */
    double sigmaXRing, sigmaYRing;
    for ( UInt iRing = 0; iRing < nR; iRing++ ) {
        
        /* Add new element */
        Rings.push_back( Ring() );
        ringIndices.push_back( iRing );

        /* X and Y axis */
        sigmaXRing = sqrt( sigmaX * sigmaX + 20*16*dH*3600.0*ringSizes[iRing] );
        sigmaYRing = sqrt( sigmaY * sigmaY + 16*dV*3600.0*ringSizes[iRing] );

        /* Assign X and Y axis to current ring */
        Rings[iRing].Assign( sigmaXRing, sigmaYRing );

        if ( semiRing ) {
            /* Build half-rings */
            
            /* Add new element */
            Rings.push_back( Ring() );
            ringIndices.push_back( iRing );
            
            /* Copy previous ring */
            Rings[iRing+1] = Rings[iRing];

            /* Iterate */
            iRing++;
        }

    }


} /* End of Cluster::Cluster */

Cluster::Cluster( const Cluster& cl )
{

    /* Constructor */
    nR          = cl.nR;
    semiRing    = cl.semiRing;
    sigmaX      = cl.sigmaX;
    sigmaY      = cl.sigmaY;
    dH          = cl.dH;
    dV          = cl.dV;
    Rings       = cl.Rings;
    ringIndices = cl.ringIndices;
    ringAreas   = cl.ringAreas;

} /* End of Cluster::Cluster */

Cluster& Cluster::operator=( const Cluster &cl )
{

    if ( &cl == this )
        return *this;

    nR       = cl.nR;
    semiRing = cl.semiRing;
    sigmaX   = cl.sigmaX;
    sigmaY   = cl.sigmaY;
    dH       = cl.dH;
    dV       = cl.dV;
    Rings    = cl.Rings;
    ringIndices = cl.ringIndices;
    ringAreas= cl.ringAreas;

    return *this;
    

} /* End of Cluster::operator= */

Cluster::~Cluster( )
{

    /* Default destructor */

} /* End of Cluster::~Cluster */
        
void Cluster::ComputeRingAreas( const Vector_2D &cellAreas, const Vector_3D &weights ) 
{

    UInt iNx, jNy, iRing;

    UInt size = 0;

#ifdef RINGS
    size = nR;
#else
    size = weights.size();
#endif /* RINGS */

    for ( iRing = 0; iRing < size; iRing++ ) {
        ringAreas.push_back( 0.0E+00 );
        for ( jNy = 0; jNy < cellAreas.size(); jNy++ ) {
            for ( iNx = 0; iNx < cellAreas[0].size(); iNx++ ) {
                if ( weights[iRing][jNy][iNx] != 0.0E+00 )
                    ringAreas[iRing] += cellAreas[jNy][iNx];
            }
        }
    }

} /* End of Cluster::ComputeRingAreas */

void Cluster::PrintRings() const
{

    for ( UInt i = 0; i < nR; i++ )
        std::cout << "Ring's horizontal and vertical axis: " << Rings[i].getHAxis() << ", " << Rings[i].getVAxis() << " [m]" << std::endl;

} /* End of Cluster::PrintRings */

void Cluster::Debug() const
{

    std::streamsize ss = std::cout.precision();

    std::cout.precision(5);

    std::cout << std::endl;
    std::cout << "**** Input Debugger ****" << std::endl;
    std::cout << "Cluster parameters: " << std::endl;
    std::cout << std::endl;
    std::setw(15);
    std::cout << "Ring number   ";
    std::setw(15);
    std::cout << "Horizontal dimension   ";
    std::setw(15);
    std::cout << "Vertical dimension";
    std::cout << std::endl;
    for ( UInt i = 0; i < nR; i++ ) {
        std::setw(7);
        std::cout << "     ";
        std::setw(15);
        std::cout << i;
        std::setw(7);
        if ( i >= 10 )
            std::cout << "     ";
        else
            std::cout << "      ";
        std::setw(15);
        std::cout << "        ";
        std::setw(15);
        std::cout << Rings[i].getHAxis();
        std::cout << "        ";
        std::setw(15);
        std::cout << "        ";
        std::setw(15);
        std::cout << Rings[i].getVAxis();
        std::cout << "        ";
        std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << std::endl;
    
    std::cout.precision(ss);

} /* End of Cluster::Debug */

/* End of Cluster.cpp */
