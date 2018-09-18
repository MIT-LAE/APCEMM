/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Mesh Program File                                                */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : Mesh.cpp                                  */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "Mesh.hpp"

Mesh::Mesh( )
{

    /* Default Constructor */

    nx = NX;
    ny = NY;

    xlim = XLIM;
    ylim = YLIM;

    hx = 2 * xlim / nx;
    hy = 2 * ylim / ny;

    for ( unsigned int i = 0; i < nx; i++ ) {
        x.push_back( (double) 0.0 );
        x[i] = i * hx - xlim + hx / 2;
    }

    for ( unsigned int j = 0; j < ny; j++ ) {
        y.push_back( (double) 0.0 );
        y[j] = j * hy - ylim + hy / 2;
    }
    
} /* End of Mesh::Mesh */

Mesh::Mesh( const Mesh &m )
{

    x = m.x;
    y = m.y;
    xlim = m.xlim;
    ylim = m.ylim;
    hx = m.hx;
    hy = m.hy;
    nx = m.nx;
    ny = m.ny;

} /* End of Mesh::Mesh */

Mesh& Mesh::operator=( const Mesh &m )
{

    if ( &m == this )
        return *this;

    x = m.x;
    y = m.y;
    xlim = m.xlim;
    ylim = m.ylim;
    hx = m.hx;
    hy = m.hy;
    nx = m.nx;
    ny = m.ny;
    return *this;

} /* End of Mesh::operator= */

Mesh::~Mesh( )
{

    /* Destructor */

} /* End of Mesh::~Mesh */

void Mesh::Ring2Mesh( Cluster &c )
{
    
    unsigned int nRing = c.nRing();

    std::vector<std::vector<bool >> v2d;
    std::vector<bool> v1d;
    v1d = std::vector<bool>( NX, 0 );
    for ( unsigned int iRing = 0; iRing < nRing; iRing++ ) {
        nCellMap.push_back( 0 );
        RingMeshMap.push_back( v2d );
        for ( unsigned int iNy = 0; iNy < NY; iNy++ ) {
            RingMeshMap[iRing].push_back( v1d );
        }
    }

    std::vector<Ring> RingV;
    RingV = c.getRings();
    double hAxis, vAxis, hAxis_in, vAxis_in;
    double xRatio, xRatio_in;
    int val;
    for ( unsigned int iRing = 0; iRing < nRing; iRing++ ) {
        val = 0;
        if ( iRing == 0 ) {
            hAxis = RingV[iRing].getHAxis();
            vAxis = RingV[iRing].getVAxis();
            for ( unsigned int iNx = 0; iNx < NX; iNx++ ) {
                xRatio = ( x[iNx]  /  hAxis ) * ( x[iNx]  /  hAxis );
                for ( unsigned int jNy = 0; jNy < NY; jNy++ ) {
                    /* If point is within the region and out of the smaller region */
                    if ( xRatio + ( y[jNy] / vAxis ) * ( y[jNy] / vAxis ) <= 1 ) {
                        RingMeshMap[iRing][jNy][iNx] = 1;
                        val++;
                    }
                }
            }
        }
        else {
            hAxis = RingV[iRing].getHAxis();
            vAxis = RingV[iRing].getVAxis();
            hAxis_in = RingV[iRing-1].getHAxis();
            vAxis_in = RingV[iRing-1].getVAxis();
            for ( unsigned int iNx = 0; iNx < NX; iNx++ ) {
                xRatio    = ( x[iNx]  /  hAxis    ) * ( x[iNx]  /  hAxis    );
                xRatio_in = ( x[iNx]  /  hAxis_in ) * ( x[iNx]  /  hAxis_in );
                for ( unsigned int jNy = 0; jNy < NY; jNy++ ) {
                    /* If point is within the region and out of the smaller region */
                    if ( ( xRatio + ( y[jNy] / vAxis ) * ( y[jNy] / vAxis ) <= 1 ) && ( xRatio_in + ( y[jNy] / vAxis_in ) * ( y[jNy] / vAxis_in ) > 1 ) ) {
                        RingMeshMap[iRing][jNy][iNx] = 1;
                        val++;
                    }
                }
            }
        }
        if ( val != 0)
            nCellMap[iRing] = val;
        else
            std::cout << "Ring " << iRing << " has no cell in it (nMap = 0)!!" << std::endl;
    }



} /* End of Mesh::Ring2Mesh */

Real_1DVector Mesh::getX( ) const
{

    return x;

} /* End of Mesh::getX */

Real_1DVector Mesh::getY( ) const
{

    return y;

} /* End of Mesh::getY */
        
RealDouble Mesh::gethx( ) const
{

    return hx;

} /* End of Mesh::gethx */

RealDouble Mesh::gethy( ) const
{

    return hy;

} /* End of Mesh::gethy */

unsigned int Mesh::getNx( ) const
{

    return nx;

} /* End of Mesh::getNx */

unsigned int Mesh::getNy( ) const
{

    return ny;

} /* End of Mesh::getNy */

std::vector<std::vector<std::vector<bool> > > Mesh::getMap( ) const
{

    return RingMeshMap;

} /* End of Mesh::getMap */

std::vector<int> Mesh::getnMap( ) const
{

    return nCellMap;

} /* End of Mesh::getnMap */

/* End of Mesh.cpp */
