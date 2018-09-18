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

/* End of Mesh.cpp */
