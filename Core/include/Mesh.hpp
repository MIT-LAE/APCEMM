/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Mesh Header File                                                 */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : Mesh.hpp                                  */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef MESH_H_INCLUDED
#define MESH_H_INCLUDED

#include <iostream>
#include <vector>

#include "Parameters.hpp"

typedef std::vector<double> Real_1DVector;
typedef double RealDouble;

class Mesh
{
    public:

        Mesh( );
        ~Mesh( );
        Mesh( const Mesh &m );
        Mesh& operator=( const Mesh &m );
        Real_1DVector getX( ) const;
        Real_1DVector getY( ) const;
        RealDouble gethx( ) const;
        RealDouble gethy( ) const;
        unsigned int getNx( ) const;
        unsigned int getNy( ) const;

    private:

        Real_1DVector x, y;
        RealDouble xlim, ylim;
        RealDouble hx, hy;
        unsigned int nx, ny;

};

#endif /* MESH_H_INCLUDED */

