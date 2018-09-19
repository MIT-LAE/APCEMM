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
#include <iomanip>
#include <vector>

#include "Parameters.hpp"
#include "Interface.hpp"
#include "Cluster.hpp"
#include "Ring.hpp"

typedef std::vector<double> Real_1DVector;

typedef double RealDouble;

class Mesh
{
    public:

        Mesh( );
        ~Mesh( );
        Mesh( const Mesh &m );
        Mesh& operator=( const Mesh &m );
        void Ring2Mesh( Cluster &c );
        Real_1DVector getX( ) const;
        Real_1DVector getY( ) const;
        RealDouble gethx( ) const;
        RealDouble gethy( ) const;
        unsigned int getNx( ) const;
        unsigned int getNy( ) const;
        std::vector<std::vector<std::vector<bool> > > getMap( ) const;
        std::vector<std::vector<std::pair<unsigned int, unsigned int> > > getList() const;
        std::vector<unsigned int> getnMap( ) const;
        void Debug() const;

    private:

        Real_1DVector x, y;
        RealDouble xlim, ylim;
        RealDouble hx, hy;
        unsigned int nx, ny;
        std::vector<unsigned int> nCellMap;
        std::vector<std::vector<std::pair<unsigned int, unsigned int> > > indList;
        std::vector<std::vector<std::vector<bool> > > RingMeshMap;

};

#endif /* MESH_H_INCLUDED */

