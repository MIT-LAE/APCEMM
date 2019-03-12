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

#include "Util/ForwardDecl.hpp"
#include "Core/Parameters.hpp"
#include "Core/Interface.hpp"
#include "Core/Cluster.hpp"
#include "Core/Ring.hpp"

class Mesh
{
    public:

        Mesh( );
        ~Mesh( );
        Mesh( const Mesh &m );
        Mesh& operator=( const Mesh &m );
        void Ring2Mesh( Cluster &c );
        const Vector_1D& x( ) const { return x_; }
        const Vector_1D& y( ) const { return y_; }
        const Vector_1D& xE( ) const { return x_e_; }
        const Vector_1D& yE( ) const { return y_e_; }
        const Vector_2D& areas( ) const { return areas_; }
        RealDouble totArea( ) const { return totArea_; }
        RealDouble hx( ) const { return hx_; }
        RealDouble hy( ) const { return hy_; }
        unsigned int Nx() const { return nx; }
        unsigned int Ny() const { return ny; }
        std::vector<std::vector<std::vector<bool> > > map( ) const { return RingMeshMap; }
        std::vector<std::vector<std::pair<unsigned int, unsigned int> > > list() const { return indList; }
        std::vector<unsigned int> nMap( ) const { return nCellMap; }
        void Debug() const;

    private:

        /* Cell center coordinates */
        Vector_1D x_, y_;

        /* Cell edges */
        Vector_1D x_e_, y_e_;

        /* Cell areas */
        Vector_2D areas_;

        /* Total area */
        RealDouble totArea_;

        RealDouble xlim, ylim_up, ylim_down;
        RealDouble hx_, hy_;
        unsigned int nx, ny;
        std::vector<unsigned int> nCellMap;
        std::vector<std::vector<std::pair<unsigned int, unsigned int> > > indList;
        std::vector<std::vector<std::vector<bool> > > RingMeshMap;

};

#endif /* MESH_H_INCLUDED */

