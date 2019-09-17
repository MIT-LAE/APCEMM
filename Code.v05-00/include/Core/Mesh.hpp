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
        void MapWeights( );
        const Vector_1D& x( ) const { return x_; }
        const Vector_1D& y( ) const { return y_; }
        const Vector_1D& xE( ) const { return x_e_; }
        const Vector_1D& yE( ) const { return y_e_; }
        const Vector_1D& dx( ) const { return dx_; }
        const Vector_1D& dy( ) const { return dy_; }
        const Vector_2D& areas( ) const { return areas_; }
        RealDouble totArea( ) const { return totArea_; }
        RealDouble hx( ) const { return hx_; }
        RealDouble hy( ) const { return hy_; }
        UInt Nx() const { return nx; }
        UInt Ny() const { return ny; }
        const Vector_3D& map( ) const { return weights; }
        const Vector_1Dui& nMap( ) const { return nCellMap; }
        const Vector_2Dui& mapIndex( ) const { return mapIndex_; }
        void Debug() const;

        Vector_3D weights;

    private:

        /* Cell center coordinates */
        Vector_1D x_, y_;

        /* Cell edges */
        Vector_1D x_e_, y_e_;

        /* Cell dimensions */
        Vector_1D dx_, dy_;

        /* Cell areas */
        Vector_2D areas_;

        /* Total area */
        RealDouble totArea_;
        /* Cell area */
        RealDouble cellArea_;

        RealDouble xlim, ylim_up, ylim_down;
        RealDouble hx_, hy_;
        UInt nx, ny;
        Vector_1Dui nCellMap;
        Vector_2Dui mapIndex_;

};

#endif /* MESH_H_INCLUDED */

