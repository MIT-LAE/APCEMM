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

#include "Util/ForwardDecl.hpp"
#include "Core/Cluster.hpp"
#include "Core/Input_Mod.hpp"

enum class MeshDomainLimitsSpec : unsigned char{
    CENTERED_LIMITS = 0, //all domain limits must be positive
    ABS_COORDS = 1 // regular x,y coords
};
class Mesh
{
    public:
        Mesh(const OptInput& optInput);
        Mesh(int nx, int ny, double xlim_left, double xlim_right, double ylim_up, double ylim_down, MeshDomainLimitsSpec limitsSpec = MeshDomainLimitsSpec::CENTERED_LIMITS);
        void Ring2Mesh( Cluster &c );
        void MapWeights( );
        inline const Vector_1D& x( ) const { return x_; }
        inline const Vector_1D& y( ) const { return y_; }
        inline const Vector_1D& xE( ) const { return x_e_; }
        inline const Vector_1D& yE( ) const { return y_e_; }
        inline const Vector_1D& dx( ) const { return dx_; }
        inline const Vector_1D& dy( ) const { return dy_; }
        inline const Vector_2D& areas( ) const { return areas_; }
        inline double hx( ) const { return hx_; }
        inline double hy( ) const { return hy_; }
        inline UInt Nx() const { return nx; }
        inline UInt Ny() const { return ny; }
        inline const Vector_3D& map( ) const { return weights; }
        inline const Vector_1Dui& nMap( ) const { return nCellMap; }
        inline const Vector_2Dui& mapIndex( ) const { return mapIndex_; }
        Vector_3D weights;

        void updateVertGrid( const Vector_1D& yE_new );

    private:

        void initCoordVectors(MeshDomainLimitsSpec limitsSpec);

        inline void calcAreas() {
            // UInt to be consistent with mesh size definition
            for ( UInt j = 0; j < ny; j++ ) {
                for ( UInt i = 0; i < nx; i++ ) {
                    areas_[j][i] = dx_[i] * dy_[j];
                }
            }
        }

        /* Cell center coordinates */
        Vector_1D x_, y_;

        /* Cell edges */
        Vector_1D x_e_, y_e_;

        /* Cell dimensions */
        Vector_1D dx_, dy_;

        /* Cell areas */
        Vector_2D areas_;

        /* Total area */
        double totArea_;
        /* Cell area */
        double cellArea_;

        double xlim_right, xlim_left, ylim_up, ylim_down;
        double hx_, hy_;
        UInt nx, ny;
        Vector_1Dui nCellMap;
        Vector_2Dui mapIndex_;

};

#endif /* MESH_H_INCLUDED */

