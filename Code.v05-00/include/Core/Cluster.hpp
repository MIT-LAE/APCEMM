/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Cluster Header File                                              */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 8/12/2018                                 */
/* File                 : Cluster.hpp                               */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef CLUSTER_H_INCLUDED
#define CLUSTER_H_INCLUDED

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include "Core/Ring.hpp"
#include "Util/ForwardDecl.hpp"

class Cluster
{

    public:

        Cluster( );
        Cluster( const UInt n, const bool sRing,                   \
                 const RealDouble sigma1, const RealDouble sigma2, \
                 const RealDouble d1, const RealDouble d2 );
        Cluster( const Cluster& cl );
        Cluster& operator=( const Cluster& cl );
        ~Cluster( );

        void ComputeRingAreas( const Vector_2D &cellAreas, const Vector_3D &weights );
        const UInt getnRing() const { return nR; }
        const bool halfRing() const { return semiRing; }
        const std::vector<Ring>& getRings() const { return Rings; }
        const std::vector<int>& getRingIndex() const { return ringIndices; }
        const Vector_1D& getRingArea( ) const { return ringAreas; }

        void PrintRings() const;
        void Debug() const;

    protected:

        UInt nR;
        bool semiRing;
        RealDouble sigmaX, sigmaY;
        RealDouble dH, dV;
        std::vector<Ring> Rings;
        std::vector<int> ringIndices;
        Vector_1D ringAreas;

};

#endif /* CLUSTER_H_INCLUDED */
