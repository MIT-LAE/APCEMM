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

class Cluster
{

    public:

        Cluster( );
        Cluster( unsigned int n, bool sRing, double sigma1, double sigma2, double d1, double d2 );
        Cluster( const Cluster& cl );
        Cluster& operator=( const Cluster& cl );
        ~Cluster( );
        void ComputeRingAreas( std::vector<std::vector<double>> cellAreas, std::vector<std::vector<std::pair<unsigned int, unsigned int>>> map );
        unsigned int getnRing() const;
        bool halfRing() const;
        std::vector<Ring> getRings() const; 
        std::vector<int> getRingIndex() const;
        std::vector<double> getRingArea( ) const;
        void PrintRings() const;
        void Debug() const;

    protected:

        unsigned int nR;
        bool semiRing;
        double sigmaX, sigmaY;
        double dH, dV;
        std::vector<Ring> Rings;
        std::vector<int> ringIndices;
        std::vector<double> ringAreas;

};

#endif /* CLUSTER_H_INCLUDED */
