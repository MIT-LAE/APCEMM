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
/* Working directory    : /home/fritzt/APCEMM-SourceCode            */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef CLUSTER_H_INCLUDED
#define CLUSTER_H_INCLUDED

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include "Ring.hpp"

class Cluster
{

    public:

        Cluster( );
        Cluster( unsigned int n, bool sRing, double sigma1, double sigma2, double d1, double d2 );
        Cluster( const Cluster& cl );
        Cluster& operator=( const Cluster& cl );
        ~Cluster( );
        unsigned int nRing() const;
        bool halfRing() const;
        std::vector<Ring> GetRings() const; 
        void PrintRings() const;
        void Debug() const;

    protected:

        unsigned int nR;
        bool semiRing;
        double sigmaX, sigmaY;
        double dH, dV;
        std::vector<Ring> Rings;

};

#endif /* CLUSTER_H_INCLUDED */
