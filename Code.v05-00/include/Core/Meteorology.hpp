/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Meteorology Header File                                          */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 11/2/2018                                 */
/* File                 : Meteorology.hpp                           */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef METEOROLOGY_H_INCLUDED
#define METEOROLOGY_H_INCLUDED

#include <iostream>
#include "Core/Mesh.hpp"
#include "Util/MetFunction.hpp"

class Meteorology
{

    public:

        Meteorology( );
        Meteorology( const bool loadFile, \
                     const Mesh &m, \
                     const double temperature_K,  \
                     const double altitude, \
                     const double LapseRate, \
                     const bool DBG );
        Meteorology( const Meteorology &met );
        ~Meteorology( );

        double alt( unsigned int j ) const { return alt_[j]; }
        double press( unsigned int j ) const { return press_[j]; }

        double temp( unsigned int j, unsigned int i) const { return temp_[j][i]; }
        double H2O( unsigned int j, unsigned int i) const { return H2O_[j][i]; }

        std::vector<std::vector<double>> Temp() const { return temp_; }
        std::vector<double> Press() const { return press_; }

        friend class Solution;

    protected:

        /* Load data from file? */
        const bool LOAD;

        /* Assume that pressure only depends on the vertical coordinate */
        std::vector<double> alt_;
        std::vector<double> press_;
        /* Temperature and humidity fields can potentially be 2D fields */
        std::vector<std::vector<double>> temp_;
        std::vector<std::vector<double>> H2O_;


};


#endif /* METEOROLOGY_H_INCLUDED */
