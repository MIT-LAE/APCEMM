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

    protected:

        /* Load data from file? */
        const bool LOAD;

        /* Assume that pressure only depends on the vertical coordinate */
        std::vector<double> alt;
        std::vector<double> press;
        /* Temperature and humidity fields can potentially be 2D fields */
        std::vector<std::vector<double>> temp;
        std::vector<std::vector<double>> H2O;


};


#endif /* METEOROLOGY_H_INCLUDED */
