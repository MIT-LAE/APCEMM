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
#include "Util/PhysConstant.hpp"
#include "Util/MetFunction.hpp"

class Meteorology
{

    public:

        Meteorology( );
        Meteorology( const bool loadFile,        \
                     const double solarTime_h,   \
                     const Mesh &m,              \
                     const double temperature_,  \
                     const double altitude_,     \
                     const double LapseRate_,    \
                     const bool cstDepth_,       \
                     const double depth_,        \
                     const bool DBG = 0 );
        Meteorology( const Meteorology &met );
        ~Meteorology( );

        void Update( const double solarTime_h, const Mesh &m, \
                     const double dTrav_x, const double dTrav_y );

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

        /* Ambient parameters */
        const double TEMPERATURE;
        const double ALTITUDE;
        const double LAPSERATE;
        const bool CSTDEPTH;
        const double DEPTH;

        double DIURNAL_AMPL;
        double DIURNAL_PHASE;

        double diurnalPert;

        /* Assume that pressure only depends on the vertical coordinate */
        std::vector<double> alt_;
        std::vector<double> press_;
        /* Temperature and humidity fields can potentially be 2D fields */
        std::vector<std::vector<double>> temp_;
        std::vector<std::vector<double>> H2O_;


};


#endif /* METEOROLOGY_H_INCLUDED */
