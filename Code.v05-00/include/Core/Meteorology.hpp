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
#include "Core/Input_Mod.hpp"
#include "Util/PhysConstant.hpp"
#include "Util/MetFunction.hpp"

class Meteorology
{

    public:

        Meteorology( );
        Meteorology( const OptInput &USERINPUT,  \
                     const double solarTime_h,   \
                     const Mesh &m,              \
                     const double temperature_,  \
                     const double relHumidity_i, \
                     const double pressure_Pa,   \
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

        /* Met input type */
        unsigned int TYPE;

        /* Ambient parameters */
        const double TEMPERATURE;
        const double PRESSURE;
        const double RHI;
        double ALTITUDE;

        /* Temperature lapse rate */
        double LAPSERATE;

        /* Diurnal temperature variations */
        double DIURNAL_AMPL;
        double DIURNAL_PHASE;
        double diurnalPert;

        /* Non-uniformity length scales and temperature perturbations */
        double DELTAT;
        double TOP, BOT, LEFT, RIGHT;

        /* Assume that pressure only depends on the vertical coordinate */
        std::vector<double> alt_;
        std::vector<double> press_;
        /* Temperature and humidity fields can potentially be 2D fields */
        std::vector<std::vector<double>> temp_;
        std::vector<std::vector<double>> H2O_;


};


#endif /* METEOROLOGY_H_INCLUDED */
