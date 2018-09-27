/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Aircraft Header File                                             */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : Aircraft.hpp                              */
/* Working directory    : /home/fritzt/APCEMM-SourceCode            */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef AIRCRAFT_H_INCLUDED
#define AIRCRAFT_H_INCLUDED

#include <string>
#include <iostream>
#include <iomanip>
#include <cstring>

#include "Engine.hpp"
#include "Vortex.hpp"

class Aircraft 
{
    public:

        Aircraft( );
        Aircraft( const char *aircraftName, double temperature_K, double pressure_Pa, double relHumidity_w );
        Aircraft( const Aircraft &ac );
        Aircraft& operator=( const Aircraft &ac ); 
        ~Aircraft( );
        void Debug( ) const;
        std::string getName() const;
        double getVFlight() const;
        double getMach() const;
        double getWingSpan() const;
        double getMTOW() const;
        double getCurrMass() const;
        Engine getEngine() const;
        double getFuelFlow() const;
        unsigned int getEngNumber() const;
        double getVortexdeltaz1() const;
        double getVortexdeltazw() const;

    protected:

        /* Aircraft name */
        std::string Name;

        /* Flight speed & mach Number */
        double vFlight_ms;
        double machNumber;

        /* Dimensions */
        double wingSpan; /* [m] */

        /* Weight */
        double MTOW; /* [kg] */
        double currMass; /* [kg] */

        /* Engine */
        Engine engine;

        /* Number of engines */
        unsigned int engNumber;

        /* Vortex */
        Vortex vortex;

    private:

};

#endif /* AIRCRAFT_H_INCLUDED */
