/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Emission Header File                                             */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : Emission.hpp                              */
/* Working directory    : /home/fritzt/APCEMM-SourceCode            */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef EMISSION_H_INCLUDED
#define EMISSION_H_INCLUDED

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include "Engine.hpp"
#include "Fuel.hpp"

class Emission
{
    public:

        Emission( );
        void Populate( const Engine &engine, const Fuel &fuel );
        ~Emission( );
        void Populate_withEngine( const Engine &engine );
        void Populate_withFuel( const Fuel &fuel );
        Emission operator+( const Emission &emission_add );
        void Debug( ) const;

        /* Gaseous species */
        double CO2;  /* [g/kg fuel] */
        double H2O;  /* [g/kg fuel] */
        double NOx;  /* [g/kg fuel] */
        double NO;   /* [g/kg fuel] */
        double NO2;  /* [g/kg fuel] */
        double HNO2; /* [g/kg fuel] */
        double SO2;  /* [g/kg fuel] */
        double CO;   /* [g/kg fuel] */
        double HC;   /* [g/kg fuel] */
        double CH4;  /* [g/kg fuel] */
        double C2H6; /* [g/kg fuel] */
        double PRPE; /* [g/kg fuel] */
        double ALK4; /* [g/kg fuel] */
        double CH2O; /* [g/kg fuel] */
        double ALD2; /* [g/kg fuel] */
        double GLYX; /* [g/kg fuel] */
        double MGLY; /* [g/kg fuel] */

        /* BC */
        double Soot; /* [g/kg fuel] */
        double SootRad;  /* [m] */

        std::string engineName;
        std::string fuelChem;

    private:

};

#endif /* EMISSION_H_INCLUDED */
