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
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef EMISSION_H_INCLUDED
#define EMISSION_H_INCLUDED

#include <string>
#include "Core/Engine.hpp"
#include "Core/Fuel.hpp"

class Emission
{
    public:

        Emission( );
        Emission( const Engine &engine, const Fuel &fuel );
        Emission( const Emission &e );
        ~Emission( );
        void Populate_withEngine( const Engine &engine );
        void Populate_withFuel( const Fuel &fuel );
        Emission& operator=( const Emission &em );
        Emission& operator+( const Emission &em );
        double getCO2( ) const;
        double getH2O( ) const;
        double getNOx( ) const;
        double getNO( ) const;
        double getNO2( ) const;
        double getHNO2( ) const;
        double getSO2( ) const;
        double getCO( ) const;
        double getHC( ) const;
        double getCH4( ) const;
        double getC2H6( ) const;
        double getPRPE( ) const;
        double getALK4( ) const;
        double getCH2O( ) const;
        double getALD2( ) const;
        double getGLYX( ) const;
        double getMGLY( ) const;
        double getSoot( ) const;
        double getSootRad( ) const;
        std::string getEngineName( ) const;
        std::string getFuelChem( ) const;
        void Debug( ) const;

    protected:
        
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
