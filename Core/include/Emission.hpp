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
#include <string>
#include <fstream>
#include <sstream>
#include "Engine.hpp"
#include "Fuel.hpp"

class Emission
{
    public:
        
        Emission( Engine engine, Fuel fuel );
        ~Emission( );
        void Populate_withEngine( Engine engine );
        void Populate_withFuel( Fuel fuel );

       /* Gaseous species */
       double EI_CO2; /* [g/kg fuel] */
       double EI_H2O; /* [g/kg fuel] */
       double EI_NOx; /* [g/kg fuel] */
       double EI_SO2; /* [g/kg fuel] */
       double EI_CO;  /* [g/kg fuel] */
       double EI_HC;  /* [g/kg fuel] */
       double EI_OH;  /* [g/kg fuel] */

       /* BC */
       double EI_Soot; /* [g/kg fuel] */
       double SootRad;  /* [m] */

    private:

};

#endif /* EMISSION_H_INCLUDED */
