/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Fuel Header File                                                 */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : Fuel.hpp                                  */
/* Working directory    : /home/fritzt/APCEMM-SourceCode            */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef FUEL_H_INCLUDED
#define FUEL_H_INCLUDED

#include <string>
#include <iostream>
#include <cstring>

class Fuel
{
    public:

        Fuel( const char *fuelName );
        ~Fuel( );
        void GetAtoms( const char *fuelChem );

        /* Atomic composition: CxHy */
        double atomC;
        double atomH;
        double atomN;
        double atomS;

        /* Sulfur content */
        double FSC; /* [ppm] */

    private:

};

#endif /* FUEL_H_INCLUDED */
