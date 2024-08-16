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
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef FUEL_H_INCLUDED
#define FUEL_H_INCLUDED

#include <string>
#include <cstring>

class Fuel
{
    public:

        Fuel( );
        Fuel( const char *fuelName );
        Fuel( const Fuel &f );
        Fuel& operator=( const Fuel &f );
        ~Fuel( );
        void getAtoms( const char *fuelChem );
        double getAtomC() const;
        double getAtomH() const;
        double getAtomN() const;
        double getAtomS() const;
        double getFSC() const;
        void setFSC( const double FSC_ );
        std::string getChemFormula() const;
        
    protected:

        /* Atomic composition: CxHy */
        double atomC;
        double atomH;
        double atomN;
        double atomS;

        /* Sulfur content */
        double FSC; /* [ppm] */

        std::string ChemFormula;

    private:

};

#endif /* FUEL_H_INCLUDED */
