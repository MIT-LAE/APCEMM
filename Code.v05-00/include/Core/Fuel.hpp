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
#include "Util/MolarWeights.hpp"

class Fuel
{
    public:

        Fuel( );
        /* EI(SO2) [g SO2/kg fuel] per ppm of fuel sulfur (mass), assuming all
           sulfur leaves as SO2: 1 ppm = 1e-3 g S/kg fuel, times M_SO2 / M_S */
        static constexpr double EI_SO2_PER_FSC = 1.0E-03 * MW_SO2 / MW_S;

        Fuel( const char *fuelName, const double FSC_ppm );
        /* Fuel whose sulfur content gives the emission index EI_SO2 [g/kg fuel] */
        static Fuel fromEISO2( const char *fuelName, const double EI_SO2 );
        Fuel( const Fuel &f );
        Fuel& operator=( const Fuel &f );
        ~Fuel( );
        void getAtoms( const char *fuelChem );
        double getAtomC() const;
        double getAtomH() const;
        double getAtomN() const;
        double getAtomS() const;
        double getFSC() const;
        double EI_SO2() const; /* [g SO2/kg fuel] */
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
