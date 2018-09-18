/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Species Header File                                              */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : Species.hpp                               */
/* Working directory    : /home/fritzt/APCEMM-SourceCode            */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef SPECIES_H_INCLUDED
#define SPECIES_H_INCLUDED

#include <iostream> 
#include <vector>

#include "Interface.hpp"
#include "Structure.hpp"
#include "Mesh.hpp"

class SpeciesArray
{

    public:

        SpeciesArray( );
        SpeciesArray( unsigned int nRing, unsigned int nTime );
        SpeciesArray( const SpeciesArray &sp );
        SpeciesArray& operator=( const SpeciesArray &sp );
        SpeciesArray& operator+( const SpeciesArray &sp );
        ~SpeciesArray( );
        void FillIn( Solution &Data, Mesh &m, unsigned int nCounter );
        unsigned int GetnRing() const;
        unsigned int GetnTime() const;

        std::vector<std::vector<double> > CO2, PPN, BrNO2, IEPOX, PMNN, N2O, N, PAN, ALK4, MAP, MPN, Cl2O2, ETP, HNO2, C3H8, RA3P, RB3P, OClO, ClNO2, ISOP, HNO4, MAOP, MP, ClOO, RP, BrCl, PP, PRPN, SO4, Br2, ETHLN, MVKN, R4P, C2H6, RIP, VRP, ATOOH, IAP, DHMOB, MOBA, MRP, N2O5, ISNOHOO, ISNP, ISOPNB, IEPOXOO, MACRNO2, ROH, MOBAOO, DIBOO, PMN, ISNOOB, INPN, H, BrNO3, PRPE, MVKOO, Cl2, ISOPND, HOBr, A3O2, PROPNN, GLYX, MAOPO2, CH4, GAOO, B3O2, ACET, MACRN, CH2OO, MGLYOO, VRO2, MGLOO, MACROO, PO2, CH3CHOO, MAN2, ISNOOA, H2O2, PRN1, ETO2, KO2, RCO3, HC5OO, GLYC, ClNO3, RIO2, R4N1, HOCl, ATO2, HNO3, ISN1, MAO3, MRO2, INO2, HAC, HC5, MGLY, ISOPNBO2, ISOPNDO2, R4O2, R4N2, BrO, RCHO, MEK, ClO, MACR, SO2, MVK, ALD2, MCO3, CH2O, H2O, Br, NO, NO3, Cl, O, O1D, O3, HO2, NO2, OH, HBr, HCl, CO, MO2, ACTA, EOH, H2, HCOOH, MOH, N2, O2, RCOOH;

    protected:

        unsigned int nRing;
        unsigned int nTime;

    private:

};

#endif /* SPECIES_H_INCLUDED */


