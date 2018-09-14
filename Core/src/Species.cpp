/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Species Program File                                             */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : Species.cpp                               */
/* Working directory    : /home/fritzt/APCEMM-SourceCode            */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "Species.hpp"

SpeciesArray::SpeciesArray( )
{
    /* Constructor */

} /* End of SpeciesArray::SpeciesArray */

void SpeciesArray::Build( unsigned int nRing, unsigned int nTime )
{
        
    for ( unsigned int i = 0; i < nTime; i++ ) {
        CO2.push_back( std::vector<double> ( nRing + 1 ) );
        PPN.push_back( std::vector<double> ( nRing + 1 ) );
        BrNO2.push_back( std::vector<double> ( nRing + 1 ) );
        IEPOX.push_back( std::vector<double> ( nRing + 1 ) );
        PMNN.push_back( std::vector<double> ( nRing + 1 ) );
        N2O.push_back( std::vector<double> ( nRing + 1 ) );
        N.push_back( std::vector<double> ( nRing + 1 ) );
        PAN.push_back( std::vector<double> ( nRing + 1 ) );
        ALK4.push_back( std::vector<double> ( nRing + 1 ) );
        MAP.push_back( std::vector<double> ( nRing + 1 ) );
        MPN.push_back( std::vector<double> ( nRing + 1 ) );
        Cl2O2.push_back( std::vector<double> ( nRing + 1 ) );
        ETP.push_back( std::vector<double> ( nRing + 1 ) );
        HNO2.push_back( std::vector<double> ( nRing + 1 ) );
        C3H8.push_back( std::vector<double> ( nRing + 1 ) );
        RA3P.push_back( std::vector<double> ( nRing + 1 ) );
        RB3P.push_back( std::vector<double> ( nRing + 1 ) );
        OClO.push_back( std::vector<double> ( nRing + 1 ) );
        ClNO2.push_back( std::vector<double> ( nRing + 1 ) );
        ISOP.push_back( std::vector<double> ( nRing + 1 ) );
        HNO4.push_back( std::vector<double> ( nRing + 1 ) );
        MAOP.push_back( std::vector<double> ( nRing + 1 ) );
        MP.push_back( std::vector<double> ( nRing + 1 ) );
        ClOO.push_back( std::vector<double> ( nRing + 1 ) );
        RP.push_back( std::vector<double> ( nRing + 1 ) );
        BrCl.push_back( std::vector<double> ( nRing + 1 ) );
        PP.push_back( std::vector<double> ( nRing + 1 ) );
        PRPN.push_back( std::vector<double> ( nRing + 1 ) );
        SO4.push_back( std::vector<double> ( nRing + 1 ) );
        Br2.push_back( std::vector<double> ( nRing + 1 ) );
        ETHLN.push_back( std::vector<double> ( nRing + 1 ) );
        MVKN.push_back( std::vector<double> ( nRing + 1 ) );
        R4P.push_back( std::vector<double> ( nRing + 1 ) );
        C2H6.push_back( std::vector<double> ( nRing + 1 ) );
        RIP.push_back( std::vector<double> ( nRing + 1 ) );
        VRP.push_back( std::vector<double> ( nRing + 1 ) );
        ATOOH.push_back( std::vector<double> ( nRing + 1 ) );
        IAP.push_back( std::vector<double> ( nRing + 1 ) );
        DHMOB.push_back( std::vector<double> ( nRing + 1 ) );
        MOBA.push_back( std::vector<double> ( nRing + 1 ) );
        MRP.push_back( std::vector<double> ( nRing + 1 ) );
        N2O5.push_back( std::vector<double> ( nRing + 1 ) );
        ISNOHOO.push_back( std::vector<double> ( nRing + 1 ) );
        ISNP.push_back( std::vector<double> ( nRing + 1 ) );
        ISOPNB.push_back( std::vector<double> ( nRing + 1 ) );
        IEPOXOO.push_back( std::vector<double> ( nRing + 1 ) );
        MACRNO2.push_back( std::vector<double> ( nRing + 1 ) );
        ROH.push_back( std::vector<double> ( nRing + 1 ) );
        MOBAOO.push_back( std::vector<double> ( nRing + 1 ) );
        DIBOO.push_back( std::vector<double> ( nRing + 1 ) );
        PMN.push_back( std::vector<double> ( nRing + 1 ) );
        ISNOOB.push_back( std::vector<double> ( nRing + 1 ) );
        INPN.push_back( std::vector<double> ( nRing + 1 ) );
        H.push_back( std::vector<double> ( nRing + 1 ) );
        BrNO3.push_back( std::vector<double> ( nRing + 1 ) );
        PRPE.push_back( std::vector<double> ( nRing + 1 ) );
        MVKOO.push_back( std::vector<double> ( nRing + 1 ) );
        Cl2.push_back( std::vector<double> ( nRing + 1 ) );
        ISOPND.push_back( std::vector<double> ( nRing + 1 ) );
        HOBr.push_back( std::vector<double> ( nRing + 1 ) );
        A3O2.push_back( std::vector<double> ( nRing + 1 ) );
        PROPNN.push_back( std::vector<double> ( nRing + 1 ) );
        GLYX.push_back( std::vector<double> ( nRing + 1 ) );
        MAOPO2.push_back( std::vector<double> ( nRing + 1 ) );
        CH4.push_back( std::vector<double> ( nRing + 1 ) );
        GAOO.push_back( std::vector<double> ( nRing + 1 ) );
        B3O2.push_back( std::vector<double> ( nRing + 1 ) );
        ACET.push_back( std::vector<double> ( nRing + 1 ) );
        MACRN.push_back( std::vector<double> ( nRing + 1 ) );
        CH2OO.push_back( std::vector<double> ( nRing + 1 ) );
        MGLYOO.push_back( std::vector<double> ( nRing + 1 ) );
        VRO2.push_back( std::vector<double> ( nRing + 1 ) );
        MGLOO.push_back( std::vector<double> ( nRing + 1 ) );
        MACROO.push_back( std::vector<double> ( nRing + 1 ) );
        PO2.push_back( std::vector<double> ( nRing + 1 ) );
        CH3CHOO.push_back( std::vector<double> ( nRing + 1 ) );
        MAN2.push_back( std::vector<double> ( nRing + 1 ) );
        ISNOOA.push_back( std::vector<double> ( nRing + 1 ) );
        H2O2.push_back( std::vector<double> ( nRing + 1 ) );
        PRN1.push_back( std::vector<double> ( nRing + 1 ) );
        ETO2.push_back( std::vector<double> ( nRing + 1 ) );
        KO2.push_back( std::vector<double> ( nRing + 1 ) );
        RCO3.push_back( std::vector<double> ( nRing + 1 ) );
        HC5OO.push_back( std::vector<double> ( nRing + 1 ) );
        GLYC.push_back( std::vector<double> ( nRing + 1 ) );
        ClNO3.push_back( std::vector<double> ( nRing + 1 ) );
        RIO2.push_back( std::vector<double> ( nRing + 1 ) );
        R4N1.push_back( std::vector<double> ( nRing + 1 ) );
        HOCl.push_back( std::vector<double> ( nRing + 1 ) );
        ATO2.push_back( std::vector<double> ( nRing + 1 ) );
        HNO3.push_back( std::vector<double> ( nRing + 1 ) );
        ISN1.push_back( std::vector<double> ( nRing + 1 ) );
        MAO3.push_back( std::vector<double> ( nRing + 1 ) );
        MRO2.push_back( std::vector<double> ( nRing + 1 ) );
        INO2.push_back( std::vector<double> ( nRing + 1 ) );
        HAC.push_back( std::vector<double> ( nRing + 1 ) );
        HC5.push_back( std::vector<double> ( nRing + 1 ) );
        MGLY.push_back( std::vector<double> ( nRing + 1 ) );
        ISOPNBO2.push_back( std::vector<double> ( nRing + 1 ) );
        ISOPNDO2.push_back( std::vector<double> ( nRing + 1 ) );
        R4O2.push_back( std::vector<double> ( nRing + 1 ) );
        R4N2.push_back( std::vector<double> ( nRing + 1 ) );
        BrO.push_back( std::vector<double> ( nRing + 1 ) );
        RCHO.push_back( std::vector<double> ( nRing + 1 ) );
        MEK.push_back( std::vector<double> ( nRing + 1 ) );
        ClO.push_back( std::vector<double> ( nRing + 1 ) );
        MACR.push_back( std::vector<double> ( nRing + 1 ) );
        SO2.push_back( std::vector<double> ( nRing + 1 ) );
        MVK.push_back( std::vector<double> ( nRing + 1 ) );
        ALD2.push_back( std::vector<double> ( nRing + 1 ) );
        MCO3.push_back( std::vector<double> ( nRing + 1 ) );
        CH2O.push_back( std::vector<double> ( nRing + 1 ) );
        H2O.push_back( std::vector<double> ( nRing + 1 ) );
        Br.push_back( std::vector<double> ( nRing + 1 ) );
        NO.push_back( std::vector<double> ( nRing + 1 ) );
        NO3.push_back( std::vector<double> ( nRing + 1 ) );
        Cl.push_back( std::vector<double> ( nRing + 1 ) );
        O.push_back( std::vector<double> ( nRing + 1 ) );
        O1D.push_back( std::vector<double> ( nRing + 1 ) );
        O3.push_back( std::vector<double> ( nRing + 1 ) );
        HO2.push_back( std::vector<double> ( nRing + 1 ) );
        NO2.push_back( std::vector<double> ( nRing + 1 ) );
        OH.push_back( std::vector<double> ( nRing + 1 ) );
        HBr.push_back( std::vector<double> ( nRing + 1 ) );
        HCl.push_back( std::vector<double> ( nRing + 1 ) );
        CO.push_back( std::vector<double> ( nRing + 1 ) );
        MO2.push_back( std::vector<double> ( nRing + 1 ) );
    }

} /* End of SpeciesArray::Build */

SpeciesArray::~SpeciesArray( )
{
    /* Destructor */

} /* End of SpeciesArray::~SpeciesArray */



