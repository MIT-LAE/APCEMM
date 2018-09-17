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
    /* Default Constructor */

} /* End of SpeciesArray::SpeciesArray */

SpeciesArray::SpeciesArray( unsigned int nRing, unsigned int nTime )
{

    /* Constructor */
        
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

} /* End of SpeciesArray::SpeciesArray */

SpeciesArray::SpeciesArray( const SpeciesArray &sp )
{

    CO2 = sp.CO2;
    PPN = sp.PPN;
    BrNO2 = sp.BrNO2;
    IEPOX = sp.IEPOX;
    PMNN = sp.PMNN;
    N2O = sp.N2O;
    N = sp.N;
    PAN = sp.PAN;
    ALK4 = sp.ALK4;
    MAP = sp.MAP;
    MPN = sp.MPN;
    Cl2O2 = sp.Cl2O2;
    ETP = sp.ETP;
    HNO2 = sp.HNO2;
    C3H8 = sp.C3H8;
    RA3P = sp.RA3P;
    RB3P = sp.RB3P;
    OClO = sp.OClO;
    ClNO2 = sp.ClNO2;
    ISOP = sp.ISOP;
    HNO4 = sp.HNO4;
    MAOP = sp.MAOP;
    MP = sp.MP;
    ClOO = sp.ClOO;
    RP = sp.RP;
    BrCl = sp.BrCl;
    PP = sp.PP;
    PRPN = sp.PRPN;
    SO4 = sp.SO4;
    Br2 = sp.Br2;
    ETHLN = sp.ETHLN;
    MVKN = sp.MVKN;
    R4P = sp.R4P;
    C2H6 = sp.C2H6;
    RIP = sp.RIP;
    VRP = sp.VRP;
    ATOOH = sp.ATOOH;
    IAP = sp.IAP;
    DHMOB = sp.DHMOB;
    MOBA = sp.MOBA;
    MRP = sp.MRP;
    N2O5 = sp.N2O5;
    ISNOHOO = sp.ISNOHOO;
    ISNP = sp.ISNP;
    ISOPNB = sp.ISOPNB;
    IEPOXOO = sp.IEPOXOO;
    MACRNO2 = sp.MACRNO2;
    ROH = sp.ROH;
    MOBAOO = sp.MOBAOO;
    DIBOO = sp.DIBOO;
    PMN = sp.PMN;
    ISNOOB = sp.ISNOOB;
    INPN = sp.INPN;
    H = sp.H;
    BrNO3 = sp.BrNO3;
    PRPE = sp.PRPE;
    MVKOO = sp.MVKOO;
    Cl2 = sp.Cl2;
    ISOPND = sp.ISOPND;
    HOBr = sp.HOBr;
    A3O2 = sp.A3O2;
    PROPNN = sp.PROPNN;
    GLYX = sp.GLYX;
    MAOPO2 = sp.MAOPO2;
    CH4 = sp.CH4;
    GAOO = sp.GAOO;
    B3O2 = sp.B3O2;
    ACET = sp.ACET;
    MACRN = sp.MACRN;
    CH2OO = sp.CH2OO;
    MGLYOO = sp.MGLYOO;
    VRO2 = sp.VRO2;
    MGLOO = sp.MGLOO;
    MACROO = sp.MACROO;
    PO2 = sp.PO2;
    CH3CHOO = sp.CH3CHOO;
    MAN2 = sp.MAN2;
    ISNOOA = sp.ISNOOA;
    H2O2 = sp.H2O2;
    PRN1 = sp.PRN1;
    ETO2 = sp.ETO2;
    KO2 = sp.KO2;
    RCO3 = sp.RCO3;
    HC5OO = sp.HC5OO;
    GLYC = sp.GLYC;
    ClNO3 = sp.ClNO3;
    RIO2 = sp.RIO2;
    R4N1 = sp.R4N1;
    HOCl = sp.HOCl;
    ATO2 = sp.ATO2;
    HNO3 = sp.HNO3;
    ISN1 = sp.ISN1;
    MAO3 = sp.MAO3;
    MRO2 = sp.MRO2;
    INO2 = sp.INO2;
    HAC = sp.HAC;
    HC5 = sp.HC5;
    MGLY = sp.MGLY;
    ISOPNBO2 = sp.ISOPNBO2;
    ISOPNDO2 = sp.ISOPNDO2;
    R4O2 = sp.R4O2;
    R4N2 = sp.R4N2;
    BrO = sp.BrO;
    RCHO = sp.RCHO;
    MEK = sp.MEK;
    ClO = sp.ClO;
    MACR = sp.MACR;
    SO2 = sp.SO2;
    MVK = sp.MVK;
    ALD2 = sp.ALD2;
    MCO3 = sp.MCO3;
    CH2O = sp.CH2O;
    H2O = sp.H2O;
    Br = sp.Br;
    NO = sp.NO;
    NO3 = sp.NO3;
    Cl = sp.Cl;
    O = sp.O;
    O1D = sp.O1D;
    O3 = sp.O3;
    HO2 = sp.HO2;
    NO2 = sp.NO2;
    OH = sp.OH;
    HBr = sp.HBr;
    HCl = sp.HCl;
    CO = sp.CO;
    MO2 = sp.MO2;
    ACTA = sp.ACTA;
    EOH = sp.EOH;
    H2 = sp.H2;
    HCOOH = sp.HCOOH;
    MOH = sp.MOH;
    N2 = sp.N2;
    O2 = sp.O2;
    RCOOH = sp.RCOOH;  

} /* End of SpeciesArray::SpeciesArray */
    
SpeciesArray& SpeciesArray::operator=( const SpeciesArray &sp )
{

    if ( &sp == this )
        return *this;
    
    CO2 = sp.CO2;
    PPN = sp.PPN;
    BrNO2 = sp.BrNO2;
    IEPOX = sp.IEPOX;
    PMNN = sp.PMNN;
    N2O = sp.N2O;
    N = sp.N;
    PAN = sp.PAN;
    ALK4 = sp.ALK4;
    MAP = sp.MAP;
    MPN = sp.MPN;
    Cl2O2 = sp.Cl2O2;
    ETP = sp.ETP;
    HNO2 = sp.HNO2;
    C3H8 = sp.C3H8;
    RA3P = sp.RA3P;
    RB3P = sp.RB3P;
    OClO = sp.OClO;
    ClNO2 = sp.ClNO2;
    ISOP = sp.ISOP;
    HNO4 = sp.HNO4;
    MAOP = sp.MAOP;
    MP = sp.MP;
    ClOO = sp.ClOO;
    RP = sp.RP;
    BrCl = sp.BrCl;
    PP = sp.PP;
    PRPN = sp.PRPN;
    SO4 = sp.SO4;
    Br2 = sp.Br2;
    ETHLN = sp.ETHLN;
    MVKN = sp.MVKN;
    R4P = sp.R4P;
    C2H6 = sp.C2H6;
    RIP = sp.RIP;
    VRP = sp.VRP;
    ATOOH = sp.ATOOH;
    IAP = sp.IAP;
    DHMOB = sp.DHMOB;
    MOBA = sp.MOBA;
    MRP = sp.MRP;
    N2O5 = sp.N2O5;
    ISNOHOO = sp.ISNOHOO;
    ISNP = sp.ISNP;
    ISOPNB = sp.ISOPNB;
    IEPOXOO = sp.IEPOXOO;
    MACRNO2 = sp.MACRNO2;
    ROH = sp.ROH;
    MOBAOO = sp.MOBAOO;
    DIBOO = sp.DIBOO;
    PMN = sp.PMN;
    ISNOOB = sp.ISNOOB;
    INPN = sp.INPN;
    H = sp.H;
    BrNO3 = sp.BrNO3;
    PRPE = sp.PRPE;
    MVKOO = sp.MVKOO;
    Cl2 = sp.Cl2;
    ISOPND = sp.ISOPND;
    HOBr = sp.HOBr;
    A3O2 = sp.A3O2;
    PROPNN = sp.PROPNN;
    GLYX = sp.GLYX;
    MAOPO2 = sp.MAOPO2;
    CH4 = sp.CH4;
    GAOO = sp.GAOO;
    B3O2 = sp.B3O2;
    ACET = sp.ACET;
    MACRN = sp.MACRN;
    CH2OO = sp.CH2OO;
    MGLYOO = sp.MGLYOO;
    VRO2 = sp.VRO2;
    MGLOO = sp.MGLOO;
    MACROO = sp.MACROO;
    PO2 = sp.PO2;
    CH3CHOO = sp.CH3CHOO;
    MAN2 = sp.MAN2;
    ISNOOA = sp.ISNOOA;
    H2O2 = sp.H2O2;
    PRN1 = sp.PRN1;
    ETO2 = sp.ETO2;
    KO2 = sp.KO2;
    RCO3 = sp.RCO3;
    HC5OO = sp.HC5OO;
    GLYC = sp.GLYC;
    ClNO3 = sp.ClNO3;
    RIO2 = sp.RIO2;
    R4N1 = sp.R4N1;
    HOCl = sp.HOCl;
    ATO2 = sp.ATO2;
    HNO3 = sp.HNO3;
    ISN1 = sp.ISN1;
    MAO3 = sp.MAO3;
    MRO2 = sp.MRO2;
    INO2 = sp.INO2;
    HAC = sp.HAC;
    HC5 = sp.HC5;
    MGLY = sp.MGLY;
    ISOPNBO2 = sp.ISOPNBO2;
    ISOPNDO2 = sp.ISOPNDO2;
    R4O2 = sp.R4O2;
    R4N2 = sp.R4N2;
    BrO = sp.BrO;
    RCHO = sp.RCHO;
    MEK = sp.MEK;
    ClO = sp.ClO;
    MACR = sp.MACR;
    SO2 = sp.SO2;
    MVK = sp.MVK;
    ALD2 = sp.ALD2;
    MCO3 = sp.MCO3;
    CH2O = sp.CH2O;
    H2O = sp.H2O;
    Br = sp.Br;
    NO = sp.NO;
    NO3 = sp.NO3;
    Cl = sp.Cl;
    O = sp.O;
    O1D = sp.O1D;
    O3 = sp.O3;
    HO2 = sp.HO2;
    NO2 = sp.NO2;
    OH = sp.OH;
    HBr = sp.HBr;
    HCl = sp.HCl;
    CO = sp.CO;
    MO2 = sp.MO2;
    ACTA = sp.ACTA;
    EOH = sp.EOH;
    H2 = sp.H2;
    HCOOH = sp.HCOOH;
    MOH = sp.MOH;
    N2 = sp.N2;
    O2 = sp.O2;
    RCOOH = sp.RCOOH;

    return *this;

} /* End of SpeciesArray::operator= */

SpeciesArray::~SpeciesArray( )
{
    /* Destructor */

} /* End of SpeciesArray::~SpeciesArray */

unsigned int SpeciesArray::GetnRing() const
{

    return nRing;

} /* End of SpeciesArray::GetnRing */

unsigned int SpeciesArray::GetnTime() const
{

    return nTime;

} /* End of SpeciesArray::GetnTime */


/* End of Species.cpp */
