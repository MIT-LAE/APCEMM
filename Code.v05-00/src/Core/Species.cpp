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
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "Core/Species.hpp"

static const RealDouble ZERO = 1.00E-50;

SpeciesArray::SpeciesArray( )
{
    /* Default Constructor */

} /* End of SpeciesArray::SpeciesArray */

SpeciesArray::SpeciesArray( const UInt nRing_, const UInt nTime_, const bool halfRing_ )
{

    /* Constructor */

    nRing = nRing_;
    nTime = nTime_;
    halfRing = halfRing_;

    Vector_1D v1d( nRing, 0.0E+00 );

    for ( UInt i = 0; i < nTime; i++ ) {
        CO2.push_back( v1d );
        PPN.push_back( v1d );
        BrNO2.push_back( v1d );
        IEPOX.push_back( v1d );
        PMNN.push_back( v1d );
        N2O.push_back( v1d );
        N.push_back( v1d );
        PAN.push_back( v1d );
        ALK4.push_back( v1d );
        MAP.push_back( v1d );
        MPN.push_back( v1d );
        Cl2O2.push_back( v1d );
        ETP.push_back( v1d );
        HNO2.push_back( v1d );
        C3H8.push_back( v1d );
        RA3P.push_back( v1d );
        RB3P.push_back( v1d );
        OClO.push_back( v1d );
        ClNO2.push_back( v1d );
        ISOP.push_back( v1d );
        HNO4.push_back( v1d );
        MAOP.push_back( v1d );
        MP.push_back( v1d );
        ClOO.push_back( v1d );
        RP.push_back( v1d );
        BrCl.push_back( v1d );
        PP.push_back( v1d );
        PRPN.push_back( v1d );
        SO4.push_back( v1d );
        Br2.push_back( v1d );
        ETHLN.push_back( v1d );
        MVKN.push_back( v1d );
        R4P.push_back( v1d );
        C2H6.push_back( v1d );
        RIP.push_back( v1d );
        VRP.push_back( v1d );
        ATOOH.push_back( v1d );
        IAP.push_back( v1d );
        DHMOB.push_back( v1d );
        MOBA.push_back( v1d );
        MRP.push_back( v1d );
        N2O5.push_back( v1d );
        ISNOHOO.push_back( v1d );
        ISNP.push_back( v1d );
        ISOPNB.push_back( v1d );
        IEPOXOO.push_back( v1d );
        MACRNO2.push_back( v1d );
        ROH.push_back( v1d );
        MOBAOO.push_back( v1d );
        DIBOO.push_back( v1d );
        PMN.push_back( v1d );
        ISNOOB.push_back( v1d );
        INPN.push_back( v1d );
        H.push_back( v1d );
        BrNO3.push_back( v1d );
        PRPE.push_back( v1d );
        MVKOO.push_back( v1d );
        Cl2.push_back( v1d );
        ISOPND.push_back( v1d );
        HOBr.push_back( v1d );
        A3O2.push_back( v1d );
        PROPNN.push_back( v1d );
        GLYX.push_back( v1d );
        MAOPO2.push_back( v1d );
        CH4.push_back( v1d );
        GAOO.push_back( v1d );
        B3O2.push_back( v1d );
        ACET.push_back( v1d );
        MACRN.push_back( v1d );
        CH2OO.push_back( v1d );
        MGLYOO.push_back( v1d );
        VRO2.push_back( v1d );
        MGLOO.push_back( v1d );
        MACROO.push_back( v1d );
        PO2.push_back( v1d );
        CH3CHOO.push_back( v1d );
        MAN2.push_back( v1d );
        ISNOOA.push_back( v1d );
        H2O2.push_back( v1d );
        PRN1.push_back( v1d );
        ETO2.push_back( v1d );
        KO2.push_back( v1d );
        RCO3.push_back( v1d );
        HC5OO.push_back( v1d );
        GLYC.push_back( v1d );
        ClNO3.push_back( v1d );
        RIO2.push_back( v1d );
        R4N1.push_back( v1d );
        HOCl.push_back( v1d );
        ATO2.push_back( v1d );
        HNO3.push_back( v1d );
        ISN1.push_back( v1d );
        MAO3.push_back( v1d );
        MRO2.push_back( v1d );
        INO2.push_back( v1d );
        HAC.push_back( v1d );
        HC5.push_back( v1d );
        MGLY.push_back( v1d );
        ISOPNBO2.push_back( v1d );
        ISOPNDO2.push_back( v1d );
        R4O2.push_back( v1d );
        R4N2.push_back( v1d );
        BrO.push_back( v1d );
        RCHO.push_back( v1d );
        MEK.push_back( v1d );
        ClO.push_back( v1d );
        MACR.push_back( v1d );
        SO2.push_back( v1d );
        MVK.push_back( v1d );
        ALD2.push_back( v1d );
        MCO3.push_back( v1d );
        CH2O.push_back( v1d );
        H2O.push_back( v1d );
        Br.push_back( v1d );
        NO.push_back( v1d );
        NO3.push_back( v1d );
        Cl.push_back( v1d );
        O.push_back( v1d );
        O1D.push_back( v1d );
        O3.push_back( v1d );
        HO2.push_back( v1d );
        NO2.push_back( v1d );
        OH.push_back( v1d );
        HBr.push_back( v1d );
        HCl.push_back( v1d );
        CO.push_back( v1d );
        MO2.push_back( v1d );

        SO4L.push_back( v1d );
        H2OL.push_back( v1d );
        HNO3L.push_back( v1d );
        HClL.push_back( v1d );
        HOClL.push_back( v1d );
        HBrL.push_back( v1d );
        HOBrL.push_back( v1d );
        H2OS.push_back( v1d );
        HNO3S.push_back( v1d );

        sootDens.push_back( v1d );
        sootRadi.push_back( v1d );
        sootArea.push_back( v1d );
        iceDens.push_back( v1d );
        iceRadi.push_back( v1d );
        iceArea.push_back( v1d );
        sulfDens.push_back( v1d );
        sulfRadi.push_back( v1d );
        sulfArea.push_back( v1d );

    }

} /* End of SpeciesArray::SpeciesArray */


SpeciesArray::SpeciesArray( const SpeciesArray &sp )
{
    
    nRing    = sp.getnRing();
    nTime    = sp.getnTime();
    halfRing = sp.gethalfRing();

    CO2      = sp.CO2;
    PPN      = sp.PPN;
    BrNO2    = sp.BrNO2;
    IEPOX    = sp.IEPOX;
    PMNN     = sp.PMNN;
    N2O      = sp.N2O;
    N        = sp.N;
    PAN      = sp.PAN;
    ALK4     = sp.ALK4;
    MAP      = sp.MAP;
    MPN      = sp.MPN;
    Cl2O2    = sp.Cl2O2;
    ETP      = sp.ETP;
    HNO2     = sp.HNO2;
    C3H8     = sp.C3H8;
    RA3P     = sp.RA3P;
    RB3P     = sp.RB3P;
    OClO     = sp.OClO;
    ClNO2    = sp.ClNO2;
    ISOP     = sp.ISOP;
    HNO4     = sp.HNO4;
    MAOP     = sp.MAOP;
    MP       = sp.MP;
    ClOO     = sp.ClOO;
    RP       = sp.RP;
    BrCl     = sp.BrCl;
    PP       = sp.PP;
    PRPN     = sp.PRPN;
    SO4      = sp.SO4;
    Br2      = sp.Br2;
    ETHLN    = sp.ETHLN;
    MVKN     = sp.MVKN;
    R4P      = sp.R4P;
    C2H6     = sp.C2H6;
    RIP      = sp.RIP;
    VRP      = sp.VRP;
    ATOOH    = sp.ATOOH;
    IAP      = sp.IAP;
    DHMOB    = sp.DHMOB;
    MOBA     = sp.MOBA;
    MRP      = sp.MRP;
    N2O5     = sp.N2O5;
    ISNOHOO  = sp.ISNOHOO;
    ISNP     = sp.ISNP;
    ISOPNB   = sp.ISOPNB;
    IEPOXOO  = sp.IEPOXOO;
    MACRNO2  = sp.MACRNO2;
    ROH      = sp.ROH;
    MOBAOO   = sp.MOBAOO;
    DIBOO    = sp.DIBOO;
    PMN      = sp.PMN;
    ISNOOB   = sp.ISNOOB;
    INPN     = sp.INPN;
    H        = sp.H;
    BrNO3    = sp.BrNO3;
    PRPE     = sp.PRPE;
    MVKOO    = sp.MVKOO;
    Cl2      = sp.Cl2;
    ISOPND   = sp.ISOPND;
    HOBr     = sp.HOBr;
    A3O2     = sp.A3O2;
    PROPNN   = sp.PROPNN;
    GLYX     = sp.GLYX;
    MAOPO2   = sp.MAOPO2;
    CH4      = sp.CH4;
    GAOO     = sp.GAOO;
    B3O2     = sp.B3O2;
    ACET     = sp.ACET;
    MACRN    = sp.MACRN;
    CH2OO    = sp.CH2OO;
    MGLYOO   = sp.MGLYOO;
    VRO2     = sp.VRO2;
    MGLOO    = sp.MGLOO;
    MACROO   = sp.MACROO;
    PO2      = sp.PO2;
    CH3CHOO  = sp.CH3CHOO;
    MAN2     = sp.MAN2;
    ISNOOA   = sp.ISNOOA;
    H2O2     = sp.H2O2;
    PRN1     = sp.PRN1;
    ETO2     = sp.ETO2;
    KO2      = sp.KO2;
    RCO3     = sp.RCO3;
    HC5OO    = sp.HC5OO;
    GLYC     = sp.GLYC;
    ClNO3    = sp.ClNO3;
    RIO2     = sp.RIO2;
    R4N1     = sp.R4N1;
    HOCl     = sp.HOCl;
    ATO2     = sp.ATO2;
    HNO3     = sp.HNO3;
    ISN1     = sp.ISN1;
    MAO3     = sp.MAO3;
    MRO2     = sp.MRO2;
    INO2     = sp.INO2;
    HAC      = sp.HAC;
    HC5      = sp.HC5;
    MGLY     = sp.MGLY;
    ISOPNBO2 = sp.ISOPNBO2;
    ISOPNDO2 = sp.ISOPNDO2;
    R4O2     = sp.R4O2;
    R4N2     = sp.R4N2;
    BrO      = sp.BrO;
    RCHO     = sp.RCHO;
    MEK      = sp.MEK;
    ClO      = sp.ClO;
    MACR     = sp.MACR;
    SO2      = sp.SO2;
    MVK      = sp.MVK;
    ALD2     = sp.ALD2;
    MCO3     = sp.MCO3;
    CH2O     = sp.CH2O;
    H2O      = sp.H2O;
    Br       = sp.Br;
    NO       = sp.NO;
    NO3      = sp.NO3;
    Cl       = sp.Cl;
    O        = sp.O;
    O1D      = sp.O1D;
    O3       = sp.O3;
    HO2      = sp.HO2;
    NO2      = sp.NO2;
    OH       = sp.OH;
    HBr      = sp.HBr;
    HCl      = sp.HCl;
    CO       = sp.CO;
    MO2      = sp.MO2;
    ACTA     = sp.ACTA;
    EOH      = sp.EOH;
    H2       = sp.H2;
    HCOOH    = sp.HCOOH;
    MOH      = sp.MOH;
    N2       = sp.N2;
    O2       = sp.O2;
    RCOOH    = sp.RCOOH;

    SO4L  = sp.SO4L;
    H2OL  = sp.H2OL;
    HNO3L = sp.HNO3L;
    HClL  = sp.HClL;
    HOClL = sp.HOClL;
    HBrL  = sp.HBrL;
    HOBrL = sp.HOBrL;

    sootDens = sp.sootDens;
    sootRadi = sp.sootRadi;
    sootArea = sp.sootArea;
    iceDens  = sp.iceDens;
    iceRadi  = sp.iceRadi;
    iceArea  = sp.iceArea;
    sulfDens = sp.sulfDens;
    sulfRadi = sp.sulfRadi;
    sulfArea = sp.sulfArea;


} /* End of SpeciesArray::SpeciesArray */
    
SpeciesArray& SpeciesArray::operator=( const SpeciesArray &sp )
{

    if ( &sp == this )
        return *this;

    nRing    = sp.getnRing();
    nTime    = sp.getnTime();
    halfRing = sp.gethalfRing();

    CO2      = sp.CO2;
    PPN      = sp.PPN;
    BrNO2    = sp.BrNO2;
    IEPOX    = sp.IEPOX;
    PMNN     = sp.PMNN;
    N2O      = sp.N2O;
    N        = sp.N;
    PAN      = sp.PAN;
    ALK4     = sp.ALK4;
    MAP      = sp.MAP;
    MPN      = sp.MPN;
    Cl2O2    = sp.Cl2O2;
    ETP      = sp.ETP;
    HNO2     = sp.HNO2;
    C3H8     = sp.C3H8;
    RA3P     = sp.RA3P;
    RB3P     = sp.RB3P;
    OClO     = sp.OClO;
    ClNO2    = sp.ClNO2;
    ISOP     = sp.ISOP;
    HNO4     = sp.HNO4;
    MAOP     = sp.MAOP;
    MP       = sp.MP;
    ClOO     = sp.ClOO;
    RP       = sp.RP;
    BrCl     = sp.BrCl;
    PP       = sp.PP;
    PRPN     = sp.PRPN;
    SO4      = sp.SO4;
    Br2      = sp.Br2;
    ETHLN    = sp.ETHLN;
    MVKN     = sp.MVKN;
    R4P      = sp.R4P;
    C2H6     = sp.C2H6;
    RIP      = sp.RIP;
    VRP      = sp.VRP;
    ATOOH    = sp.ATOOH;
    IAP      = sp.IAP;
    DHMOB    = sp.DHMOB;
    MOBA     = sp.MOBA;
    MRP      = sp.MRP;
    N2O5     = sp.N2O5;
    ISNOHOO  = sp.ISNOHOO;
    ISNP     = sp.ISNP;
    ISOPNB   = sp.ISOPNB;
    IEPOXOO  = sp.IEPOXOO;
    MACRNO2  = sp.MACRNO2;
    ROH      = sp.ROH;
    MOBAOO   = sp.MOBAOO;
    DIBOO    = sp.DIBOO;
    PMN      = sp.PMN;
    ISNOOB   = sp.ISNOOB;
    INPN     = sp.INPN;
    H        = sp.H;
    BrNO3    = sp.BrNO3;
    PRPE     = sp.PRPE;
    MVKOO    = sp.MVKOO;
    Cl2      = sp.Cl2;
    ISOPND   = sp.ISOPND;
    HOBr     = sp.HOBr;
    A3O2     = sp.A3O2;
    PROPNN   = sp.PROPNN;
    GLYX     = sp.GLYX;
    MAOPO2   = sp.MAOPO2;
    CH4      = sp.CH4;
    GAOO     = sp.GAOO;
    B3O2     = sp.B3O2;
    ACET     = sp.ACET;
    MACRN    = sp.MACRN;
    CH2OO    = sp.CH2OO;
    MGLYOO   = sp.MGLYOO;
    VRO2     = sp.VRO2;
    MGLOO    = sp.MGLOO;
    MACROO   = sp.MACROO;
    PO2      = sp.PO2;
    CH3CHOO  = sp.CH3CHOO;
    MAN2     = sp.MAN2;
    ISNOOA   = sp.ISNOOA;
    H2O2     = sp.H2O2;
    PRN1     = sp.PRN1;
    ETO2     = sp.ETO2;
    KO2      = sp.KO2;
    RCO3     = sp.RCO3;
    HC5OO    = sp.HC5OO;
    GLYC     = sp.GLYC;
    ClNO3    = sp.ClNO3;
    RIO2     = sp.RIO2;
    R4N1     = sp.R4N1;
    HOCl     = sp.HOCl;
    ATO2     = sp.ATO2;
    HNO3     = sp.HNO3;
    ISN1     = sp.ISN1;
    MAO3     = sp.MAO3;
    MRO2     = sp.MRO2;
    INO2     = sp.INO2;
    HAC      = sp.HAC;
    HC5      = sp.HC5;
    MGLY     = sp.MGLY;
    ISOPNBO2 = sp.ISOPNBO2;
    ISOPNDO2 = sp.ISOPNDO2;
    R4O2     = sp.R4O2;
    R4N2     = sp.R4N2;
    BrO      = sp.BrO;
    RCHO     = sp.RCHO;
    MEK      = sp.MEK;
    ClO      = sp.ClO;
    MACR     = sp.MACR;
    SO2      = sp.SO2;
    MVK      = sp.MVK;
    ALD2     = sp.ALD2;
    MCO3     = sp.MCO3;
    CH2O     = sp.CH2O;
    H2O      = sp.H2O;
    Br       = sp.Br;
    NO       = sp.NO;
    NO3      = sp.NO3;
    Cl       = sp.Cl;
    O        = sp.O;
    O1D      = sp.O1D;
    O3       = sp.O3;
    HO2      = sp.HO2;
    NO2      = sp.NO2;
    OH       = sp.OH;
    HBr      = sp.HBr;
    HCl      = sp.HCl;
    CO       = sp.CO;
    MO2      = sp.MO2;
    ACTA     = sp.ACTA;
    EOH      = sp.EOH;
    H2       = sp.H2;
    HCOOH    = sp.HCOOH;
    MOH      = sp.MOH;
    N2       = sp.N2;
    O2       = sp.O2;
    RCOOH    = sp.RCOOH;
    
    SO4L  = sp.SO4L;
    H2OL  = sp.H2OL;
    HNO3L = sp.HNO3L;
    HClL  = sp.HClL;
    HOClL = sp.HOClL;
    HBrL  = sp.HBrL;
    HOBrL = sp.HOBrL;

    sootDens = sp.sootDens;
    sootRadi = sp.sootRadi;
    sootArea = sp.sootArea;
    iceDens  = sp.iceDens;
    iceRadi  = sp.iceRadi;
    iceArea  = sp.iceArea;
    sulfDens = sp.sulfDens;
    sulfRadi = sp.sulfRadi;
    sulfArea = sp.sulfArea;
    return *this;

} /* End of SpeciesArray::operator= */

SpeciesArray& SpeciesArray::operator+( const SpeciesArray &sp )
{

    if ( nRing != sp.getnRing() ) {
        std::cout << "Can't perform + on SpeciesArray: nRing exception: " << nRing << " != " << sp.getnRing() << std::endl;
        return *this;
    }
    
    if ( nTime != sp.getnTime() ) {
        std::cout << "Can't perform + on SpeciesArray: nTime exception: " << nTime << " != " << sp.getnTime() << std::endl;
        return *this;
    }

    for ( UInt iRing = 0; iRing < nRing; iRing++ ) {
        for ( UInt iTime = 0; iTime < nTime; iTime++ ) {

            CO2[iTime][iRing]      += sp.CO2[iTime][iRing];
            PPN[iTime][iRing]      += sp.PPN[iTime][iRing];
            BrNO2[iTime][iRing]    += sp.BrNO2[iTime][iRing];
            IEPOX[iTime][iRing]    += sp.IEPOX[iTime][iRing];
            PMNN[iTime][iRing]     += sp.PMNN[iTime][iRing];
            N2O[iTime][iRing]      += sp.N2O[iTime][iRing];
            N[iTime][iRing]        += sp.N[iTime][iRing];
            PAN[iTime][iRing]      += sp.PAN[iTime][iRing];
            ALK4[iTime][iRing]     += sp.ALK4[iTime][iRing];
            MAP[iTime][iRing]      += sp.MAP[iTime][iRing];
            MPN[iTime][iRing]      += sp.MPN[iTime][iRing];
            Cl2O2[iTime][iRing]    += sp.Cl2O2[iTime][iRing];
            ETP[iTime][iRing]      += sp.ETP[iTime][iRing];
            HNO2[iTime][iRing]     += sp.HNO2[iTime][iRing];
            C3H8[iTime][iRing]     += sp.C3H8[iTime][iRing];
            RA3P[iTime][iRing]     += sp.RA3P[iTime][iRing];
            RB3P[iTime][iRing]     += sp.RB3P[iTime][iRing];
            OClO[iTime][iRing]     += sp.OClO[iTime][iRing];
            ClNO2[iTime][iRing]    += sp.ClNO2[iTime][iRing];
            ISOP[iTime][iRing]     += sp.ISOP[iTime][iRing];
            HNO4[iTime][iRing]     += sp.HNO4[iTime][iRing];
            MAOP[iTime][iRing]     += sp.MAOP[iTime][iRing];
            MP[iTime][iRing]       += sp.MP[iTime][iRing];
            ClOO[iTime][iRing]     += sp.ClOO[iTime][iRing];
            RP[iTime][iRing]       += sp.RP[iTime][iRing];
            BrCl[iTime][iRing]     += sp.BrCl[iTime][iRing];
            PP[iTime][iRing]       += sp.PP[iTime][iRing];
            PRPN[iTime][iRing]     += sp.PRPN[iTime][iRing];
            SO4[iTime][iRing]      += sp.SO4[iTime][iRing];
            Br2[iTime][iRing]      += sp.Br2[iTime][iRing];
            ETHLN[iTime][iRing]    += sp.ETHLN[iTime][iRing];
            MVKN[iTime][iRing]     += sp.MVKN[iTime][iRing];
            R4P[iTime][iRing]      += sp.R4P[iTime][iRing];
            C2H6[iTime][iRing]     += sp.C2H6[iTime][iRing];
            RIP[iTime][iRing]      += sp.RIP[iTime][iRing];
            VRP[iTime][iRing]      += sp.VRP[iTime][iRing];
            ATOOH[iTime][iRing]    += sp.ATOOH[iTime][iRing];
            IAP[iTime][iRing]      += sp.IAP[iTime][iRing];
            DHMOB[iTime][iRing]    += sp.DHMOB[iTime][iRing];
            MOBA[iTime][iRing]     += sp.MOBA[iTime][iRing];
            MRP[iTime][iRing]      += sp.MRP[iTime][iRing];
            N2O5[iTime][iRing]     += sp.N2O5[iTime][iRing];
            ISNOHOO[iTime][iRing]  += sp.ISNOHOO[iTime][iRing];
            ISNP[iTime][iRing]     += sp.ISNP[iTime][iRing];
            ISOPNB[iTime][iRing]   += sp.ISOPNB[iTime][iRing];
            IEPOXOO[iTime][iRing]  += sp.IEPOXOO[iTime][iRing];
            MACRNO2[iTime][iRing]  += sp.MACRNO2[iTime][iRing];
            ROH[iTime][iRing]      += sp.ROH[iTime][iRing];
            MOBAOO[iTime][iRing]   += sp.MOBAOO[iTime][iRing];
            DIBOO[iTime][iRing]    += sp.DIBOO[iTime][iRing];
            PMN[iTime][iRing]      += sp.PMN[iTime][iRing];
            ISNOOB[iTime][iRing]   += sp.ISNOOB[iTime][iRing];
            INPN[iTime][iRing]     += sp.INPN[iTime][iRing];
            H[iTime][iRing]        += sp.H[iTime][iRing];
            BrNO3[iTime][iRing]    += sp.BrNO3[iTime][iRing];
            PRPE[iTime][iRing]     += sp.PRPE[iTime][iRing];
            MVKOO[iTime][iRing]    += sp.MVKOO[iTime][iRing];
            Cl2[iTime][iRing]      += sp.Cl2[iTime][iRing];
            ISOPND[iTime][iRing]   += sp.ISOPND[iTime][iRing];
            HOBr[iTime][iRing]     += sp.HOBr[iTime][iRing];
            A3O2[iTime][iRing]     += sp.A3O2[iTime][iRing];
            PROPNN[iTime][iRing]   += sp.PROPNN[iTime][iRing];
            GLYX[iTime][iRing]     += sp.GLYX[iTime][iRing];
            MAOPO2[iTime][iRing]   += sp.MAOPO2[iTime][iRing];
            CH4[iTime][iRing]      += sp.CH4[iTime][iRing];
            GAOO[iTime][iRing]     += sp.GAOO[iTime][iRing];
            B3O2[iTime][iRing]     += sp.B3O2[iTime][iRing];
            ACET[iTime][iRing]     += sp.ACET[iTime][iRing];
            MACRN[iTime][iRing]    += sp.MACRN[iTime][iRing];
            CH2OO[iTime][iRing]    += sp.CH2OO[iTime][iRing];
            MGLYOO[iTime][iRing]   += sp.MGLYOO[iTime][iRing];
            VRO2[iTime][iRing]     += sp.VRO2[iTime][iRing];
            MGLOO[iTime][iRing]    += sp.MGLOO[iTime][iRing];
            MACROO[iTime][iRing]   += sp.MACROO[iTime][iRing];
            PO2[iTime][iRing]      += sp.PO2[iTime][iRing];
            CH3CHOO[iTime][iRing]  += sp.CH3CHOO[iTime][iRing];
            MAN2[iTime][iRing]     += sp.MAN2[iTime][iRing];
            ISNOOA[iTime][iRing]   += sp.ISNOOA[iTime][iRing];
            H2O2[iTime][iRing]     += sp.H2O2[iTime][iRing];
            PRN1[iTime][iRing]     += sp.PRN1[iTime][iRing];
            ETO2[iTime][iRing]     += sp.ETO2[iTime][iRing];
            KO2[iTime][iRing]      += sp.KO2[iTime][iRing];
            RCO3[iTime][iRing]     += sp.RCO3[iTime][iRing];
            HC5OO[iTime][iRing]    += sp.HC5OO[iTime][iRing];
            GLYC[iTime][iRing]     += sp.GLYC[iTime][iRing];
            ClNO3[iTime][iRing]    += sp.ClNO3[iTime][iRing];
            RIO2[iTime][iRing]     += sp.RIO2[iTime][iRing];
            R4N1[iTime][iRing]     += sp.R4N1[iTime][iRing];
            HOCl[iTime][iRing]     += sp.HOCl[iTime][iRing];
            ATO2[iTime][iRing]     += sp.ATO2[iTime][iRing];
            HNO3[iTime][iRing]     += sp.HNO3[iTime][iRing];
            ISN1[iTime][iRing]     += sp.ISN1[iTime][iRing];
            MAO3[iTime][iRing]     += sp.MAO3[iTime][iRing];
            MRO2[iTime][iRing]     += sp.MRO2[iTime][iRing];
            INO2[iTime][iRing]     += sp.INO2[iTime][iRing];
            HAC[iTime][iRing]      += sp.HAC[iTime][iRing];
            HC5[iTime][iRing]      += sp.HC5[iTime][iRing];
            MGLY[iTime][iRing]     += sp.MGLY[iTime][iRing];
            ISOPNBO2[iTime][iRing] += sp.ISOPNBO2[iTime][iRing];
            ISOPNDO2[iTime][iRing] += sp.ISOPNDO2[iTime][iRing];
            R4O2[iTime][iRing]     += sp.R4O2[iTime][iRing];
            R4N2[iTime][iRing]     += sp.R4N2[iTime][iRing];
            BrO[iTime][iRing]      += sp.BrO[iTime][iRing];
            RCHO[iTime][iRing]     += sp.RCHO[iTime][iRing];
            MEK[iTime][iRing]      += sp.MEK[iTime][iRing];
            ClO[iTime][iRing]      += sp.ClO[iTime][iRing];
            MACR[iTime][iRing]     += sp.MACR[iTime][iRing];
            SO2[iTime][iRing]      += sp.SO2[iTime][iRing];
            MVK[iTime][iRing]      += sp.MVK[iTime][iRing];
            ALD2[iTime][iRing]     += sp.ALD2[iTime][iRing];
            MCO3[iTime][iRing]     += sp.MCO3[iTime][iRing];
            CH2O[iTime][iRing]     += sp.CH2O[iTime][iRing];
            H2O[iTime][iRing]      += sp.H2O[iTime][iRing];
            Br[iTime][iRing]       += sp.Br[iTime][iRing];
            NO[iTime][iRing]       += sp.NO[iTime][iRing];
            NO3[iTime][iRing]      += sp.NO3[iTime][iRing];
            Cl[iTime][iRing]       += sp.Cl[iTime][iRing];
            O[iTime][iRing]        += sp.O[iTime][iRing];
            O1D[iTime][iRing]      += sp.O1D[iTime][iRing];
            O3[iTime][iRing]       += sp.O3[iTime][iRing];
            HO2[iTime][iRing]      += sp.HO2[iTime][iRing];
            NO2[iTime][iRing]      += sp.NO2[iTime][iRing];
            OH[iTime][iRing]       += sp.OH[iTime][iRing];
            HBr[iTime][iRing]      += sp.HBr[iTime][iRing];
            HCl[iTime][iRing]      += sp.HCl[iTime][iRing];
            CO[iTime][iRing]       += sp.CO[iTime][iRing];
            MO2[iTime][iRing]      += sp.MO2[iTime][iRing];
            ACTA     += sp.ACTA;
            EOH      += sp.EOH;
            H2       += sp.H2;
            HCOOH    += sp.HCOOH;
            MOH      += sp.MOH;
            N2       += sp.N2;
            O2       += sp.O2;
            RCOOH    += sp.RCOOH;
   
            SO4L[iTime][iRing]  += sp.SO4L[iTime][iRing];
            H2OL[iTime][iRing]  += sp.H2OL[iTime][iRing];
            HNO3L[iTime][iRing] += sp.HNO3L[iTime][iRing];
            HClL[iTime][iRing]  += sp.HClL[iTime][iRing];
            HOClL[iTime][iRing] += sp.HOClL[iTime][iRing];
            HBrL[iTime][iRing]  += sp.HBrL[iTime][iRing];
            HOBrL[iTime][iRing] += sp.HOBrL[iTime][iRing];
            H2OS[iTime][iRing]  += sp.H2OS[iTime][iRing];
            HNO3S[iTime][iRing] += sp.HNO3S[iTime][iRing];

            sootDens[iTime][iRing] += sp.sootDens[iTime][iRing];
            sootRadi[iTime][iRing] += sp.sootRadi[iTime][iRing];
            sootArea[iTime][iRing] += sp.sootArea[iTime][iRing];
            iceDens[iTime][iRing]  += sp.iceDens[iTime][iRing];
            iceRadi[iTime][iRing]  += sp.iceRadi[iTime][iRing];
            iceArea[iTime][iRing]  += sp.iceArea[iTime][iRing];
            sulfDens[iTime][iRing] += sp.sulfDens[iTime][iRing];
            sulfRadi[iTime][iRing] += sp.sulfRadi[iTime][iRing];
            sulfArea[iTime][iRing] += sp.sulfArea[iTime][iRing];

        }
    }
    return *this;

} /* End of SpeciesArray::operator+ */

SpeciesArray& SpeciesArray::operator-( const SpeciesArray &sp )
{

    if ( nRing != sp.getnRing() ) {
        std::cout << "Can't perform + on SpeciesArray: nRing exception: " << nRing << " != " << sp.nRing << std::endl;
        return *this;
    }
    
    if ( nTime != sp.getnTime() ) {
        std::cout << "Can't perform + on SpeciesArray: nTime exception: " << nTime << " != " << sp.nTime << std::endl;
        return *this;
    }

    for ( UInt iRing = 0; iRing < nRing; iRing++ ) {
        for ( UInt iTime = 0; iTime < nTime; iTime++ ) {

            CO2[iTime][iRing]      -= sp.CO2[iTime][iRing];
            PPN[iTime][iRing]      -= sp.PPN[iTime][iRing];
            BrNO2[iTime][iRing]    -= sp.BrNO2[iTime][iRing];
            IEPOX[iTime][iRing]    -= sp.IEPOX[iTime][iRing];
            PMNN[iTime][iRing]     -= sp.PMNN[iTime][iRing];
            N2O[iTime][iRing]      -= sp.N2O[iTime][iRing];
            N[iTime][iRing]        -= sp.N[iTime][iRing];
            PAN[iTime][iRing]      -= sp.PAN[iTime][iRing];
            ALK4[iTime][iRing]     -= sp.ALK4[iTime][iRing];
            MAP[iTime][iRing]      -= sp.MAP[iTime][iRing];
            MPN[iTime][iRing]      -= sp.MPN[iTime][iRing];
            Cl2O2[iTime][iRing]    -= sp.Cl2O2[iTime][iRing];
            ETP[iTime][iRing]      -= sp.ETP[iTime][iRing];
            HNO2[iTime][iRing]     -= sp.HNO2[iTime][iRing];
            C3H8[iTime][iRing]     -= sp.C3H8[iTime][iRing];
            RA3P[iTime][iRing]     -= sp.RA3P[iTime][iRing];
            RB3P[iTime][iRing]     -= sp.RB3P[iTime][iRing];
            OClO[iTime][iRing]     -= sp.OClO[iTime][iRing];
            ClNO2[iTime][iRing]    -= sp.ClNO2[iTime][iRing];
            ISOP[iTime][iRing]     -= sp.ISOP[iTime][iRing];
            HNO4[iTime][iRing]     -= sp.HNO4[iTime][iRing];
            MAOP[iTime][iRing]     -= sp.MAOP[iTime][iRing];
            MP[iTime][iRing]       -= sp.MP[iTime][iRing];
            ClOO[iTime][iRing]     -= sp.ClOO[iTime][iRing];
            RP[iTime][iRing]       -= sp.RP[iTime][iRing];
            BrCl[iTime][iRing]     -= sp.BrCl[iTime][iRing];
            PP[iTime][iRing]       -= sp.PP[iTime][iRing];
            PRPN[iTime][iRing]     -= sp.PRPN[iTime][iRing];
            SO4[iTime][iRing]      -= sp.SO4[iTime][iRing];
            Br2[iTime][iRing]      -= sp.Br2[iTime][iRing];
            ETHLN[iTime][iRing]    -= sp.ETHLN[iTime][iRing];
            MVKN[iTime][iRing]     -= sp.MVKN[iTime][iRing];
            R4P[iTime][iRing]      -= sp.R4P[iTime][iRing];
            C2H6[iTime][iRing]     -= sp.C2H6[iTime][iRing];
            RIP[iTime][iRing]      -= sp.RIP[iTime][iRing];
            VRP[iTime][iRing]      -= sp.VRP[iTime][iRing];
            ATOOH[iTime][iRing]    -= sp.ATOOH[iTime][iRing];
            IAP[iTime][iRing]      -= sp.IAP[iTime][iRing];
            DHMOB[iTime][iRing]    -= sp.DHMOB[iTime][iRing];
            MOBA[iTime][iRing]     -= sp.MOBA[iTime][iRing];
            MRP[iTime][iRing]      -= sp.MRP[iTime][iRing];
            N2O5[iTime][iRing]     -= sp.N2O5[iTime][iRing];
            ISNOHOO[iTime][iRing]  -= sp.ISNOHOO[iTime][iRing];
            ISNP[iTime][iRing]     -= sp.ISNP[iTime][iRing];
            ISOPNB[iTime][iRing]   -= sp.ISOPNB[iTime][iRing];
            IEPOXOO[iTime][iRing]  -= sp.IEPOXOO[iTime][iRing];
            MACRNO2[iTime][iRing]  -= sp.MACRNO2[iTime][iRing];
            ROH[iTime][iRing]      -= sp.ROH[iTime][iRing];
            MOBAOO[iTime][iRing]   -= sp.MOBAOO[iTime][iRing];
            DIBOO[iTime][iRing]    -= sp.DIBOO[iTime][iRing];
            PMN[iTime][iRing]      -= sp.PMN[iTime][iRing];
            ISNOOB[iTime][iRing]   -= sp.ISNOOB[iTime][iRing];
            INPN[iTime][iRing]     -= sp.INPN[iTime][iRing];
            H[iTime][iRing]        -= sp.H[iTime][iRing];
            BrNO3[iTime][iRing]    -= sp.BrNO3[iTime][iRing];
            PRPE[iTime][iRing]     -= sp.PRPE[iTime][iRing];
            MVKOO[iTime][iRing]    -= sp.MVKOO[iTime][iRing];
            Cl2[iTime][iRing]      -= sp.Cl2[iTime][iRing];
            ISOPND[iTime][iRing]   -= sp.ISOPND[iTime][iRing];
            HOBr[iTime][iRing]     -= sp.HOBr[iTime][iRing];
            A3O2[iTime][iRing]     -= sp.A3O2[iTime][iRing];
            PROPNN[iTime][iRing]   -= sp.PROPNN[iTime][iRing];
            GLYX[iTime][iRing]     -= sp.GLYX[iTime][iRing];
            MAOPO2[iTime][iRing]   -= sp.MAOPO2[iTime][iRing];
            CH4[iTime][iRing]      -= sp.CH4[iTime][iRing];
            GAOO[iTime][iRing]     -= sp.GAOO[iTime][iRing];
            B3O2[iTime][iRing]     -= sp.B3O2[iTime][iRing];
            ACET[iTime][iRing]     -= sp.ACET[iTime][iRing];
            MACRN[iTime][iRing]    -= sp.MACRN[iTime][iRing];
            CH2OO[iTime][iRing]    -= sp.CH2OO[iTime][iRing];
            MGLYOO[iTime][iRing]   -= sp.MGLYOO[iTime][iRing];
            VRO2[iTime][iRing]     -= sp.VRO2[iTime][iRing];
            MGLOO[iTime][iRing]    -= sp.MGLOO[iTime][iRing];
            MACROO[iTime][iRing]   -= sp.MACROO[iTime][iRing];
            PO2[iTime][iRing]      -= sp.PO2[iTime][iRing];
            CH3CHOO[iTime][iRing]  -= sp.CH3CHOO[iTime][iRing];
            MAN2[iTime][iRing]     -= sp.MAN2[iTime][iRing];
            ISNOOA[iTime][iRing]   -= sp.ISNOOA[iTime][iRing];
            H2O2[iTime][iRing]     -= sp.H2O2[iTime][iRing];
            PRN1[iTime][iRing]     -= sp.PRN1[iTime][iRing];
            ETO2[iTime][iRing]     -= sp.ETO2[iTime][iRing];
            KO2[iTime][iRing]      -= sp.KO2[iTime][iRing];
            RCO3[iTime][iRing]     -= sp.RCO3[iTime][iRing];
            HC5OO[iTime][iRing]    -= sp.HC5OO[iTime][iRing];
            GLYC[iTime][iRing]     -= sp.GLYC[iTime][iRing];
            ClNO3[iTime][iRing]    -= sp.ClNO3[iTime][iRing];
            RIO2[iTime][iRing]     -= sp.RIO2[iTime][iRing];
            R4N1[iTime][iRing]     -= sp.R4N1[iTime][iRing];
            HOCl[iTime][iRing]     -= sp.HOCl[iTime][iRing];
            ATO2[iTime][iRing]     -= sp.ATO2[iTime][iRing];
            HNO3[iTime][iRing]     -= sp.HNO3[iTime][iRing];
            ISN1[iTime][iRing]     -= sp.ISN1[iTime][iRing];
            MAO3[iTime][iRing]     -= sp.MAO3[iTime][iRing];
            MRO2[iTime][iRing]     -= sp.MRO2[iTime][iRing];
            INO2[iTime][iRing]     -= sp.INO2[iTime][iRing];
            HAC[iTime][iRing]      -= sp.HAC[iTime][iRing];
            HC5[iTime][iRing]      -= sp.HC5[iTime][iRing];
            MGLY[iTime][iRing]     -= sp.MGLY[iTime][iRing];
            ISOPNBO2[iTime][iRing] -= sp.ISOPNBO2[iTime][iRing];
            ISOPNDO2[iTime][iRing] -= sp.ISOPNDO2[iTime][iRing];
            R4O2[iTime][iRing]     -= sp.R4O2[iTime][iRing];
            R4N2[iTime][iRing]     -= sp.R4N2[iTime][iRing];
            BrO[iTime][iRing]      -= sp.BrO[iTime][iRing];
            RCHO[iTime][iRing]     -= sp.RCHO[iTime][iRing];
            MEK[iTime][iRing]      -= sp.MEK[iTime][iRing];
            ClO[iTime][iRing]      -= sp.ClO[iTime][iRing];
            MACR[iTime][iRing]     -= sp.MACR[iTime][iRing];
            SO2[iTime][iRing]      -= sp.SO2[iTime][iRing];
            MVK[iTime][iRing]      -= sp.MVK[iTime][iRing];
            ALD2[iTime][iRing]     -= sp.ALD2[iTime][iRing];
            MCO3[iTime][iRing]     -= sp.MCO3[iTime][iRing];
            CH2O[iTime][iRing]     -= sp.CH2O[iTime][iRing];
            H2O[iTime][iRing]      -= sp.H2O[iTime][iRing];
            Br[iTime][iRing]       -= sp.Br[iTime][iRing];
            NO[iTime][iRing]       -= sp.NO[iTime][iRing];
            NO3[iTime][iRing]      -= sp.NO3[iTime][iRing];
            Cl[iTime][iRing]       -= sp.Cl[iTime][iRing];
            O[iTime][iRing]        -= sp.O[iTime][iRing];
            O1D[iTime][iRing]      -= sp.O1D[iTime][iRing];
            O3[iTime][iRing]       -= sp.O3[iTime][iRing];
            HO2[iTime][iRing]      -= sp.HO2[iTime][iRing];
            NO2[iTime][iRing]      -= sp.NO2[iTime][iRing];
            OH[iTime][iRing]       -= sp.OH[iTime][iRing];
            HBr[iTime][iRing]      -= sp.HBr[iTime][iRing];
            HCl[iTime][iRing]      -= sp.HCl[iTime][iRing];
            CO[iTime][iRing]       -= sp.CO[iTime][iRing];
            MO2[iTime][iRing]      -= sp.MO2[iTime][iRing];
            ACTA                   -= sp.ACTA;
            EOH                    -= sp.EOH;
            H2                     -= sp.H2;
            HCOOH                  -= sp.HCOOH;
            MOH                    -= sp.MOH;
            N2                     -= sp.N2;
            O2                     -= sp.O2;
            RCOOH                  -= sp.RCOOH;
            
            SO4L[iTime][iRing]     -= sp.SO4L[iTime][iRing];
            H2OL[iTime][iRing]     -= sp.H2OL[iTime][iRing];
            HNO3L[iTime][iRing]    -= sp.HNO3L[iTime][iRing];
            HClL[iTime][iRing]     -= sp.HClL[iTime][iRing];
            HOClL[iTime][iRing]    -= sp.HOClL[iTime][iRing];
            HBrL[iTime][iRing]     -= sp.HBrL[iTime][iRing];
            HOBrL[iTime][iRing]    -= sp.HOBrL[iTime][iRing];
            H2OS[iTime][iRing]     -= sp.H2OS[iTime][iRing];
            HNO3S[iTime][iRing]    -= sp.HNO3S[iTime][iRing];
            
            sootDens[iTime][iRing] -= sp.sootDens[iTime][iRing];
            sootRadi[iTime][iRing] -= sp.sootRadi[iTime][iRing];
            sootArea[iTime][iRing] -= sp.sootArea[iTime][iRing];
            iceDens[iTime][iRing]  -= sp.iceDens[iTime][iRing];
            iceRadi[iTime][iRing]  -= sp.iceRadi[iTime][iRing];
            iceArea[iTime][iRing]  -= sp.iceArea[iTime][iRing];
            sulfDens[iTime][iRing] -= sp.sulfDens[iTime][iRing];
            sulfRadi[iTime][iRing] -= sp.sulfRadi[iTime][iRing];
            sulfArea[iTime][iRing] -= sp.sulfArea[iTime][iRing];

        }
    }
    return *this;

} /* End of SpeciesArray::operator- */

SpeciesArray::~SpeciesArray( )
{
    /* Destructor */

} /* End of SpeciesArray::~SpeciesArray */

void SpeciesArray::FillIn( Solution &Data,             \
                           const Vector_3D &weights,   \
                           UInt nCounter )
{

    Vector_1D properties_LA( 4, 0.0E+00 );
    Vector_1D properties_PA( 4, 0.0E+00 );

    UInt jNy, iNx, iRing;
    RealDouble w = 0.0E+00;
    RealDouble totW = 0.0E+00;

    for ( iRing = 0; iRing < nRing; iRing++ ) {

        /* Precompute total weights */
        totW = 0.0E+00;
        for ( jNy = 0; jNy < Data.CO2.size(); jNy++ ) {
            for ( iNx = 0; iNx < Data.CO2[0].size(); iNx++ )
                totW += weights[iRing][jNy][iNx];
        }
        for ( jNy = 0; jNy < Data.CO2.size(); jNy++ ) {
            for ( iNx = 0; iNx < Data.CO2[0].size(); iNx++ ) {

                w = weights[iRing][jNy][iNx] / totW;

                CO2[nCounter][iRing]      += Data.CO2[jNy][iNx]      * w;
                PPN[nCounter][iRing]      += Data.PPN[jNy][iNx]      * w;
                BrNO2[nCounter][iRing]    += Data.BrNO2[jNy][iNx]    * w;
                IEPOX[nCounter][iRing]    += Data.IEPOX[jNy][iNx]    * w;
                PMNN[nCounter][iRing]     += Data.PMNN[jNy][iNx]     * w;
                N2O[nCounter][iRing]      += Data.N2O[jNy][iNx]      * w;
                N[nCounter][iRing]        += Data.N[jNy][iNx]        * w;
                PAN[nCounter][iRing]      += Data.PAN[jNy][iNx]      * w;
                ALK4[nCounter][iRing]     += Data.ALK4[jNy][iNx]     * w;
                MAP[nCounter][iRing]      += Data.MAP[jNy][iNx]      * w;
                MPN[nCounter][iRing]      += Data.MPN[jNy][iNx]      * w;
                Cl2O2[nCounter][iRing]    += Data.Cl2O2[jNy][iNx]    * w;
                ETP[nCounter][iRing]      += Data.ETP[jNy][iNx]      * w;
                HNO2[nCounter][iRing]     += Data.HNO2[jNy][iNx]     * w;
                C3H8[nCounter][iRing]     += Data.C3H8[jNy][iNx]     * w;
                RA3P[nCounter][iRing]     += Data.RA3P[jNy][iNx]     * w;
                RB3P[nCounter][iRing]     += Data.RB3P[jNy][iNx]     * w;
                OClO[nCounter][iRing]     += Data.OClO[jNy][iNx]     * w;
                ClNO2[nCounter][iRing]    += Data.ClNO2[jNy][iNx]    * w;
                ISOP[nCounter][iRing]     += Data.ISOP[jNy][iNx]     * w;
                HNO4[nCounter][iRing]     += Data.HNO4[jNy][iNx]     * w;
                MAOP[nCounter][iRing]     += Data.MAOP[jNy][iNx]     * w;
                MP[nCounter][iRing]       += Data.MP[jNy][iNx]       * w;
                ClOO[nCounter][iRing]     += Data.ClOO[jNy][iNx]     * w;
                RP[nCounter][iRing]       += Data.RP[jNy][iNx]       * w;
                BrCl[nCounter][iRing]     += Data.BrCl[jNy][iNx]     * w;
                PP[nCounter][iRing]       += Data.PP[jNy][iNx]       * w;
                PRPN[nCounter][iRing]     += Data.PRPN[jNy][iNx]     * w;
                SO4[nCounter][iRing]      += Data.SO4[jNy][iNx]      * w;
                Br2[nCounter][iRing]      += Data.Br2[jNy][iNx]      * w;
                ETHLN[nCounter][iRing]    += Data.ETHLN[jNy][iNx]    * w;
                MVKN[nCounter][iRing]     += Data.MVKN[jNy][iNx]     * w;
                R4P[nCounter][iRing]      += Data.R4P[jNy][iNx]      * w;
                C2H6[nCounter][iRing]     += Data.C2H6[jNy][iNx]     * w;
                RIP[nCounter][iRing]      += Data.RIP[jNy][iNx]      * w;
                VRP[nCounter][iRing]      += Data.VRP[jNy][iNx]      * w;
                ATOOH[nCounter][iRing]    += Data.ATOOH[jNy][iNx]    * w;
                IAP[nCounter][iRing]      += Data.IAP[jNy][iNx]      * w;
                DHMOB[nCounter][iRing]    += Data.DHMOB[jNy][iNx]    * w;
                MOBA[nCounter][iRing]     += Data.MOBA[jNy][iNx]     * w;
                MRP[nCounter][iRing]      += Data.MRP[jNy][iNx]      * w;
                N2O5[nCounter][iRing]     += Data.N2O5[jNy][iNx]     * w;
                ISNOHOO[nCounter][iRing]  += Data.ISNOHOO[jNy][iNx]  * w;
                ISNP[nCounter][iRing]     += Data.ISNP[jNy][iNx]     * w;
                ISOPNB[nCounter][iRing]   += Data.ISOPNB[jNy][iNx]   * w;
                IEPOXOO[nCounter][iRing]  += Data.IEPOXOO[jNy][iNx]  * w;
                MACRNO2[nCounter][iRing]  += Data.MACRNO2[jNy][iNx]  * w;
                ROH[nCounter][iRing]      += Data.ROH[jNy][iNx]      * w;
                MOBAOO[nCounter][iRing]   += Data.MOBAOO[jNy][iNx]   * w;
                DIBOO[nCounter][iRing]    += Data.DIBOO[jNy][iNx]    * w;
                PMN[nCounter][iRing]      += Data.PMN[jNy][iNx]      * w;
                ISNOOB[nCounter][iRing]   += Data.ISNOOB[jNy][iNx]   * w;
                INPN[nCounter][iRing]     += Data.INPN[jNy][iNx]     * w;
                H[nCounter][iRing]        += Data.H[jNy][iNx]        * w;
                BrNO3[nCounter][iRing]    += Data.BrNO3[jNy][iNx]    * w;
                PRPE[nCounter][iRing]     += Data.PRPE[jNy][iNx]     * w;
                MVKOO[nCounter][iRing]    += Data.MVKOO[jNy][iNx]    * w;
                Cl2[nCounter][iRing]      += Data.Cl2[jNy][iNx]      * w;
                ISOPND[nCounter][iRing]   += Data.ISOPND[jNy][iNx]   * w;
                HOBr[nCounter][iRing]     += Data.HOBr[jNy][iNx]     * w;
                A3O2[nCounter][iRing]     += Data.A3O2[jNy][iNx]     * w;
                PROPNN[nCounter][iRing]   += Data.PROPNN[jNy][iNx]   * w;
                GLYX[nCounter][iRing]     += Data.GLYX[jNy][iNx]     * w;
                MAOPO2[nCounter][iRing]   += Data.MAOPO2[jNy][iNx]   * w;
                CH4[nCounter][iRing]      += Data.CH4[jNy][iNx]      * w;
                GAOO[nCounter][iRing]     += Data.GAOO[jNy][iNx]     * w;
                B3O2[nCounter][iRing]     += Data.B3O2[jNy][iNx]     * w;
                ACET[nCounter][iRing]     += Data.ACET[jNy][iNx]     * w;
                MACRN[nCounter][iRing]    += Data.MACRN[jNy][iNx]    * w;
                CH2OO[nCounter][iRing]    += Data.CH2OO[jNy][iNx]    * w;
                MGLYOO[nCounter][iRing]   += Data.MGLYOO[jNy][iNx]   * w;
                VRO2[nCounter][iRing]     += Data.VRO2[jNy][iNx]     * w;
                MGLOO[nCounter][iRing]    += Data.MGLOO[jNy][iNx]    * w;
                MACROO[nCounter][iRing]   += Data.MACROO[jNy][iNx]   * w;
                PO2[nCounter][iRing]      += Data.PO2[jNy][iNx]      * w;
                CH3CHOO[nCounter][iRing]  += Data.CH3CHOO[jNy][iNx]  * w;
                MAN2[nCounter][iRing]     += Data.MAN2[jNy][iNx]     * w;
                ISNOOA[nCounter][iRing]   += Data.ISNOOA[jNy][iNx]   * w;
                H2O2[nCounter][iRing]     += Data.H2O2[jNy][iNx]     * w;
                PRN1[nCounter][iRing]     += Data.PRN1[jNy][iNx]     * w;
                ETO2[nCounter][iRing]     += Data.ETO2[jNy][iNx]     * w;
                KO2[nCounter][iRing]      += Data.KO2[jNy][iNx]      * w;
                RCO3[nCounter][iRing]     += Data.RCO3[jNy][iNx]     * w;
                HC5OO[nCounter][iRing]    += Data.HC5OO[jNy][iNx]    * w;
                GLYC[nCounter][iRing]     += Data.GLYC[jNy][iNx]     * w;
                ClNO3[nCounter][iRing]    += Data.ClNO3[jNy][iNx]    * w;
                RIO2[nCounter][iRing]     += Data.RIO2[jNy][iNx]     * w;
                R4N1[nCounter][iRing]     += Data.R4N1[jNy][iNx]     * w;
                HOCl[nCounter][iRing]     += Data.HOCl[jNy][iNx]     * w;
                ATO2[nCounter][iRing]     += Data.ATO2[jNy][iNx]     * w;
                HNO3[nCounter][iRing]     += Data.HNO3[jNy][iNx]     * w;
                ISN1[nCounter][iRing]     += Data.ISN1[jNy][iNx]     * w;
                MAO3[nCounter][iRing]     += Data.MAO3[jNy][iNx]     * w;
                MRO2[nCounter][iRing]     += Data.MRO2[jNy][iNx]     * w;
                INO2[nCounter][iRing]     += Data.INO2[jNy][iNx]     * w;
                HAC[nCounter][iRing]      += Data.HAC[jNy][iNx]      * w;
                HC5[nCounter][iRing]      += Data.HC5[jNy][iNx]      * w;
                MGLY[nCounter][iRing]     += Data.MGLY[jNy][iNx]     * w;
                ISOPNBO2[nCounter][iRing] += Data.ISOPNBO2[jNy][iNx] * w;
                ISOPNDO2[nCounter][iRing] += Data.ISOPNDO2[jNy][iNx] * w;
                R4O2[nCounter][iRing]     += Data.R4O2[jNy][iNx]     * w;
                R4N2[nCounter][iRing]     += Data.R4N2[jNy][iNx]     * w;
                BrO[nCounter][iRing]      += Data.BrO[jNy][iNx]      * w;
                RCHO[nCounter][iRing]     += Data.RCHO[jNy][iNx]     * w;
                MEK[nCounter][iRing]      += Data.MEK[jNy][iNx]      * w;
                ClO[nCounter][iRing]      += Data.ClO[jNy][iNx]      * w;
                MACR[nCounter][iRing]     += Data.MACR[jNy][iNx]     * w;
                SO2[nCounter][iRing]      += Data.SO2[jNy][iNx]      * w;
                MVK[nCounter][iRing]      += Data.MVK[jNy][iNx]      * w;
                ALD2[nCounter][iRing]     += Data.ALD2[jNy][iNx]     * w;
                MCO3[nCounter][iRing]     += Data.MCO3[jNy][iNx]     * w;
                CH2O[nCounter][iRing]     += Data.CH2O[jNy][iNx]     * w;
                H2O[nCounter][iRing]      += Data.H2O[jNy][iNx]      * w;
                Br[nCounter][iRing]       += Data.Br[jNy][iNx]       * w;
                NO[nCounter][iRing]       += Data.NO[jNy][iNx]       * w;
                NO3[nCounter][iRing]      += Data.NO3[jNy][iNx]      * w;
                Cl[nCounter][iRing]       += Data.Cl[jNy][iNx]       * w;
                O[nCounter][iRing]        += Data.O[jNy][iNx]        * w;
                O1D[nCounter][iRing]      += Data.O1D[jNy][iNx]      * w;
                O3[nCounter][iRing]       += Data.O3[jNy][iNx]       * w;
                HO2[nCounter][iRing]      += Data.HO2[jNy][iNx]      * w;
                NO2[nCounter][iRing]      += Data.NO2[jNy][iNx]      * w;
                OH[nCounter][iRing]       += Data.OH[jNy][iNx]       * w;
                HBr[nCounter][iRing]      += Data.HBr[jNy][iNx]      * w;
                HCl[nCounter][iRing]      += Data.HCl[jNy][iNx]      * w;
                CO[nCounter][iRing]       += Data.CO[jNy][iNx]       * w;
                MO2[nCounter][iRing]      += Data.MO2[jNy][iNx]      * w;

                SO4L[nCounter][iRing]     += Data.SO4L[jNy][iNx]     * w;
                H2OL[nCounter][iRing]     += Data.H2OL[jNy][iNx]     * w;
                HNO3L[nCounter][iRing]    += Data.HNO3L[jNy][iNx]    * w;
                HClL[nCounter][iRing]     += Data.HClL[jNy][iNx]     * w;
                HOClL[nCounter][iRing]    += Data.HOClL[jNy][iNx]    * w;
                HBrL[nCounter][iRing]     += Data.HBrL[jNy][iNx]     * w;
                HOBrL[nCounter][iRing]    += Data.HOBrL[jNy][iNx]    * w;
                H2OS[nCounter][iRing]     += Data.H2OS[jNy][iNx]     * w;
                HNO3S[nCounter][iRing]    += Data.HNO3S[jNy][iNx]    * w;

                sootDens[nCounter][iRing] += Data.sootDens[jNy][iNx] * w;
                sootRadi[nCounter][iRing] += Data.sootRadi[jNy][iNx] * w;
                sootArea[nCounter][iRing] += Data.sootArea[jNy][iNx] * w;

            }
        }

        properties_LA = Data.liquidAerosol.Average( weights[iRing], \
                                                    totW );
        sulfDens[nCounter][iRing] = properties_LA[0];
        sulfRadi[nCounter][iRing] = properties_LA[1];
        sulfArea[nCounter][iRing] = properties_LA[2];

        properties_PA = Data.solidAerosol.Average( weights[iRing], \
                                                   totW );
        iceDens[nCounter][iRing] = properties_PA[0];
        iceRadi[nCounter][iRing] = properties_PA[1];
        iceArea[nCounter][iRing] = properties_PA[2];

    } 

    ACTA  = Data.ACTA[0][0];
    EOH   = Data.EOH[0][0];
    H2    = Data.H2[0][0];
    HCOOH = Data.HCOOH[0][0];
    MOH   = Data.MOH[0][0];
    N2    = Data.N2[0][0];
    O2    = Data.O2[0][0];
    RCOOH = Data.RCOOH[0][0];

} /* End of SpeciesArray::FillIn */

void SpeciesArray::FillIn( RealDouble varArray[], UInt iTime, UInt iRing )
{

    /* Ensure positiveness */
    for ( UInt i = 0; i < NVAR; i++ ) {
        if ( varArray[i] <= 0.0 ) {
            varArray[i] = ZERO;
        }
    }

    CO2[iTime][iRing]      = varArray[  0];
    PPN[iTime][iRing]      = varArray[  1];
    BrNO2[iTime][iRing]    = varArray[  2];
    IEPOX[iTime][iRing]    = varArray[  3];
    PMNN[iTime][iRing]     = varArray[  4];
    N2O[iTime][iRing]      = varArray[  5];
    N[iTime][iRing]        = varArray[  6];
    PAN[iTime][iRing]      = varArray[  7];
    ALK4[iTime][iRing]     = varArray[  8];
    MAP[iTime][iRing]      = varArray[  9];
    MPN[iTime][iRing]      = varArray[ 10];
    Cl2O2[iTime][iRing]    = varArray[ 11];
    ETP[iTime][iRing]      = varArray[ 12];
    HNO2[iTime][iRing]     = varArray[ 13];
    C3H8[iTime][iRing]     = varArray[ 14];
    RA3P[iTime][iRing]     = varArray[ 15];
    RB3P[iTime][iRing]     = varArray[ 16];
    OClO[iTime][iRing]     = varArray[ 17];
    ClNO2[iTime][iRing]    = varArray[ 18];
    ISOP[iTime][iRing]     = varArray[ 19];
    HNO4[iTime][iRing]     = varArray[ 20];
    MAOP[iTime][iRing]     = varArray[ 21];
    MP[iTime][iRing]       = varArray[ 22];
    ClOO[iTime][iRing]     = varArray[ 23];
    RP[iTime][iRing]       = varArray[ 24];
    BrCl[iTime][iRing]     = varArray[ 25];
    PP[iTime][iRing]       = varArray[ 26];
    PRPN[iTime][iRing]     = varArray[ 27];
    SO4[iTime][iRing]      = varArray[ 28];
    Br2[iTime][iRing]      = varArray[ 29];
    ETHLN[iTime][iRing]    = varArray[ 30];
    MVKN[iTime][iRing]     = varArray[ 31];
    R4P[iTime][iRing]      = varArray[ 32];
    C2H6[iTime][iRing]     = varArray[ 33];
    RIP[iTime][iRing]      = varArray[ 34];
    VRP[iTime][iRing]      = varArray[ 35];
    ATOOH[iTime][iRing]    = varArray[ 36];
    IAP[iTime][iRing]      = varArray[ 37];
    DHMOB[iTime][iRing]    = varArray[ 38];
    MOBA[iTime][iRing]     = varArray[ 39];
    MRP[iTime][iRing]      = varArray[ 40];
    N2O5[iTime][iRing]     = varArray[ 41];
    ISNOHOO[iTime][iRing]  = varArray[ 42];
    ISNP[iTime][iRing]     = varArray[ 43];
    ISOPNB[iTime][iRing]   = varArray[ 44];
    IEPOXOO[iTime][iRing]  = varArray[ 45];
    MACRNO2[iTime][iRing]  = varArray[ 46];
    ROH[iTime][iRing]      = varArray[ 47];
    MOBAOO[iTime][iRing]   = varArray[ 48];
    DIBOO[iTime][iRing]    = varArray[ 49];
    PMN[iTime][iRing]      = varArray[ 50];
    ISNOOB[iTime][iRing]   = varArray[ 51];
    INPN[iTime][iRing]     = varArray[ 52];
    H[iTime][iRing]        = varArray[ 53];
    BrNO3[iTime][iRing]    = varArray[ 54];
    PRPE[iTime][iRing]     = varArray[ 55];
    MVKOO[iTime][iRing]    = varArray[ 56];
    Cl2[iTime][iRing]      = varArray[ 57];
    ISOPND[iTime][iRing]   = varArray[ 58];
    HOBr[iTime][iRing]     = varArray[ 59];
    A3O2[iTime][iRing]     = varArray[ 60];
    PROPNN[iTime][iRing]   = varArray[ 61];
    GLYX[iTime][iRing]     = varArray[ 62];
    MAOPO2[iTime][iRing]   = varArray[ 63];
    CH4[iTime][iRing]      = varArray[ 64];
    GAOO[iTime][iRing]     = varArray[ 65];
    B3O2[iTime][iRing]     = varArray[ 66];
    ACET[iTime][iRing]     = varArray[ 67];
    MACRN[iTime][iRing]    = varArray[ 68];
    CH2OO[iTime][iRing]    = varArray[ 69];
    MGLYOO[iTime][iRing]   = varArray[ 70];
    VRO2[iTime][iRing]     = varArray[ 71];
    MGLOO[iTime][iRing]    = varArray[ 72];
    MACROO[iTime][iRing]   = varArray[ 73];
    PO2[iTime][iRing]      = varArray[ 74];
    CH3CHOO[iTime][iRing]  = varArray[ 75];
    MAN2[iTime][iRing]     = varArray[ 76];
    ISNOOA[iTime][iRing]   = varArray[ 77];
    H2O2[iTime][iRing]     = varArray[ 78];
    PRN1[iTime][iRing]     = varArray[ 79];
    ETO2[iTime][iRing]     = varArray[ 80];
    KO2[iTime][iRing]      = varArray[ 81];
    RCO3[iTime][iRing]     = varArray[ 82];
    HC5OO[iTime][iRing]    = varArray[ 83];
    GLYC[iTime][iRing]     = varArray[ 84];
    ClNO3[iTime][iRing]    = varArray[ 85];
    RIO2[iTime][iRing]     = varArray[ 86];
    R4N1[iTime][iRing]     = varArray[ 87];
    HOCl[iTime][iRing]     = varArray[ 88];
    ATO2[iTime][iRing]     = varArray[ 89];
    HNO3[iTime][iRing]     = varArray[ 90];
    ISN1[iTime][iRing]     = varArray[ 91];
    MAO3[iTime][iRing]     = varArray[ 92];
    MRO2[iTime][iRing]     = varArray[ 93];
    INO2[iTime][iRing]     = varArray[ 94];
    HAC[iTime][iRing]      = varArray[ 95];
    HC5[iTime][iRing]      = varArray[ 96];
    MGLY[iTime][iRing]     = varArray[ 97];
    ISOPNBO2[iTime][iRing] = varArray[ 98];
    ISOPNDO2[iTime][iRing] = varArray[ 99];
    R4O2[iTime][iRing]     = varArray[100];
    R4N2[iTime][iRing]     = varArray[101];
    BrO[iTime][iRing]      = varArray[102];
    RCHO[iTime][iRing]     = varArray[103];
    MEK[iTime][iRing]      = varArray[104];
    ClO[iTime][iRing]      = varArray[105];
    MACR[iTime][iRing]     = varArray[106];
    SO2[iTime][iRing]      = varArray[107];
    MVK[iTime][iRing]      = varArray[108];
    ALD2[iTime][iRing]     = varArray[109];
    MCO3[iTime][iRing]     = varArray[110];
    CH2O[iTime][iRing]     = varArray[111];
    H2O[iTime][iRing]      = varArray[112];
    Br[iTime][iRing]       = varArray[113];
    NO[iTime][iRing]       = varArray[114];
    NO3[iTime][iRing]      = varArray[115];
    Cl[iTime][iRing]       = varArray[116];
    O[iTime][iRing]        = varArray[117];
    O1D[iTime][iRing]      = varArray[118];
    O3[iTime][iRing]       = varArray[119];
    HO2[iTime][iRing]      = varArray[120];
    NO2[iTime][iRing]      = varArray[121];
    OH[iTime][iRing]       = varArray[122];
    HBr[iTime][iRing]      = varArray[123];
    HCl[iTime][iRing]      = varArray[124];
    CO[iTime][iRing]       = varArray[125];
    MO2[iTime][iRing]      = varArray[126];


} /* End of SpeciesArray::FillIn */

void SpeciesArray::getData( RealDouble varArray[], RealDouble fixArray[], UInt iTime, UInt iRing )
{

    varArray[  0] = CO2[iTime][iRing];
    varArray[  1] = PPN[iTime][iRing];
    varArray[  2] = BrNO2[iTime][iRing];
    varArray[  3] = IEPOX[iTime][iRing];
    varArray[  4] = PMNN[iTime][iRing];
    varArray[  5] = N2O[iTime][iRing];
    varArray[  6] = N[iTime][iRing];
    varArray[  7] = PAN[iTime][iRing];
    varArray[  8] = ALK4[iTime][iRing];
    varArray[  9] = MAP[iTime][iRing];
    varArray[ 10] = MPN[iTime][iRing];
    varArray[ 11] = Cl2O2[iTime][iRing];
    varArray[ 12] = ETP[iTime][iRing];
    varArray[ 13] = HNO2[iTime][iRing];
    varArray[ 14] = C3H8[iTime][iRing];
    varArray[ 15] = RA3P[iTime][iRing];
    varArray[ 16] = RB3P[iTime][iRing];
    varArray[ 17] = OClO[iTime][iRing];
    varArray[ 18] = ClNO2[iTime][iRing];
    varArray[ 19] = ISOP[iTime][iRing];
    varArray[ 20] = HNO4[iTime][iRing];
    varArray[ 21] = MAOP[iTime][iRing];
    varArray[ 22] = MP[iTime][iRing];
    varArray[ 23] = ClOO[iTime][iRing];
    varArray[ 24] = RP[iTime][iRing];
    varArray[ 25] = BrCl[iTime][iRing];
    varArray[ 26] = PP[iTime][iRing];
    varArray[ 27] = PRPN[iTime][iRing];
    varArray[ 28] = SO4[iTime][iRing];
    varArray[ 29] = Br2[iTime][iRing];
    varArray[ 30] = ETHLN[iTime][iRing];
    varArray[ 31] = MVKN[iTime][iRing];
    varArray[ 32] = R4P[iTime][iRing];
    varArray[ 33] = C2H6[iTime][iRing];
    varArray[ 34] = RIP[iTime][iRing];
    varArray[ 35] = VRP[iTime][iRing];
    varArray[ 36] = ATOOH[iTime][iRing];
    varArray[ 37] = IAP[iTime][iRing];
    varArray[ 38] = DHMOB[iTime][iRing];
    varArray[ 39] = MOBA[iTime][iRing];
    varArray[ 40] = MRP[iTime][iRing];
    varArray[ 41] = N2O5[iTime][iRing];
    varArray[ 42] = ISNOHOO[iTime][iRing];
    varArray[ 43] = ISNP[iTime][iRing];
    varArray[ 44] = ISOPNB[iTime][iRing];
    varArray[ 45] = IEPOXOO[iTime][iRing];
    varArray[ 46] = MACRNO2[iTime][iRing];
    varArray[ 47] = ROH[iTime][iRing];
    varArray[ 48] = MOBAOO[iTime][iRing];
    varArray[ 49] = DIBOO[iTime][iRing];
    varArray[ 50] = PMN[iTime][iRing];
    varArray[ 51] = ISNOOB[iTime][iRing];
    varArray[ 52] = INPN[iTime][iRing];
    varArray[ 53] = H[iTime][iRing];
    varArray[ 54] = BrNO3[iTime][iRing];
    varArray[ 55] = PRPE[iTime][iRing];
    varArray[ 56] = MVKOO[iTime][iRing];
    varArray[ 57] = Cl2[iTime][iRing];
    varArray[ 58] = ISOPND[iTime][iRing];
    varArray[ 59] = HOBr[iTime][iRing];
    varArray[ 60] = A3O2[iTime][iRing];
    varArray[ 61] = PROPNN[iTime][iRing];
    varArray[ 62] = GLYX[iTime][iRing];
    varArray[ 63] = MAOPO2[iTime][iRing];
    varArray[ 64] = CH4[iTime][iRing];
    varArray[ 65] = GAOO[iTime][iRing];
    varArray[ 66] = B3O2[iTime][iRing];
    varArray[ 67] = ACET[iTime][iRing];
    varArray[ 68] = MACRN[iTime][iRing];
    varArray[ 69] = CH2OO[iTime][iRing];
    varArray[ 70] = MGLYOO[iTime][iRing];
    varArray[ 71] = VRO2[iTime][iRing];
    varArray[ 72] = MGLOO[iTime][iRing];
    varArray[ 73] = MACROO[iTime][iRing];
    varArray[ 74] = PO2[iTime][iRing];
    varArray[ 75] = CH3CHOO[iTime][iRing];
    varArray[ 76] = MAN2[iTime][iRing];
    varArray[ 77] = ISNOOA[iTime][iRing];
    varArray[ 78] = H2O2[iTime][iRing];
    varArray[ 79] = PRN1[iTime][iRing];
    varArray[ 80] = ETO2[iTime][iRing];
    varArray[ 81] = KO2[iTime][iRing];
    varArray[ 82] = RCO3[iTime][iRing];
    varArray[ 83] = HC5OO[iTime][iRing];
    varArray[ 84] = GLYC[iTime][iRing];
    varArray[ 85] = ClNO3[iTime][iRing];
    varArray[ 86] = RIO2[iTime][iRing];
    varArray[ 87] = R4N1[iTime][iRing];
    varArray[ 88] = HOCl[iTime][iRing];
    varArray[ 89] = ATO2[iTime][iRing];
    varArray[ 90] = HNO3[iTime][iRing];
    varArray[ 91] = ISN1[iTime][iRing];
    varArray[ 92] = MAO3[iTime][iRing];
    varArray[ 93] = MRO2[iTime][iRing];
    varArray[ 94] = INO2[iTime][iRing];
    varArray[ 95] = HAC[iTime][iRing];
    varArray[ 96] = HC5[iTime][iRing];
    varArray[ 97] = MGLY[iTime][iRing];
    varArray[ 98] = ISOPNBO2[iTime][iRing];
    varArray[ 99] = ISOPNDO2[iTime][iRing];
    varArray[100] = R4O2[iTime][iRing];
    varArray[101] = R4N2[iTime][iRing];
    varArray[102] = BrO[iTime][iRing];
    varArray[103] = RCHO[iTime][iRing];
    varArray[104] = MEK[iTime][iRing];
    varArray[105] = ClO[iTime][iRing];
    varArray[106] = MACR[iTime][iRing];
    varArray[107] = SO2[iTime][iRing];
    varArray[108] = MVK[iTime][iRing];
    varArray[109] = ALD2[iTime][iRing];
    varArray[110] = MCO3[iTime][iRing];
    varArray[111] = CH2O[iTime][iRing];
    varArray[112] = H2O[iTime][iRing];
    varArray[113] = Br[iTime][iRing];
    varArray[114] = NO[iTime][iRing];
    varArray[115] = NO3[iTime][iRing];
    varArray[116] = Cl[iTime][iRing];
    varArray[117] = O[iTime][iRing];
    varArray[118] = O1D[iTime][iRing];
    varArray[119] = O3[iTime][iRing];
    varArray[120] = HO2[iTime][iRing];
    varArray[121] = NO2[iTime][iRing];
    varArray[122] = OH[iTime][iRing];
    varArray[123] = HBr[iTime][iRing];
    varArray[124] = HCl[iTime][iRing];
    varArray[125] = CO[iTime][iRing];
    varArray[126] = MO2[iTime][iRing];
    fixArray[  0] = ACTA;
    fixArray[  1] = EOH;
    fixArray[  2] = H2;
    fixArray[  3] = HCOOH;
    fixArray[  4] = MOH;
    fixArray[  5] = N2;
    fixArray[  6] = O2;
    fixArray[  7] = RCOOH;

    /* Ensure positiveness */
    for ( UInt i = 0; i < NVAR; i++ ) {
        if ( varArray[i] <= 0.0 ) {
            varArray[i] = ZERO;
        }
    }

} /* End of SpeciesArray::getData */

Vector_1D SpeciesArray::RingAverage( const Vector_1D ringArea, const RealDouble totArea, \
                                     const UInt iNt ) const
{

    Vector_1D ringAverage( NVAR, 0.0E+00 );
    UInt iRing;
    RealDouble area;

    for ( iRing = 0; iRing < nRing; iRing++ ) {
        area = ringArea[iRing] / totArea;
        ringAverage[ind_CO2]      += CO2[iNt][iRing]      * area;
        ringAverage[ind_PPN]      += PPN[iNt][iRing]      * area;
        ringAverage[ind_BrNO2]    += BrNO2[iNt][iRing]    * area;
        ringAverage[ind_IEPOX]    += IEPOX[iNt][iRing]    * area;
        ringAverage[ind_PMNN]     += PMNN[iNt][iRing]     * area;
        ringAverage[ind_N2O]      += N2O[iNt][iRing]      * area;
        ringAverage[ind_N]        += N[iNt][iRing]        * area;
        ringAverage[ind_PAN]      += PAN[iNt][iRing]      * area;
        ringAverage[ind_ALK4]     += ALK4[iNt][iRing]     * area;
        ringAverage[ind_MAP]      += MAP[iNt][iRing]      * area;
        ringAverage[ind_MPN]      += MPN[iNt][iRing]      * area;
        ringAverage[ind_Cl2O2]    += Cl2O2[iNt][iRing]    * area;
        ringAverage[ind_ETP]      += ETP[iNt][iRing]      * area;
        ringAverage[ind_HNO2]     += HNO2[iNt][iRing]     * area;
        ringAverage[ind_C3H8]     += C3H8[iNt][iRing]     * area;
        ringAverage[ind_RA3P]     += RA3P[iNt][iRing]     * area;
        ringAverage[ind_RB3P]     += RB3P[iNt][iRing]     * area;
        ringAverage[ind_OClO]     += OClO[iNt][iRing]     * area;
        ringAverage[ind_ClNO2]    += ClNO2[iNt][iRing]    * area;
        ringAverage[ind_ISOP]     += ISOP[iNt][iRing]     * area;
        ringAverage[ind_HNO4]     += HNO4[iNt][iRing]     * area;
        ringAverage[ind_MAOP]     += MAOP[iNt][iRing]     * area;
        ringAverage[ind_MP]       += MP[iNt][iRing]       * area;
        ringAverage[ind_ClOO]     += ClOO[iNt][iRing]     * area;
        ringAverage[ind_RP]       += RP[iNt][iRing]       * area;
        ringAverage[ind_BrCl]     += BrCl[iNt][iRing]     * area;
        ringAverage[ind_PP]       += PP[iNt][iRing]       * area;
        ringAverage[ind_PRPN]     += PRPN[iNt][iRing]     * area;
        ringAverage[ind_SO4]      += SO4[iNt][iRing]      * area;
        ringAverage[ind_Br2]      += Br2[iNt][iRing]      * area;
        ringAverage[ind_ETHLN]    += ETHLN[iNt][iRing]    * area;
        ringAverage[ind_MVKN]     += MVKN[iNt][iRing]     * area;
        ringAverage[ind_R4P]      += R4P[iNt][iRing]      * area;
        ringAverage[ind_C2H6]     += C2H6[iNt][iRing]     * area;
        ringAverage[ind_RIP]      += RIP[iNt][iRing]      * area;
        ringAverage[ind_VRP]      += VRP[iNt][iRing]      * area;
        ringAverage[ind_ATOOH]    += ATOOH[iNt][iRing]    * area;
        ringAverage[ind_IAP]      += IAP[iNt][iRing]      * area;
        ringAverage[ind_DHMOB]    += DHMOB[iNt][iRing]    * area;
        ringAverage[ind_MOBA]     += MOBA[iNt][iRing]     * area;
        ringAverage[ind_MRP]      += MRP[iNt][iRing]      * area;
        ringAverage[ind_N2O5]     += N2O5[iNt][iRing]     * area;
        ringAverage[ind_ISNOHOO]  += ISNOHOO[iNt][iRing]  * area;
        ringAverage[ind_ISNP]     += ISNP[iNt][iRing]     * area;
        ringAverage[ind_ISOPNB]   += ISOPNB[iNt][iRing]   * area;
        ringAverage[ind_IEPOXOO]  += IEPOXOO[iNt][iRing]  * area;
        ringAverage[ind_MACRNO2]  += MACRNO2[iNt][iRing]  * area;
        ringAverage[ind_ROH]      += ROH[iNt][iRing]      * area;
        ringAverage[ind_MOBAOO]   += MOBAOO[iNt][iRing]   * area;
        ringAverage[ind_DIBOO]    += DIBOO[iNt][iRing]    * area;
        ringAverage[ind_PMN]      += PMN[iNt][iRing]      * area;
        ringAverage[ind_ISNOOB]   += ISNOOB[iNt][iRing]   * area;
        ringAverage[ind_INPN]     += INPN[iNt][iRing]     * area;
        ringAverage[ind_H]        += H[iNt][iRing]        * area;
        ringAverage[ind_BrNO3]    += BrNO3[iNt][iRing]    * area;
        ringAverage[ind_PRPE]     += PRPE[iNt][iRing]     * area;
        ringAverage[ind_MVKOO]    += MVKOO[iNt][iRing]    * area;
        ringAverage[ind_Cl2]      += Cl2[iNt][iRing]      * area;
        ringAverage[ind_ISOPND]   += ISOPND[iNt][iRing]   * area;
        ringAverage[ind_HOBr]     += HOBr[iNt][iRing]     * area;
        ringAverage[ind_A3O2]     += A3O2[iNt][iRing]     * area;
        ringAverage[ind_PROPNN]   += PROPNN[iNt][iRing]   * area;
        ringAverage[ind_GLYX]     += GLYX[iNt][iRing]     * area;
        ringAverage[ind_MAOPO2]   += MAOPO2[iNt][iRing]   * area;
        ringAverage[ind_CH4]      += CH4[iNt][iRing]      * area;
        ringAverage[ind_GAOO]     += GAOO[iNt][iRing]     * area;
        ringAverage[ind_B3O2]     += B3O2[iNt][iRing]     * area;
        ringAverage[ind_ACET]     += ACET[iNt][iRing]     * area;
        ringAverage[ind_MACRN]    += MACRN[iNt][iRing]    * area;
        ringAverage[ind_CH2OO]    += CH2OO[iNt][iRing]    * area;
        ringAverage[ind_MGLYOO]   += MGLYOO[iNt][iRing]   * area;
        ringAverage[ind_VRO2]     += VRO2[iNt][iRing]     * area;
        ringAverage[ind_MGLOO]    += MGLOO[iNt][iRing]    * area;
        ringAverage[ind_MACROO]   += MACROO[iNt][iRing]   * area;
        ringAverage[ind_PO2]      += PO2[iNt][iRing]      * area;
        ringAverage[ind_CH3CHOO]  += CH3CHOO[iNt][iRing]  * area;
        ringAverage[ind_MAN2]     += MAN2[iNt][iRing]     * area;
        ringAverage[ind_ISNOOA]   += ISNOOA[iNt][iRing]   * area;
        ringAverage[ind_H2O2]     += H2O2[iNt][iRing]     * area;
        ringAverage[ind_PRN1]     += PRN1[iNt][iRing]     * area;
        ringAverage[ind_ETO2]     += ETO2[iNt][iRing]     * area;
        ringAverage[ind_KO2]      += KO2[iNt][iRing]      * area;
        ringAverage[ind_RCO3]     += RCO3[iNt][iRing]     * area;
        ringAverage[ind_HC5OO]    += HC5OO[iNt][iRing]    * area;
        ringAverage[ind_GLYC]     += GLYC[iNt][iRing]     * area;
        ringAverage[ind_ClNO3]    += ClNO3[iNt][iRing]    * area;
        ringAverage[ind_RIO2]     += RIO2[iNt][iRing]     * area;
        ringAverage[ind_R4N1]     += R4N1[iNt][iRing]     * area;
        ringAverage[ind_HOCl]     += HOCl[iNt][iRing]     * area;
        ringAverage[ind_ATO2]     += ATO2[iNt][iRing]     * area;
        ringAverage[ind_HNO3]     += HNO3[iNt][iRing]     * area;
        ringAverage[ind_ISN1]     += ISN1[iNt][iRing]     * area;
        ringAverage[ind_MAO3]     += MAO3[iNt][iRing]     * area;
        ringAverage[ind_MRO2]     += MRO2[iNt][iRing]     * area;
        ringAverage[ind_INO2]     += INO2[iNt][iRing]     * area;
        ringAverage[ind_HAC]      += HAC[iNt][iRing]      * area;
        ringAverage[ind_HC5]      += HC5[iNt][iRing]      * area;
        ringAverage[ind_MGLY]     += MGLY[iNt][iRing]     * area;
        ringAverage[ind_ISOPNBO2] += ISOPNBO2[iNt][iRing] * area;
        ringAverage[ind_ISOPNDO2] += ISOPNDO2[iNt][iRing] * area;
        ringAverage[ind_R4O2]     += R4O2[iNt][iRing]     * area;
        ringAverage[ind_R4N2]     += R4N2[iNt][iRing]     * area;
        ringAverage[ind_BrO]      += BrO[iNt][iRing]      * area;
        ringAverage[ind_RCHO]     += RCHO[iNt][iRing]     * area;
        ringAverage[ind_MEK]      += MEK[iNt][iRing]      * area;
        ringAverage[ind_ClO]      += ClO[iNt][iRing]      * area;
        ringAverage[ind_MACR]     += MACR[iNt][iRing]     * area;
        ringAverage[ind_SO2]      += SO2[iNt][iRing]      * area;
        ringAverage[ind_MVK]      += MVK[iNt][iRing]      * area;
        ringAverage[ind_ALD2]     += ALD2[iNt][iRing]     * area;
        ringAverage[ind_MCO3]     += MCO3[iNt][iRing]     * area;
        ringAverage[ind_CH2O]     += CH2O[iNt][iRing]     * area;
        ringAverage[ind_H2O]      += H2O[iNt][iRing]      * area;
        ringAverage[ind_Br]       += Br[iNt][iRing]       * area;
        ringAverage[ind_NO]       += NO[iNt][iRing]       * area;
        ringAverage[ind_NO3]      += NO3[iNt][iRing]      * area;
        ringAverage[ind_Cl]       += Cl[iNt][iRing]       * area;
        ringAverage[ind_O]        += O[iNt][iRing]        * area;
        ringAverage[ind_O1D]      += O1D[iNt][iRing]      * area;
        ringAverage[ind_O3]       += O3[iNt][iRing]       * area;
        ringAverage[ind_HO2]      += HO2[iNt][iRing]      * area;
        ringAverage[ind_NO2]      += NO2[iNt][iRing]      * area;
        ringAverage[ind_OH]       += OH[iNt][iRing]       * area;
        ringAverage[ind_HBr]      += HBr[iNt][iRing]      * area;
        ringAverage[ind_HCl]      += HCl[iNt][iRing]      * area;
        ringAverage[ind_CO]       += CO[iNt][iRing]       * area;
        ringAverage[ind_MO2]      += MO2[iNt][iRing]      * area;
    }

    return ringAverage;

} /* End of SpeciesArray::RingAverage */

Vector_2D SpeciesArray::RingAverage( const Vector_1D ringArea, \
                                     const RealDouble totArea ) const
{

    Vector_2D ringAverage( nTime, Vector_1D( NVAR, 0.0E+00 ) );
    UInt iRing, iTime;
    RealDouble area;

    for ( iRing = 0; iRing < nRing; iRing++ ) {
        area = ringArea[iRing] / totArea;
        for ( iTime = 0; iTime < nTime; iTime++ ) {
            ringAverage[iTime][ind_CO2]      += CO2[iTime][iRing]      * area;
            ringAverage[iTime][ind_PPN]      += PPN[iTime][iRing]      * area;
            ringAverage[iTime][ind_BrNO2]    += BrNO2[iTime][iRing]    * area;
            ringAverage[iTime][ind_IEPOX]    += IEPOX[iTime][iRing]    * area;
            ringAverage[iTime][ind_PMNN]     += PMNN[iTime][iRing]     * area;
            ringAverage[iTime][ind_N2O]      += N2O[iTime][iRing]      * area;
            ringAverage[iTime][ind_N]        += N[iTime][iRing]        * area;
            ringAverage[iTime][ind_PAN]      += PAN[iTime][iRing]      * area;
            ringAverage[iTime][ind_ALK4]     += ALK4[iTime][iRing]     * area;
            ringAverage[iTime][ind_MAP]      += MAP[iTime][iRing]      * area;
            ringAverage[iTime][ind_MPN]      += MPN[iTime][iRing]      * area;
            ringAverage[iTime][ind_Cl2O2]    += Cl2O2[iTime][iRing]    * area;
            ringAverage[iTime][ind_ETP]      += ETP[iTime][iRing]      * area;
            ringAverage[iTime][ind_HNO2]     += HNO2[iTime][iRing]     * area;
            ringAverage[iTime][ind_C3H8]     += C3H8[iTime][iRing]     * area;
            ringAverage[iTime][ind_RA3P]     += RA3P[iTime][iRing]     * area;
            ringAverage[iTime][ind_RB3P]     += RB3P[iTime][iRing]     * area;
            ringAverage[iTime][ind_OClO]     += OClO[iTime][iRing]     * area;
            ringAverage[iTime][ind_ClNO2]    += ClNO2[iTime][iRing]    * area;
            ringAverage[iTime][ind_ISOP]     += ISOP[iTime][iRing]     * area;
            ringAverage[iTime][ind_HNO4]     += HNO4[iTime][iRing]     * area;
            ringAverage[iTime][ind_MAOP]     += MAOP[iTime][iRing]     * area;
            ringAverage[iTime][ind_MP]       += MP[iTime][iRing]       * area;
            ringAverage[iTime][ind_ClOO]     += ClOO[iTime][iRing]     * area;
            ringAverage[iTime][ind_RP]       += RP[iTime][iRing]       * area;
            ringAverage[iTime][ind_BrCl]     += BrCl[iTime][iRing]     * area;
            ringAverage[iTime][ind_PP]       += PP[iTime][iRing]       * area;
            ringAverage[iTime][ind_PRPN]     += PRPN[iTime][iRing]     * area;
            ringAverage[iTime][ind_SO4]      += SO4[iTime][iRing]      * area;
            ringAverage[iTime][ind_Br2]      += Br2[iTime][iRing]      * area;
            ringAverage[iTime][ind_ETHLN]    += ETHLN[iTime][iRing]    * area;
            ringAverage[iTime][ind_MVKN]     += MVKN[iTime][iRing]     * area;
            ringAverage[iTime][ind_R4P]      += R4P[iTime][iRing]      * area;
            ringAverage[iTime][ind_C2H6]     += C2H6[iTime][iRing]     * area;
            ringAverage[iTime][ind_RIP]      += RIP[iTime][iRing]      * area;
            ringAverage[iTime][ind_VRP]      += VRP[iTime][iRing]      * area;
            ringAverage[iTime][ind_ATOOH]    += ATOOH[iTime][iRing]    * area;
            ringAverage[iTime][ind_IAP]      += IAP[iTime][iRing]      * area;
            ringAverage[iTime][ind_DHMOB]    += DHMOB[iTime][iRing]    * area;
            ringAverage[iTime][ind_MOBA]     += MOBA[iTime][iRing]     * area;
            ringAverage[iTime][ind_MRP]      += MRP[iTime][iRing]      * area;
            ringAverage[iTime][ind_N2O5]     += N2O5[iTime][iRing]     * area;
            ringAverage[iTime][ind_ISNOHOO]  += ISNOHOO[iTime][iRing]  * area;
            ringAverage[iTime][ind_ISNP]     += ISNP[iTime][iRing]     * area;
            ringAverage[iTime][ind_ISOPNB]   += ISOPNB[iTime][iRing]   * area;
            ringAverage[iTime][ind_IEPOXOO]  += IEPOXOO[iTime][iRing]  * area;
            ringAverage[iTime][ind_MACRNO2]  += MACRNO2[iTime][iRing]  * area;
            ringAverage[iTime][ind_ROH]      += ROH[iTime][iRing]      * area;
            ringAverage[iTime][ind_MOBAOO]   += MOBAOO[iTime][iRing]   * area;
            ringAverage[iTime][ind_DIBOO]    += DIBOO[iTime][iRing]    * area;
            ringAverage[iTime][ind_PMN]      += PMN[iTime][iRing]      * area;
            ringAverage[iTime][ind_ISNOOB]   += ISNOOB[iTime][iRing]   * area;
            ringAverage[iTime][ind_INPN]     += INPN[iTime][iRing]     * area;
            ringAverage[iTime][ind_H]        += H[iTime][iRing]        * area;
            ringAverage[iTime][ind_BrNO3]    += BrNO3[iTime][iRing]    * area;
            ringAverage[iTime][ind_PRPE]     += PRPE[iTime][iRing]     * area;
            ringAverage[iTime][ind_MVKOO]    += MVKOO[iTime][iRing]    * area;
            ringAverage[iTime][ind_Cl2]      += Cl2[iTime][iRing]      * area;
            ringAverage[iTime][ind_ISOPND]   += ISOPND[iTime][iRing]   * area;
            ringAverage[iTime][ind_HOBr]     += HOBr[iTime][iRing]     * area;
            ringAverage[iTime][ind_A3O2]     += A3O2[iTime][iRing]     * area;
            ringAverage[iTime][ind_PROPNN]   += PROPNN[iTime][iRing]   * area;
            ringAverage[iTime][ind_GLYX]     += GLYX[iTime][iRing]     * area;
            ringAverage[iTime][ind_MAOPO2]   += MAOPO2[iTime][iRing]   * area;
            ringAverage[iTime][ind_CH4]      += CH4[iTime][iRing]      * area;
            ringAverage[iTime][ind_GAOO]     += GAOO[iTime][iRing]     * area;
            ringAverage[iTime][ind_B3O2]     += B3O2[iTime][iRing]     * area;
            ringAverage[iTime][ind_ACET]     += ACET[iTime][iRing]     * area;
            ringAverage[iTime][ind_MACRN]    += MACRN[iTime][iRing]    * area;
            ringAverage[iTime][ind_CH2OO]    += CH2OO[iTime][iRing]    * area;
            ringAverage[iTime][ind_MGLYOO]   += MGLYOO[iTime][iRing]   * area;
            ringAverage[iTime][ind_VRO2]     += VRO2[iTime][iRing]     * area;
            ringAverage[iTime][ind_MGLOO]    += MGLOO[iTime][iRing]    * area;
            ringAverage[iTime][ind_MACROO]   += MACROO[iTime][iRing]   * area;
            ringAverage[iTime][ind_PO2]      += PO2[iTime][iRing]      * area;
            ringAverage[iTime][ind_CH3CHOO]  += CH3CHOO[iTime][iRing]  * area;
            ringAverage[iTime][ind_MAN2]     += MAN2[iTime][iRing]     * area;
            ringAverage[iTime][ind_ISNOOA]   += ISNOOA[iTime][iRing]   * area;
            ringAverage[iTime][ind_H2O2]     += H2O2[iTime][iRing]     * area;
            ringAverage[iTime][ind_PRN1]     += PRN1[iTime][iRing]     * area;
            ringAverage[iTime][ind_ETO2]     += ETO2[iTime][iRing]     * area;
            ringAverage[iTime][ind_KO2]      += KO2[iTime][iRing]      * area;
            ringAverage[iTime][ind_RCO3]     += RCO3[iTime][iRing]     * area;
            ringAverage[iTime][ind_HC5OO]    += HC5OO[iTime][iRing]    * area;
            ringAverage[iTime][ind_GLYC]     += GLYC[iTime][iRing]     * area;
            ringAverage[iTime][ind_ClNO3]    += ClNO3[iTime][iRing]    * area;
            ringAverage[iTime][ind_RIO2]     += RIO2[iTime][iRing]     * area;
            ringAverage[iTime][ind_R4N1]     += R4N1[iTime][iRing]     * area;
            ringAverage[iTime][ind_HOCl]     += HOCl[iTime][iRing]     * area;
            ringAverage[iTime][ind_ATO2]     += ATO2[iTime][iRing]     * area;
            ringAverage[iTime][ind_HNO3]     += HNO3[iTime][iRing]     * area;
            ringAverage[iTime][ind_ISN1]     += ISN1[iTime][iRing]     * area;
            ringAverage[iTime][ind_MAO3]     += MAO3[iTime][iRing]     * area;
            ringAverage[iTime][ind_MRO2]     += MRO2[iTime][iRing]     * area;
            ringAverage[iTime][ind_INO2]     += INO2[iTime][iRing]     * area;
            ringAverage[iTime][ind_HAC]      += HAC[iTime][iRing]      * area;
            ringAverage[iTime][ind_HC5]      += HC5[iTime][iRing]      * area;
            ringAverage[iTime][ind_MGLY]     += MGLY[iTime][iRing]     * area;
            ringAverage[iTime][ind_ISOPNBO2] += ISOPNBO2[iTime][iRing] * area;
            ringAverage[iTime][ind_ISOPNDO2] += ISOPNDO2[iTime][iRing] * area;
            ringAverage[iTime][ind_R4O2]     += R4O2[iTime][iRing]     * area;
            ringAverage[iTime][ind_R4N2]     += R4N2[iTime][iRing]     * area;
            ringAverage[iTime][ind_BrO]      += BrO[iTime][iRing]      * area;
            ringAverage[iTime][ind_RCHO]     += RCHO[iTime][iRing]     * area;
            ringAverage[iTime][ind_MEK]      += MEK[iTime][iRing]      * area;
            ringAverage[iTime][ind_ClO]      += ClO[iTime][iRing]      * area;
            ringAverage[iTime][ind_MACR]     += MACR[iTime][iRing]     * area;
            ringAverage[iTime][ind_SO2]      += SO2[iTime][iRing]      * area;
            ringAverage[iTime][ind_MVK]      += MVK[iTime][iRing]      * area;
            ringAverage[iTime][ind_ALD2]     += ALD2[iTime][iRing]     * area;
            ringAverage[iTime][ind_MCO3]     += MCO3[iTime][iRing]     * area;
            ringAverage[iTime][ind_CH2O]     += CH2O[iTime][iRing]     * area;
            ringAverage[iTime][ind_H2O]      += H2O[iTime][iRing]      * area;
            ringAverage[iTime][ind_Br]       += Br[iTime][iRing]       * area;
            ringAverage[iTime][ind_NO]       += NO[iTime][iRing]       * area;
            ringAverage[iTime][ind_NO3]      += NO3[iTime][iRing]      * area;
            ringAverage[iTime][ind_Cl]       += Cl[iTime][iRing]       * area;
            ringAverage[iTime][ind_O]        += O[iTime][iRing]        * area;
            ringAverage[iTime][ind_O1D]      += O1D[iTime][iRing]      * area;
            ringAverage[iTime][ind_O3]       += O3[iTime][iRing]       * area;
            ringAverage[iTime][ind_HO2]      += HO2[iTime][iRing]      * area;
            ringAverage[iTime][ind_NO2]      += NO2[iTime][iRing]      * area;
            ringAverage[iTime][ind_OH]       += OH[iTime][iRing]       * area;
            ringAverage[iTime][ind_HBr]      += HBr[iTime][iRing]      * area;
            ringAverage[iTime][ind_HCl]      += HCl[iTime][iRing]      * area;
            ringAverage[iTime][ind_CO]       += CO[iTime][iRing]       * area;
            ringAverage[iTime][ind_MO2]      += MO2[iTime][iRing]      * area;
        }
    }

    return ringAverage;

} /* End of SpeciesArray::RingAverage */

/* End of Species.cpp */
