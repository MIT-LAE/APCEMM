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

#include "Species.hpp"

static const double ZERO = 1.00E-50;

SpeciesArray::SpeciesArray( )
{
    /* Default Constructor */

} /* End of SpeciesArray::SpeciesArray */

SpeciesArray::SpeciesArray( unsigned int nRing_, unsigned int nTime_, bool halfRing_ )
{

    /* Constructor */

    nRing = nRing_;
    nTime = nTime_;
    halfRing = halfRing_;

    std::vector<double> v1d = std::vector<double>( nRing, 0.0 );

    for ( unsigned int i = 0; i < nTime; i++ ) {
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

    for ( unsigned int iRing = 0; iRing < nRing; iRing++ ) {
        for ( unsigned int iTime = 0; iTime < nTime; iTime++ ) {

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

    for ( unsigned int iRing = 0; iRing < nRing; iRing++ ) {
        for ( unsigned int iTime = 0; iTime < nTime; iTime++ ) {

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
            ACTA     -= sp.ACTA;
            EOH      -= sp.EOH;
            H2       -= sp.H2;
            HCOOH    -= sp.HCOOH;
            MOH      -= sp.MOH;
            N2       -= sp.N2;
            O2       -= sp.O2;
            RCOOH    -= sp.RCOOH;

        }
    }
    return *this;

} /* End of SpeciesArray::operator- */

SpeciesArray::~SpeciesArray( )
{
    /* Destructor */

} /* End of SpeciesArray::~SpeciesArray */

void SpeciesArray::FillIn( Solution &Data, Mesh &m, unsigned int nCounter )
{

    std::vector<std::vector<std::vector<bool> > > map = m.getMap();
    std::vector<unsigned int> nMap = m.getnMap();
    std::vector<std::vector<std::pair<unsigned int, unsigned int>>> indList = m.getList();

    unsigned int i, j;
    unsigned int nCell;
    for ( unsigned int iRing = 0; iRing < nRing; iRing++ ) {
        for ( unsigned int iList = 0; iList < indList[iRing].size(); iList++ ) {
            i = indList[iRing][iList].first;
            j = indList[iRing][iList].second;
        
            /* Sum concentrations over rings */
                
            CO2[nCounter][iRing]      += Data.CO2[j][i];
            PPN[nCounter][iRing]      += Data.PPN[j][i];
            BrNO2[nCounter][iRing]    += Data.BrNO2[j][i];
            IEPOX[nCounter][iRing]    += Data.IEPOX[j][i];
            PMNN[nCounter][iRing]     += Data.PMNN[j][i];
            N2O[nCounter][iRing]      += Data.N2O[j][i];
            N[nCounter][iRing]        += Data.N[j][i];
            PAN[nCounter][iRing]      += Data.PAN[j][i];
            ALK4[nCounter][iRing]     += Data.ALK4[j][i];
            MAP[nCounter][iRing]      += Data.MAP[j][i];
            MPN[nCounter][iRing]      += Data.MPN[j][i];
            Cl2O2[nCounter][iRing]    += Data.Cl2O2[j][i];
            ETP[nCounter][iRing]      += Data.ETP[j][i];
            HNO2[nCounter][iRing]     += Data.HNO2[j][i];
            C3H8[nCounter][iRing]     += Data.C3H8[j][i];
            RA3P[nCounter][iRing]     += Data.RA3P[j][i];
            RB3P[nCounter][iRing]     += Data.RB3P[j][i];
            OClO[nCounter][iRing]     += Data.OClO[j][i];
            ClNO2[nCounter][iRing]    += Data.ClNO2[j][i];
            ISOP[nCounter][iRing]     += Data.ISOP[j][i];
            HNO4[nCounter][iRing]     += Data.HNO4[j][i];
            MAOP[nCounter][iRing]     += Data.MAOP[j][i];
            MP[nCounter][iRing]       += Data.MP[j][i];
            ClOO[nCounter][iRing]     += Data.ClOO[j][i];
            RP[nCounter][iRing]       += Data.RP[j][i];
            BrCl[nCounter][iRing]     += Data.BrCl[j][i];
            PP[nCounter][iRing]       += Data.PP[j][i];
            PRPN[nCounter][iRing]     += Data.PRPN[j][i];
            SO4[nCounter][iRing]      += Data.SO4[j][i];
            Br2[nCounter][iRing]      += Data.Br2[j][i];
            ETHLN[nCounter][iRing]    += Data.ETHLN[j][i];
            MVKN[nCounter][iRing]     += Data.MVKN[j][i];
            R4P[nCounter][iRing]      += Data.R4P[j][i];
            C2H6[nCounter][iRing]     += Data.C2H6[j][i];
            RIP[nCounter][iRing]      += Data.RIP[j][i];
            VRP[nCounter][iRing]      += Data.VRP[j][i];
            ATOOH[nCounter][iRing]    += Data.ATOOH[j][i];
            IAP[nCounter][iRing]      += Data.IAP[j][i];
            DHMOB[nCounter][iRing]    += Data.DHMOB[j][i];
            MOBA[nCounter][iRing]     += Data.MOBA[j][i];
            MRP[nCounter][iRing]      += Data.MRP[j][i];
            N2O5[nCounter][iRing]     += Data.N2O5[j][i];
            ISNOHOO[nCounter][iRing]  += Data.ISNOHOO[j][i];
            ISNP[nCounter][iRing]     += Data.ISNP[j][i];
            ISOPNB[nCounter][iRing]   += Data.ISOPNB[j][i];
            IEPOXOO[nCounter][iRing]  += Data.IEPOXOO[j][i];
            MACRNO2[nCounter][iRing]  += Data.MACRNO2[j][i];
            ROH[nCounter][iRing]      += Data.ROH[j][i];
            MOBAOO[nCounter][iRing]   += Data.MOBAOO[j][i];
            DIBOO[nCounter][iRing]    += Data.DIBOO[j][i];
            PMN[nCounter][iRing]      += Data.PMN[j][i];
            ISNOOB[nCounter][iRing]   += Data.ISNOOB[j][i];
            INPN[nCounter][iRing]     += Data.INPN[j][i];
            H[nCounter][iRing]        += Data.H[j][i];
            BrNO3[nCounter][iRing]    += Data.BrNO3[j][i];
            PRPE[nCounter][iRing]     += Data.PRPE[j][i];
            MVKOO[nCounter][iRing]    += Data.MVKOO[j][i];
            Cl2[nCounter][iRing]      += Data.Cl2[j][i];
            ISOPND[nCounter][iRing]   += Data.ISOPND[j][i];
            HOBr[nCounter][iRing]     += Data.HOBr[j][i];
            A3O2[nCounter][iRing]     += Data.A3O2[j][i];
            PROPNN[nCounter][iRing]   += Data.PROPNN[j][i];
            GLYX[nCounter][iRing]     += Data.GLYX[j][i];
            MAOPO2[nCounter][iRing]   += Data.MAOPO2[j][i];
            CH4[nCounter][iRing]      += Data.CH4[j][i];
            GAOO[nCounter][iRing]     += Data.GAOO[j][i];
            B3O2[nCounter][iRing]     += Data.B3O2[j][i];
            ACET[nCounter][iRing]     += Data.ACET[j][i];
            MACRN[nCounter][iRing]    += Data.MACRN[j][i];
            CH2OO[nCounter][iRing]    += Data.CH2OO[j][i];
            MGLYOO[nCounter][iRing]   += Data.MGLYOO[j][i];
            VRO2[nCounter][iRing]     += Data.VRO2[j][i];
            MGLOO[nCounter][iRing]    += Data.MGLOO[j][i];
            MACROO[nCounter][iRing]   += Data.MACROO[j][i];
            PO2[nCounter][iRing]      += Data.PO2[j][i];
            CH3CHOO[nCounter][iRing]  += Data.CH3CHOO[j][i];
            MAN2[nCounter][iRing]     += Data.MAN2[j][i];
            ISNOOA[nCounter][iRing]   += Data.ISNOOA[j][i];
            H2O2[nCounter][iRing]     += Data.H2O2[j][i];
            PRN1[nCounter][iRing]     += Data.PRN1[j][i];
            ETO2[nCounter][iRing]     += Data.ETO2[j][i];
            KO2[nCounter][iRing]      += Data.KO2[j][i];
            RCO3[nCounter][iRing]     += Data.RCO3[j][i];
            HC5OO[nCounter][iRing]    += Data.HC5OO[j][i];
            GLYC[nCounter][iRing]     += Data.GLYC[j][i];
            ClNO3[nCounter][iRing]    += Data.ClNO3[j][i];
            RIO2[nCounter][iRing]     += Data.RIO2[j][i];
            R4N1[nCounter][iRing]     += Data.R4N1[j][i];
            HOCl[nCounter][iRing]     += Data.HOCl[j][i];
            ATO2[nCounter][iRing]     += Data.ATO2[j][i];
            HNO3[nCounter][iRing]     += Data.HNO3[j][i];
            ISN1[nCounter][iRing]     += Data.ISN1[j][i];
            MAO3[nCounter][iRing]     += Data.MAO3[j][i];
            MRO2[nCounter][iRing]     += Data.MRO2[j][i];
            INO2[nCounter][iRing]     += Data.INO2[j][i];
            HAC[nCounter][iRing]      += Data.HAC[j][i];
            HC5[nCounter][iRing]      += Data.HC5[j][i];
            MGLY[nCounter][iRing]     += Data.MGLY[j][i];
            ISOPNBO2[nCounter][iRing] += Data.ISOPNBO2[j][i];
            ISOPNDO2[nCounter][iRing] += Data.ISOPNDO2[j][i];
            R4O2[nCounter][iRing]     += Data.R4O2[j][i];
            R4N2[nCounter][iRing]     += Data.R4N2[j][i];
            BrO[nCounter][iRing]      += Data.BrO[j][i];
            RCHO[nCounter][iRing]     += Data.RCHO[j][i];
            MEK[nCounter][iRing]      += Data.MEK[j][i];
            ClO[nCounter][iRing]      += Data.ClO[j][i];
            MACR[nCounter][iRing]     += Data.MACR[j][i];
            SO2[nCounter][iRing]      += Data.SO2[j][i];
            MVK[nCounter][iRing]      += Data.MVK[j][i];
            ALD2[nCounter][iRing]     += Data.ALD2[j][i];
            MCO3[nCounter][iRing]     += Data.MCO3[j][i];
            CH2O[nCounter][iRing]     += Data.CH2O[j][i];
            H2O[nCounter][iRing]      += Data.H2O[j][i];
            Br[nCounter][iRing]       += Data.Br[j][i];
            NO[nCounter][iRing]       += Data.NO[j][i];
            NO3[nCounter][iRing]      += Data.NO3[j][i];
            Cl[nCounter][iRing]       += Data.Cl[j][i];
            O[nCounter][iRing]        += Data.O[j][i];
            O1D[nCounter][iRing]      += Data.O1D[j][i];
            O3[nCounter][iRing]       += Data.O3[j][i];
            HO2[nCounter][iRing]      += Data.HO2[j][i];
            NO2[nCounter][iRing]      += Data.NO2[j][i];
            OH[nCounter][iRing]       += Data.OH[j][i];
            HBr[nCounter][iRing]      += Data.HBr[j][i];
            HCl[nCounter][iRing]      += Data.HCl[j][i];
            CO[nCounter][iRing]       += Data.CO[j][i];
            MO2[nCounter][iRing]      += Data.MO2[j][i];
        }
        
        nCell = indList[iRing].size();
    
        /* Divide by number of cells in each ring */
        CO2[nCounter][iRing]      /= nCell;
        PPN[nCounter][iRing]      /= nCell;
        BrNO2[nCounter][iRing]    /= nCell;
        IEPOX[nCounter][iRing]    /= nCell;
        PMNN[nCounter][iRing]     /= nCell;
        N2O[nCounter][iRing]      /= nCell;
        N[nCounter][iRing]        /= nCell;
        PAN[nCounter][iRing]      /= nCell;
        ALK4[nCounter][iRing]     /= nCell;
        MAP[nCounter][iRing]      /= nCell;
        MPN[nCounter][iRing]      /= nCell;
        Cl2O2[nCounter][iRing]    /= nCell;
        ETP[nCounter][iRing]      /= nCell;
        HNO2[nCounter][iRing]     /= nCell;
        C3H8[nCounter][iRing]     /= nCell;
        RA3P[nCounter][iRing]     /= nCell;
        RB3P[nCounter][iRing]     /= nCell;
        OClO[nCounter][iRing]     /= nCell;
        ClNO2[nCounter][iRing]    /= nCell;
        ISOP[nCounter][iRing]     /= nCell;
        HNO4[nCounter][iRing]     /= nCell;
        MAOP[nCounter][iRing]     /= nCell;
        MP[nCounter][iRing]       /= nCell;
        ClOO[nCounter][iRing]     /= nCell;
        RP[nCounter][iRing]       /= nCell;
        BrCl[nCounter][iRing]     /= nCell;
        PP[nCounter][iRing]       /= nCell;
        PRPN[nCounter][iRing]     /= nCell;
        SO4[nCounter][iRing]      /= nCell;
        Br2[nCounter][iRing]      /= nCell;
        ETHLN[nCounter][iRing]    /= nCell;
        MVKN[nCounter][iRing]     /= nCell;
        R4P[nCounter][iRing]      /= nCell;
        C2H6[nCounter][iRing]     /= nCell;
        RIP[nCounter][iRing]      /= nCell;
        VRP[nCounter][iRing]      /= nCell;
        ATOOH[nCounter][iRing]    /= nCell;
        IAP[nCounter][iRing]      /= nCell;
        DHMOB[nCounter][iRing]    /= nCell;
        MOBA[nCounter][iRing]     /= nCell;
        MRP[nCounter][iRing]      /= nCell;
        N2O5[nCounter][iRing]     /= nCell;
        ISNOHOO[nCounter][iRing]  /= nCell;
        ISNP[nCounter][iRing]     /= nCell;
        ISOPNB[nCounter][iRing]   /= nCell;
        IEPOXOO[nCounter][iRing]  /= nCell;
        MACRNO2[nCounter][iRing]  /= nCell;
        ROH[nCounter][iRing]      /= nCell;
        MOBAOO[nCounter][iRing]   /= nCell;
        DIBOO[nCounter][iRing]    /= nCell;
        PMN[nCounter][iRing]      /= nCell;
        ISNOOB[nCounter][iRing]   /= nCell;
        INPN[nCounter][iRing]     /= nCell;
        H[nCounter][iRing]        /= nCell;
        BrNO3[nCounter][iRing]    /= nCell;
        PRPE[nCounter][iRing]     /= nCell;
        MVKOO[nCounter][iRing]    /= nCell;
        Cl2[nCounter][iRing]      /= nCell;
        ISOPND[nCounter][iRing]   /= nCell;
        HOBr[nCounter][iRing]     /= nCell;
        A3O2[nCounter][iRing]     /= nCell;
        PROPNN[nCounter][iRing]   /= nCell;
        GLYX[nCounter][iRing]     /= nCell;
        MAOPO2[nCounter][iRing]   /= nCell;
        CH4[nCounter][iRing]      /= nCell;
        GAOO[nCounter][iRing]     /= nCell;
        B3O2[nCounter][iRing]     /= nCell;
        ACET[nCounter][iRing]     /= nCell;
        MACRN[nCounter][iRing]    /= nCell;
        CH2OO[nCounter][iRing]    /= nCell;
        MGLYOO[nCounter][iRing]   /= nCell;
        VRO2[nCounter][iRing]     /= nCell;
        MGLOO[nCounter][iRing]    /= nCell;
        MACROO[nCounter][iRing]   /= nCell;
        PO2[nCounter][iRing]      /= nCell;
        CH3CHOO[nCounter][iRing]  /= nCell;
        MAN2[nCounter][iRing]     /= nCell;
        ISNOOA[nCounter][iRing]   /= nCell;
        H2O2[nCounter][iRing]     /= nCell;
        PRN1[nCounter][iRing]     /= nCell;
        ETO2[nCounter][iRing]     /= nCell;
        KO2[nCounter][iRing]      /= nCell;
        RCO3[nCounter][iRing]     /= nCell;
        HC5OO[nCounter][iRing]    /= nCell;
        GLYC[nCounter][iRing]     /= nCell;
        ClNO3[nCounter][iRing]    /= nCell;
        RIO2[nCounter][iRing]     /= nCell;
        R4N1[nCounter][iRing]     /= nCell;
        HOCl[nCounter][iRing]     /= nCell;
        ATO2[nCounter][iRing]     /= nCell;
        HNO3[nCounter][iRing]     /= nCell;
        ISN1[nCounter][iRing]     /= nCell;
        MAO3[nCounter][iRing]     /= nCell;
        MRO2[nCounter][iRing]     /= nCell;
        INO2[nCounter][iRing]     /= nCell;
        HAC[nCounter][iRing]      /= nCell;
        HC5[nCounter][iRing]      /= nCell;
        MGLY[nCounter][iRing]     /= nCell;
        ISOPNBO2[nCounter][iRing] /= nCell;
        ISOPNDO2[nCounter][iRing] /= nCell;
        R4O2[nCounter][iRing]     /= nCell;
        R4N2[nCounter][iRing]     /= nCell;
        BrO[nCounter][iRing]      /= nCell;
        RCHO[nCounter][iRing]     /= nCell;
        MEK[nCounter][iRing]      /= nCell;
        ClO[nCounter][iRing]      /= nCell;
        MACR[nCounter][iRing]     /= nCell;
        SO2[nCounter][iRing]      /= nCell;
        MVK[nCounter][iRing]      /= nCell;
        ALD2[nCounter][iRing]     /= nCell;
        MCO3[nCounter][iRing]     /= nCell;
        CH2O[nCounter][iRing]     /= nCell;
        H2O[nCounter][iRing]      /= nCell;
        Br[nCounter][iRing]       /= nCell;
        NO[nCounter][iRing]       /= nCell;
        NO3[nCounter][iRing]      /= nCell;
        Cl[nCounter][iRing]       /= nCell;
        O[nCounter][iRing]        /= nCell;
        O1D[nCounter][iRing]      /= nCell;
        O3[nCounter][iRing]       /= nCell;
        HO2[nCounter][iRing]      /= nCell;
        NO2[nCounter][iRing]      /= nCell;
        OH[nCounter][iRing]       /= nCell;
        HBr[nCounter][iRing]      /= nCell;
        HCl[nCounter][iRing]      /= nCell;
        CO[nCounter][iRing]       /= nCell;
        MO2[nCounter][iRing]      /= nCell;

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

void SpeciesArray::FillIn( double varArray[], unsigned int iTime, unsigned int iRing )
{

    /* Ensure positiveness */
    for ( unsigned int i = 0; i < N_VAR; i++ ) {
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

void SpeciesArray::getData( double varArray[], double fixArray[], unsigned int iTime, unsigned int iRing )
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
    for ( unsigned int i = 0; i < N_VAR; i++ ) {
        if ( varArray[i] <= 0.0 ) {
            varArray[i] = ZERO;
        }
    }

} /* End of SpeciesArray::getData */

unsigned int SpeciesArray::getnRing() const
{

    return nRing;

} /* End of SpeciesArray::getnRing */

unsigned int SpeciesArray::getnTime() const
{

    return nTime;

} /* End of SpeciesArray::getnTime */

bool SpeciesArray::gethalfRing() const
{

    return halfRing;

} /* End of SpeciesArray::gethalfRing */

/* End of Species.cpp */
