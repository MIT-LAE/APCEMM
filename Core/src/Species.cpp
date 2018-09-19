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

SpeciesArray::SpeciesArray( unsigned int nRing_, unsigned int nTime_ )
{

    /* Constructor */

    nRing = nRing_;
    nTime = nTime_;

    std::vector<double> v1d;
    v1d = std::vector<double>( nRing + 1, 0.0 );

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
                
            CO2[nCounter][iRing]   += Data.CO2[j][i];
            PPN[nCounter][iRing]   += Data.PPN[j][i];
            BrNO2[nCounter][iRing] += Data.BrNO2[j][i];
            IEPOX[nCounter][iRing] += Data.IEPOX[j][i];
            PMNN[nCounter][iRing]  += Data.PMNN[j][i];
            N2O[nCounter][iRing]   += Data.N2O[j][i];
            N[nCounter][iRing]     += Data.N[j][i];
            PAN[nCounter][iRing]   += Data.PAN[j][i];
            ALK4[nCounter][iRing]  += Data.ALK4[j][i];
            MAP[nCounter][iRing]   += Data.MAP[j][i];
            MPN[nCounter][iRing]   += Data.MPN[j][i];
            Cl2O2[nCounter][iRing] += Data.Cl2O2[j][i];
            ETP[nCounter][iRing]   += Data.ETP[j][i];
            HNO2[nCounter][iRing]  += Data.HNO2[j][i];
            C3H8[nCounter][iRing]  += Data.C3H8[j][i];
            RA3P[nCounter][iRing]  += Data.RA3P[j][i];
            RB3P[nCounter][iRing]  += Data.RB3P[j][i];
            OClO[nCounter][iRing]  += Data.OClO[j][i];
            ClNO2[nCounter][iRing] += Data.ClNO2[j][i];
            ISOP[nCounter][iRing]  += Data.ISOP[j][i];
            HNO4[nCounter][iRing]  += Data.HNO4[j][i];
            MAOP[nCounter][iRing]  += Data.MAOP[j][i];
            MP[nCounter][iRing]    += Data.MP[j][i];
            ClOO[nCounter][iRing]  += Data.ClOO[j][i];
            RP[nCounter][iRing]    += Data.RP[j][i];
            BrCl[nCounter][iRing]  += Data.BrCl[j][i];
            PP[nCounter][iRing]    += Data.PP[j][i];
            PRPN[nCounter][iRing]  += Data.PRPN[j][i];
            SO4[nCounter][iRing]   += Data.SO4[j][i];
            Br2[nCounter][iRing]   += Data.Br2[j][i];
            ETHLN[nCounter][iRing] += Data.ETHLN[j][i];
            MVKN[nCounter][iRing]  += Data.MVKN[j][i];
            R4P[nCounter][iRing]   += Data.R4P[j][i];
            C2H6[nCounter][iRing]  += Data.C2H6[j][i];
            RIP[nCounter][iRing]   += Data.RIP[j][i];
            VRP[nCounter][iRing]   += Data.VRP[j][i];
            ATOOH[nCounter][iRing] += Data.ATOOH[j][i];
            IAP[nCounter][iRing]   += Data.IAP[j][i];
            DHMOB[nCounter][iRing] += Data.DHMOB[j][i];
            MOBA[nCounter][iRing]  += Data.MOBA[j][i];
            MRP[nCounter][iRing]   += Data.MRP[j][i];
            N2O5[nCounter][iRing]  += Data.N2O5[j][i];
            ISNOHOO[nCounter][iRing]+= Data.ISNOHOO[j][i];
            ISNP[nCounter][iRing]  += Data.ISNP[j][i];
            ISOPNB[nCounter][iRing]+= Data.ISOPNB[j][i];
            IEPOXOO[nCounter][iRing]+= Data.IEPOXOO[j][i];
            MACRNO2[nCounter][iRing]+= Data.MACRNO2[j][i];
            ROH[nCounter][iRing]   += Data.ROH[j][i];
            MOBAOO[nCounter][iRing]+= Data.MOBAOO[j][i];
            DIBOO[nCounter][iRing] += Data.DIBOO[j][i];
            PMN[nCounter][iRing]   += Data.PMN[j][i];
            ISNOOB[nCounter][iRing]+= Data.ISNOOB[j][i];
            INPN[nCounter][iRing]  += Data.INPN[j][i];
            H[nCounter][iRing]     += Data.H[j][i];
            BrNO3[nCounter][iRing] += Data.BrNO3[j][i];
            PRPE[nCounter][iRing]  += Data.PRPE[j][i];
            MVKOO[nCounter][iRing] += Data.MVKOO[j][i];
            Cl2[nCounter][iRing]   += Data.Cl2[j][i];
            ISOPND[nCounter][iRing]+= Data.ISOPND[j][i];
            HOBr[nCounter][iRing]  += Data.HOBr[j][i];
            A3O2[nCounter][iRing]  += Data.A3O2[j][i];
            PROPNN[nCounter][iRing]+= Data.PROPNN[j][i];
            GLYX[nCounter][iRing]  += Data.GLYX[j][i];
            MAOPO2[nCounter][iRing]+= Data.MAOPO2[j][i];
            CH4[nCounter][iRing]   += Data.CH4[j][i];
            GAOO[nCounter][iRing]  += Data.GAOO[j][i];
            B3O2[nCounter][iRing]  += Data.B3O2[j][i];
            ACET[nCounter][iRing]  += Data.ACET[j][i];
            MACRN[nCounter][iRing] += Data.MACRN[j][i];
            CH2OO[nCounter][iRing] += Data.CH2OO[j][i];
            MGLYOO[nCounter][iRing]+= Data.MGLYOO[j][i];
            VRO2[nCounter][iRing]  += Data.VRO2[j][i];
            MGLOO[nCounter][iRing] += Data.MGLOO[j][i];
            MACROO[nCounter][iRing]+= Data.MACROO[j][i];
            PO2[nCounter][iRing]   += Data.PO2[j][i];
            CH3CHOO[nCounter][iRing]+= Data.CH3CHOO[j][i];
            MAN2[nCounter][iRing]  += Data.MAN2[j][i];
            ISNOOA[nCounter][iRing]+= Data.ISNOOA[j][i];
            H2O2[nCounter][iRing]  += Data.H2O2[j][i];
            PRN1[nCounter][iRing]  += Data.PRN1[j][i];
            ETO2[nCounter][iRing]  += Data.ETO2[j][i];
            KO2[nCounter][iRing]   += Data.KO2[j][i];
            RCO3[nCounter][iRing]  += Data.RCO3[j][i];
            HC5OO[nCounter][iRing] += Data.HC5OO[j][i];
            GLYC[nCounter][iRing]  += Data.GLYC[j][i];
            ClNO3[nCounter][iRing] += Data.ClNO3[j][i];
            RIO2[nCounter][iRing]  += Data.RIO2[j][i];
            R4N1[nCounter][iRing]  += Data.R4N1[j][i];
            HOCl[nCounter][iRing]  += Data.HOCl[j][i];
            ATO2[nCounter][iRing]  += Data.ATO2[j][i];
            HNO3[nCounter][iRing]  += Data.HNO3[j][i];
            ISN1[nCounter][iRing]  += Data.ISN1[j][i];
            MAO3[nCounter][iRing]  += Data.MAO3[j][i];
            MRO2[nCounter][iRing]  += Data.MRO2[j][i];
            INO2[nCounter][iRing]  += Data.INO2[j][i];
            HAC[nCounter][iRing]   += Data.HAC[j][i];
            HC5[nCounter][iRing]   += Data.HC5[j][i];
            MGLY[nCounter][iRing]  += Data.MGLY[j][i];
            ISOPNBO2[nCounter][iRing]+= Data.ISOPNBO2[j][i];
            ISOPNDO2[nCounter][iRing]+= Data.ISOPNDO2[j][i];
            R4O2[nCounter][iRing]  += Data.R4O2[j][i];
            R4N2[nCounter][iRing]  += Data.R4N2[j][i];
            BrO[nCounter][iRing]   += Data.BrO[j][i];
            RCHO[nCounter][iRing]  += Data.RCHO[j][i];
            MEK[nCounter][iRing]   += Data.MEK[j][i];
            ClO[nCounter][iRing]   += Data.ClO[j][i];
            MACR[nCounter][iRing]  += Data.MACR[j][i];
            SO2[nCounter][iRing]   += Data.SO2[j][i];
            MVK[nCounter][iRing]   += Data.MVK[j][i];
            ALD2[nCounter][iRing]  += Data.ALD2[j][i];
            MCO3[nCounter][iRing]  += Data.MCO3[j][i];
            CH2O[nCounter][iRing]  += Data.CH2O[j][i];
            H2O[nCounter][iRing]   += Data.H2O[j][i];
            Br[nCounter][iRing]    += Data.Br[j][i];
            NO[nCounter][iRing]    += Data.NO[j][i];
            NO3[nCounter][iRing]   += Data.NO3[j][i];
            Cl[nCounter][iRing]    += Data.Cl[j][i];
            O[nCounter][iRing]     += Data.O[j][i];
            O1D[nCounter][iRing]   += Data.O1D[j][i];
            O3[nCounter][iRing]    += Data.O3[j][i];
            HO2[nCounter][iRing]   += Data.HO2[j][i];
            NO2[nCounter][iRing]   += Data.NO2[j][i];
            OH[nCounter][iRing]    += Data.OH[j][i];
            HBr[nCounter][iRing]   += Data.HBr[j][i];
            HCl[nCounter][iRing]   += Data.HCl[j][i];
            CO[nCounter][iRing]    += Data.CO[j][i];
            MO2[nCounter][iRing]   += Data.MO2[j][i];
        }
        
        nCell = indList[iRing].size();
    
        /* Divide by number of cells in each ring */
        CO2[nCounter][iRing]   /= nCell;
        PPN[nCounter][iRing]   /= nCell;
        BrNO2[nCounter][iRing] /= nCell;
        IEPOX[nCounter][iRing] /= nCell;
        PMNN[nCounter][iRing]  /= nCell;
        N2O[nCounter][iRing]   /= nCell;
        N[nCounter][iRing]     /= nCell;
        PAN[nCounter][iRing]   /= nCell;
        ALK4[nCounter][iRing]  /= nCell;
        MAP[nCounter][iRing]   /= nCell;
        MPN[nCounter][iRing]   /= nCell;
        Cl2O2[nCounter][iRing] /= nCell;
        ETP[nCounter][iRing]   /= nCell;
        HNO2[nCounter][iRing]  /= nCell;
        C3H8[nCounter][iRing]  /= nCell;
        RA3P[nCounter][iRing]  /= nCell;
        RB3P[nCounter][iRing]  /= nCell;
        OClO[nCounter][iRing]  /= nCell;
        ClNO2[nCounter][iRing] /= nCell;
        ISOP[nCounter][iRing]  /= nCell;
        HNO4[nCounter][iRing]  /= nCell;
        MAOP[nCounter][iRing]  /= nCell;
        MP[nCounter][iRing]    /= nCell;
        ClOO[nCounter][iRing]  /= nCell;
        RP[nCounter][iRing]    /= nCell;
        BrCl[nCounter][iRing]  /= nCell;
        PP[nCounter][iRing]    /= nCell;
        PRPN[nCounter][iRing]  /= nCell;
        SO4[nCounter][iRing]   /= nCell;
        Br2[nCounter][iRing]   /= nCell;
        ETHLN[nCounter][iRing] /= nCell;
        MVKN[nCounter][iRing]  /= nCell;
        R4P[nCounter][iRing]   /= nCell;
        C2H6[nCounter][iRing]  /= nCell;
        RIP[nCounter][iRing]   /= nCell;
        VRP[nCounter][iRing]   /= nCell;
        ATOOH[nCounter][iRing] /= nCell;
        IAP[nCounter][iRing]   /= nCell;
        DHMOB[nCounter][iRing] /= nCell;
        MOBA[nCounter][iRing]  /= nCell;
        MRP[nCounter][iRing]   /= nCell;
        N2O5[nCounter][iRing]  /= nCell;
        ISNOHOO[nCounter][iRing]/= nCell;
        ISNP[nCounter][iRing]  /= nCell;
        ISOPNB[nCounter][iRing]/= nCell;
        IEPOXOO[nCounter][iRing]/= nCell;
        MACRNO2[nCounter][iRing]/= nCell;
        ROH[nCounter][iRing]   /= nCell;
        MOBAOO[nCounter][iRing]/= nCell;
        DIBOO[nCounter][iRing] /= nCell;
        PMN[nCounter][iRing]   /= nCell;
        ISNOOB[nCounter][iRing]/= nCell;
        INPN[nCounter][iRing]  /= nCell;
        H[nCounter][iRing]     /= nCell;
        BrNO3[nCounter][iRing] /= nCell;
        PRPE[nCounter][iRing]  /= nCell;
        MVKOO[nCounter][iRing] /= nCell;
        Cl2[nCounter][iRing]   /= nCell;
        ISOPND[nCounter][iRing]/= nCell;
        HOBr[nCounter][iRing]  /= nCell;
        A3O2[nCounter][iRing]  /= nCell;
        PROPNN[nCounter][iRing]/= nCell;
        GLYX[nCounter][iRing]  /= nCell;
        MAOPO2[nCounter][iRing]/= nCell;
        CH4[nCounter][iRing]   /= nCell;
        GAOO[nCounter][iRing]  /= nCell;
        B3O2[nCounter][iRing]  /= nCell;
        ACET[nCounter][iRing]  /= nCell;
        MACRN[nCounter][iRing] /= nCell;
        CH2OO[nCounter][iRing] /= nCell;
        MGLYOO[nCounter][iRing]/= nCell;
        VRO2[nCounter][iRing]  /= nCell;
        MGLOO[nCounter][iRing] /= nCell;
        MACROO[nCounter][iRing]/= nCell;
        PO2[nCounter][iRing]   /= nCell;
        CH3CHOO[nCounter][iRing]/= nCell;
        MAN2[nCounter][iRing]  /= nCell;
        ISNOOA[nCounter][iRing]/= nCell;
        H2O2[nCounter][iRing]  /= nCell;
        PRN1[nCounter][iRing]  /= nCell;
        ETO2[nCounter][iRing]  /= nCell;
        KO2[nCounter][iRing]   /= nCell;
        RCO3[nCounter][iRing]  /= nCell;
        HC5OO[nCounter][iRing] /= nCell;
        GLYC[nCounter][iRing]  /= nCell;
        ClNO3[nCounter][iRing] /= nCell;
        RIO2[nCounter][iRing]  /= nCell;
        R4N1[nCounter][iRing]  /= nCell;
        HOCl[nCounter][iRing]  /= nCell;
        ATO2[nCounter][iRing]  /= nCell;
        HNO3[nCounter][iRing]  /= nCell;
        ISN1[nCounter][iRing]  /= nCell;
        MAO3[nCounter][iRing]  /= nCell;
        MRO2[nCounter][iRing]  /= nCell;
        INO2[nCounter][iRing]  /= nCell;
        HAC[nCounter][iRing]   /= nCell;
        HC5[nCounter][iRing]   /= nCell;
        MGLY[nCounter][iRing]  /= nCell;
        ISOPNBO2[nCounter][iRing]/= nCell;
        ISOPNDO2[nCounter][iRing]/= nCell;
        R4O2[nCounter][iRing]  /= nCell;
        R4N2[nCounter][iRing]  /= nCell;
        BrO[nCounter][iRing]   /= nCell;
        RCHO[nCounter][iRing]  /= nCell;
        MEK[nCounter][iRing]   /= nCell;
        ClO[nCounter][iRing]   /= nCell;
        MACR[nCounter][iRing]  /= nCell;
        SO2[nCounter][iRing]   /= nCell;
        MVK[nCounter][iRing]   /= nCell;
        ALD2[nCounter][iRing]  /= nCell;
        MCO3[nCounter][iRing]  /= nCell;
        CH2O[nCounter][iRing]  /= nCell;
        H2O[nCounter][iRing]   /= nCell;
        Br[nCounter][iRing]    /= nCell;
        NO[nCounter][iRing]    /= nCell;
        NO3[nCounter][iRing]   /= nCell;
        Cl[nCounter][iRing]    /= nCell;
        O[nCounter][iRing]     /= nCell;
        O1D[nCounter][iRing]   /= nCell;
        O3[nCounter][iRing]    /= nCell;
        HO2[nCounter][iRing]   /= nCell;
        NO2[nCounter][iRing]   /= nCell;
        OH[nCounter][iRing]    /= nCell;
        HBr[nCounter][iRing]   /= nCell;
        HCl[nCounter][iRing]   /= nCell;
        CO[nCounter][iRing]    /= nCell;
        MO2[nCounter][iRing]   /= nCell;

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

unsigned int SpeciesArray::GetnRing() const
{

    return nRing;

} /* End of SpeciesArray::GetnRing */

unsigned int SpeciesArray::GetnTime() const
{

    return nTime;

} /* End of SpeciesArray::GetnTime */


/* End of Species.cpp */
