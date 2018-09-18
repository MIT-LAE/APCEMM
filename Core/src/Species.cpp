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
    std::vector<int> nMap = m.getnMap();
    unsigned int nx_max, ny_max;
    unsigned int fMult = 1;

#if Y_SYMMETRY
    nx_max = std::ceil(NX/2);
    fMult *= 2;
#else
    nx_max = NX;
#endif /* Y_SYMMETRY */

#if X_SYMMETRY
    ny_max = std::ceil(NY/2);
    fMult *= 2;
#else
    ny_max = NY;
#endif /* X_SYMMETRY */

    for ( unsigned int iRing = 0; iRing < nRing; iRing++ ) {
        for ( unsigned int iNx = 0; iNx < nx_max; iNx++ ) {
            for ( unsigned int jNy = 0; jNy < ny_max; jNy++ ) {
                /* Sum concentrations over rings */
                CO2[nCounter][iRing]   += Data.CO2[jNy][iNx] * map[iRing][jNy][iNx];
                PPN[nCounter][iRing]   += Data.PPN[jNy][iNx] * map[iRing][jNy][iNx];
                BrNO2[nCounter][iRing] += Data.BrNO2[jNy][iNx] * map[iRing][jNy][iNx];
                IEPOX[nCounter][iRing] += Data.IEPOX[jNy][iNx] * map[iRing][jNy][iNx];
                PMNN[nCounter][iRing]  += Data.PMNN[jNy][iNx] * map[iRing][jNy][iNx];
                N2O[nCounter][iRing]   += Data.N2O[jNy][iNx] * map[iRing][jNy][iNx];
                N[nCounter][iRing]     += Data.N[jNy][iNx] * map[iRing][jNy][iNx];
                PAN[nCounter][iRing]   += Data.PAN[jNy][iNx] * map[iRing][jNy][iNx];
                ALK4[nCounter][iRing]  += Data.ALK4[jNy][iNx] * map[iRing][jNy][iNx];
                MAP[nCounter][iRing]   += Data.MAP[jNy][iNx] * map[iRing][jNy][iNx];
                MPN[nCounter][iRing]   += Data.MPN[jNy][iNx] * map[iRing][jNy][iNx];
                Cl2O2[nCounter][iRing] += Data.Cl2O2[jNy][iNx] * map[iRing][jNy][iNx];
                ETP[nCounter][iRing]   += Data.ETP[jNy][iNx] * map[iRing][jNy][iNx];
                HNO2[nCounter][iRing]  += Data.HNO2[jNy][iNx] * map[iRing][jNy][iNx];
                C3H8[nCounter][iRing]  += Data.C3H8[jNy][iNx] * map[iRing][jNy][iNx];
                RA3P[nCounter][iRing]  += Data.RA3P[jNy][iNx] * map[iRing][jNy][iNx];
                RB3P[nCounter][iRing]  += Data.RB3P[jNy][iNx] * map[iRing][jNy][iNx];
                OClO[nCounter][iRing]  += Data.OClO[jNy][iNx] * map[iRing][jNy][iNx];
                ClNO2[nCounter][iRing] += Data.ClNO2[jNy][iNx] * map[iRing][jNy][iNx];
                ISOP[nCounter][iRing]  += Data.ISOP[jNy][iNx] * map[iRing][jNy][iNx];
                HNO4[nCounter][iRing]  += Data.HNO4[jNy][iNx] * map[iRing][jNy][iNx];
                MAOP[nCounter][iRing]  += Data.MAOP[jNy][iNx] * map[iRing][jNy][iNx];
                MP[nCounter][iRing]    += Data.MP[jNy][iNx] * map[iRing][jNy][iNx];
                ClOO[nCounter][iRing]  += Data.ClOO[jNy][iNx] * map[iRing][jNy][iNx];
                RP[nCounter][iRing]    += Data.RP[jNy][iNx] * map[iRing][jNy][iNx];
                BrCl[nCounter][iRing]  += Data.BrCl[jNy][iNx] * map[iRing][jNy][iNx];
                PP[nCounter][iRing]    += Data.PP[jNy][iNx] * map[iRing][jNy][iNx];
                PRPN[nCounter][iRing]  += Data.PRPN[jNy][iNx] * map[iRing][jNy][iNx];
                SO4[nCounter][iRing]   += Data.SO4[jNy][iNx] * map[iRing][jNy][iNx];
                Br2[nCounter][iRing]   += Data.Br2[jNy][iNx] * map[iRing][jNy][iNx];
                ETHLN[nCounter][iRing] += Data.ETHLN[jNy][iNx] * map[iRing][jNy][iNx];
                MVKN[nCounter][iRing]  += Data.MVKN[jNy][iNx] * map[iRing][jNy][iNx];
                R4P[nCounter][iRing]   += Data.R4P[jNy][iNx] * map[iRing][jNy][iNx];
                C2H6[nCounter][iRing]  += Data.C2H6[jNy][iNx] * map[iRing][jNy][iNx];
                RIP[nCounter][iRing]   += Data.RIP[jNy][iNx] * map[iRing][jNy][iNx];
                VRP[nCounter][iRing]   += Data.VRP[jNy][iNx] * map[iRing][jNy][iNx];
                ATOOH[nCounter][iRing] += Data.ATOOH[jNy][iNx] * map[iRing][jNy][iNx];
                IAP[nCounter][iRing]   += Data.IAP[jNy][iNx] * map[iRing][jNy][iNx];
                DHMOB[nCounter][iRing] += Data.DHMOB[jNy][iNx] * map[iRing][jNy][iNx];
                MOBA[nCounter][iRing]  += Data.MOBA[jNy][iNx] * map[iRing][jNy][iNx];
                MRP[nCounter][iRing]   += Data.MRP[jNy][iNx] * map[iRing][jNy][iNx];
                N2O5[nCounter][iRing]  += Data.N2O5[jNy][iNx] * map[iRing][jNy][iNx];
                ISNOHOO[nCounter][iRing]+= Data.ISNOHOO[jNy][iNx] * map[iRing][jNy][iNx];
                ISNP[nCounter][iRing]  += Data.ISNP[jNy][iNx] * map[iRing][jNy][iNx];
                ISOPNB[nCounter][iRing]+= Data.ISOPNB[jNy][iNx] * map[iRing][jNy][iNx];
                IEPOXOO[nCounter][iRing]+= Data.IEPOXOO[jNy][iNx] * map[iRing][jNy][iNx];
                MACRNO2[nCounter][iRing]+= Data.MACRNO2[jNy][iNx] * map[iRing][jNy][iNx];
                ROH[nCounter][iRing]   += Data.ROH[jNy][iNx] * map[iRing][jNy][iNx];
                MOBAOO[nCounter][iRing]+= Data.MOBAOO[jNy][iNx] * map[iRing][jNy][iNx];
                DIBOO[nCounter][iRing] += Data.DIBOO[jNy][iNx] * map[iRing][jNy][iNx];
                PMN[nCounter][iRing]   += Data.PMN[jNy][iNx] * map[iRing][jNy][iNx];
                ISNOOB[nCounter][iRing]+= Data.ISNOOB[jNy][iNx] * map[iRing][jNy][iNx];
                INPN[nCounter][iRing]  += Data.INPN[jNy][iNx] * map[iRing][jNy][iNx];
                H[nCounter][iRing]     += Data.H[jNy][iNx] * map[iRing][jNy][iNx];
                BrNO3[nCounter][iRing] += Data.BrNO3[jNy][iNx] * map[iRing][jNy][iNx];
                PRPE[nCounter][iRing]  += Data.PRPE[jNy][iNx] * map[iRing][jNy][iNx];
                MVKOO[nCounter][iRing] += Data.MVKOO[jNy][iNx] * map[iRing][jNy][iNx];
                Cl2[nCounter][iRing]   += Data.Cl2[jNy][iNx] * map[iRing][jNy][iNx];
                ISOPND[nCounter][iRing]+= Data.ISOPND[jNy][iNx] * map[iRing][jNy][iNx];
                HOBr[nCounter][iRing]  += Data.HOBr[jNy][iNx] * map[iRing][jNy][iNx];
                A3O2[nCounter][iRing]  += Data.A3O2[jNy][iNx] * map[iRing][jNy][iNx];
                PROPNN[nCounter][iRing]+= Data.PROPNN[jNy][iNx] * map[iRing][jNy][iNx];
                GLYX[nCounter][iRing]  += Data.GLYX[jNy][iNx] * map[iRing][jNy][iNx];
                MAOPO2[nCounter][iRing]+= Data.MAOPO2[jNy][iNx] * map[iRing][jNy][iNx];
                CH4[nCounter][iRing]   += Data.CH4[jNy][iNx] * map[iRing][jNy][iNx];
                GAOO[nCounter][iRing]  += Data.GAOO[jNy][iNx] * map[iRing][jNy][iNx];
                B3O2[nCounter][iRing]  += Data.B3O2[jNy][iNx] * map[iRing][jNy][iNx];
                ACET[nCounter][iRing]  += Data.ACET[jNy][iNx] * map[iRing][jNy][iNx];
                MACRN[nCounter][iRing] += Data.MACRN[jNy][iNx] * map[iRing][jNy][iNx];
                CH2OO[nCounter][iRing] += Data.CH2OO[jNy][iNx] * map[iRing][jNy][iNx];
                MGLYOO[nCounter][iRing]+= Data.MGLYOO[jNy][iNx] * map[iRing][jNy][iNx];
                VRO2[nCounter][iRing]  += Data.VRO2[jNy][iNx] * map[iRing][jNy][iNx];
                MGLOO[nCounter][iRing] += Data.MGLOO[jNy][iNx] * map[iRing][jNy][iNx];
                MACROO[nCounter][iRing]+= Data.MACROO[jNy][iNx] * map[iRing][jNy][iNx];
                PO2[nCounter][iRing]   += Data.PO2[jNy][iNx] * map[iRing][jNy][iNx];
                CH3CHOO[nCounter][iRing]+= Data.CH3CHOO[jNy][iNx] * map[iRing][jNy][iNx];
                MAN2[nCounter][iRing]  += Data.MAN2[jNy][iNx] * map[iRing][jNy][iNx];
                ISNOOA[nCounter][iRing]+= Data.ISNOOA[jNy][iNx] * map[iRing][jNy][iNx];
                H2O2[nCounter][iRing]  += Data.H2O2[jNy][iNx] * map[iRing][jNy][iNx];
                PRN1[nCounter][iRing]  += Data.PRN1[jNy][iNx] * map[iRing][jNy][iNx];
                ETO2[nCounter][iRing]  += Data.ETO2[jNy][iNx] * map[iRing][jNy][iNx];
                KO2[nCounter][iRing]   += Data.KO2[jNy][iNx] * map[iRing][jNy][iNx];
                RCO3[nCounter][iRing]  += Data.RCO3[jNy][iNx] * map[iRing][jNy][iNx];
                HC5OO[nCounter][iRing] += Data.HC5OO[jNy][iNx] * map[iRing][jNy][iNx];
                GLYC[nCounter][iRing]  += Data.GLYC[jNy][iNx] * map[iRing][jNy][iNx];
                ClNO3[nCounter][iRing] += Data.ClNO3[jNy][iNx] * map[iRing][jNy][iNx];
                RIO2[nCounter][iRing]  += Data.RIO2[jNy][iNx] * map[iRing][jNy][iNx];
                R4N1[nCounter][iRing]  += Data.R4N1[jNy][iNx] * map[iRing][jNy][iNx];
                HOCl[nCounter][iRing]  += Data.HOCl[jNy][iNx] * map[iRing][jNy][iNx];
                ATO2[nCounter][iRing]  += Data.ATO2[jNy][iNx] * map[iRing][jNy][iNx];
                HNO3[nCounter][iRing]  += Data.HNO3[jNy][iNx] * map[iRing][jNy][iNx];
                ISN1[nCounter][iRing]  += Data.ISN1[jNy][iNx] * map[iRing][jNy][iNx];
                MAO3[nCounter][iRing]  += Data.MAO3[jNy][iNx] * map[iRing][jNy][iNx];
                MRO2[nCounter][iRing]  += Data.MRO2[jNy][iNx] * map[iRing][jNy][iNx];
                INO2[nCounter][iRing]  += Data.INO2[jNy][iNx] * map[iRing][jNy][iNx];
                HAC[nCounter][iRing]   += Data.HAC[jNy][iNx] * map[iRing][jNy][iNx];
                HC5[nCounter][iRing]   += Data.HC5[jNy][iNx] * map[iRing][jNy][iNx];
                MGLY[nCounter][iRing]  += Data.MGLY[jNy][iNx] * map[iRing][jNy][iNx];
                ISOPNBO2[nCounter][iRing]+= Data.ISOPNBO2[jNy][iNx] * map[iRing][jNy][iNx];
                ISOPNDO2[nCounter][iRing]+= Data.ISOPNDO2[jNy][iNx] * map[iRing][jNy][iNx];
                R4O2[nCounter][iRing]  += Data.R4O2[jNy][iNx] * map[iRing][jNy][iNx];
                R4N2[nCounter][iRing]  += Data.R4N2[jNy][iNx] * map[iRing][jNy][iNx];
                BrO[nCounter][iRing]   += Data.BrO[jNy][iNx] * map[iRing][jNy][iNx];
                RCHO[nCounter][iRing]  += Data.RCHO[jNy][iNx] * map[iRing][jNy][iNx];
                MEK[nCounter][iRing]   += Data.MEK[jNy][iNx] * map[iRing][jNy][iNx];
                ClO[nCounter][iRing]   += Data.ClO[jNy][iNx] * map[iRing][jNy][iNx];
                MACR[nCounter][iRing]  += Data.MACR[jNy][iNx] * map[iRing][jNy][iNx];
                SO2[nCounter][iRing]   += Data.SO2[jNy][iNx] * map[iRing][jNy][iNx];
                MVK[nCounter][iRing]   += Data.MVK[jNy][iNx] * map[iRing][jNy][iNx];
                ALD2[nCounter][iRing]  += Data.ALD2[jNy][iNx] * map[iRing][jNy][iNx];
                MCO3[nCounter][iRing]  += Data.MCO3[jNy][iNx] * map[iRing][jNy][iNx];
                CH2O[nCounter][iRing]  += Data.CH2O[jNy][iNx] * map[iRing][jNy][iNx];
                H2O[nCounter][iRing]   += Data.H2O[jNy][iNx] * map[iRing][jNy][iNx];
                Br[nCounter][iRing]    += Data.Br[jNy][iNx] * map[iRing][jNy][iNx];
                NO[nCounter][iRing]    += Data.NO[jNy][iNx] * map[iRing][jNy][iNx];
                NO3[nCounter][iRing]   += Data.NO3[jNy][iNx] * map[iRing][jNy][iNx];
                Cl[nCounter][iRing]    += Data.Cl[jNy][iNx] * map[iRing][jNy][iNx];
                O[nCounter][iRing]     += Data.O[jNy][iNx] * map[iRing][jNy][iNx];
                O1D[nCounter][iRing]   += Data.O1D[jNy][iNx] * map[iRing][jNy][iNx];
                O3[nCounter][iRing]    += Data.O3[jNy][iNx] * map[iRing][jNy][iNx];
                HO2[nCounter][iRing]   += Data.HO2[jNy][iNx] * map[iRing][jNy][iNx];
                NO2[nCounter][iRing]   += Data.NO2[jNy][iNx] * map[iRing][jNy][iNx];
                OH[nCounter][iRing]    += Data.OH[jNy][iNx] * map[iRing][jNy][iNx];
                HBr[nCounter][iRing]   += Data.HBr[jNy][iNx] * map[iRing][jNy][iNx];
                HCl[nCounter][iRing]   += Data.HCl[jNy][iNx] * map[iRing][jNy][iNx];
                CO[nCounter][iRing]    += Data.CO[jNy][iNx] * map[iRing][jNy][iNx];
                MO2[nCounter][iRing]   += Data.MO2[jNy][iNx] * map[iRing][jNy][iNx];

            }
        }

        CO2[nCounter][iRing]   /= ( nMap[iRing] / fMult );
        PPN[nCounter][iRing]   /= ( nMap[iRing] / fMult );
        BrNO2[nCounter][iRing] /= ( nMap[iRing] / fMult );
        IEPOX[nCounter][iRing] /= ( nMap[iRing] / fMult );
        PMNN[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        N2O[nCounter][iRing]   /= ( nMap[iRing] / fMult );
        N[nCounter][iRing]     /= ( nMap[iRing] / fMult );
        PAN[nCounter][iRing]   /= ( nMap[iRing] / fMult );
        ALK4[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        MAP[nCounter][iRing]   /= ( nMap[iRing] / fMult );
        MPN[nCounter][iRing]   /= ( nMap[iRing] / fMult );
        Cl2O2[nCounter][iRing] /= ( nMap[iRing] / fMult );
        ETP[nCounter][iRing]   /= ( nMap[iRing] / fMult );
        HNO2[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        C3H8[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        RA3P[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        RB3P[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        OClO[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        ClNO2[nCounter][iRing] /= ( nMap[iRing] / fMult );
        ISOP[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        HNO4[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        MAOP[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        MP[nCounter][iRing]    /= ( nMap[iRing] / fMult );
        ClOO[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        RP[nCounter][iRing]    /= ( nMap[iRing] / fMult );
        BrCl[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        PP[nCounter][iRing]    /= ( nMap[iRing] / fMult );
        PRPN[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        SO4[nCounter][iRing]   /= ( nMap[iRing] / fMult );
        Br2[nCounter][iRing]   /= ( nMap[iRing] / fMult );
        ETHLN[nCounter][iRing] /= ( nMap[iRing] / fMult );
        MVKN[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        R4P[nCounter][iRing]   /= ( nMap[iRing] / fMult );
        C2H6[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        RIP[nCounter][iRing]   /= ( nMap[iRing] / fMult );
        VRP[nCounter][iRing]   /= ( nMap[iRing] / fMult );
        ATOOH[nCounter][iRing] /= ( nMap[iRing] / fMult );
        IAP[nCounter][iRing]   /= ( nMap[iRing] / fMult );
        DHMOB[nCounter][iRing] /= ( nMap[iRing] / fMult );
        MOBA[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        MRP[nCounter][iRing]   /= ( nMap[iRing] / fMult );
        N2O5[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        ISNOHOO[nCounter][iRing]/= ( nMap[iRing] / fMult );
        ISNP[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        ISOPNB[nCounter][iRing]/= ( nMap[iRing] / fMult );
        IEPOXOO[nCounter][iRing]/= ( nMap[iRing] / fMult );
        MACRNO2[nCounter][iRing]/= ( nMap[iRing] / fMult );
        ROH[nCounter][iRing]   /= ( nMap[iRing] / fMult );
        MOBAOO[nCounter][iRing]/= ( nMap[iRing] / fMult );
        DIBOO[nCounter][iRing] /= ( nMap[iRing] / fMult );
        PMN[nCounter][iRing]   /= ( nMap[iRing] / fMult );
        ISNOOB[nCounter][iRing]/= ( nMap[iRing] / fMult );
        INPN[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        H[nCounter][iRing]     /= ( nMap[iRing] / fMult );
        BrNO3[nCounter][iRing] /= ( nMap[iRing] / fMult );
        PRPE[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        MVKOO[nCounter][iRing] /= ( nMap[iRing] / fMult );
        Cl2[nCounter][iRing]   /= ( nMap[iRing] / fMult );
        ISOPND[nCounter][iRing]/= ( nMap[iRing] / fMult );
        HOBr[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        A3O2[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        PROPNN[nCounter][iRing]/= ( nMap[iRing] / fMult );
        GLYX[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        MAOPO2[nCounter][iRing]/= ( nMap[iRing] / fMult );
        CH4[nCounter][iRing]   /= ( nMap[iRing] / fMult );
        GAOO[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        B3O2[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        ACET[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        MACRN[nCounter][iRing] /= ( nMap[iRing] / fMult );
        CH2OO[nCounter][iRing] /= ( nMap[iRing] / fMult );
        MGLYOO[nCounter][iRing]/= ( nMap[iRing] / fMult );
        VRO2[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        MGLOO[nCounter][iRing] /= ( nMap[iRing] / fMult );
        MACROO[nCounter][iRing]/= ( nMap[iRing] / fMult );
        PO2[nCounter][iRing]   /= ( nMap[iRing] / fMult );
        CH3CHOO[nCounter][iRing]/= ( nMap[iRing] / fMult );
        MAN2[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        ISNOOA[nCounter][iRing]/= ( nMap[iRing] / fMult );
        H2O2[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        PRN1[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        ETO2[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        KO2[nCounter][iRing]   /= ( nMap[iRing] / fMult );
        RCO3[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        HC5OO[nCounter][iRing] /= ( nMap[iRing] / fMult );
        GLYC[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        ClNO3[nCounter][iRing] /= ( nMap[iRing] / fMult );
        RIO2[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        R4N1[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        HOCl[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        ATO2[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        HNO3[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        ISN1[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        MAO3[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        MRO2[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        INO2[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        HAC[nCounter][iRing]   /= ( nMap[iRing] / fMult );
        HC5[nCounter][iRing]   /= ( nMap[iRing] / fMult );
        MGLY[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        ISOPNBO2[nCounter][iRing]/= ( nMap[iRing] / fMult );
        ISOPNDO2[nCounter][iRing]/= ( nMap[iRing] / fMult );
        R4O2[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        R4N2[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        BrO[nCounter][iRing]   /= ( nMap[iRing] / fMult );
        RCHO[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        MEK[nCounter][iRing]   /= ( nMap[iRing] / fMult );
        ClO[nCounter][iRing]   /= ( nMap[iRing] / fMult );
        MACR[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        SO2[nCounter][iRing]   /= ( nMap[iRing] / fMult );
        MVK[nCounter][iRing]   /= ( nMap[iRing] / fMult );
        ALD2[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        MCO3[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        CH2O[nCounter][iRing]  /= ( nMap[iRing] / fMult );
        H2O[nCounter][iRing]   /= ( nMap[iRing] / fMult );
        Br[nCounter][iRing]    /= ( nMap[iRing] / fMult );
        NO[nCounter][iRing]    /= ( nMap[iRing] / fMult );
        NO3[nCounter][iRing]   /= ( nMap[iRing] / fMult );
        Cl[nCounter][iRing]    /= ( nMap[iRing] / fMult );
        O[nCounter][iRing]     /= ( nMap[iRing] / fMult );
        O1D[nCounter][iRing]   /= ( nMap[iRing] / fMult );
        O3[nCounter][iRing]    /= ( nMap[iRing] / fMult );
        HO2[nCounter][iRing]   /= ( nMap[iRing] / fMult );
        NO2[nCounter][iRing]   /= ( nMap[iRing] / fMult );
        OH[nCounter][iRing]    /= ( nMap[iRing] / fMult );
        HBr[nCounter][iRing]   /= ( nMap[iRing] / fMult );
        HCl[nCounter][iRing]   /= ( nMap[iRing] / fMult );
        CO[nCounter][iRing]    /= ( nMap[iRing] / fMult );
        MO2[nCounter][iRing]   /= ( nMap[iRing] / fMult );

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
