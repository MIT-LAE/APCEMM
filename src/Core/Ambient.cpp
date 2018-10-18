/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Ambient Program File                                             */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : Ambient.cpp                               */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "Core/Ambient.hpp"

Ambient::Ambient( )
{

    /* Default Constructor */

} /* End of Ambient::Ambient */

Ambient::Ambient( unsigned int nTime_, std::vector<double> ambientVector, std::vector<std::vector<double> > aerVector )
{

    /* Constructor */

    nTime = nTime_;

    CO2.assign      ( nTime, ambientVector[  0] );
    PPN.assign      ( nTime, ambientVector[  1] );
    BrNO2.assign    ( nTime, ambientVector[  2] );
    IEPOX.assign    ( nTime, ambientVector[  3] );
    PMNN.assign     ( nTime, ambientVector[  4] );
    N2O.assign      ( nTime, ambientVector[  5] );
    N.assign        ( nTime, ambientVector[  6] );
    PAN.assign      ( nTime, ambientVector[  7] );
    ALK4.assign     ( nTime, ambientVector[  8] );
    MAP.assign      ( nTime, ambientVector[  9] );
    MPN.assign      ( nTime, ambientVector[ 10] );
    Cl2O2.assign    ( nTime, ambientVector[ 11] );
    ETP.assign      ( nTime, ambientVector[ 12] );
    HNO2.assign     ( nTime, ambientVector[ 13] );
    C3H8.assign     ( nTime, ambientVector[ 14] );
    RA3P.assign     ( nTime, ambientVector[ 15] );
    RB3P.assign     ( nTime, ambientVector[ 16] );
    OClO.assign     ( nTime, ambientVector[ 17] );
    ClNO2.assign    ( nTime, ambientVector[ 18] );
    ISOP.assign     ( nTime, ambientVector[ 19] );
    HNO4.assign     ( nTime, ambientVector[ 20] );
    MAOP.assign     ( nTime, ambientVector[ 21] );
    MP.assign       ( nTime, ambientVector[ 22] );
    ClOO.assign     ( nTime, ambientVector[ 23] );
    RP.assign       ( nTime, ambientVector[ 24] );
    BrCl.assign     ( nTime, ambientVector[ 25] );
    PP.assign       ( nTime, ambientVector[ 26] );
    PRPN.assign     ( nTime, ambientVector[ 27] );
    SO4.assign      ( nTime, ambientVector[ 28] );
    Br2.assign      ( nTime, ambientVector[ 29] );
    ETHLN.assign    ( nTime, ambientVector[ 30] );
    MVKN.assign     ( nTime, ambientVector[ 31] );
    R4P.assign      ( nTime, ambientVector[ 32] );
    C2H6.assign     ( nTime, ambientVector[ 33] );
    RIP.assign      ( nTime, ambientVector[ 34] );
    VRP.assign      ( nTime, ambientVector[ 35] );
    ATOOH.assign    ( nTime, ambientVector[ 36] );
    IAP.assign      ( nTime, ambientVector[ 37] );
    DHMOB.assign    ( nTime, ambientVector[ 38] );
    MOBA.assign     ( nTime, ambientVector[ 39] );
    MRP.assign      ( nTime, ambientVector[ 40] );
    N2O5.assign     ( nTime, ambientVector[ 41] );
    ISNOHOO.assign  ( nTime, ambientVector[ 42] );
    ISNP.assign     ( nTime, ambientVector[ 43] );
    ISOPNB.assign   ( nTime, ambientVector[ 44] );
    IEPOXOO.assign  ( nTime, ambientVector[ 45] );
    MACRNO2.assign  ( nTime, ambientVector[ 46] );
    ROH.assign      ( nTime, ambientVector[ 47] );
    MOBAOO.assign   ( nTime, ambientVector[ 48] );
    DIBOO.assign    ( nTime, ambientVector[ 49] );
    PMN.assign      ( nTime, ambientVector[ 50] );
    ISNOOB.assign   ( nTime, ambientVector[ 51] );
    INPN.assign     ( nTime, ambientVector[ 52] );
    H.assign        ( nTime, ambientVector[ 53] );
    BrNO3.assign    ( nTime, ambientVector[ 54] );
    PRPE.assign     ( nTime, ambientVector[ 55] );
    MVKOO.assign    ( nTime, ambientVector[ 56] );
    Cl2.assign      ( nTime, ambientVector[ 57] );
    ISOPND.assign   ( nTime, ambientVector[ 58] );
    HOBr.assign     ( nTime, ambientVector[ 59] );
    A3O2.assign     ( nTime, ambientVector[ 60] );
    PROPNN.assign   ( nTime, ambientVector[ 61] );
    GLYX.assign     ( nTime, ambientVector[ 62] );
    MAOPO2.assign   ( nTime, ambientVector[ 63] );
    CH4.assign      ( nTime, ambientVector[ 64] );
    GAOO.assign     ( nTime, ambientVector[ 65] );
    B3O2.assign     ( nTime, ambientVector[ 66] );
    ACET.assign     ( nTime, ambientVector[ 67] );
    MACRN.assign    ( nTime, ambientVector[ 68] );
    CH2OO.assign    ( nTime, ambientVector[ 69] );
    MGLYOO.assign   ( nTime, ambientVector[ 70] );
    VRO2.assign     ( nTime, ambientVector[ 71] );
    MGLOO.assign    ( nTime, ambientVector[ 72] );
    MACROO.assign   ( nTime, ambientVector[ 73] );
    PO2.assign      ( nTime, ambientVector[ 74] );
    CH3CHOO.assign  ( nTime, ambientVector[ 75] );
    MAN2.assign     ( nTime, ambientVector[ 76] );
    ISNOOA.assign   ( nTime, ambientVector[ 77] );
    H2O2.assign     ( nTime, ambientVector[ 78] );
    PRN1.assign     ( nTime, ambientVector[ 79] );
    ETO2.assign     ( nTime, ambientVector[ 80] );
    KO2.assign      ( nTime, ambientVector[ 81] );
    RCO3.assign     ( nTime, ambientVector[ 82] );
    HC5OO.assign    ( nTime, ambientVector[ 83] );
    GLYC.assign     ( nTime, ambientVector[ 84] );
    ClNO3.assign    ( nTime, ambientVector[ 85] );
    RIO2.assign     ( nTime, ambientVector[ 86] );
    R4N1.assign     ( nTime, ambientVector[ 87] );
    HOCl.assign     ( nTime, ambientVector[ 88] );
    ATO2.assign     ( nTime, ambientVector[ 89] );
    HNO3.assign     ( nTime, ambientVector[ 90] );
    ISN1.assign     ( nTime, ambientVector[ 91] );
    MAO3.assign     ( nTime, ambientVector[ 92] );
    MRO2.assign     ( nTime, ambientVector[ 93] );
    INO2.assign     ( nTime, ambientVector[ 94] );
    HAC.assign      ( nTime, ambientVector[ 95] );
    HC5.assign      ( nTime, ambientVector[ 96] );
    MGLY.assign     ( nTime, ambientVector[ 97] );
    ISOPNBO2.assign ( nTime, ambientVector[ 98] );
    ISOPNDO2.assign ( nTime, ambientVector[ 99] );
    R4O2.assign     ( nTime, ambientVector[100] );
    R4N2.assign     ( nTime, ambientVector[101] );
    BrO.assign      ( nTime, ambientVector[102] );
    RCHO.assign     ( nTime, ambientVector[103] );
    MEK.assign      ( nTime, ambientVector[104] );
    ClO.assign      ( nTime, ambientVector[105] );
    MACR.assign     ( nTime, ambientVector[106] );
    SO2.assign      ( nTime, ambientVector[107] );
    MVK.assign      ( nTime, ambientVector[108] );
    ALD2.assign     ( nTime, ambientVector[109] );
    MCO3.assign     ( nTime, ambientVector[110] );
    CH2O.assign     ( nTime, ambientVector[111] );
    H2O.assign      ( nTime, ambientVector[112] );
    Br.assign       ( nTime, ambientVector[113] );
    NO.assign       ( nTime, ambientVector[114] );
    NO3.assign      ( nTime, ambientVector[115] );
    Cl.assign       ( nTime, ambientVector[116] );
    O.assign        ( nTime, ambientVector[117] );
    O1D.assign      ( nTime, ambientVector[118] );
    O3.assign       ( nTime, ambientVector[119] );
    HO2.assign      ( nTime, ambientVector[120] );
    NO2.assign      ( nTime, ambientVector[121] );
    OH.assign       ( nTime, ambientVector[122] );
    HBr.assign      ( nTime, ambientVector[123] );
    HCl.assign      ( nTime, ambientVector[124] );
    CO.assign       ( nTime, ambientVector[125] );
    MO2.assign      ( nTime, ambientVector[126] );

    ACTA     = ambientVector[127];
    EOH      = ambientVector[128];
    H2       = ambientVector[129];
    HCOOH    = ambientVector[130];
    MOH      = ambientVector[131];
    N2       = ambientVector[132];
    O2       = ambientVector[133];
    RCOOH    = ambientVector[134];

    sootDens.assign( nTime, aerVector[  0][0] );
    sootRadi.assign( nTime, aerVector[  0][1] );
    iceDens.assign ( nTime, aerVector[  1][0] );
    iceRadi.assign ( nTime, aerVector[  1][1] );
    sulfDens.assign( nTime, aerVector[  2][0] );
    sulfRadi.assign( nTime, aerVector[  2][1] );

} /* End of Ambient::Ambient */

Ambient::~Ambient( )
{

    /* Destructor */

} /* End of Ambient::~Ambient */ 


Ambient::Ambient( const Ambient &a )
{

    /* Copy size */
    nTime = a.nTime;
    
    /* Copy arrays */
    CO2      = a.CO2      ;
    PPN      = a.PPN      ;
    BrNO2    = a.BrNO2    ;
    IEPOX    = a.IEPOX    ;
    PMNN     = a.PMNN     ;
    N2O      = a.N2O      ;
    N        = a.N        ;
    PAN      = a.PAN      ;
    ALK4     = a.ALK4     ;
    MAP      = a.MAP      ;
    MPN      = a.MPN      ;
    Cl2O2    = a.Cl2O2    ;
    ETP      = a.ETP      ;
    HNO2     = a.HNO2     ;
    C3H8     = a.C3H8     ;
    RA3P     = a.RA3P     ;
    RB3P     = a.RB3P     ;
    OClO     = a.OClO     ;
    ClNO2    = a.ClNO2    ;
    ISOP     = a.ISOP     ;
    HNO4     = a.HNO4     ;
    MAOP     = a.MAOP     ;
    MP       = a.MP       ;
    ClOO     = a.ClOO     ;
    RP       = a.RP       ;
    BrCl     = a.BrCl     ;
    PP       = a.PP       ;
    PRPN     = a.PRPN     ;
    SO4      = a.SO4      ;
    Br2      = a.Br2      ;
    ETHLN    = a.ETHLN    ;
    MVKN     = a.MVKN     ;
    R4P      = a.R4P      ;
    C2H6     = a.C2H6     ;
    RIP      = a.RIP      ;
    VRP      = a.VRP      ;
    ATOOH    = a.ATOOH    ;
    IAP      = a.IAP      ;
    DHMOB    = a.DHMOB    ;
    MOBA     = a.MOBA     ;
    MRP      = a.MRP      ;
    N2O5     = a.N2O5     ;
    ISNOHOO  = a.ISNOHOO  ;
    ISNP     = a.ISNP     ;
    ISOPNB   = a.ISOPNB   ;
    IEPOXOO  = a.IEPOXOO  ;
    MACRNO2  = a.MACRNO2  ;
    ROH      = a.ROH      ;
    MOBAOO   = a.MOBAOO   ;
    DIBOO    = a.DIBOO    ;
    PMN      = a.PMN      ;
    ISNOOB   = a.ISNOOB   ;
    INPN     = a.INPN     ;
    H        = a.H        ;
    BrNO3    = a.BrNO3    ;
    PRPE     = a.PRPE     ;
    MVKOO    = a.MVKOO    ;
    Cl2      = a.Cl2      ;
    ISOPND   = a.ISOPND   ;
    HOBr     = a.HOBr     ;
    A3O2     = a.A3O2     ;
    PROPNN   = a.PROPNN   ;
    GLYX     = a.GLYX     ;
    MAOPO2   = a.MAOPO2   ;
    CH4      = a.CH4      ;
    GAOO     = a.GAOO     ;
    B3O2     = a.B3O2     ;
    ACET     = a.ACET     ;
    MACRN    = a.MACRN    ;
    CH2OO    = a.CH2OO    ;
    MGLYOO   = a.MGLYOO   ;
    VRO2     = a.VRO2     ;
    MGLOO    = a.MGLOO    ;
    MACROO   = a.MACROO   ;
    PO2      = a.PO2      ;
    CH3CHOO  = a.CH3CHOO  ;
    MAN2     = a.MAN2     ;
    ISNOOA   = a.ISNOOA   ;
    H2O2     = a.H2O2     ;
    PRN1     = a.PRN1     ;
    ETO2     = a.ETO2     ;
    KO2      = a.KO2      ;
    RCO3     = a.RCO3     ;
    HC5OO    = a.HC5OO    ;
    GLYC     = a.GLYC     ;
    ClNO3    = a.ClNO3    ;
    RIO2     = a.RIO2     ;
    R4N1     = a.R4N1     ;
    HOCl     = a.HOCl     ;
    ATO2     = a.ATO2     ;
    HNO3     = a.HNO3     ;
    ISN1     = a.ISN1     ;
    MAO3     = a.MAO3     ;
    MRO2     = a.MRO2     ;
    INO2     = a.INO2     ;
    HAC      = a.HAC      ;
    HC5      = a.HC5      ;
    MGLY     = a.MGLY     ;
    ISOPNBO2 = a.ISOPNBO2 ;
    ISOPNDO2 = a.ISOPNDO2 ;
    R4O2     = a.R4O2     ;
    R4N2     = a.R4N2     ;
    BrO      = a.BrO      ;
    RCHO     = a.RCHO     ;
    MEK      = a.MEK      ;
    ClO      = a.ClO      ;
    MACR     = a.MACR     ;
    SO2      = a.SO2      ;
    MVK      = a.MVK      ;
    ALD2     = a.ALD2     ;
    MCO3     = a.MCO3     ;
    CH2O     = a.CH2O     ;
    H2O      = a.H2O      ;
    Br       = a.Br       ;
    NO       = a.NO       ;
    NO3      = a.NO3      ;
    Cl       = a.Cl       ;
    O        = a.O        ;
    O1D      = a.O1D      ;
    O3       = a.O3       ;
    HO2      = a.HO2      ;
    NO2      = a.NO2      ;
    OH       = a.OH       ;
    HBr      = a.HBr      ;
    HCl      = a.HCl      ;
    CO       = a.CO       ;
    MO2      = a.MO2      ;
    ACTA     = a.ACTA     ;
    EOH      = a.EOH      ;
    H2       = a.H2       ;
    HCOOH    = a.HCOOH    ;
    MOH      = a.MOH      ;
    N2       = a.N2       ;
    O2       = a.O2       ;
    RCOOH    = a.RCOOH    ;

} /* End of Ambient::Ambient */

Ambient& Ambient::operator=( const Ambient &a )
{

    if ( &a == this )
        return *this;

    /* Assign size */
    nTime = a.nTime;
    
    /* Assign arrays */
    CO2      = a.CO2      ;
    PPN      = a.PPN      ;
    BrNO2    = a.BrNO2    ;
    IEPOX    = a.IEPOX    ;
    PMNN     = a.PMNN     ;
    N2O      = a.N2O      ;
    N        = a.N        ;
    PAN      = a.PAN      ;
    ALK4     = a.ALK4     ;
    MAP      = a.MAP      ;
    MPN      = a.MPN      ;
    Cl2O2    = a.Cl2O2    ;
    ETP      = a.ETP      ;
    HNO2     = a.HNO2     ;
    C3H8     = a.C3H8     ;
    RA3P     = a.RA3P     ;
    RB3P     = a.RB3P     ;
    OClO     = a.OClO     ;
    ClNO2    = a.ClNO2    ;
    ISOP     = a.ISOP     ;
    HNO4     = a.HNO4     ;
    MAOP     = a.MAOP     ;
    MP       = a.MP       ;
    ClOO     = a.ClOO     ;
    RP       = a.RP       ;
    BrCl     = a.BrCl     ;
    PP       = a.PP       ;
    PRPN     = a.PRPN     ;
    SO4      = a.SO4      ;
    Br2      = a.Br2      ;
    ETHLN    = a.ETHLN    ;
    MVKN     = a.MVKN     ;
    R4P      = a.R4P      ;
    C2H6     = a.C2H6     ;
    RIP      = a.RIP      ;
    VRP      = a.VRP      ;
    ATOOH    = a.ATOOH    ;
    IAP      = a.IAP      ;
    DHMOB    = a.DHMOB    ;
    MOBA     = a.MOBA     ;
    MRP      = a.MRP      ;
    N2O5     = a.N2O5     ;
    ISNOHOO  = a.ISNOHOO  ;
    ISNP     = a.ISNP     ;
    ISOPNB   = a.ISOPNB   ;
    IEPOXOO  = a.IEPOXOO  ;
    MACRNO2  = a.MACRNO2  ;
    ROH      = a.ROH      ;
    MOBAOO   = a.MOBAOO   ;
    DIBOO    = a.DIBOO    ;
    PMN      = a.PMN      ;
    ISNOOB   = a.ISNOOB   ;
    INPN     = a.INPN     ;
    H        = a.H        ;
    BrNO3    = a.BrNO3    ;
    PRPE     = a.PRPE     ;
    MVKOO    = a.MVKOO    ;
    Cl2      = a.Cl2      ;
    ISOPND   = a.ISOPND   ;
    HOBr     = a.HOBr     ;
    A3O2     = a.A3O2     ;
    PROPNN   = a.PROPNN   ;
    GLYX     = a.GLYX     ;
    MAOPO2   = a.MAOPO2   ;
    CH4      = a.CH4      ;
    GAOO     = a.GAOO     ;
    B3O2     = a.B3O2     ;
    ACET     = a.ACET     ;
    MACRN    = a.MACRN    ;
    CH2OO    = a.CH2OO    ;
    MGLYOO   = a.MGLYOO   ;
    VRO2     = a.VRO2     ;
    MGLOO    = a.MGLOO    ;
    MACROO   = a.MACROO   ;
    PO2      = a.PO2      ;
    CH3CHOO  = a.CH3CHOO  ;
    MAN2     = a.MAN2     ;
    ISNOOA   = a.ISNOOA   ;
    H2O2     = a.H2O2     ;
    PRN1     = a.PRN1     ;
    ETO2     = a.ETO2     ;
    KO2      = a.KO2      ;
    RCO3     = a.RCO3     ;
    HC5OO    = a.HC5OO    ;
    GLYC     = a.GLYC     ;
    ClNO3    = a.ClNO3    ;
    RIO2     = a.RIO2     ;
    R4N1     = a.R4N1     ;
    HOCl     = a.HOCl     ;
    ATO2     = a.ATO2     ;
    HNO3     = a.HNO3     ;
    ISN1     = a.ISN1     ;
    MAO3     = a.MAO3     ;
    MRO2     = a.MRO2     ;
    INO2     = a.INO2     ;
    HAC      = a.HAC      ;
    HC5      = a.HC5      ;
    MGLY     = a.MGLY     ;
    ISOPNBO2 = a.ISOPNBO2 ;
    ISOPNDO2 = a.ISOPNDO2 ;
    R4O2     = a.R4O2     ;
    R4N2     = a.R4N2     ;
    BrO      = a.BrO      ;
    RCHO     = a.RCHO     ;
    MEK      = a.MEK      ;
    ClO      = a.ClO      ;
    MACR     = a.MACR     ;
    SO2      = a.SO2      ;
    MVK      = a.MVK      ;
    ALD2     = a.ALD2     ;
    MCO3     = a.MCO3     ;
    CH2O     = a.CH2O     ;
    H2O      = a.H2O      ;
    Br       = a.Br       ;
    NO       = a.NO       ;
    NO3      = a.NO3      ;
    Cl       = a.Cl       ;
    O        = a.O        ;
    O1D      = a.O1D      ;
    O3       = a.O3       ;
    HO2      = a.HO2      ;
    NO2      = a.NO2      ;
    OH       = a.OH       ;
    HBr      = a.HBr      ;
    HCl      = a.HCl      ;
    CO       = a.CO       ;
    MO2      = a.MO2      ;
    ACTA     = a.ACTA     ;
    EOH      = a.EOH      ;
    H2       = a.H2       ;
    HCOOH    = a.HCOOH    ;
    MOH      = a.MOH      ;
    N2       = a.N2       ;
    O2       = a.O2       ;
    RCOOH    = a.RCOOH    ;
    return *this;

} /* End of Ambient::operator= */

void Ambient::getData( double varArray[], double fixArray[], double aerArray[][2], unsigned int iTime ) const
{

    varArray[  0] = CO2[iTime];
    varArray[  1] = PPN[iTime];
    varArray[  2] = BrNO2[iTime];
    varArray[  3] = IEPOX[iTime];
    varArray[  4] = PMNN[iTime];
    varArray[  5] = N2O[iTime];
    varArray[  6] = N[iTime];
    varArray[  7] = PAN[iTime];
    varArray[  8] = ALK4[iTime];
    varArray[  9] = MAP[iTime];
    varArray[ 10] = MPN[iTime];
    varArray[ 11] = Cl2O2[iTime];
    varArray[ 12] = ETP[iTime];
    varArray[ 13] = HNO2[iTime];
    varArray[ 14] = C3H8[iTime];
    varArray[ 15] = RA3P[iTime];
    varArray[ 16] = RB3P[iTime];
    varArray[ 17] = OClO[iTime];
    varArray[ 18] = ClNO2[iTime];
    varArray[ 19] = ISOP[iTime];
    varArray[ 20] = HNO4[iTime];
    varArray[ 21] = MAOP[iTime];
    varArray[ 22] = MP[iTime];
    varArray[ 23] = ClOO[iTime];
    varArray[ 24] = RP[iTime];
    varArray[ 25] = BrCl[iTime];
    varArray[ 26] = PP[iTime];
    varArray[ 27] = PRPN[iTime];
    varArray[ 28] = SO4[iTime];
    varArray[ 29] = Br2[iTime];
    varArray[ 30] = ETHLN[iTime];
    varArray[ 31] = MVKN[iTime];
    varArray[ 32] = R4P[iTime];
    varArray[ 33] = C2H6[iTime];
    varArray[ 34] = RIP[iTime];
    varArray[ 35] = VRP[iTime];
    varArray[ 36] = ATOOH[iTime];
    varArray[ 37] = IAP[iTime];
    varArray[ 38] = DHMOB[iTime];
    varArray[ 39] = MOBA[iTime];
    varArray[ 40] = MRP[iTime];
    varArray[ 41] = N2O5[iTime];
    varArray[ 42] = ISNOHOO[iTime];
    varArray[ 43] = ISNP[iTime];
    varArray[ 44] = ISOPNB[iTime];
    varArray[ 45] = IEPOXOO[iTime];
    varArray[ 46] = MACRNO2[iTime];
    varArray[ 47] = ROH[iTime];
    varArray[ 48] = MOBAOO[iTime];
    varArray[ 49] = DIBOO[iTime];
    varArray[ 50] = PMN[iTime];
    varArray[ 51] = ISNOOB[iTime];
    varArray[ 52] = INPN[iTime];
    varArray[ 53] = H[iTime];
    varArray[ 54] = BrNO3[iTime];
    varArray[ 55] = PRPE[iTime];
    varArray[ 56] = MVKOO[iTime];
    varArray[ 57] = Cl2[iTime];
    varArray[ 58] = ISOPND[iTime];
    varArray[ 59] = HOBr[iTime];
    varArray[ 60] = A3O2[iTime];
    varArray[ 61] = PROPNN[iTime];
    varArray[ 62] = GLYX[iTime];
    varArray[ 63] = MAOPO2[iTime];
    varArray[ 64] = CH4[iTime];
    varArray[ 65] = GAOO[iTime];
    varArray[ 66] = B3O2[iTime];
    varArray[ 67] = ACET[iTime];
    varArray[ 68] = MACRN[iTime];
    varArray[ 69] = CH2OO[iTime];
    varArray[ 70] = MGLYOO[iTime];
    varArray[ 71] = VRO2[iTime];
    varArray[ 72] = MGLOO[iTime];
    varArray[ 73] = MACROO[iTime];
    varArray[ 74] = PO2[iTime];
    varArray[ 75] = CH3CHOO[iTime];
    varArray[ 76] = MAN2[iTime];
    varArray[ 77] = ISNOOA[iTime];
    varArray[ 78] = H2O2[iTime];
    varArray[ 79] = PRN1[iTime];
    varArray[ 80] = ETO2[iTime];
    varArray[ 81] = KO2[iTime];
    varArray[ 82] = RCO3[iTime];
    varArray[ 83] = HC5OO[iTime];
    varArray[ 84] = GLYC[iTime];
    varArray[ 85] = ClNO3[iTime];
    varArray[ 86] = RIO2[iTime];
    varArray[ 87] = R4N1[iTime];
    varArray[ 88] = HOCl[iTime];
    varArray[ 89] = ATO2[iTime];
    varArray[ 90] = HNO3[iTime];
    varArray[ 91] = ISN1[iTime];
    varArray[ 92] = MAO3[iTime];
    varArray[ 93] = MRO2[iTime];
    varArray[ 94] = INO2[iTime];
    varArray[ 95] = HAC[iTime];
    varArray[ 96] = HC5[iTime];
    varArray[ 97] = MGLY[iTime];
    varArray[ 98] = ISOPNBO2[iTime];
    varArray[ 99] = ISOPNDO2[iTime];
    varArray[100] = R4O2[iTime];
    varArray[101] = R4N2[iTime];
    varArray[102] = BrO[iTime];
    varArray[103] = RCHO[iTime];
    varArray[104] = MEK[iTime];
    varArray[105] = ClO[iTime];
    varArray[106] = MACR[iTime];
    varArray[107] = SO2[iTime];
    varArray[108] = MVK[iTime];
    varArray[109] = ALD2[iTime];
    varArray[110] = MCO3[iTime];
    varArray[111] = CH2O[iTime];
    varArray[112] = H2O[iTime];
    varArray[113] = Br[iTime];
    varArray[114] = NO[iTime];
    varArray[115] = NO3[iTime];
    varArray[116] = Cl[iTime];
    varArray[117] = O[iTime];
    varArray[118] = O1D[iTime];
    varArray[119] = O3[iTime];
    varArray[120] = HO2[iTime];
    varArray[121] = NO2[iTime];
    varArray[122] = OH[iTime];
    varArray[123] = HBr[iTime];
    varArray[124] = HCl[iTime];
    varArray[125] = CO[iTime];
    varArray[126] = MO2[iTime];
    fixArray[  0] = ACTA;
    fixArray[  1] = EOH;
    fixArray[  2] = H2;
    fixArray[  3] = HCOOH;
    fixArray[  4] = MOH;
    fixArray[  5] = N2;
    fixArray[  6] = O2;
    fixArray[  7] = RCOOH;
   
    aerArray[  0][0] = sootDens[iTime];
    aerArray[  0][1] = sootRadi[iTime];
    aerArray[  1][0] = iceDens[iTime];
    aerArray[  1][1] = iceRadi[iTime];
    aerArray[  2][0] = sulfDens[iTime];
    aerArray[  2][1] = sulfRadi[iTime];

    /* Ensure positiveness */
    for ( unsigned int i = 0; i < NVAR; i++ ) {
        if ( varArray[i] <= 0.0 ) {
            varArray[i] = 1.0E-50;
        }
    }

} /* End of Ambient::getData */

void Ambient::FillIn( double varArray[], unsigned int iTime )
{

    /* Ensure positiveness */
    for ( unsigned int i = 0; i < NVAR; i++ ) {
        if ( varArray[i] <= 0.0 ) {
            varArray[i] = 1.0E-50;
        }
    }
    
    CO2[iTime]      = varArray[  0];
    PPN[iTime]      = varArray[  1];
    BrNO2[iTime]    = varArray[  2];
    IEPOX[iTime]    = varArray[  3];
    PMNN[iTime]     = varArray[  4];
    N2O[iTime]      = varArray[  5];
    N[iTime]        = varArray[  6];
    PAN[iTime]      = varArray[  7];
    ALK4[iTime]     = varArray[  8];
    MAP[iTime]      = varArray[  9];
    MPN[iTime]      = varArray[ 10];
    Cl2O2[iTime]    = varArray[ 11];
    ETP[iTime]      = varArray[ 12];
    HNO2[iTime]     = varArray[ 13];
    C3H8[iTime]     = varArray[ 14];
    RA3P[iTime]     = varArray[ 15];
    RB3P[iTime]     = varArray[ 16];
    OClO[iTime]     = varArray[ 17];
    ClNO2[iTime]    = varArray[ 18];
    ISOP[iTime]     = varArray[ 19];
    HNO4[iTime]     = varArray[ 20];
    MAOP[iTime]     = varArray[ 21];
    MP[iTime]       = varArray[ 22];
    ClOO[iTime]     = varArray[ 23];
    RP[iTime]       = varArray[ 24];
    BrCl[iTime]     = varArray[ 25];
    PP[iTime]       = varArray[ 26];
    PRPN[iTime]     = varArray[ 27];
    SO4[iTime]      = varArray[ 28];
    Br2[iTime]      = varArray[ 29];
    ETHLN[iTime]    = varArray[ 30];
    MVKN[iTime]     = varArray[ 31];
    R4P[iTime]      = varArray[ 32];
    C2H6[iTime]     = varArray[ 33];
    RIP[iTime]      = varArray[ 34];
    VRP[iTime]      = varArray[ 35];
    ATOOH[iTime]    = varArray[ 36];
    IAP[iTime]      = varArray[ 37];
    DHMOB[iTime]    = varArray[ 38];
    MOBA[iTime]     = varArray[ 39];
    MRP[iTime]      = varArray[ 40];
    N2O5[iTime]     = varArray[ 41];
    ISNOHOO[iTime]  = varArray[ 42];
    ISNP[iTime]     = varArray[ 43];
    ISOPNB[iTime]   = varArray[ 44];
    IEPOXOO[iTime]  = varArray[ 45];
    MACRNO2[iTime]  = varArray[ 46];
    ROH[iTime]      = varArray[ 47];
    MOBAOO[iTime]   = varArray[ 48];
    DIBOO[iTime]    = varArray[ 49];
    PMN[iTime]      = varArray[ 50];
    ISNOOB[iTime]   = varArray[ 51];
    INPN[iTime]     = varArray[ 52];
    H[iTime]        = varArray[ 53];
    BrNO3[iTime]    = varArray[ 54];
    PRPE[iTime]     = varArray[ 55];
    MVKOO[iTime]    = varArray[ 56];
    Cl2[iTime]      = varArray[ 57];
    ISOPND[iTime]   = varArray[ 58];
    HOBr[iTime]     = varArray[ 59];
    A3O2[iTime]     = varArray[ 60];
    PROPNN[iTime]   = varArray[ 61];
    GLYX[iTime]     = varArray[ 62];
    MAOPO2[iTime]   = varArray[ 63];
    CH4[iTime]      = varArray[ 64];
    GAOO[iTime]     = varArray[ 65];
    B3O2[iTime]     = varArray[ 66];
    ACET[iTime]     = varArray[ 67];
    MACRN[iTime]    = varArray[ 68];
    CH2OO[iTime]    = varArray[ 69];
    MGLYOO[iTime]   = varArray[ 70];
    VRO2[iTime]     = varArray[ 71];
    MGLOO[iTime]    = varArray[ 72];
    MACROO[iTime]   = varArray[ 73];
    PO2[iTime]      = varArray[ 74];
    CH3CHOO[iTime]  = varArray[ 75];
    MAN2[iTime]     = varArray[ 76];
    ISNOOA[iTime]   = varArray[ 77];
    H2O2[iTime]     = varArray[ 78];
    PRN1[iTime]     = varArray[ 79];
    ETO2[iTime]     = varArray[ 80];
    KO2[iTime]      = varArray[ 81];
    RCO3[iTime]     = varArray[ 82];
    HC5OO[iTime]    = varArray[ 83];
    GLYC[iTime]     = varArray[ 84];
    ClNO3[iTime]    = varArray[ 85];
    RIO2[iTime]     = varArray[ 86];
    R4N1[iTime]     = varArray[ 87];
    HOCl[iTime]     = varArray[ 88];
    ATO2[iTime]     = varArray[ 89];
    HNO3[iTime]     = varArray[ 90];
    ISN1[iTime]     = varArray[ 91];
    MAO3[iTime]     = varArray[ 92];
    MRO2[iTime]     = varArray[ 93];
    INO2[iTime]     = varArray[ 94];
    HAC[iTime]      = varArray[ 95];
    HC5[iTime]      = varArray[ 96];
    MGLY[iTime]     = varArray[ 97];
    ISOPNBO2[iTime] = varArray[ 98];
    ISOPNDO2[iTime] = varArray[ 99];
    R4O2[iTime]     = varArray[100];
    R4N2[iTime]     = varArray[101];
    BrO[iTime]      = varArray[102];
    RCHO[iTime]     = varArray[103];
    MEK[iTime]      = varArray[104];
    ClO[iTime]      = varArray[105];
    MACR[iTime]     = varArray[106];
    SO2[iTime]      = varArray[107];
    MVK[iTime]      = varArray[108];
    ALD2[iTime]     = varArray[109];
    MCO3[iTime]     = varArray[110];
    CH2O[iTime]     = varArray[111];
    H2O[iTime]      = varArray[112];
    Br[iTime]       = varArray[113];
    NO[iTime]       = varArray[114];
    NO3[iTime]      = varArray[115];
    Cl[iTime]       = varArray[116];
    O[iTime]        = varArray[117];
    O1D[iTime]      = varArray[118];
    O3[iTime]       = varArray[119];
    HO2[iTime]      = varArray[120];
    NO2[iTime]      = varArray[121];
    OH[iTime]       = varArray[122];
    HBr[iTime]      = varArray[123];
    HCl[iTime]      = varArray[124];
    CO[iTime]       = varArray[125];
    MO2[iTime]      = varArray[126];
    
} /* End of Ambient::FillIn */

unsigned int Ambient::getnTime() const
{

    return nTime;

} /* End of Ambient::getnTime */


/* End of Ambient.cpp */
