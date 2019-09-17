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

Ambient::Ambient( UInt nTime_, Vector_1D ambientVector, \
                  Vector_2D aerVector, Vector_1D liqVector )
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

    SO4L.assign  ( nTime, liqVector[0] );
    H2OL.assign  ( nTime, liqVector[1] );
    HNO3L.assign ( nTime, liqVector[2] );
    HClL.assign  ( nTime, liqVector[3] );
    HOClL.assign ( nTime, liqVector[4] );
    HBrL.assign  ( nTime, liqVector[5] );
    HOBrL.assign ( nTime, liqVector[6] );
    H2OS.assign  ( nTime, liqVector[7] );
    HNO3S.assign ( nTime, liqVector[8] );
   
    NIT.assign   ( nTime, 0.0 );

    SO4T.assign  ( nTime, ambientVector[ 28] + liqVector[0] );
    
    sootDens.assign( nTime, aerVector[  0][0] );
    sootRadi.assign( nTime, aerVector[  0][1] );
    sootArea.assign( nTime, aerVector[  0][2] );
    iceDens.assign ( nTime, aerVector[  1][0] );
    iceRadi.assign ( nTime, aerVector[  1][1] );
    iceArea.assign ( nTime, aerVector[  1][2] );
    sulfDens.assign( nTime, aerVector[  2][0] );
    sulfRadi.assign( nTime, aerVector[  2][1] );
    sulfArea.assign( nTime, aerVector[  1][2] );

    cosSZA.assign( nTime - 1, 0.0 );

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

    SO4L     = a.SO4L     ;
    H2OL     = a.H2OL     ;
    HNO3L    = a.HNO3L    ;
    HClL     = a.HClL     ;
    HOClL    = a.HOClL    ;
    HBrL     = a.HBrL     ;
    HOBrL    = a.HOBrL    ;
    H2OS     = a.H2OS     ;
    HNO3S    = a.HNO3S    ;
    NIT      = a.NIT      ;
    SO4T     = a.SO4T     ;

    sootDens = a.sootDens ;
    sootRadi = a.sootRadi ;
    iceDens  = a.iceDens  ;
    sulfDens = a.sulfDens ;
    sulfRadi = a.sulfRadi ;

    cosSZA = a.cosSZA;

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
    
    SO4L     = a.SO4L     ;
    H2OL     = a.H2OL     ;
    HNO3L    = a.HNO3L    ;
    HClL     = a.HClL     ;
    HOClL    = a.HOClL    ;
    HBrL     = a.HBrL     ;
    HOBrL    = a.HOBrL    ;
    H2OS     = a.H2OS     ;
    HNO3S    = a.HNO3S    ;
    NIT      = a.NIT      ;
    SO4T     = a.SO4T     ;

    sootDens = a.sootDens ;
    sootRadi = a.sootRadi ;
    iceDens  = a.iceDens  ;
    sulfDens = a.sulfDens ;
    sulfRadi = a.sulfRadi ;

    cosSZA = a.cosSZA;

    return *this;

} /* End of Ambient::operator= */

void Ambient::getData( RealDouble aerArray[][2], UInt iTime ) const
{

    VAR[  0] = CO2[iTime];
    VAR[  1] = PPN[iTime];
    VAR[  2] = BrNO2[iTime];
    VAR[  3] = IEPOX[iTime];
    VAR[  4] = PMNN[iTime];
    VAR[  5] = N2O[iTime];
    VAR[  6] = N[iTime];
    VAR[  7] = PAN[iTime];
    VAR[  8] = ALK4[iTime];
    VAR[  9] = MAP[iTime];
    VAR[ 10] = MPN[iTime];
    VAR[ 11] = Cl2O2[iTime];
    VAR[ 12] = ETP[iTime];
    VAR[ 13] = HNO2[iTime];
    VAR[ 14] = C3H8[iTime];
    VAR[ 15] = RA3P[iTime];
    VAR[ 16] = RB3P[iTime];
    VAR[ 17] = OClO[iTime];
    VAR[ 18] = ClNO2[iTime];
    VAR[ 19] = ISOP[iTime];
    VAR[ 20] = HNO4[iTime];
    VAR[ 21] = MAOP[iTime];
    VAR[ 22] = MP[iTime];
    VAR[ 23] = ClOO[iTime];
    VAR[ 24] = RP[iTime];
    VAR[ 25] = BrCl[iTime];
    VAR[ 26] = PP[iTime];
    VAR[ 27] = PRPN[iTime];
    VAR[ 28] = SO4[iTime];
    VAR[ 29] = Br2[iTime];
    VAR[ 30] = ETHLN[iTime];
    VAR[ 31] = MVKN[iTime];
    VAR[ 32] = R4P[iTime];
    VAR[ 33] = C2H6[iTime];
    VAR[ 34] = RIP[iTime];
    VAR[ 35] = VRP[iTime];
    VAR[ 36] = ATOOH[iTime];
    VAR[ 37] = IAP[iTime];
    VAR[ 38] = DHMOB[iTime];
    VAR[ 39] = MOBA[iTime];
    VAR[ 40] = MRP[iTime];
    VAR[ 41] = N2O5[iTime];
    VAR[ 42] = ISNOHOO[iTime];
    VAR[ 43] = ISNP[iTime];
    VAR[ 44] = ISOPNB[iTime];
    VAR[ 45] = IEPOXOO[iTime];
    VAR[ 46] = MACRNO2[iTime];
    VAR[ 47] = ROH[iTime];
    VAR[ 48] = MOBAOO[iTime];
    VAR[ 49] = DIBOO[iTime];
    VAR[ 50] = PMN[iTime];
    VAR[ 51] = ISNOOB[iTime];
    VAR[ 52] = INPN[iTime];
    VAR[ 53] = H[iTime];
    VAR[ 54] = BrNO3[iTime];
    VAR[ 55] = PRPE[iTime];
    VAR[ 56] = MVKOO[iTime];
    VAR[ 57] = Cl2[iTime];
    VAR[ 58] = ISOPND[iTime];
    VAR[ 59] = HOBr[iTime];
    VAR[ 60] = A3O2[iTime];
    VAR[ 61] = PROPNN[iTime];
    VAR[ 62] = GLYX[iTime];
    VAR[ 63] = MAOPO2[iTime];
    VAR[ 64] = CH4[iTime];
    VAR[ 65] = GAOO[iTime];
    VAR[ 66] = B3O2[iTime];
    VAR[ 67] = ACET[iTime];
    VAR[ 68] = MACRN[iTime];
    VAR[ 69] = CH2OO[iTime];
    VAR[ 70] = MGLYOO[iTime];
    VAR[ 71] = VRO2[iTime];
    VAR[ 72] = MGLOO[iTime];
    VAR[ 73] = MACROO[iTime];
    VAR[ 74] = PO2[iTime];
    VAR[ 75] = CH3CHOO[iTime];
    VAR[ 76] = MAN2[iTime];
    VAR[ 77] = ISNOOA[iTime];
    VAR[ 78] = H2O2[iTime];
    VAR[ 79] = PRN1[iTime];
    VAR[ 80] = ETO2[iTime];
    VAR[ 81] = KO2[iTime];
    VAR[ 82] = RCO3[iTime];
    VAR[ 83] = HC5OO[iTime];
    VAR[ 84] = GLYC[iTime];
    VAR[ 85] = ClNO3[iTime];
    VAR[ 86] = RIO2[iTime];
    VAR[ 87] = R4N1[iTime];
    VAR[ 88] = HOCl[iTime];
    VAR[ 89] = ATO2[iTime];
    VAR[ 90] = HNO3[iTime];
    VAR[ 91] = ISN1[iTime];
    VAR[ 92] = MAO3[iTime];
    VAR[ 93] = MRO2[iTime];
    VAR[ 94] = INO2[iTime];
    VAR[ 95] = HAC[iTime];
    VAR[ 96] = HC5[iTime];
    VAR[ 97] = MGLY[iTime];
    VAR[ 98] = ISOPNBO2[iTime];
    VAR[ 99] = ISOPNDO2[iTime];
    VAR[100] = R4O2[iTime];
    VAR[101] = R4N2[iTime];
    VAR[102] = BrO[iTime];
    VAR[103] = RCHO[iTime];
    VAR[104] = MEK[iTime];
    VAR[105] = ClO[iTime];
    VAR[106] = MACR[iTime];
    VAR[107] = SO2[iTime];
    VAR[108] = MVK[iTime];
    VAR[109] = ALD2[iTime];
    VAR[110] = MCO3[iTime];
    VAR[111] = CH2O[iTime];
    VAR[112] = H2O[iTime];
    VAR[113] = Br[iTime];
    VAR[114] = NO[iTime];
    VAR[115] = NO3[iTime];
    VAR[116] = Cl[iTime];
    VAR[117] = O[iTime];
    VAR[118] = O1D[iTime];
    VAR[119] = O3[iTime];
    VAR[120] = HO2[iTime];
    VAR[121] = NO2[iTime];
    VAR[122] = OH[iTime];
    VAR[123] = HBr[iTime];
    VAR[124] = HCl[iTime];
    VAR[125] = CO[iTime];
    VAR[126] = MO2[iTime];
    FIX[  0] = ACTA;
    FIX[  1] = EOH;
    FIX[  2] = H2;
    FIX[  3] = HCOOH;
    FIX[  4] = MOH;
    FIX[  5] = N2;
    FIX[  6] = O2;
    FIX[  7] = RCOOH;
   
    aerArray[  0][0] = sootDens[iTime];
    aerArray[  0][1] = sootRadi[iTime];
    aerArray[  1][0] = iceDens[iTime];
    aerArray[  1][1] = iceRadi[iTime];
    aerArray[  2][0] = sulfDens[iTime];
    aerArray[  2][1] = sulfRadi[iTime];

    /* Ensure positiveness */
    for ( UInt i = 0; i < NVAR; i++ ) {
        if ( VAR[i] <= 0.0 ) {
            VAR[i] = 1.0E-50;
        }
    }

} /* End of Ambient::getData */

void Ambient::FillIn( UInt iTime )
{

    /* Ensure positiveness */
    for ( UInt i = 0; i < NVAR; i++ ) {
        if ( VAR[i] <= 0.0 ) {
            VAR[i] = 1.0E-50;
        }
    }
    
    CO2[iTime]      = VAR[  0];
    PPN[iTime]      = VAR[  1];
    BrNO2[iTime]    = VAR[  2];
    IEPOX[iTime]    = VAR[  3];
    PMNN[iTime]     = VAR[  4];
    N2O[iTime]      = VAR[  5];
    N[iTime]        = VAR[  6];
    PAN[iTime]      = VAR[  7];
    ALK4[iTime]     = VAR[  8];
    MAP[iTime]      = VAR[  9];
    MPN[iTime]      = VAR[ 10];
    Cl2O2[iTime]    = VAR[ 11];
    ETP[iTime]      = VAR[ 12];
    HNO2[iTime]     = VAR[ 13];
    C3H8[iTime]     = VAR[ 14];
    RA3P[iTime]     = VAR[ 15];
    RB3P[iTime]     = VAR[ 16];
    OClO[iTime]     = VAR[ 17];
    ClNO2[iTime]    = VAR[ 18];
    ISOP[iTime]     = VAR[ 19];
    HNO4[iTime]     = VAR[ 20];
    MAOP[iTime]     = VAR[ 21];
    MP[iTime]       = VAR[ 22];
    ClOO[iTime]     = VAR[ 23];
    RP[iTime]       = VAR[ 24];
    BrCl[iTime]     = VAR[ 25];
    PP[iTime]       = VAR[ 26];
    PRPN[iTime]     = VAR[ 27];
    SO4[iTime]      = VAR[ 28];
    Br2[iTime]      = VAR[ 29];
    ETHLN[iTime]    = VAR[ 30];
    MVKN[iTime]     = VAR[ 31];
    R4P[iTime]      = VAR[ 32];
    C2H6[iTime]     = VAR[ 33];
    RIP[iTime]      = VAR[ 34];
    VRP[iTime]      = VAR[ 35];
    ATOOH[iTime]    = VAR[ 36];
    IAP[iTime]      = VAR[ 37];
    DHMOB[iTime]    = VAR[ 38];
    MOBA[iTime]     = VAR[ 39];
    MRP[iTime]      = VAR[ 40];
    N2O5[iTime]     = VAR[ 41];
    ISNOHOO[iTime]  = VAR[ 42];
    ISNP[iTime]     = VAR[ 43];
    ISOPNB[iTime]   = VAR[ 44];
    IEPOXOO[iTime]  = VAR[ 45];
    MACRNO2[iTime]  = VAR[ 46];
    ROH[iTime]      = VAR[ 47];
    MOBAOO[iTime]   = VAR[ 48];
    DIBOO[iTime]    = VAR[ 49];
    PMN[iTime]      = VAR[ 50];
    ISNOOB[iTime]   = VAR[ 51];
    INPN[iTime]     = VAR[ 52];
    H[iTime]        = VAR[ 53];
    BrNO3[iTime]    = VAR[ 54];
    PRPE[iTime]     = VAR[ 55];
    MVKOO[iTime]    = VAR[ 56];
    Cl2[iTime]      = VAR[ 57];
    ISOPND[iTime]   = VAR[ 58];
    HOBr[iTime]     = VAR[ 59];
    A3O2[iTime]     = VAR[ 60];
    PROPNN[iTime]   = VAR[ 61];
    GLYX[iTime]     = VAR[ 62];
    MAOPO2[iTime]   = VAR[ 63];
    CH4[iTime]      = VAR[ 64];
    GAOO[iTime]     = VAR[ 65];
    B3O2[iTime]     = VAR[ 66];
    ACET[iTime]     = VAR[ 67];
    MACRN[iTime]    = VAR[ 68];
    CH2OO[iTime]    = VAR[ 69];
    MGLYOO[iTime]   = VAR[ 70];
    VRO2[iTime]     = VAR[ 71];
    MGLOO[iTime]    = VAR[ 72];
    MACROO[iTime]   = VAR[ 73];
    PO2[iTime]      = VAR[ 74];
    CH3CHOO[iTime]  = VAR[ 75];
    MAN2[iTime]     = VAR[ 76];
    ISNOOA[iTime]   = VAR[ 77];
    H2O2[iTime]     = VAR[ 78];
    PRN1[iTime]     = VAR[ 79];
    ETO2[iTime]     = VAR[ 80];
    KO2[iTime]      = VAR[ 81];
    RCO3[iTime]     = VAR[ 82];
    HC5OO[iTime]    = VAR[ 83];
    GLYC[iTime]     = VAR[ 84];
    ClNO3[iTime]    = VAR[ 85];
    RIO2[iTime]     = VAR[ 86];
    R4N1[iTime]     = VAR[ 87];
    HOCl[iTime]     = VAR[ 88];
    ATO2[iTime]     = VAR[ 89];
    HNO3[iTime]     = VAR[ 90];
    ISN1[iTime]     = VAR[ 91];
    MAO3[iTime]     = VAR[ 92];
    MRO2[iTime]     = VAR[ 93];
    INO2[iTime]     = VAR[ 94];
    HAC[iTime]      = VAR[ 95];
    HC5[iTime]      = VAR[ 96];
    MGLY[iTime]     = VAR[ 97];
    ISOPNBO2[iTime] = VAR[ 98];
    ISOPNDO2[iTime] = VAR[ 99];
    R4O2[iTime]     = VAR[100];
    R4N2[iTime]     = VAR[101];
    BrO[iTime]      = VAR[102];
    RCHO[iTime]     = VAR[103];
    MEK[iTime]      = VAR[104];
    ClO[iTime]      = VAR[105];
    MACR[iTime]     = VAR[106];
    SO2[iTime]      = VAR[107];
    MVK[iTime]      = VAR[108];
    ALD2[iTime]     = VAR[109];
    MCO3[iTime]     = VAR[110];
    CH2O[iTime]     = VAR[111];
    H2O[iTime]      = VAR[112];
    Br[iTime]       = VAR[113];
    NO[iTime]       = VAR[114];
    NO3[iTime]      = VAR[115];
    Cl[iTime]       = VAR[116];
    O[iTime]        = VAR[117];
    O1D[iTime]      = VAR[118];
    O3[iTime]       = VAR[119];
    HO2[iTime]      = VAR[120];
    NO2[iTime]      = VAR[121];
    OH[iTime]       = VAR[122];
    HBr[iTime]      = VAR[123];
    HCl[iTime]      = VAR[124];
    CO[iTime]       = VAR[125];
    MO2[iTime]      = VAR[126];
    
} /* End of Ambient::FillIn */

UInt Ambient::getnTime() const
{

    return nTime;

} /* End of Ambient::getnTime */


/* End of Ambient.cpp */
