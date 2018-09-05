/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Struct_SetShape Program File                                     */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : Struct_SetShape.cpp                       */
/* Working directory    : /home/fritzt/APCEMM-SourceCode            */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "Structure.h"
#include "Parameters.h"

void Struct_SetShape( Solution &sol )
{
    /* Gaseous species */
    sol.SetShape( sol.CO2, NY, NX );
    sol.SetShape( sol.PPN, NY, NX );
    sol.SetShape( sol.BrNO2, NY, NX );
    sol.SetShape( sol.IEPOX, NY, NX );
    sol.SetShape( sol.PMNN, NY, NX );
    sol.SetShape( sol.N2O, NY, NX );
    sol.SetShape( sol.N, NY, NX );
    sol.SetShape( sol.PAN, NY, NX );
    sol.SetShape( sol.ALK4, NY, NX );
    sol.SetShape( sol.MAP, NY, NX );
    sol.SetShape( sol.MPN, NY, NX );
    sol.SetShape( sol.Cl2O2, NY, NX );
    sol.SetShape( sol.ETP, NY, NX );
    sol.SetShape( sol.HNO2, NY, NX );
    sol.SetShape( sol.C3H8, NY, NX );
    sol.SetShape( sol.RA3P, NY, NX );
    sol.SetShape( sol.RB3P, NY, NX );
    sol.SetShape( sol.OClO, NY, NX );
    sol.SetShape( sol.ClNO2, NY, NX );
    sol.SetShape( sol.ISOP, NY, NX );
    sol.SetShape( sol.HNO4, NY, NX );
    sol.SetShape( sol.MAOP, NY, NX );
    sol.SetShape( sol.MP, NY, NX );
    sol.SetShape( sol.ClOO, NY, NX );
    sol.SetShape( sol.RP, NY, NX );
    sol.SetShape( sol.BrCl, NY, NX );
    sol.SetShape( sol.PP, NY, NX );
    sol.SetShape( sol.PRPN, NY, NX );
    sol.SetShape( sol.SO4, NY, NX );
    sol.SetShape( sol.Br2, NY, NX );
    sol.SetShape( sol.ETHLN, NY, NX );
    sol.SetShape( sol.MVKN, NY, NX );
    sol.SetShape( sol.R4P, NY, NX );
    sol.SetShape( sol.C2H6, NY, NX );
    sol.SetShape( sol.RIP, NY, NX );
    sol.SetShape( sol.VRP, NY, NX );
    sol.SetShape( sol.ATOOH, NY, NX );
    sol.SetShape( sol.IAP, NY, NX );
    sol.SetShape( sol.DHMOB, NY, NX );
    sol.SetShape( sol.MOBA, NY, NX );
    sol.SetShape( sol.MRP, NY, NX );
    sol.SetShape( sol.N2O5, NY, NX );
    sol.SetShape( sol.ISNOHOO, NY, NX );
    sol.SetShape( sol.ISNP, NY, NX );
    sol.SetShape( sol.ISOPNB, NY, NX );
    sol.SetShape( sol.IEPOXOO, NY, NX );
    sol.SetShape( sol.MACRNO2, NY, NX );
    sol.SetShape( sol.ROH, NY, NX );
    sol.SetShape( sol.MOBAOO, NY, NX );
    sol.SetShape( sol.DIBOO, NY, NX );
    sol.SetShape( sol.PMN, NY, NX );
    sol.SetShape( sol.ISNOOB, NY, NX );
    sol.SetShape( sol.INPN, NY, NX );
    sol.SetShape( sol.H, NY, NX );
    sol.SetShape( sol.BrNO3, NY, NX );
    sol.SetShape( sol.PRPE, NY, NX );
    sol.SetShape( sol.MVKOO, NY, NX );
    sol.SetShape( sol.Cl2, NY, NX );
    sol.SetShape( sol.ISOPND, NY, NX );
    sol.SetShape( sol.HOBr, NY, NX );
    sol.SetShape( sol.A3O2, NY, NX );
    sol.SetShape( sol.PROPNN, NY, NX );
    sol.SetShape( sol.GLYX, NY, NX );
    sol.SetShape( sol.MAOPO2, NY, NX );
    sol.SetShape( sol.CH4, NY, NX );
    sol.SetShape( sol.GAOO, NY, NX );
    sol.SetShape( sol.B3O2, NY, NX );
    sol.SetShape( sol.ACET, NY, NX );
    sol.SetShape( sol.MACRN, NY, NX );
    sol.SetShape( sol.CH2OO, NY, NX );
    sol.SetShape( sol.MGLYOO, NY, NX );
    sol.SetShape( sol.VRO2, NY, NX );
    sol.SetShape( sol.MGLOO, NY, NX );
    sol.SetShape( sol.MACROO, NY, NX );
    sol.SetShape( sol.PO2, NY, NX );
    sol.SetShape( sol.CH3CHOO, NY, NX );
    sol.SetShape( sol.MAN2, NY, NX );
    sol.SetShape( sol.ISNOOA, NY, NX );
    sol.SetShape( sol.H2O2, NY, NX );
    sol.SetShape( sol.PRN1, NY, NX );
    sol.SetShape( sol.ETO2, NY, NX );
    sol.SetShape( sol.KO2, NY, NX );
    sol.SetShape( sol.RCO3, NY, NX );
    sol.SetShape( sol.HC5OO, NY, NX );
    sol.SetShape( sol.GLYC, NY, NX );
    sol.SetShape( sol.ClNO3, NY, NX );
    sol.SetShape( sol.RIO2, NY, NX );
    sol.SetShape( sol.R4N1, NY, NX );
    sol.SetShape( sol.HOCl, NY, NX );
    sol.SetShape( sol.ATO2, NY, NX );
    sol.SetShape( sol.HNO3, NY, NX );
    sol.SetShape( sol.ISN1, NY, NX );
    sol.SetShape( sol.MAO3, NY, NX );
    sol.SetShape( sol.MRO2, NY, NX );
    sol.SetShape( sol.INO2, NY, NX );
    sol.SetShape( sol.HAC, NY, NX );
    sol.SetShape( sol.HC5, NY, NX );
    sol.SetShape( sol.MGLY, NY, NX );
    sol.SetShape( sol.ISOPNBO2, NY, NX );
    sol.SetShape( sol.ISOPNDO2, NY, NX );
    sol.SetShape( sol.R4O2, NY, NX );
    sol.SetShape( sol.R4N2, NY, NX );
    sol.SetShape( sol.BrO, NY, NX );
    sol.SetShape( sol.RCHO, NY, NX );
    sol.SetShape( sol.MEK, NY, NX );
    sol.SetShape( sol.ClO, NY, NX );
    sol.SetShape( sol.MACR, NY, NX );
    sol.SetShape( sol.SO2, NY, NX );
    sol.SetShape( sol.MVK, NY, NX );
    sol.SetShape( sol.ALD2, NY, NX );
    sol.SetShape( sol.MCO3, NY, NX );
    sol.SetShape( sol.CH2O, NY, NX );
    sol.SetShape( sol.H2O, NY, NX );
    sol.SetShape( sol.Br, NY, NX );
    sol.SetShape( sol.NO, NY, NX );
    sol.SetShape( sol.NO3, NY, NX );
    sol.SetShape( sol.Cl, NY, NX );
    sol.SetShape( sol.O, NY, NX );
    sol.SetShape( sol.O1D, NY, NX );
    sol.SetShape( sol.O3, NY, NX );
    sol.SetShape( sol.HO2, NY, NX );
    sol.SetShape( sol.NO2, NY, NX );
    sol.SetShape( sol.OH, NY, NX );
    sol.SetShape( sol.HBr, NY, NX );
    sol.SetShape( sol.HCl, NY, NX );
    sol.SetShape( sol.CO, NY, NX );
    sol.SetShape( sol.MO2, NY, NX );
    sol.SetShape( sol.ACTA, NY, NX );
    sol.SetShape( sol.EOH, NY, NX );
    sol.SetShape( sol.H2, NY, NX );
    sol.SetShape( sol.HCOOH, NY, NX );
    sol.SetShape( sol.MOH, NY, NX );
    sol.SetShape( sol.N2, NY, NX );
    sol.SetShape( sol.O2, NY, NX );
    sol.SetShape( sol.RCOOH, NY, NX );

    /* Liquid/solid species */
    sol.SetShape( sol.SO4L, NY, NX );
    sol.SetShape( sol.H2OL, NY, NX );
    sol.SetShape( sol.H2OS, NY, NX );
    sol.SetShape( sol.HNO3L, NY, NX );
    sol.SetShape( sol.HNO3S, NY, NX );
    sol.SetShape( sol.HClL, NY, NX );
    sol.SetShape( sol.HOClL, NY, NX );
    sol.SetShape( sol.HBrL, NY, NX );
    sol.SetShape( sol.HOBrL, NY, NX );

    /* Aerosols */
    sol.SetShape( sol.Soot, NY, NX );
    sol.SetShape( sol.SLA, NY, NX );
    sol.SetShape( sol.SPA, NY, NX );

}
