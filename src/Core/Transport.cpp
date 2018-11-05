/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Transport Program File                                           */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : Transport.cpp                             */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <iostream>

#include "Core/Structure.hpp"
#include "SANDS/Solver.hpp"

void Transport( Solution& data, Solver& GasPhase )
{
        
    const bool isInputReal = 1;

    /* Gaseous species */
    GasPhase.Solve( data.CO2      , isInputReal );
    GasPhase.Solve( data.PPN      , isInputReal );
    GasPhase.Solve( data.BrNO2    , isInputReal );
    GasPhase.Solve( data.IEPOX    , isInputReal );
    GasPhase.Solve( data.PMNN     , isInputReal );
    GasPhase.Solve( data.N2O      , isInputReal );
    GasPhase.Solve( data.N        , isInputReal );
    GasPhase.Solve( data.PAN      , isInputReal );
    GasPhase.Solve( data.ALK4     , isInputReal );
    GasPhase.Solve( data.MAP      , isInputReal );
    GasPhase.Solve( data.MPN      , isInputReal );
    GasPhase.Solve( data.Cl2O2    , isInputReal );
    GasPhase.Solve( data.ETP      , isInputReal );
    GasPhase.Solve( data.HNO2     , isInputReal );
    GasPhase.Solve( data.C3H8     , isInputReal );
    GasPhase.Solve( data.RA3P     , isInputReal );
    GasPhase.Solve( data.RB3P     , isInputReal );
    GasPhase.Solve( data.OClO     , isInputReal );
    GasPhase.Solve( data.ClNO2    , isInputReal );
    GasPhase.Solve( data.ISOP     , isInputReal );
    GasPhase.Solve( data.HNO4     , isInputReal );
    GasPhase.Solve( data.MAOP     , isInputReal );
    GasPhase.Solve( data.MP       , isInputReal );
    GasPhase.Solve( data.ClOO     , isInputReal );
    GasPhase.Solve( data.RP       , isInputReal );
    GasPhase.Solve( data.BrCl     , isInputReal );
    GasPhase.Solve( data.PP       , isInputReal );
    GasPhase.Solve( data.PRPN     , isInputReal );
    GasPhase.Solve( data.SO4      , isInputReal );
    GasPhase.Solve( data.Br2      , isInputReal );
    GasPhase.Solve( data.ETHLN    , isInputReal );
    GasPhase.Solve( data.MVKN     , isInputReal );
    GasPhase.Solve( data.R4P      , isInputReal );
    GasPhase.Solve( data.C2H6     , isInputReal );
    GasPhase.Solve( data.RIP      , isInputReal );
    GasPhase.Solve( data.VRP      , isInputReal );
    GasPhase.Solve( data.ATOOH    , isInputReal );
    GasPhase.Solve( data.IAP      , isInputReal );
    GasPhase.Solve( data.DHMOB    , isInputReal );
    GasPhase.Solve( data.MOBA     , isInputReal );
    GasPhase.Solve( data.MRP      , isInputReal );
    GasPhase.Solve( data.N2O5     , isInputReal );
    GasPhase.Solve( data.ISNOHOO  , isInputReal );
    GasPhase.Solve( data.ISNP     , isInputReal );
    GasPhase.Solve( data.ISOPNB   , isInputReal );
    GasPhase.Solve( data.IEPOXOO  , isInputReal );
    GasPhase.Solve( data.MACRNO2  , isInputReal );
    GasPhase.Solve( data.ROH      , isInputReal );
    GasPhase.Solve( data.MOBAOO   , isInputReal );
    GasPhase.Solve( data.DIBOO    , isInputReal );
    GasPhase.Solve( data.PMN      , isInputReal );
    GasPhase.Solve( data.ISNOOB   , isInputReal );
    GasPhase.Solve( data.INPN     , isInputReal );
    GasPhase.Solve( data.H        , isInputReal );
    GasPhase.Solve( data.BrNO3    , isInputReal );
    GasPhase.Solve( data.PRPE     , isInputReal );
    GasPhase.Solve( data.MVKOO    , isInputReal );
    GasPhase.Solve( data.Cl2      , isInputReal );
    GasPhase.Solve( data.ISOPND   , isInputReal );
    GasPhase.Solve( data.HOBr     , isInputReal );
    GasPhase.Solve( data.A3O2     , isInputReal );
    GasPhase.Solve( data.PROPNN   , isInputReal );
    GasPhase.Solve( data.GLYX     , isInputReal );
    GasPhase.Solve( data.MAOPO2   , isInputReal );
    GasPhase.Solve( data.CH4      , isInputReal );
    GasPhase.Solve( data.GAOO     , isInputReal );
    GasPhase.Solve( data.B3O2     , isInputReal );
    GasPhase.Solve( data.ACET     , isInputReal );
    GasPhase.Solve( data.MACRN    , isInputReal );
    GasPhase.Solve( data.CH2OO    , isInputReal );
    GasPhase.Solve( data.MGLYOO   , isInputReal );
    GasPhase.Solve( data.VRO2     , isInputReal );
    GasPhase.Solve( data.MGLOO    , isInputReal );
    GasPhase.Solve( data.MACROO   , isInputReal );
    GasPhase.Solve( data.PO2      , isInputReal );
    GasPhase.Solve( data.CH3CHOO  , isInputReal );
    GasPhase.Solve( data.MAN2     , isInputReal );
    GasPhase.Solve( data.ISNOOA   , isInputReal );
    GasPhase.Solve( data.H2O2     , isInputReal );
    GasPhase.Solve( data.PRN1     , isInputReal );
    GasPhase.Solve( data.ETO2     , isInputReal );
    GasPhase.Solve( data.KO2      , isInputReal );
    GasPhase.Solve( data.RCO3     , isInputReal );
    GasPhase.Solve( data.HC5OO    , isInputReal );
    GasPhase.Solve( data.GLYC     , isInputReal );
    GasPhase.Solve( data.ClNO3    , isInputReal );
    GasPhase.Solve( data.RIO2     , isInputReal );
    GasPhase.Solve( data.R4N1     , isInputReal );
    GasPhase.Solve( data.HOCl     , isInputReal );
    GasPhase.Solve( data.ATO2     , isInputReal );
    GasPhase.Solve( data.HNO3     , isInputReal );
    GasPhase.Solve( data.ISN1     , isInputReal );
    GasPhase.Solve( data.MAO3     , isInputReal );
    GasPhase.Solve( data.MRO2     , isInputReal );
    GasPhase.Solve( data.INO2     , isInputReal );
    GasPhase.Solve( data.HAC      , isInputReal );
    GasPhase.Solve( data.HC5      , isInputReal );
    GasPhase.Solve( data.MGLY     , isInputReal );
    GasPhase.Solve( data.ISOPNBO2 , isInputReal );
    GasPhase.Solve( data.ISOPNDO2 , isInputReal );
    GasPhase.Solve( data.R4O2     , isInputReal );
    GasPhase.Solve( data.R4N2     , isInputReal );
    GasPhase.Solve( data.BrO      , isInputReal );
    GasPhase.Solve( data.RCHO     , isInputReal );
    GasPhase.Solve( data.MEK      , isInputReal );
    GasPhase.Solve( data.ClO      , isInputReal );
    GasPhase.Solve( data.MACR     , isInputReal );
    GasPhase.Solve( data.SO2      , isInputReal );
    GasPhase.Solve( data.MVK      , isInputReal );
    GasPhase.Solve( data.ALD2     , isInputReal );
    GasPhase.Solve( data.MCO3     , isInputReal );
    GasPhase.Solve( data.CH2O     , isInputReal );
    GasPhase.Solve( data.H2O      , isInputReal );
    GasPhase.Solve( data.Br       , isInputReal );
    GasPhase.Solve( data.NO       , isInputReal );
    GasPhase.Solve( data.NO3      , isInputReal );
    GasPhase.Solve( data.Cl       , isInputReal );
    GasPhase.Solve( data.O        , isInputReal );
    GasPhase.Solve( data.O1D      , isInputReal );
    GasPhase.Solve( data.O3       , isInputReal );
    GasPhase.Solve( data.HO2      , isInputReal );
    GasPhase.Solve( data.NO2      , isInputReal );
    GasPhase.Solve( data.OH       , isInputReal );
    GasPhase.Solve( data.HBr      , isInputReal );
    GasPhase.Solve( data.HCl      , isInputReal );
    GasPhase.Solve( data.CO       , isInputReal );
    GasPhase.Solve( data.MO2      , isInputReal );

} /* End of CallSolver */

