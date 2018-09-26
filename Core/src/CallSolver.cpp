/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* CallSolver Program File                                          */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : CallSolver.cpp                            */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <iostream>

#include "Structure.hpp"
#include "../../SANDS/include/Solver.hpp"

void CallSolver( Solution& data, Solver& solv )
{
        
    const bool isInputReal = 1;

    /* Gaseous species */
    solv.Solve( data.CO2      , isInputReal );
    solv.Solve( data.PPN      , isInputReal );
    solv.Solve( data.BrNO2    , isInputReal );
    solv.Solve( data.IEPOX    , isInputReal );
    solv.Solve( data.PMNN     , isInputReal );
    solv.Solve( data.N2O      , isInputReal );
    solv.Solve( data.N        , isInputReal );
    solv.Solve( data.PAN      , isInputReal );
    solv.Solve( data.ALK4     , isInputReal );
    solv.Solve( data.MAP      , isInputReal );
    solv.Solve( data.MPN      , isInputReal );
    solv.Solve( data.Cl2O2    , isInputReal );
    solv.Solve( data.ETP      , isInputReal );
    solv.Solve( data.HNO2     , isInputReal );
    solv.Solve( data.C3H8     , isInputReal );
    solv.Solve( data.RA3P     , isInputReal );
    solv.Solve( data.RB3P     , isInputReal );
    solv.Solve( data.OClO     , isInputReal );
    solv.Solve( data.ClNO2    , isInputReal );
    solv.Solve( data.ISOP     , isInputReal );
    solv.Solve( data.HNO4     , isInputReal );
    solv.Solve( data.MAOP     , isInputReal );
    solv.Solve( data.MP       , isInputReal );
    solv.Solve( data.ClOO     , isInputReal );
    solv.Solve( data.RP       , isInputReal );
    solv.Solve( data.BrCl     , isInputReal );
    solv.Solve( data.PP       , isInputReal );
    solv.Solve( data.PRPN     , isInputReal );
    solv.Solve( data.SO4      , isInputReal );
    solv.Solve( data.Br2      , isInputReal );
    solv.Solve( data.ETHLN    , isInputReal );
    solv.Solve( data.MVKN     , isInputReal );
    solv.Solve( data.R4P      , isInputReal );
    solv.Solve( data.C2H6     , isInputReal );
    solv.Solve( data.RIP      , isInputReal );
    solv.Solve( data.VRP      , isInputReal );
    solv.Solve( data.ATOOH    , isInputReal );
    solv.Solve( data.IAP      , isInputReal );
    solv.Solve( data.DHMOB    , isInputReal );
    solv.Solve( data.MOBA     , isInputReal );
    solv.Solve( data.MRP      , isInputReal );
    solv.Solve( data.N2O5     , isInputReal );
    solv.Solve( data.ISNOHOO  , isInputReal );
    solv.Solve( data.ISNP     , isInputReal );
    solv.Solve( data.ISOPNB   , isInputReal );
    solv.Solve( data.IEPOXOO  , isInputReal );
    solv.Solve( data.MACRNO2  , isInputReal );
    solv.Solve( data.ROH      , isInputReal );
    solv.Solve( data.MOBAOO   , isInputReal );
    solv.Solve( data.DIBOO    , isInputReal );
    solv.Solve( data.PMN      , isInputReal );
    solv.Solve( data.ISNOOB   , isInputReal );
    solv.Solve( data.INPN     , isInputReal );
    solv.Solve( data.H        , isInputReal );
    solv.Solve( data.BrNO3    , isInputReal );
    solv.Solve( data.PRPE     , isInputReal );
    solv.Solve( data.MVKOO    , isInputReal );
    solv.Solve( data.Cl2      , isInputReal );
    solv.Solve( data.ISOPND   , isInputReal );
    solv.Solve( data.HOBr     , isInputReal );
    solv.Solve( data.A3O2     , isInputReal );
    solv.Solve( data.PROPNN   , isInputReal );
    solv.Solve( data.GLYX     , isInputReal );
    solv.Solve( data.MAOPO2   , isInputReal );
    solv.Solve( data.CH4      , isInputReal );
    solv.Solve( data.GAOO     , isInputReal );
    solv.Solve( data.B3O2     , isInputReal );
    solv.Solve( data.ACET     , isInputReal );
    solv.Solve( data.MACRN    , isInputReal );
    solv.Solve( data.CH2OO    , isInputReal );
    solv.Solve( data.MGLYOO   , isInputReal );
    solv.Solve( data.VRO2     , isInputReal );
    solv.Solve( data.MGLOO    , isInputReal );
    solv.Solve( data.MACROO   , isInputReal );
    solv.Solve( data.PO2      , isInputReal );
    solv.Solve( data.CH3CHOO  , isInputReal );
    solv.Solve( data.MAN2     , isInputReal );
    solv.Solve( data.ISNOOA   , isInputReal );
    solv.Solve( data.H2O2     , isInputReal );
    solv.Solve( data.PRN1     , isInputReal );
    solv.Solve( data.ETO2     , isInputReal );
    solv.Solve( data.KO2      , isInputReal );
    solv.Solve( data.RCO3     , isInputReal );
    solv.Solve( data.HC5OO    , isInputReal );
    solv.Solve( data.GLYC     , isInputReal );
    solv.Solve( data.ClNO3    , isInputReal );
    solv.Solve( data.RIO2     , isInputReal );
    solv.Solve( data.R4N1     , isInputReal );
    solv.Solve( data.HOCl     , isInputReal );
    solv.Solve( data.ATO2     , isInputReal );
    solv.Solve( data.HNO3     , isInputReal );
    solv.Solve( data.ISN1     , isInputReal );
    solv.Solve( data.MAO3     , isInputReal );
    solv.Solve( data.MRO2     , isInputReal );
    solv.Solve( data.INO2     , isInputReal );
    solv.Solve( data.HAC      , isInputReal );
    solv.Solve( data.HC5      , isInputReal );
    solv.Solve( data.MGLY     , isInputReal );
    solv.Solve( data.ISOPNBO2 , isInputReal );
    solv.Solve( data.ISOPNDO2 , isInputReal );
    solv.Solve( data.R4O2     , isInputReal );
    solv.Solve( data.R4N2     , isInputReal );
    solv.Solve( data.BrO      , isInputReal );
    solv.Solve( data.RCHO     , isInputReal );
    solv.Solve( data.MEK      , isInputReal );
    solv.Solve( data.ClO      , isInputReal );
    solv.Solve( data.MACR     , isInputReal );
    solv.Solve( data.SO2      , isInputReal );
    solv.Solve( data.MVK      , isInputReal );
    solv.Solve( data.ALD2     , isInputReal );
    solv.Solve( data.MCO3     , isInputReal );
    solv.Solve( data.CH2O     , isInputReal );
    solv.Solve( data.H2O      , isInputReal );
    solv.Solve( data.Br       , isInputReal );
    solv.Solve( data.NO       , isInputReal );
    solv.Solve( data.NO3      , isInputReal );
    solv.Solve( data.Cl       , isInputReal );
    solv.Solve( data.O        , isInputReal );
    solv.Solve( data.O1D      , isInputReal );
    solv.Solve( data.O3       , isInputReal );
    solv.Solve( data.HO2      , isInputReal );
    solv.Solve( data.NO2      , isInputReal );
    solv.Solve( data.OH       , isInputReal );
    solv.Solve( data.HBr      , isInputReal );
    solv.Solve( data.HCl      , isInputReal );
    solv.Solve( data.CO       , isInputReal );
    solv.Solve( data.MO2      , isInputReal );

} /* End of CallSolver */

