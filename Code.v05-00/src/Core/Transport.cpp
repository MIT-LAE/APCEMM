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
#include <fftw3.h>

#include "Util/ForwardDecl.hpp"
#include "Core/Structure.hpp"
#include "SANDS/Solver.hpp"

void Transport( Solution& data, SANDS::Solver& Solver, \
                const Vector_2D &cellAreas )
{
        
    /* Gaseous species */
    Solver.Run( data.CO2     , cellAreas );
    Solver.Run( data.PPN     , cellAreas );
    Solver.Run( data.BrNO2   , cellAreas );
    Solver.Run( data.IEPOX   , cellAreas );
    Solver.Run( data.PMNN    , cellAreas );
    Solver.Run( data.N2O     , cellAreas );
    Solver.Run( data.N       , cellAreas );
    Solver.Run( data.PAN     , cellAreas );
    Solver.Run( data.ALK4    , cellAreas );
    Solver.Run( data.MAP     , cellAreas );
    Solver.Run( data.MPN     , cellAreas );
    Solver.Run( data.Cl2O2   , cellAreas );
    Solver.Run( data.ETP     , cellAreas );
    Solver.Run( data.HNO2    , cellAreas );
    Solver.Run( data.C3H8    , cellAreas );
    Solver.Run( data.RA3P    , cellAreas );
    Solver.Run( data.RB3P    , cellAreas );
    Solver.Run( data.OClO    , cellAreas );
    Solver.Run( data.ClNO2   , cellAreas );
    Solver.Run( data.ISOP    , cellAreas );
    Solver.Run( data.HNO4    , cellAreas );
    Solver.Run( data.MAOP    , cellAreas );
    Solver.Run( data.MP      , cellAreas );
    Solver.Run( data.ClOO    , cellAreas );
    Solver.Run( data.RP      , cellAreas );
    Solver.Run( data.BrCl    , cellAreas );
    Solver.Run( data.PP      , cellAreas );
    Solver.Run( data.PRPN    , cellAreas );
    Solver.Run( data.SO4     , cellAreas );
    Solver.Run( data.Br2     , cellAreas );
    Solver.Run( data.ETHLN   , cellAreas );
    Solver.Run( data.MVKN    , cellAreas );
    Solver.Run( data.R4P     , cellAreas );
    Solver.Run( data.C2H6    , cellAreas );
    Solver.Run( data.RIP     , cellAreas );
    Solver.Run( data.VRP     , cellAreas );
    Solver.Run( data.ATOOH   , cellAreas );
    Solver.Run( data.IAP     , cellAreas );
    Solver.Run( data.DHMOB   , cellAreas );
    Solver.Run( data.MOBA    , cellAreas );
    Solver.Run( data.MRP     , cellAreas );
    Solver.Run( data.N2O5    , cellAreas );
    Solver.Run( data.ISNOHOO , cellAreas );
    Solver.Run( data.ISNP    , cellAreas );
    Solver.Run( data.ISOPNB  , cellAreas );
    Solver.Run( data.IEPOXOO , cellAreas );
    Solver.Run( data.MACRNO2 , cellAreas );
    Solver.Run( data.ROH     , cellAreas );
    Solver.Run( data.MOBAOO  , cellAreas );
    Solver.Run( data.DIBOO   , cellAreas );
    Solver.Run( data.PMN     , cellAreas );
    Solver.Run( data.ISNOOB  , cellAreas );
    Solver.Run( data.INPN    , cellAreas );
    Solver.Run( data.H       , cellAreas );
    Solver.Run( data.BrNO3   , cellAreas );
    Solver.Run( data.PRPE    , cellAreas );
    Solver.Run( data.MVKOO   , cellAreas );
    Solver.Run( data.Cl2     , cellAreas );
    Solver.Run( data.ISOPND  , cellAreas );
    Solver.Run( data.HOBr    , cellAreas );
    Solver.Run( data.A3O2    , cellAreas );
    Solver.Run( data.PROPNN  , cellAreas );
    Solver.Run( data.GLYX    , cellAreas );
    Solver.Run( data.MAOPO2  , cellAreas );
    Solver.Run( data.CH4     , cellAreas );
    Solver.Run( data.GAOO    , cellAreas );
    Solver.Run( data.B3O2    , cellAreas );
    Solver.Run( data.ACET    , cellAreas );
    Solver.Run( data.MACRN   , cellAreas );
    Solver.Run( data.CH2OO   , cellAreas );
    Solver.Run( data.MGLYOO  , cellAreas );
    Solver.Run( data.VRO2    , cellAreas );
    Solver.Run( data.MGLOO   , cellAreas );
    Solver.Run( data.MACROO  , cellAreas );
    Solver.Run( data.PO2     , cellAreas );
    Solver.Run( data.CH3CHOO , cellAreas );
    Solver.Run( data.MAN2    , cellAreas );
    Solver.Run( data.ISNOOA  , cellAreas );
    Solver.Run( data.H2O2    , cellAreas );
    Solver.Run( data.PRN1    , cellAreas );
    Solver.Run( data.ETO2    , cellAreas );
    Solver.Run( data.KO2     , cellAreas );
    Solver.Run( data.RCO3    , cellAreas );
    Solver.Run( data.HC5OO   , cellAreas );
    Solver.Run( data.GLYC    , cellAreas );
    Solver.Run( data.ClNO3   , cellAreas );
    Solver.Run( data.RIO2    , cellAreas );
    Solver.Run( data.R4N1    , cellAreas );
    Solver.Run( data.HOCl    , cellAreas );
    Solver.Run( data.ATO2    , cellAreas );
    Solver.Run( data.HNO3    , cellAreas );
    Solver.Run( data.ISN1    , cellAreas );
    Solver.Run( data.MAO3    , cellAreas );
    Solver.Run( data.MRO2    , cellAreas );
    Solver.Run( data.INO2    , cellAreas );
    Solver.Run( data.HAC     , cellAreas );
    Solver.Run( data.HC5     , cellAreas );
    Solver.Run( data.MGLY    , cellAreas );
    Solver.Run( data.ISOPNBO2, cellAreas );
    Solver.Run( data.ISOPNDO2, cellAreas );
    Solver.Run( data.R4O2    , cellAreas );
    Solver.Run( data.R4N2    , cellAreas );
    Solver.Run( data.BrO     , cellAreas );
    Solver.Run( data.RCHO    , cellAreas );
    Solver.Run( data.MEK     , cellAreas );
    Solver.Run( data.ClO     , cellAreas );
    Solver.Run( data.MACR    , cellAreas );
    Solver.Run( data.SO2     , cellAreas );
    Solver.Run( data.MVK     , cellAreas );
    Solver.Run( data.ALD2    , cellAreas );
    Solver.Run( data.MCO3    , cellAreas );
    Solver.Run( data.CH2O    , cellAreas );
    Solver.Run( data.H2O     , cellAreas, 1 );
    Solver.Run( data.Br      , cellAreas );
    Solver.Run( data.NO      , cellAreas );
    Solver.Run( data.NO3     , cellAreas );
    Solver.Run( data.Cl      , cellAreas );
    Solver.Run( data.O       , cellAreas );
    Solver.Run( data.O1D     , cellAreas );
    Solver.Run( data.O3      , cellAreas );
    Solver.Run( data.HO2     , cellAreas );
    Solver.Run( data.NO2     , cellAreas );
    Solver.Run( data.OH      , cellAreas );
    Solver.Run( data.HBr     , cellAreas );
    Solver.Run( data.HCl     , cellAreas );
    Solver.Run( data.CO      , cellAreas );
    Solver.Run( data.MO2     , cellAreas );

} /* End of Transport */

