#include <iostream>
#include <vector>
#include <complex>

#include "Structure.h"

void DiffusionSolver( std::vector<std::vector<double> >& vect, \
                      std::vector<std::vector<double> >& diffFactor, \
                      std::vector<std::vector<std::complex<double> > >& advFactor, \
                      const char* fileName_FFTW, \
                      const bool realInput );

void CallSpeciesDiffusion( Solution& Data, \
                           std::vector<std::vector<double> >& diffFactor, \
                           std::vector<std::vector<std::complex<double> > >& advFactor, \
                           const char* fileName_FFTW)
{
        
    const bool isInputReal = 1;


    DiffusionSolver( Data.O3, diffFactor, advFactor, fileName_FFTW, isInputReal );
    /* Gaseous species */
    DiffusionSolver( Data.CO2  , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.PPN  , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.BrNO2, diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.IEPOX, diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.PMNN , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.N2O  , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.N    , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.PAN  , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.ALK4 , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.MAP  , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.MPN  , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.Cl2O2, diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.ETP  , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.HNO2 , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.C3H8 , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.RA3P , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.RB3P , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.OClO , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.ClNO2, diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.ISOP , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.HNO4 , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.MAOP , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.MP   , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.ClOO , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.RP   , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.BrCl , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.PP   , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.PRPN , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.SO4  , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.Br2  , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.ETHLN, diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.MVKN , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.R4P  , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.C2H6 , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.RIP  , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.VRP  , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.ATOOH, diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.IAP  , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.DHMOB, diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.MOBA , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.MRP  , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.N2O5 , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.ISNOHOO, diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.ISNP , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.ISOPNB, diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.IEPOXOO, diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.MACRNO2, diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.ROH  , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.MOBAOO, diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.DIBOO, diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.PMN  , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.ISNOOB, diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.INPN , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.H    , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.BrNO3, diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.PRPE , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.MVKOO, diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.Cl2  , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.ISOPND, diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.HOBr , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.A3O2 , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.PROPNN, diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.GLYX , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.MAOPO2, diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.CH4  , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.GAOO , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.B3O2 , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.ACET , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.MACRN, diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.CH2OO, diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.MGLYOO, diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.VRO2 , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.MGLOO, diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.MACROO, diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.PO2  , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.CH3CHOO, diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.MAN2 , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.ISNOOA, diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.H2O2 , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.PRN1 , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.ETO2 , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.KO2  , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.RCO3 , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.HC5OO, diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.GLYC , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.ClNO3, diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.RIO2 , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.R4N1 , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.HOCl , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.ATO2 , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.HNO3 , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.ISN1 , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.MAO3 , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.MRO2 , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.INO2 , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.HAC  , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.HC5  , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.MGLY , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.ISOPNBO2, diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.ISOPNDO2, diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.R4O2 , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.R4N2 , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.BrO  , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.RCHO , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.MEK  , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.ClO  , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.MACR , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.SO2  , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.MVK  , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.ALD2 , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.MCO3 , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.CH2O , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.H2O  , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.Br   , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.NO   , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.NO3  , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.Cl   , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.O    , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.O1D  , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.O3   , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.HO2  , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.NO2  , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.OH   , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.HBr  , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.HCl  , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.CO   , diffFactor, advFactor, fileName_FFTW, isInputReal );
    DiffusionSolver( Data.MO2  , diffFactor, advFactor, fileName_FFTW, isInputReal );

} /* End of CallSpeciesDiffusion */

