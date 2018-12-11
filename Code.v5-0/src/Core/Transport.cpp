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

#include "Core/Structure.hpp"
#include "SANDS/Solver.hpp"

void Transport( Solution& data, SANDS::Solver& Solver )
{
        
    /* Gaseous species */
    Solver.Run( data.CO2      );
    Solver.Run( data.PPN      );
    Solver.Run( data.BrNO2    );
    Solver.Run( data.IEPOX    );
    Solver.Run( data.PMNN     );
    Solver.Run( data.N2O      );
    Solver.Run( data.N        );
    Solver.Run( data.PAN      );
    Solver.Run( data.ALK4     );
    Solver.Run( data.MAP      );
    Solver.Run( data.MPN      );
    Solver.Run( data.Cl2O2    );
    Solver.Run( data.ETP      );
    Solver.Run( data.HNO2     );
    Solver.Run( data.C3H8     );
    Solver.Run( data.RA3P     );
    Solver.Run( data.RB3P     );
    Solver.Run( data.OClO     );
    Solver.Run( data.ClNO2    );
    Solver.Run( data.ISOP     );
    Solver.Run( data.HNO4     );
    Solver.Run( data.MAOP     );
    Solver.Run( data.MP       );
    Solver.Run( data.ClOO     );
    Solver.Run( data.RP       );
    Solver.Run( data.BrCl     );
    Solver.Run( data.PP       );
    Solver.Run( data.PRPN     );
    Solver.Run( data.SO4      );
    Solver.Run( data.Br2      );
    Solver.Run( data.ETHLN    );
    Solver.Run( data.MVKN     );
    Solver.Run( data.R4P      );
    Solver.Run( data.C2H6     );
    Solver.Run( data.RIP      );
    Solver.Run( data.VRP      );
    Solver.Run( data.ATOOH    );
    Solver.Run( data.IAP      );
    Solver.Run( data.DHMOB    );
    Solver.Run( data.MOBA     );
    Solver.Run( data.MRP      );
    Solver.Run( data.N2O5     );
    Solver.Run( data.ISNOHOO  );
    Solver.Run( data.ISNP     );
    Solver.Run( data.ISOPNB   );
    Solver.Run( data.IEPOXOO  );
    Solver.Run( data.MACRNO2  );
    Solver.Run( data.ROH      );
    Solver.Run( data.MOBAOO   );
    Solver.Run( data.DIBOO    );
    Solver.Run( data.PMN      );
    Solver.Run( data.ISNOOB   );
    Solver.Run( data.INPN     );
    Solver.Run( data.H        );
    Solver.Run( data.BrNO3    );
    Solver.Run( data.PRPE     );
    Solver.Run( data.MVKOO    );
    Solver.Run( data.Cl2      );
    Solver.Run( data.ISOPND   );
    Solver.Run( data.HOBr     );
    Solver.Run( data.A3O2     );
    Solver.Run( data.PROPNN   );
    Solver.Run( data.GLYX     );
    Solver.Run( data.MAOPO2   );
    Solver.Run( data.CH4      );
    Solver.Run( data.GAOO     );
    Solver.Run( data.B3O2     );
    Solver.Run( data.ACET     );
    Solver.Run( data.MACRN    );
    Solver.Run( data.CH2OO    );
    Solver.Run( data.MGLYOO   );
    Solver.Run( data.VRO2     );
    Solver.Run( data.MGLOO    );
    Solver.Run( data.MACROO   );
    Solver.Run( data.PO2      );
    Solver.Run( data.CH3CHOO  );
    Solver.Run( data.MAN2     );
    Solver.Run( data.ISNOOA   );
    Solver.Run( data.H2O2     );
    Solver.Run( data.PRN1     );
    Solver.Run( data.ETO2     );
    Solver.Run( data.KO2      );
    Solver.Run( data.RCO3     );
    Solver.Run( data.HC5OO    );
    Solver.Run( data.GLYC     );
    Solver.Run( data.ClNO3    );
    Solver.Run( data.RIO2     );
    Solver.Run( data.R4N1     );
    Solver.Run( data.HOCl     );
    Solver.Run( data.ATO2     );
    Solver.Run( data.HNO3     );
    Solver.Run( data.ISN1     );
    Solver.Run( data.MAO3     );
    Solver.Run( data.MRO2     );
    Solver.Run( data.INO2     );
    Solver.Run( data.HAC      );
    Solver.Run( data.HC5      );
    Solver.Run( data.MGLY     );
    Solver.Run( data.ISOPNBO2 );
    Solver.Run( data.ISOPNDO2 );
    Solver.Run( data.R4O2     );
    Solver.Run( data.R4N2     );
    Solver.Run( data.BrO      );
    Solver.Run( data.RCHO     );
    Solver.Run( data.MEK      );
    Solver.Run( data.ClO      );
    Solver.Run( data.MACR     );
    Solver.Run( data.SO2      );
    Solver.Run( data.MVK      );
    Solver.Run( data.ALD2     );
    Solver.Run( data.MCO3     );
    Solver.Run( data.CH2O     );
    Solver.Run( data.H2O      );
    Solver.Run( data.Br       );
    Solver.Run( data.NO       );
    Solver.Run( data.NO3      );
    Solver.Run( data.Cl       );
    Solver.Run( data.O        );
    Solver.Run( data.O1D      );
    Solver.Run( data.O3       );
    Solver.Run( data.HO2      );
    Solver.Run( data.NO2      );
    Solver.Run( data.OH       );
    Solver.Run( data.HBr      );
    Solver.Run( data.HCl      );
    Solver.Run( data.CO       );
    Solver.Run( data.MO2      );

} /* End of Transport */

