/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Struct_SetToAmbient Program File                                 */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : Struct_SetToAmbient.cpp                   */
/* Working directory    : /home/fritzt/APCEMM-SourceCode            */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "Structure.h"
#include "Parameters.h"
#include "PhysConstant.h"

double pSat_H2Ol( double T );

void Struct_SetToAmbient( Solution &sol, unsigned int n_x, unsigned int n_y, \
                          double airDens, double temperature, double relHum )
{
    double amb_Value[N_SPC];
    char const *fileName("./src/Ambient.txt");
    std::ifstream file;

    file.open( fileName );

    if ( file.is_open() )
    {
        std::cout << "Reading ambient data from file: " << fileName << std::endl;
        std::string line;
        unsigned int i = 0;

        while ( std::getline( file, line ) ) {
            if ( line[0] != '#' ) {
                std::istringstream iss(line);
                iss >> amb_Value[i];
                i++;
            }
        }
        file.close();
    }
    else
    {
        std::string const currFile("Struct_SetToAmbient.cpp");
        std::cout << "ERROR: In " << currFile << ": Can't read (" << fileName << ")" << std::endl;
    }

    /* Gaseous species */
    sol.SetShape( sol.CO2  , n_x, n_y, amb_Value[  0] * airDens );
    sol.SetShape( sol.PPN  , n_x, n_y, amb_Value[  1] * airDens );
    sol.SetShape( sol.BrNO2, n_x, n_y, amb_Value[  2] * airDens );
    sol.SetShape( sol.IEPOX, n_x, n_y, amb_Value[  3] * airDens );
    sol.SetShape( sol.PMNN , n_x, n_y, amb_Value[  4] * airDens );
    sol.SetShape( sol.N2O  , n_x, n_y, amb_Value[  5] * airDens );
    sol.SetShape( sol.N    , n_x, n_y, amb_Value[  6] * airDens );
    sol.SetShape( sol.PAN  , n_x, n_y, amb_Value[  7] * airDens );
    sol.SetShape( sol.ALK4 , n_x, n_y, amb_Value[  8] * airDens );
    sol.SetShape( sol.MAP  , n_x, n_y, amb_Value[  9] * airDens );
    sol.SetShape( sol.MPN  , n_x, n_y, amb_Value[ 10] * airDens );
    sol.SetShape( sol.Cl2O2, n_x, n_y, amb_Value[ 11] * airDens );
    sol.SetShape( sol.ETP  , n_x, n_y, amb_Value[ 12] * airDens );
    sol.SetShape( sol.HNO2 , n_x, n_y, amb_Value[ 13] * airDens );
    sol.SetShape( sol.C3H8 , n_x, n_y, amb_Value[ 14] * airDens );
    sol.SetShape( sol.RA3P , n_x, n_y, amb_Value[ 15] * airDens );
    sol.SetShape( sol.RB3P , n_x, n_y, amb_Value[ 16] * airDens );
    sol.SetShape( sol.OClO , n_x, n_y, amb_Value[ 17] * airDens );
    sol.SetShape( sol.ClNO2, n_x, n_y, amb_Value[ 18] * airDens );
    sol.SetShape( sol.ISOP , n_x, n_y, amb_Value[ 19] * airDens );
    sol.SetShape( sol.HNO4 , n_x, n_y, amb_Value[ 20] * airDens );
    sol.SetShape( sol.MAOP , n_x, n_y, amb_Value[ 21] * airDens );
    sol.SetShape( sol.MP   , n_x, n_y, amb_Value[ 22] * airDens );
    sol.SetShape( sol.ClOO , n_x, n_y, amb_Value[ 23] * airDens );
    sol.SetShape( sol.RP   , n_x, n_y, amb_Value[ 24] * airDens );
    sol.SetShape( sol.BrCl , n_x, n_y, amb_Value[ 25] * airDens );
    sol.SetShape( sol.PP   , n_x, n_y, amb_Value[ 26] * airDens );
    sol.SetShape( sol.PRPN , n_x, n_y, amb_Value[ 27] * airDens );
    sol.SetShape( sol.SO4  , n_x, n_y, amb_Value[ 28] * airDens );
    sol.SetShape( sol.Br2  , n_x, n_y, amb_Value[ 29] * airDens );
    sol.SetShape( sol.ETHLN, n_x, n_y, amb_Value[ 30] * airDens );
    sol.SetShape( sol.MVKN , n_x, n_y, amb_Value[ 31] * airDens );
    sol.SetShape( sol.R4P  , n_x, n_y, amb_Value[ 32] * airDens );
    sol.SetShape( sol.C2H6 , n_x, n_y, amb_Value[ 33] * airDens );
    sol.SetShape( sol.RIP  , n_x, n_y, amb_Value[ 34] * airDens );
    sol.SetShape( sol.VRP  , n_x, n_y, amb_Value[ 35] * airDens );
    sol.SetShape( sol.ATOOH, n_x, n_y, amb_Value[ 36] * airDens );
    sol.SetShape( sol.IAP  , n_x, n_y, amb_Value[ 37] * airDens );
    sol.SetShape( sol.DHMOB, n_x, n_y, amb_Value[ 38] * airDens );
    sol.SetShape( sol.MOBA , n_x, n_y, amb_Value[ 39] * airDens );
    sol.SetShape( sol.MRP  , n_x, n_y, amb_Value[ 40] * airDens );
    sol.SetShape( sol.N2O5 , n_x, n_y, amb_Value[ 41] * airDens );
    sol.SetShape( sol.ISNOHOO, n_x, n_y, amb_Value[ 42] * airDens );
    sol.SetShape( sol.ISNP , n_x, n_y, amb_Value[ 43] * airDens );
    sol.SetShape( sol.ISOPNB, n_x, n_y, amb_Value[ 44] * airDens );
    sol.SetShape( sol.IEPOXOO, n_x, n_y, amb_Value[ 45] * airDens );
    sol.SetShape( sol.MACRNO2, n_x, n_y, amb_Value[ 46] * airDens );
    sol.SetShape( sol.ROH  , n_x, n_y, amb_Value[ 47] * airDens );
    sol.SetShape( sol.MOBAOO, n_x, n_y, amb_Value[ 48] * airDens );
    sol.SetShape( sol.DIBOO, n_x, n_y, amb_Value[ 49] * airDens );
    sol.SetShape( sol.PMN  , n_x, n_y, amb_Value[ 50] * airDens );
    sol.SetShape( sol.ISNOOB, n_x, n_y, amb_Value[ 51] * airDens );
    sol.SetShape( sol.INPN , n_x, n_y, amb_Value[ 52] * airDens );
    sol.SetShape( sol.H    , n_x, n_y, amb_Value[ 53] * airDens );
    sol.SetShape( sol.BrNO3, n_x, n_y, amb_Value[ 54] * airDens );
    sol.SetShape( sol.PRPE , n_x, n_y, amb_Value[ 55] * airDens );
    sol.SetShape( sol.MVKOO, n_x, n_y, amb_Value[ 56] * airDens );
    sol.SetShape( sol.Cl2  , n_x, n_y, amb_Value[ 57] * airDens );
    sol.SetShape( sol.ISOPND, n_x, n_y, amb_Value[ 58] * airDens );
    sol.SetShape( sol.HOBr , n_x, n_y, amb_Value[ 59] * airDens );
    sol.SetShape( sol.A3O2 , n_x, n_y, amb_Value[ 60] * airDens );
    sol.SetShape( sol.PROPNN, n_x, n_y, amb_Value[ 61] * airDens );
    sol.SetShape( sol.GLYX , n_x, n_y, amb_Value[ 62] * airDens );
    sol.SetShape( sol.MAOPO2, n_x, n_y, amb_Value[ 63] * airDens );
    sol.SetShape( sol.CH4  , n_x, n_y, amb_Value[ 64] * airDens );
    sol.SetShape( sol.GAOO , n_x, n_y, amb_Value[ 65] * airDens );
    sol.SetShape( sol.B3O2 , n_x, n_y, amb_Value[ 66] * airDens );
    sol.SetShape( sol.ACET , n_x, n_y, amb_Value[ 67] * airDens );
    sol.SetShape( sol.MACRN, n_x, n_y, amb_Value[ 68] * airDens );
    sol.SetShape( sol.CH2OO, n_x, n_y, amb_Value[ 69] * airDens );
    sol.SetShape( sol.MGLYOO, n_x, n_y, amb_Value[ 70] * airDens );
    sol.SetShape( sol.VRO2 , n_x, n_y, amb_Value[ 71] * airDens );
    sol.SetShape( sol.MGLOO, n_x, n_y, amb_Value[ 72] * airDens );
    sol.SetShape( sol.MACROO, n_x, n_y, amb_Value[ 73] * airDens );
    sol.SetShape( sol.PO2  , n_x, n_y, amb_Value[ 74] * airDens );
    sol.SetShape( sol.CH3CHOO, n_x, n_y, amb_Value[ 75] * airDens );
    sol.SetShape( sol.MAN2 , n_x, n_y, amb_Value[ 76] * airDens );
    sol.SetShape( sol.ISNOOA, n_x, n_y, amb_Value[ 77] * airDens );
    sol.SetShape( sol.H2O2 , n_x, n_y, amb_Value[ 78] * airDens );
    sol.SetShape( sol.PRN1 , n_x, n_y, amb_Value[ 79] * airDens );
    sol.SetShape( sol.ETO2 , n_x, n_y, amb_Value[ 80] * airDens );
    sol.SetShape( sol.KO2  , n_x, n_y, amb_Value[ 81] * airDens );
    sol.SetShape( sol.RCO3 , n_x, n_y, amb_Value[ 82] * airDens );
    sol.SetShape( sol.HC5OO, n_x, n_y, amb_Value[ 83] * airDens );
    sol.SetShape( sol.GLYC , n_x, n_y, amb_Value[ 84] * airDens );
    sol.SetShape( sol.ClNO3, n_x, n_y, amb_Value[ 85] * airDens );
    sol.SetShape( sol.RIO2 , n_x, n_y, amb_Value[ 86] * airDens );
    sol.SetShape( sol.R4N1 , n_x, n_y, amb_Value[ 87] * airDens );
    sol.SetShape( sol.HOCl , n_x, n_y, amb_Value[ 88] * airDens );
    sol.SetShape( sol.ATO2 , n_x, n_y, amb_Value[ 89] * airDens );
    sol.SetShape( sol.HNO3 , n_x, n_y, amb_Value[ 90] * airDens );
    sol.SetShape( sol.ISN1 , n_x, n_y, amb_Value[ 91] * airDens );
    sol.SetShape( sol.MAO3 , n_x, n_y, amb_Value[ 92] * airDens );
    sol.SetShape( sol.MRO2 , n_x, n_y, amb_Value[ 93] * airDens );
    sol.SetShape( sol.INO2 , n_x, n_y, amb_Value[ 94] * airDens );
    sol.SetShape( sol.HAC  , n_x, n_y, amb_Value[ 95] * airDens );
    sol.SetShape( sol.HC5  , n_x, n_y, amb_Value[ 96] * airDens );
    sol.SetShape( sol.MGLY , n_x, n_y, amb_Value[ 97] * airDens );
    sol.SetShape( sol.ISOPNBO2, n_x, n_y, amb_Value[ 98] * airDens );
    sol.SetShape( sol.ISOPNDO2, n_x, n_y, amb_Value[ 99] * airDens );
    sol.SetShape( sol.R4O2 , n_x, n_y, amb_Value[100] * airDens );
    sol.SetShape( sol.R4N2 , n_x, n_y, amb_Value[101] * airDens );
    sol.SetShape( sol.BrO  , n_x, n_y, amb_Value[102] * airDens );
    sol.SetShape( sol.RCHO , n_x, n_y, amb_Value[103] * airDens );
    sol.SetShape( sol.MEK  , n_x, n_y, amb_Value[104] * airDens );
    sol.SetShape( sol.ClO  , n_x, n_y, amb_Value[105] * airDens );
    sol.SetShape( sol.MACR , n_x, n_y, amb_Value[106] * airDens );
    sol.SetShape( sol.SO2  , n_x, n_y, amb_Value[107] * airDens );
    sol.SetShape( sol.MVK  , n_x, n_y, amb_Value[108] * airDens );
    sol.SetShape( sol.ALD2 , n_x, n_y, amb_Value[109] * airDens );
    sol.SetShape( sol.MCO3 , n_x, n_y, amb_Value[110] * airDens );
    sol.SetShape( sol.CH2O , n_x, n_y, amb_Value[111] * airDens );
    sol.SetShape( sol.H2O  , n_x, n_y, amb_Value[112] * airDens );
    sol.SetShape( sol.Br   , n_x, n_y, amb_Value[113] * airDens );
    sol.SetShape( sol.NO   , n_x, n_y, amb_Value[114] * airDens );
    sol.SetShape( sol.NO3  , n_x, n_y, amb_Value[115] * airDens );
    sol.SetShape( sol.Cl   , n_x, n_y, amb_Value[116] * airDens );
    sol.SetShape( sol.O    , n_x, n_y, amb_Value[117] * airDens );
    sol.SetShape( sol.O1D  , n_x, n_y, amb_Value[118] * airDens );
    sol.SetShape( sol.O3   , n_x, n_y, amb_Value[119] * airDens );
    sol.SetShape( sol.HO2  , n_x, n_y, amb_Value[120] * airDens );
    sol.SetShape( sol.NO2  , n_x, n_y, amb_Value[121] * airDens );
    sol.SetShape( sol.OH   , n_x, n_y, amb_Value[122] * airDens );
    sol.SetShape( sol.HBr  , n_x, n_y, amb_Value[123] * airDens );
    sol.SetShape( sol.HCl  , n_x, n_y, amb_Value[124] * airDens );
    sol.SetShape( sol.CO   , n_x, n_y, amb_Value[125] * airDens );
    sol.SetShape( sol.MO2  , n_x, n_y, amb_Value[126] * airDens );
    sol.SetShape( sol.ACTA , n_x, n_y, amb_Value[127] * airDens );
    sol.SetShape( sol.EOH  , n_x, n_y, amb_Value[128] * airDens );
    sol.SetShape( sol.H2   , n_x, n_y, amb_Value[129] * airDens );
    sol.SetShape( sol.HCOOH, n_x, n_y, amb_Value[130] * airDens );
    sol.SetShape( sol.MOH  , n_x, n_y, amb_Value[131] * airDens );
    sol.SetShape( sol.N2   , n_x, n_y, amb_Value[132] * airDens );
    sol.SetShape( sol.O2   , n_x, n_y, amb_Value[133] * airDens );
    sol.SetShape( sol.RCOOH, n_x, n_y, amb_Value[134] * airDens );

    /* Setting up ambient water vapor */
    for ( unsigned int i = 0; i < sol.H2O.size(); i++ ) {
        for ( unsigned int j = 0; j < sol.H2O[0].size(); j++ ) {
            sol.H2O[i][j] = (relHum/((double) 100.0) * \
                             pSat_H2Ol( temperature ) / ( kB * temperature )) / 1.00E+06;
            /* RH_w = x_H2O*P/Psat_H2Ol = [H2O]*1E6 * kB*T/Psat_H2Ol */

        }
    }

    /* Liquid/solid species */
    sol.SetShape( sol.SO4L , n_x, n_y, (double) 0.0 );
    sol.SetShape( sol.H2OL , n_x, n_y, (double) 0.0 );
    sol.SetShape( sol.H2OS , n_x, n_y, (double) 0.0 );
    sol.SetShape( sol.HNO3L, n_x, n_y, (double) 0.0 );
    sol.SetShape( sol.HNO3S, n_x, n_y, (double) 0.0 );
    sol.SetShape( sol.HClL , n_x, n_y, (double) 0.0 );
    sol.SetShape( sol.HOClL, n_x, n_y, (double) 0.0 );
    sol.SetShape( sol.HBrL , n_x, n_y, (double) 0.0 );
    sol.SetShape( sol.HOBrL, n_x, n_y, (double) 0.0 );

    /* Aerosols */
    sol.SetShape( sol.Soot, n_x, n_y, (double) 0.0 );
    sol.SetShape( sol.SLA , n_x, n_y, (double) 0.0 );
    sol.SetShape( sol.SPA , n_x, n_y, (double) 0.0 );

}


