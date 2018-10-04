/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Structure Header File                                            */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : Structure.hpp                             */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef STRUCTURE_H_INCLUDED
#define STRUCTURE_H_INCLUDED

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>

#include "Core/Parameters.hpp"
#include "Util/PhysConstant.hpp"
#include "Util/PhysFunction.hpp"
#include "Core/Emission.hpp"
#include "Core/Aircraft.hpp"
#include "Core/Engine.hpp"

class Solution
{
    public:

        Solution( const unsigned int nVar, const unsigned int n_x, const unsigned int n_y );
        ~Solution();
        void Clear( std::vector<std::vector<double> >& vector_2D );
        void SetShape( std::vector<std::vector<double> >& vector_2D, unsigned int n_x, unsigned int n_y, double value = 0.0 );
        void SetToValue( std::vector<std::vector<double> >& vector_2D, double value = 0.0 );
        void Print( std::vector<std::vector<double> >& vector_2D, unsigned int i_max = 1, unsigned int j_max = 1 );
        void Initialize( char const *fileName, double temperature, double airDens, double relHum );
        void getData( double varArray[], double fixArray[], unsigned int i = 0, unsigned int j = 0 );
        void applyData( double varArray[], unsigned int i = 0, unsigned int j = 0 );
        void applyRing( double varArray[], double tempArray[], std::vector<std::vector<std::pair<unsigned int, unsigned int>>> mapRing2Mesh, unsigned int iRing );
        void applyAmbient( double varArray[], std::vector<std::vector<std::pair<unsigned int, unsigned int>>> mapRing2Mesh, unsigned int ambIndex );
        void addEmission( const Emission &EI, const Aircraft &ac, std::vector<std::vector<std::pair<unsigned int, unsigned int>>> &map, std::vector<std::vector<double>> cellAreas, bool halfRing );
        std::vector<double> getAmbient( ) const;
        unsigned int getNx() const;
        unsigned int getNy() const;
        void Debug( double airDens );

        /* Gaseous species */
        std::vector<std::vector<double> > CO2, PPN, BrNO2, IEPOX, PMNN, N2O, N, PAN,\
                                          ALK4, MAP, MPN, Cl2O2, ETP, HNO2, C3H8, RA3P,\
                                          RB3P, OClO, ClNO2, ISOP, HNO4, MAOP, MP, ClOO,\
                                          RP, BrCl, PP, PRPN, SO4, Br2, ETHLN, MVKN,\
                                          R4P, C2H6, RIP, VRP, ATOOH, IAP, DHMOB, MOBA,\
                                          MRP, N2O5, ISNOHOO, ISNP, ISOPNB, IEPOXOO, MACRNO2, ROH,\
                                          MOBAOO, DIBOO, PMN, ISNOOB, INPN, H, BrNO3, PRPE,\
                                          MVKOO, Cl2, ISOPND, HOBr, A3O2, PROPNN, GLYX, MAOPO2,\
                                          CH4, GAOO, B3O2, ACET, MACRN, CH2OO, MGLYOO, VRO2,\
                                          MGLOO, MACROO, PO2, CH3CHOO, MAN2, ISNOOA, H2O2, PRN1,\
                                          ETO2, KO2, RCO3, HC5OO, GLYC, ClNO3, RIO2, R4N1,\
                                          HOCl, ATO2, HNO3, ISN1, MAO3, MRO2, INO2, HAC,\
                                          HC5, MGLY, ISOPNBO2, ISOPNDO2, R4O2, R4N2, BrO, RCHO,\
                                          MEK, ClO, MACR, SO2, MVK, ALD2, MCO3, CH2O,\
                                          H2O, Br, NO, NO3, Cl, O, O1D, O3,\
                                          HO2, NO2, OH, HBr, HCl, CO, MO2, ACTA,\
                                          EOH, H2, HCOOH, MOH, N2, O2, RCOOH;

        /* Liquid species */
        std::vector<std::vector<double> > SO4L, H2OL, HNO3L, HClL, HOClL, HBrL, HOBrL;

        /* Solid species */
        std::vector<std::vector<double> > H2OS, HNO3S;

        /* Aerosols */
        std::vector<std::vector<double> > Soot, SLA, SPA;

    private:

        const unsigned int nVariables;
        const unsigned int size_x;
        const unsigned int size_y;

};

#endif /* STRUCTURE_H_INCLUDED */
