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
#include <cmath> 

#include "Core/Interface.hpp"
#include "Core/Parameters.hpp"
#include "Util/ForwardDecl.hpp"
#include "Core/Input.hpp"
#include "KPP/KPP_Parameters.h"
#include "KPP/KPP_Global.h"
#include "KPP/KPP.hpp"
#include "Util/PhysConstant.hpp"
#include "Util/PhysFunction.hpp"
#include "Core/Emission.hpp"
#include "Core/Aircraft.hpp"
#include "Core/Engine.hpp"
#include "Core/LiquidAer.hpp"
#include "Core/SZA.hpp"
#include "AIM/Aerosol.hpp"
#include "Core/Meteorology.hpp"

class Solution
{
    public:

        Solution();

        ~Solution();

        void Clear( Vector_2D& vector_2D );

        void SetShape( Vector_2D& vector_2D, \
                       const UInt n_x,       \
                       const UInt n_y,       \
                       const RealDouble value = 0.0 );

        void SetToValue( Vector_2D& vector_2D, \
                         const RealDouble value = 0.0 );

        void Print( const Vector_2D& vector_2D, \
                    const UInt i_max = 1,       \
                    const UInt j_max = 1 ) const;

        void Initialize( char const *fileName,     \
                         const Input &input,       \
                         const RealDouble airDens, \
                         const Meteorology &met,   \
                         const bool DBG );

        void getData( const UInt i = 0, \
                      const UInt j = 0 );

        void applyData( const UInt i = 0, \
                        const UInt j = 0 );

        void applyRing( RealDouble tempArray[],        \
                        const Vector_2Dui &mapIndices, \
                        const UInt iRing );

        void applyAmbient( const Vector_2Dui &mapIndices, \
                           const UInt iRing );

        void addEmission( const Emission &EI, const Aircraft &AC,            \
                          const Mesh &m,                                     \
                          bool halfRing,                                     \
                          const RealDouble temperature, bool set2Saturation, \
                          AIM::Aerosol liqAer, AIM::Aerosol iceAer,          \
                          const RealDouble Soot_Den,                         \
                          const Meteorology &met, const RealDouble areaPlume );

        Vector_1D getAmbient( ) const;
        Vector_1D getLiqSpecies() const;
        Vector_2D getAerosol( ) const;
        Vector_1D getAerosolDens( ) const;
        Vector_1D getAerosolRadi( ) const;
        Vector_1D getAerosolArea( ) const;

        void getAerosolProp( RealDouble ( &radi )[4], \
                             RealDouble ( &area )[4], \
                             RealDouble &IWC,         \
                             const Vector_2D &weights ) const;

        int SpinUp( Vector_1D &amb_Value,       \
                    const Input &input,         \
                    const RealDouble airDens,   \
                    const RealDouble startTime, \
                    const bool DGB = 0 );

        UInt Nx() const { return size_x; };
        UInt Ny() const { return size_y; };
        void Debug( const RealDouble airDens );

        /* Gaseous species */
        Vector_2D CO2, PPN, BrNO2, IEPOX, PMNN, N2O, N, PAN,               \
                  ALK4, MAP, MPN, Cl2O2, ETP, HNO2, C3H8, RA3P,            \
                  RB3P, OClO, ClNO2, ISOP, HNO4, MAOP, MP, ClOO,           \
                  RP, BrCl, PP, PRPN, SO4, Br2, ETHLN, MVKN,               \
                  R4P, C2H6, RIP, VRP, ATOOH, IAP, DHMOB, MOBA,            \
                  MRP, N2O5, ISNOHOO, ISNP, ISOPNB, IEPOXOO, MACRNO2, ROH, \
                  MOBAOO, DIBOO, PMN, ISNOOB, INPN, H, BrNO3, PRPE,        \
                  MVKOO, Cl2, ISOPND, HOBr, A3O2, PROPNN, GLYX, MAOPO2,    \
                  CH4, GAOO, B3O2, ACET, MACRN, CH2OO, MGLYOO, VRO2,       \
                  MGLOO, MACROO, PO2, CH3CHOO, MAN2, ISNOOA, H2O2, PRN1,   \
                  ETO2, KO2, RCO3, HC5OO, GLYC, ClNO3, RIO2, R4N1,         \
                  HOCl, ATO2, HNO3, ISN1, MAO3, MRO2, INO2, HAC,           \
                  HC5, MGLY, ISOPNBO2, ISOPNDO2, R4O2, R4N2, BrO, RCHO,    \
                  MEK, ClO, MACR, SO2, MVK, ALD2, MCO3, CH2O,              \
                  H2O, Br, NO, NO3, Cl, O, O1D, O3,                        \
                  HO2, NO2, OH, HBr, HCl, CO, MO2, ACTA,                   \
                  EOH, H2, HCOOH, MOH, N2, O2, RCOOH;

        /* Liquid species */
        Vector_2D SO4L, H2OL, HNO3L, HClL, HOClL, HBrL, HOBrL;

        /* Solid species */
        Vector_2D H2OS, HNO3S;

        Vector_2D NIT, NAT;

        /* Tracers */
        Vector_2D SO4T;

        /* Aerosols */
        Vector_2D sootDens, sootRadi, sootArea;

        AIM::Grid_Aerosol liquidAerosol, solidAerosol;

        UInt nBin_PA;
        UInt nBin_LA;

//        Vector_1D LA_rE;
//        Vector_1D LA_rJ;
//        Vector_1D LA_vJ;
//
//        Vector_1D PA_rE;
//        Vector_1D PA_rJ;
//        Vector_1D PA_vJ;

        RealDouble LA_nDens, LA_rEff, LA_SAD;
        RealDouble PA_nDens, PA_rEff, PA_SAD;

        AIM::Coagulation LA_Kernel, PA_Kernel;

        Vector_1D KHETI_SLA;
        Vector_1D AERFRAC, SOLIDFRAC;
        UInt STATE_PSC; 

    private:

        const UInt nVariables;
        const UInt nAer;
        const UInt size_x;
        const UInt size_y;

};

#endif /* STRUCTURE_H_INCLUDED */
