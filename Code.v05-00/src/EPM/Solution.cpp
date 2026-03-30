/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Structure Program File                                           */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : Structure.cpp                             */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <fstream>

#include "KPP/KPP.hpp"
#include "KPP/KPP_Parameters.h"
#include "Core/LiquidAer.hpp"
#include "Core/Parameters.hpp"
#include "Core/SZA.hpp"
#include "EPM/Solution.hpp"
#include "Util/PhysConstant.hpp"

using physConst::PI, physConst::kB;

namespace EPM {

Solution::Solution(const OptInput& optInput) :
        liquidAerosol(),
        solidAerosol(),
        nVariables(NSPEC),
        nAer(N_AER),
        size_x(optInput.ADV_GRID_NX),
        size_y(optInput.ADV_GRID_NY),
        reducedSize(0)
{
    VectorUtils::set_shape(sootDens, size_x, size_y);
    VectorUtils::set_shape(sootRadi, size_x, size_y);
    VectorUtils::set_shape(sootArea, size_x, size_y);
}

void Solution::Initialize(std::string fileName,
                          const Input &input,
                          const double airDens,
                          const Meteorology &met,
                          const OptInput &Input_Opt,
                          Vector_1D &varSpeciesArray,
                          const bool DBG) {
    Vector_1D amb_Value(NSPECALL, 0.0);
    Vector_2D aer_Value(nAer, Vector_1D(2, 0.0));

    /* Read input background conditions */
    readInputBackgroundConditions(input, amb_Value, aer_Value, fileName);

    const double AMBIENT_VALID_TIME = 8.0; //hours
    SpinUp(amb_Value, input, airDens, AMBIENT_VALID_TIME, varSpeciesArray);

    /* Enforce pre-defined values? *
     * Read input defined values for background concentrations */

    /* Split NOx as NO and NO2 using the current NO/NO2 ratio */

    /* Inputs are in ppb */

    Vector_1D amb_Value_old = amb_Value;

    setAmbientConcentrations(input, amb_Value);

      /* Initialize and allocate space for species */
    initializeSpeciesH2O(input, Input_Opt, amb_Value, airDens, met);

    Vector_1D stratData{ Species[ind_SO4][0][0],   Species[ind_HNO3][0][0],  \
                         Species[ind_HCl][0][0],   Species[ind_HOCl][0][0],  \
                         Species[ind_HBr][0][0],   Species[ind_HOBr][0][0],  \
                         Species[ind_H2O][0][0],   Species[ind_ClNO3][0][0], \
                         Species[ind_BrNO3][0][0], Species[ind_NIT][0][0],   \
                         Species[ind_NAT][0][0] };

    KHETI_SLA.assign( 11, 0.0 );
    AERFRAC.assign( 7, 0.0 );
    SOLIDFRAC.assign( 7, 0.0 );

    Vector_1D RAD, RHO, KG, NDENS, SAD;
    RAD.assign( 2, 0.0 );
    RHO.assign( 2, 0.0 );
    KG.assign( 2, 0.0 );
    NDENS.assign( 2, 0.0 );
    SAD.assign( 2, 0.0 );

    const double boxArea = (Input_Opt.ADV_GRID_XLIM_LEFT + Input_Opt.ADV_GRID_XLIM_RIGHT) * (Input_Opt.ADV_GRID_YLIM_UP + Input_Opt.ADV_GRID_YLIM_DOWN);
    STATE_PSC = STRAT_AER( input.temperature_K(), input.pressure_Pa(), airDens,  \
                           input.latitude_deg(), stratData,                      \
                           boxArea, KHETI_SLA, SOLIDFRAC,                        \
                           AERFRAC, RAD, RHO, KG, NDENS, SAD, Input_Opt.ADV_TROPOPAUSE_PRESSURE, DBG );

    setSpeciesValues(AERFRAC, SOLIDFRAC, stratData);

    /* Aerosols */
    /* Assume that soot particles are monodisperse */
    VectorUtils::set_value(sootDens, aer_Value[0][0]);
    VectorUtils::set_value(sootRadi, aer_Value[0][1]);
    VectorUtils::set_value(sootArea, 4.0 / double(3.0) * PI * aer_Value[0][0] * aer_Value[0][1] * aer_Value[0][1] * aer_Value[0][1]);
    
    double la_r_hig_low = LA_R_HIG/LA_R_LOW;
    nBin_LA = std::floor(1 + log(la_r_hig_low * la_r_hig_low * la_r_hig_low) / log(LA_VRAT));

    Vector_1D LA_rE( nBin_LA + 1, 0.0 ); /* Bin edges in m */
    Vector_1D LA_rJ( nBin_LA    , 0.0 ); /* Bin center radius in m */
    Vector_1D LA_vJ( nBin_LA    , 0.0 ); /* Bin volume centers in m^3 */

    const double LA_RRAT = cbrt( LA_VRAT );
    LA_rE[0] = LA_R_LOW;
    for ( UInt iBin_LA = 1; iBin_LA < nBin_LA + 1; iBin_LA++ )
        LA_rE[iBin_LA] = LA_rE[iBin_LA-1] * LA_RRAT; /* [m] */

    for ( UInt iBin_LA = 0; iBin_LA < nBin_LA; iBin_LA++ ) {
        LA_rJ[iBin_LA] = 0.5 * ( LA_rE[iBin_LA] + LA_rE[iBin_LA+1] );                       /* [m] */
        LA_vJ[iBin_LA] = 4.0 / double(3.0) * PI * \
                         ( LA_rE[iBin_LA] * LA_rE[iBin_LA] * LA_rE[iBin_LA] \
                         + LA_rE[iBin_LA+1] * LA_rE[iBin_LA+1] * LA_rE[iBin_LA+1] ) * 0.5;  /* [m^3] */
    }

    LA_nDens = NDENS[1] * 1.00E-06; /* [#/cm^3]      */
    LA_rEff  = RAD[1]   * 1.00E+09; /* [nm]          */
    LA_SAD   = SAD[1]   * 1.00E+06; /* [\mum^2/cm^3] */

    if ( LA_nDens > 0.0E+00 && LA_rEff > 0 && LA_SAD > 0) {
        /* For a lognormal distribution:
         * r_eff = r_m * exp( 5/2 * ln(S)^2 )
         * A     = 4\pi N0 r_m^2 * exp ( 2 * ln(S)^2 )
         * A/r_eff^2 = 4\pi N0 * exp( - 3 * ln(S)^2 )
         *
         * ln(S) = sqrt(-1/3*ln(A/(4\pi r_eff^2 * N0)));
         * r_m = r_eff * exp( -5/2 * ln(S)^2 ); */

        const double sLA = sqrt( - 1.0 / (3.0) * log(SAD[1]/(4.0 * PI * RAD[1] * RAD[1] * NDENS[1] ) ) );
        const double rLA = std::max( RAD[1] * exp( - 2.5 * sLA * sLA ), 1.5 * LA_R_LOW );
        const AIM::Grid_Aerosol LAAerosol( size_x, size_y, LA_rJ, LA_rE, LA_nDens, rLA, exp(sLA), "lognormal" );

        liquidAerosol = LAAerosol;
    }
    else {
        nBin_LA = 2;
        //dumb hardcoded Grid_Aerosol default constructor
    }
    double pa_r_hig_low = PA_R_HIG/PA_R_LOW;
    nBin_PA = std::floor( 1 + log( pa_r_hig_low * pa_r_hig_low* pa_r_hig_low ) / log( PA_VRAT ) );

    Vector_1D PA_rE( nBin_PA + 1, 0.0 ); /* Bin edges in m */
    Vector_1D PA_rJ( nBin_PA    , 0.0 ); /* Bin center radius in m */
    Vector_1D PA_vJ( nBin_PA    , 0.0 ); /* Bin volume centers in m^3 */

    const double PA_RRAT = cbrt( PA_VRAT );
    PA_rE[0] = PA_R_LOW;
    for ( UInt iBin_PA = 1; iBin_PA < nBin_PA + 1; iBin_PA++ )
        PA_rE[iBin_PA] = PA_rE[iBin_PA-1] * PA_RRAT;                                        /* [m]   */

    for ( UInt iBin_PA = 0; iBin_PA < nBin_PA; iBin_PA++ ) {
        PA_rJ[iBin_PA] = 0.5 * ( PA_rE[iBin_PA] + PA_rE[iBin_PA+1] );                       /* [m]   */
        PA_vJ[iBin_PA] = 4.0 / double(3.0) * PI * \
                         ( PA_rE[iBin_PA] * PA_rE[iBin_PA] * PA_rE[iBin_PA] \
                         + PA_rE[iBin_PA+1] * PA_rE[iBin_PA+1] * PA_rE[iBin_PA+1] ) * 0.5;  /* [m^3] */
    }

    /* TODO: Figure out what PA_nDens and LA_nDens are. They seem to be only used at the original timestep and never updated? -MX */
    PA_nDens = NDENS[0] * 1.00E-06; /* [#/cm^3]      */
    PA_rEff  = RAD[0]   * 1.00E+09; /* [nm]          */
    PA_SAD   = SAD[0]   * 1.00E+06; /* [\mum^2/cm^3] */

    if ( PA_nDens >= 0.0E+00 ) {
        const double expsPA = 1.6;
        const double rPA = std::max( RAD[0] * exp( - 2.5 * log(expsPA) * log(expsPA) ), 1.5 * PA_R_LOW );
        AIM::Grid_Aerosol PAAerosol( size_x, size_y, PA_rJ, PA_rE, PA_nDens, rPA, expsPA, "lognormal" );

        solidAerosol = PAAerosol;
    }
} /* End of Solution::Initialize */

void Solution::processInputBackgroundLine(std::istream &s, Vector_1D &amb_Value, Vector_2D &aer_Value) {
  std::string line;
  UInt i = 0;

  while (std::getline(s, line) && i < nVariables + nAer) {
    if (line.length() > 0 && line != "\r" && line != "\n" && line[0] != '#') {
      std::istringstream iss_line(line);
      if (i < nVariables) {
        iss_line >> amb_Value[i];
      } else if (i >= nVariables && i < nVariables + nAer) {
        iss_line >> aer_Value[i - nVariables][0];
        std::getline(s, line);
        std::istringstream iss_line2(line);
        iss_line2 >> aer_Value[i - nVariables][1];
      }
      i++;
    }
  }
}

// Read default ambient conditions from CMake-generated include file.
const std::string default_ambient =
#include "Defaults/Ambient.hpp"
;

void Solution::readInputBackgroundConditions(const Input& input, Vector_1D& amb_Value, Vector_2D& aer_Value, std::string fileName){
    if (fileName == "=DEFAULT=") {
        // Use default ambient conditions.
        std::istringstream iss(default_ambient);
        processInputBackgroundLine(iss, amb_Value, aer_Value);
    } else {
        std::ifstream file(fileName);
        if (!file) {
            std::cout << "ERROR: In Structure::Initialize: Can't read (" << fileName
                      << ")" << std::endl;
            exit(-1);
        }

        processInputBackgroundLine(file, amb_Value, aer_Value);
    }
}

void Solution::setAmbientConcentrations(const Input& input, Vector_1D& amb_Value){
    if ( input.backgNOx() > 0.0E+00 ) {
        const double NONO2rat = amb_Value[ind_NO]/amb_Value[ind_NO2];
        /* NOx = NO + NO2 = NO2 * ( r + 1 )
         * NO2 = NOx / ( r + 1 );
         * NO  = NOx - NO2; */
        amb_Value[ind_NO2] = input.backgNOx() / ( NONO2rat + 1 );
        amb_Value[ind_NO]  = input.backgNOx() - amb_Value[ind_NO2];

        /* Convert to mixing ratio */
        amb_Value[ind_NO2] /= 1.0E+09;
        amb_Value[ind_NO]  /= 1.0E+09;
    }

    if ( input.backgHNO3() > 0.0E+00 ) {
        /* Convert from ppb to mixing ratio */
        amb_Value[ind_HNO3] = input.backgHNO3() / 1.0E+09;
    }

    if ( input.backgO3() > 0.0E+00 ) {
        /* Convert from ppb to mixing ratio */
        amb_Value[ind_O3] = input.backgO3() / 1.0E+09;
    }

    if ( input.backgCO() > 0.0E+00 ) {
        /* Convert from ppb to mixing ratio */
        amb_Value[ind_CO] = input.backgCO() / 1.0E+09;
    }

    if ( input.backgCH4() > 0.0E+00 ) {
        /* Convert from ppb to mixing ratio */
        amb_Value[ind_CH4] = input.backgCH4() / 1.0E+09;
    }

    if ( input.backgSO2() > 0.0E+00 ) {
        /* Convert from ppb to mixing ratio */
        amb_Value[ind_SO2] = input.backgSO2() / 1.0E+09;
    }

}

void Solution::initializeSpeciesH2O(const Input& input, const OptInput& Input_Opt, Vector_1D& amb_Value, const double airDens, const Meteorology& met){
    UInt actualX = size_x;
    UInt actualY = size_y;
    if ( !Input_Opt.CHEMISTRY_CHEMISTRY ) {
        actualX = 1;
        actualY = 1;

        reducedSize = 1;
    }

    Vector_2D tmpArray( size_y, Vector_1D( size_x, 0.0E+00 ) );
    Vector_2D tmpArray_Reduced( actualY, Vector_1D( actualX, 0.0E+00 ) );

    for ( UInt N = 0; N < NSPECALL; N++ ) {
        if ( ( N == ind_H2O      ) || \
             ( N == ind_H2Omet   ) || \
             ( N == ind_H2Oplume ) || \
             ( N == ind_H2OL     ) || \
             ( N == ind_H2OS     ) ) {
            VectorUtils::set_shape( tmpArray, size_x, size_y, amb_Value[N] * airDens );
            Species.push_back( tmpArray );
        } else {
            VectorUtils::set_shape( tmpArray_Reduced, actualX, actualY, amb_Value[N] * airDens );
            Species.push_back( tmpArray_Reduced );
        }
    }

    if ( Input_Opt.MET_LOADMET ) {
        /* Use meteorological input? */

        //TODO: Fix this insanely wasteful copy. Probably not happening without significant refactoring everything.
        Species[ind_H2Omet] = met.H2O_field();
        /* Update H2O */
        for ( UInt i = 0; i < size_x; i++ ) {
            for ( UInt j = 0; j < size_y; j++ ) {
                Species[ind_H2O][j][i] = Species[ind_H2Omet][j][i] \
                                       + Species[ind_H2Oplume][j][i];
            }
        }
    } else {
        /* Else use user-defined H2O profile */
        double H2Oval = (input.relHumidity_w()/((double) 100.0) * \
                          physFunc::pSat_H2Ol( input.temperature_K() ) / ( kB * input.temperature_K() )) / 1.00E+06;
        for ( UInt i = 0; i < size_x; i++ ) {
            for ( UInt j = 0; j < size_y; j++ ) {
                //H2O[j][i] = H2Oval;
                Species[ind_H2Omet][j][i] = H2Oval;
                /* RH_w = x_H2O * P / Psat_H2Ol(T) = [H2O](#/cm3) * 1E6 * kB * T / Psat_H2Ol(T) */
                Species[ind_H2O][j][i] = Species[ind_H2Omet][j][i] + Species[ind_H2Oplume][j][i];
            }
        }
    }

}

void Solution::setSpeciesValues(Vector_1D& AERFRAC,  Vector_1D& SOLIDFRAC, const Vector_1D& stratData){
    /* Liquid/solid species */
    VectorUtils::set_value( Species[ind_SO4L], (double) AERFRAC[0]                          * stratData[0] );
    VectorUtils::set_value( Species[ind_SO4] , (double) ( 1.0 - AERFRAC[0] )                * stratData[0] );

    AERFRAC[6] = 0.0E+00;
    SOLIDFRAC[6] = 0.0E+00;
    VectorUtils::set_value( Species[ind_H2OL] , (double) AERFRAC[6]                          * stratData[6] );
    VectorUtils::set_value( Species[ind_H2OS] , (double) SOLIDFRAC[6]                        * stratData[6] );
    /* Do not overwrite H2O!! */
    //VectorUtils::set_value( Species[ind_H2O]  , (double) ( 1.0 - AERFRAC[6] - SOLIDFRAC[6] ) * stratData[6] );

    VectorUtils::set_value( Species[ind_HNO3L], (double) AERFRAC[1]                          * stratData[1] );
    VectorUtils::set_value( Species[ind_HNO3S], (double) SOLIDFRAC[1]                        * stratData[1] );
    VectorUtils::set_value( Species[ind_HNO3] , (double) ( 1.0 - AERFRAC[1] - SOLIDFRAC[1] ) * stratData[1] );

    VectorUtils::set_value( Species[ind_HClL] , (double) AERFRAC[2]                          * stratData[2] );
    VectorUtils::set_value( Species[ind_HCl]  , (double) ( 1.0 - AERFRAC[2] )                * stratData[2] );

    VectorUtils::set_value( Species[ind_HOClL], (double) AERFRAC[3]                          * stratData[3] );
    VectorUtils::set_value( Species[ind_HOCl] , (double) ( 1.0 - AERFRAC[3] )                * stratData[3] );

    VectorUtils::set_value( Species[ind_HBrL] , (double) AERFRAC[4]                          * stratData[4] );
    VectorUtils::set_value( Species[ind_HBr]  , (double) ( 1.0 - AERFRAC[4] )                * stratData[4] );

    VectorUtils::set_value( Species[ind_HOBrL], (double) AERFRAC[5]                          * stratData[5] );
    VectorUtils::set_value( Species[ind_HOBr] , (double) ( 1.0 - AERFRAC[5] )                * stratData[5] );

    VectorUtils::set_value( Species[ind_NIT]  , (double) stratData[ 9] );
    VectorUtils::set_value( Species[ind_NAT]  , (double) stratData[10] );
}


void Solution::getData(Vector_1D &varSpeciesArray, const UInt i, const UInt j) {
    for ( UInt N = 0; N < NVAR; N++ ) {
        if ( ( N == ind_H2O      ) || \
                ( N == ind_H2Omet   ) || \
                ( N == ind_H2Oplume ) || \
                ( N == ind_H2OL     ) || \
                ( N == ind_H2OS     ) ) {
            varSpeciesArray[N] = Species[N][j][i];
        }
        else {
            varSpeciesArray[N] = Species[N][0][0];
        }
    }
} /* End of Solution::getData */

/* 
No functions check the return code of this function, temp fix is to return void
*/
void Solution::SpinUp(Vector_1D &amb_Value, const Input &input, const double airDens,
                      const double startTime, Vector_1D &varSpeciesArray, const bool DBG) {
    /* Chemistry timestep
     * DT_CHEM               = 10 mins */
    const double DT_CHEM = 10.0 * 60.0;
    double curr_Time_s   = startTime * 3600.0;

    /* Integrate chemistry from startTime until endTime
     * Make sure than endTime is greater than startTime,
     * if not integrate until next day at the same time */
    double RunUntil = input.emissionTime();

//    /* If emission time corresponds to data from file exit here */
//    if ( startTime == RunUntil )
//        return 1;

//    if ( RunUntil < startTime )
//        RunUntil += 24.0;
    RunUntil += 24.0;

    /* Convert to seconds */
    RunUntil *= 3600.0;

    /* Allocate arrays for KPP */
    int IERR = 0;
    double STEPMIN = (double)0.0;

    double RTOL[NVAR];
    double ATOL[NVAR];

    for( UInt i = 0; i < NVAR; i++ ) {
        RTOL[i] = KPP_RTOLS;
        ATOL[i] = KPP_ATOLS;
    }


    /* Initialize arrays */
    for ( UInt iVar = 0; iVar < NVAR; iVar++ )
        varSpeciesArray[iVar] = amb_Value[iVar] * airDens;

    Vector_1D fixSpeciesArray(NFIX);
    for ( UInt iFix = 0; iFix < NFIX; iFix++ )
        fixSpeciesArray[iFix] = amb_Value[NVAR+iFix] * airDens;

    /* Define sun parameters */
    /* FIXME: We don't need this on the heap. It's a goddamn local variable. */
    SZA sun( input.latitude_deg(), input.emissionDOY() );

    if ( DBG )
        std::cout << "\n Running spin-up from " << curr_Time_s / 3600.0 << " to " << RunUntil / 3600.0 << " [hr]\n";

    while ( curr_Time_s < RunUntil ) {

        /* Compute the cosize of solar zenith angle midway through the integration step */
        sun.Update( curr_Time_s + DT_CHEM/2 );

        for ( UInt iPhotol = 0; iPhotol < NPHOTOL; iPhotol++ )
            PHOTOL[iPhotol] = 0.0E+00;

        if ( sun.CSZA > 0.0E+00 )
            Update_JRates( PHOTOL, sun.CSZA );

        if ( DBG ) {
            std::cout << "\n DEBUG : (In SpinUp)\n";
            for ( UInt iPhotol = 0; iPhotol < NPHOTOL; iPhotol++ )
                std::cout << "         PHOTOL[" << iPhotol << "] = " << PHOTOL[iPhotol] << "\n";
        }

        /* Update reaction rates */
        for ( UInt iReact = 0; iReact < NREACT; iReact++ )
            RCONST[iReact] = 0.0E+00;

        Update_RCONST( input.temperature_K(), input.pressure_Pa(), airDens, varSpeciesArray[ind_H2O] );

        /* ~~~~~~~~~~~~~~~~~~~~~~~~ */
        /* ~~~~~ Integration ~~~~~~ */
        /* ~~~~~~~~~~~~~~~~~~~~~~~~ */

        IERR = INTEGRATE(varSpeciesArray.data(), fixSpeciesArray.data(),
                         curr_Time_s, curr_Time_s + DT_CHEM, ATOL, RTOL, STEPMIN);

        if ( IERR < 0 ) {
            /* Integration failed */

            std::cout << " SpinUp Integration failed";
            #ifdef OMP
                std::cout << " on " << omp_get_thread_num();
            #endif /* OMP */
            std::cout << " at time t = " << curr_Time_s/3600.0 << "\n";

            if ( DBG ) {
                std::cout << " ~~~ Printing reaction rates:\n";
                for ( UInt iReact = 0; iReact < NREACT; iReact++ ) {
                    std::cout << "Reaction " << iReact << ": " << RCONST[iReact] << " [molec/cm^3/s]\n";
                }
                std::cout << " ~~~ Printing concentrations:\n";
                for ( UInt iSpec = 0; iSpec < NVAR; iSpec++ ) {
                    std::cout << "Species " << iSpec << ": " << varSpeciesArray[iSpec]/airDens*1.0E+09 << " [ppb]\n";
                }
            }

            // return KPP_FAIL;
        }

        curr_Time_s += DT_CHEM;

    }

    for ( UInt iVar = 0; iVar < NVAR; iVar++ )
        amb_Value[iVar] = varSpeciesArray[iVar] / airDens;

    // return IERR;

} /* End of Solution::SpinUp */

}
