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

#include "KPP/KPP.hpp"
#include "KPP/KPP_Parameters.h"
#include "Core/LiquidAer.hpp"
#include "Core/Parameters.hpp"
#include "Core/SZA.hpp"
#include "Core/Structure.hpp"
#include "Util/PhysConstant.hpp"

Solution::Solution(const OptInput& optInput) : \
        liquidAerosol( ), 
        solidAerosol( ),
        LA_Kernel( ), 
        PA_Kernel( ),
        nVariables( NSPEC ), 
        nAer( N_AER ), 
        size_x( optInput.ADV_GRID_NX ), 
        size_y( optInput.ADV_GRID_NY ),
        reducedSize( 0 )
{

    /* Constructor */

} /* End of Solution::Solution */

Solution::~Solution()
{

    /* Destructor */

} /* End of Solution::~Solution */

/* FIXME: This function is not unique to this specific class. Move it to some utility functions namespace or sth idk */
void Solution::Clear( Vector_2D& vector_2D )
{

    for ( UInt i = 0; i < vector_2D.size(); i++ ) {
        vector_2D[i].clear();
    }
    vector_2D.clear();

} /* End of Solution::Clear */

/* FIXME: see above comment on clear() */
void Solution::SetShape( Vector_2D& vector_2D, \
                         const UInt n_x,       \
                         const UInt n_y,       \
                         const double value )
{

    Clear( vector_2D );

    /* Dimensions are transposed! */
    for ( UInt i = 0; i < n_y; i++ ) {
        vector_2D.push_back( Vector_1D( n_x, value ) );
    }

} /* End of Solution::SetShape */

/* FIXME: see above comment on clear() */
void Solution::SetToValue( Vector_2D& vector_2D, \
                           const double value )
{

    for ( UInt i = 0; i < vector_2D.size(); i++ ) {
        for ( UInt j = 0; j < vector_2D[0].size(); j++ ) {
            vector_2D[i][j] = value;
        }
    }

} /* End of Solution::SetToValue */

/* FIXME: see above comment on clear() */
void Solution::Print( const Vector_2D& vector_2D, \
                      const UInt i_max,           \
                      const UInt j_max ) const
{

    for ( UInt i = 0; i < i_max; i++ ) {
        for ( UInt j = 0; j < j_max; j++ ) {
            std::cout << vector_2D[i][j] << ", ";
        }
        std::cout << std::endl;
    }

} /* End of Solution::Print */

void Solution::Initialize( char const *fileName,      \
                           const Input &input,        \
                           const double airDens,  \
                           const Meteorology &met,    \
                           const OptInput &Input_Opt, \
                           double* varSpeciesArray, double* fixSpeciesArray, 
                           const bool DBG )
{

    Vector_1D amb_Value(NSPECALL, 0.0);
    Vector_2D aer_Value(nAer, Vector_1D(2, 0.0));

    /* Read input background conditions */
    readInputBackgroundConditions(input, amb_Value, aer_Value, fileName);

    const double AMBIENT_VALID_TIME = 8.0; //hours
    SpinUp( amb_Value, input, airDens, AMBIENT_VALID_TIME, varSpeciesArray, fixSpeciesArray );

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
    SetShape( sootDens , size_x, size_y, (double) aer_Value[  0][0] );
    SetShape( sootRadi , size_x, size_y, (double) aer_Value[  0][1] );
    SetShape( sootArea , size_x, size_y, (double) 4.0 / double(3.0) * physConst::PI * aer_Value[  0][0] * aer_Value[  0][1] * aer_Value[  0][1] * aer_Value[  0][1] );


    nBin_LA = std::floor( 1 + log( pow( (LA_R_HIG/LA_R_LOW), 3.0 ) ) / log( LA_VRAT ) );

    Vector_1D LA_rE( nBin_LA + 1, 0.0 ); /* Bin edges in m */
    Vector_1D LA_rJ( nBin_LA    , 0.0 ); /* Bin center radius in m */
    Vector_1D LA_vJ( nBin_LA    , 0.0 ); /* Bin volume centers in m^3 */

    const double LA_RRAT = pow( LA_VRAT, 1.0 / double(3.0) );
    LA_rE[0] = LA_R_LOW;
    for ( UInt iBin_LA = 1; iBin_LA < nBin_LA + 1; iBin_LA++ )                             
        LA_rE[iBin_LA] = LA_rE[iBin_LA-1] * LA_RRAT; /* [m] */

    for ( UInt iBin_LA = 0; iBin_LA < nBin_LA; iBin_LA++ ) {
        LA_rJ[iBin_LA] = 0.5 * ( LA_rE[iBin_LA] + LA_rE[iBin_LA+1] );                       /* [m] */
        LA_vJ[iBin_LA] = 4.0 / double(3.0) * physConst::PI * \
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

        const double sLA = sqrt( - 1.0 / (3.0) * log(SAD[1]/(4.0 * physConst::PI * RAD[1] * RAD[1] * NDENS[1] ) ) );
        const double rLA = std::max( RAD[1] * exp( - 2.5 * sLA * sLA ), 1.5 * LA_R_LOW );
        const AIM::Grid_Aerosol LAAerosol( size_x, size_y, LA_rJ, LA_rE, LA_nDens, rLA, exp(sLA), "lognormal" );

        liquidAerosol = LAAerosol;
    }
    else {
        nBin_LA = 2;
        //dumb hardcoded Grid_Aerosol default constructor
    }
    const AIM::Coagulation kernel1( "liquid", LA_rJ, LA_vJ, physConst::RHO_SULF, \
                                    input.temperature_K(), input.pressure_Pa() );

    LA_Kernel = kernel1;

    nBin_PA = std::floor( 1 + log( pow( (PA_R_HIG/PA_R_LOW), 3.0 ) ) / log( PA_VRAT ) );

    Vector_1D PA_rE( nBin_PA + 1, 0.0 ); /* Bin edges in m */
    Vector_1D PA_rJ( nBin_PA    , 0.0 ); /* Bin center radius in m */
    Vector_1D PA_vJ( nBin_PA    , 0.0 ); /* Bin volume centers in m^3 */

    const double PA_RRAT = pow( PA_VRAT, 1.0 / double(3.0) );
    PA_rE[0] = PA_R_LOW;
    for ( UInt iBin_PA = 1; iBin_PA < nBin_PA + 1; iBin_PA++ )
        PA_rE[iBin_PA] = PA_rE[iBin_PA-1] * PA_RRAT;                                        /* [m]   */

    for ( UInt iBin_PA = 0; iBin_PA < nBin_PA; iBin_PA++ ) {
        PA_rJ[iBin_PA] = 0.5 * ( PA_rE[iBin_PA] + PA_rE[iBin_PA+1] );                       /* [m]   */
        PA_vJ[iBin_PA] = 4.0 / double(3.0) * physConst::PI * \
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

    const AIM::Coagulation kernel2( "ice", PA_rJ, PA_vJ, physConst::RHO_ICE, \
                                    input.temperature_K(), input.pressure_Pa() );

    PA_Kernel = kernel2;

} /* End of Solution::Initialize */

void Solution::readInputBackgroundConditions(const Input& input, Vector_1D& amb_Value, Vector_2D& aer_Value, const char* fileName){
    std::ifstream file;
    file.open( fileName );
    if ( file.is_open() ) {
        std::string line;
        UInt i = 0;

        while ( ( std::getline( file, line ) ) && ( i < nVariables + nAer ) ) {
            if ( ( line.length() > 0 ) && ( line != "\r" ) && ( line != "\n" ) && ( line[0] != '#' ) ) {
                std::istringstream iss(line);
                if ( i < nVariables ) {
                    iss >> amb_Value[i];
                }
                else if ( ( i >= nVariables ) && ( i < nVariables + nAer ) ) {
                    iss >> aer_Value[i - nVariables][0];
                    std::getline( file, line );
                    std::istringstream iss(line);
                    iss >> aer_Value[i - nVariables][1];
                }
                i++;
            }
        }
        file.close();
    }
    else {
        std::string const currFunc("Structure::Initialize");
        std::cout << "ERROR: In " << currFunc << ": Can't read (" << fileName << ")" << std::endl;
        exit(-1);
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
            SetShape( tmpArray, size_x, size_y, amb_Value[N] * airDens );
            Species.push_back( tmpArray );
        } else {
            SetShape( tmpArray_Reduced, actualX, actualY, amb_Value[N] * airDens );
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
                          physFunc::pSat_H2Ol( input.temperature_K() ) / ( physConst::kB * input.temperature_K() )) / 1.00E+06;
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
    SetToValue( Species[ind_SO4L], (double) AERFRAC[0]                          * stratData[0] );
    SetToValue( Species[ind_SO4] , (double) ( 1.0 - AERFRAC[0] )                * stratData[0] );

    AERFRAC[6] = 0.0E+00;
    SOLIDFRAC[6] = 0.0E+00;
    SetToValue( Species[ind_H2OL] , (double) AERFRAC[6]                          * stratData[6] );
    SetToValue( Species[ind_H2OS] , (double) SOLIDFRAC[6]                        * stratData[6] );
    /* Do not overwrite H2O!! */
    //SetToValue( Species[ind_H2O]  , (double) ( 1.0 - AERFRAC[6] - SOLIDFRAC[6] ) * stratData[6] );

    SetToValue( Species[ind_HNO3L], (double) AERFRAC[1]                          * stratData[1] );
    SetToValue( Species[ind_HNO3S], (double) SOLIDFRAC[1]                        * stratData[1] );
    SetToValue( Species[ind_HNO3] , (double) ( 1.0 - AERFRAC[1] - SOLIDFRAC[1] ) * stratData[1] );

    SetToValue( Species[ind_HClL] , (double) AERFRAC[2]                          * stratData[2] );
    SetToValue( Species[ind_HCl]  , (double) ( 1.0 - AERFRAC[2] )                * stratData[2] );

    SetToValue( Species[ind_HOClL], (double) AERFRAC[3]                          * stratData[3] );
    SetToValue( Species[ind_HOCl] , (double) ( 1.0 - AERFRAC[3] )                * stratData[3] );

    SetToValue( Species[ind_HBrL] , (double) AERFRAC[4]                          * stratData[4] );
    SetToValue( Species[ind_HBr]  , (double) ( 1.0 - AERFRAC[4] )                * stratData[4] );

    SetToValue( Species[ind_HOBrL], (double) AERFRAC[5]                          * stratData[5] );
    SetToValue( Species[ind_HOBr] , (double) ( 1.0 - AERFRAC[5] )                * stratData[5] );

    SetToValue( Species[ind_NIT]  , (double) stratData[ 9] );
    SetToValue( Species[ind_NAT]  , (double) stratData[10] );
}


void Solution::getData(  double* varSpeciesArray, double* fixSpeciesArray, const UInt i, \
                        const UInt j, \
	                const bool CHEMISTRY )
{

    for ( UInt N = 0; N < NVAR; N++ ) {
	if ( CHEMISTRY )
            varSpeciesArray[N] = Species[N][j][i];
	else {
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
    }

    for ( UInt N = 0; N < NFIX; N++ ) {
	if ( CHEMISTRY )
            fixSpeciesArray[N] = Species[N+NVAR][j][i];
	else {
            if ( ( N == ind_H2O      ) || \
                 ( N == ind_H2Omet   ) || \
                 ( N == ind_H2Oplume ) || \
                 ( N == ind_H2OL     ) || \
                 ( N == ind_H2OS     ) ) {
                fixSpeciesArray[N] = Species[N+NVAR][j][i];
            }
	    else {
                fixSpeciesArray[N] = Species[N+NVAR][0][0];
	    }
	
	}

    }

} /* End of Solution::getData */

void Solution::applyData( const double* varSpeciesArray, const UInt i, \
                          const UInt j )
{

    for ( UInt N = 0; N < NVAR; N++ )
        Species[N][j][i] = varSpeciesArray[N];

} /* End of Solution::applyData */

void Solution::applyRing( const double* varSpeciesArray,
                          double tempArray[],        \
                          const Vector_2Dui &mapIndices, \
                          const UInt iRing )
{

    UInt iNx = 0;
    UInt jNy = 0;
    UInt N   = 0;

    for ( jNy = 0; jNy < size_y; jNy++ ) {
        for ( iNx = 0; iNx < size_x; iNx++ ) {
            if ( mapIndices[jNy][iNx] == iRing ) {

                for ( N = 0; N < NVAR; N++ ) {
                    if ( ( N == ind_N ) || ( N == ind_O ) || ( N == ind_O1D ) ) {
                        /* Special handlings! */
                        Species[ind_N][jNy][iNx]     = varSpeciesArray[ind_N];
                        Species[ind_O][jNy][iNx]     = varSpeciesArray[ind_O];
                        Species[ind_O1D][jNy][iNx]   = varSpeciesArray[ind_O1D];
                    } else
                        Species[N][jNy][iNx] *= varSpeciesArray[N] / tempArray[N];
                }

            }
        }
    }

} /* End of Solution::applyRing */

void Solution::applyAmbient( const double* varSpeciesArray, const Vector_2Dui &mapIndices, \
                             const UInt iRing )
{

    UInt iNx = 0;
    UInt jNy = 0;
    UInt N   = 0;

    for ( jNy = 0; jNy < size_y; jNy++ ) {
        for ( iNx = 0; iNx < size_x; iNx++ ) {
            if ( mapIndices[jNy][iNx] == iRing ) {
                for ( N = 0; N < NVAR; N++ )
                    Species[N][jNy][iNx] = varSpeciesArray[N];
            }
        }
    }

} /* End of Solution::applyAmbient */

/* TODO !!!! : Write an addEmission that doesn't have to deal with the archaic RINGS framework.
Though first we need to figure out what this function does... hahaha */
void Solution::addEmission( const Emission &EI, const Aircraft &AC,        \
                            const Mesh &m,                                 \
                            bool halfRing,                                 \
                            const double temperature, bool set2Saturation, \
                            AIM::Aerosol liqAer, AIM::Aerosol iceAer,      \
                            const double Soot_Den,                         \
                            const Meteorology &met, const double areaPlume )
{
    /*  TODO: Release as Gaussian instead of top-hat? 
        Also, there has to be a better way add to liquidAerosol and solidAerosol than passing by value and wasting memory,
        but will need proper unit tests to bother with.*/

    UInt innerRing;
    UInt iNx, jNy;
    double w;
    double E_CO2, E_H2O, E_NO, E_NO2, E_HNO2, E_SO2, E_CO, E_CH4, E_C2H6, \
               E_PRPE, E_ALK4, E_CH2O, E_ALD2, E_GLYX, E_MGLY;
    // double E_Soot;
    const double rad = EI.getSootRad();
    const double fuelPerDist = AC.FuelFlow() / AC.VFlight();

    /* Unit check:  [kg/m]   =   [kg fuel/s] /     [m/s] */
    E_CO2  = EI.getCO2()  / ( MW_CO2  * 1.0E+03 ) * fuelPerDist * physConst::Na;
    /*     = [g(CO2)/kg f]/ ( [kg/mol]* [g/kg]  ) * [kg fuel/m] * [molec/mol]
     *     = [molec/m]
     */
    E_NO   = EI.getNO()   / ( MW_NO   * 1.0E+03 ) * fuelPerDist * physConst::Na;
    /*     = [g(NO)/kg f] / ( g(NO)/mol         ) * [kg fuel/m] * [molec/mol]
     *     = [molec/m] */
    E_NO2  = EI.getNO2()  / ( MW_NO2  * 1.0E+03 ) * fuelPerDist * physConst::Na;
    E_HNO2 = EI.getHNO2() / ( MW_HNO2 * 1.0E+03 ) * fuelPerDist * physConst::Na;
    E_CO   = EI.getCO()   / ( MW_CO   * 1.0E+03 ) * fuelPerDist * physConst::Na;
    E_CH4  = EI.getCH4()  / ( MW_CH4  * 1.0E+03 ) * fuelPerDist * physConst::Na;
    E_C2H6 = EI.getC2H6() / ( MW_C2H6 * 1.0E+03 ) * fuelPerDist * physConst::Na;
    E_PRPE = EI.getPRPE() / ( MW_PRPE * 1.0E+03 ) * fuelPerDist * physConst::Na;
    E_ALK4 = EI.getALK4() / ( MW_ALK4 * 1.0E+03 ) * fuelPerDist * physConst::Na;
    E_CH2O = EI.getCH2O() / ( MW_CH2O * 1.0E+03 ) * fuelPerDist * physConst::Na;
    E_ALD2 = EI.getALD2() / ( MW_ALD2 * 1.0E+03 ) * fuelPerDist * physConst::Na;
    E_GLYX = EI.getGLYX() / ( MW_GLYX * 1.0E+03 ) * fuelPerDist * physConst::Na;
    E_MGLY = EI.getMGLY() / ( MW_MGLY * 1.0E+03 ) * fuelPerDist * physConst::Na;
    if ( !set2Saturation ) {
        E_H2O  = EI.getH2O()  / ( MW_H2O  * 1.0E+03 ) * fuelPerDist * physConst::Na;
    }
    E_SO2  = ( 1.0 - SO2TOSO4 ) * \
             EI.getSO2()  / ( MW_SO2  * 1.0E+03 ) * fuelPerDist * physConst::Na;

    // E_Soot = EI.getSoot() / ( 4.0 / 3.0 * physConst::PI * physConst::RHO_SOOT * 1.00E+03 * rad * rad * rad ) * fuelPerDist;
    /*     = [g_soot/kg_fuel]/ (                        * [kg_soot/m^3]       * [g/kg]   * [m^3]           ) * [kg_fuel/m]
     *     = [part/m]
     */

    Vector_2D cellAreas  = m.areas();
    Vector_3D weights    = m.weights;
    Vector_1Dui nCellMap = m.nMap();
    UInt nCell   = 0;

    /* Apply area scaling for aerosols before releasing into the grid
     * *Aer are initially in #/cm^3
     * Scaling them the initial plume area such that they are now in #*m^2/cm^3
     * When applying them to the grid, they are scaled by the total ring area
     * in Grid_Aerosol::addPDF */
    liqAer.scalePdf( areaPlume );
    iceAer.scalePdf( areaPlume );


    if ( !halfRing ) {
        /* Full rings */
        innerRing = 0;
        nCell     = nCellMap[innerRing];
        for ( jNy = 0; jNy < size_y; jNy++ ) {
            for ( iNx = 0; iNx < size_x; iNx++ ) {

                Species[ind_H2Omet][jNy][iNx] = met.H2O(jNy, iNx);

                w = weights[innerRing][jNy][iNx];

                /* Initially weights are either 0 or 1 */
                if ( w != 0.0E+00 ) {

                    if ( !reducedSize ) {
                        /* Emissions will only be added if chemistry is turned on! */

                        /* Each of the following arrays are in [# / cm^3], where
                         * # can either be expressed in molecules for gaseous
                         * species or in number of particles for aerosols */
                        Species[ind_CO2][jNy][iNx]  += ( E_CO2  * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                        Species[ind_NO][jNy][iNx]   += ( E_NO   * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                        Species[ind_NO2][jNy][iNx]  += ( E_NO2  * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                        Species[ind_HNO2][jNy][iNx] += ( E_HNO2 * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                        Species[ind_CO][jNy][iNx]   += ( E_CO   * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                        Species[ind_CH4][jNy][iNx]  += ( E_CH4  * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                        Species[ind_C2H6][jNy][iNx] += ( E_C2H6 * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                        Species[ind_PRPE][jNy][iNx] += ( E_PRPE * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                        Species[ind_ALK4][jNy][iNx] += ( E_ALK4 * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                        Species[ind_CH2O][jNy][iNx] += ( E_CH2O * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                        Species[ind_ALD2][jNy][iNx] += ( E_ALD2 * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                        Species[ind_GLYX][jNy][iNx] += ( E_GLYX * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                        Species[ind_MGLY][jNy][iNx] += ( E_MGLY * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                        Species[ind_SO2][jNy][iNx]  += ( E_SO2  * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                    }

                    if ( set2Saturation ) {
                        /* If supersaturated, then set water vapor to saturation and no bare soot particles
                         * as they are all covered with ice */
                        Species[ind_H2Oplume][jNy][iNx] += physFunc::pSat_H2Os( met.temp(jNy, iNx) ) / ( physConst::kB * met.temp(jNy, iNx) * 1.00E+06 ) - Species[ind_H2Omet][jNy][iNx]; /* [molec / cm^3] */
                    } else {
                        /* If subsaturated, then emit water and soot */
                        Species[ind_H2Oplume][jNy][iNx]      += ( E_H2O * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) ); /* [molec / cm^3] */
                        sootDens[jNy][iNx] += ( Soot_Den * areaPlume / ( nCell * cellAreas[jNy][iNx] ) );
                        sootRadi[jNy][iNx] = rad;
                        sootArea[jNy][iNx] = 4.0 * physConst::PI * rad * rad * sootDens[jNy][iNx];
                    }

                }

            }

        }

        if ( ( std::isfinite(iceAer.Moment()) ) && ( iceAer.Moment() > 0.0E+00 ) )
            solidAerosol.addPDF( iceAer, weights[innerRing], cellAreas, nCell );
        if ( ( std::isfinite(liqAer.Moment()) ) && ( liqAer.Moment() > 0.0E+00 ) )
            liquidAerosol.addPDF( liqAer, weights[innerRing], cellAreas, nCell );


    } else {
        /* Half rings */

        for ( innerRing = 0; innerRing <= 1; innerRing++ ) {
            /* Concentrations should be identical whether it is half-rings or
             * full rings. Therefore, make sure that nCell is doubled when
             * using half-rings! */
            nCell     = 2.0 * nCellMap[innerRing];
            for ( jNy = 0; jNy < size_y; jNy++ ) {
                for ( iNx = 0; iNx < size_x; iNx++ ) {

                    Species[ind_H2Omet][jNy][iNx] = met.H2O(jNy, iNx);

                    w = weights[innerRing][jNy][iNx];

                    /* Initially weights are either 0 or 1 */
                    if ( w != 0.0E+00 ) {
                        if ( !reducedSize ) {
                            /* Emissions will only be added if chemistry is turned on! */

                            /* Each of the following arrays are in [# / cm^3], where
                             * # can either be expressed in molecules for gaseous
                             * species or in number of particles for aerosols */
                            Species[ind_CO2][jNy][iNx]  += ( E_CO2  * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                            Species[ind_NO][jNy][iNx]   += ( E_NO   * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                            Species[ind_NO2][jNy][iNx]  += ( E_NO2  * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                            Species[ind_HNO2][jNy][iNx] += ( E_HNO2 * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                            Species[ind_CO][jNy][iNx]   += ( E_CO   * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                            Species[ind_CH4][jNy][iNx]  += ( E_CH4  * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                            Species[ind_C2H6][jNy][iNx] += ( E_C2H6 * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                            Species[ind_PRPE][jNy][iNx] += ( E_PRPE * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                            Species[ind_ALK4][jNy][iNx] += ( E_ALK4 * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                            Species[ind_CH2O][jNy][iNx] += ( E_CH2O * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                            Species[ind_ALD2][jNy][iNx] += ( E_ALD2 * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                            Species[ind_GLYX][jNy][iNx] += ( E_GLYX * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                            Species[ind_MGLY][jNy][iNx] += ( E_MGLY * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                            Species[ind_SO2][jNy][iNx]  += ( E_SO2  * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) );
                        }

                        if ( set2Saturation ) {
                            /* If supersaturated, then set water vapor to saturation and no
                             * bare soot particles as they are all covered with ice */
                            Species[ind_H2Oplume][jNy][iNx] += physFunc::pSat_H2Os( met.temp(jNy, iNx) ) / ( physConst::kB * met.temp(jNy, iNx) * 1.00E+06 ) - Species[ind_H2Omet][jNy][iNx]; /* [molec / cm^3] */
                        } else {
                            /* If subsaturated, then emit water and soot */
                            Species[ind_H2Oplume][jNy][iNx]      += ( E_H2O * 1.0E-06 / ( nCell * cellAreas[jNy][iNx] ) ); /* [molec / cm^3] */
                            sootDens[jNy][iNx] += ( Soot_Den * areaPlume / ( nCell * cellAreas[jNy][iNx] ) );
                            sootRadi[jNy][iNx] = rad;
                            sootArea[jNy][iNx] = 4.0 * physConst::PI * rad * rad * sootDens[jNy][iNx];
                        }

                    }
                }
            }
            
	    if ( ( std::isfinite(iceAer.Moment()) ) && ( iceAer.Moment() > 0.0E+00 ) )
                solidAerosol.addPDF( iceAer, weights[innerRing], cellAreas, nCell );
            if ( ( std::isfinite(liqAer.Moment()) ) && ( liqAer.Moment() > 0.0E+00 ) )
                liquidAerosol.addPDF( liqAer, weights[innerRing], cellAreas, nCell );

        }

    }

    for ( jNy = 0; jNy < size_y; jNy++ ) {
        for ( iNx = 0; iNx < size_x; iNx++ ) {
            Species[ind_H2O][jNy][iNx] = Species[ind_H2Omet][jNy][iNx] + Species[ind_H2Oplume][jNy][iNx];
        }
    }

} /* End of Solution::addEmission */

Vector_1D Solution::getAmbient() const
{

    Vector_1D ambVector( NSPECREACT, 0.0 );

    for ( UInt iVar = 0; iVar < NSPECREACT; iVar++ )
        ambVector[iVar] = Species[iVar][0][0];

    return ambVector;

} /* End of Solution::getAmbient */

Vector_1D Solution::getLiqSpecies( ) const
{

    Vector_1D liqAerVector( 9, 0.0 );

    liqAerVector[0] = Species[ind_SO4L][0][0];
    liqAerVector[1] = Species[ind_H2OL][0][0];
    liqAerVector[2] = Species[ind_HNO3L][0][0];
    liqAerVector[3] = Species[ind_HClL][0][0];
    liqAerVector[4] = Species[ind_HOClL][0][0];
    liqAerVector[5] = Species[ind_HBrL][0][0];
    liqAerVector[6] = Species[ind_HOBrL][0][0];
    liqAerVector[7] = Species[ind_H2OS][0][0];
    liqAerVector[8] = Species[ind_HNO3S][0][0];

    return liqAerVector;

} /* End of Solution::getLiqSpecies */

Vector_2D Solution::getAerosol( ) const
{

    Vector_2D aerVector( nAer, Vector_1D( 3, 0.0 ) );

    aerVector[  0][0] = sootDens[0][0];
    aerVector[  0][1] = sootRadi[0][0];
    aerVector[  0][2] = sootArea[0][0];
    aerVector[  1][0] = PA_nDens;
    aerVector[  1][1] = PA_rEff;
    aerVector[  1][2] = PA_SAD;
    aerVector[  2][0] = LA_nDens;
    aerVector[  2][1] = LA_rEff;
    aerVector[  2][2] = LA_SAD;

    return aerVector;

} /* End of Solution::getAerosol */

Vector_1D Solution::getAerosolDens( ) const
{

    Vector_1D aerVector( nAer, 0.0 );

    aerVector[  0] = sootDens[0][0];
    aerVector[  1] = PA_nDens;
    aerVector[  2] = LA_nDens;

    return aerVector;

} /* End of Solution::getAerosolDens */

Vector_1D Solution::getAerosolRadi( ) const
{

    Vector_1D aerVector( nAer, 0.0 );

    aerVector[  0] = sootRadi[0][0];
    aerVector[  1] = PA_rEff;
    aerVector[  2] = LA_rEff;

    return aerVector;

} /* End of Solution::getAerosolRadi */

Vector_1D Solution::getAerosolArea( ) const
{

    Vector_1D aerVector( nAer, 0.0 );

    aerVector[  0] = sootArea[0][0];
    aerVector[  1] = PA_SAD;
    aerVector[  2] = LA_SAD;

    return aerVector;

} /* End of Solution::getAerosolArea */

void Solution::getAerosolProp( double ( &radi )[4], \
                               double ( &area )[4], \
                               double &IWC,         \
                               const Vector_2D &weights ) const
{

    UInt jNy, iNx;

    Vector_1D aerosolProp( 4, 0.0E+00 );
    double totalWeight = 0.0E+00;

    for ( jNy = 0; jNy < size_y; jNy++ ) {
        for ( iNx = 0; iNx < size_x; iNx++ )
            totalWeight += weights[jNy][iNx];
    }

    /* Compute aerosol microphysical properties for ice/NAT */
    aerosolProp = solidAerosol.Average( weights, totalWeight );

    radi[0] = aerosolProp[1];                      /* [m]        */
    area[0] = aerosolProp[2];                      /* [m^2/cm^3] */
    IWC     = physConst::RHO_ICE * aerosolProp[3]; /* [kg/cm^3] */

    /* Compute aerosol microphysical properties for stratospheric liquid aerosols */
    aerosolProp = liquidAerosol.Average( weights, totalWeight );

    radi[1] = aerosolProp[1]; /* [m]        */
    area[1] = aerosolProp[2]; /* [m^2/cm^3] */

    /* Compute aerosol microphysical properties for tropospheric sulfates (near ground) */

    radi[2] = 0.0E+00;
    area[2] = 0.0E+00;

    /* Compute aerosol microphysical properties for BC particles */

    radi[3] = 0.0E+00;
    area[3] = 0.0E+00;

    for ( jNy = 0; jNy < size_y; jNy++ ) {
        for ( iNx = 0; iNx < size_x; iNx++ ) {
            area[3] += sootDens[jNy][iNx] * 4.0 * physConst::PI * \
                       sootRadi[jNy][iNx] * sootRadi[jNy][iNx] *  \
                       weights[jNy][iNx] / totalWeight;
            radi[3] += sootRadi[jNy][iNx] * \
                       weights[jNy][iNx] / totalWeight;
        }
    }

} /* End of Solution::getAerosolProp */

/* 
No functions check the return code of this function, temp fix is to return void
*/
void Solution::SpinUp( Vector_1D &amb_Value,       \
                      const Input &input,         \
                      const double airDens,   \
                      const double startTime, \
                      double* varSpeciesArray, double* fixSpeciesArray, const bool DBG )
{

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

        IERR = INTEGRATE( varSpeciesArray, fixSpeciesArray, curr_Time_s, curr_Time_s + DT_CHEM, \
                          ATOL, RTOL, STEPMIN );

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

void Solution::Debug( const double airDens )
{

    UInt iNx, jNy;
    iNx = 0;
    jNy = 0;

    std::cout << std::endl;
    std::cout << "**** Input Debugger ****" << std::endl;
    std::cout << "Background concentrations: " << std::endl;
    std::cout << std::endl;
    std::cout << "   Species " << "    Value" << " Units " << std::endl;

    std::cout << std::setw(9);
    std::cout << "CO2" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_CO2][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "PPN" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_PPN][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "BrNO2" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_BrNO2][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "IEPOX" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_IEPOX][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "PMNN" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_PMNN][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "N2O" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_N2O][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "N" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_N][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "PAN" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_PAN][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ALK4" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_ALK4][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MAP" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_MAP][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MPN" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_MPN][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "Cl2O2" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_Cl2O2][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ETP" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_ETP][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "HNO2" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_HNO2][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "C3H8" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_C3H8][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "RA3P" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_RA3P][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "RB3P" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_RB3P][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "OClO" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_OClO][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ClNO2" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_ClNO2][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ISOP" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_ISOP][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "HNO4" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_HNO4][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MAOP" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_MAOP][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MP" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_MP][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ClOO" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_ClOO][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "RP" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_RP][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "BrCl" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_BrCl][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "PP" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_PP][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "PRPN" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_PRPN][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "SO4" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_SO4][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "Br2" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_Br2][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ETHLN" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_ETHLN][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MVKN" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_MVKN][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "R4P" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_R4P][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "C2H6" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_C2H6][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "RIP" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_RIP][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "VRP" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_VRP][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ATOOH" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_ATOOH][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "IAP" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_IAP][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "DHMOB" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_DHMOB][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MOBA" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_MOBA][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MRP" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_MRP][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "N2O5" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_N2O5][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ISNOHOO" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_ISNOHOO][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ISNP" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_ISNP][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ISOPNB" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_ISOPNB][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "IEPOXOO" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_IEPOXOO][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MACRNO2" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_MACRNO2][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ROH" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_ROH][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MOBAOO" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_MOBAOO][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "DIBOO" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_DIBOO][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "PMN" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_PMN][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ISNOOB" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_ISNOOB][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "INPN" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_INPN][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "H" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_H][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "BrNO3" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_BrNO3][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "PRPE" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_PRPE][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MVKOO" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_MVKOO][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "Cl2" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_Cl2][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ISOPND" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_ISOPND][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "HOBr" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_HOBr][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "A3O2" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_A3O2][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "PROPNN" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_PROPNN][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "GLYX" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_GLYX][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MAOPO2" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_MAOPO2][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "CH4" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_CH4][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "GAOO" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_GAOO][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "B3O2" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_B3O2][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ACET" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_ACET][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MACRN" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_MACRN][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "CH2OO" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_CH2OO][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MGLYOO" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_MGLYOO][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "VRO2" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_VRO2][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MGLOO" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_MGLOO][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MACROO" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_MACROO][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "PO2" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_PO2][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "CH3CHOO" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_CH3CHOO][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MAN2" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_MAN2][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ISNOOA" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_ISNOOA][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "H2O2" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_H2O2][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "PRN1" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_PRN1][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ETO2" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_ETO2][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "KO2" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_KO2][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "RCO3" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_RCO3][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "HC5OO" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_HC5OO][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "GLYC" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_GLYC][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ClNO3" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_ClNO3][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "RIO2" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_RIO2][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "R4N1" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_R4N1][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "HOCl" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_HOCl][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ATO2" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_ATO2][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "HNO3" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_HNO3][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ISN1" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_ISN1][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MAO3" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_MAO3][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MRO2" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_MRO2][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "INO2" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_INO2][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "HAC" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_HAC][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "HC5" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_HC5][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MGLY" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_MGLY][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ISOPNBO2" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_ISOPNBO2][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ISOPNDO2" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_ISOPNDO2][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "R4O2" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_R4O2][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "R4N2" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_R4N2][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "BrO" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_BrO][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "RCHO" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_RCHO][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MEK" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_MEK][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ClO" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_ClO][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MACR" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_MACR][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "SO2" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_SO2][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MVK" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_MVK][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ALD2" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_ALD2][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MCO3" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_MCO3][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "CH2O" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_CH2O][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "H2O" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_H2O][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "Br" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_Br][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "NO" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_NO][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "NO3" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_NO3][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "Cl" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_Cl][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "O" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_O][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "O1D" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_O1D][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "O3" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_O3][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "HO2" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_HO2][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "NO2" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_NO2][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "OH" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_OH][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "HBr" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_HBr][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "HCl" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_HCl][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "CO" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_CO][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MO2" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_MO2][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "ACTA" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_ACTA][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "EOH" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_EOH][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "H2" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_H2][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "HCOOH" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_HCOOH][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "MOH" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_MOH][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "N2" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_N2][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "O2" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_O2][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::setw(9);
    std::cout << "RCOOH" << ": ";
    std::cout << std::setw(9);
    std::cout << Species[ind_RCOOH][jNy][iNx]/airDens*1.0E+09 << " [ppb]" << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

} /* End of Solution::Debug */

/* End of Structure.cpp */
