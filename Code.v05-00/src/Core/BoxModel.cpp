/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* BoxModel Program File                                            */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 3/18/2018                                 */
/* File                 : BoxModel.cpp                              */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


static int SUCCESS     =  1;

/* STL includes */
#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>
#include <algorithm>
#include <complex>
#include <ctime>
#include <sys/stat.h>
#ifdef OMP
    #include "omp.h"
#endif /* OMP */

#include "Util/ForwardDecl.hpp"
#include "Core/Input_Mod.hpp"
#include "Core/Parameters.hpp"
#include "Core/Interface.hpp"
#include "Core/Input.hpp"
#include "Core/Structure.hpp"
#include "Core/Monitor.hpp"
#include "KPP/KPP.hpp"
#include "KPP/KPP_Parameters.h"
#include "KPP/KPP_Global.h"
#include "Core/SZA.hpp"
#include "Core/Ambient.hpp"
#include "Core/Fuel.hpp"
#include "Core/Engine.hpp"
#include "Core/Aircraft.hpp"
#include "Core/Emission.hpp"

#include "Core/Save.hpp"
static int SAVE_FAIL    = -2;

Vector_1D BuildTime( const double tStart, const double tEnd,    \
                     const double sunRise, const double sunSet, \
                     const double DYN_DT );
int BoxModel( const OptInput &Input_Opt, const Input &input )
{

    bool printDEBUG = 0;
    int isSaved_Box = 0;

#ifdef DEBUG

    std::cout << "\n DEBUG is turned ON!\n\n";
    printDEBUG = 1;

#endif /* DEBUG */
    
    /* ======================================================================= */
    /* --- Input options from the SIMULATION MENU ---------------------------- */
    /* ======================================================================= */
    
    const bool BUILD_LUT      = Input_Opt.SIMULATION_PARAMETER_SWEEP;
    /* Filename for background conditions */
    const char* BACKG_FILENAME= Input_Opt.SIMULATION_INPUT_BACKG_COND.c_str();
    
    /* ======================================================================= */
    /* ---- Input options from the CHEMISTRY MENU ---------------------------- */
    /* ======================================================================= */

    const bool HETCHEM        = Input_Opt.CHEMISTRY_HETCHEM;
    /* Folder for photolysis rates */
    const char* JRATE_FOLDER  = Input_Opt.CHEMISTRY_JRATE_FOLDER.c_str();
   
    /* Make sure that same timesteps as plume model are used */

    const double CHEMISTRY_DT = Input_Opt.CHEMISTRY_TIMESTEP;

    /* ======================================================================= */
    /* ---- Input options from the TRANSPORT MENU ---------------------------- */
    /* ======================================================================= */

    const double TRANSPORT_DT = Input_Opt.TRANSPORT_TIMESTEP;

    /* ======================================================================= */
    /* ---- Input options from the TIMESERIES MENU --------------------------- */
    /* ======================================================================= */

    const std::vector<int> TS_SPEC_LIST = Input_Opt.TS_SPECIES;

    /* Define dynamic timestep in s */
    double DYN_DT;

    /* If either TRANSPORT or CHEMISTRY is set to 0, then pick the non-zero 
     * timestep.
     * If both are non-zero, then pick the smallest timestep.
     *
     * Both TRANSPORT_DT and CHEMISTRY_DT are expressed in minutes! */

    if ( ( TRANSPORT_DT == 0.0E+00 ) || ( CHEMISTRY_DT == 0.0E+00 ) )
        DYN_DT = std::max( TRANSPORT_DT, CHEMISTRY_DT ) * 60.0;
    else if ( ( TRANSPORT_DT > 0.0E+00 ) && ( CHEMISTRY_DT > 0.0E+00 ) )
        DYN_DT = std::min( TRANSPORT_DT, CHEMISTRY_DT ) * 60.0;
    else {
        std::cout << " Invalid option when setting the dynamic timestep. Abort!" << std::endl;
        std::cout << " TRANSPORT_DT = " << TRANSPORT_DT << " min" << std::endl;
        std::cout << " CHEMISTRY_DT = " << CHEMISTRY_DT << " min" << std::endl;
        exit(1);
    }

    const double BOX_AREA = ( Input_Opt.ADV_GRID_XLIM_LEFT + Input_Opt.ADV_GRID_XLIM_RIGHT ) * ( Input_Opt.ADV_GRID_YLIM_UP + Input_Opt.ADV_GRID_YLIM_DOWN );

    /* Assign parameters */
    
    double temperature_K = input.temperature_K();
    double pressure_Pa   = input.pressure_Pa();
    double relHumidity_w = input.relHumidity_w();
    double relHumidity;
    double AerosolArea[NAERO];
    double AerosolRadi[NAERO];
    double IWC = 0;

    /* Compute relative humidity w.r.t ice */
    double relHumidity_i = relHumidity_w * physFunc::pSat_H2Ol( temperature_K )\
                                             / physFunc::pSat_H2Os( temperature_K );

    int IERR;

    /* ======================================================================= */
    /* ----------------------------------------------------------------------- */
    /* -------------- SOLAR ZENITH ANGLE ------------------------------------- */
    /* ----------------------------------------------------------------------- */
    /* ======================================================================= */

    /* Define sun parameters, this include sunrise and sunset hours and updates
     * the local solar zenith angle. */
    SZA *sun = new SZA( input.latitude_deg(), input.emissionDOY() ); 

    /* ======================================================================= */
    /* ----------------------------------------------------------------------- */
    /* ------------------------------ TIMESTEPS ------------------------------ */
    /* ----------------------------------------------------------------------- */
    /* ======================================================================= */

    /*  
     *  - tEmission is the local emission time expressed in hours 
     *  (between 0.0 and 24.0)
     *  - tInitial is the local time at which the simulation starts in hours
     *  - simulationTime represents the simulation time (in hours) (now read from
     *    input file)
     *  - tFinal corresponds to the final time of the simulation expressed in hours
     */ 

    /* Define emission and simulation time */
    const double tEmission_h = input.emissionTime();                 /* [hr] */
    const double tInitial_h  = tEmission_h;                          /* [hr] */
    const double tFinal_h    = tInitial_h + input.simulationTime();  /* [hr] */
    const double tInitial_s  = tInitial_h * 3600.0;                  /* [s] */
    const double tFinal_s    = tFinal_h   * 3600.0;                  /* [s] */

    /* Current time in [s] */
    double curr_Time_s = tInitial_s; /* [s] */
    /* Time step in [s] */
    double dt = 0;                   /* [s] */

    /* Create time array */

    /* Vector of time in [s] */
    const Vector_1D timeArray = BuildTime ( tInitial_s, tFinal_s, 3600.0*sun->sunRise, 3600.0*sun->sunSet, DYN_DT );

    /* Time counter [-] */
    UInt nTime = 0;

    /* ======================================================================= */
    /* ----------------------------------------------------------------------- */
    /* ------------------------ BACKGROUND CONDITIONS ------------------------ */
    /* ----------------------------------------------------------------------- */
    /* ======================================================================= */
    
    /* Declare solution structure */
    Solution Data(Input_Opt);

    /* Compute airDens from pressure and temperature */
    double airDens = pressure_Pa / ( physConst::kB   * temperature_K ) * 1.00E-06;
    /* [molec/cm3] = [Pa = J/m3] / ([J/K]            * [K]           ) * [m3/cm3] */
   
    /* Read ambient concentrations */
    Vector_1D amb_Value(NSPEC, 0.0);
    Vector_2D aer_Value(N_AER, Vector_1D(2, 0.0));
    std::ifstream file;

    file.open( BACKG_FILENAME );

    if ( file.is_open() ) {
        std::string line;
        UInt i = 0;

        while ( ( std::getline( file, line ) ) && ( i < NSPEC + N_AER ) ) {
            if ( ( line.length() > 0 ) && ( line != "\r" ) && ( line != "\n" ) && ( line[0] != '#' ) ) {
                std::istringstream iss(line);
                if ( i < NSPEC ) {
                    iss >> amb_Value[i];
                }
                else if ( ( i >= NSPEC ) && ( i < NSPEC + N_AER ) ) {
                    iss >> aer_Value[i - NSPEC][0];
                    std::getline( file, line );
                    std::istringstream iss(line);
                    iss >> aer_Value[i - NSPEC][1];
                }
                i++;
            }
        }
        file.close();
    }
    else {
        std::cout << "ERROR: Can't read (" << BACKG_FILENAME << ")" << std::endl;
        exit(-1);
    }

    Data.SpinUp( amb_Value, input, airDens, \
            /* Time for which ambient file is valid in hr */ (const double) 8.0 );

    /* Enforce pre-defined values? *
     * Read input defined values for background concentrations */

    /* Split NOx as NO and NO2 using the current NO/NO2 ratio */

    /* Inputs are in ppb */

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

    for ( UInt iSpec = 0; iSpec < amb_Value.size(); iSpec++ )
        amb_Value[iSpec] *= airDens;

    amb_Value[ind_H2O] = input.relHumidity_w() / ((double) 100.0) * \
                         physFunc::pSat_H2Ol( input.temperature_K() ) / ( physConst::kB * input.temperature_K() ) / 1.00E+06;

    /* Create ambient struture */
    Ambient ambientData( timeArray.size(), amb_Value, aer_Value, Vector_1D( 9, 0.0E+00 ) );

    /* ======================================================================= */
    /* ----------------------------------------------------------------------- */
    /* ----------------------------- EMISSIONS ------------------------------- */
    /* ----------------------------------------------------------------------- */
    /* ======================================================================= */

    /* Emission
     * The emissions is a combination of 
     * engine-fuel characteristics.
     * - CO2, H2O and FSC are fuel characteristics
     * - NOx, CO, HC and Soot are engine dependent.
     * An aircraft is paired with its engine.
     */

    /* Define fuel */
    char const *ChemFormula("C12H24");
    Fuel JetA( ChemFormula );

    /* Define aircraft */
    char const *aircraftName("B747-800");
    double aircraftMass = input.aircraftMass();
    std::string engineInputFile = Input_Opt.SIMULATION_INPUT_ENG_EI;
    Aircraft aircraft( aircraftName, engineInputFile, aircraftMass, \
                       temperature_K, pressure_Pa, \
                       relHumidity_w, input.nBV() );

    if ( BUILD_LUT ) {
        aircraft.setEI_NOx( input.EI_NOx() );
        aircraft.setEI_CO( input.EI_CO() );
        aircraft.setEI_HC( input.EI_HC() );
        aircraft.setEI_Soot( input.EI_Soot() );
        aircraft.setSootRad( input.sootRad() );
        aircraft.setFuelFlow( input.fuelFlow() );
        JetA.setFSC( input.EI_SO2() * (double) 500.0 );
    }

    /* Print AC Debug? */
    if ( DEBUG_AC_INPUT )
        aircraft.Debug();

    /* Aggregate emissions from engine and fuel characteristics */
    const Emission EI( aircraft.engine(), JetA );

    /* Print Emission Debug? */
    if ( DEBUG_EI_INPUT )
        EI.Debug();

    /* Add emissions */
    double E_CO2, E_H2O, E_NO, E_NO2, E_HNO2, E_SO2, E_CO, E_CH4, E_C2H6, E_PRPE, E_ALK4, E_CH2O, E_ALD2, E_GLYX, E_MGLY;
    double E_Soot;
    const double rad = EI.getSootRad();
    const double fuelPerDist = aircraft.FuelFlow() / aircraft.VFlight();
    /* Unit check:  [kg/m]   =   [kg fuel/s]    /     [m/s] */
    E_CO2  = EI.getCO2()  / ( MW_CO2  * 1.0E+03 ) * fuelPerDist * physConst::Na / BOX_AREA * 1.0E-06;
    /*     = [g/kg fuel]  / ( [kg/mol]* [g/kg]  ) * [kg fuel/m] * [molec/mol]   / [m^2]    * [m^3/cm^3]
     *     = [molec/cm^3]
     */
    E_NO   = EI.getNO()   / ( MW_NO   * 1.0E+03 ) * fuelPerDist * physConst::Na / BOX_AREA * 1.0E-06;
    E_NO2  = EI.getNO2()  / ( MW_NO2  * 1.0E+03 ) * fuelPerDist * physConst::Na / BOX_AREA * 1.0E-06;
    E_HNO2 = EI.getHNO2() / ( MW_HNO2 * 1.0E+03 ) * fuelPerDist * physConst::Na / BOX_AREA * 1.0E-06;
    E_CO   = EI.getCO()   / ( MW_CO   * 1.0E+03 ) * fuelPerDist * physConst::Na / BOX_AREA * 1.0E-06;
    E_CH4  = EI.getCH4()  / ( MW_CH4  * 1.0E+03 ) * fuelPerDist * physConst::Na / BOX_AREA * 1.0E-06;
    E_C2H6 = EI.getC2H6() / ( MW_C2H6 * 1.0E+03 ) * fuelPerDist * physConst::Na / BOX_AREA * 1.0E-06;
    E_PRPE = EI.getPRPE() / ( MW_PRPE * 1.0E+03 ) * fuelPerDist * physConst::Na / BOX_AREA * 1.0E-06;
    E_ALK4 = EI.getALK4() / ( MW_ALK4 * 1.0E+03 ) * fuelPerDist * physConst::Na / BOX_AREA * 1.0E-06;
    E_CH2O = EI.getCH2O() / ( MW_CH2O * 1.0E+03 ) * fuelPerDist * physConst::Na / BOX_AREA * 1.0E-06;
    E_ALD2 = EI.getALD2() / ( MW_ALD2 * 1.0E+03 ) * fuelPerDist * physConst::Na / BOX_AREA * 1.0E-06;
    E_GLYX = EI.getGLYX() / ( MW_GLYX * 1.0E+03 ) * fuelPerDist * physConst::Na / BOX_AREA * 1.0E-06;
    E_MGLY = EI.getMGLY() / ( MW_MGLY * 1.0E+03 ) * fuelPerDist * physConst::Na / BOX_AREA * 1.0E-06;
    E_H2O  = EI.getH2O()  / ( MW_H2O  * 1.0E+03 ) * fuelPerDist * physConst::Na / BOX_AREA * 1.0E-06;
    E_SO2  = ( 1.0 - SO2TOSO4 ) * \
             EI.getSO2()  / ( MW_SO2  * 1.0E+03 ) * fuelPerDist * physConst::Na / BOX_AREA * 1.0E-06;

    E_Soot = EI.getSoot() / ( 4.0 / 3.0 * physConst::PI * physConst::RHO_SOOT * 1.00E+03 * rad * rad * rad ) * fuelPerDist / BOX_AREA * 1.0E-06;
    /*     = [g_soot/kg_fuel]/ (                        * [kg_soot/m^3]       * [g/kg]   * [m^3]           ) * [kg_fuel/m] / [m^2] * [m^3/cm^3]
     *     = [part/cm^3]
     */

    ambientData.Species[ind_CO2][0]      += E_CO2;
    ambientData.Species[ind_NO][0]       += E_NO;
    ambientData.Species[ind_NO2][0]      += E_NO2;
    ambientData.Species[ind_HNO2][0]     += E_HNO2;
    ambientData.Species[ind_CO][0]       += E_CO;
    ambientData.Species[ind_CH4][0]      += E_CH4;
    ambientData.Species[ind_C2H6][0]     += E_C2H6;
    ambientData.Species[ind_PRPE][0]     += E_PRPE;
    ambientData.Species[ind_ALK4][0]     += E_ALK4;
    ambientData.Species[ind_CH2O][0]     += E_CH2O;
    ambientData.Species[ind_ALD2][0]     += E_ALD2;
    ambientData.Species[ind_GLYX][0]     += E_GLYX;
    ambientData.Species[ind_MGLY][0]     += E_MGLY;
    ambientData.Species[ind_H2O][0]      += E_H2O;
    ambientData.Species[ind_SO2][0]      += E_SO2;
    ambientData.sootDens[0] += E_Soot;

    
    /* ======================================================================= */
    /* ----------------------------------------------------------------------- */
    /* ----------------------------- CHEMISTRY ------------------------------- */
    /* ----------------------------------------------------------------------- */
    /* ======================================================================= */

    /* Allocate arrays for KPP */

    double STEPMIN = (double)0.0;

    double RTOL[NVAR];
    double ATOL[NVAR];

    /* Allocate photolysis rate array */
    double jRate[NPHOTOL];

    for( UInt i = 0; i < NVAR; i++ ) {
        RTOL[i] = KPP_RTOLS; 
        ATOL[i] = KPP_ATOLS; 
    }

    /* aerArray stores all the number concentrations of aerosols */
    double aerArray[N_AER][2];

    /* Ambient chemistry */
    ambientData.getData( aerArray, nTime );

    
    /* ======================================================================= */
    /* ----------------------------------------------------------------------- */
    /* ------------------------ TIME LOOP STARTS HERE ------------------------ */
    /* ----------------------------------------------------------------------- */
    /* ======================================================================= */

    while ( curr_Time_s < tFinal_s ) {
        if ( printDEBUG ) {
            /* Print message */
            std::cout << "\n";
            std::cout << "\n - Time step: " << nTime + 1 << " out of " << timeArray.size();
            #ifdef OMP
                std::cout << " (thread: " << omp_get_thread_num() << ")";
            #endif /* OMP */
            std::cout << "\n -> Solar time: " << std::fmod( curr_Time_s/3600.0, 24.0 ) << " [hr]" << std::endl;
        }

        /* ======================================================================= */
        /* ----------------------------------------------------------------------- */
        /* --------------------------- UPDATE TIMESTEP --------------------------- */
        /* ----------------------------------------------------------------------- */
        /* ======================================================================= */

        /* Compute time step */
        dt = timeArray[nTime+1] - timeArray[nTime];

        /* ======================================================================= */
        /* ----------------------------------------------------------------------- */
        /* ----------- UPDATE SOLAR ZENITH ANGLE AND PHOTOLYSIS RATES ------------ */
        /* ----------------------------------------------------------------------- */
        /* ======================================================================= */

        /* Compute the cosize of solar zenith angle midway through the integration step */
        sun->Update( curr_Time_s + dt/2 );

        /* Store cosine of solar zenith angle */
        ambientData.cosSZA[nTime] = sun->CSZA;

        if ( printDEBUG ) {
            std::cout << "\n DEBUG : \n";
            std::cout << "         CSZA = " << sun->CSZA << "\n";
        }

        /* Reset photolysis rates */
        for ( UInt iPhotol = 0; iPhotol < NPHOTOL; iPhotol++ )
            jRate[iPhotol] = 0.0E+00;

        /* If daytime, update photolysis rates */
        if ( sun->CSZA > 0.0E+00 )
            Update_JRates( jRate, sun->CSZA );

        if ( printDEBUG ) {
            std::cout << "\n DEBUG : \n";
            for ( UInt iPhotol = 0; iPhotol < NPHOTOL; iPhotol++ )
                std::cout << "         PHOTOL[" << iPhotol << "] = " << jRate[iPhotol] << "\n";
        }


        /* ======================================================================= */
        /* ----------------------------------------------------------------------- */
        /* ------------------------------- RUN KPP ------------------------------- */
        /* ------------------------ The Kinetics Pre-Processor ------------------- */
        /* ----------------------------------------------------------------------- */
        /* ======================================================================= */

        /* Ambient chemistry */
        ambientData.getData( aerArray, nTime );

        /* ~~~~~~~~~~~~~~~~~~~~~~~~ */
        /* ~~~~ Chemical rates ~~~~ */
        /* ~~~~~~~~~~~~~~~~~~~~~~~~ */

        /* Zero-out reaction rate */
        for ( UInt iReact = 0; iReact < NREACT; iReact++ )
            RCONST[iReact] = 0.0E+00;

        /* Update photolysis rates */
        for ( UInt iPhotol = 0; iPhotol < NPHOTOL; iPhotol++ )
            PHOTOL[iPhotol] = jRate[iPhotol];

        /* Update reaction rates */
        Update_RCONST( temperature_K, pressure_Pa, airDens, VAR[ind_H2O] );

        /* ~~~~~~~~~~~~~~~~~~~~~~~~ */
        /* ~~~~~ Integration ~~~~~~ */
        /* ~~~~~~~~~~~~~~~~~~~~~~~~ */

        IERR = INTEGRATE( VAR, curr_Time_s, curr_Time_s + dt, \
                          ATOL, RTOL, STEPMIN );

        if ( IERR < 0 ) {
            /* Integration failed */

                std::cout << "Integration failed";
                #ifdef OMP
                    std::cout << " on " << omp_get_thread_num();
                #endif /* OMP */
                std::cout << " for ambient conditions at time t = " << curr_Time_s/3600.0 << " ( nTime = " << nTime << " )\n";

            if ( printDEBUG ) {
                std::cout << " ~~~ Printing reaction rates:\n";
                for ( UInt iReact = 0; iReact < NREACT; iReact++ ) {
                    std::cout << "Reaction " << iReact << ": " << RCONST[iReact] << " [molec/cm^3/s]\n";
                }
                std::cout << " ~~~ Printing concentrations:\n";
                for ( UInt iSpec = 0; iSpec < NVAR; iSpec++ ) {
                    std::cout << "Species " << iSpec << ": " << VAR[iSpec]/airDens*1.0E+09 << " [ppb]\n";
                }
            }

            /* Clear dynamically allocated variable(s) */
            if ( sun != NULL )
                sun->~SZA();
            return IERR;
        }

        ambientData.FillIn( nTime + 1 );
        
        curr_Time_s += dt;
        nTime++;
    }

    /* ======================================================================= */
    /* ----------------------------------------------------------------------- */
    /* ------------------------- TIME LOOP ENDS HERE ------------------------- */
    /* ----------------------------------------------------------------------- */
    /* ======================================================================= */

    
    #pragma omp critical
    {
        isSaved_Box = output::Write_Box( input.fileName_BOX2char(), \
                                         TS_SPEC_LIST,              \
                                         ambientData,               \
                                         timeArray,                 \
                                         input,                     \
                                         airDens, relHumidity_i );
    }
    if ( isSaved_Box == output::SAVE_FAILURE ) {
        std::cout << " Saving to ring-averaged concentrations to file failed...\n";
        return SAVE_FAIL;
    }

    return SUCCESS;

} /* End of BoxModel */

/* End of BoxModel.cpp */
