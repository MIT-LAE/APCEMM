/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* PlumeModel Program File                                          */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : PlumeModel.cpp                            */
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
#include <fftw3.h>
#ifdef OMP
    #include "omp.h"
#endif /* OMP */

#include "Util/ForwardDecl.hpp"
#include "Core/Input_Mod.hpp"
#include "Core/Parameters.hpp"
#include "Core/Interface.hpp"
#include "Core/Input.hpp"
#include "Core/Monitor.hpp"
#include "SANDS/Solver.hpp"
#include "AIM/Aerosol.hpp"
#include "AIM/Coagulation.hpp"
#include "AIM/Settling.hpp"
#include "EPM/Integrate.hpp"
#include "KPP/KPP.hpp"
#include "KPP/KPP_Parameters.h"
#include "KPP/KPP_Global.h"
#include "Core/SZA.hpp"
#include "Core/Mesh.hpp"
#include "Core/Meteorology.hpp"
#include "Core/Structure.hpp"
#include "Core/Ambient.hpp"
#include "Core/Fuel.hpp"
#include "Core/Engine.hpp"
#include "Core/Aircraft.hpp"
#include "Core/Emission.hpp"
#include "Core/ReadJRates.hpp"

/* For RINGS */
#include "Core/Cluster.hpp"
#include "Core/Species.hpp"

/* For DIAGNOSTIC */
#include "Core/Diag_Mod.hpp"

#ifdef TIME_IT
    #include "Core/Timer.hpp"
#endif /* TIME_IT */

#include "Core/Save.hpp"
int isSaved = 1;
static int SAVE_FAIL   = -2;

RealDouble C[NSPEC];             /* Concentration of all species */
RealDouble * VAR = &C[0];        /* Concentration of variable species (global) */
RealDouble * FIX = &C[NVAR];     /* Concentration of fixed species (global) */

RealDouble RCONST[NREACT];       /* Rate constants (global) */
RealDouble NOON_JRATES[NPHOTOL]; /* Noon-time photolysis rates (global) */
RealDouble PHOTOL[NPHOTOL];      /* Photolysis rates (global) */
RealDouble HET[NSPEC][3];        /* Heterogeneous chemistry rates (global) */

RealDouble TIME;                 /* Current integration time (global) */

/* Require this for adjoint integration */
RealDouble SZA_CST[3];

int BoxModel( const OptInput &Input_Opt, const Input &input );
void DiffParam( const RealDouble time, RealDouble &d_x, RealDouble &d_y, \
                const RealDouble D_X, const RealDouble D_Y );
void AdvGlobal( const RealDouble time, const RealDouble T_UPDRAFT, \
                const RealDouble V_UPDRAFT,                        \
                RealDouble &v_x, RealDouble &v_y,                  \
                RealDouble &dTrav_x, RealDouble &dTrav_y );
Vector_1D BuildTime( const RealDouble tStart, const RealDouble tEnd,    \
                     const RealDouble sunRise, const RealDouble sunSet, \
                     const RealDouble DYN_DT );


int PlumeModel( OptInput &Input_Opt, const Input &input )
{

    bool printDEBUG = 1;

#ifdef DEBUG

    std::cout << "\n DEBUG is turned ON!\n\n";
    printDEBUG = 1;

#endif /* DEBUG */

    /* ======================================================================= */
    /* --- Input options from the SIMULATION MENU ---------------------------- */
    /* ======================================================================= */

    const bool RUN_BOXMODEL   = Input_Opt.SIMULATION_BOXMODEL;
    const bool BUILD_LUT      = Input_Opt.SIMULATION_PARAMETER_SWEEP;
    const bool SAVE_FORWARD   = Input_Opt.SIMULATION_SAVE_FORWARD;
    const bool ADJOINT        = Input_Opt.SIMULATION_ADJOINT;
    const char* BACKG_FILENAME= Input_Opt.SIMULATION_INPUT_BACKG_COND.c_str();
    const bool THREADED_FFT   = Input_Opt.SIMULATION_THREADED_FFT;
    const bool USE_WISDOM     = Input_Opt.SIMULATION_USE_FFTW_WISDOM;
    const char* FFTW_DIR      = Input_Opt.SIMULATION_DIRECTORY_W_WRITE_PERMISSION.c_str();

    /* ======================================================================= */
    /* ---- Input options from the TRANSPORT MENU ---------------------------- */
    /* ======================================================================= */

    const bool TRANSPORT          = Input_Opt.TRANSPORT_TRANSPORT;
    const RealDouble TRANSPORT_DT = Input_Opt.TRANSPORT_TIMESTEP;

    #ifdef RINGS
        /* The RINGS option requires that negative values are filled with
         * positive values. Otherwise, chemistry spits out garbage.
         * It's better (and safer) to always let this option turned on. */
        const bool FILLNEG        = 1;
    #else
        const bool FILLNEG        = Input_Opt.TRANSPORT_FILL;
    #endif /* RINGS */

    const bool FLUX_CORRECTION    = Input_Opt.TRANSPORT_PART_FLUX;
    const bool UPDRAFT            = Input_Opt.TRANSPORT_UPDRAFT;
    const RealDouble UPDRAFT_TIME = Input_Opt.TRANSPORT_UPDRAFT_TIMESCALE;
    const RealDouble UPDRAFT_VEL  = Input_Opt.TRANSPORT_UPDRAFT_VELOCITY;

    /* ======================================================================= */
    /* ---- Input options from the CHEMISTRY MENU ---------------------------- */
    /* ======================================================================= */

    const bool CHEMISTRY          = Input_Opt.CHEMISTRY_CHEMISTRY;
    const RealDouble CHEMISTRY_DT = Input_Opt.CHEMISTRY_TIMESTEP;
    const bool HETCHEM            = Input_Opt.CHEMISTRY_HETCHEM;
    const char* JRATE_FOLDER      = Input_Opt.CHEMISTRY_JRATE_FOLDER.c_str();

    /* ======================================================================= */
    /* ---- Input options from the AEROSOL MENU ------------------------------ */
    /* ======================================================================= */

    const bool GRAVSETTLING  = Input_Opt.AEROSOL_GRAVSETTLING;
    const bool ICE_COAG      = Input_Opt.AEROSOL_COAGULATION_SOLID;
    const bool LIQ_COAG      = Input_Opt.AEROSOL_COAGULATION_LIQUID;
    const RealDouble COAG_DT = Input_Opt.AEROSOL_COAGULATION_TIMESTEP;
    const bool ICE_GROWTH    = Input_Opt.AEROSOL_ICE_GROWTH;

    /* ======================================================================= */
    /* ---- Input options from the METEOROLOGY MENU -------------------------- */
    /* ======================================================================= */

    /* Input options from the METEOROLOGY MENU are read in the Meteorology
     * subroutine */

    /* ======================================================================= */
    /* ---- Input options from the DIAGNOSTIC MENU --------------------------- */
    /* ======================================================================= */

    const char* DIAG_FILENAME = Input_Opt.DIAG_FILENAME.c_str();

    /* ======================================================================= */
    /* ---- Input options from the TIMESERIES MENU --------------------------- */
    /* ======================================================================= */

    std::string TS_FOLDER = "";

    /* If timeseries output is desired for multiple cases at a time, uncomment
     * the following lines */
//    TS_FOLDER += "Case" + std::to_string(input.Case());
//    TS_FOLDER += "/";

    std::string TS_FILE1, TS_FILE2;
    const bool TS_SPEC                  = Input_Opt.TS_SPEC;
    TS_FILE1                            = TS_FOLDER + Input_Opt.TS_FILENAME;
    const char* TS_SPEC_FILENAME        = TS_FILE1.c_str();
    const std::vector<int> TS_SPEC_LIST = Input_Opt.TS_SPECIES;
    const RealDouble TS_FREQ            = Input_Opt.TS_FREQ;

    const bool TS_AERO                  = Input_Opt.TS_AERO;
    TS_FILE2                            = TS_FOLDER + Input_Opt.TS_AERO_FILENAME;
    const char* TS_AERO_FILENAME        = TS_FILE2.c_str();
    const std::vector<int> TS_AERO_LIST = Input_Opt.TS_AEROSOL;
    const RealDouble TS_AERO_FREQ       = Input_Opt.TS_AERO_FREQ;

    if ( TS_SPEC )
        std::cout << "\n Saving TS files to: " << TS_SPEC_FILENAME << std::endl;

    if ( TS_AERO )
        std::cout << "\n Saving TS_AERO files to: " << TS_AERO_FILENAME << std::endl;

    if ( ( TS_SPEC || TS_AERO ) && ( TS_FOLDER.compare("") != 0 ) ) {

        /* Create output directory for timeseries */
        struct stat sb;
        if ( !( stat( TS_FOLDER.c_str(), &sb ) == 0 \
                    && S_ISDIR(sb.st_mode) ) ) {

            /* Create directory */
            const int dir_err = \
                    mkdir( TS_FOLDER.c_str(), \
                            S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );

            if ( dir_err == -1 ) {
                std::cout << " Could not create directory: ";
                std::cout << TS_FOLDER << std::endl;
                std::cout << " You may not have write permission" << std::endl;
                exit(1);
            }
        }
    }

    /* ======================================================================= */
    /* ---- Input options from the PROD & LOSS MENU -------------------------- */
    /* ======================================================================= */

    /* TODO: Implement PL rates */
    const bool SAVE_PL   = Input_Opt.PL_PL;
    const bool SAVE_O3PL = Input_Opt.PL_O3;

    /* Define dynamic timestep in s */
    RealDouble DYN_DT;

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

    /* If DIAG_OUTPUT is turned on, make sure that output timestep is a
     * multiple of the dynamic timestep */
    if ( ( TS_SPEC ) && ( TS_FREQ > 0.0E+00 ) && \
         ( std::fmod( TS_FREQ / DYN_DT * 60, 1 ) != 0 ) ) {
        std::cout << " Timeseries frequency should be a multiple of the dynamic timestep!" << std::endl;
        std::cout << " Output might be compromised!" << std::endl;
    }

    /* Assign parameters */

    RealDouble temperature_K = input.temperature_K();
    RealDouble pressure_Pa   = input.pressure_Pa();
    RealDouble relHumidity_w = input.relHumidity_w();

    /* Compute relative humidity w.r.t ice */
    RealDouble relHumidity_i = relHumidity_w * physFunc::pSat_H2Ol( temperature_K )\
                                             / physFunc::pSat_H2Os( temperature_K );

    /* Grid indices */
    UInt iNx = 0;
    UInt jNy = 0;

    /* Species index */
    UInt N   = 0;

    int IERR;

#ifdef TIME_IT

    Timer Stopwatch, Stopwatch_cumul;
    unsigned long SANDS_clock, KPP_clock;
    unsigned long SANDS_clock_cumul = 0;
    unsigned long KPP_clock_cumul = 0;
    unsigned long clock_cumul = 0;
    bool reset = 1;

#endif /* TIME_IT */

#if ( NOy_MASS_CHECK )

    RealDouble mass_Ambient_NOy, mass_Emitted_NOy;

    #ifdef RINGS
        RealDouble mass_Emitted_NOy_Rings;
    #endif /* RINGS */

#endif /* NOy_MASS_CHECK */

#if ( CO2_MASS_CHECK )

    RealDouble mass_Ambient_CO2, mass_Emitted_CO2;

    #ifdef RINGS
        RealDouble mass_Emitted_CO2_Rings;
    #endif /* RINGS */

#endif /* CO2_MASS_CHECK */

#if ( H2O_MASS_CHECK )

    /* Conversion factor from ice volume [m^3] to [molecules] */
    const RealDouble UNITCONVERSION = physConst::RHO_ICE / MW_H2O * physConst::Na;
    /* Unit check: [kg/m^3] / [kg/mol] * [molec/mol] = [molec/m^3] */

    RealDouble mass_Ambient_H2O, mass_H2O;

    Vector_2D totIceVol;

#endif /* H2O_MASS_CHECK */

    /* ======================================================================= */
    /* ----------------------------------------------------------------------- */
    /* --------------------------------- MESH -------------------------------- */
    /* ----------------------------------------------------------------------- */
    /* ======================================================================= */

    Mesh m;
    const Vector_1D xE = m.xE();
    const Vector_1D yE = m.yE();

    /* Get cell areas */
    const Vector_2D cellAreas = m.areas();

    /* ======================================================================= */
    /* ----------------------------------------------------------------------- */
    /* -------------- SOLAR ZENITH ANGLE + PHOTOLYSIS RATES ------------------ */
    /* ----------------------------------------------------------------------- */
    /* ======================================================================= */

    /* Allocate photolysis rate array */
    RealDouble jRate[NPHOTOL];

    /* Define sun parameters, this include sunrise and sunset hours and updates
     * the local solar zenith angle. */
    SZA *sun = new SZA( input.latitude_deg(), input.emissionDOY() );

    /* Initialize noon time photolysis rates
     * The data is after the quantum yield has been applied and represents
     * the value of the photolysis rates at 12:00 (noon) locally.
     * The photolysis rates at any given time are obtained by multiplying
     * those by the cosine of the solar zenith angle, when positive. */
    for ( UInt iPhotol = 0; iPhotol < NPHOTOL; iPhotol++ )
        NOON_JRATES[iPhotol] = 0.0E+00;

    /* Allocating noon-time photolysis rates. */

    if ( CHEMISTRY ) {
        #pragma omp critical
        {
            ReadJRates( JRATE_FOLDER,  \
                input.emissionMonth(), \
                input.emissionDay(),   \
                input.longitude_deg(), \
                input.latitude_deg(),  \
                pressure_Pa/100.0,     \
                NOON_JRATES );
        }

        if ( printDEBUG ) {
            std::cout << "\n DEBUG : \n";
            for ( UInt iPhotol = 0; iPhotol < NPHOTOL; iPhotol++ ) {
                std::cout << "         NOON_JRATES[" << iPhotol << "] = ";
                std::cout << NOON_JRATES[iPhotol] << "\n";
            }
        }
    }

    /* Run box model */
    if ( RUN_BOXMODEL )
        BoxModel( Input_Opt, input );


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
    const RealDouble tEmission_h = input.emissionTime();                 /* [hr] */
    const RealDouble tInitial_h  = tEmission_h;                          /* [hr] */
    const RealDouble tFinal_h    = tInitial_h + input.simulationTime();  /* [hr] */
    const RealDouble tInitial_s  = tInitial_h * 3600.0;                  /* [s] */
    const RealDouble tFinal_s    = tFinal_h   * 3600.0;                  /* [s] */

    /* Current time in [s] */
    RealDouble curr_Time_s = tInitial_s; /* [s] */
    /* Time step in [s] */
    RealDouble dt = 0;                   /* [s] */

    /* Create time array */

    /* Vector of time in [s] */
    const Vector_1D timeArray = BuildTime ( tInitial_s, tFinal_s, 3600.0*sun->sunRise, 3600.0*sun->sunSet, DYN_DT );

    /* Time counter [-] */
    UInt nTime = 0;

    bool LAST_STEP                    = 0;
    bool ITS_TIME_FOR_LIQ_COAGULATION = 0;
    bool ITS_TIME_FOR_ICE_COAGULATION = 0;
    bool ITS_TIME_FOR_ICE_GROWTH      = 0;
    RealDouble lastTimeLiqCoag        = curr_Time_s;
    RealDouble lastTimeIceCoag        = curr_Time_s;
    RealDouble lastTimeIceGrowth      = curr_Time_s;
    RealDouble dtLiqCoag              = 0.0E+00;
    RealDouble dtIceCoag              = 0.0E+00;
    RealDouble dtIceGrowth            = 0.0E+00;

    /* ======================================================================= */
    /* ----------------------------------------------------------------------- */
    /* ----------------------------- METEOROLOGY ----------------------------- */
    /* ----------------------------------------------------------------------- */
    /* ======================================================================= */

    Meteorology Met( Input_Opt, curr_Time_s / 3600.0, m,        \
                     temperature_K, pressure_Pa, relHumidity_i, \
                     printDEBUG );
    if ( Input_Opt.MET_LOADMET && Input_Opt.MET_LOADTEMP ) {
        temperature_K = Met.temp_user;
        relHumidity_i = relHumidity_w * physFunc::pSat_H2Ol( temperature_K )\
                                      / physFunc::pSat_H2Os( temperature_K );
    }
    if ( Input_Opt.MET_LOADMET && Input_Opt.MET_LOADH2O ) {
        relHumidity_w = Met.RHw_user;
        Input_Opt.MET_DEPTH = Met.satdepth_user;
        relHumidity_i = relHumidity_w * physFunc::pSat_H2Ol( temperature_K )\
                                      / physFunc::pSat_H2Os( temperature_K );
    }
    std::cout << "Temperature      = " << temperature_K << " K" << std::endl;
    std::cout << "Rel. humidity    = " << relHumidity_w << " %" << std::endl;
    std::cout << "Saturation depth = " << Input_Opt.MET_DEPTH << " m" << std::endl;

    /* ======================================================================= */
    /* ----------------------------------------------------------------------- */
    /* ------------------------ BACKGROUND CONDITIONS ------------------------ */
    /* ----------------------------------------------------------------------- */
    /* ======================================================================= */

    /* Declare solution structure */
    Solution Data;

    /* Compute airDens from pressure and temperature */
    RealDouble airDens = pressure_Pa / ( physConst::kB   * temperature_K ) * 1.00E-06;
    /*     [molec/cm3] = [Pa = J/m3] / ([J/K]            * [K]           ) * [m3/cm3] */

    /* Set solution arrays to ambient data */
    Data.Initialize( BACKG_FILENAME,      \
                     input, airDens, Met, \
                     Input_Opt,           \
                     printDEBUG );

    /* Print Background Debug? */
    if ( DEBUG_BG_INPUT )
        Data.Debug( airDens );

    /* Create ambient struture */
    Ambient ambientData( timeArray.size(), Data.getAmbient(), Data.getAerosol(), Data.getLiqSpecies() );

    /* ======================================================================= */
    /* ----------------------------------------------------------------------- */
    /* --------------------------- TRANSPORT SOLVER -------------------------- */
    /* ----------------------------------------------------------------------- */
    /* ======================================================================= */

    /* Allocate horizontal and vertical diffusion parameters */
    RealDouble d_x, d_y;

    /* Allocate horizontal and vertical advection parameters */
    /* These correspond to domain-wide advection velocities (updraft, downdraft) */
    RealDouble vGlob_x, vGlob_y;

    /* Allocate horizontal and vertical distance traveled */
    RealDouble dTrav_x, dTrav_y;

    /* Allocate steady-state diffusion parameters */
    RealDouble D_X, D_Y;

    /* Allocate shear */
    RealDouble shear;
    int LASTINDEX_SHEAR;

    /* Initialize */
    d_x = 0;       /* [m2/s] */
    d_y = 0;       /* [m2/s] */
    vGlob_x = 0;   /* [m/s]  */
    vGlob_y = 0;   /* [m/s]  */
    dTrav_x = 0;   /* [m]    */
    dTrav_y = 0;   /* [m]    */

    D_X   = input.horizDiff(); /* [m^2/s] */
    D_Y   = input.vertiDiff(); /* [m^2/s] */
    shear = input.shear();     /* [1/s] */
    if ( shear >= 0.0E+00 )
        LASTINDEX_SHEAR = NX-1;
    else
        LASTINDEX_SHEAR = 0;

    /* Fill with? */
    const RealDouble fillWith = 0.0E+00;

    /* Allocate Solvers */
    SANDS::Solver Solver;
    #pragma omp critical
    {
        std::cout << "\n Initializing solver..." << std::endl;
        Solver.Initialize( /* Use threaded FFT?    */ THREADED_FFT, \
                           /* Use FFTW wisdom?     */ USE_WISDOM,   \
                           /* FFTW Directory       */ FFTW_DIR,     \
                           /* Fill negative values */ FILLNEG,      \
                           /* Fill with this value */ fillWith );
        std::cout << "\n Initialization complete..." << std::endl;
    }


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
    RealDouble aircraftMass = input.aircraftMass();
    Aircraft aircraft( aircraftName, aircraftMass, \
                       temperature_K, pressure_Pa, relHumidity_w );

    /* Overwrite engine conditions with input parameters */
    aircraft.setEI_NOx( input.EI_NOx() );
    aircraft.setEI_CO( input.EI_CO() );
    aircraft.setEI_HC( input.EI_HC() );
    aircraft.setEI_Soot( input.EI_Soot() );
    aircraft.setSootRad( input.sootRad() );
    aircraft.setVFlight( input.flightSpeed(), temperature_K );
    aircraft.setEngNumber( input.numEngines() );
    aircraft.setWingspan( input.wingspan() );
    aircraft.setFuelFlow( input.fuelFlow() );
    JetA.setFSC( input.EI_SO2() * (RealDouble) 500.0 );

    /* Print AC Debug? */
    if ( DEBUG_AC_INPUT )
        aircraft.Debug();

    /* Aggregate emissions from engine and fuel characteristics */
    const Emission EI( aircraft.engine(), JetA );

    /* Print Emission Debug? */
    if ( DEBUG_EI_INPUT )
        EI.Debug();

    /* ======================================================================= */
    /* ----------------------------------------------------------------------- */
    /* ----------------------------- CHEMISTRY ------------------------------- */
    /* ----------------------------------------------------------------------- */
    /* ======================================================================= */

    /* Allocate arrays for KPP */

    RealDouble STEPMIN = (RealDouble)0.0;

    /* Allocate tolerances */
    RealDouble RTOL[NVAR];
    RealDouble ATOL[NVAR];

    /* Allocate RealDoubles to store RH and IWC */
    RealDouble relHumidity, IWC;

    /* Initialize tolerances */
    for( UInt i = 0; i < NVAR; i++ ) {
        RTOL[i] = KPP_RTOLS;
        ATOL[i] = KPP_ATOLS;
    }


    /* aerArray stores all the number concentrations of aerosols */
    RealDouble aerArray[N_AER][2];

    /* Ambient chemistry */
    ambientData.getData( aerArray, nTime );

    /* ======================================================================= */
    /* ----------------------------------------------------------------------- */
    /* -------------------------- EARLY MICROPHYSICS ------------------------- */
    /* ----------------------------------------------------------------------- */
    /* ======================================================================= */
    std::cout << "\n Starting EPM..." << std::endl;

    RealDouble Ice_rad, Ice_den, Soot_den, H2O_mol, SO4g_mol, SO4l_mol;
    RealDouble areaPlume;
    RealDouble Ab0 = input.bypassArea();
    RealDouble Tc0 = input.coreExitTemp();
    AIM::Aerosol liquidAer, iceAer;
    EPM::Integrate( temperature_K, pressure_Pa, relHumidity_w, VAR, FIX, \
                    aerArray, aircraft, EI, Ice_rad, Ice_den, Soot_den,  \
                    H2O_mol, SO4g_mol, SO4l_mol, liquidAer, iceAer, areaPlume, \
                    Ab0, Tc0 );

    /* Compute initial plume area.
     * If 2 engines, we assume that after 3 mins, the two plumes haven't fully mixed yet and result in a total
     * area of 2 * the area computed for one engine
     * If 3 or more engines, we assume that the plumes originating from the same wing have mixed. */

    areaPlume *= 2;
    if ( aircraft.EngNumber() != 2 ) {
        Ice_den  *= aircraft.EngNumber() / 2.0;
        liquidAer.scalePdf( aircraft.EngNumber() / 2.0 );
        iceAer.scalePdf( aircraft.EngNumber() / 2.0 );
        Soot_den *= aircraft.EngNumber() / 2.0;
    }

    if ( iceAer.Moment() != 0 ) {
        /* Apply ice particle vortex losses using parameterization from
         * large-eddy simulations */
        /* TODO: Change Input_Opt.MET_DEPTH to actual depth from meteorology and not just
         * user-specified input */
        const RealDouble iceNumFrac = aircraft.VortexLosses( EI.getSoot(),    \
                                                             EI.getSootRad(), \
                                                             Input_Opt.MET_DEPTH );
        if ( iceNumFrac <= 0.00E+00 && !CHEMISTRY ) {
            std::cout << "EndSim: vortex sinking" << std::endl;
            exit(0);
        }
        iceAer.scalePdf( iceNumFrac );
    }

    const RealDouble semiYaxis = 0.5 * aircraft.deltaz1();
    const RealDouble semiXaxis = areaPlume / ( physConst::PI * semiYaxis );


    /* Liquid aerosol considerations */
    UInt LA_MICROPHYSICS;

    /* Do we have emitted sulfate aerosols? */
    if ( liquidAer.Moment() != 0 )
        /* If yes, then microphysics in all grid cells */
        LA_MICROPHYSICS = 2;
    else {
        /* If no, then do we have background liquid aerosols? */
        if ( Data.LA_nDens != 0 )
            /* If yes, then all grid cells' microphysics will be the same */
            LA_MICROPHYSICS = 1;
        else
            /* If no, then we have no liquid particles */
            LA_MICROPHYSICS = 0;
    }

    /* Transport for liquid aerosols? */
    const bool TRANSPORT_LA = ( LA_MICROPHYSICS == 2 );

    /* Solid aerosol considerations */
    UInt PA_MICROPHYSICS;

    /* Do we have a contrail? */
    if ( Ice_den != 0 )
        /* If yes, then microphysics in all grid cells */
        PA_MICROPHYSICS = 2;
    else {
        /* If no, then do we have background solid aerosols? */
        if ( Data.PA_nDens != 0 )
            /* If yes, then all grid cells' microphysics will be the same */
            PA_MICROPHYSICS = 1;
        else
            /* If no, then we have no solid particles */
            PA_MICROPHYSICS = 0;
    }

    /* Transport for solid aerosols? */
    const bool TRANSPORT_PA = ( PA_MICROPHYSICS == 2 );

    Vector_1D vFall( Data.nBin_PA, 0.0E+00 );
    if ( TRANSPORT_PA && GRAVSETTLING ) {
        /* Compute settling velocities */
        vFall = AIM::SettlingVelocity( Data.solidAerosol.getBinCenters(), \
                                       temperature_K, pressure_Pa );
    }

#ifdef RINGS

    /* ======================================================================= */
    /* ----------------------------------------------------------------------- */
    /* ---------------------------- RING STRUCTURE --------------------------- */
    /* ----------------------------------------------------------------------- */
    /* ======================================================================= */

    /* Ring index */
    UInt iRing = 0;

    /* Create cluster of rings */
    //Cluster ringCluster( NRING, ( relHumidity_i > 100.0 ) );
    Cluster ringCluster( NRING, 0 );

    /* Number of rings */
    const UInt nRing = ringCluster.getnRing();

    /* Print Ring Debug? */
    if ( DEBUG_RINGS )
        ringCluster.Debug();

    /* Allocate species-ring vector */
    SpeciesArray ringData( nRing, timeArray.size(), ringCluster.halfRing() );

    /* Compute Grid to Ring mapping */
    m.Ring2Mesh( ringCluster );
    Vector_2Dui mapIndices = m.mapIndex();

    /* Print ring to mesh mapping? */
    if ( DEBUG_MAPPING )
        m.Debug();

    /* Compute ring areas
     * Note: The rings are only affected by shear and NOT diffusion.
     * When shear is applied to a N-D potato, it does NOT modify its area. Think of
     * shear as advection in infinitesimal layers, each having a different velocity.
     * We can thus compute the ring areas initially, once and for all. */
    ringCluster.ComputeRingAreas( cellAreas, m.weights );
    const Vector_1D ringArea = ringCluster.getRingArea();
    const RealDouble totArea = std::accumulate( ringArea.begin(), ringArea.end(), 0 );

    /* Add emission into the grid */
    Data.addEmission( EI, aircraft, m, ringCluster.halfRing(),  \
                      temperature_K, ( relHumidity_i > 100.0 ), \
                      liquidAer, iceAer, Soot_den, Met, areaPlume );
    /* Fill in variables species for initial time */
    ringData.FillIn( Data, m.weights, nTime );

    /* Allocate an additional array for KPP */
    RealDouble tempArray[NVAR];

    /* Allocate a ring's relative humidity */
    RealDouble AerosolArea[NAERO];
    RealDouble AerosolRadi[NAERO];

    /* Otherwise we do not have a ring structure and chemistry is solved on
     * the grid */

#else

    /* Initialization at the grid scale level */

    /* Use the most inner ring to initialize at the grid-scale level */
    Cluster ringCluster( NRING, 0 );

    /* Compute Grid to Ring mapping */
    m.Ring2Mesh( ringCluster );

    /* Compute ring areas
     * Note: The rings are only affected by shear and NOT diffusion.
     * When shear is applied to a N-D potato, it does NOT modify its area. Think of
     * shear as advection in infinitesimal layers, each having a different velocity.
     * We can thus compute the ring areas initially, once and for all. */
    ringCluster.ComputeRingAreas( cellAreas, m.weights );
    const Vector_1D ringArea = ringCluster.getRingArea();

    /* Add emission into the grid */
    Data.addEmission( EI, aircraft, m, ringCluster.halfRing(),  \
                      temperature_K, ( relHumidity_i > 100.0 ), \
                      liquidAer, iceAer, Soot_den, Met, areaPlume );

#endif /* RINGS */

    /* ======================================================================= */
    /* ----------------------------------------------------------------------- */
    /* ---------------- OUTPUT RUN CHARACTERISTICS TO LOG FILE --------------- */
    /* ----------------------------------------------------------------------- */
    /* ======================================================================= */

    #pragma omp critical
    {
        #ifdef OMP
            std::cout << "\n\n ## ON THREAD: " << omp_get_thread_num() << "\n ##";
        #endif /* OMP */
        const UInt coutPrecision = 5;
        const UInt txtWidth      = coutPrecision + 2;
        std::cout << std::setprecision(coutPrecision);
        std::cout << "\n ## ATMOSPHERIC COND.:";
        std::cout << "\n ##\n";
        std::cout << " ## - Temperature: " << std::setw(txtWidth) << temperature_K          << " [    K]\n";
        std::cout << " ## - Pressure   : " << std::setw(txtWidth) << pressure_Pa * 1.00E-02 << " [  hPa]\n";
        std::cout << " ## - Rel. Hum. I: " << std::setw(txtWidth) << relHumidity_i          << " [    %]\n";
        std::cout << " ## - Horiz. Diff: " << std::setw(txtWidth) << input.horizDiff()      << " [ m2/s]\n";
        std::cout << " ## - Verti. Diff: " << std::setw(txtWidth) << input.vertiDiff()      << " [ m2/s]\n";
        std::cout << " ## - Shear      : " << std::setw(txtWidth) << input.shear()          << " [  1/s]\n";
        std::cout << " ## - Latitude   : " << std::setw(txtWidth) << input.latitude_deg()   << " [  deg]\n";
        std::cout << " ## - Longitude  : " << std::setw(txtWidth) << input.longitude_deg()  << " [  deg]\n";
        std::cout << " ## - Max CSZA   : " << std::setw(txtWidth) << sun->CSZA_max          << " [  -  ]\n";
        std::cout << " ## - Emiss. time: " << std::setw(txtWidth) << input.emissionTime()   << " [  hrs]\n";
        std::cout << " ## - Emiss. day : " << std::setw(txtWidth-3) << input.emissionMonth() << "/" << input.emissionDay() << "\n";

        std::cout << "\n ## EMISSIONS:";
        std::cout << "\n ##\n";
        std::cout << " ## - E_CO2 = " << std::setw(txtWidth+3) << EI.getCO2() * aircraft.FuelFlow() / aircraft.VFlight()           << " [kg(CO2)/km]"\
            " ( EI  = " << std::setw(txtWidth) << EI.getCO2() * 1.00E-03 << " [kg/kg_fuel]      )\n";
        std::cout << " ## - E_CO  = " << std::setw(txtWidth+3) << EI.getCO()  * aircraft.FuelFlow() / aircraft.VFlight() * 1.0E+03 << " [ g(CO) /km]"\
            " ( EI  = " << std::setw(txtWidth) << EI.getCO()             << " [ g/kg_fuel]      )\n";
        std::cout << " ## - E_CH4 = " << std::setw(txtWidth+3) << EI.getCH4() * aircraft.FuelFlow() / aircraft.VFlight() * 1.0E+06 << " [mg(CH4)/km]"\
            " ( EI  = " << std::setw(txtWidth) << EI.getCH4() * 1.00E+03 << " [mg/kg_fuel]      )\n";
        std::cout << " ## - E_NOx = " << std::setw(txtWidth+3) << ( EI.getNO() / MW_NO + EI.getNO2() / MW_NO2 + EI.getHNO2() / MW_HNO2 ) * MW_NO2 \
                                                  * aircraft.FuelFlow() / aircraft.VFlight() * 1.0E+03 << " [ g(NO2)/km]"\
            " ( EI  = " << std::setw(txtWidth) << EI.getNOx()            << " [ g(NO2)/kg_fuel] )\n";
        std::cout << " ## - E_NOx = " << std::setw(txtWidth+3) << ( EI.getNO() / MW_NO + EI.getNO2() / MW_NO2 + EI.getHNO2() / MW_HNO2 ) * MW_N \
                                                  * aircraft.FuelFlow() / aircraft.VFlight() * 1.0E+03 << " [ g(N)  /km]\n";
        std::cout << " ## - E_SO2 = " << std::setw(txtWidth+3) << EI.getSO2() * aircraft.FuelFlow() / aircraft.VFlight() * 1.0E+03 << " [ g(SO2)/km]"\
            " ( EI  = " << std::setw(txtWidth) << EI.getSO2()            << " [ g/kg_fuel]      )\n";
        std::cout << " ##                                   ( FSC = " << std::setw(txtWidth) << JetA.getFSC() << " [-]               )\n";
        std::cout << " ## - E_Soo = " << std::setw(txtWidth+3) << EI.getSoot() * aircraft.FuelFlow() / aircraft.VFlight() * 1.0E+03 << " [ g(Soo)/km]"\
            " ( EI  = " << std::setw(txtWidth) << EI.getSoot()* 1.00E+03 << " [mg/kg_fuel]      )\n";
        std::cout << " ## - E_Soo = " << std::setw(txtWidth+3) << EI.getSoot() * aircraft.FuelFlow() / aircraft.VFlight() * 1.0E+03 / ( 4.0 / 3.0 * physConst::PI * physConst::RHO_SOOT * 1.00E+03 * EI.getSootRad() * EI.getSootRad() * EI.getSootRad() ) << " [ #(Soo)/km]"\
            " ( GMD = " << std::setw(txtWidth) << 2.0 * EI.getSootRad() * 1.0E+09 << " [nm]              )\n";
        std::cout << " ## - Fflow = " << std::setw(txtWidth+3) << aircraft.FuelFlow() << " [      kg/s]\n";
        std::cout << " ## - AMass = " << std::setw(txtWidth+3) << input.aircraftMass() << " [      kg  ]\n";

        std::cout << "\n ## AEROSOLS:";
        std::cout << "\n ##\n";
        std::cout << " ## - LA : " << std::setw(txtWidth+3) << liquidAer.Moment() << " [#/cm^3], \n";
        std::cout << " ##        " << std::setw(txtWidth+3) << liquidAer.EffRadius() * 1.0E+09 << " [nm], \n";
        std::cout << " ##        " << std::setw(txtWidth+3) << liquidAer.Moment(2) * 1.0E+12 << " [mum^2/cm^3] \n";
        std::cout << " ##\n";
        std::cout << " ## - PA : " << std::setw(txtWidth+3) << iceAer.Moment() << " [#/cm^3], \n";
        if ( iceAer.Moment(2) > 0 )
            std::cout << " ##        " << std::setw(txtWidth+3) << iceAer.EffRadius() * 1.0E+06 << " [mum], \n";
        else
            std::cout << " ##        " << std::setw(txtWidth+3) << 0.0E+00 << " [mum], \n";
        std::cout << " ##        " << std::setw(txtWidth+3) << iceAer.Moment(2) * 1.0E+12 << " [mum^2/cm^3] \n";

        std::cout << "\n ## BACKG COND.:";
        std::cout << "\n ##\n";
        std::cout << " ## - NOx  = " << std::setw(txtWidth) << ( VAR[ind_NO] + VAR[ind_NO2] ) / airDens * 1.0E+12 << " [ppt]\n";
        std::cout << " ## - HNO3 = " << std::setw(txtWidth) << ( VAR[ind_HNO3] ) / airDens * 1.0E+12 << " [ppt]\n";
        std::cout << " ## - O3   = " << std::setw(txtWidth) << ( VAR[ind_O3] )   / airDens * 1.0E+09 << " [ppb]\n";
        std::cout << " ## - CO   = " << std::setw(txtWidth) << ( VAR[ind_CO] )   / airDens * 1.0E+09 << " [ppb]\n";
        std::cout << " ##\n";
        std::cout << " ## - LA : " << std::setw(txtWidth+3) << Data.LA_nDens << " [#/cm^3], \n";
        std::cout << " ##        " << std::setw(txtWidth+3) << Data.LA_rEff  << " [nm], \n";
        std::cout << " ##        " << std::setw(txtWidth+3) << Data.LA_SAD   << " [mum^2/cm^3] \n";
        std::cout << " ##\n";
        std::cout << " ## - PA : " << std::setw(txtWidth+3) << Data.PA_nDens << " [#/cm^3], \n";
        std::cout << " ##        " << std::setw(txtWidth+3) << Data.PA_rEff * 1.00E-03  << " [mum], \n";
        std::cout << " ##        " << std::setw(txtWidth+3) << Data.PA_SAD   << " [mum^2/cm^3] \n";

        /* Force flush */
        std::cout << std::endl;

    }

    /* Timeseries diagnostics */
    if ( TS_SPEC ) {
        std::cout << "Saving chemistry" << std::endl;
        int hh = (int) (curr_Time_s - timeArray[0])/3600;
        int mm = (int) (curr_Time_s - timeArray[0])/60   - 60 * hh;
        int ss = (int) (curr_Time_s - timeArray[0])      - 60 * ( mm + 60 * hh );
        Diag_TS_Chem( TS_SPEC_FILENAME, TS_SPEC_LIST, hh, mm, ss, \
                      Data, m );
    }

    if ( TS_AERO ) {
        std::cout << "Saving aerosol" << std::endl;
        int hh = (int) (curr_Time_s - timeArray[0])/3600;
        int mm = (int) (curr_Time_s - timeArray[0])/60   - 60 * hh;
        int ss = (int) (curr_Time_s - timeArray[0])      - 60 * ( mm + 60 * hh );
        Diag_TS_Phys( TS_AERO_FILENAME, TS_AERO_LIST, hh, mm, ss, \
                      Data, m, Met );
        float totalIceParticles = Data.solidAerosol.TotalNumber_sum( cellAreas );
        float totalIceMass = Data.solidAerosol.TotalIceMass_sum( cellAreas );
        if ( totalIceParticles <= 1.00E+1 && totalIceMass <= 1.00E-5 && !CHEMISTRY ) {
            std::cout << "EndSim: no particles remain" << std::endl;
            exit(0);
        }
    }

    /* Prod & loss diagnostics */

#ifdef RINGS

    /* Rates before chemistry is performed.
     * Chemistry is performed NT-1 times, the size is thus:
     * (NT-1) x NRING x NFAM
     * and
     * (NT-1) x NFAM for ambient conditions */

    UInt NFAM_ = NFAM;
    if ( !SAVE_PL && SAVE_O3PL )
        NFAM_ = 2;
    else if ( !SAVE_PL && !SAVE_O3PL )
        NFAM_ = 0;

    Vector_3D plumeRates( timeArray.size() - 1, Vector_2D( NRING, Vector_1D( NFAM_, 0.0E+00 ) ) );
    Vector_2D ambientRates( timeArray.size() - 1, Vector_1D( NFAM_, 0.0E+00 ) );

#else

    // TODO!!
    if ( SAVE_PL ) {

        /* If chemistry is performed at the grid cell level, then the
         * rates are stored as:
         * NY x NX x NFAM in [molec/cm^3/s]
         * into netCDF files at a frequency specified by the input file */

        int hh = (int) (curr_Time_s - timeArray[0])/3600;
        int mm = (int) (curr_Time_s - timeArray[0])/60   - 60 * hh;
        int ss = (int) (curr_Time_s - timeArray[0])      - 60 * ( mm + 60 * hh );
//        Diag_PL( "PL_hhmmss.nc", hh, mm, ss, Data, m );

    } else {
        if ( SAVE_O3PL ) {
        }
    }

#endif /* RINGS */

    /* ======================================================================= */
    /* ----------------------------------------------------------------------- */
    /* ------------------------ TIME LOOP STARTS HERE ------------------------ */
    /* ----------------------------------------------------------------------- */
    /* ======================================================================= */

#ifdef TIME_IT

    Stopwatch_cumul.Start( );

#endif /* TIME_IT */
    
    //std::cout << curr_Time_s < tFinal_s << std::endl;
    while ( curr_Time_s < tFinal_s ) {
        
        if ( printDEBUG || 1 ) {
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
        LAST_STEP = ( curr_Time_s + dt >= tFinal_s );
        Solver.UpdateTimeStep( dt );
 
        /* ======================================================================= */
        /* ----------------------------------------------------------------------- */
        /* ---------------------- UPDATE TRANSPORT PARAMETERS -------------------- */
        /* ----------------------------------------------------------------------- */
        /* ======================================================================= */

        /* Compute diffusion parameters */
        /* Is transport turned on? */
        if ( TRANSPORT ) {

            /* Compute diffusion parameters at mid time step */

            /* d_x: horizontal diffusion coefficient [m^2/s]
             * d_y: vertical diffusion coefficient [m^2/s]
             */

            DiffParam( curr_Time_s - tInitial_s + dt/2.0, d_x, d_y, D_X, D_Y );
        }
        else {

            /* If diffusion is turned off, set diffusion parameters to 0 */

            d_x = 0.0;
            d_y = 0.0;

        }
        
        /* Compute advection parameters */
        /* Is plume updraft on? */
        if ( UPDRAFT ) {

            /* Compute global advection velocities at mid time step */

            /* vGlob_x > 0 means left, < 0 means right [m/s]
             * vGlob_y > 0 means upwards, < 0 means downwards [m/s]
             * dTrav_x: distance traveled on the x-axis through advection [m]
             * dTrav_y: distance traveled on the y-axis through advection [m]
             */

            AdvGlobal( curr_Time_s - tInitial_s + dt/2.0, UPDRAFT_TIME, UPDRAFT_VEL, \
                       vGlob_x, vGlob_y, dTrav_x, dTrav_y );
        }
        else {

            /* If advection is turned off, set advection parameters to 0 */

            vGlob_x = 0;
            vGlob_y = 0;

        }

        /* Update diffusion and advection arrays */
        Solver.UpdateDiff ( d_x, d_y );
        /* Assume no plume advection */
        Solver.UpdateAdv  ( 0.0E+00, 0.0E+00 );
        /* Microphysics settling is considered for each bin independently */
        /* Update shear */
        Solver.UpdateShear( shear, m.y() );

        /* ======================================================================= */
        /* ----------------------------------------------------------------------- */
        /* ------------------------------- RUN SANDS ----------------------------- */
        /* ---------------- Spectral Advection aNd Diffusion Solver -------------- */
        /* ----------------------------------------------------------------------- */
        /* ======================================================================= */

#ifdef TIME_IT

        Stopwatch.Start( reset );

#endif /* TIME_IT */

        if ( TRANSPORT ) {

            if ( CHEMISTRY ) {
                /* Advection and diffusion of gas phase species */
                for ( N = 0; N < NVAR; N++ ) {
                    if ( N == ind_H2O )
                        Solver.Run( Data.Species[ind_H2Oplume], cellAreas, 1 );
                    else
                        Solver.Run( Data.Species[N], cellAreas );
                }
            } else {
                /* Advection and diffusion of condensable species */
                /* Advection and diffusion of plume affected H2O */
                Solver.Run( Data.Species[ind_H2Oplume], cellAreas, -1 );

            }

            /* Update H2O */
            for ( jNy = 0; jNy < NY; jNy++ ) {
                for ( iNx = 0; iNx < NX; iNx++ ) {
                    Data.Species[ind_H2O][jNy][iNx] = Data.Species[ind_H2Omet][jNy][iNx] + Data.Species[ind_H2Oplume][jNy][iNx];
                }
            }

            /* Advection and diffusion for aerosol particles */
            Solver.Run( Data.sootDens, cellAreas );
            /* Monodisperse assumption for soot particles */
            Solver.Run( Data.sootRadi, cellAreas );
            Solver.Run( Data.sootArea, cellAreas );

            /* We assume that sulfate aerosols do not settle */
            if ( TRANSPORT_LA ) {
                /* Transport of liquid aerosols */
                for ( UInt iBin_LA = 0; iBin_LA < Data.nBin_LA; iBin_LA++ )
                    Solver.Run( Data.liquidAerosol.pdf[iBin_LA], cellAreas );
            }

            if ( TRANSPORT_PA ) {
                /* Transport of solid aerosols */

                /* Ice volume per bin (NBIN x NY x NX) in [m^3/cm^3 air] */
                Vector_3D iceVolume = Data.solidAerosol.Volume();

                for ( UInt iBin_PA = 0; iBin_PA < Data.nBin_PA; iBin_PA++ ) {
                    /* Transport particle number and volume for each bin and
                     * recompute centers of each bin for each grid cell
                     * accordingly */
                    Solver.UpdateAdv ( 0.0E+00, vFall[iBin_PA] );

                    Solver.Run( Data.solidAerosol.pdf[iBin_PA], cellAreas, -1 );
                    Solver.Run( iceVolume[iBin_PA], cellAreas, -1 );

                }

                if ( FLUX_CORRECTION ) {

                    /* Limit flux of ice particles through top boundary */

#pragma omp parallel for                        \
                    if      ( !PARALLEL_CASES ) \
                    default ( shared          ) \
                    private ( iNx, jNy        ) \
                    schedule( dynamic, 1      )
                    for ( jNy = 0; jNy < NY; jNy++ ) {
                        if ( ( yE[jNy] > YLIM_UP - 200.0 ) && ( yE[jNy] > 400.0 ) ) {
                            for ( iNx = 0; iNx < NX; iNx++ ) {
                                for ( UInt iBin_PA = 0; iBin_PA < Data.nBin_PA; iBin_PA++ ) {
                                    Data.solidAerosol.pdf[iBin_PA][jNy][iNx] = 0.0E+00;
                                    iceVolume[iBin_PA][jNy][iNx] = 0.0E+00;
                                }
                                Data.Species[ind_H2O][jNy][iNx] = Data.Species[ind_H2O][jNy][LASTINDEX_SHEAR];
                            }
                        }
                    }

                    /* Limit flux of ice particles through left and right boundary */

#pragma omp parallel for                        \
                    if      ( !PARALLEL_CASES ) \
                    default ( shared          ) \
                    private ( iNx, jNy        ) \
                    schedule( dynamic, 1      )
                    for ( iNx = 0; iNx < NX; iNx++ ) {
#ifndef XLIM
                        if ( ( xE[iNx] < -XLIM_LEFT + 5.0E+03 ) || ( xE[iNx] > XLIM_RIGHT - 5.0E+03 ) ) {
#else
                        if ( ( xE[iNx] < -XLIM + 5.0E+03 ) || ( xE[iNx] > XLIM - 5.0E+03 ) ) {
#endif
                            for ( jNy = 0; jNy < NY; jNy++ ) {
                                for ( UInt iBin_PA = 0; iBin_PA < Data.nBin_PA; iBin_PA++ ) {
                                    Data.solidAerosol.pdf[iBin_PA][jNy][iNx] = 0.0E+00;
                                    iceVolume[iBin_PA][jNy][iNx] = 0.0E+00;
                                }
                                Data.Species[ind_H2O][jNy][iNx] = Data.Species[ind_H2O][NY-1][iNx];
                            }
                        }
                    }

                } /* FLUX_CORRECTION */

                /* Update centers of each bin */
                Data.solidAerosol.UpdateCenters( iceVolume, Data.solidAerosol.pdf );

            }

#ifdef RINGS

            /* If using rings and shear is non zero, then stretch rings to capture the
             * asymmetric expansion of the plume */
            if ( shear != 0.0E+00 ) {

                /* Rings do NOT get diffused, nor advected. They only get distorted
                 * through shear */

                /* Update diffusion and advection arrays */
                Solver.UpdateDiff ( 0.0E+00, 0.0E+00 );
                /* Assume no plume advection */
                Solver.UpdateAdv  ( 0.0E+00, 0.0E+00 );
                /* Update shear */
                Solver.UpdateShear( shear, m.y() );

                /* Do not apply any filling option: -1 */
                for ( iRing = 0; iRing < nRing + 1; iRing++ )
                    Solver.Run( m.weights[iRing], cellAreas, -1 );

                /* Recompute the map to mesh mapping, i.e. for each grid cell,
                 * find the corresponding ring */
                m.MapWeights();

                mapIndices = m.mapIndex();

            }

#endif /* RINGS */

        }

#ifdef TIME_IT

        Stopwatch.Stop( );
        SANDS_clock = Stopwatch.Elapsed( );

#endif /* TIME_IT */

        /* ======================================================================= */
        /* ----------------------------------------------------------------------- */
        /* -------------------------- UPDATE METEOROLOGY ------------------------- */
        /* ----------------------------------------------------------------------- */
        /* ======================================================================= */

        /* Met only matters for contrail evolution */
        if ( TRANSPORT_PA ) {
            /* Update met fields at mid time step */
            Met.Update( ( curr_Time_s + dt/2 ) / 3600.0, m, dTrav_x, dTrav_y );

            /* Update H2O */
            for ( jNy = 0; jNy < NY; jNy++ ) {
                for ( iNx = 0; iNx < NX; iNx++ ) {
                    Data.Species[ind_H2Omet][jNy][iNx] = Met.H2O(jNy,iNx);
                    Data.Species[ind_H2O][jNy][iNx] = Data.Species[ind_H2Omet][jNy][iNx] + Data.Species[ind_H2Oplume][jNy][iNx];
                }
            }
        }

        /* ======================================================================= */
        /* ----------------------------------------------------------------------- */
        /* ----------- UPDATE SOLAR ZENITH ANGLE AND PHOTOLYSIS RATES ------------ */
        /* ----------------------------------------------------------------------- */
        /* ======================================================================= */

        /* Compute the cosine of solar zenith angle midway through the integration step */
        sun->Update( curr_Time_s + dt/2 );

        /* Store cosine of solar zenith angle */
        ambientData.cosSZA[nTime] = sun->CSZA;

        if ( printDEBUG ) {
            std::cout << "\n DEBUG : \n";
            std::cout << "         CSZA = " << sun->CSZA << "\n";
        }

        /* Reset photolysis rates */
#pragma omp parallel for       \
        if ( !PARALLEL_CASES ) \
        default ( shared     ) \
        schedule( dynamic, 1 )
        for ( UInt iPhotol = 0; iPhotol < NPHOTOL; iPhotol++ )
            jRate[iPhotol] = 0.0E+00;

        /* If daytime, update photolysis rates */
        if ( sun->CSZA > 0.0E+00 )
            Update_JRates( jRate, sun->CSZA );

        if ( printDEBUG ) {
            std::cout << "\n DEBUG : \n";
            for ( UInt iPhotol = 0; iPhotol < NPHOTOL; iPhotol++ )
                std::cout << "         PHOTOL[" << iPhotol << "] = " << PHOTOL[iPhotol] << "\n";
        }


        /* ======================================================================= */
        /* ----------------------------------------------------------------------- */
        /* ------------------------------- RUN KPP ------------------------------- */
        /* ------------------------ The Kinetics Pre-Processor ------------------- */
        /* ----------------------------------------------------------------------- */
        /* ======================================================================= */


#ifdef TIME_IT

        Stopwatch.Start( reset );

#endif /* TIME_IT */

        /* Are we solving the chemistry in a ring structure? */
        #ifdef RINGS

            /* Fill in variables species for current time */
            ringData.FillIn( Data, m.weights, nTime + 1 );

            /* Is chemistry turned on? */
            if ( CHEMISTRY ) {

                /* In-ring chemistry */
                for ( iRing = 0; iRing < nRing; iRing++ ) {

                    /* Convert ring structure to KPP inputs (VAR and FIX) */
                    ringData.getData( nTime + 1, iRing );

                    for ( UInt iSpec = 0; iSpec < NVAR; iSpec++ )
                        tempArray[iSpec] = VAR[iSpec];

                    /* ===================================================== */
                    /* ================= Chemical rates ==================== */
                    /* ===================================================== */

                    /* Update heterogeneous chemistry reaction rates */
                    if ( HETCHEM ) {

                        for ( UInt iSpec = 0; iSpec < NSPEC; iSpec++ ) {
                            HET[iSpec][0] = 0.0E+00;
                            HET[iSpec][1] = 0.0E+00;
                            HET[iSpec][2] = 0.0E+00;
                        }

                        Data.getAerosolProp( AerosolRadi, AerosolArea, IWC,  \
                                             m.weights[iRing] );

                        relHumidity = VAR[ind_H2O] * \
                                      physConst::kB * temperature_K * 1.00E+06 / \
                                      physFunc::pSat_H2Ol( temperature_K );
                        GC_SETHET( temperature_K, pressure_Pa, airDens, relHumidity, \
                                   Data.STATE_PSC, VAR, AerosolArea, AerosolRadi, IWC, &(Data.KHETI_SLA[0]) );

                        if ( printDEBUG ) {
                            std::cout << "\n DEBUG :  Heterogeneous chemistry rates (Ring:  " << iRing << ")\n";
                            std::cout << "       :  Aerosol properties\n";
                            std::cout << "       :  Radius ice/NAT    = " << AerosolRadi[0] * 1.0E+06 << " [mum]\n";
                            std::cout << "       :  Radius strat. liq = " << AerosolRadi[1] * 1.0E+09 << " [nm]\n";
                            std::cout << "       :  Radius trop. sulf = " << AerosolRadi[2] * 1.0E+09 << " [nm]\n";
                            std::cout << "       :  Radius soot part. = " << AerosolRadi[3] * 1.0E+09 << " [nm]\n";
                            std::cout << "       :  Area ice/NAT      = " << AerosolArea[0] * 1.0E+12 << " [mum^2/cm^3]\n";
                            std::cout << "       :  Area strat. liq   = " << AerosolArea[1] * 1.0E+12 << " [mum^2/cm^3]\n";
                            std::cout << "       :  Area trop. sulf   = " << AerosolArea[2] * 1.0E+12 << " [mum^2/cm^3]\n";
                            std::cout << "       :  Area soot part.   = " << AerosolArea[3] * 1.0E+12 << " [mum^2/cm^3]\n";
                            std::cout << "       :  HET[ind_HO2][0]   = " << HET[ind_HO2][0]          << " [molec/cm^3/s]\n";
                            std::cout << "       :  HET[ind_NO2][0]   = " << HET[ind_NO2][0]          << " [molec/cm^3/s]\n";
                            std::cout << "       :  HET[ind_NO3][0]   = " << HET[ind_NO3][0]          << " [molec/cm^3/s]\n";
                            std::cout << "       :  HET[ind_N2O5][0]  = " << HET[ind_N2O5][0]         << " [molec/cm^3/s]\n";
                            std::cout << "       :  HET[ind_BrNO3][0] = " << HET[ind_BrNO3][0]        << " [molec/cm^3/s]\n";
                            std::cout << "       :  HET[ind_HOBr][0]  = " << HET[ind_HOBr][0]         << " [molec/cm^3/s]\n";
                            std::cout << "       :  HET[ind_HBr][0]   = " << HET[ind_HBr][0]          << " [molec/cm^3/s]\n";
                            std::cout << "       :  HET[ind_HOBr][1]  = " << HET[ind_HOBr][1]         << " [molec/cm^3/s]\n";
                            std::cout << "       :  PSC Rates:\n";
                            std::cout << "       :  HET[ind_N2O5][1]  = " << HET[ind_N2O5][1]         << " [molec/cm^3/s]\n";
                            std::cout << "       :  HET[ind_ClNO3][0] = " << HET[ind_ClNO3][0]        << " [molec/cm^3/s]\n";
                            std::cout << "       :  HET[ind_ClNO3][1] = " << HET[ind_ClNO3][1]        << " [molec/cm^3/s]\n";
                            std::cout << "       :  HET[ind_ClNO3][2] = " << HET[ind_ClNO3][2]        << " [molec/cm^3/s]\n";
                            std::cout << "       :  HET[ind_BrNO3][1] = " << HET[ind_BrNO3][1]        << " [molec/cm^3/s]\n";
                            std::cout << "       :  HET[ind_HOCl][0]  = " << HET[ind_HOCl][0]         << " [molec/cm^3/s]\n";
                            std::cout << "       :  HET[ind_HOCl][1]  = " << HET[ind_HOCl][1]         << " [molec/cm^3/s]\n";
                            std::cout << "       :  HET[ind_HOBr][2]  = " << HET[ind_HOBr][2]         << " [molec/cm^3/s]\n";
                        }
                    }

                    /* Zero-out reaction rate */
                    for ( UInt iReact = 0; iReact < NREACT; iReact++ )
                        RCONST[iReact] = 0.0E+00;

                    /* Update photolysis rates */
                    for ( UInt iPhotol = 0; iPhotol < NPHOTOL; iPhotol++ )
                        PHOTOL[iPhotol] = jRate[iPhotol];

                    /* Update reaction rates */
                    Update_RCONST( temperature_K, pressure_Pa, airDens, VAR[ind_H2O] );

                    if ( SAVE_PL ) {

                        RealDouble familyRate[NFAM];

                        for ( UInt iFam = 0; iFam < NFAM; iFam++ )
                            familyRate[iFam] = 0.0E+00;

                        /* If chemistry is performed within each rings, then the rates
                         * are stored as:
                         * NRING x (NT-1) x NFAM in [molec/cm^3/s]
                         * into the "forward" output file at a frequency specified by
                         * the input file "input.apcemm" */

                        /* Compute family rates */
                        ComputeFamilies( VAR, FIX, RCONST, familyRate );

                        for ( UInt iFam = 0; iFam < NFAM; iFam++ )
                            plumeRates[nTime][iRing][iFam] = familyRate[iFam];

                    } else {

                        if ( SAVE_O3PL ) {

                            RealDouble familyRate[NFAM];

                            for ( UInt iFam = 0; iFam < NFAM; iFam++ )
                                familyRate[iFam] = 0.0E+00;

                            /* If chemistry is performed within each rings, then the rates
                             * are stored as:
                             * NRING x (NT-1) x NFAM in [molec/cm^3/s]
                             * into the "forward" output file at a frequency specified by
                             * the input file "input.apcemm" */

                            /* Compute family rates */
                            ComputeFamilies( VAR, FIX, RCONST, familyRate );

                            for ( UInt iFam = 0; iFam < 2; iFam++ )
                                plumeRates[nTime][iRing][iFam] = familyRate[iFam];

                        }
                    }

                    /* ===================================================== */
                    /* ============== Chemical integration ================= */
                    /* ===================================================== */

                    IERR = INTEGRATE( VAR, curr_Time_s, curr_Time_s + dt, \
                                      ATOL, RTOL, STEPMIN );

                    if ( IERR < 0 ) {
                        /* Integration failed */

                        std::cout << "Integration failed";
                        #ifdef OMP
                            std::cout << " on " << omp_get_thread_num();
                        #endif /* OMP */
                        std::cout << " for ring = " << iRing << " at time t = " << curr_Time_s/3600.0 << " ( nTime = " << nTime << " )\n";

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
                        return KPP_FAIL;
                    }

                    ringData.FillIn( nTime + 1, iRing );

                    Data.applyRing( tempArray, mapIndices, iRing );

                }

                /* Ambient chemistry */
                ambientData.getData( aerArray, nTime );

                for ( UInt iSpec = 0; iSpec < NVAR; iSpec++ )
                    tempArray[iSpec] = VAR[iSpec];

                /* ========================================================= */
                /* =================== Chemical rates ====================== */
                /* ========================================================= */

                /* Update heterogeneous chemistry reaction rates */
                if ( HETCHEM ) {

                    for ( UInt iSpec = 0; iSpec < NSPEC; iSpec++ ) {
                        HET[iSpec][0] = 0.0E+00;
                        HET[iSpec][1] = 0.0E+00;
                        HET[iSpec][2] = 0.0E+00;
                    }

                    relHumidity = VAR[ind_H2O] * \
                                  physConst::kB * temperature_K * 1.00E+06 / \
                                  physFunc::pSat_H2Ol( temperature_K );
                    GC_SETHET( temperature_K, pressure_Pa, airDens, relHumidity, \
                               Data.STATE_PSC, VAR, AerosolArea, AerosolRadi, IWC, &(Data.KHETI_SLA[0]) );

                    if ( printDEBUG ) {
                        std::cout << "\n DEBUG :   Heterogeneous chemistry rates (Ambient)\n";
                        std::cout << "       :   Aerosol properties\n";
                        std::cout << "       :   Radius ice/NAT    = " << AerosolRadi[0] * 1.0E+06 << " [mum]\n";
                        std::cout << "       :   Radius strat. liq = " << AerosolRadi[1] * 1.0E+09 << " [nm]\n";
                        std::cout << "       :   Radius trop. sulf = " << AerosolRadi[2] * 1.0E+09 << " [nm]\n";
                        std::cout << "       :   Radius soot part. = " << AerosolRadi[3] * 1.0E+09 << " [nm]\n";
                        std::cout << "       :   Area ice/NAT      = " << AerosolArea[0] * 1.0E+12 << " [mum^2/cm^3]\n";
                        std::cout << "       :   Area strat. liq   = " << AerosolArea[1] * 1.0E+12 << " [mum^2/cm^3]\n";
                        std::cout << "       :   Area trop. sulf   = " << AerosolArea[2] * 1.0E+12 << " [mum^2/cm^3]\n";
                        std::cout << "       :   Area soot part.   = " << AerosolArea[3] * 1.0E+12 << " [mum^2/cm^3]\n";
                        std::cout << "       :   HET[ind_HO2][0]   = " << HET[ind_HO2][0]          << " [molec/cm^3/s]\n";
                        std::cout << "       :   HET[ind_NO2][0]   = " << HET[ind_NO2][0]          << " [molec/cm^3/s]\n";
                        std::cout << "       :   HET[ind_NO3][0]   = " << HET[ind_NO3][0]          << " [molec/cm^3/s]\n";
                        std::cout << "       :   HET[ind_N2O5][0]  = " << HET[ind_N2O5][0]         << " [molec/cm^3/s]\n";
                        std::cout << "       :   HET[ind_BrNO3][0] = " << HET[ind_BrNO3][0]        << " [molec/cm^3/s]\n";
                        std::cout << "       :   HET[ind_HOBr][0]  = " << HET[ind_HOBr][0]         << " [molec/cm^3/s]\n";
                        std::cout << "       :   HET[ind_HBr][0]   = " << HET[ind_HBr][0]          << " [molec/cm^3/s]\n";
                        std::cout << "       :   HET[ind_HOBr][1]  = " << HET[ind_HOBr][1]         << " [molec/cm^3/s]\n";
                        std::cout << "       :   PSC Rates:\n";
                        std::cout << "       :   HET[ind_N2O5][1]  = " << HET[ind_N2O5][1]         << " [molec/cm^3/s]\n";
                        std::cout << "       :   HET[ind_ClNO3][0] = " << HET[ind_ClNO3][0]        << " [molec/cm^3/s]\n";
                        std::cout << "       :   HET[ind_ClNO3][1] = " << HET[ind_ClNO3][1]        << " [molec/cm^3/s]\n";
                        std::cout << "       :   HET[ind_ClNO3][2] = " << HET[ind_ClNO3][2]        << " [molec/cm^3/s]\n";
                        std::cout << "       :   HET[ind_BrNO3][1] = " << HET[ind_BrNO3][1]        << " [molec/cm^3/s]\n";
                        std::cout << "       :   HET[ind_HOCl][0]  = " << HET[ind_HOCl][0]         << " [molec/cm^3/s]\n";
                        std::cout << "       :   HET[ind_HOCl][1]  = " << HET[ind_HOCl][1]         << " [molec/cm^3/s]\n";
                        std::cout << "       :   HET[ind_HOBr][2]  = " << HET[ind_HOBr][2]         << " [molec/cm^3/s]\n";
                    }
                }

                /* Zero-out reaction rate */
                for ( UInt iReact = 0; iReact < NREACT; iReact++ )
                    RCONST[iReact] = 0.0E+00;

                /* Update photolysis rates */
                for ( UInt iPhotol = 0; iPhotol < NPHOTOL; iPhotol++ )
                    PHOTOL[iPhotol] = jRate[iPhotol];

                /* Update reaction rates */
                Update_RCONST( temperature_K, pressure_Pa, airDens, VAR[ind_H2O] );

                if ( SAVE_PL ) {

                    RealDouble familyRate[NFAM];

                    for ( UInt iFam = 0; iFam < NFAM; iFam++ )
                        familyRate[iFam] = 0.0E+00;

                    /* If chemistry is performed within each rings, then the rates
                     * are stored as:
                     * NRING x (NT-1) x NFAM in [molec/cm^3/s]
                     * into the "forward" output file at a frequency specified by
                     * the input file "input.apcemm" */

                    /* Compute family rates */
                    ComputeFamilies( VAR, FIX, RCONST, familyRate );

                    for ( UInt iFam = 0; iFam < NFAM; iFam++ )
                        ambientRates[nTime][iFam] = familyRate[iFam];

                } else {

                    if ( SAVE_O3PL ) {

                        RealDouble familyRate[NFAM];

                        for ( UInt iFam = 0; iFam < NFAM; iFam++ )
                            familyRate[iFam] = 0.0E+00;

                        /* If chemistry is performed within each rings, then the rates
                         * are stored as:
                         * NRING x (NT-1) x NFAM in [molec/cm^3/s]
                         * into the "forward" output file at a frequency specified by
                         * the input file "input.apcemm" */

                        /* Compute family rates */
                        ComputeFamilies( VAR, FIX, RCONST, familyRate );

                        for ( UInt iFam = 0; iFam < 2; iFam++ )
                            ambientRates[nTime][iFam] = familyRate[iFam];

                    }
                }

                /* ========================================================= */
                /* ================ Chemical integration =================== */
                /* ========================================================= */

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
                    return KPP_FAIL;
                }

                ambientData.FillIn( nTime + 1 );

                Data.applyRing( tempArray, mapIndices, iRing );

            }

        #else

            /* Otherwise solve chemistry on the grid */

            /* Is chemistry turned on? */
            if ( CHEMISTRY ) {

                Vector_2D iceVolume_ = Data.solidAerosol.TotalVolume();

#pragma omp parallel for                             \
                if       ( !PARALLEL_CASES         ) \
                default ( shared                   ) \
                private ( iNx, jNy                 ) \
                private ( relHumidity, IWC         ) \
                schedule( dynamic, 1               )
                for ( iNx = 0; iNx < NX; iNx++ ) {
                    for ( jNy = 0; jNy < NY; jNy++ ) {

                        RealDouble AerosolArea[NAERO];
                        RealDouble AerosolRadi[NAERO];

                        /* Convert data structure to KPP inputs (VAR and FIX) */
                        Data.getData( iNx, jNy );

                        /* ================================================= */
                        /* =============== Chemical rates ================== */
                        /* ================================================= */

                        /* Update heterogeneous chemistry reaction rates */
                        if ( HETCHEM ) {

                            for ( UInt iSpec = 0; iSpec < NSPEC; iSpec++ ) {
                                HET[iSpec][0] = 0.0E+00;
                                HET[iSpec][1] = 0.0E+00;
                                HET[iSpec][2] = 0.0E+00;
                            }

                            relHumidity = VAR[ind_H2O] * \
                                          physConst::kB * Met.temp(jNy,iNx) * 1.00E+06 / \
                                          physFunc::pSat_H2Ol( Met.temp(jNy,iNx) );

                            /* Ice/NAT */
                            AerosolArea[0] = Data.solidAerosol.Moment( 2, jNy, iNx );
                            AerosolRadi[0] = std::max( std::min( Data.solidAerosol.Radius( jNy, iNx ), 1.00E-04 ), 1.00E-10 );

                            /* Stratospheric liquid aerosols */
                            AerosolArea[1] = Data.liquidAerosol.Moment( 2, jNy, iNx );
                            AerosolRadi[1] = std::max( std::min( Data.liquidAerosol.Radius( jNy, iNx ), 1.00E-06 ), 1.00E-10 );

                            /* Tropospheric aerosols.
                             * Zero it out */
                            AerosolArea[2] = 0.0E+00;
                            AerosolRadi[2] = 1.0E-07;

                            /* Black carbon */
                            AerosolArea[3] = Data.sootArea[jNy][iNx];
                            AerosolRadi[3] = std::max( std::min( Data.sootRadi[jNy][iNx], 1.00E-08 ), 1.00E-10 );

                            IWC            = Data.solidAerosol.Moment( 3, jNy, iNx ) \
                                           * physConst::RHO_ICE; /* [kg/cm^3] */

                            GC_SETHET( Met.temp(jNy,iNx), Met.press(jNy), \
                                       Met.airDens(jNy,iNx), relHumidity, \
                                       Data.STATE_PSC, VAR, AerosolArea,  \
                                       AerosolRadi, IWC, &(Data.KHETI_SLA[0]) );
                        }

                        /* Zero-out reaction rate */
                        for ( UInt iReact = 0; iReact < NREACT; iReact++ )
                            RCONST[iReact] = 0.0E+00;

                        /* Update photolysis rates */
                        for ( UInt iPhotol = 0; iPhotol < NPHOTOL; iPhotol++ )
                            PHOTOL[iPhotol] = jRate[iPhotol];

                        /* Update reaction rates */
                        Update_RCONST( Met.temp(jNy,iNx), Met.press(jNy), \
                                       Met.airDens(jNy,iNx), VAR[ind_H2O] );

                        /* ================================================= */
                        /* ============= Chemical integration ============== */
                        /* ================================================= */

                        IERR = INTEGRATE( VAR, curr_Time_s, curr_Time_s + dt, \
                                          ATOL, RTOL, STEPMIN );

                        if ( IERR < 0 ) {
                            /* Integration failed */

                            std::cout << "Integration failed";
                            #ifdef OMP
                                std::cout << " on " << omp_get_thread_num();
                            #endif /* OMP */
                            std::cout << " for grid cell = (" << jNy << ", " << iNx << ") at time t = " << curr_Time_s/3600.0 << " ( nTime = " << nTime << " )\n";

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

//                                /* Tweak would be to define a bool "stop", initialize to False and then set it as True if IERR < 0 */
//                                { /* Clear dynamically allocated variable(s) */
//                                    if ( sun != NULL ) sun->~SZA();
//                                    return KPP_FAIL;
//                                }
                        }

                        /* Convert KPP output back to data structure */
                        Data.applyData( iNx, jNy );

                    }
                }

                RealDouble AerosolArea[NAERO];
                RealDouble AerosolRadi[NAERO];

                /* Ambient chemistry */
                ambientData.getData( aerArray, nTime );

                /* ========================================================= */
                /* ==================== Chemical rates ===================== */
                /* ========================================================= */

                /* Update heterogeneous chemistry reaction rates */
                if ( HETCHEM ) {

                    for ( UInt iSpec = 0; iSpec < NSPEC; iSpec++ ) {
                        HET[iSpec][0] = 0.0E+00;
                        HET[iSpec][1] = 0.0E+00;
                        HET[iSpec][2] = 0.0E+00;
                    }

                    relHumidity = VAR[ind_H2O] * \
                                  physConst::kB * temperature_K * 1.00E+06 / \
                                  physFunc::pSat_H2Ol( temperature_K );
                    GC_SETHET( temperature_K, pressure_Pa, airDens, relHumidity, \
                               Data.STATE_PSC, VAR, AerosolArea, AerosolRadi, IWC, &(Data.KHETI_SLA[0]) );
                }

                /* Zero-out reaction rate */
                for ( UInt iReact = 0; iReact < NREACT; iReact++ )
                    RCONST[iReact] = 0.0E+00;

                /* Update photolysis rates */
                for ( UInt iPhotol = 0; iPhotol < NPHOTOL; iPhotol++ )
                    PHOTOL[iPhotol] = jRate[iPhotol];

                /* Update reaction rates */
                Update_RCONST( temperature_K, pressure_Pa, airDens, VAR[ind_H2O] );

                /* ========================================================= */
                /* ================= Chemical integration ================== */
                /* ========================================================= */

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
                    return KPP_FAIL;
                }

                ambientData.FillIn( nTime + 1 );
            }

        #endif /* RINGS */

#ifdef TIME_IT

        Stopwatch.Stop( );
        KPP_clock = Stopwatch.Elapsed( );

#endif /* TIME_IT */

        /* ======================================================================= */
        /* ----------------------------------------------------------------------- */
        /* ------------------------ AEROSOL MICROPHYSICS ------------------------- */
        /* ---------------------------- COAGULATION ------------------------------ */
        /* ----------------------------------------------------------------------- */
        /* ======================================================================= */

        ITS_TIME_FOR_LIQ_COAGULATION = ( ( ( curr_Time_s + dt - lastTimeLiqCoag ) >= COAG_DT * 60.0 ) || LAST_STEP );
        /* Liquid aerosol coagulation */
        if ( ITS_TIME_FOR_LIQ_COAGULATION && LIQ_COAG ) {
            dtLiqCoag = ( curr_Time_s + dt - lastTimeLiqCoag );
            if ( printDEBUG )
                std::cout << "\n DEBUG (Liquid Coagulation): Current time: " << ( curr_Time_s + dt - tInitial_s ) / 3600.0 << " hr. Last coagulation event was at: " << ( lastTimeLiqCoag - tInitial_s ) / 3600.0 << " hr. Running for " << dtLiqCoag << " s\n";

            lastTimeLiqCoag = curr_Time_s + dt;
            /* If shear = 0, take advantage of the symmetry around the Y-axis */
            Data.liquidAerosol.Coagulate( dtLiqCoag, Data.LA_Kernel, LA_MICROPHYSICS, ( shear == 0.0E+00 ) && ( XLIM_LEFT == XLIM_RIGHT ) );
        }

        ITS_TIME_FOR_ICE_COAGULATION = ( ( ( curr_Time_s + dt - lastTimeIceCoag ) >= COAG_DT * 60.0 ) || LAST_STEP );
        /* Solid aerosol coagulation */
        if ( ITS_TIME_FOR_ICE_COAGULATION && ICE_COAG ) {
            dtIceCoag = ( curr_Time_s + dt - lastTimeIceCoag );
            if ( printDEBUG )
                std::cout << "\n DEBUG (Solid Coagulation): Current time: " << ( curr_Time_s + dt - tInitial_s ) / 3600.0 << " hr. Last coagulation event was at: " << ( lastTimeIceCoag - tInitial_s ) / 3600.0 << " hr. Running for " << dtIceCoag << " s\n";

            lastTimeIceCoag = curr_Time_s + dt;
            /* If shear = 0, take advantage of the symmetry around the Y-axis */
            Data.solidAerosol.Coagulate ( dtIceCoag, Data.PA_Kernel, PA_MICROPHYSICS, ( shear == 0.0E+00 ) && ( XLIM_LEFT == XLIM_RIGHT ) );
        }

        /* ======================================================================= */
        /* ----------------------------------------------------------------------- */
        /* ------------------------ AEROSOL MICROPHYSICS ------------------------- */
        /* ------------------------- ICE CRYSTAL GROWTH -------------------------- */
        /* ----------------------------------------------------------------------- */
        /* ======================================================================= */

        /* TODO: For now perform growth at every time step */
        ITS_TIME_FOR_ICE_GROWTH = ( ( ( curr_Time_s - lastTimeIceGrowth ) >= 0 ) || LAST_STEP );
        /* Solid aerosol growth */
        if ( ITS_TIME_FOR_ICE_GROWTH && ICE_GROWTH ) {
            dtIceGrowth = ( curr_Time_s + dt - lastTimeIceGrowth );
            if ( printDEBUG )
                std::cout << "\n DEBUG (Solid Aerosol Growth): Current time: " << ( curr_Time_s - tInitial_s ) / 3600.0 << " hr. Last growth event was at: " << ( lastTimeIceGrowth - tInitial_s ) / 3600.0 << " hr. Running for " << dtIceGrowth << " s\n";

            lastTimeIceGrowth = curr_Time_s + dt;
            /* If shear = 0, take advantage of the symmetry around the Y-axis */
            Data.solidAerosol.Grow( dtIceGrowth, Data.Species[ind_H2O], Met.Temp(), Met.Press(), PA_MICROPHYSICS, ( shear == 0.0E+00 ) && ( XLIM_LEFT == XLIM_RIGHT ) );
        }

        /* ======================================================================= */
        /* ----------------------------------------------------------------------- */
        /* ------------------------ PERFORM MASS CHECKS -------------------------- */
        /* ----------------------------------------------------------------------- */
        /* ======================================================================= */

#if ( NOy_MASS_CHECK )

        /* Compute ambient concentrations */
        mass_Ambient_NOy = ambientData.Species[ind_NO][nTime+1]     + ambientData.Species[ind_NO2][nTime+1]   \
                         + ambientData.Species[ind_NO3][nTime+1]    + ambientData.Species[ind_HNO2][nTime+1]  \
                         + ambientData.Species[ind_HNO3][nTime+1]   + ambientData.Species[ind_HNO4][nTime+1]  \`
                         + 2*ambientData.Species[ind_N2O5][nTime+1] + ambientData.Species[ind_PAN][nTime+1]   \
                         + ambientData.Species[ind_MPN][nTime+1]    + ambientData.Species[ind_N][nTime+1]     \
                         + ambientData.Species[ind_PROPNN][nTime+1] + ambientData.Species[ind_BrNO2][nTime+1] \
                         + ambientData.Species[ind_BrNO3][nTime+1]  + ambientData.Species[ind_ClNO2][nTime+1] \
                         + ambientData.Species[ind_ClNO3][nTime+1]  + ambientData.Species[ind_PPN][nTime+1]   \
                         + ambientData.Species[ind_PRPN][nTime+1]   + ambientData.Species[ind_R4N1][nTime+1]  \
                         + ambientData.Species[ind_PRN1][nTime+1]   + ambientData.Species[ind_R4N2][nTime+1]  \
                         + 2*ambientData.Species[ind_N2O][nTime+1];

        /* Compute emitted */
        mass_Emitted_NOy = 0;
#pragma omp parallel for                \
        default  ( shared             ) \
        private  ( iNx, jNy           ) \
        reduction( +:mass_Emitted_NOy ) \
        schedule ( dynamic, 1         ) \
        if       ( !PARALLEL_CASES    )
        for ( iNx = 0; iNx < NX; iNx++ ) {
            for ( jNy = 0; jNy < NY; jNy++ ) {
                mass_Emitted_NOy += ( Data.Species[ind_NO][jNy][iNx]     + Data.Species[ind_NO2][jNy][iNx]   \
                                    + Data.Species[ind_NO3][jNy][iNx]    + Data.Species[ind_HNO2][jNy][iNx]  \
                                    + Data.Species[ind_HNO3][jNy][iNx]   + Data.Species[ind_HNO4][jNy][iNx]  \
                                    + 2*Data.Species[ind_N2O5][jNy][iNx] + Data.Species[ind_PAN][jNy][iNx]   \
                                    + Data.Species[ind_MPN][jNy][iNx]    + Data.Species[ind_N][jNy][iNx]     \
                                    + Data.Species[ind_PROPNN][jNy][iNx] + Data.Species[ind_BrNO2][jNy][iNx] \
                                    + Data.Species[ind_BrNO3][jNy][iNx]  + Data.Species[ind_ClNO2][jNy][iNx] \
                                    + Data.Species[ind_ClNO3][jNy][iNx]  + Data.Species[ind_PPN][jNy][iNx]   \
                                    + Data.Species[ind_PRPN][jNy][iNx]   + Data.Species[ind_R4N1][jNy][iNx]  \
                                    + Data.Species[ind_PRN1][jNy][iNx]   + Data.Species[ind_R4N2][jNy][iNx]  \
                                    + 2*Data.Species[ind_N2O][jNy][iNx]                         \
                                    - mass_Ambient_NOy ) * cellAreas[jNy][iNx];
            }
        }

        /* Print to console */
        std::cout << "\n\n    " << " *** NOy mass check: ";
        std::cout << "\n    "   << " ~~> Emitted NOy: ";
        std::cout << std::setw(6) << mass_Emitted_NOy * 1.0E+06 / physConst::Na * MW_N * 1.0E+06 << " [g(N)/km] ";
        /*                           [molec/cm3 * m2] * [m3/cm3]/ [molec/mole]  * [kg/mole]*[g/kg*m/km] = [g/km] */
        std::cout << std::endl;

#ifdef RINGS

        mass_Emitted_NOy_Rings = 0;
        for ( iRing = 0; iRing < nRing; iRing++ ) {
            mass_Emitted_NOy_Rings += ( ringData.NO[nTime+1][iRing]     \
                                      + ringData.NO2[nTime+1][iRing]    \
                                      + ringData.NO3[nTime+1][iRing]    \
                                      + ringData.HNO2[nTime+1][iRing]   \
                                      + ringData.HNO3[nTime+1][iRing]   \
                                      + ringData.HNO4[nTime+1][iRing]   \
                                      + 2*ringData.N2O5[nTime+1][iRing] \
                                      + ringData.PAN[nTime+1][iRing]    \
                                      + ringData.MPN[nTime+1][iRing]    \
                                      + ringData.N[nTime+1][iRing]      \
                                      + ringData.PROPNN[nTime+1][iRing] \
                                      + ringData.BrNO2[nTime+1][iRing]  \
                                      + ringData.BrNO3[nTime+1][iRing]  \
                                      + ringData.ClNO2[nTime+1][iRing]  \
                                      + ringData.ClNO3[nTime+1][iRing]  \
                                      + ringData.PPN[nTime+1][iRing]    \
                                      + ringData.PRPN[nTime+1][iRing]   \
                                      + ringData.R4N1[nTime+1][iRing]   \
                                      + ringData.PRN1[nTime+1][iRing]   \
                                      + ringData.R4N2[nTime+1][iRing]   \
                                      + 2*ringData.N2O[nTime+1][iRing]  \
                                      - mass_Ambient_NOy ) * ringArea[iRing];
        }
        /* How much of this emitted mass is still in the rings? FR = Fraction in rings */
        std::cout << "(FR: " << 100 * mass_Emitted_NOy_Rings / mass_Emitted_NOy << " %)";

#endif /* RINGS */

#endif /* NOy_MASS_CHECK */

#if ( CO2_MASS_CHECK )

        /* CO2 is not an exactly conserved quantity because of the oxidation CO and other compounds (unless chemistry is turned off) */

        mass_Ambient_CO2 = ambientData.Species[ind_CO2][nTime+1];

        /* Compute emitted */
        mass_Emitted_CO2 = 0;
#pragma omp parallel for                \
        default  ( shared             ) \
        private  ( iNx, jNy           ) \
        reduction( +:mass_Emitted_CO2 ) \
        schedule ( dynamic, 1         ) \
        if       ( !PARALLEL_CASES    )
        for ( iNx = 0; iNx < NX; iNx++ ) {
            for ( jNy = 0; jNy < NY; jNy++ ) {
                mass_Emitted_CO2 += ( Data.Species[ind_CO2][jNy][iNx] \
                                    - mass_Ambient_CO2 ) * cellAreas[jNy][iNx];
            }
        }

        std::cout << "\n\n    " << " *** CO2 mass check: ";

        std::cout << "\n    " << " ~~> Emitted CO2: ";
        std::cout << std::setw(6) << mass_Emitted_CO2 * 1.0E+06 / physConst::Na * MW_CO2 * 1.0E+03 << " [kg/km]   ";
        /*                           [molec/cm3 * m2] * [m3/cm3]/ [molec/mole]  *[kg/mole]*[m/km] = [kg/km] */
        std::cout << std::endl;

#ifdef RINGS
            mass_Emitted_CO2_Rings = 0;
            for ( iRing = 0; iRing < nRing; iRing++ ) {
                mass_Emitted_CO2_Rings += ( ringData.Species[ind_CO2][nTime+1][iRing] \
                                          - mass_Ambient_CO2 ) * ringArea[iRing];
            }
            /* How much of this emitted mass is still in the rings? FR = Fraction in rings */
            std::cout << "(FR: " << 100 * mass_Emitted_CO2_Rings / mass_Emitted_CO2 << " %)\n";

#endif /* RINGS */

#endif /* CO2_MASS_CHECK */

#if ( H2O_MASS_CHECK )

        mass_Ambient_H2O = ambientData.Species[ind_H2O][nTime+1];
        totIceVol = Data.solidAerosol.TotalVolume( );

        /* Compute total water */
        mass_H2O = 0;
#pragma omp parallel for             \
        default  ( shared          ) \
        private  ( iNx, jNy        ) \
        reduction( +:mass_H2O      ) \
        schedule ( dynamic, 1      ) \
        if       ( !PARALLEL_CASES )
        for ( iNx = 0; iNx < NX; iNx++ ) {
            for ( jNy = 0; jNy < NY; jNy++ ) {
                mass_H2O += ( ( Data.Species[ind_H2O][jNy][iNx] - mass_Ambient_H2O ) \
                            + totIceVol[jNy][iNx] * UNITCONVERSION ) * \
                            cellAreas[jNy][iNx];
            }
        }

        std::cout << "\n\n    " << " *** H2O mass check: ";
        std::cout << "\n    " << " ~~> H2O: ";
        std::cout << std::setw(6) << mass_H2O      * 1.0E+06 / physConst::Na * MW_H2O * 1.0E+03 << " [kg/km]   ";
        /*                          [molec/cm3*m2] * [m3/cm3]/ [molec/mole]  *[kg/mole]*[m/km] = [kg/km] */
        std::cout << std::endl;

#endif /* H2O_MASS_CHECK */

#ifdef TIME_IT

        SANDS_clock_cumul += SANDS_clock;
        KPP_clock_cumul   += KPP_clock;
        std::cout << "\n    " << " *** Clock breakdown: ";
        std::cout << "\n    " << " *** ----------------- ";
        std::cout << "\n    " << " *** Total: " << SANDS_clock + KPP_clock << " [ms]";
        std::cout << " ( SANDS: ";
        std::cout << 100 * ( SANDS_clock / RealDouble( SANDS_clock + KPP_clock ) ) << "%";
        std::cout << ", KPP: " << 100 * ( KPP_clock / RealDouble( SANDS_clock + KPP_clock ) ) << "% )";

#endif /* TIME_IT */

        curr_Time_s += dt;
        nTime++;

        /* Timeseries diagnostics */
        if ( TS_SPEC && \
           (( TS_FREQ == 0 ) || \
            ( std::fmod((curr_Time_s - timeArray[0])/60.0, TS_FREQ) == 0.0E+00 )) ) {
           int hh = (int) (curr_Time_s - timeArray[0])/3600;
           int mm = (int) (curr_Time_s - timeArray[0])/60   - 60 * hh;
           int ss = (int) (curr_Time_s - timeArray[0])      - 60 * ( mm + 60 * hh );
           Diag_TS_Chem( TS_SPEC_FILENAME, TS_SPEC_LIST, hh, mm, ss, \
                         Data, m );
        }

        if ( TS_AERO && \
           (( TS_AERO_FREQ == 0 ) || \
            ( std::fmod((curr_Time_s - timeArray[0])/60.0, TS_AERO_FREQ) == 0.0E+00 )) ) {
            int hh = (int) (curr_Time_s - timeArray[0])/3600;
            int mm = (int) (curr_Time_s - timeArray[0])/60   - 60 * hh;
            int ss = (int) (curr_Time_s - timeArray[0])      - 60 * ( mm + 60 * hh );
            Diag_TS_Phys( TS_AERO_FILENAME, TS_AERO_LIST, hh, mm, ss, \
                          Data, m, Met );    
            float totalIceParticles = Data.solidAerosol.TotalNumber_sum( cellAreas );
            float totalIceMass = Data.solidAerosol.TotalIceMass_sum( cellAreas );
            if ( totalIceParticles <= 1.00E+1 && totalIceMass <= 1.00E-5 && !CHEMISTRY ) {
                std::cout << "EndSim: no particles remain" << std::endl;
                std::cout << "# ice particles: " << totalIceParticles << std::endl;
                std::cout << "Total ice mass [g]: " << totalIceMass << std::endl;
                exit(0);
            }
        }

    }
 
    /* ===================================================================== */
    /* --------------------------------------------------------------------- */
    /* ------------------------ TIME LOOP ENDS HERE ------------------------ */
    /* --------------------------------------------------------------------- */
    /* ===================================================================== */


#ifdef TIME_IT

    Stopwatch_cumul.Stop( );
    clock_cumul = Stopwatch_cumul.Elapsed( );

    std::cout << "\n";
    std::cout << " ** Final clock breakdown: " << "\n";

    std::cout << " ** -> SANDS: ";
    std::cout << std::setw(6) <<SANDS_clock_cumul / RealDouble(1000);
    std::cout << " [s] (" << 100 * SANDS_clock_cumul / RealDouble(clock_cumul) << " %)";
    std::cout << std::endl;

    std::cout << " ** -> KPP  : ";
    std::cout << std::setw(6) << KPP_clock_cumul / RealDouble(1000) << " [s] ";
    std::cout << "(" << 100 * KPP_clock_cumul / RealDouble(clock_cumul) << " %)";
    std::cout << std::endl;

    std::cout << " ** -> Rem. : ";
    std::cout << std::setw(6) << ( clock_cumul - SANDS_clock_cumul - KPP_clock_cumul ) / RealDouble(1000) << " [s] ";
    std::cout << "(" << 100 * ( clock_cumul - SANDS_clock_cumul - KPP_clock_cumul ) / RealDouble(clock_cumul) << " %)" << "\n";
    std::cout << std::endl;

    std::cout << " ** ----------------- " << "\n";
    std::cout << " ** Total   : ";
    std::cout << std::setw(6) << clock_cumul / RealDouble(1000) << " [s]" << "\n";
    std::cout << std::endl;

#endif /* TIME_IT */


#ifdef RINGS

#pragma omp critical
    {
        if ( SAVE_FORWARD ) {
            isSaved = output::Write( input.fileName2char(),               \
                                     Input_Opt,                           \
                                     TS_SPEC_LIST,                        \
                                     ringData,                            \
                                     ambientData,                         \
                                     ringCluster,                         \
                                     timeArray,                           \
                                     input,                               \
                                     airDens, relHumidity_i,              \
                                     sun->sunRise, sun->sunSet,           \
                                     plumeRates, ambientRates );
        }
    }

    if ( isSaved == output::SAVE_FAILURE ) {
        std::cout << " Saving to ring-averaged concentrations to file failed...\n";
        return SAVE_FAIL;
    }
#endif /* RINGS */

    /* ======================================================================= */
    /* ----------------------------------------------------------------------- */
    /* -------------------------- PERFORM ADJOINT ---------------------------- */
    /* ----------------------------------------------------------------------- */
    /* ======================================================================= */

#ifdef RINGS

    if ( ADJOINT ) {
        #ifdef OMP
            #pragma omp critical
            { std::cout << "\n\n ## ON THREAD " << omp_get_thread_num() << ": Starting adjoint calculation...\n"; }
        #else
            std::cout << "\n\n Starting adjoint calculation...\n";
        #endif /* OMP */

        std::copy(sun->CSZA_Vector.begin(), sun->CSZA_Vector.end(), SZA_CST);
        RealDouble VAR_OPT[NVAR];
        RealDouble METRIC;

        const Vector_1D initBackg = ringData.RingAverage( ringArea, totArea, 0 );
        const Vector_1D finalPlume = ringData.RingAverage( ringArea, totArea, timeArray.size() - 1 );

        IERR = KPP_Main_ADJ( &(finalPlume)[0], &(initBackg)[0],   \
                             temperature_K, pressure_Pa, airDens, \
                             &(timeArray)[0], timeArray.size(),   \
                             KPPADJ_RTOLS, KPPADJ_ATOLS,          \
                             /* Output */ VAR_OPT,                \
                             /* Output metric */ &METRIC,         \
                             /* Debug? */ DEBUG_ADJOINT,          \
                             /* 2nd try? */ 0 );

        if ( IERR < 0 ) {
            /* Integration succeeded but convergence was poor. Try again with
             * new initial direction */

            RealDouble VAR_OPT2[NVAR];
            RealDouble METRIC2;
            IERR = KPP_Main_ADJ( &(finalPlume)[0], &(initBackg)[0],   \
                                 temperature_K, pressure_Pa, airDens, \
                                 &(timeArray)[0], timeArray.size(),   \
                                 KPPADJ_RTOLS, KPPADJ_ATOLS,          \
                                 /* Output */ VAR_OPT2,               \
                                 /* Output metric */ &METRIC2,        \
                                 /* Debug? */ DEBUG_ADJOINT,          \
                                 /* 2nd try? */ 1 );

            /* If 2nd try provides a better metric, store the new optimized
             * array; otherwise stick to the original one. */
            if ( METRIC2 < METRIC ) {
                for ( UInt iSpec = 0; iSpec < NVAR; iSpec++ )
                    VAR_OPT[iSpec] = VAR_OPT2[iSpec];
            } else {
                #ifdef OMP
                    #pragma omp critical
                    { std::cout << "\n ## ON THREAD " << omp_get_thread_num() << ": Sticking to original solution.\n"; }
                #else
                    std::cout << "\n Sticking to original optimization solution.\n";
                #endif /* OMP */
            }
            /* This should be changed eventually */
            IERR = 0;

        }

        if ( IERR < 0 ) {
            /* Adjoint integration failed */
            return KPPADJ_FAIL;
        }

        /* Apply optimized initial conditions to concentration array */
        for ( UInt iSpec = 0; iSpec < NVAR; iSpec++ )
            VAR[iSpec] = VAR_OPT[iSpec];

        /* Create ambient struture */
        Ambient adjointData( timeArray.size(), Data.getAmbient(), Data.getAerosol(), Data.getLiqSpecies() );
        adjointData.FillIn( 0 );


        /* Perform forward integration with optimized initial conditions */
        const int NADJ = NVAR;
        RealDouble Y_adj[NADJ][NVAR];

        /* Allocate tolerances */

        RealDouble RTOL_adj[NADJ][NVAR];
        RealDouble ATOL_adj[NADJ][NVAR];



        /* ---- TOLERANCES ---------------------- */

        for( UInt i = 0; i < NVAR; i++ ) {
            RTOL[i] = KPPADJ_RTOLS;
            ATOL[i] = KPPADJ_ATOLS;
        }

        /* Tolerances for calculating adjoints are
         * used for controlling adjoint truncation
         * error and for solving the linear adjoint
         * equations by iterations.
         * Note: Adjoints typically span many orders
         * of magnitude and a careful tuning of
         * ATOL_adj may be necessary */

        for( UInt i = 0; i < NADJ; i++ ) {
            for( UInt j = 0; j < NVAR; j++ ) {
                RTOL_adj[i][j] = 1.0e-5;
                ATOL_adj[i][j] = 1.0e-10;
            }
        }

        /* ICNTRL , RCNTRL  = Adjoint settings
         * ISTATUS, RSTATUS = Adjoint statistics */

        RealDouble RCNTRL[20], RSTATUS[20];
        int ICNTRL[20], ISTATUS[20];

        /* Default control options */
        for( UInt i = 0; i < 20; i++ ) {
            ICNTRL[i] = 0;
            RCNTRL[i] = (RealDouble)0.0;
            ISTATUS[i] = 0;
            RSTATUS[i] = (RealDouble)0.0;
        }

        ICNTRL[6] = 1;

        /* Declare aerosol quantities */
        RealDouble relHumidity;
        RealDouble AerosolArea[NAERO];
        RealDouble AerosolRadi[NAERO];
        RealDouble IWC = 0;

        for ( UInt iAero = 0; iAero < NAERO; iAero++ ) {
            AerosolArea[iAero] = 0.0E+00;
            AerosolRadi[iAero] = 0.0E+00;
        }

        /* Reinitialize time */
        curr_Time_s = tInitial_s; /* [s] */
        nTime = 0;

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
                    std::cout << " ( on thread " << omp_get_thread_num() << " )";
                #endif /* OMP */
                std::cout << "\n -> Solar time: " << std::fmod( curr_Time_s/3600.0, 24.0 ) << " [hr]";
            }

            dt = timeArray[nTime+1] - timeArray[nTime];

            /* ======================================================================= */
            /* ----------------------------------------------------------------------- */
            /* ----------- UPDATE SOLAR ZENITH ANGLE AND PHOTOLYSIS RATES ------------ */
            /* ----------------------------------------------------------------------- */
            /* ======================================================================= */

            /* Compute the cosine of solar zenith angle midway through the integration step */
            sun->Update( curr_Time_s + dt/2 );

            /* Store cosine of solar zenith angle */
            adjointData.cosSZA[nTime] = sun->CSZA;

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
                    std::cout << "         PHOTOL[" << iPhotol << "] = " << PHOTOL[iPhotol] << "\n";
            }

            /* ============================================================= */
            /* ------------------------------------------------------------- */
            /* -------------------------- RUN KPP -------------------------- */
            /* ------------------- The Kinetics Pre-Processor -------------- */
            /* ------------------------------------------------------------- */
            /* ============================================================= */


            /* Ambient chemistry */
            adjointData.getData( aerArray, nTime );

            /* ============================================================= */
            /* ====================== Chemical rates ======================= */
            /* ============================================================= */

            /* Update heterogeneous chemistry reaction rates */
            if ( HETCHEM ) {

                for ( UInt iSpec = 0; iSpec < NSPEC; iSpec++ ) {
                    HET[iSpec][0] = 0.0E+00;
                    HET[iSpec][1] = 0.0E+00;
                    HET[iSpec][2] = 0.0E+00;
                }

                relHumidity = VAR[ind_H2O] * \
                              physConst::kB * temperature_K * 1.00E+06 / \
                              physFunc::pSat_H2Ol( temperature_K );
                GC_SETHET( temperature_K, pressure_Pa, airDens, relHumidity, \
                           Data.STATE_PSC, VAR, AerosolArea, AerosolRadi, IWC, &(Data.KHETI_SLA[0]) );

                if ( printDEBUG ) {
                    std::cout << "\n DEBUG :   Heterogeneous chemistry rates (Ambient)\n";
                    std::cout << "       :   Aerosol properties\n";
                    std::cout << "       :   Radius ice/NAT    = " << AerosolRadi[0] * 1.0E+06 << " [mum]\n";
                    std::cout << "       :   Radius strat. liq = " << AerosolRadi[1] * 1.0E+09 << " [nm]\n";
                    std::cout << "       :   Radius trop. sulf = " << AerosolRadi[2] * 1.0E+09 << " [nm]\n";
                    std::cout << "       :   Radius soot part. = " << AerosolRadi[3] * 1.0E+09 << " [nm]\n";
                    std::cout << "       :   Area ice/NAT      = " << AerosolArea[0] * 1.0E+12 << " [mum^2/cm^3]\n";
                    std::cout << "       :   Area strat. liq   = " << AerosolArea[1] * 1.0E+12 << " [mum^2/cm^3]\n";
                    std::cout << "       :   Area trop. sulf   = " << AerosolArea[2] * 1.0E+12 << " [mum^2/cm^3]\n";
                    std::cout << "       :   Area soot part.   = " << AerosolArea[3] * 1.0E+12 << " [mum^2/cm^3]\n";
                    std::cout << "       :   HET[ind_HO2][0]   = " << HET[ind_HO2][0]          << " [molec/cm^3/s]\n";
                    std::cout << "       :   HET[ind_NO2][0]   = " << HET[ind_NO2][0]          << " [molec/cm^3/s]\n";
                    std::cout << "       :   HET[ind_NO3][0]   = " << HET[ind_NO3][0]          << " [molec/cm^3/s]\n";
                    std::cout << "       :   HET[ind_N2O5][0]  = " << HET[ind_N2O5][0]         << " [molec/cm^3/s]\n";
                    std::cout << "       :   HET[ind_BrNO3][0] = " << HET[ind_BrNO3][0]        << " [molec/cm^3/s]\n";
                    std::cout << "       :   HET[ind_HOBr][0]  = " << HET[ind_HOBr][0]         << " [molec/cm^3/s]\n";
                    std::cout << "       :   HET[ind_HBr][0]   = " << HET[ind_HBr][0]          << " [molec/cm^3/s]\n";
                    std::cout << "       :   HET[ind_HOBr][1]  = " << HET[ind_HOBr][1]         << " [molec/cm^3/s]\n";
                    std::cout << "       :   PSC Rates:\n";
                    std::cout << "       :   HET[ind_N2O5][1]  = " << HET[ind_N2O5][1]         << " [molec/cm^3/s]\n";
                    std::cout << "       :   HET[ind_ClNO3][0] = " << HET[ind_ClNO3][0]        << " [molec/cm^3/s]\n";
                    std::cout << "       :   HET[ind_ClNO3][1] = " << HET[ind_ClNO3][1]        << " [molec/cm^3/s]\n";
                    std::cout << "       :   HET[ind_ClNO3][2] = " << HET[ind_ClNO3][2]        << " [molec/cm^3/s]\n";
                    std::cout << "       :   HET[ind_BrNO3][1] = " << HET[ind_BrNO3][1]        << " [molec/cm^3/s]\n";
                    std::cout << "       :   HET[ind_HOCl][0]  = " << HET[ind_HOCl][0]         << " [molec/cm^3/s]\n";
                    std::cout << "       :   HET[ind_HOCl][1]  = " << HET[ind_HOCl][1]         << " [molec/cm^3/s]\n";
                    std::cout << "       :   HET[ind_HOBr][2]  = " << HET[ind_HOBr][2]         << " [molec/cm^3/s]\n";
                }
            }

            /* Zero-out reaction rate */
            for ( UInt iReact = 0; iReact < NREACT; iReact++ )
                RCONST[iReact] = 0.0E+00;

            /* Update photolysis rates */
            for ( UInt iPhotol = 0; iPhotol < NPHOTOL; iPhotol++ )
                PHOTOL[iPhotol] = jRate[iPhotol];

            /* Update reaction rates */
            Update_RCONST( temperature_K, pressure_Pa, airDens, VAR[ind_H2O] );

            /* ============================================================= */
            /* =================== Chemical integration ==================== */
            /* ============================================================= */

            IERR = INTEGRATE_ADJ( NADJ, VAR, Y_adj, timeArray[nTime], timeArray[nTime+1], ATOL_adj, RTOL_adj, ATOL, RTOL, ICNTRL, RCNTRL, ISTATUS, RSTATUS, STEPMIN );

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
                return KPP_FAIL;
            }

            adjointData.FillIn( nTime + 1 );

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
            isSaved = output::Write_Adjoint( input.fileName_ADJ2char(), \
                                             TS_SPEC_LIST,              \
                                             ringData, ambientData,     \
                                             adjointData,               \
                                             ringArea, totArea,         \
                                             timeArray,                 \
                                             input,                     \
                                             airDens, relHumidity_i );
        }
        if ( isSaved == output::SAVE_FAILURE ) {
            std::cout << " Saving to adjoint data to file failed...\n";
            return SAVE_FAIL;
        }

    }

#endif /* RINGS */

    /* Clear dynamically allocated variable(s) */
    if ( sun != NULL )
        sun->~SZA();

    return SUCCESS;

} /* End of PlumeModel */

/* End of PlumeModel.cpp */
