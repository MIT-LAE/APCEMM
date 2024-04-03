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
#include <memory>
#include <sys/stat.h>
#include <fftw3.h>
#ifdef OMP
    #include "omp.h"
#endif /* OMP */

#include "APCEMM.h"
#include "Util/ForwardDecl.hpp"
#include "Core/Input_Mod.hpp"
#include "Core/Parameters.hpp"
#include "Core/Interface.hpp"
#include "Core/Input.hpp"
#include "Core/Monitor.hpp"
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
#include "Core/Fuel.hpp"
#include "Core/Engine.hpp"
#include "Core/Aircraft.hpp"
#include "Core/Emission.hpp"
#include "Core/ReadJRates.hpp"
#include "Core/TimestepVarsWrapper.hpp"
#include "Core/MPMSimVarsWrapper.hpp"
#include "FVM_ANDS/FVM_Solver.hpp"
#include "Util/PlumeModelUtils.hpp"
#include "Util/VectorUtils.hpp"

/* For RINGS */
#include "Core/Cluster.hpp"
#include "Core/Species.hpp"

/* For DIAGNOSTIC */
#include "Core/Diag_Mod.hpp"

#include "Core/Status.hpp"

/* (URGENT) FIXME: Whoever works on chemistry, REFACTOR THESE REMAINING GLOBAL VARIABLES by refactoring!
            Using non-const global variables is the biggest software engineering antipattern ever */
int isSaved = 1;
static int SAVE_FAIL   = -2;

double RCONST[NREACT];       /* Rate constants (global) */
double NOON_JRATES[NPHOTOL]; /* Noon-time photolysis rates (global) */
double PHOTOL[NPHOTOL];      /* Photolysis rates (global) */
double HET[NSPEC][3];        /* Heterogeneous chemistry rates (global) */

double TIME;                 /* Current integration time (global) */

/* Require this for adjoint integration */
double SZA_CST[3];

double totalH2OMass(const Solution& Data, const Vector_2D& cellAreas){

    double totIceMass = Data.solidAerosol.TotalIceMass_sum(cellAreas); //kg/m
    std::cout << std::setprecision(8);

    /* Compute total water */
    double mass_H2O = 0;
    int nx = cellAreas[0].size();
    int ny = cellAreas.size();
    for (int iNx = 0; iNx < nx; iNx++ ) {
        for (int jNy = 0; jNy < ny; jNy++ ) {
            mass_H2O += Data.Species[ind_H2O][jNy][iNx] * 1.0e6 * cellAreas[jNy][iNx] / physConst::Na * MW_H2O;
                        // #/cm3 * cm3/m3 * m2 = #/m
                        // #/m * mol/# * kg/mol = kg/m

        }
    }
    return totIceMass + mass_H2O;
        
}
SimStatus PlumeModel( OptInput &Input_Opt, const Input &input )
{
    auto start = std::chrono::high_resolution_clock::now();

    double C[NSPEC];             /* Concentration of all species */
    double * VAR = &C[0];        /* Concentration of variable species */
    double * FIX = &C[NVAR];     /* Concentration of fixed species */

    omp_set_num_threads(Input_Opt.SIMULATION_OMP_NUM_THREADS);
    bool printDEBUG = false;

#ifdef DEBUG

    std::cout << "\n DEBUG is turned ON!\n\n";
    printDEBUG = true;

#endif /* DEBUG */

    MPMSimVarsWrapper simVars = MPMSimVarsWrapper(input, Input_Opt);

    if ( simVars.TS_SPEC )
        std::cout << "\n Saving TS files to: " << simVars.TS_SPEC_FILEPATH << std::endl;

    if ( simVars.TS_AERO )
        std::cout << "\n Saving TS_AERO files to: " << simVars.TS_AERO_FILEPATH << std::endl;
    
    if ( ( simVars.TS_SPEC || simVars.TS_AERO ) && ( simVars.TS_FOLDER.compare("") != 0 ) ) {
        std::cout << "Creating directory" << std::endl;
        std::cout << simVars.TS_FOLDER << std::endl;
        /* Create output directory for timeseries */
        struct stat sb;
        if ( !( stat( simVars.TS_FOLDER.c_str(), &sb ) == 0 \
                    && S_ISDIR(sb.st_mode) ) ) {

            /* Create directory */
            const int dir_err = \
                    mkdir( simVars.TS_FOLDER.c_str(), \
                            S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );

            if ( dir_err == -1 ) {
                std::cout << " Could not create directory: ";
                std::cout << simVars.TS_FOLDER << std::endl;
                std::cout << " You may not have write permission" << std::endl;
                exit(1);
            }
        }
    }

    /* Grid indices */
    UInt iNx = 0;
    UInt jNy = 0;

    /* Species index */
    UInt N   = 0;

    int IERR;

#if ( NOy_MASS_CHECK )

    double mass_Ambient_NOy, mass_Emitted_NOy;

#endif /* NOy_MASS_CHECK */

#if ( CO2_MASS_CHECK )

    double mass_Ambient_CO2, mass_Emitted_CO2;

#endif /* CO2_MASS_CHECK */

#if ( H2O_MASS_CHECK )

    /* Conversion factor from ice volume [m^3] to [molecules] */
    const double UNITCONVERSION = physConst::RHO_ICE / MW_H2O * physConst::Na;
    /* Unit check: [kg/m^3] / [kg/mol] * [molec/mol] = [molec/m^3] */

    double mass_Ambient_H2O, mass_H2O;

    Vector_2D totIceVol;

#endif /* H2O_MASS_CHECK */

    /* ======================================================================= */
    /* ----------------------------------------------------------------------- */
    /* --------------------------------- MESH -------------------------------- */
    /* ----------------------------------------------------------------------- */
    /* ======================================================================= */

    Mesh m(Input_Opt);
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
    double jRate[NPHOTOL];

    /* Define sun parameters, this include sunrise and sunset hours and updates
     * the local solar zenith angle. */
    std::unique_ptr<SZA> sun = std::make_unique<SZA>(input.latitude_deg(), input.emissionDOY() );

    /* Initialize noon time photolysis rates
     * The data is after the quantum yield has been applied and represents
     * the value of the photolysis rates at 12:00 (noon) locally.
     * The photolysis rates at any given time are obtained by multiplying
     * those by the cosine of the solar zenith angle, when positive. */
    for ( UInt iPhotol = 0; iPhotol < NPHOTOL; iPhotol++ )
        NOON_JRATES[iPhotol] = 0.0E+00;

    /* Allocating noon-time photolysis rates. */

    if ( simVars.CHEMISTRY ) {
        #pragma omp critical
        {
            ReadJRates( simVars.JRATE_FOLDER.c_str(),  \
                input.emissionMonth(), \
                input.emissionDay(),   \
                input.longitude_deg(), \
                input.latitude_deg(),  \
                simVars.pressure_Pa/100.0,     \
                NOON_JRATES );
        }

    }

    /* ======================================================================= */
    /* ----------------------------------------------------------------------- */
    /* ------------------------------ TIMESTEPS ------------------------------ */
    /* ----------------------------------------------------------------------- */
    /* ======================================================================= */

    TimestepVarsWrapper timestepVars(input, Input_Opt);
    timestepVars.setTimeArray(PlumeModelUtils::BuildTime ( timestepVars.tInitial_s, timestepVars.tFinal_s, 3600.0*sun->sunRise, 3600.0*sun->sunSet, timestepVars.dt ));

    /* ======================================================================= */
    /* ----------------------------------------------------------------------- */
    /* ----------------------------- METEOROLOGY ----------------------------- */
    /* ----------------------------------------------------------------------- */
    /* ======================================================================= */

    AmbientMetParams ambMetParams;
    ambMetParams.solarTime_h = timestepVars.curr_Time_s / 3600.0;
    ambMetParams.rhi = simVars.relHumidity_i;
    ambMetParams.temp_K = simVars.temperature_K; 
    ambMetParams.press_Pa = simVars.pressure_Pa;
    ambMetParams.shear = input.shear();
    
    Meteorology Met( Input_Opt, ambMetParams, m.y(), m.yE());

    if ( Input_Opt.MET_LOADMET && Input_Opt.MET_LOADTEMP ) {
        simVars.temperature_K = Met.tempRef();
    }
    if ( Input_Opt.MET_LOADMET && Input_Opt.MET_LOADRH ) {
        simVars.relHumidity_w = Met.rhwRef();
        Input_Opt.MET_DEPTH = Met.satdepthUser();
        simVars.relHumidity_i = simVars.relHumidity_w * physFunc::pSat_H2Ol( simVars.temperature_K )\
                                      / physFunc::pSat_H2Os( simVars.temperature_K );
    }
    std::cout << "Temperature      = " << simVars.temperature_K << " K" << std::endl;
    std::cout << "Rel. humidity    = " << simVars.relHumidity_w << " %" << std::endl;
    std::cout << "Saturation depth = " << Input_Opt.MET_DEPTH << " m" << std::endl;

    /* ======================================================================= */
    /* ----------------------------------------------------------------------- */
    /* ------------------------ BACKGROUND CONDITIONS ------------------------ */
    /* ----------------------------------------------------------------------- */
    /* ======================================================================= */
    /* Declare solution structure */
    Solution Data(Input_Opt);

    /* Compute airDens from pressure and temperature */
    double airDens = simVars.pressure_Pa / ( physConst::kB   * simVars.temperature_K ) * 1.00E-06;
    /*     [molec/cm3] = [Pa = J/m3] / ([J/K]            * [K]           ) * [m3/cm3] */

    /* Set solution arrays to ambient data */
    Data.Initialize( simVars.BACKG_FILENAME.c_str(),      \
                     input, airDens, Met, \
                     Input_Opt, VAR, FIX, printDEBUG );


    /* Print Background Debug? */
    if ( DEBUG_BG_INPUT )
        Data.Debug( airDens );

    /* ======================================================================= */
    /* ----------------------------------------------------------------------- */
    /* --------------------------- TRANSPORT SOLVER -------------------------- */
    /* ----------------------------------------------------------------------- */
    /* ======================================================================= */

    /* Allocate horizontal and vertical diffusion parameters */
    double d_x, d_y;

    /* Allocate horizontal and vertical advection parameters */
    /* These correspond to domain-wide advection velocities (updraft, downdraft) */
    double vGlob_x, vGlob_y;

    /* Allocate horizontal and vertical distance traveled */
    double dTrav_x, dTrav_y;

    /* Allocate steady-state diffusion parameters */
    double D_X, D_Y;

    /* Allocate shear */
    double shear;
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
    shear = Met.shear();     /* [1/s] */
    if ( shear >= 0.0E+00 )
        LASTINDEX_SHEAR = Input_Opt.ADV_GRID_NX - 1;
    else
        LASTINDEX_SHEAR = 0;

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
    std::string engineInputFilePath = Input_Opt.SIMULATION_INPUT_ENG_EI;
    Aircraft aircraft = Aircraft(input, engineInputFilePath);

    /* Multiply by 500 since it gets multiplied by 1/500  
    * within the Emission object ... */ 
    JetA.setFSC( input.EI_SO2() * (double) 500.0 );
    
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

    double STEPMIN = (double)0.0;

    /* Allocate tolerances */
    double RTOL[NVAR];
    double ATOL[NVAR];

    /* Allocate doubles to store RH and IWC */
    double relHumidity, IWC;

    /* Initialize tolerances */
    for( UInt i = 0; i < NVAR; i++ ) {
        RTOL[i] = KPP_RTOLS;
        ATOL[i] = KPP_ATOLS;
    }
    
    /* Assign */
    UInt i_0 = std::floor( xE[0]/(xE[0]-xE[1]) ); //index i where x = 0
    UInt j_0 = std::floor( yE[0]/(yE[0]-yE[1]) ); //index j where y = 0
    Data.getData( VAR, FIX, i_0, j_0, simVars.CHEMISTRY );

    /* ======================================================================= */
    /* ----------------------------------------------------------------------- */
    /* -------------------------- EARLY MICROPHYSICS ------------------------- */
    /* ----------------------------------------------------------------------- */
    /* ======================================================================= */
    std::cout << "\n Starting EPM..." << std::endl;

    double Ice_rad, Ice_den, Soot_den, H2O_mol, SO4g_mol, SO4l_mol;
    double areaPlume;
    double Ab0 = input.bypassArea();
    double Tc0 = input.coreExitTemp();
    AIM::Aerosol liquidAer, iceAer;

    Vector_2D aerArray = Data.getAerosol();
    SimStatus EPM_RC = EPM::Integrate( simVars.temperature_K, simVars.pressure_Pa, simVars.relHumidity_w, VAR, \
                                 aerArray, aircraft, EI, Ice_rad, Ice_den, Soot_den,  \
                                 H2O_mol, SO4g_mol, SO4l_mol, liquidAer, iceAer, areaPlume, \
        		             Ab0, Tc0, simVars.CHEMISTRY, Input_Opt.ADV_AMBIENT_LAPSERATE, input.fileName_micro() );

    if((!simVars.CHEMISTRY) && (EPM_RC != SimStatus::EPMSuccess)) {
        return EPM_RC;
    }

    /* Compute initial plume area and scale initial ice aerosol properties based on number engines.
     * Note that EPM results are for ONLY ONE ENGINE.
     * If 2 engines, we assume that after 3 mins, the two plumes haven't fully mixed yet and result in a total
     * area of 2 * the area computed for one engine
     * If 3 or more engines, we assume that the plumes originating from the same wing have mixed. */

    areaPlume *= 2.0;
    if ( aircraft.EngNumber() != 2 ) {
        Ice_den  *= aircraft.EngNumber() / 2.0; //Scale densities by this factor to account for the plume size already doubling.
        Soot_den *= aircraft.EngNumber() / 2.0;
        liquidAer.scalePdf( aircraft.EngNumber()); //EPM only runs for one engine, so scale aerosol pdfs by engine number.
        iceAer.scalePdf( aircraft.EngNumber());
    }

    if ( iceAer.Moment() != 0 ) {
        /* Apply ice particle vortex losses using parameterization from
         * large-eddy simulations */
        const double iceNumFrac = aircraft.VortexLosses( EI.getSoot(),    \
                                                             EI.getSootRad(), \
                                                             Input_Opt.MET_DEPTH );

        if ( iceNumFrac <= 0.00E+00 && !simVars.CHEMISTRY ) {
            std::cout << "EndSim: vortex sinking" << std::endl;
            //exit(0);
            return SimStatus::NoSurvivalVortex;
        }
        iceAer.scalePdf( iceNumFrac );
    }

    const double semiYaxis = 0.5 * aircraft.deltaz1();
    const double semiXaxis = areaPlume / ( physConst::PI * semiYaxis );

    //Initializes/Calculates TRANSPORT_PA, TRANSPORT_LA
    simVars.initMicrophysicsVars(Data, liquidAer, Ice_den);

    Vector_1D vFall( Data.nBin_PA, 0.0E+00 );
    if ( simVars.TRANSPORT_PA() && simVars.GRAVSETTLING ) {
        /* Compute settling velocities */
        vFall = AIM::SettlingVelocity( Data.solidAerosol.getBinCenters(), \
                                       simVars.temperature_K, simVars.pressure_Pa );
    }

    /* Initialization at the grid scale level */

    /* Use the most inner ring to initialize at the grid-scale level */

    //Initializes Cluster with NRING rings and sets "ringCluster.halfRing() / semiring to false."
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
                      simVars.temperature_K, ( simVars.relHumidity_i > 100.0 ), \
                      liquidAer, iceAer, Soot_den, Met, areaPlume );


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
        std::cout << " ## - Temperature: " << std::setw(txtWidth) << simVars.temperature_K          << " [    K]\n";
        std::cout << " ## - Pressure   : " << std::setw(txtWidth) << simVars.pressure_Pa * 1.00E-02 << " [  hPa]\n";
        std::cout << " ## - Rel. Hum. I: " << std::setw(txtWidth) << simVars.relHumidity_i          << " [    %]\n";
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
    std::cout << "EPM Ice Particles: " << iceAer.Moment(0) * 1e6 * areaPlume << std::endl;
    /* Timeseries diagnostics */
    //This segfaults if TS_SPEC is true and chemistry is not enabled. -Michael
    if ( simVars.TS_SPEC ) {
        std::cout << "Saving chemistry" << std::endl;
        int hh = (int) (timestepVars.curr_Time_s - timestepVars.timeArray[0])/3600;
        int mm = (int) (timestepVars.curr_Time_s - timestepVars.timeArray[0])/60   - 60 * hh;
        int ss = (int) (timestepVars.curr_Time_s - timestepVars.timeArray[0])      - 60 * ( mm + 60 * hh );
        Diag::Diag_TS_Chem( simVars.TS_SPEC_FILEPATH.c_str(), simVars.TS_SPEC_LIST, hh, mm, ss, \
                      Data, m );
    }

    if ( simVars.TS_AERO ) {
        std::cout << "Saving aerosol" << std::endl;
        int hh = (int) (timestepVars.curr_Time_s - timestepVars.timeArray[0])/3600;
        int mm = (int) (timestepVars.curr_Time_s - timestepVars.timeArray[0])/60   - 60 * hh;
        int ss = (int) (timestepVars.curr_Time_s - timestepVars.timeArray[0])      - 60 * ( mm + 60 * hh );
        Diag::Diag_TS_Phys( simVars.TS_AERO_FILEPATH.c_str(), hh, mm, ss, \
                      Data.solidAerosol, Data.Species[ind_H2O], 
                      m.x(), m.y(), m.xE(), m.yE(), Met);

        float totalIceParticles = Data.solidAerosol.TotalNumber_sum( cellAreas );
        float totalIceMass = Data.solidAerosol.TotalIceMass_sum( cellAreas );
	std::cout << "# particles: " << totalIceParticles << ", ice mass: " << totalIceMass << std::endl;
        if ( totalIceParticles <= 1.00E+6 && totalIceMass <= 1.00E-2 && !simVars.CHEMISTRY ) {
            std::cout << "EndSim: no particles remain" << std::endl;
            return SimStatus::Complete;
        }
    }

    /* Prod & loss diagnostics */

    // TODO!!
    if ( simVars.SAVE_PL ) {

        /* If chemistry is performed at the grid cell level, then the
         * rates are stored as:
         * NY x NX x NFAM in [molec/cm^3/s]
         * into netCDF files at a frequency specified by the input file */

        int hh = (int) (timestepVars.curr_Time_s - timestepVars.timeArray[0])/3600;
        int mm = (int) (timestepVars.curr_Time_s - timestepVars.timeArray[0])/60   - 60 * hh;
        int ss = (int) (timestepVars.curr_Time_s - timestepVars.timeArray[0])      - 60 * ( mm + 60 * hh );
//        Diag_PL( "PL_hhmmss.nc", hh, mm, ss, Data, m );

    } else {
        if ( simVars.SAVE_O3PL ) {
        }
    }
        /* Fill with? */
    const double fillWith = 0.0E+00;

    /* Allocate Solvers */
    // SANDS::Solver Solver;
    // #pragma omp critical
    // {
    //     std::cout << "\n Initializing solver..." << std::endl;
    //     Solver.Initialize( /* Use threaded FFT?    */ simVars.THREADED_FFT, \
    //                        /* Use FFTW wisdom?     */ simVars.USE_WISDOM,   \
    //                        /* FFTW Directory       */ simVars.FFTW_DIR.c_str(),     \
    //                        /* Fill negative values */ simVars.FILLNEG,      \
    //                        /* Fill with this value */ fillWith );
    //     std::cout << "\n Initialization complete..." << std::endl;
    // }
    const FVM_ANDS::AdvDiffParams fvmSolverInitParams(0, 0, shear, D_X, D_Y, timestepVars.TRANSPORT_DT);
    const FVM_ANDS::BoundaryConditions H2O_BOUNDARY_COND = FVM_ANDS::bcFrom2DVector(Data.Species[ind_H2Oplume]);
    const FVM_ANDS::BoundaryConditions ZERO_BOUNDARY_COND = FVM_ANDS::bcFrom2DVector(Data.Species[ind_H2Oplume], true);

    std::vector<std::unique_ptr<FVM_ANDS::FVM_Solver>> fvmSolversVec;

    for(int i = 0; i < omp_get_max_threads(); i++) {
        fvmSolversVec.push_back(std::make_unique<FVM_ANDS::FVM_Solver>(fvmSolverInitParams, m.x(), m.y(), ZERO_BOUNDARY_COND, FVM_ANDS::std2dVec_to_eigenVec(Data.Species[ind_H2O]) ));
    }
    FVM_ANDS::FVM_Solver& fvmSolver = *fvmSolversVec[0];

    /* ======================================================================= */
    /* ----------------------------------------------------------------------- */
    /* ------------------------ TIME LOOP STARTS HERE ------------------------ */
    /* ----------------------------------------------------------------------- */
    /* ======================================================================= */
    
    bool early_stop;
    early_stop = false;
    while ( (timestepVars.curr_Time_s < timestepVars.tFinal_s) && (!early_stop) ) {
        /* Print message */
        std::cout << "\n";
        std::cout << "\n - Time step: " << timestepVars.nTime + 1 << " out of " << timestepVars.timeArray.size();
        #ifdef OMP
            std::cout << " (thread: " << omp_get_thread_num() << ")";
        #endif /* OMP */
        std::cout << "\n -> Solar time: " << std::fmod( timestepVars.curr_Time_s/3600.0, 24.0 ) << " [hr]" << std::endl;
        
        /* ======================================================================= */
        /* --------------------------- UPDATE TIMESTEP --------------------------- */
        /* ======================================================================= */

        /* Compute time step */
        timestepVars.checkLastStep();
	
        //PRE-SANDS H2O MASS CHECK
        std::cout << "Pre-transport total H2O Mass: " << totalH2OMass(Data, cellAreas) << std::endl;


        /* ======================================================================= */
        /* --------------------------- RUN TRANSPORT ----------------------------- */
        /* ======================================================================= */

        if (simVars.TRANSPORT && (timestepVars.nTime == 0 || timestepVars.checkTimeForTransport())) {

            std::cout << "Running transport..." << std::endl;
            timestepVars.lastTimeTransport = timestepVars.curr_Time_s + timestepVars.dt;
            //Update Transport Parameters. Computing at current time + dt_transport / 2
            PlumeModelUtils::DiffParam( timestepVars.curr_Time_s - timestepVars.tInitial_s + timestepVars.TRANSPORT_DT / 2.0, d_x, d_y, D_X, D_Y );
            /* Compute advection parameters */
            /* Is plume updraft on? */
            if ( simVars.UPDRAFT ) {

                /* Compute global advection velocities at mid time step */

                /* vGlob_x > 0 means left, < 0 means right [m/s]
                * vGlob_y > 0 means upwards, < 0 means downwards [m/s]
                * dTrav_x: distance traveled on the x-axis through advection [m]
                * dTrav_y: distance traveled on the y-axis through advection [m]
                */

                PlumeModelUtils::AdvGlobal( timestepVars.curr_Time_s - timestepVars.tInitial_s + timestepVars.TRANSPORT_DT / 2.0, simVars.UPDRAFT_TIME, simVars.UPDRAFT_VEL, \
                        vGlob_x, vGlob_y, dTrav_x, dTrav_y );
            }

            for(auto& solver: fvmSolversVec) {
                solver->updateTimestep(timestepVars.TRANSPORT_DT);
                solver->updateDiffusion(d_x, d_y);
                solver->updateAdvection(0, 0, shear);
            }

            /* TODO: If anyone ever decides to re-enable chemistry, please paralleize these transport solves at the NVAR level!
                     See the code in the if (simVars.TRANSPORT_PA) section for reference.
                     This is much more efficient than going through each bin and parallel-solving each bin
                     because the diffusion cannot be parallelized otherwise since it's Gauss-Seidel, and this 
                     also eliminates all the headaches and overhead of things like thread communication, false sharing, etc
                     when looping over a very long array.
                     -Michael*/
            if ( simVars.CHEMISTRY ) {
                /* Advection and diffusion of gas phase species */
                //Figure out what BCs to use with this later once chemistry is re-enabled. -Michael
                for ( N = 0; N < NVAR; N++ ) {
                    if ( N == ind_H2O )
                        fvmSolver.operatorSplitSolve2DVec( Data.Species[ind_H2O], H2O_BOUNDARY_COND);
                    else
                        /* Advection and diffusion of condensable species */
                        fvmSolver.operatorSplitSolve2DVec( Data.Species[N], ZERO_BOUNDARY_COND);
                }
            } 
            else {
                /* Advection and diffusion of H2O */
                for ( jNy = 0; jNy < Input_Opt.ADV_GRID_NY; jNy++ ) {
                    for ( iNx = 0; iNx < Input_Opt.ADV_GRID_NX; iNx++ ) {
                        Data.Species[ind_H2Oplume][jNy][iNx] = Data.Species[ind_H2O][jNy][iNx] - Data.Species[ind_H2Omet][jNy][iNx];
 
                    }
                }
                fvmSolver.operatorSplitSolve2DVec(Data.Species[ind_H2Oplume], H2O_BOUNDARY_COND, true);
            }
            /* Update H2O */
            for ( jNy = 0; jNy < Input_Opt.ADV_GRID_NY; jNy++ ) {
                for ( iNx = 0; iNx < Input_Opt.ADV_GRID_NX; iNx++ ) {
                    Data.Species[ind_H2O][jNy][iNx] = Data.Species[ind_H2Omet][jNy][iNx] + Data.Species[ind_H2Oplume][jNy][iNx];
                }
            }
            // /* Advection and diffusion for aerosol particles */
            // fvmSolver.operatorSplitSolve2DVec(Data.sootDens, ZERO_BOUNDARY_COND, true);

            // /* Monodisperse assumption for soot particles */
            // fvmSolver.operatorSplitSolve2DVec(Data.sootRadi, ZERO_BOUNDARY_COND, true);
            // fvmSolver.operatorSplitSolve2DVec(Data.sootArea, ZERO_BOUNDARY_COND, true);

            /* We assume that sulfate aerosols do not settle */
            /* TODO: Make this loop parallelize on the iBin_LA level, see if TRANSPORT_PA() loop below for reference.
                     Currently commenting out because this has no impact on contrail ice crystal behavior/impact
                     - Michael */
            // if ( simVars.TRANSPORT_LA() ) {
            //     /* Transport of liquid aerosols */
            //     for ( UInt iBin_LA = 0; iBin_LA < Data.nBin_LA; iBin_LA++ ){
            //         fvmSolver.operatorSplitSolve2DVec(Data.liquidAerosol.pdf[iBin_LA], ZERO_BOUNDARY_COND);
            //     }
            // }

            if ( simVars.TRANSPORT_PA() ) {
                /* Transport of solid aerosols */

                #pragma omp parallel for if(!PARALLEL_CASES) default(shared)
                for ( int iBin_PA = 0; iBin_PA < Data.nBin_PA; iBin_PA++ ) {
                    /* Transport particle number and volume for each bin and
                     * recompute centers of each bin for each grid cell
                     * accordingly */
                    int threadID = omp_get_thread_num();
                    FVM_ANDS::FVM_Solver& solver = *fvmSolversVec[threadID];
                    solver.updateAdvection(0, -vFall[iBin_PA], shear);

                    //passing in "false" to the "parallelAdvection" param to not spawn more threads
                    solver.operatorSplitSolve2DVec(Data.solidAerosol.getPDF_nonConstRef()[iBin_PA], ZERO_BOUNDARY_COND, false);

                }

                /* Check how much particle number and mass change before/after flux correction */
                timestepVars.totalIceParticles_before = Data.solidAerosol.TotalNumber_sum( cellAreas );
                timestepVars.totalIceMass_before = Data.solidAerosol.TotalIceMass_sum( cellAreas );
                if ( timestepVars.nTime == 0 ) timestepVars.totalIceMass_initial = timestepVars.totalIceMass_before;
                if ( timestepVars.nTime == 0 ) timestepVars.totalIceParticles_initial = timestepVars.totalIceParticles_before;

                /* Update met fields at mid time step */
                /* Update H2O */
                for ( jNy = 0; jNy < Input_Opt.ADV_GRID_NY; jNy++ ) {
                    for ( iNx = 0; iNx < Input_Opt.ADV_GRID_NX; iNx++ ) {
                        Data.Species[ind_H2Omet][jNy][iNx] = Met.H2O(jNy,iNx);
                        Data.Species[ind_H2O][jNy][iNx] = Data.Species[ind_H2Omet][jNy][iNx] + Data.Species[ind_H2Oplume][jNy][iNx];
                    }
                }
                Met.Update( timestepVars.TRANSPORT_DT, ( timestepVars.curr_Time_s + timestepVars.TRANSPORT_DT / 2 ) / 3600.0, ( timestepVars.curr_Time_s + timestepVars.TRANSPORT_DT / 2 - timestepVars.timeArray[0] ) / 3600, dTrav_x, dTrav_y );
                shear = Met.shear(); //shear needs to be updated because using local variable here...
            }
        } //TRANSPORT_PA
        std::cout << "Post-transport total H2O Mass: " << totalH2OMass(Data, cellAreas) << std::endl;

        /* ======================================================================= */
        /* -------------------- Apply Temp. Perturbation ------------------------- */
        /* ======================================================================= */

        if (simVars.TEMP_PERTURB && (timestepVars.nTime == 0 || timestepVars.checkTimeForTempPerturb())){
            std::cout << "Running temp. perturb..." << std::endl;
            Met.updateTempPerturb();
            timestepVars.lastTimeTempPerturb = timestepVars.curr_Time_s + timestepVars.dt;
        }

        /* ======================================================================= */
        /* ----------- UPDATE SOLAR ZENITH ANGLE AND PHOTOLYSIS RATES ------------ */
        /* ======================================================================= */

        /* Compute the cosine of solar zenith angle midway through the integration step */
        sun->Update( timestepVars.curr_Time_s + timestepVars.dt/2 );

        /* Reset photolysis rates */
        #pragma omp parallel for \
        if ( !PARALLEL_CASES ) \
        default ( shared     ) \
        schedule( dynamic, 1 )
        for ( UInt iPhotol = 0; iPhotol < NPHOTOL; iPhotol++ )
            jRate[iPhotol] = 0.0E+00;

        /* If daytime, update photolysis rates */
        if ( sun->CSZA > 0.0E+00 )
            Update_JRates( jRate, sun->CSZA );

        /* ======================================================================= */
        /* ------------------------------- RUN KPP ------------------------------- */
        /* ------------------------ The Kinetics Pre-Processor ------------------- */
        /* ======================================================================= */

        /* Is chemistry turned on? */
        if ( simVars.CHEMISTRY && timestepVars.checkTimeForChem()) {
            std::cout << "Running chemistry..." << std::endl;
            timestepVars.lastTimeChem = timestepVars.curr_Time_s + timestepVars.dt;
            Vector_2D iceVolume_ = Data.solidAerosol.TotalVolume();

            #pragma omp parallel for             \
            if      ( !PARALLEL_CASES         ) \
            default ( shared                   ) \
            private ( iNx, jNy                 ) \
            private ( relHumidity, IWC         ) \
            schedule( dynamic, 1               )
            for ( jNy = 0; jNy < Input_Opt.ADV_GRID_NY; jNy++ ) {
                for ( iNx = 0; iNx < Input_Opt.ADV_GRID_NX; iNx++ ) {

                    double AerosolArea[NAERO];
                    double AerosolRadi[NAERO];

                    /* Convert data structure to KPP inputs (VAR and FIX) */
                    Data.getData( VAR, FIX, iNx, jNy, simVars.CHEMISTRY );

                    /* ================================================= */
                    /* =============== Chemical rates ================== */
                    /* ================================================= */

                    /* Update heterogeneous chemistry reaction rates */
                    if ( simVars.HETCHEM ) {

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
                                    Met.airMolecDens(jNy,iNx), relHumidity, \
                                    Data.STATE_PSC, VAR, AerosolArea,  \
                                    AerosolRadi, IWC, &(Data.KHETI_SLA[0]), Input_Opt.ADV_TROPOPAUSE_PRESSURE);
                    }

                    /* Zero-out reaction rate */
                    for ( UInt iReact = 0; iReact < NREACT; iReact++ )
                        RCONST[iReact] = 0.0E+00;

                    /* Update photolysis rates */
                    for ( UInt iPhotol = 0; iPhotol < NPHOTOL; iPhotol++ )
                        PHOTOL[iPhotol] = jRate[iPhotol];

                    /* Update reaction rates */
                    Update_RCONST( Met.temp(jNy,iNx), Met.press(jNy), \
                                    Met.airMolecDens(jNy,iNx), VAR[ind_H2O] );

                    /* ================================================= */
                    /* ============= Chemical integration ============== */
                    /* ================================================= */

                    IERR = INTEGRATE( VAR, FIX, timestepVars.curr_Time_s, timestepVars.curr_Time_s + timestepVars.dt, \
                                        ATOL, RTOL, STEPMIN );

                    if ( IERR < 0 ) {
                        /* Integration failed */

                        std::cout << "Integration failed";
                        #ifdef OMP
                            std::cout << " on " << omp_get_thread_num();
                        #endif /* OMP */
                        std::cout << " for grid cell = (" << jNy << ", " << iNx << ") at time t = " << timestepVars.curr_Time_s/3600.0 << " ( nTime = " << timestepVars.nTime << " )\n";

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
                    Data.applyData(VAR, iNx, jNy );

                }
            }

            double AerosolArea[NAERO];
            double AerosolRadi[NAERO];

            /* Retrieve chemistry results into VAR and FIX */
            Data.getData(VAR, FIX);

            /* ========================================================= */
            /* ==================== Chemical rates ===================== */
            /* ========================================================= */

            /* Update heterogeneous chemistry reaction rates */
            if ( simVars.HETCHEM ) {

                for ( UInt iSpec = 0; iSpec < NSPEC; iSpec++ ) {
                    HET[iSpec][0] = 0.0E+00;
                    HET[iSpec][1] = 0.0E+00;
                    HET[iSpec][2] = 0.0E+00;
                }

                relHumidity = VAR[ind_H2O] * \
                                physConst::kB * simVars.temperature_K * 1.00E+06 / \
                                physFunc::pSat_H2Ol( simVars.temperature_K );
                GC_SETHET( simVars.temperature_K, simVars.pressure_Pa, airDens, relHumidity, \
                            Data.STATE_PSC, VAR, AerosolArea, AerosolRadi, IWC, &(Data.KHETI_SLA[0]), Input_Opt.ADV_TROPOPAUSE_PRESSURE );
            }

            /* Zero-out reaction rate */
            for ( UInt iReact = 0; iReact < NREACT; iReact++ )
                RCONST[iReact] = 0.0E+00;

            /* Update photolysis rates */
            for ( UInt iPhotol = 0; iPhotol < NPHOTOL; iPhotol++ )
                PHOTOL[iPhotol] = jRate[iPhotol];

            /* Update reaction rates */
            Update_RCONST( simVars.temperature_K, simVars.pressure_Pa, airDens, VAR[ind_H2O] );

            /* ========================================================= */
            /* ================= Chemical integration ================== */
            /* ========================================================= */

            IERR = INTEGRATE( VAR, FIX, timestepVars.curr_Time_s, timestepVars.curr_Time_s + timestepVars.dt, \
                                ATOL, RTOL, STEPMIN );

            if ( IERR < 0 ) {
                /* Integration failed */

                std::cout << "Integration failed";
                #ifdef OMP
                    std::cout << " on " << omp_get_thread_num();
                #endif /* OMP */
                std::cout << " for ambient conditions at time t = " << timestepVars.curr_Time_s/3600.0 << " ( nTime = " << timestepVars.nTime << " )\n";
                // return KPP_FAIL;
                return SimStatus::Failed;
            }

        }

        /* ======================================================================= */
        /* ------------------------ AEROSOL MICROPHYSICS ------------------------- */
        /* ---------------------------- COAGULATION ------------------------------ */
        /* ======================================================================= */

        //PRE-COAGULATION H2O MASS CHECK
        //std::cout << "Pre-Coagulation total H2O Mass: " << totalH2OMass(Data, cellAreas) << std::endl;

        if ( printDEBUG ) printf("DEBUG: T%06d -> %12.2e particles\n", 90, Data.solidAerosol.TotalNumber_sum( cellAreas ));

        /* Liquid aerosol coagulation */
        if (simVars.LIQ_COAG && timestepVars.checkTimeForLiqCoag()) {
            std::cout << "Running liquid coagulation..." << std::endl;
            timestepVars.lastTimeLiqCoag = timestepVars.curr_Time_s + timestepVars.dt;
            /* If shear = 0, take advantage of the symmetry around the Y-axis */
            Data.liquidAerosol.Coagulate( timestepVars.COAG_DT, Data.LA_Kernel, simVars.LA_MICROPHYSICS(), ( shear == 0.0E+00 ) && ( Input_Opt.ADV_GRID_XLIM_LEFT == Input_Opt.ADV_GRID_XLIM_RIGHT ) );
        }

        /* Solid aerosol coagulation */
        if (simVars.ICE_COAG && timestepVars.checkTimeForIceCoag()) {
            std::cout << "Running solid coagulation..." << std::endl;
            timestepVars.lastTimeIceCoag = timestepVars.curr_Time_s + timestepVars.dt;
            /* If shear = 0, take advantage of the symmetry around the Y-axis */
            Data.solidAerosol.Coagulate ( timestepVars.COAG_DT, Data.PA_Kernel, simVars.PA_MICROPHYSICS(), ( shear == 0.0E+00 ) && ( Input_Opt.ADV_GRID_XLIM_LEFT == Input_Opt.ADV_GRID_XLIM_RIGHT ) );
        }

        /* ======================================================================= */
        /* ------------------------ AEROSOL MICROPHYSICS ------------------------- */
        /* ------------------------- ICE CRYSTAL GROWTH -------------------------- */
        /* ======================================================================= */

        //PRE-ICE GROWTH H2O MASS CHECK
        std::cout << "Pre-Ice Growth total H2O Mass: " << totalH2OMass(Data, cellAreas) << std::endl;

        /* Solid aerosol growth */
        if (simVars.ICE_GROWTH && timestepVars.checkTimeForIceGrowth()) {
            std::cout << "Running ice growth..." << std::endl;
            timestepVars.lastTimeIceGrowth = timestepVars.curr_Time_s + timestepVars.dt;
            /* If shear = 0, take advantage of the symmetry around the Y-axis */
            Data.solidAerosol.Grow( timestepVars.ICE_GROWTH_DT, Data.Species[ind_H2O], Met.Temp(), Met.Press(), simVars.PA_MICROPHYSICS(), ( shear == 0.0E+00 ) && (  Input_Opt.ADV_GRID_XLIM_LEFT == Input_Opt.ADV_GRID_XLIM_RIGHT ) );
            std::cout<<"Ice Mass: " << Data.solidAerosol.TotalIceMass_sum(cellAreas)<<std::endl;
        }

        // If we do not perform chemistry, let's abort the simulation if there's no ice mass left
        if ( !simVars.CHEMISTRY && timestepVars.totalIceParticles_initial > 0.0E+00 ) {
            timestepVars.totalIceParticles_after = Data.solidAerosol.TotalNumber_sum( cellAreas );
            if ( timestepVars.totalIceParticles_after / timestepVars.totalIceParticles_initial < timestepVars.ABORT_THRESHOLD ) {
                std::cout << "Only " << 100 * timestepVars.totalIceParticles_after / timestepVars.totalIceParticles_initial << " % of initial (post-vortex sinking) ice particles remaining. Let's abort!" << std::endl;
                early_stop = true;
            }
        }
        //FINAL H2O MASS CHECK
        std::cout << "Final total H2O Mass: " << totalH2OMass(Data, cellAreas) << std::endl;

        timestepVars.curr_Time_s += timestepVars.dt;
        timestepVars.nTime++;

        /* Timeseries diagnostics */
        if ( simVars.TS_SPEC && \
           (( simVars.TS_FREQ == 0 ) || \
            ( std::fmod((timestepVars.curr_Time_s - timestepVars.timeArray[0])/60.0, simVars.TS_FREQ) == 0.0E+00 )) ) {
           int hh = (int) (timestepVars.curr_Time_s - timestepVars.timeArray[0])/3600;
           int mm = (int) (timestepVars.curr_Time_s - timestepVars.timeArray[0])/60   - 60 * hh;
           int ss = (int) (timestepVars.curr_Time_s - timestepVars.timeArray[0])      - 60 * ( mm + 60 * hh );
           Diag::Diag_TS_Chem( simVars.TS_SPEC_FILEPATH.c_str(), simVars.TS_SPEC_LIST, hh, mm, ss, \
                         Data, m );
        }

        if ( simVars.TS_AERO && \
           (( simVars.TS_AERO_FREQ == 0 ) || \
            ( std::fmod((timestepVars.curr_Time_s - timestepVars.timeArray[0])/60.0, simVars.TS_AERO_FREQ) == 0.0E+00 )) ) {
            int hh = (int) (timestepVars.curr_Time_s - timestepVars.timeArray[0])/3600;
            int mm = (int) (timestepVars.curr_Time_s - timestepVars.timeArray[0])/60   - 60 * hh;
            int ss = (int) (timestepVars.curr_Time_s - timestepVars.timeArray[0])      - 60 * ( mm + 60 * hh );
	        std::cout << "part lost=" << timestepVars.totPart_lost << ", ice lost=" << timestepVars.totIce_lost << std::endl;

            Diag::Diag_TS_Phys( simVars.TS_AERO_FILEPATH.c_str(), hh, mm, ss, 
                          Data.solidAerosol, Data.Species[ind_H2O], 
                          m.x(), m.y(), m.xE(), m.yE(), Met );    
            float totalIceParticles = Data.solidAerosol.TotalNumber_sum( cellAreas );
            float totalIceMass = Data.solidAerosol.TotalIceMass_sum( cellAreas );
	        std::cout << "# particles: " << totalIceParticles << ", ice mass: " << totalIceMass << std::endl;
            if ( totalIceParticles <= 1.00E+6 && totalIceMass <= 1.00E-2 && !simVars.CHEMISTRY ) {
                std::cout << "EndSim: no particles remain" << std::endl;
                std::cout << "# ice particles: " << totalIceParticles << std::endl;
                std::cout << "Total ice mass [g]: " << totalIceMass << std::endl;
                //exit(0);
                early_stop = true;
            }
        }
        if(early_stop) {
            break;
        }
    }
 
    SimStatus status;
    if (early_stop){
        status = SimStatus::Complete;
    } else {
        status = SimStatus::Incomplete;
    }
    /* ===================================================================== */
    /* ------------------------ TIME LOOP ENDS HERE ------------------------ */
    /* ===================================================================== */

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
    std::cout << "APCEMM LAGRID Plume Model Run Finished! Run time: " << duration.count() << "ms" << std::endl;
    return status;

} /* End of PlumeModel */

/* End of PlumeModel.cpp */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
