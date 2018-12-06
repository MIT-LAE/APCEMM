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
#include "Core/Parameters.hpp"
#include "Core/Interface.hpp"
#include "Core/Input.hpp"
#include "Core/Monitor.hpp"
#include "SANDS/Solver.hpp"
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

#if ( RINGS )
    #include "Core/Cluster.hpp"
    #include "Core/Species.hpp"
#endif /* RINGS */

#if ( TIME_IT )
    #include "Core/Timer.hpp"
#endif /* TIME_IT */

#if ( SAVE_FORWARD || SAVE_ADJOINT || SAVE_PA_MICROPHYS || SAVE_LA_MICROPHYS )
    #include "Core/Save.hpp"
    int isSaved = 1;
    static int SAVE_FAIL   = -2;
#endif /* SAVE_FORWARD || SAVE_ADJOINT || SAVE_PA_MICROPHYS || SAVE_LA_MICROPHYS */

#if ( CHEMISTRY )
    double C[NSPEC];          /* Concentration of all species */
    double * VAR = &C[0];     /* Concentration of variable species (global) */
    double * FIX = &C[NVAR];  /* Concentration of fixed species (global) */

    double RCONST[NREACT];    /* Rate constants (global) */
    double PHOTOL[NPHOTOL];   /* Photolysis rates (global) */
    double HET[NSPEC][3];     /* Heterogeneous chemistry rates (global) */

    double TIME;              /* Current integration time (global) */
#endif /* CHEMISTRY */

#if ( ADJOINT )
    double SZA_CST[3];
#endif /* ADJOINT */


void DiffParam( double time, double &d_x, double &d_y );
void AdvGlobal( double time, double &v_x, double &v_y, double &dTrav_x, double &dTrav_y );
std::vector<double> BuildTime( double tStart, double tFinal, \
               double sunRise, double sunSet );
double UpdateTime( double time, double tStart, \
                   double sunRise, double sunSet );
void Transport( Solution& Data, SANDS::Solver& Solver );


int PlumeModel( const Input &input )
{

bool printDEBUG = 0;

#ifdef DEBUG

    std::cout << "\n DEBUG is turned ON!\n\n";
    printDEBUG = 1;

#endif /* DEBUG */
    
    RealDouble temperature_K = input.temperature_K();
    RealDouble pressure_Pa   = input.pressure_Pa();
    RealDouble relHumidity_w = input.relHumidity_w();

    const RealDouble longitude_deg = input.longitude_deg();
    const RealDouble latitude_deg  = input.latitude_deg();

    const UInt dayGMT = input.dayGMT();
    const RealDouble emissionTime = input.emissionTime();

    /* Grid indices */
    unsigned int iNx = 0;
    unsigned int jNy = 0;
    
    int IERR;
    bool LAST_STEP = 0;
    bool ITS_TIME_FOR_LIQ_COAGULATION = 0;
    double lastTimeLiqCoag, dtLiqCoag;
    bool ITS_TIME_FOR_ICE_COAGULATION = 0;
    double lastTimeIceCoag, dtIceCoag;

#if ( TIME_IT )

    Timer Stopwatch, Stopwatch_cumul;
    unsigned long SANDS_clock, KPP_clock;
    unsigned long SANDS_clock_cumul = 0;
    unsigned long KPP_clock_cumul = 0;
    unsigned long clock_cumul = 0;
    bool reset = 1;
    
#endif /* TIME_IT */

#if ( NOy_MASS_CHECK )

    double mass_Ambient_NOy, mass_Emitted_NOy;

    #if ( RINGS )

    double mass_Emitted_NOy_Rings;

    #endif /* RINGS */

#endif /* NOy_MASS_CHECK */

#if ( CO2_MASS_CHECK )

    double mass_Ambient_CO2, mass_Emitted_CO2;

    #if ( RINGS )

    double mass_Emitted_CO2_Rings;

    #endif /* RINGS */

#endif /* CO2_MASS_CHECK */

    /* Compute relative humidity w.r.t ice */
    double relHumidity_i = relHumidity_w * physFunc::pSat_H2Ol( temperature_K )\
                                         / physFunc::pSat_H2Os( temperature_K );

    /* Define sun parameters */
    SZA *sun = new SZA( latitude_deg, dayGMT );    


    /** ~~~~~~~~~~~~~~~~~ **/
    /**        Mesh       **/
    /** ~~~~~~~~~~~~~~~~~ **/
    
    Mesh m;

    /* Get cell areas */
    const std::vector<std::vector<double>> cellAreas = m.areas();

    /** ~~~~~~~~~~~~~~~~~ **/
    /**        Time       **/
    /** ~~~~~~~~~~~~~~~~~ **/

    /*  
     *  - tEmission is the local emission time expressed in hours 
     *  (between 0.0 and 24.0)
     *  - tInitial is the local time at which the simulation starts in hours
     *  - TSIMUL represents the simulation time (in hours)
     *  - tFinal corresponds to the final time of the simulation expressed in hours
     */ 

    /* Define emission and simulation time */
    const double tEmission_h = emissionTime;        /* [hr] */
    const double tInitial_h  = tEmission_h;         /* [hr] */
    const double tFinal_h    = tInitial_h + TSIMUL; /* [hr] */
    const double tInitial_s  = tInitial_h * 3600.0; /* [s] */
    const double tFinal_s    = tFinal_h   * 3600.0; /* [s] */

    /* Current time in [s] */
    double curr_Time_s = tInitial_s; /* [s] */
    /* Time step in [s] */
    double dt = 0;                   /* [s] */

    /* Create time array */

    /* Vector of time in [s] */
    const std::vector<double> timeArray = BuildTime ( tInitial_s, tFinal_s, 3600.0*sun->sunRise, 3600.0*sun->sunSet );

    /* Time counter [-] */
    unsigned int nTime = 0;
    

    /** ~~~~~~~~~~~~~~~~~~ **/
    /**     Meteorology    **/
    /** ~~~~~~~~~~~~~~~~~~ **/

    Meteorology Met( LOAD_MET, m, temperature_K, 11.2E+03, -3.0E-03, printDEBUG );
    
    /** ~~~~~~~~~~~~~~~~~~ **/
    /**     Background     **/
    /** ~~~~~~~~~~~~~~~~~~ **/

    /* Declare solution structure */
    Solution Data;

    /* Compute airDens from pressure and temperature */
    double airDens = pressure_Pa / ( physConst::kB   * temperature_K ) * 1.00E-06;
    /* [molec/cm3] = [Pa = J/m3] / ([J/K]            * [K]           ) * [m3/cm3] */

    /* Set solution arrays to ambient data */
    Data.Initialize( AMBFILE, input, airDens, Met, printDEBUG );

    /* Print Background Debug? */
    if ( DEBUG_BG_INPUT )
        Data.Debug( airDens );

    /* Create ambient struture */
    Ambient ambientData( timeArray.size(), Data.getAmbient(), Data.getAerosol(), Data.getLiqSpecies() );

#if ( SAVE_LA_MICROPHYS )
    
    bool ITS_TIME_TO_SAVE_LA_OUTPUT = 0;
    std::vector<double> saveTime_LA;
    std::vector<std::vector<std::vector<std::vector<double>>>> saveOutput_LA( 1, std::vector<std::vector<std::vector<double>>>( Data.nBin_LA, std::vector<std::vector<double>>( NY, std::vector<double>( NX, 0.0E+00 ))));

#endif /* SAVE_LA_MICROPHYS */

#if ( SAVE_PA_MICROPHYS )

    bool ITS_TIME_TO_SAVE_PA_OUTPUT = 0;
    std::vector<double> saveTime_PA;
    std::vector<std::vector<std::vector<std::vector<double>>>> saveOutput_PA( 1, std::vector<std::vector<std::vector<double>>>( Data.nBin_PA, std::vector<std::vector<double>>( NY, std::vector<double>( NX, 0.0E+00 ))));

#endif /* SAVE_PA_MICROPHYS */

   
    /** ~~~~~~~~~~~~~~~~~ **/
    /**      Solver       **/
    /** ~~~~~~~~~~~~~~~~~ **/
    
    /* Allocate horizontal and vertical diffusion parameters */
    double d_x, d_y;

    /* Allocate horizontal and vertical advection parameters */
    double vGlob_x, vGlob_y; /* These correspond to domain-wide advection velocities (updraft, downdraft) */

    /* Allocate horizontal and vertical distance traveled */
    double dTrav_x, dTrav_y;

    /* Fill negative values? */
    const bool fillNegValues = 1;
    /* Fill with? */
    const double fillWith = 0.0E+00;

    /* Allocate Solvers */
    SANDS::Solver Solver;
    #pragma omp critical
    {
        Solver.Initialize( fillNegValues, fillWith );
    }
    

    /** ~~~~~~~~~~~~~~~~~ **/
    /**     Emissions     **/
    /** ~~~~~~~~~~~~~~~~~ **/

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
    Aircraft aircraft( aircraftName, temperature_K, pressure_Pa, relHumidity_w );

#if ( APCEMM_LUT )

    aircraft.setEI_NOx( input.EI_NOx() );
    aircraft.setEI_CO( input.EI_CO() );
    aircraft.setEI_HC( input.EI_HC() );
    aircraft.setEI_Soot( input.EI_Soot() );
    aircraft.setSootRad( input.sootRad() );
    aircraft.setFuelFlow( input.fuelFlow() );
    JetA.setFSC( input.EI_SO2() * (double) 500.0 );

#endif /* APCEMM_LUT */

    /* Print AC Debug? */
    if ( DEBUG_AC_INPUT )
        aircraft.Debug();

    /* Aggregate emissions from engine and fuel characteristics */
    const Emission EI( aircraft.getEngine(), JetA );

    /* Print Emission Debug? */
    if ( DEBUG_EI_INPUT )
        EI.Debug(); 

    
    /** ~~~~~~~~~~~~~~~~~ **/
    /**     Chemistry     **/
    /** ~~~~~~~~~~~~~~~~~ **/

    /* Allocate arrays for KPP */

    double STEPMIN = (double)0.0;

    double RTOL[NVAR];
    double ATOL[NVAR];

    for( unsigned int i = 0; i < NVAR; i++ ) {
        RTOL[i] = KPP_RTOLS; 
        ATOL[i] = KPP_ATOLS; 
    }

    /* aerArray stores all the number concentrations of aerosols */
    double aerArray[N_AER][2];
    
    /* Ambient chemistry */
    ambientData.getData( VAR, FIX, aerArray, nTime );

    
    /** ~~~~~~~~~~~~~~~~~~~~~~~ **/
    /**    Early Microphysics   **/
    /** ~~~~~~~~~~~~~~~~~~~~~~~ **/

    double Ice_rad, Ice_den, Soot_den, H2O_mol, SO4g_mol, SO4l_mol;
    double areaPlume; 
    AIM::Aerosol liquidAer, iceAer;
    EPM::Integrate( temperature_K, pressure_Pa, relHumidity_w, VAR, FIX, aerArray, aircraft, EI, \
                    Ice_rad, Ice_den, Soot_den, H2O_mol, SO4g_mol, SO4l_mol, liquidAer, iceAer, areaPlume );

    /* Compute initial plume area.
     * If 2 engines, we assume that after 3 mins, the two plumes haven't fully mixed yet and result in a total
     * area of 2 * the area computed for one engine
     * If 3 or more engines, we assume that the plumes originating from the same wing have mixed. */

    areaPlume *= 2;
    if ( aircraft.getEngNumber() != 2 ) {
        Ice_den  *= aircraft.getEngNumber() / 2.0;
        liquidAer.scalePdf( aircraft.getEngNumber() / 2.0 );
        iceAer.scalePdf( aircraft.getEngNumber() / 2.0 );
        Soot_den *= aircraft.getEngNumber() / 2.0;
    }

    double semiYaxis = 0.5*aircraft.getVortexdeltaz1();
    double semiXaxis = areaPlume/(physConst::PI*0.5*aircraft.getVortexdeltaz1());
   

    /* Liquid aerosol considerations */
    unsigned int LA_MICROPHYSICS;

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
    unsigned int PA_MICROPHYSICS;

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

    std::vector<double> vFall( Data.nBin_PA, 0.0E+00 );
    if ( TRANSPORT_PA ) {
        /* Compute settling velocities */
        vFall = AIM::SettlingVelocity( Data.solidAerosol.getBinCenters(), \
                                       temperature_K, pressure_Pa );
    }
    

    /** ~~~~~~~~~~~~~~~~~ **/
    /**      Rings?       **/
    /** ~~~~~~~~~~~~~~~~~ **/

#if ( RINGS )
        
    /** ~~~~~~~~~~~~~~~~~~ **/
    /**  Cluster of rings  **/
    /** ~~~~~~~~~~~~~~~~~~ **/

    /* Ring index */
    unsigned int iRing = 0;

    /* Create cluster of rings */
    Cluster ringCluster( NRING, ( relHumidity_i > 100.0 ), semiXaxis, semiYaxis, 0.0, 0.0 );
    
    /* Number of rings */
    const unsigned int nRing = ringCluster.getnRing();

    /* Print Ring Debug? */
    if ( DEBUG_RINGS )
        ringCluster.Debug();

    /* Allocate species-ring vector */
    SpeciesArray ringSpecies( nRing, timeArray.size(), ringCluster.halfRing() );

    /* Compute Grid to Ring mapping */        
    m.Ring2Mesh( ringCluster );
   
    /* Get mapping */
    const std::vector<std::vector<std::pair<unsigned int, unsigned int>>> mapRing2Mesh = m.list();
    
    /* Print ring to mesh mapping? */
    if ( DEBUG_MAPPING )
        m.Debug();
    
    /* Compute ring areas */
    ringCluster.ComputeRingAreas( cellAreas, mapRing2Mesh );
    const std::vector<double> ringArea = ringCluster.getRingArea();
    const double totArea = std::accumulate( ringArea.begin(), ringArea.end(), 0 );

    /* Add emission into the grid */
    Data.addEmission( EI, aircraft, mapRing2Mesh, cellAreas, ringCluster.halfRing(), temperature_K, ( relHumidity_i > 100.0 ), \
                      liquidAer, iceAer, Soot_den * areaPlume / ringArea[0] ); 
   
    /* Fill in variables species for initial time */
    ringSpecies.FillIn( Data, m, nTime );

    /* Allocate an additional array for KPP */
    double tempArray[NVAR];

    /* Allocate a ring's relative humidity */
    double relHumidity_Ring;
 
    /* Otherwise we do not have a ring structure and chemistry is solved on the grid */

#endif /* RINGS */

    /* Output run characteristics to log file/console */

#pragma omp critical
    {

        #ifdef OMP
            std::cout << "\n\n ## ON THREAD: " << omp_get_thread_num() << "\n ##";
        #endif /* OMP */
        const unsigned int coutPrecision = 5;
        const unsigned int txtWidth      = coutPrecision + 2;
        std::cout << std::setprecision(coutPrecision);
        std::cout << "\n ## ATMOSPHERIC COND.:";
        std::cout << "\n ##\n";
        std::cout << " ## - Temperature: " << std::setw(txtWidth) << temperature_K          << " [  K]\n";
        std::cout << " ## - Pressure   : " << std::setw(txtWidth) << pressure_Pa * 1.00E-02 << " [hPa]\n";
        std::cout << " ## - Rel. Hum. I: " << std::setw(txtWidth) << relHumidity_i          << " [  %]\n";
        std::cout << " ## - Latitude   : " << std::setw(txtWidth) << latitude_deg           << " [deg]\n";
        std::cout << " ## - Max CSZA   : " << std::setw(txtWidth) << sun->CSZA_max          << " [ - ]\n";

        std::cout << "\n ## EMISSIONS:";
        std::cout << "\n ##\n";
        std::cout << " ## - E_CO2 = " << std::setw(txtWidth+3) << EI.getCO2() * aircraft.getFuelFlow() / aircraft.getVFlight()           << " [kg(CO2)/km]"\
            " ( EI  = " << std::setw(txtWidth) << EI.getCO2() * 1.00E-03 << " [kg/kg_fuel] )\n";
        std::cout << " ## - E_CO  = " << std::setw(txtWidth+3) << EI.getCO()  * aircraft.getFuelFlow() / aircraft.getVFlight() * 1.0E+03 << " [ g(CO) /km]"\
            " ( EI  = " << std::setw(txtWidth) << EI.getCO()             << " [ g/kg_fuel] )\n";
        std::cout << " ## - E_CH4 = " << std::setw(txtWidth+3) << EI.getCH4() * aircraft.getFuelFlow() / aircraft.getVFlight() * 1.0E+06 << " [mg(CH4)/km]"\
            " ( EI  = " << std::setw(txtWidth) << EI.getCH4() * 1.00E+03 << " [mg/kg_fuel] )\n";
        std::cout << " ## - E_NOx = " << std::setw(txtWidth+3) << ( EI.getNO() / MW_NO + EI.getNO2() / MW_NO2 + EI.getHNO2() / MW_HNO2 ) * MW_N \
                                                  * aircraft.getFuelFlow() / aircraft.getVFlight() * 1.0E+03 << " [ g(N)  /km]"\
            " ( EI  = " << std::setw(txtWidth) << EI.getNOx()            << " [ g/kg_fuel] )\n";
        std::cout << " ## - E_SO2 = " << std::setw(txtWidth+3) << EI.getSO2() * aircraft.getFuelFlow() / aircraft.getVFlight() * 1.0E+03 << " [ g(SO2)/km]"\
            " ( EI  = " << std::setw(txtWidth) << EI.getSO2()            << " [ g/kg_fuel] )\n";
        std::cout << " ##                                   ( FSC = " << std::setw(txtWidth) << JetA.getFSC() << " [-]          )\n";
        std::cout << " ## - E_Soo = " << std::setw(txtWidth+3) << EI.getSoot() * aircraft.getFuelFlow() / aircraft.getVFlight() * 1.0E+03 << " [ g(Soo)/km]"\
            " ( EI  = " << std::setw(txtWidth) << EI.getSoot()* 1.00E+03 << " [mg/kg_fuel] )\n";
        std::cout << " ## - E_Soo = " << std::setw(txtWidth+3) << EI.getSoot() * aircraft.getFuelFlow() / aircraft.getVFlight() * 1.0E+03 / ( 4.0 / 3.0 * physConst::PI * physConst::RHO_SOOT * 1.00E+03 * EI.getSootRad() * EI.getSootRad() * EI.getSootRad() ) << " [ #(Soo)/km]"\
            " ( GMD = " << std::setw(txtWidth) << 2.0 * EI.getSootRad() * 1.0E+09 << " [nm]         )\n";

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

    }


    double frac_gSO4 = 0.0E+00;
    std::vector<double> aerosolProp( 3, 0.0E+00 );
    double AerosolArea[NAERO];
    double AerosolRadi[NAERO];
    double IWC = 0;
    double kheti_sla[11];

    lastTimeLiqCoag = curr_Time_s;
    lastTimeIceCoag = curr_Time_s;

#if ( SAVE_LA_MICROPHYS )
    
    saveOutput_LA[0] = Data.liquidAerosol.pdf;
    saveTime_LA.push_back( curr_Time_s );

#endif /* SAVE_LA_MICROPHYS */

#if ( SAVE_PA_MICROPHYS )

    saveOutput_PA[0] = Data.solidAerosol.pdf;
    saveTime_PA.push_back( curr_Time_s );

#endif /* SAVE_PA_MICROPHYS */



    /** ~~~~~~~~~~~~~~~~~~~~~~~~~~ **/
    /**         Time Loop          **/
    /** ~~~~~~~~~~~~~~~~~~~~~~~~~~ **/

#if ( TIME_IT )

    Stopwatch_cumul.Start( );

#endif /* TIME_IT */

    while ( curr_Time_s < tFinal_s ) {

        if ( printDEBUG ) {
            /* Print message */
            std::cout << "\n";
            std::cout << "\n - Time step: " << nTime + 1 << " out of " << timeArray.size();
            #ifdef OMP
                std::cout << " (thread: " << omp_get_thread_num() << ")";
            #endif /* OMP */
            std::cout << "\n -> Solar time: " << std::fmod( curr_Time_s/3600.0, 24.0 ) << " [hr]";
        }

        /** ~~~~~~~~~~~~~~~~~~~~~~~~~~ **/
        /**      Update Time Step      **/
        /** ~~~~~~~~~~~~~~~~~~~~~~~~~~ **/
        
        /* Compute time step */
        dt = timeArray[nTime+1] - timeArray[nTime];
        LAST_STEP = ( curr_Time_s + dt >= tFinal_s );

        Solver.UpdateTimeStep( dt );

        /** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ **/
        /**     Advection & Diffusion    **/
        /** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ **/
    
        /* Compute diffusion parameters */
        /* Is diffusion turned on? */
        if ( DIFFUSION ) {

            /* d_x: horizontal diffusion coefficient [m^2/s]
             * d_y: vertical diffusion coefficient [m^2/s]
             */

            DiffParam( curr_Time_s - tInitial_s, d_x, d_y );

        }
        else {

            /* If diffusion is turned off, set diffusion parameters to 0 */

            d_x = 0.0;
            d_y = 0.0;

        }

        /* Compute advection parameters */
        /* Is advection turned on? */
        if ( ADVECTION ) {

            /* Compute global advection velocities */

            /* vGlob_x > 0 means left, < 0 means right [m/s]
             * vGlob_y > 0 means upwards, < 0 means downwards [m/s]
             * dTrav_x: distance traveled on the x-axis through advection [m]
             * dTrav_y: distance traveled on the y-axis through advection [m]
             */
            
            AdvGlobal( curr_Time_s - tInitial_s, vGlob_x, vGlob_y, dTrav_x, dTrav_y ); 
            
        }
        else {

            /* If advection is turned off, set advection parameters to 0 */
            
            vGlob_x = 0;
            vGlob_y = 0;

        }
   
        /* Update diffusion and advection arrays */
        Solver.UpdateDiff( d_x, d_y );
        /* Assume no plume advection */
        Solver.UpdateAdv ( 0.0E+00, 0.0E+00 );
        /* Microphysics settling is considered for each bin independently */

        
        /** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ **/
        /** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ **/
        /** ~~~~~~~~~~~~~~~~~~~ SANDS ~~~~~~~~~~~~~~~~~~~~~ **/
        /** ~~~ Spectral Advection aNd Diffusion Solver ~~~ **/
        /** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ **/
        /** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ **/
        
#if ( TIME_IT )

        Stopwatch.Start( reset );

#endif /* TIME_IT */

#if ( DIFFUSION || ADVECTION )
        /* Advection and diffusion for gas phase species */
        Transport( Data, Solver ); //, plan_FFT, plan_IFFT );
        
        /* Advection and diffusion for aerosol particles */
        Solver.UpdateAdv ( 0.0E+00, 0.0E+00 );
        Solver.Run( Data.sootDens ); //, plan_FFT, plan_IFFT );
        /* Monodisperse assumption for soot particles */
        Solver.Run( Data.sootRadi ); //, plan_FFT, plan_IFFT );
        Solver.Run( Data.sootArea ); //, plan_FFT, plan_IFFT );
        
        /* We assume that sulfate aerosols do not settle */
        if ( TRANSPORT_LA ) {
            /* Transport of liquid aerosols */
            for ( unsigned int iBin_LA = 0; iBin_LA < Data.nBin_LA; iBin_LA++ )
                Solver.Run( Data.liquidAerosol.pdf[iBin_LA] ); //, plan_FFT, plan_IFFT );
        }

        if ( TRANSPORT_PA ) {
            /* Transport of solid aerosols */
            for ( unsigned int iBin_PA = 0; iBin_PA < Data.nBin_PA; iBin_PA++ ) {
                Solver.UpdateAdv ( 0.0E+00, vFall[iBin_PA] );
                Solver.Run( Data.solidAerosol.pdf[iBin_PA] ); //, plan_FFT, plan_IFFT );
            }
        }

#endif /* ( DIFFUSION || ADVECTION ) */
            
#if ( TIME_IT )
    
        Stopwatch.Stop( );
        SANDS_clock = Stopwatch.Elapsed( );

#endif /* TIME_IT */
        
        /* Update temperature field and pressure at new location */
        /*
         * To be implemented
         * Use dTrav_x and dTrav_y to update the temperature and pressure
         */

        /** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ **/
        /** ~~~~~~~~~~~~~~~ SO4 partitioning ~~~~~~~~~~~~~~ **/
        /** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ **/

        /* Compute SO4_l fraction */
        for ( iNx = 0; iNx < NX; iNx++ ) {
            for ( jNy = 0; jNy < NY; jNy++ ) {
                frac_gSO4 = H2SO4_GASFRAC( temperature_K, Data.SO4[jNy][iNx] );
                Data.SO4L[jNy][iNx] = ( 1.0 - frac_gSO4 ) * Data.SO4T[jNy][iNx];
                Data.SO4[jNy][iNx]  = Data.SO4T[jNy][iNx] - Data.SO4L[jNy][iNx];
            }
        }

        /** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ **/
        /** ~~~~ Update cosine of solar zenight angle ~~~~~ **/
        /** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ **/

        /* Compute the cosize of solar zenith angle midway through the integration step */
        sun->Update( curr_Time_s + dt/2 );

        /* Store cosine of solar zenith angle */
        ambientData.cosSZA[nTime+1] = sun->CSZA;

        if ( printDEBUG ) {
            std::cout << "\n DEBUG : \n";
            std::cout << "         CSZA = " << sun->CSZA << "\n";
        }


        /** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ **/
        /** ~~~~~~~~~~~ Update photolysis rates ~~~~~~~~~~~ **/
        /** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ **/

        for ( unsigned int iPhotol = 0; iPhotol < NPHOTOL; iPhotol++ )
            PHOTOL[iPhotol] = 0.0E+00;

        if ( sun->CSZA > 0.0E+00 )
            Read_JRates( PHOTOL, sun->CSZA );

        if ( printDEBUG ) {
            std::cout << "\n DEBUG : \n";
            for ( unsigned int iPhotol = 0; iPhotol < NPHOTOL; iPhotol++ )
                std::cout << "         PHOTOL[" << iPhotol << "] = " << PHOTOL[iPhotol] << "\n";
        }


        /** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ **/
        /** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ **/
        /** ~~~~~~~~~~~~~~~~~~~~ KPP ~~~~~~~~~~~~~~~~~~~~~~ **/
        /** ~~~~~~~~~ The Kinetics Pre-Processor ~~~~~~~~~~ **/
        /** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ **/
        /** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ **/
        

#if ( TIME_IT )
       
        Stopwatch.Start( reset );

#endif /* TIME_IT */
        
        /* Are we solving the chemistry in a ring structure? */
#if ( RINGS )
       
        /* Fill in variables species for current time */
        ringSpecies.FillIn( Data, m, nTime + 1 );

        /* Is chemistry turned on? */
        if ( CHEMISTRY ) {
        
            /* In-ring chemistry */
            for ( iRing = 0; iRing < nRing ; iRing++ ) {

                /* Convert ring structure to KPP inputs (VAR and FIX) */
                ringSpecies.getData( VAR, FIX, nTime + 1, iRing );

                for ( unsigned int iSpec = 0; iSpec < NVAR; iSpec++ )
                    tempArray[iSpec] = VAR[iSpec];

                /* ~~~~~~~~~~~~~~~~~~~~~~~~ */
                /* ~~~~ Chemical rates ~~~~ */
                /* ~~~~~~~~~~~~~~~~~~~~~~~~ */

                /* Update heterogeneous chemistry reaction rates */
                if ( HETCHEMISTRY ) {
                    
                    for ( unsigned int iSpec = 0; iSpec < NSPEC; iSpec++ ) {
                        HET[iSpec][0] = 0.0E+00;
                        HET[iSpec][1] = 0.0E+00;
                        HET[iSpec][2] = 0.0E+00;
                    }

                    Data.getAerosolProp( AerosolRadi, AerosolArea, IWC, mapRing2Mesh[iRing] );

                    relHumidity_Ring = VAR[ind_H2O] * \
                                       physConst::kB * temperature_K * 1.00E+06 / \
                                       physFunc::pSat_H2Ol( temperature_K );
                    GC_SETHET( temperature_K, pressure_Pa, airDens, relHumidity_Ring, \
                               Data.STATE_PSC, VAR, AerosolArea, AerosolRadi, IWC, kheti_sla );

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

                /* Update reaction rates */
                for ( unsigned int iReact = 0; iReact < NREACT; iReact++ )
                    RCONST[iReact] = 0.0E+00;

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
                    std::cout << " for ring = " << iRing << " at time t = " << curr_Time_s/3600.0 << " ( nTime = " << nTime << " )\n";

                    if ( printDEBUG ) {
                        std::cout << " ~~~ Printing reaction rates:\n";
                        for ( unsigned int iReact = 0; iReact < NREACT; iReact++ ) {
                            std::cout << "Reaction " << iReact << ": " << RCONST[iReact] << " [molec/cm^3/s]\n";
                        }
                        std::cout << " ~~~ Printing concentrations:\n";
                        for ( unsigned int iSpec = 0; iSpec < NVAR; iSpec++ ) {
                            std::cout << "Species " << iSpec << ": " << VAR[iSpec]/airDens*1.0E+09 << " [ppb]\n";
                        }
                    }
                    
                    delete sun; sun = NULL;
                    return KPP_FAIL;
                }
                
                ringSpecies.FillIn( VAR, nTime + 1, iRing );
            
                Data.applyRing( VAR, tempArray, mapRing2Mesh, iRing );

            }

            /* Ambient chemistry */
            ambientData.getData( VAR, FIX, aerArray, nTime );

            /* ~~~~~~~~~~~~~~~~~~~~~~~~ */
            /* ~~~~ Chemical rates ~~~~ */
            /* ~~~~~~~~~~~~~~~~~~~~~~~~ */

            /* Update heterogeneous chemistry reaction rates */
            if ( HETCHEMISTRY ) {

                for ( unsigned int iSpec = 0; iSpec < NSPEC; iSpec++ ) {
                    HET[iSpec][0] = 0.0E+00;
                    HET[iSpec][1] = 0.0E+00;
                    HET[iSpec][2] = 0.0E+00;
                }

                relHumidity_Ring = VAR[ind_H2O] * \
                                   physConst::kB * temperature_K * 1.00E+06 / \
                                   physFunc::pSat_H2Ol( temperature_K );
                GC_SETHET( temperature_K, pressure_Pa, airDens, relHumidity_Ring, \
                           Data.STATE_PSC, VAR, AerosolArea, AerosolRadi, IWC, kheti_sla );

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

            /* Update reaction rates */
            for ( unsigned int iReact = 0; iReact < NREACT; iReact++ )
                RCONST[iReact] = 0.0E+00;

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
                    for ( unsigned int iReact = 0; iReact < NREACT; iReact++ ) {
                        std::cout << "Reaction " << iReact << ": " << RCONST[iReact] << " [molec/cm^3/s]\n";
                    }
                    std::cout << " ~~~ Printing concentrations:\n";
                    for ( unsigned int iSpec = 0; iSpec < NVAR; iSpec++ ) {
                        std::cout << "Species " << iSpec << ": " << VAR[iSpec]/airDens*1.0E+09 << " [ppb]\n";
                    }
                }

                delete sun; sun = NULL;
                return KPP_FAIL;
            }

            ambientData.FillIn( VAR, nTime + 1 );

            Data.applyAmbient( VAR, mapRing2Mesh, nRing );

        }


        /* Otherwise solve chemistry on the grid */
#else

        /* Is chemistry turned on? */
        if ( CHEMISTRY ) {

            for ( iNx = 0; iNx < NX; iNx++ ) {
                for ( jNy = 0; jNy < NY; jNy++ ) {

                    /* Convert data structure to KPP inputs (VAR and FIX) */
                    Data.getData( VAR, FIX, iNx, jNy );

                    /* ~~~~~~~~~~~~~~~~~~~~~~~~ */
                    /* ~~~~ Chemical rates ~~~~ */
                    /* ~~~~~~~~~~~~~~~~~~~~~~~~ */

                    /* Update heterogeneous chemistry reaction rates */
                    if ( HETCHEMISTRY ) {
                        
                        for ( unsigned int iSpec = 0; iSpec < NSPEC; iSpec++ ) {
                            HET[iSpec][0] = 0.0E+00;
                            HET[iSpec][1] = 0.0E+00;
                            HET[iSpec][2] = 0.0E+00;
                        }

                        relHumidity_Ring = VAR[ind_H2O] * \
                                           physConst::kB * met.temp[jNy][iNx] * 1.00E+06 / \
                                           physFunc::pSat_H2Ol( met.temp[jNy][iNx] );
                        GC_SETHET( met.temp[jNy][iNx], met.press[jNy], airDens, relHumidity_Ring, \
                                   Data.STATE_PSC, VAR, AerosolArea, AerosolRadi, IWC, kheti_sla );
                    }

                    /* Update reaction rates */
                    for ( unsigned int iReact = 0; iReact < NREACT; iReact++ )
                        RCONST[iReact] = 0.0E+00;

                    Update_RCONST( met.temp[jNy][iNx], met.press[jNy], airDens, VAR[ind_H2O] );

                    /* ~~~~~~~~~~~~~~~~~~~~~~~~ */
                    /* ~~~~~ Integration ~~~~~~ */
                    /* ~~~~~~~~~~~~~~~~~~~~~~~~ */

                    IERR = KPP_Main( VAR, FIX, curr_Time_s, dt, \
                                     KPP_RTOLS, KPP_ATOLS );
                    
                    if ( IERR < 0 ) {
                        /* Integration failed */

                        std::cout << "Integration failed";
                        #ifdef OMP
                            std::cout << " on " << omp_get_thread_num();
                        #endif /* OMP */
                        std::cout << " for grid cell = (" << jNy << ", " << iNx << ") at time t = " << curr_Time_s/3600.0 << " ( nTime = " << nTime << " )\n";
                    
                        if ( printDEBUG ) {
                            std::cout << " ~~~ Printing reaction rates:\n";
                            for ( unsigned int iReact = 0; iReact < NREACT; iReact++ ) {
                                std::cout << "Reaction " << iReact << ": " << RCONST[iReact] << " [molec/cm^3/s]\n";
                            }
                            std::cout << " ~~~ Printing concentrations:\n";
                            for ( unsigned int iSpec = 0; iSpec < NVAR; iSpec++ ) {
                                std::cout << "Species " << iSpec << ": " << VAR[iSpec]/airDens*1.0E+09 << " [ppb]\n";
                            }
                        }

                        delete sun; sun = NULL;
                        return KPP_FAIL;
                    }

                    /* Convert KPP output back to data structure */
                    Data.applyData( VAR, iNx, jNy );
                }
            }
            
            /* Ambient chemistry */
            ambientData.getData( VAR, FIX, aerArray, nTime );
            
            /* ~~~~~~~~~~~~~~~~~~~~~~~~ */
            /* ~~~~ Chemical rates ~~~~ */
            /* ~~~~~~~~~~~~~~~~~~~~~~~~ */

            /* Update heterogeneous chemistry reaction rates */
            if ( HETCHEMISTRY ) {

                for ( unsigned int iSpec = 0; iSpec < NSPEC; iSpec++ ) {
                    HET[iSpec][0] = 0.0E+00;
                    HET[iSpec][1] = 0.0E+00;
                    HET[iSpec][2] = 0.0E+00;
                }

                relHumidity_Ring = VAR[ind_H2O] * \
                                   physConst::kB * met.temp[jNy][iNx] * 1.00E+06 / \
                                   physFunc::pSat_H2Ol( met.temp[jNy][iNx] );
                GC_SETHET( met.temp[jNy][iNx], met.press[jNy], airDens, relHumidity_Ring, \
                           Data.STATE_PSC, VAR, AerosolArea, AerosolRadi, IWC, kheti_sla );
            }


            /* Update reaction rates */
            for ( unsigned int iReact = 0; iReact < NREACT; iReact++ )
                RCONST[iReact] = 0.0E+00;
            
            Update_RCONST( met.temp[jNy][iNx], met.press[jNy], airDens, VAR[ind_H2O], sun->CSZA );

            /* ~~~~~~~~~~~~~~~~~~~~~~~~ */
            /* ~~~~~~ Integration ~~~~~ */
            /* ~~~~~~~~~~~~~~~~~~~~~~~~ */

            IERR = KPP_Main( VAR, FIX, curr_Time_s, dt, \
                             KPP_RTOLS, KPP_ATOLS );

            if ( IERR < 0 ) {
                /* Integration failed */

                std::cout << "Integration failed";
                #ifdef OMP
                    std::cout << " on " << omp_get_thread_num();
                #endif /* OMP */
                std::cout << " for ambient conditions at time t = " << curr_Time_s/3600.0 << " ( nTime = " << nTime << " )\n";

                if ( printDEBUG ) {
                    std::cout << " ~~~ Printing reaction rates:\n";
                    for ( unsigned int iReact = 0; iReact < NREACT; iReact++ ) {
                        std::cout << "Reaction " << iReact << ": " << RCONST[iReact] << " [molec/cm^3/s]\n";
                    }
                    std::cout << " ~~~ Printing concentrations:\n";
                    for ( unsigned int iSpec = 0; iSpec < NVAR; iSpec++ ) {
                        std::cout << "Species " << iSpec << ": " << VAR[iSpec]/airDens*1.0E+09 << " [ppb]\n";
                    }
                }

                delete sun; sun = NULL;
                return KPP_FAIL;
            }

            ambientData.FillIn( VAR, nTime + 1 );

        }

#endif /* RINGS */


#if ( TIME_IT )
        
        Stopwatch.Stop( );
        KPP_clock = Stopwatch.Elapsed( );

#endif /* TIME_IT */


        /** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ **/
        /** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ **/
        /** ~~~~~~~~~~~~~ Aerosol Microphysics ~~~~~~~~~~~~ **/
        /** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ **/
        /** ~~~~~~~~~~~~~~~~~ Coagulation ~~~~~~~~~~~~~~~~~ **/
        /** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ **/
        /** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ **/
    
        ITS_TIME_FOR_LIQ_COAGULATION = ( ( ( curr_Time_s - lastTimeLiqCoag ) >= LIQCOAG_TSTEP ) || LAST_STEP );
        /* Liquid aerosol coagulation */
        if ( ITS_TIME_FOR_LIQ_COAGULATION && LIQ_MICROPHYSICS ) {
            dtLiqCoag = ( curr_Time_s - lastTimeLiqCoag );
            if ( printDEBUG )
                std::cout << "\n DEBUG (Liquid Coagulation): Current time: " << ( curr_Time_s - tInitial_s ) / 3600.0 << " hr. Last coagulation event was at: " << ( lastTimeLiqCoag - tInitial_s ) / 3600.0 << " hr. Running for " << dtLiqCoag << " s\n";

            lastTimeLiqCoag = curr_Time_s;
            /* Here we assume that the sulfate aerosol fields are symmetric around the X and Y axis */
            Data.liquidAerosol.Coagulate( dtLiqCoag, Data.LA_Kernel, LA_MICROPHYSICS, (unsigned int) 2 );
        }

        ITS_TIME_FOR_ICE_COAGULATION = ( ( ( curr_Time_s - lastTimeIceCoag ) >= ICECOAG_TSTEP ) || LAST_STEP );
        /* Solid aerosol coagulation */
        if ( ITS_TIME_FOR_ICE_COAGULATION && ICE_MICROPHYSICS ) {
            dtIceCoag = ( curr_Time_s - lastTimeIceCoag );
            if ( printDEBUG )
                std::cout << "\n DEBUG (Solid Coagulation): Current time: " << ( curr_Time_s - tInitial_s ) / 3600.0 << " hr. Last coagulation event was at: " << ( lastTimeIceCoag - tInitial_s ) / 3600.0 << " hr. Running for " << dtIceCoag << " s\n";

            lastTimeIceCoag = curr_Time_s;
            /* Here we assume that the solid aerosol fields are symmetric around the X axis */
            Data.solidAerosol.Coagulate ( dtIceCoag, Data.PA_Kernel, PA_MICROPHYSICS, (unsigned int) ( 2 - ( relHumidity_i > 100.0 ) ) );
        }
        
        /** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ **/
        /** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ **/
        /** ~~~~~~~~~~~~~~~~~~~ Growth ~~~~~~~~~~~~~~~~~~~~ **/
        /** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ **/
        /** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ **/


#if ( SAVE_LA_MICROPHYS )
    
        ITS_TIME_TO_SAVE_LA_OUTPUT = ( ( ( curr_Time_s - saveTime_LA.back() ) >= SAVE_LA_DT ) || LAST_STEP );
        /* Save liquid aerosol at current time */
        if ( ITS_TIME_TO_SAVE_LA_OUTPUT ) {
            if ( printDEBUG )
                std::cout << "\n DEBUG (Save Liquid Aerosols): Current time: " << ( curr_Time_s - tInitial_s ) / 3600.0 << " hr. Last time liquid aerosols were saved: " << ( saveTime_LA.back() - tInitial_s ) / 3600.0 << " hr\n";

            if ( LAST_STEP )
                saveTime_LA.push_back( curr_Time_s + dt );
            else
                saveTime_LA.push_back( curr_Time_s );
            saveOutput_LA.push_back( Data.liquidAerosol.pdf );
        }

#endif /* SAVE_LA_MICROPHYS */

#if ( SAVE_PA_MICROPHYS )
    
        ITS_TIME_TO_SAVE_PA_OUTPUT = ( ( ( curr_Time_s - saveTime_PA.back() ) >= SAVE_PA_DT ) || LAST_STEP );
        /* Save solid aerosol at current time */
        if ( ITS_TIME_TO_SAVE_PA_OUTPUT ) {
            if ( printDEBUG )
                std::cout << "\n DEBUG (Save Solid Aerosols): Current time: " << ( curr_Time_s - tInitial_s ) / 3600.0 << " hr. Last time solid aerosols were saved: " << ( saveTime_PA.back() - tInitial_s ) / 3600.0 << " hr\n";

            if ( LAST_STEP )
                saveTime_PA.push_back( curr_Time_s + dt );
            else
                saveTime_PA.push_back( curr_Time_s );
            saveOutput_PA.push_back( Data.solidAerosol.pdf );
        }

#endif /* SAVE_PA_MICROPHYS */


#if ( NOy_MASS_CHECK )

        /* Compute ambient concentrations */
        mass_Ambient_NOy = ambientData.NO[nTime+1] + ambientData.NO2[nTime+1] + ambientData.NO3[nTime+1] +  ambientData.HNO2[nTime+1] \
                         + ambientData.HNO3[nTime+1] + ambientData.HNO4[nTime+1] + 2*ambientData.N2O5[nTime+1] + ambientData.PAN[nTime+1] \
                         + ambientData.MPN[nTime+1] + ambientData.N[nTime+1] + ambientData.PROPNN[nTime+1] + ambientData.BrNO2[nTime+1] \
                         + ambientData.BrNO3[nTime+1] + ambientData.ClNO2[nTime+1] + ambientData.ClNO3[nTime+1] + ambientData.PPN[nTime+1] \
                         + ambientData.PRPN[nTime+1] + ambientData.R4N1[nTime+1] + ambientData.PRN1[nTime+1] + ambientData.R4N2[nTime+1] \
                         + 2*ambientData.N2O[nTime+1];

        /* Compute emitted */
        mass_Emitted_NOy = 0;
        for ( iNx = 0; iNx < NX; iNx++ ) {
            for ( jNy = 0; jNy < NY; jNy++ ) {
                mass_Emitted_NOy += ( Data.NO[jNy][iNx] + Data.NO2[jNy][iNx] + Data.NO3[jNy][iNx] + Data.HNO2[jNy][iNx] \
                                    + Data.HNO3[jNy][iNx] + Data.HNO4[jNy][iNx] + 2*Data.N2O5[jNy][iNx] + Data.PAN[jNy][iNx] \
                                    + Data.MPN[jNy][iNx] + Data.N[jNy][iNx] + Data.PROPNN[jNy][iNx] + Data.BrNO2[jNy][iNx] \
                                    + Data.BrNO3[jNy][iNx] + Data.ClNO2[jNy][iNx] + Data.ClNO3[jNy][iNx]  + Data.PPN[jNy][iNx] \
                                    + Data.PRPN[jNy][iNx] + Data.R4N1[jNy][iNx] + Data.PRN1[jNy][iNx] + Data.R4N2[jNy][iNx] \
                                    + 2*Data.N2O[jNy][iNx] \
                                    - mass_Ambient_NOy ) * cellAreas[jNy][iNx];
            }
        }

        /* Print to console */
        std::cout << "\n\n    " << " *** NOy mass check: ";
        std::cout << "\n    " << " ~~> Emitted NOy: " << std::setw(6) << mass_Emitted_NOy * 1.0E+06 / physConst::Na * MW_N * 1.0E+06 << " [g(N)/km] ";
        /*                                                               [molec/cm3 * m2] * [m3/cm3]/ [molec/mole]  * [kg/mole]*[g/kg*m/km] = [g/km] */
    
    #if ( RINGS ) 
        
        mass_Emitted_NOy_Rings = 0;
        for ( iRing = 0; iRing < nRing; iRing++ ) {
            mass_Emitted_NOy_Rings += ( ringSpecies.NO[nTime+1][iRing] + ringSpecies.NO2[nTime+1][iRing] + ringSpecies.NO3[nTime+1][iRing] \
                                      + ringSpecies.HNO2[nTime+1][iRing] + ringSpecies.HNO3[nTime+1][iRing] + ringSpecies.HNO4[nTime+1][iRing] \
                                      + 2*ringSpecies.N2O5[nTime+1][iRing] + ringSpecies.PAN[nTime+1][iRing] + ringSpecies.MPN[nTime+1][iRing] \
                                      + ringSpecies.N[nTime+1][iRing] + ringSpecies.PROPNN[nTime+1][iRing] + ringSpecies.BrNO2[nTime+1][iRing] \
                                      + ringSpecies.BrNO3[nTime+1][iRing] + ringSpecies.ClNO2[nTime+1][iRing] + ringSpecies.ClNO3[nTime+1][iRing] \
                                      + ringSpecies.PPN[nTime+1][iRing] + ringSpecies.PRPN[nTime+1][iRing] + ringSpecies.R4N1[nTime+1][iRing] \
                                      + ringSpecies.PRN1[nTime+1][iRing] + ringSpecies.R4N2[nTime+1][iRing] + 2*ringSpecies.N2O[nTime+1][iRing] \
                                      - mass_Ambient_NOy ) * ringArea[iRing]; 
        }
        /* How much of this emitted mass is still in the rings? FR = Fraction in rings */
        std::cout << "(FR: " << 100 * mass_Emitted_NOy_Rings / mass_Emitted_NOy << " %)";

    #endif /* RINGS */

#endif /* NOy_MASS_CHECK */

#if ( CO2_MASS_CHECK )

        /* CO2 is not an exactly conserved quantity because of the oxidation CO and other compounds (unless chemistry is turned off) */

        mass_Ambient_CO2 = ambientData.CO2[nTime+1];
        
        /* Compute emitted */
        mass_Emitted_CO2 = 0;
        for ( iNx = 0; iNx < NX; iNx++ ) {
            for ( jNy = 0; jNy < NY; jNy++ ) {
                mass_Emitted_CO2 += ( Data.CO2[jNy][iNx] - mass_Ambient_CO2 ) * cellAreas[jNy][iNx];
            }
        }

        std::cout << "\n\n    " << " *** CO2 mass check: ";

        std::cout << "\n    " << " ~~> Emitted CO2: " << std::setw(6) << mass_Emitted_CO2 * 1.0E+06 / physConst::Na * MW_CO2 * 1.0E+03 << " [kg/km]   ";
        /*                                                               [molec/cm3 * m2] * [m3/cm3]/ [molec/mole]  *[kg/mole]*[m/km] = [kg/km] */

    #if ( RINGS ) 
        
        mass_Emitted_CO2_Rings = 0;
        for ( iRing = 0; iRing < nRing; iRing++ ) {
            mass_Emitted_CO2_Rings += ( ringSpecies.CO2[nTime+1][iRing] - mass_Ambient_CO2 ) * ringArea[iRing]; 
        }
        /* How much of this emitted mass is still in the rings? FR = Fraction in rings */
        std::cout << "(FR: " << 100 * mass_Emitted_CO2_Rings / mass_Emitted_CO2 << " %)\n";

    #endif /* RINGS */

#endif /* CO2_MASS_CHECK */

#if ( TIME_IT )

        SANDS_clock_cumul += SANDS_clock;
        KPP_clock_cumul   += KPP_clock;
        std::cout << "\n    " << " *** Clock breakdown: ";
        std::cout << "\n    " << " *** ----------------- ";
        std::cout << "\n    " << " *** Total: " << SANDS_clock + KPP_clock << " [ms]";
        std::cout << " ( SANDS: " << 100 * ( SANDS_clock / double( SANDS_clock + KPP_clock ) ) << "% , KPP: " << 100 * ( KPP_clock / double( SANDS_clock + KPP_clock ) ) << "% )";

#endif /* TIME_IT */

        curr_Time_s += dt;
        nTime++;

    }

#if ( TIME_IT )

    Stopwatch_cumul.Stop( );
    clock_cumul = Stopwatch_cumul.Elapsed( );

    std::cout << "\n";
    std::cout << " ** Final clock breakdown: " << "\n";
    
    std::cout << " ** -> SANDS: ";
    std::cout << std::setw(6) <<SANDS_clock_cumul / double(1000) << " [s] (" << 100 * SANDS_clock_cumul / double(clock_cumul) << " %)" << "\n";
    
    std::cout << " ** -> KPP  : ";
    std::cout << std::setw(6) << KPP_clock_cumul / double(1000) << " [s] (" << 100 * KPP_clock_cumul / double(clock_cumul) << " %)" << "\n";
    
    std::cout << " ** -> Rem. : ";
    std::cout << std::setw(6) << ( clock_cumul - SANDS_clock_cumul - KPP_clock_cumul ) / double(1000) << " [s] (" << 100 * ( clock_cumul - SANDS_clock_cumul - KPP_clock_cumul ) / double(clock_cumul) << " %)" << "\n";
    
    std::cout << " ** ----------------- " << "\n";
    std::cout << " ** Total   : ";
    std::cout << std::setw(6) << clock_cumul / double(1000) << " [s]" << "\n";
    std::cout << "\n";

#endif /* TIME_IT */


    #if ( SAVE_FORWARD )

    { 
        isSaved = output::Write( input.fileName2char(),               \
                                 ringSpecies,                         \
                                 ambientData,                         \
                                 ringCluster,                         \
                                 timeArray,                           \
                                 input,                               \
                                 airDens, relHumidity_i,              \
                                 sun->sunRise, sun->sunSet );
        if ( isSaved == output::SAVE_FAILURE ) {
            std::cout << " Saving to ring-averaged concentrations to file failed...\n";
            return SAVE_FAIL;
        }
    }

    #endif /* SAVE_FORWARD */

    #if ( SAVE_LA_MICROPHYS )
    
    isSaved = output::Write_MicroPhys( OUT_FILE_LA, saveOutput_LA, saveTime_LA, \
                                       Data.liquidAerosol.getBinCenters(), \
                                       m.x(), m.y(), temperature_K, pressure_Pa, 0.0, \
                                       relHumidity_w, relHumidity_i );
    if ( isSaved == output::SAVE_FAILURE ) {
        std::cout << " Saving liquid aerosol's properties failed...\n";
        return SAVE_FAIL;
    }

    #endif /* SAVE_LA_MICROPHYS */

    #if ( SAVE_PA_MICROPHYS )

    isSaved = output::Write_MicroPhys( OUT_FILE_PA, saveOutput_PA, saveTime_PA, \
                                       Data.solidAerosol.getBinCenters(), \
                                       m.x(), m.y(), temperature_K, pressure_Pa, 0.0, \
                                       relHumidity_w, relHumidity_i );
    if ( isSaved == output::SAVE_FAILURE ) {
        std::cout << " Saving solid aerosol's properties failed...\n";
        return SAVE_FAIL;
    }

    #endif /* SAVE_PA_MICROPHYS */


    #if ( RINGS && ADJOINT )

        #ifdef OMP
            #pragma omp critical 
            { std::cout << "\n\n ## ON THREAD " << omp_get_thread_num() << ": Starting adjoint calculation...\n"; }
        #else
            std::cout << "\n\n Starting adjoint calculation...\n";
        #endif /* OMP */

    std::copy(sun->CSZA_Vector.begin(), sun->CSZA_Vector.end(), SZA_CST);
    double VAR_OPT[NVAR];

    const std::vector<double> initBackg = ringSpecies.RingAverage( ringArea, totArea, 0 );
    const std::vector<double> finalPlume = ringSpecies.RingAverage( ringArea, totArea, timeArray.size() - 1 );

    IERR = KPP_Main_ADJ( &(finalPlume)[0], &(initBackg)[0],   \
                         temperature_K, pressure_Pa, airDens, \
                         &(timeArray)[0], timeArray.size(),   \
                         KPPADJ_RTOLS, KPPADJ_ATOLS,          \
                         /* Output */ VAR_OPT,                \
                         /* Debug? */ DEBUG_ADJOINT );

    if ( IERR < 0 ) {
        /* Adjoint integration failed */
        return KPPADJ_FAIL;
    }
    
    /* Create ambient struture */
    Ambient adjointData( timeArray.size(), Data.getAmbient(), Data.getAerosol(), Data.getLiqSpecies() );
    adjointData.FillIn( VAR_OPT, 0 );


    /* Perform forward integration with optimized initial conditions */
    const int NADJ = NVAR;
    double Y_adj[NADJ][NVAR], ATOL_adj[NADJ][NVAR], RTOL_adj[NADJ][NVAR];



    /* ---- TOLERANCES ---------------------- */

    for( unsigned int i = 0; i < NVAR; i++ ) {
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

    for( unsigned int i = 0; i < NADJ; i++ ) {
        for( unsigned int j = 0; j < NVAR; j++ ) {
            RTOL_adj[i][j] = 1.0e-5;
            ATOL_adj[i][j] = 1.0e-10;
        }
    }
    
    /* ICNTRL , RCNTRL  = Adjoint settings
     * ISTATUS, RSTATUS = Adjoint statistics */

    double RCNTRL[20], RSTATUS[20];
    int ICNTRL[20], ISTATUS[20];

    /* Default control options */
    for( unsigned int i = 0; i < 20; i++ ) {
        ICNTRL[i] = 0;
        RCNTRL[i] = (double)0.0;
        ISTATUS[i] = 0;
        RSTATUS[i] = (double)0.0;
    }

    ICNTRL[6] = 1;

    /* Reinitialize time */
    curr_Time_s = tInitial_s; /* [s] */
    nTime = 0;

    /* Time loop starts here */

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

        dt = timeArray[nTime+1] - timeArray[nTime]; //UpdateTime( curr_Time_s, tInitial_s, 3600.0*sun->sunRise, 3600.0*sun->sunSet );

        /** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ **/
        /** ~~~~ Update cosine of solar zenight angle ~~~~~ **/
        /** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ **/

        /* Compute the cosize of solar zenith angle midway through the integration step */
        sun->Update( curr_Time_s + dt/2 );

        /* Store cosine of solar zenith angle */
        adjointData.cosSZA[nTime+1] = sun->CSZA;

        if ( printDEBUG ) {
            std::cout << "\n DEBUG : \n";
            std::cout << "         CSZA = " << sun->CSZA << "\n";
        }


        /** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ **/
        /** ~~~~~~~~~~~ Update photolysis rates ~~~~~~~~~~~ **/
        /** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ **/

        for ( unsigned int iPhotol = 0; iPhotol < NPHOTOL; iPhotol++ )
            PHOTOL[iPhotol] = 0.0E+00;

        if ( sun->CSZA > 0.0E+00 )
            Read_JRates( PHOTOL, sun->CSZA );

        if ( printDEBUG ) {
            std::cout << "\n DEBUG : \n";
            for ( unsigned int iPhotol = 0; iPhotol < NPHOTOL; iPhotol++ )
                std::cout << "         PHOTOL[" << iPhotol << "] = " << PHOTOL[iPhotol] << "\n";
        }


        /** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ **/
        /** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ **/
        /** ~~~~~~~~~~~~~~~~~~~~ KPP ~~~~~~~~~~~~~~~~~~~~~~ **/
        /** ~~~~~~~~~ The Kinetics Pre-Processor ~~~~~~~~~~ **/
        /** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ **/
        /** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ **/


        /* Ambient chemistry */
        adjointData.getData( VAR, FIX, aerArray, nTime );

        /* ~~~~~~~~~~~~~~~~~~~~~~~~ */
        /* ~~~~ Chemical rates ~~~~ */
        /* ~~~~~~~~~~~~~~~~~~~~~~~~ */

        /* Update heterogeneous chemistry reaction rates */
        if ( HETCHEMISTRY ) {

            for ( unsigned int iSpec = 0; iSpec < NSPEC; iSpec++ ) {
                HET[iSpec][0] = 0.0E+00;
                HET[iSpec][1] = 0.0E+00;
                HET[iSpec][2] = 0.0E+00;
            }

            relHumidity_Ring = VAR[ind_H2O] * \
                               physConst::kB * temperature_K * 1.00E+06 / \
                               physFunc::pSat_H2Ol( temperature_K );
            GC_SETHET( temperature_K, pressure_Pa, airDens, relHumidity_Ring, \
                       Data.STATE_PSC, VAR, AerosolArea, AerosolRadi, IWC, kheti_sla );

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

        /* Update reaction rates */
        for ( unsigned int iReact = 0; iReact < NREACT; iReact++ )
            RCONST[iReact] = 0.0E+00;

        Update_RCONST( temperature_K, pressure_Pa, airDens, VAR[ind_H2O] );

        /* ~~~~~~~~~~~~~~~~~~~~~~~~ */
        /* ~~~~~ Integration ~~~~~~ */
        /* ~~~~~~~~~~~~~~~~~~~~~~~~ */
 
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
                for ( unsigned int iReact = 0; iReact < NREACT; iReact++ ) {
                    std::cout << "Reaction " << iReact << ": " << RCONST[iReact] << " [molec/cm^3/s]\n";
                }
                std::cout << " ~~~ Printing concentrations:\n";
                for ( unsigned int iSpec = 0; iSpec < NVAR; iSpec++ ) {
                    std::cout << "Species " << iSpec << ": " << VAR[iSpec]/airDens*1.0E+09 << " [ppb]\n";
                }
            }

            delete sun; sun = NULL;
            return KPP_FAIL;
        }

        adjointData.FillIn( VAR, nTime + 1 );
        
        curr_Time_s += dt;
        nTime++;
        
    }

    #if ( SAVE_ADJOINT )
    {
        #pragma omp critical
        {
            isSaved = output::Write_Adjoint( input.fileName_ADJ2char(), \
                                             ringSpecies, ambientData,  \
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

    #endif /* SAVE_ADJOINT */

    #endif /* RINGS */

    /* Clear dynamically allocated variable(s) */
    if ( sun != NULL )
        sun->~SZA();
    
    return SUCCESS;

} /* End of PlumeModel */

/* End of PlumeModel.cpp */
