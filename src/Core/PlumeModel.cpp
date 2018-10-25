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

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <complex>
#include <ctime>
#include "omp.h"

#include "Core/Parameters.hpp"
#include "Core/Interface.hpp"
#include "Core/Monitor.hpp"
#include "SANDS/Solver.hpp"
#include "AIM/Coagulation.hpp"
#include "EPM/Integrate.hpp"
#include "KPP/KPP.hpp"
#include "KPP/KPP_Parameters.h"
#include "KPP/KPP_Global.h"
#if ( RINGS )
    #include "Core/Cluster.hpp"
    #include "Core/Species.hpp"
#endif /* RINGS */
#include "Core/SZA.hpp"
#include "Core/Mesh.hpp"
#include "Core/Structure.hpp"
#include "Core/Ambient.hpp"
#include "Core/Fuel.hpp"
#include "Core/Engine.hpp"
#include "Core/Aircraft.hpp"
#include "Core/Emission.hpp"
#if ( TIME_IT )
    #include "Core/Timer.hpp"
#endif /* TIME_IT */
#if ( SAVE_OUTPUT )
    #include "Core/Save.hpp"
#endif /* SAVE_OUTPUT */

static int SUCCESS   =  1;
static int KPP_FAIL  = -1;
static int SAVE_FAIL = -2;

typedef std::complex<double> Complex;
typedef std::vector<double> Real_1DVector;
typedef std::vector<Real_1DVector> Real_2DVector;

void DiffParam( double time, double &d_x, double &d_y );
void AdvParam( double time, double &v_x, double &v_y );
void AdvGlobal( double time, double &v_x, double &v_y, double &dTrav_x, double &dTrav_y );
std::vector<double> BuildTime( double tStart, double tFinal, \
               double sunRise, double sunSet );
double UpdateTime( double time, double tStart, \
                   double sunRise, double sunSet );
void CallSolver( Solution& Data, Solver& SANDS );


int PlumeModel( double temperature_K, double pressure_Pa, \
                double relHumidity_w, double longitude_deg, \
                double latitude_deg )
{

    const bool DBG = 0;

    /* Grid indices */
    unsigned int iNx = 0;
    unsigned int jNy = 0;

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

    int IERR;

    /* Compute relative humidity w.r.t ice */
    double relHumidity_i = relHumidity_w * physFunc::pSat_H2Ol( temperature_K )\
                                         / physFunc::pSat_H2Os( temperature_K );

//    /* If ICE_MICROPHYSICS is turned on and the domain is supersaturated, break the x-symmetry */
//#if ( ICE_MICROPHYSICS && X_SYMMETRY )
//
//    /* Is supersaturated? */
//    if ( relHumidity_i > 100.0 ) {
//
//#undef X_SYMMETRY
//#define X_SYMMETRY                  0     /* Set x-symmetry to 0 */
//
//    }
//#endif /* ICE_MICROPHYSICS && X_SYMMETRY */


    const unsigned int dayGMT(81);

    /* Define sun parameters */
    SZA sun( latitude_deg, dayGMT );
    


    /** ~~~~~~~~~~~~~~~~~ **/
    /**        Mesh       **/
    /** ~~~~~~~~~~~~~~~~~ **/
    
    Mesh m;

    /* Get cell areas */
    const std::vector<std::vector<double>> cellAreas = m.getAreas();

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
    const double tEmission_h = std::fmod(8.0, 24.0); /* [hr] */
    const double tInitial_h  = tEmission_h;          /* [hr] */
    const double tFinal_h    = tInitial_h + TSIMUL;  /* [hr] */
    const double tInitial_s  = tInitial_h * 3600.0;  /* [s] */
    const double tFinal_s    = tFinal_h   * 3600.0;  /* [s] */

    /* Current time in [s] */
    double curr_Time_s = tInitial_s; /* [s] */
    /* Time step in [s] */
    double dt = 0;                   /* [s] */

    /* Create time array */

    /* Vector of time in [s] */
    const std::vector<double> timeArray = BuildTime ( tInitial_s, tFinal_s, 3600.0*sun.sunRise, 3600.0*sun.sunSet );

    /* Time counter [-] */
    unsigned int nTime = 0;
    

    
    /** ~~~~~~~~~~~~~~~~~ **/
    /**     Background    **/
    /** ~~~~~~~~~~~~~~~~~ **/

    /* Declare solution structure */
    Solution Data;

    /* Compute airDens from pressure and temperature */
    double airDens = pressure_Pa / ( physConst::kB   * temperature_K ) * 1.00E-06;
    /* [molec/cm3] = [Pa = J/m3] / ([J/K]            * [K]           ) * [m3/cm3] */

    /* Set solution arrays to ambient data */
    Data.Initialize( AMBFILE, temperature_K, pressure_Pa, airDens, relHumidity_w, latitude_deg, DBG );

    /* Print Background Debug? */
    if ( DEBUG_BG_INPUT || DBG )
        Data.Debug( airDens );

    /* Create ambient struture */
    Ambient ambientData( timeArray.size(), Data.getAmbient(), Data.getAerosol(), Data.getLiqSpecies() );
   
    /** ~~~~~~~~~~~~~~~~~ **/
    /**      Solver       **/
    /** ~~~~~~~~~~~~~~~~~ **/
    
    /* Allocate horizontal and vertical diffusion parameters */
    double d_x, d_y;

    /* Allocate horizontal and vertical advection parameters */
    double v_x, v_y; /* These can vary from one species/particle to another */
    double vGlob_x, vGlob_y; /* These correspond to domain-wide advection velocities (updraft, downdraft) */

    /* Allocate horizontal and vertical distance traveled */
    double dTrav_x, dTrav_y;

    /* Fill negative values? */
    const bool fillNegValues = 1;
    /* Fill with? */
    const double fillWith = 0.0;

    /* Allocate Solver */
    Solver SANDS_Solver( fillNegValues, fillWith );
    
    /* Run FFTW_Wisdom? */
    if ( FFTW_WISDOM ) {
        std::cout << "FFTW_Wisdom..." << "\n";
        SANDS_Solver.Wisdom( Data.CO2 );
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
    const Fuel JetA( ChemFormula );

    /* Define aircraft */
    char const *aircraftName("B747-800");
    const Aircraft aircraft( aircraftName, temperature_K, pressure_Pa, relHumidity_w );

    /* Print AC Debug? */
    if ( DEBUG_AC_INPUT || DBG )
        aircraft.Debug();

    /* Aggregate emissions from engine and fuel characteristics */
    const Emission EI( aircraft.getEngine(), JetA );

    /* Print Emission Debug? */
    if ( DEBUG_EI_INPUT || DBG )
        EI.Debug(); 

    
    
    /** ~~~~~~~~~~~~~~~~~ **/
    /**     Chemistry     **/
    /** ~~~~~~~~~~~~~~~~~ **/

    /* Allocate arrays for KPP */
    /* varArray stores all the concentrations of variable species */
    double varArray[NVAR];
    
    /* fixArray stores all the concentrations of fixed species */
    double fixArray[NFIX];
   
    /* aerArray stores all the number concentrations of aerosols */
    double aerArray[N_AER][2];

    /* Ambient chemistry */
    ambientData.getData( varArray, fixArray, aerArray, nTime );

    
    
    /** ~~~~~~~~~~~~~~~~~~~~~~~ **/
    /**    Early Microphysics   **/
    /** ~~~~~~~~~~~~~~~~~~~~~~~ **/

    double Ice_rad, Ice_den, Soot_den, H2O_mol, SO4g_mol, SO4l_mol;
    double areaPlume; 
    AIM::Aerosol liquidAer, iceAer;
    EPM::Integrate( temperature_K, pressure_Pa, relHumidity_w, varArray, fixArray, aerArray, aircraft, EI, \
                    Ice_rad, Ice_den, Soot_den, H2O_mol, SO4g_mol, SO4l_mol, liquidAer, iceAer, areaPlume );

    /* Compute initial plume area.
     * If 2 engines, we assume that after 3 mins, the two plume haven't fully mixed yet and result in a total
     * area of 2 * the area computed for one engine
     * If 3 or more engines, we assume that the plumes originating from the same wing have mixed. */

    areaPlume *= 2;
    double semiYaxis = 0.5*aircraft.getVortexdeltaz1();
    double semiXaxis = areaPlume/(physConst::PI*0.5*aircraft.getVortexdeltaz1());
   

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
    if ( DEBUG_RINGS | DBG )
        ringCluster.Debug();

    /* Allocate species-ring vector */
    SpeciesArray ringSpecies( nRing, timeArray.size(), ringCluster.halfRing() );

    /* Compute Grid to Ring mapping */        
    m.Ring2Mesh( ringCluster );
   
    /* Get mapping */
    const std::vector<std::vector<std::pair<unsigned int, unsigned int>>> mapRing2Mesh = m.getList();
    
    /* Print ring to mesh mapping? */
    if ( DEBUG_MAPPING || DBG )
        m.Debug();
    
    /* Compute ring areas */
    ringCluster.ComputeRingAreas( cellAreas, mapRing2Mesh );
    const std::vector<double> ringArea = ringCluster.getRingArea();

    /* Add emission into the grid */
    Data.addEmission( EI, aircraft, mapRing2Mesh, cellAreas, ringCluster.halfRing(), temperature_K, ( relHumidity_i > 100.0 ), \
                      liquidAer, iceAer ); 
   
    /* Fill in variables species for initial time */
    ringSpecies.FillIn( Data, m, nTime );

    /* Allocate an additional array for KPP */
    double tempArray[NVAR];

    /* Allocate a ring's relative humidity */
    double relHumidity_Ring;
 
    /* Otherwise we do not have a ring structure and chemistry is solved on the grid */

#endif /* RINGS */

    bool IS_PSC = 0;
    double frac_gSO4 = 0.0E+00;
    double areaHET[NAERO];
    double radiHET[NAERO];
    double iwcHET = 0;
    double kheti_sla[11];

#pragma omp critical
    {
        std::cout << "\n\n ## ON THREAD: " << omp_get_thread_num() << "\n ##";

        const unsigned int coutPrecision = 5;
        const unsigned int txtWidth      = coutPrecision + 2;
        std::cout << std::setprecision(coutPrecision);
        std::cout << "\n ## ATMOSPHERIC COND.:";
        std::cout << "\n ##\n";
        std::cout << " ## - Temperature: " << std::setw(txtWidth) << temperature_K          << " [  K]\n";
        std::cout << " ## - Pressure   : " << std::setw(txtWidth) << pressure_Pa * 1.00E-02 << " [hPa]\n";
        std::cout << " ## - Rel. Hum. I: " << std::setw(txtWidth) << relHumidity_i          << " [  %]\n";
        std::cout << " ## - Latitude   : " << std::setw(txtWidth) << latitude_deg           << " [deg]\n";
        std::cout << " ## - Max CSZA   : " << std::setw(txtWidth) << sun.CSZA_max           << " [ - ]\n";

        std::cout << "\n ## EMISSIONS:";
        std::cout << "\n ##\n";
        std::cout << " ## - E_CO2 = " << std::setw(txtWidth) << EI.getCO2() * aircraft.getFuelFlow() / aircraft.getVFlight()           << " [kg(CO2)/km]"\
            " ( EI = " << std::setw(txtWidth) << EI.getCO2() * 1.00E-03 << " [kg/kg_fuel] )\n";
        std::cout << " ## - E_CO  = " << std::setw(txtWidth) << EI.getCO()  * aircraft.getFuelFlow() / aircraft.getVFlight() * 1.0E+03 << " [ g(CO) /km]"\
            " ( EI = " << std::setw(txtWidth) << EI.getCO()             << " [ g/kg_fuel] )\n";
        std::cout << " ## - E_CH4 = " << std::setw(txtWidth) << EI.getCH4() * aircraft.getFuelFlow() / aircraft.getVFlight() * 1.0E+06 << " [mg(CH4)/km]"\
            " ( EI = " << std::setw(txtWidth) << EI.getCH4() * 1.00E+03 << " [mg/kg_fuel] )\n";
        std::cout << " ## - E_NOx = " << std::setw(txtWidth) << ( EI.getNO() / MW_NO + EI.getNO2() / MW_NO2 + EI.getHNO2() / MW_HNO2 ) * MW_N \
                                                  * aircraft.getFuelFlow() / aircraft.getVFlight() * 1.0E+03 << " [ g(N)  /km]"\
            " ( EI = " << std::setw(txtWidth) << EI.getNOx()            << " [ g/kg_fuel] )\n";
        std::cout << " ## - E_SO2 = " << std::setw(txtWidth) << EI.getSO2() * aircraft.getFuelFlow() / aircraft.getVFlight() * 1.0E+03 << " [ g(SO2)/km]"\
            " ( EI = " << std::setw(txtWidth) << EI.getSO2()            << " [ g/kg_fuel], \n";
        std::cout << " ##                                 FSC = " << std::setw(txtWidth) << JetA.getFSC() << " [-]          )\n";
        std::cout << " ## - E_Soo = " << std::setw(txtWidth) << EI.getSoot()* aircraft.getFuelFlow() / aircraft.getVFlight() * 1.0E+03 << " [ g(Soo)/km]"\
            " ( EI = " << std::setw(txtWidth) << EI.getSoot()* 1.00E+03 << " [mg/kg_fuel] )\n";
        std::cout << " ## - rSoot = " << std::setw(txtWidth) << EI.getSootRad() * 1.0E+09 << " [nm] \n";

        std::cout << "\n ## AEROSOLS:";
        std::cout << "\n ##\n";
        std::cout << " ## - LA : " << std::setw(txtWidth+3) << liquidAer.Moment() << " [#/cm^3], \n";
        std::cout << " ##        " << std::setw(txtWidth+3) << liquidAer.getEffRadius() * 1.0E+09 << " [nm], \n";
        std::cout << " ##        " << std::setw(txtWidth+3) << liquidAer.Moment(2) * 1.0E+12 << " [mum^2/cm^3] \n";

        std::cout << "\n ## BACKG COND.:";
        std::cout << "\n ##\n";
        std::cout << " ## - NOx = " << std::setw(txtWidth) << ( varArray[ind_NO] + varArray[ind_NO2] ) / airDens * 1.0E+12 << " [ppt]\n";
        std::cout << " ## - O3  = " << std::setw(txtWidth) << ( varArray[ind_O3] ) / airDens * 1.0E+09 << " [ppb]\n";
        std::cout << " ## - CO  = " << std::setw(txtWidth) << ( varArray[ind_CO] ) / airDens * 1.0E+09 << " [ppb]\n";
        std::cout << " ##\n";
        std::cout << " ## - LA : " << std::setw(txtWidth+3) << Data.LA_nDens << " [#/cm^3], \n";
        std::cout << " ##        " << std::setw(txtWidth+3) << Data.LA_rEff  << " [nm], \n";
        std::cout << " ##        " << std::setw(txtWidth+3) << Data.LA_SAD   << " [mum^2/cm^3] \n";
        std::cout << " ##\n";
        std::cout << " ## - PA : " << std::setw(txtWidth+3) << Data.PA_nDens << " [#/cm^3], \n";
        std::cout << " ##        " << std::setw(txtWidth+3) << Data.PA_rEff  << " [nm], \n";
        std::cout << " ##        " << std::setw(txtWidth+3) << Data.PA_SAD   << " [mum^2/cm^3] \n";

    }

    /** ~~~~~~~~~~~~~~~~~~~~~~~~~~ **/
    /**         Time Loop          **/
    /** ~~~~~~~~~~~~~~~~~~~~~~~~~~ **/

#if ( TIME_IT )

    Stopwatch_cumul.Start( );

#endif /* TIME_IT */

    while ( curr_Time_s < tFinal_s ) {

        /* Print message */
        std::cout << "\n";
        std::cout << "\n - Time step: " << nTime << " out of " << timeArray.size();
        std::cout << "\n -> Solar time: " << std::fmod( curr_Time_s/3600.0, 24.0 ) << " [hr]";

        /** ~~~~~~~~~~~~~~~~~~~~~~~~~~ **/
        /**      Update Time Step      **/
        /** ~~~~~~~~~~~~~~~~~~~~~~~~~~ **/
        
        /* Compute time step */
        dt = UpdateTime( curr_Time_s, tInitial_s, 3600.0*sun.sunRise, 3600.0*sun.sunSet );
        SANDS_Solver.UpdateTimeStep( dt );

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
            

            AdvParam( curr_Time_s - tInitial_s, v_x, v_y ); 

        }
        else {

            /* If advection is turned off, set advection parameters to 0 */
            
            vGlob_x = 0;
            vGlob_y = 0;

        }
   
        /* Update diffusion and advection arrays */
        SANDS_Solver.UpdateDiff( d_x, d_y );
        SANDS_Solver.UpdateAdv ( v_x, v_y );

        
        /** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ **/
        /** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ **/
        /** ~~~~~~~~~~~~~~~~~~~ SANDS ~~~~~~~~~~~~~~~~~~~~~ **/
        /** ~~~ Spectral Advection aNd Diffusion Solver ~~~ **/
        /** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ **/
        /** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ **/
        
#if ( TIME_IT )

        Stopwatch.Start( reset );

#endif /* TIME_IT */

#pragma omp critical /* Not sure why omp critical is needed here, otherwise leads to segmentation faults... */
        {
            CallSolver( Data, SANDS_Solver );
        }

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

        /* Compute SUN */
        sun.Update( curr_Time_s );

        /* Store cosine of solar zenith angle */
        ambientData.cosSZA[nTime] = sun.CSZA;

        if ( DBG ) {
            std::cout << "\n DEBUG : \n";
            std::cout << "         CSZA = " << sun.CSZA << "\n";
        }


        /** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ **/
        /** ~~~~~~~~~~~ Update photolysis rates ~~~~~~~~~~~ **/
        /** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ **/

        for ( unsigned int iPhotol = 0; iPhotol < NPHOTOL; iPhotol++ )
            PHOTOL[iPhotol] = 0.0E+00;

        if ( sun.CSZA > 0.0E+00 )
            Read_JRates( PHOTOL, sun.CSZA );

        if ( DBG ) {
            std::cout << "\n DEBUG : \n";
            for ( unsigned int iPhotol = 0; iPhotol < NPHOTOL; iPhotol++ )
                std::cout << "         PHOTOL[" << iPhotol << "] = " << PHOTOL[iPhotol] << "\n";
        }


        /** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ **/
        /** ~~~~~~~~~~~~~~~~~~~~ KPP ~~~~~~~~~~~~~~~~~~~~~~ **/
        /** ~~~~~~~~~ The Kinetics Pre-Processor ~~~~~~~~~~ **/
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

                /* Convert ring structure to KPP inputs (varArray and fixArray) */
                ringSpecies.getData( varArray, fixArray, nTime + 1, iRing );

                for ( unsigned int iSpec = 0; iSpec < NVAR; iSpec++ )
                    tempArray[iSpec] = varArray[iSpec];

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

                    relHumidity_Ring = varArray[ind_H2O] * \
                                       physConst::kB * temperature_K * 1.00E+06 / \
                                       physFunc::pSat_H2Ol( temperature_K );
                    GC_SETHET( temperature_K, pressure_Pa, airDens, relHumidity_Ring, \
                               IS_PSC, varArray, areaHET, radiHET, iwcHET, kheti_sla );
                }
             
                /* Update reaction rates */
                for ( unsigned int iReact = 0; iReact < NREACT; iReact++ )
                    RCONST[iReact] = 0.0E+00;

                Update_RCONST( temperature_K, pressure_Pa, airDens, varArray[ind_H2O] );

                /* ~~~~~~~~~~~~~~~~~~~~~~~~ */
                /* ~~~~~ Integration ~~~~~~ */
                /* ~~~~~~~~~~~~~~~~~~~~~~~~ */

                IERR = KPP_Main( varArray, fixArray, curr_Time_s, dt, \
                                 KPP_RTOLS, KPP_ATOLS );

                if ( IERR < 0 ) {
                    /* Integration failed */

                    std::cout << "Integration failed for ring = " << iRing << " at time t = " << curr_Time_s << "( nTime = " << nTime << " )\n";
                    std::cout << " ~~~ Printing reaction rates:\n";
                    for ( unsigned int iReact = 0; iReact < NREACT; iReact++ ) {
                        std::cout << "Reaction " << iReact << ": " << RCONST[iReact] << " [molec/cm^3/s]\n";
                    }
                    std::cout << " ~~~ Printing concentrations:\n";
                    for ( unsigned int iSpec = 0; iSpec < NVAR; iSpec++ ) {
                        std::cout << "Species " << iSpec << ": " << varArray[iSpec]/airDens*1.0E+09 << " [ppb]\n";
                    }

                    return KPP_FAIL;
                }
                
                ringSpecies.FillIn( varArray, nTime + 1, iRing );
            
                Data.applyRing( varArray, tempArray, mapRing2Mesh, iRing );

            }

            /* Ambient chemistry */
            ambientData.getData( varArray, fixArray, aerArray, nTime );

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

                relHumidity_Ring = varArray[ind_H2O] * \
                                   physConst::kB * temperature_K * 1.00E+06 / \
                                   physFunc::pSat_H2Ol( temperature_K );
                GC_SETHET( temperature_K, pressure_Pa, airDens, relHumidity_Ring, \
                           IS_PSC, varArray, areaHET, radiHET, iwcHET, kheti_sla );
            }

            /* Update reaction rates */
            for ( unsigned int iReact = 0; iReact < NREACT; iReact++ )
                RCONST[iReact] = 0.0E+00;

            Update_RCONST( temperature_K, pressure_Pa, airDens, varArray[ind_H2O] );

            /* ~~~~~~~~~~~~~~~~~~~~~~~~ */
            /* ~~~~~ Integration ~~~~~~ */
            /* ~~~~~~~~~~~~~~~~~~~~~~~~ */

            IERR = KPP_Main( varArray, fixArray, curr_Time_s, dt, \
                             KPP_RTOLS, KPP_ATOLS );

            if ( IERR < 0 ) {
                /* Integration failed */

                std::cout << "Integration failed for ring = " << iRing << " at time t = " << curr_Time_s << "( nTime = " << nTime << " )\n";
                std::cout << " ~~~ Printing reaction rates:\n";
                for ( unsigned int iReact = 0; iReact < NREACT; iReact++ ) {
                    std::cout << "Reaction " << iReact << ": " << RCONST[iReact] << " [molec/cm^3/s]\n";
                }
                std::cout << " ~~~ Printing concentrations:\n";
                for ( unsigned int iSpec = 0; iSpec < NVAR; iSpec++ ) {
                    std::cout << "Species " << iSpec << ": " << varArray[iSpec]/airDens*1.0E+09 << " [ppb]\n";
                }

                return KPP_FAIL;
            }

            ambientData.FillIn( varArray, nTime + 1 );
                
            Data.applyAmbient( varArray, mapRing2Mesh, nRing );

        }


        /* Otherwise solve chemistry on the grid */
#else

        /* Is chemistry turned on? */
        if ( CHEMISTRY ) {

            for ( iNx = 0; iNx < NX; iNx++ ) {
                for ( jNy = 0; jNy < NY; jNy++ ) {

                    /* Convert data structure to KPP inputs (varArray and fixArray) */
                    Data.getData( varArray, fixArray, iNx, jNy );

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

                        relHumidity_Ring = varArray[ind_H2O] * \
                                           physConst::kB * temperature_K * 1.00E+06 / \
                                           physFunc::pSat_H2Ol( temperature_K );
                        GC_SETHET( temperature_K, pressure_Pa, airDens, relHumidity_Ring, \
                                   IS_PSC, varArray, areaHET, radiHET, iwcHET, kheti_sla );
                    }

                    /* Update reaction rates */
                    for ( unsigned int iReact = 0; iReact < NREACT; iReact++ )
                        RCONST[iReact] = 0.0E+00;

                    Update_RCONST( temperature_K, pressure_Pa, airDens, varArray[ind_H2O] );

                    /* ~~~~~~~~~~~~~~~~~~~~~~~~ */
                    /* ~~~~~ Integration ~~~~~~ */
                    /* ~~~~~~~~~~~~~~~~~~~~~~~~ */

                    IERR = KPP_Main( varArray, fixArray, curr_Time_s, dt, \
                                     KPP_RTOLS, KPP_ATOLS );
                    
                    if ( IERR < 0 ) {
                        /* Integration failed */

                        std::cout << "Integration failed for ring = " << iRing << " at time t = " << curr_Time_s << "( nTime = " << nTime << " )\n";
                        std::cout << " ~~~ Printing reaction rates:\n";
                        for ( unsigned int iReact = 0; iReact < NREACT; iReact++ ) {
                            std::cout << "Reaction " << iReact << ": " << RCONST[iReact] << " [molec/cm^3/s]\n";
                        }
                        std::cout << " ~~~ Printing concentrations:\n";
                        for ( unsigned int iSpec = 0; iSpec < NVAR; iSpec++ ) {
                            std::cout << "Species " << iSpec << ": " << varArray[iSpec]/airDens*1.0E+09 << " [ppb]\n";
                        }

                        return KPP_FAIL;
                    }

                    /* Convert KPP output back to data structure */
                    Data.applyData( varArray, iNx, jNy );
                }
            }
            
            /* Ambient chemistry */
            ambientData.getData( varArray, fixArray, aerArray, nTime );
            
            /* ~~~~~~~~~~~~~~~~~~~~~~~~ */
            /* ~~~~ Chemical rates ~~~~ */
            /* ~~~~~~~~~~~~~~~~~~~~~~~~ */

            /* Update reaction rates */
            for ( unsigned int iReact = 0; iReact < NREACT; iReact++ )
                RCONST[iReact] = 0.0E+00;
            
            Update_RCONST( temperature_K, pressure_Pa, airDens, varArray[ind_H2O], sun.CSZA );

            /* ~~~~~~~~~~~~~~~~~~~~~~~~ */
            /* ~~~~~~ Integration ~~~~~ */
            /* ~~~~~~~~~~~~~~~~~~~~~~~~ */

            IERR = KPP_Main( varArray, fixArray, curr_Time_s, dt, \
                             KPP_RTOLS, KPP_ATOLS );

            if ( IERR < 0 ) {
                /* Integration failed */

                std::cout << "Integration failed for ring = " << iRing << " at time t = " << curr_Time_s << "( nTime = " << nTime << " )\n";
                std::cout << " ~~~ Printing reaction rates:\n";
                for ( unsigned int iReact = 0; iReact < NREACT; iReact++ ) {
                    std::cout << "Reaction " << iReact << ": " << RCONST[iReact] << " [molec/cm^3/s]\n";
                }
                std::cout << " ~~~ Printing concentrations:\n";
                for ( unsigned int iSpec = 0; iSpec < NVAR; iSpec++ ) {
                    std::cout << "Species " << iSpec << ": " << varArray[iSpec]/airDens*1.0E+09 << " [ppb]\n";
                }

                return KPP_FAIL;
            }

            ambientData.FillIn( varArray, nTime + 1 );

        }

#endif /* RINGS */

        
#if ( TIME_IT )
        
        Stopwatch.Stop( );
        KPP_clock = Stopwatch.Elapsed( );

#endif /* TIME_IT */

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
    
#if ( SAVE_OUTPUT )

    int isSaved = output::Write( ringSpecies, ambientData, ringCluster, timeArray, temperature_K, pressure_Pa, airDens, relHumidity_w, relHumidity_i, longitude_deg, latitude_deg, sun.sunRise, sun.sunSet );
    if ( isSaved == output::SAVE_FAILURE ) {
        std::cout << "Saving file failed...\n";
        return SAVE_FAIL;
    }

#endif /* SAVE_OUTPUT */


    return SUCCESS;

} /* End of PlumeModel */

