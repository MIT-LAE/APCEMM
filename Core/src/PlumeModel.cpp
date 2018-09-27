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
/* Working directory    : /home/fritzt/APCEMM-SourceCode            */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <complex>
#include <ctime>

#include "Parameters.hpp"
#include "PhysConstant.hpp"
#include "Interface.hpp"
#include "Monitor.hpp"
#include "SANDS_Solver.hpp"
#if ( RINGS )
    #include "Cluster.hpp"
    #include "Species.hpp"
#endif /* RINGS */
#include "Mesh.hpp"
#include "Structure.hpp"
#include "Ambient.hpp"
#include "Fuel.hpp"
#include "Engine.hpp"
#include "Aircraft.hpp"
#include "Emission.hpp"
#if ( TIME_IT )
    #include "Timer.hpp"
#endif /* TIME_IT */
#if ( SAVE_OUTPUT )
    #include "Save.hpp"
#endif /* SAVE_OUTPUT */

typedef std::complex<double> Complex;
typedef std::vector<double> Real_1DVector;
typedef std::vector<Real_1DVector> Real_2DVector;

void SZA( double latitude_deg, int dayGMT,\
          double &sunRise, double &sunSet,\
          double &SZASINLAT, double &SZACOSLAT,\
          double &SZASINDEC, double &SZACOSDEC );
void DiffParam( double time, double &d_x, double &d_y );
void AdvParam( double time, double &v_x, double &v_y );
void AdvGlobal( double time, double &v_x, double &v_y, double &dTrav_x, double &dTrav_y );
std::vector<double> BuildTime( double tStart, double tFinal, \
               double sunRise, double sunSet );
double UpdateTime( double time, double tStart, \
                   double sunRise, double sunSet );
double pSat_H2Ol( double T );
double pSat_H2Os( double T );
void CallSolver( Solution& Data, Solver& SANDS );

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
int KPP_Main( double varArray[], double fixArray[], double currentT, double dt, \
              double airDens, double temperature, double pressure, \
              double sinLAT, double cosLAT, double sinDEC, double cosDEC, \
              double rtols, double atols );
#ifdef __cplusplus
}
#endif /* __cplusplus */


int PlumeModel( double temperature_K, double pressure_Pa, \
                double relHumidity_w, double longitude_deg, \
                double latitude_deg )
{

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
    double relHumidity_i = relHumidity_w * pSat_H2Ol( temperature_K )\
                                         / pSat_H2Os( temperature_K );

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


    unsigned int dayGMT(81);
    double sunRise, sunSet; /* sun rise and sun set in hours */
    double SZASINLAT, SZACOSLAT, SZASINDEC, SZACOSDEC;

    /* Define sun parameters */
    SZA( latitude_deg, dayGMT,\
         sunRise, sunSet,\
         SZASINLAT, SZACOSLAT, SZASINDEC, SZACOSDEC );
    


    /** ~~~~~~~~~~~~~~~~~ **/
    /**        Mesh       **/
    /** ~~~~~~~~~~~~~~~~~ **/
    
    Mesh m;

    /* Get cell areas */
    std::vector<std::vector<double>> cellAreas;
    cellAreas = m.getAreas();


    /** ~~~~~~~~~~~~~~~~~ **/
    /**        Time       **/
    /** ~~~~~~~~~~~~~~~~~ **/

    /*  
     *  - tEmission is the local emission time expressed in hours 
     *  (between 0.0 and 24.0)
     *  - TSTART corresponds to the time at which the simulation starts
     *  after emission (set to 0). It is expressed in hours
     *  - tInitial is the local time at which the simulation starts in hours
     *  - TSIMUL represents the simulation time (in hours)
     *  - tFinal corresponds to the final time of the simulation expressed in hours
     */ 

    /* Define emission and simulation time */
    double tEmission_h = std::fmod(8.0,24.0);  /* [hr] */
    double tInitial_h  = tEmission_h + TSTART; /* [hr] */
    double tFinal_h    = tInitial_h + TSIMUL;  /* [hr] */
    double tInitial_s  = tInitial_h * 3600.0;  /* [s] */
    double tFinal_s    = tFinal_h   * 3600.0;  /* [s] */

    /* Current time in [s] */
    double curr_Time_s = tInitial_s; /* [s] */
    /* Time step in [s] */
    double dt = 0;                   /* [s] */

    /* Create time array */

    /* Vector of time in [s] */
    std::vector<double> timeArray;
    timeArray = BuildTime ( tInitial_s, tFinal_s, 3600.0*sunRise, 3600.0*sunSet );
    /* Time counter [-] */
    unsigned int nTime = 0;
    

    
    /** ~~~~~~~~~~~~~~~~~ **/
    /**     Background    **/
    /** ~~~~~~~~~~~~~~~~~ **/
    
    /* Assign solution structure */
    Solution Data(N_SPC, NX, NY);

    /* Compute airDens from pressure and temperature */
    double airDens = pressure_Pa / ( kB   * temperature_K ) * 1.00E-06;
    /* [molec/cm3] = [Pa = J/m3] / ([J/K] * [K])            * [m3/cm3] */

    /* Set solution arrays to ambient data */
    Data.Initialize( AMBFILE, temperature_K, airDens, relHumidity_w );

    /* Print Background Debug? */
    if ( DEBUG_BG_INPUT )
        Data.Debug( airDens );

    /* Create ambient struture */
    Ambient ambientData( timeArray.size(), Data.getAmbient() );
    
    
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
    bool fillNegValues = 1;
    /* Fill with? */
    double fillWith = 0.0;

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
    /* varArray stores all the concentrations of variable species */
    double varArray[N_VAR];
    /* fixArray stores all the concentrations of fixed species */
    double fixArray[N_FIX];



    /** ~~~~~~~~~~~~~~~~~ **/
    /**      Rings?       **/
    /** ~~~~~~~~~~~~~~~~~ **/

#if ( RINGS )
        
    /** ~~~~~~~~~~~~~~~~~~ **/
    /**  Cluster of rings  **/
    /** ~~~~~~~~~~~~~~~~~~ **/

    /* Ring index */
    unsigned int iRing = 0;

    /* Number of rings */
    unsigned int nRing;

    /* Create cluster of rings */
    Cluster ringCluster( NRING, ( relHumidity_i > 100.0 ), 0.0, 0.0, 0.0, 0.0 );
    nRing = ringCluster.getnRing();

    /* Print Ring Debug? */
    if ( DEBUG_RINGS )
        ringCluster.Debug();

    /* Allocate species-ring vector */
    SpeciesArray ringSpecies( ringCluster.getnRing(), timeArray.size(), ringCluster.halfRing() );

    /* Compute Grid to Ring mapping */        
    m.Ring2Mesh( ringCluster );
   
    /* Get mapping */
    std::vector<std::vector<std::pair<unsigned int, unsigned int>>> mapRing2Mesh;
    mapRing2Mesh = m.getList();
    
    /* Print ring to mesh mapping? */
    if ( DEBUG_MAPPING )
        m.Debug();
    
    /* Compute ring areas */
    ringCluster.ComputeRingAreas( cellAreas, mapRing2Mesh );
    std::vector<double> ringArea = ringCluster.getRingArea();

    /* Add emission into the grid */
    Data.addEmission( EI, aircraft, mapRing2Mesh, cellAreas, ringCluster.halfRing() ); 
   
    /* Fill in variables species for initial time */
    ringSpecies.FillIn( Data, m, nTime );

    /* Allocate an additional array for KPP */
    double tempArray[N_VAR];

 
    /* Otherwise we do not have a ring structure and chemistry is solved on the grid */
#else
    
    /* Grid indices */
    unsigned int iNx = 0;
    unsigned int jNy = 0;

#endif /* RINGS */
   

    std::cout << "\n\n *** Time loop starts now ***";

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
        dt = UpdateTime( curr_Time_s, tInitial_s, 3600.0*sunRise, 3600.0*sunSet );
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

        CallSolver( Data, SANDS_Solver );

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

                std::copy( varArray, varArray + N_VAR, tempArray );

                /* Call KPP */
                KPP_Main( varArray, fixArray, curr_Time_s, dt, \
                      airDens, temperature_K, pressure_Pa, \
                      SZASINLAT, SZACOSLAT, SZASINDEC, SZACOSDEC, \
                      RTOLS, ATOLS );
                
                ringSpecies.FillIn( varArray, nTime + 1, iRing );
   
                Data.applyRing( varArray, tempArray, mapRing2Mesh, iRing );
                
            }

            /* Ambient chemistry */
            ambientData.getData( varArray, fixArray, nTime );

            /* Call KPP */
            KPP_Main( varArray, fixArray, curr_Time_s, dt, \
                  airDens, temperature_K, pressure_Pa, \
                  SZASINLAT, SZACOSLAT, SZASINDEC, SZACOSDEC, \
                  RTOLS, ATOLS );

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


                    /* Call KPP */
                    KPP_Main( varArray, fixArray, curr_Time_s, dt, \
                          airDens, temperature_K, pressure_Pa, \
                          SZASINLAT, SZACOSLAT, SZASINDEC, SZACOSDEC, \
                          RTOLS, ATOLS );

                    /* Convert KPP output back to data structure */
                    Data.applyData( varArray, iNx, jNy );
                }
            }

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
        for ( unsigned int iNx = 0; iNx < NX; iNx++ ) {
            for ( unsigned int jNy = 0; jNy < NY; jNy++ ) {
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
        std::cout << "\n    " << " ~~> Emitted NOy: " << std::setw(6) << mass_Emitted_NOy * 1.0E+06 / Na * MW_N * 1.0E+06 << " [g(N)/km] ";
        /*                                                           [molec/cm3 * m2] * [m3/cm3]/[molec/mole]*[kg/mole]*[g/kg*m/km] = [g/km] */
    
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
        for ( unsigned int iNx = 0; iNx < NX; iNx++ ) {
            for ( unsigned int jNy = 0; jNy < NY; jNy++ ) {
                mass_Emitted_CO2 += ( Data.CO2[jNy][iNx] - mass_Ambient_CO2 ) * cellAreas[jNy][iNx];
            }
        }

        std::cout << "\n\n    " << " *** CO2 mass check: ";

        std::cout << "\n    " << " ~~> Emitted CO2: " << std::setw(6) << mass_Emitted_CO2 * 1.0E+06 / Na * MW_CO2 * 1.0E+03 << " [kg/km]   ";
        /*                                                           [molec/cm3 * m2] * [m3/cm3]/[molec/mole]*[kg/mole]*[m/km] = [kg/km] */

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

    int isSaved = output::Write( ringSpecies, ambientData, ringCluster, timeArray, temperature_K, pressure_Pa, airDens, relHumidity_w, relHumidity_i, longitude_deg, latitude_deg, sunRise, sunSet );
    if ( isSaved == output::SAVE_FAILURE ) {
        std::cout << "Saving file failed...\n";
        return -2;
    }

#endif /* SAVE_OUTPUT */


    return 1;

} /* End of PlumeModel */

