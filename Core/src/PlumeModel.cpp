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
#include <complex>

#include <ctime> // for clock();

#include "Parameters.hpp"
#include "PhysConstant.hpp"
#include "Interface.hpp"
#include "Monitor.hpp"
#include "SANDS_Solver.hpp"
#if RINGS == 1
    #include "Cluster.hpp"
    #include "Species.hpp"
#endif
#include "Structure.hpp"
#include "Fuel.hpp"
#include "Engine.hpp"
#include "Aircraft.hpp"
#include "Emission.hpp"

typedef std::complex<double> Complex;
typedef std::vector<double> Real_1DVector;
typedef std::vector<Real_1DVector> Real_2DVector;

void BuildMesh( double *x, double *y, \
                double const xlim, double const ylim, \
                unsigned int const nx, unsigned int const ny );
void SZA( double latitude_deg, int dayGMT,\
          double &sunRise, double &sunSet,\
          double &SZASINLAT, double &SZACOSLAT,\
          double &SZASINDEC, double &SZACOSDEC );
void DiffParam( double time, double &d_H, double &d_V );
std::vector<double> BuildTime( double tStart, double tFinal, \
               double sunRise, double sunSet );
double pSat_H2Ol( double T );
double pSat_H2Os( double T );
double UpdateTime( double time, double tStart, \
                   double sunRise, double sunSet );
void CallSolver( Solution& Data, Solver& SANDS );

#ifdef __cplusplus
extern "C" {
#endif
int KPP_Main( double varArray[], double fixArray[], double currentT, double dt, \
              double airDens, double temperature, double pressure, \
              double sinLAT, double cosLAT, double sinDEC, double cosDEC, \
              double rtols, double atols );
#ifdef __cplusplus
}
#endif


int PlumeModel( double temperature_K, double pressure_Pa, \
                double relHumidity_w, double longitude_deg, \
                double latitude_deg )
{

    /* For clock */
    int start_s, stop_s;

    double relHumidity_i = relHumidity_w * pSat_H2Ol( temperature_K )\
                                         / pSat_H2Os( temperature_K );

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
    double x[NX], y[NY];
    BuildMesh(  x,  y, XLIM, YLIM, NX, NY );

    bool fillNegValues = 1;
    double fillWith = 0.0;
    Solver SANDS_Solver( fillNegValues, fillWith );



    /** ~~~~~~~~~~~~~~~~~ **/
    /**        Time       **/
    /** ~~~~~~~~~~~~~~~~~ **/
    /* Define emission and simulation time */
    double tEmission = std::fmod(8.0,24.0); /* [hr] */
    double tInitial = tEmission + TSTART;   /* [hr] */
    double tFinal = tInitial + TSIMUL;      /* [hr] */

    double curr_Time = 3600.0*tInitial; /* [s] */

    /* Create time array */
    std::vector<double> timeArray;
    int nTime = 0;
    timeArray = BuildTime ( 3600.0*tInitial, 3600.0*tFinal, 3600.0*sunRise, 3600.0*sunSet );
    


    /* Rings? */
    if ( RINGS ) {
        /** ~~~~~~~~~~~~~~~~~~ **/
        /**  Cluster of rings  **/
        /** ~~~~~~~~~~~~~~~~~~ **/
        
        /* Create cluster of rings */
        Cluster ringCluster( NRING, ( relHumidity_i > 100.0 ), 0.0, 0.0, 0.0, 0.0 );
      
        /* Print Ring Debug? */
        if ( DEBUG_RINGS )
            ringCluster.Debug();

        /* Allocate species-ring vector */
        SpeciesArray ringSpecies( ringCluster.nRing(), timeArray.size() );

    }
   

    /** ~~~~~~~~~~~~~~~~~ **/
    /**     Emissions     **/
    /** ~~~~~~~~~~~~~~~~~ **/
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
    Emission EI( aircraft.getEngine(), JetA );

    /* Print Emission Debug? */
    if ( DEBUG_EI_INPUT )
        EI.Debug(); 



    /** ~~~~~~~~~~~~~~~~~ **/
    /**     Background    **/
    /** ~~~~~~~~~~~~~~~~~ **/
    /* Assign solution structure */
    Solution Data(N_SPC, NX, NY);

    /* Compute airDens from pressure and temperature */
    double airDens = pressure_Pa / ( kB   * temperature_K ) * 1.00E-06;
    /* [molec/cm3] = [Pa = J/m3] / ([J/K] * [K])            * [m3/cm3] */

    char const *fileName("data/Ambient.txt");
    /* Set solution arrays to ambient data */
    Data.Initialize( fileName, temperature_K, airDens, relHumidity_w );

    /* Print Background Debug? */
    if ( DEBUG_BG_INPUT )
        Data.Debug( airDens );

    std::cout << "NY: " << Data.O3.size() << std::endl;
    std::cout << "NX: " << Data.O3[0].size() << std::endl;

    /* Define initial time step and diffusion and advection arrays */
    double dt;
    dt = UpdateTime( curr_Time, 3600.0*tInitial, 3600.0*sunRise, 3600.0*sunSet );
    SANDS_Solver.UpdateTimeStep( dt );
    
    /** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ **/
    /**     Advection & Diffusion    **/
    /** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ **/
    /* Compute diffusion parameters */
    double d_x, d_y;
    if ( DIFFUSION )
        DiffParam( curr_Time - 3600.0*tInitial, d_x, d_y );
    else {
        d_x = 0.0;
        d_y = 0.0;
    }

    /* Compute advection parameters */
    double v_x, v_y;
    if ( ADVECTION ) {
        v_x = 0; // > 0 means left, < 0 means right
        v_y = 0; // > 0 means upwards, < 0 means downwards
    }
    else {
        v_x = 0;
        v_y = 0;
    }
   
    SANDS_Solver.UpdateDiff( d_x, d_y );
    SANDS_Solver.UpdateAdv( v_x, v_y );

    /* Run FFTW_Wisdom */
    if ( ( nTime == 0 ) && ( FFTW_WISDOM ) ) {
        int start_wisdom, stop_wisdom;
        start_wisdom = clock();
        std::cout << "FFTW_Wisdom..." << std::endl;
        SANDS_Solver.Wisdom( Data.CO2 );
        stop_wisdom = clock();
        std::cout << "time: " << (stop_wisdom-start_wisdom)/double(CLOCKS_PER_SEC) << " [s]" << std::endl;
    }

    double sum = 0;

    double BackG = Data.O3[0][0];
    for ( int i = 0; i < 2; i++ ) {
        if ( i == 0 ) {
            for ( int k = 0; k < NX; k++ ) {
                for ( int l = 0; l < NY; l++ )
                    sum += (Data.O3[l][k] - BackG);
            }
            Data.O3[NY/2][NX/2] += BackG;
            Data.O3[NY/2-1][NX/2] += BackG;
            Data.O3[NY/2][NX/2-1] += BackG;
            Data.O3[NY/2-1][NX/2-1] += BackG;
            sum = 0;
            for ( int k = 0; k < NX; k++ ) {
                for ( int l = 0; l < NY; l++ )
                    sum += (Data.O3[l][k] - BackG);
            }
            std::cout << sum/airDens*1E9 << std::endl;
        }
        start_s = clock();
        
        CallSolver( Data, SANDS_Solver );

        stop_s = clock();
        std::cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << " [ms], " << std::endl;

        if ( i > -1 ) {
            sum = 0;
            for ( int k = 0; k < NX; k++ ) {
                for ( int l = 0; l < NY; l++ )
                    sum += (Data.O3[l][k] - BackG);
            }
//            std::cout << sum/airDens*1E9 << std::endl;
        }
            
        for ( int k = -3; k < 3; k++ ) {
            for ( int l = -3; l < 3; l++ ) 
                std::cout << (Data.O3[NY/2+k][NX/2+l])/airDens*1E9 << ", "; 
            std::cout << std::endl;
        }

    }


    unsigned int iNx = 0; 
    unsigned int jNy = 0;
    double varArray[N_VAR];
    double fixArray[N_FIX];
    Data.GetData( varArray, fixArray, iNx, jNy );
    
    start_s = clock();
    for ( int i = 0; i < 100; i++ ) {
        KPP_Main( varArray, fixArray, curr_Time, 300, \
                  airDens, temperature_K, pressure_Pa, \
                  SZASINLAT, SZACOSLAT, SZASINDEC, SZACOSDEC, \
                  RTOLS, ATOLS );

    }
    stop_s = clock();

    std::cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << " [ms], " << std::endl;
 
 
 
 
 
    return 1;

} /* End of PlumeModel */

void SZA( double latitude_deg, int dayGMT,\
          double &sunRise, double &sunSet,\
          double &SZASINLAT, double &SZACOSLAT,\
          double &SZASINDEC, double &SZACOSDEC )
{
    double const A0 = 0.006918;
    double const A1 = 0.399912;
    double const A2 = 0.006758;
    double const A3 = 0.002697;
    double const B1 = 0.070257;
    double const B2 = 0.000907;
    double const B3 = 0.000148;

    const double PI = 3.141592653589793238460; /* \pi */

    double r_SZA = 2*PI*(floor(dayGMT) - 1)/365.0;

    double DEC = A0 - A1*cos(1*r_SZA) + B1*sin(1*r_SZA)\
                    - A2*cos(2*r_SZA) + B2*sin(2*r_SZA)\
                    - A3*cos(3*r_SZA) + B3*sin(3*r_SZA);

    SZASINLAT = std::sin(latitude_deg*PI/180);
    SZACOSLAT = std::cos(latitude_deg*PI/180);
    SZASINDEC = std::sin(DEC);
    SZACOSDEC = std::cos(DEC);

    sunRise = std::max((12.0 - 180.0/(PI*15.0)*acos(-(SZASINLAT * SZASINDEC)\
                                              / (SZACOSLAT * SZACOSDEC))), 0.0);
    sunSet  = std::min((12.0 + 180.0/(PI*15.0)*acos(-(SZASINLAT * SZASINDEC)\
                                              / (SZACOSLAT * SZACOSDEC))), 24.0);

} /* End of SZA */


