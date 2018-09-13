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
#include "Structure.hpp"
#include "Fuel.hpp"
#include "Engine.hpp"

#ifdef __cplusplus
extern "C"{
#endif

#include "cheader.h"

#ifdef __cplusplus
}
#endif

typedef std::complex<double> Complex;

void BuildMesh( double *x, double *y, \
                double const xlim, double const ylim, \
                unsigned int const nx, unsigned int const ny );
void BuildFreq( double *kx, double *ky, \
                double *kxx, double *kyy, \
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
void SaveWisdomFile( std::vector<std::vector<double> >& Data, const char* fileName );
void Assign_diffFactor( std::vector<std::vector<double> >& diffFactor, \
                        double d_x, double d_y, \
                        double kxx[], double kyy[], double dt );
void Assign_advFactor( std::vector<std::vector<Complex > >& advFactor, \
                        double v_x, double v_y, \
                        double kx[], double ky[], double dt );
void SPCDiffusion( Solution& Data, \
                   std::vector<std::vector<double> >& diffFactor, \
                   std::vector<std::vector<Complex > >& advFactor, \
                   const char* fileName_FFTW );


int PlumeModel( double temperature_K, double pressure_Pa, \
                 double relHumidity_w, double longitude_deg, \
                 double latitude_deg )
{

    /* For clock */
    int start_s, stop_s;

    double x[NX], y[NY];
    double kx[NX], ky[NY], kxx[NX], kyy[NY];
    
    BuildMesh(  x,  y, XLIM, YLIM, NX, NY );
    BuildFreq( kx, ky, kxx, kyy, XLIM, YLIM, NX, NY );

    double relHumidity_i = relHumidity_w * pSat_H2Ol( temperature_K )\
                                         / pSat_H2Os( temperature_K );

    unsigned int dayGMT(81);
    double sunRise, sunSet; /* sun rise and sun set in hours */
    double SZASINLAT, SZACOSLAT, SZASINDEC, SZACOSDEC;

    SZA( latitude_deg, dayGMT,\
         sunRise, sunSet,\
         SZASINLAT, SZACOSLAT, SZASINDEC, SZACOSDEC );

    double tEmission = std::fmod(8.0,24.0);
    double tInitial = tEmission + TSTART;
    double tFinal = tInitial + TSIMUL;
    double curr_Time = 3600.0*tInitial;

    std::vector<double> timeArray;
    int nTime = 0;

    /* Create time array */
    timeArray = BuildTime ( 3600.0*tInitial, 3600.0*tFinal, 3600.0*sunRise, 3600.0*sunSet );
    //nTime = timeArray.size();

    /* Define fuel */
    char const *ChemFormula("C12H24");
    Fuel JetA( ChemFormula );

    /* Define engine */
    char const *EngineFile("./data/ENG_EI.txt");
    char const *EngineName("GE90-90B");
    std::cout << "Engine name: " << EngineName << std::endl;
    double vFlight = 250; 
    Engine eng( EngineFile, EngineName, temperature_K, pressure_Pa, relHumidity_w, vFlight );

    /* Assign solution structure */
    Solution Data(N_SPC, NX, NY);

    /* Compute airDens from pressure and temperature */
    double airDens = pressure_Pa / ( kB * temperature_K ) * 1.00E-06;
    /* [molec/cm3] = [Pa = J/m3] / ([J/K] * [K])         * [m3/cm3] */

    char const *fileName("data/Ambient.txt");
    /* Set solution arrays to ambient data */
    Data.Initialize( fileName, temperature_K, airDens, relHumidity_w );

    /* 2D Diffusion and advection arrays */
    std::vector<std::vector<double> > diffFactor;
    std::vector<std::vector<Complex > > advFactor;
    for ( unsigned int i = 0; i < NY; i++ ) {
        diffFactor.push_back( std::vector<double>( NX ) );
        advFactor .push_back( std::vector<Complex >( NX ) );
    }

    std::cout << "NY: " << Data.O3.size() << std::endl;
    std::cout << "NX: " << Data.O3[0].size() << std::endl;

    /* Compute diffusion parameters */
    double d_x, d_y;
    if ( DIFFUSION == 1 )
        DiffParam( curr_Time - 3600.0*tInitial, d_x, d_y );
    else {
        d_x = 0.0;
        d_y = 0.0;
    }

    /* Compute advection parameters */
    double v_x, v_y;
    if ( ADVECTION == 1 ) {
        v_x = 0; // > 0 means left, < 0 means right
        v_y = 0; // > 0 means upwards, < 0 means downwards
    }
    else {
        v_x = 0;
        v_y = 0;
    }

    /* Define initial time step and diffusion and advection arrays */
    double dt;
    dt = UpdateTime( curr_Time, 3600.0*tInitial, 3600.0*sunRise, 3600.0*sunSet );
    Assign_diffFactor( diffFactor, d_x, d_y, kxx, kyy, dt );
    Assign_advFactor ( advFactor , v_x, v_y, kx , ky , dt );
    std::cout << "Diffusion param: " << d_x << ", " << d_y << " [m^2/s]" << std::endl;
    std::cout << "Advection param: " << v_x << ", " << v_y << " [m/s]" << std::endl;
    
    /* Run FFTW_Wisdom */
    const char *fileName_FFTW("data/FFTW_Wisdom.txt");
    if ( ( nTime == 0 ) && ( FFTW_WISDOM ) ) {
        int start_wisdom, stop_wisdom;
        start_wisdom = clock();
        std::cout << "FFTW_Wisdom..." << std::endl;
        SaveWisdomFile( Data.CO2, fileName_FFTW );
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
        
        SPCDiffusion( Data, diffFactor, advFactor, fileName_FFTW );

        stop_s = clock();

        std::cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << " [ms], ";
        if ( i > -1 ) {
            sum = 0;
            for ( int k = 0; k < NX; k++ ) {
                for ( int l = 0; l < NY; l++ )
                    sum += (Data.O3[l][k] - BackG);
            }
            std::cout << sum/airDens*1E9 << std::endl;
        }
    }


    unsigned int iNx = 0; 
    unsigned int jNy = 0;
    double varArray[N_VAR];
    double fixArray[N_FIX];
    Data.GetData( varArray, fixArray, iNx, jNy );
    std::cout << varArray[ind_CO2]/airDens*1E6 << std::endl;
    std::cout << varArray[ind_O3]/airDens*1E9 << std::endl;
    
    start_s = clock();
    for ( int i = 0; i < 100; i++ ) {
        KPP_Main( varArray, fixArray, curr_Time, 300, \
                  airDens, temperature_K, pressure_Pa, \
                  SZASINLAT, SZACOSLAT, SZASINDEC, SZACOSDEC, \
                  RTOLS, ATOLS );

    }
    stop_s = clock();

    std::cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << " [ms], " << std::endl;
    std::cout << varArray[ind_CO2]/airDens*1E6 << std::endl;
    std::cout << varArray[ind_O3]/airDens*1E9 << std::endl;
 
 
 
 
 
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

void Assign_diffFactor( std::vector<std::vector<double> >& diffFactor, double d_x, double d_y, \
                        double kxx[], double kyy[], double dt )
{

    for ( unsigned int i = 0; i < NX; i++ ) {
        for ( unsigned int j = 0; j < NY; j++ )
            diffFactor[j][i] = exp( dt * ( d_x * kxx[i] + d_y * kyy[j] ) );
    }

} /* End of Assign_diffFactor */

void Assign_advFactor( std::vector<std::vector<Complex > >& advFactor, double v_x, double v_y, \
                        double kx[], double ky[], double dt )
{

    Complex _1j (0.0, 1.0);
    for ( unsigned int i = 0; i < NX; i++ ) {
        for ( unsigned int j = 0; j < NY; j++ )
            advFactor[j][i] = exp( _1j * dt * ( v_x * kx[i] + v_y * ky[j] ) ); 
    }

} /* End of Assign_advFactor */

