/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Main Program File                                                */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : Main.cpp                                  */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <string>
#include <vector>
#include "omp.h"

void PrintMessage( bool doPrint );
std::vector<std::vector<double> > ReadParameters( );
int PlumeModel( double temperature_K, double pressure_Pa, \
                double relHumidity_w, double longitude_deg, \
                double latitude_deg );

int main( int , char* [] )
{

    const bool doPrint = false;
    std::vector<std::vector<double> > parameters;
    unsigned int nCases;

    /* Print welcome message */
    PrintMessage( doPrint );

    /* Read in parameters */
    parameters = ReadParameters();

    /* Number of cases */
    nCases  = parameters[0].size();

    std::cout << "Running model for " << nCases << " cases on " << omp_get_num_procs() << " procs\n";

//    double temperature_K;
//    double pressure_Pa;
//    double relHumidity_w; 
//    double longitude_deg;
//    double latitude_deg;

#pragma omp parallel for
    for ( unsigned int i = 0; i < nCases; i++ ) {

#pragma omp critical
        {
            std::cout << "Running case " << i << " on thread " << omp_get_thread_num() << "\n";
        }

        double temperature_K = parameters[0][i];
        double pressure_Pa   = parameters[1][i];
        double relHumidity_w = parameters[2][i];
        double longitude_deg = parameters[3][i];
        double latitude_deg  = parameters[4][i];

        const unsigned int model = 1;
        /*
        * model = 0 -> Box     Model
        * model = 1 -> Plume   Model (APCEMM)
        * model = 2 -> Adjoint Model
        * model = 3 -> Box + Plume Model
        */

        int iERR = 0;

        switch (model) {

            /* Box Model */
            case 0:

                std::cout << "Not implemented yet\n";
                break;

            /* Plume Model (APCEMM) */
            case 1:

                iERR = PlumeModel( temperature_K, pressure_Pa, \
                                   relHumidity_w, longitude_deg, \
                                   latitude_deg );

                break;

            /* Adjoint Model */
            case 2:

                std::cout << "Not implemented yet\n";
                break;

            case 3:

                std::cout << "Not implemented yet\n";
                break;

            default:

                std::cout << "Wrong input for model\n";
                std::cout << "model = " << model << "\n";
                std::cout << "Value should be between 0 and 3\n";
                break;
                
        }

        if ( iERR < 0 ) {
            std::cout.precision(3);
            std::cout << "APCEMM Case: " << i << " failed." << "\n";
            std::cout << "Error: " << iERR << "\n";
            std::cout << std::fixed;
            std::cout << std::setprecision(3);
            std::cout << "T   : " << std::setw(8) << temperature_K << " [K]\n";
            std::cout << "P   : " << std::setw(8) << pressure_Pa/((double) 100.0) << " [hPa]\n";
            std::cout << "RH_w: " << std::setw(8) << relHumidity_w << " [%]\n";
            std::cout << "LON : " << std::setw(8) << longitude_deg << " [deg]\n";
            std::cout << "LAT : " << std::setw(8) << latitude_deg << " [deg]\n";
        }
        else {
            std::cout << "APCEMM Case: " << i << " completed.\n";
        }

    }

    return 0;

}
