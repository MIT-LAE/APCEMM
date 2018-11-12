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
int PlumeModel( const unsigned int iCase,                   \
                double temperature_K, double pressure_Pa,   \
                double relHumidity_w, double longitude_deg, \
                double latitude_deg,                        \
                const std::vector<double> &inputEmission );

int main( int , char* [] )
{
            
    std::vector<std::vector<double> > parameters;
    unsigned int iCase, nCases;

    #pragma omp master
        {
        
            const bool doPrint = false;

            /* Print welcome message */
            PrintMessage( doPrint );

            /* Read in parameters */
            parameters = ReadParameters();

            /* Number of cases */
            nCases  = parameters[0].size();

        }

    /* Synchronize the threads */
    #pragma omp barrier

    #pragma omp single
        {

            std::cout << "\n Running model for " << nCases << " case(s)\n"; 
            std::cout << " Number of processors: " << omp_get_num_procs() << "\n";
            std::cout << " Number of threads   : " << omp_get_max_threads() << "\n";

        }

    #pragma omp parallel for schedule(dynamic, 1) shared(parameters,nCases)
        for ( iCase = 0; iCase < nCases; iCase++ ) {

            double temperature_K = parameters[0][iCase];
            double pressure_Pa   = parameters[1][iCase];
            double relHumidity_w = parameters[2][iCase];
            double longitude_deg = parameters[3][iCase];
            double latitude_deg  = parameters[4][iCase];
            std::vector<double> emissionIndices( 6, 0.0E+00 );
            emissionIndices[0] = parameters[5][iCase];
            emissionIndices[1] = parameters[6][iCase];
            emissionIndices[2] = parameters[7][iCase];
            emissionIndices[3] = parameters[8][iCase];
            emissionIndices[4] = parameters[9][iCase];
            emissionIndices[5] = parameters[10][iCase];

        #pragma omp critical
            { std::cout << "-> Running case " << iCase << " on thread " << omp_get_thread_num() << "\n"; }


            const unsigned int model = 1;

            /*
             * model = 0 -> Box Model
             *
             * model = 1 -> Plume Model
             *                   + 
             *              Adjoint Model
             *
             * model = 2 -> Adjoint Model
             *
             * model = 3 -> Box Model
             *                  +
             *              Plume Model
             */

            int iERR = 0;

            switch (model) {

                /* Box Model */
                case 0:

                    std::cout << "Not implemented yet\n";
                    break;

                /* Plume Model (APCEMM) */
                case 1:

                    iERR = PlumeModel( iCase, temperature_K, pressure_Pa, \
                                       relHumidity_w, longitude_deg,      \
                                       latitude_deg, emissionIndices );
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

            #pragma omp critical 
            {

                if ( iERR < 0 ) {
                    std::cout.precision(3);
                    std::cout << "\n APCEMM Case: " << iCase << " failed on thread " << omp_get_thread_num() << ".\n";
                    std::cout << " Error: " << iERR << "\n";
                    std::cout << std::fixed;
                    std::cout << std::setprecision(3);
                    std::cout << " T   : " << std::setw(8) << temperature_K << " [K]\n";
                    std::cout << " P   : " << std::setw(8) << pressure_Pa/((double) 100.0) << " [hPa]\n";
                    std::cout << " RH_w: " << std::setw(8) << relHumidity_w << " [%]\n";
                    std::cout << " LON : " << std::setw(8) << longitude_deg << " [deg]\n";
                    std::cout << " LAT : " << std::setw(8) << latitude_deg << " [deg]\n";
                }
                else {
                    std::cout << " APCEMM Case: " << iCase << " completed.\n";
                }

            }

        }

    return 0;

}
