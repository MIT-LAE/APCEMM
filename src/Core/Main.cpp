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

#include "Core/Parameters.hpp"
#include "Core/Input.hpp"

void PrintMessage( bool doPrint );
std::vector<std::vector<double> > ReadParameters( );
int PlumeModel( const unsigned int iCase, \
                const Input &inputCase );

int main( int , char* [] )
{

    std::vector<std::vector<double> > parameters;
    unsigned int iCase, nCases;
    unsigned int jCase;
    const unsigned int iOFFSET = 0;
    
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
        if ( nCases > 1 )
            std::cout << "\n Running model for " << nCases << " cases on " << omp_get_num_procs() << " processors\n"; 
        else
            std::cout << "\n Running model for " << nCases << " case on " << omp_get_num_procs() << " processors\n"; 
    }

    /* Synchronize the threads */
    #pragma omp barrier

    #pragma omp parallel for schedule(dynamic, 1) shared(parameters, nCases)
    for ( iCase = 0; iCase < nCases; iCase++ ) {

        const Input inputCase( iCase, parameters );

        #pragma omp critical
        { std::cout << "-> Running case " << iCase << " on thread " << omp_get_thread_num() << "\n"; }


        int iERR = 0;

        switch (model) {

            /* Box Model */
            case 0:

                std::cout << "Not implemented yet\n";
                break;

            /* Plume Model (APCEMM) */
            case 1:

                jCase = iOFFSET + iCase;
                iERR = PlumeModel( jCase, inputCase );
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
                std::cout << " T   : " << std::setw(8) << inputCase.temperature_K() << " [K]\n";
                std::cout << " P   : " << std::setw(8) << inputCase.pressure_Pa()/((double) 100.0) << " [hPa]\n";
                std::cout << " RH_w: " << std::setw(8) << inputCase.relHumidity_w() << " [%]\n";
                std::cout << " LON : " << std::setw(8) << inputCase.longitude_deg() << " [deg]\n";
                std::cout << " LAT : " << std::setw(8) << inputCase.latitude_deg() << " [deg]\n";
            }
            else { std::cout << " APCEMM Case: " << iCase << " completed.\n"; }
        }

    }
    
    return 0;

} /* End of Main */

void PrintMessage( bool doPrint )
{
    std::string welcomeMessage, authorsMessage;

    welcomeMessage = "\n\n\
            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\
            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\
            ~~~~~~~~                                  ~~~~~~~~\n\
            ~~~~~~~~             APCEMM               ~~~~~~~~\n\
            ~~~~~~~~                                  ~~~~~~~~\n\
            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\
            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\0";

/*
                            _\
                           | \ \
                          _|  \____________________________\
                         /    o  o  o  o  o  o  o  o  o  |_\ \
                         \_________________________________/ \
                                         /    /\
                                        /    /\
                                       /    /\
                                      /____/\0
*/

    authorsMessage = "\
            \n   Version: 5.0\n\
            \n   Author: Thibaud M. Fritz (fritzt@mit.edu),\
            \n     with contributions from:\
            \n         - Sebastian D. Eastham (seastham@mit.edu),\
            \n         - Raymond L. Speth (speth@mit.edu).\n\
            \n   Corresponding author: Sebastian D. Eastham (seastham@mit.edu)\n\
            \n   This project was funded by NASA and developed at \
            \n   the Laboratory for Aviation and the Environment,\
            \n   Department of Aeronautics and Astronautics\
            \n   Massachusetts Institute of Technology\
            \n   Cambridge, MA, USA\n\0";

    std::cout << welcomeMessage << std::endl;

    if ( doPrint )
        std::cout << authorsMessage << std::endl;

} /* End of PrintMessage */

/* End of Main.cpp */

