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
#include <sys/stat.h>

#include "Core/Interface.hpp"
#include "Core/Parameters.hpp"
#include "Core/Input.hpp"

static int DIR_FAIL = -9;

void PrintMessage( bool doPrint );
std::vector<std::vector<double> > ReadParameters( );
int PlumeModel( const Input &inputCase );

inline bool exist( const std::string &name )
{

    struct stat buffer;
    return (stat (name.c_str(), &buffer) == 0 );

} /* End of exist */

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
        
    struct stat sb;
    if (!(stat( OUT_PATH.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode))) {
        const int dir_err = mkdir( OUT_PATH.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
        if ( dir_err == -1 ) {
            std::cout << " Could not create directory: " << OUT_PATH << "\n";
            return DIR_FAIL; 
        }
    }


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

        jCase = iOFFSET + iCase;

        std::string fullPath, fullPath_ADJ;

        std::stringstream ss, ss_ADJ;
        ss << std::setw(5) << std::setfill('0') << jCase;
        std::string file = "APCEMM_Case_" + ss.str();
        ss_ADJ << std::setw(5) << std::setfill('0') << jCase;
        std::string file_ADJ = "APCEMM_ADJ_Case_" + ss_ADJ.str();
        if ( OUT_PATH.back() == '/' ) {
            fullPath = OUT_PATH + file;
            fullPath_ADJ = OUT_PATH + file_ADJ;
        } else {
            fullPath = OUT_PATH + '/' + file;
            fullPath_ADJ = OUT_PATH + '/' + file_ADJ;
        }
        fullPath = fullPath + ".nc";
        fullPath_ADJ = fullPath_ADJ + ".nc";

        bool fileExist = 0;
        
        #pragma omp critical 
        { fileExist = exist( fullPath_ADJ ); }

        if ( !fileExist || REBUILD ) {
            
            const Input inputCase( iCase, parameters, fullPath, fullPath_ADJ );

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

                    iERR = PlumeModel( inputCase );
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
