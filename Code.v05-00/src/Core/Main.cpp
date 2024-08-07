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
#include <fstream>
#include <cstdio>
#include <ctime>
#include <filesystem>
#include <unistd.h>
#include <limits.h>
#include <sys/stat.h>
#include <YamlInputReader/YamlInputReader.hpp>
#include "Core/Interface.hpp"
#include "Core/Parameters.hpp"
#include "Core/Input.hpp"
#include "Core/LAGRIDPlumeModel.hpp"
#include "Core/Status.hpp"

static int DIR_FAIL = -9;

void CreateREADME( const std::string folder, const std::string fileName, \
                   const std::string purpose );
void CreateStatusOutput(const std::string folder, const int caseNumber, const SimStatus status);
int PlumeModel( OptInput &Input_Opt, const Input &inputCase );

inline bool exist( const std::string &name )
{

    struct stat buffer;
    return (stat (name.c_str(), &buffer) == 0 );

} /* End of exist */

int main( int argc, char* argv[])
{

    std::vector<std::unordered_map<std::string, double> > parameters;
    unsigned int iCase, nCases;
    const unsigned int iOFFSET = 0;
    
    const unsigned int model = 1;

    #ifdef DEBUG
        std::cout << "-------- DEBUG is enabled --------" << std::endl;
    #endif

    // Set the seed once at the top-level
    setSeed();

    /* Declaring the Input Option object for use in APCEMM */
    OptInput Input_Opt; // Input Option object

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
    if(argc < 2){
        std::cout << "No Input File Detected!" << std::endl;
        std::cout << "Exiting ... " << std::endl;
        return 1;
    }
    else if(argc > 2){
        std::cout << "Unexpected Input: " << argv[2] << std::endl;
        std::cout << "Exiting ... " << std::endl;
        return 1;
    }

    #pragma omp master
    {
        std::string FILESEP = "/";
        std::string FILENAME = argv[1];
        std::filesystem::path INPUT_FILE_PATH(FILENAME);
        INPUT_FILE_PATH = std::filesystem::canonical(INPUT_FILE_PATH);

        YamlInputReader::readYamlInputFile( Input_Opt, INPUT_FILE_PATH.generic_string() );

        /* Collect parameters and create cases */
        parameters = YamlInputReader::generateCases( Input_Opt );

        /* Number of cases */
        nCases  = parameters.size();
        
        /* Create output directory */
        struct stat sb;
        if ( !( stat( Input_Opt.SIMULATION_OUTPUT_FOLDER.c_str(), &sb) == 0 \
                    && S_ISDIR(sb.st_mode) ) ) {

            /* Create directory */
            const int dir_err = \
                    mkdir( Input_Opt.SIMULATION_OUTPUT_FOLDER.c_str(), \
                            S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );

            if ( dir_err == -1 ) {
                std::cout << " Could not create directory: ";
                std::cout << Input_Opt.SIMULATION_OUTPUT_FOLDER << std::endl;
                std::cout << " You may not have write permission" << std::endl;
                exit(1);
            }
            
            /* Create README */
            const std::string description = "";
            CreateREADME( Input_Opt.SIMULATION_OUTPUT_FOLDER, "README", description );

        }
    } /* master CPU */

    /* ====================================================================== */
    /* ---- Synchronize the threads ----------------------------------------- */
    /* ====================================================================== */
    // #pragma omp barrier

    // /* Print number of cases considered */
    // #pragma omp single
    // {
    //     #ifdef OMP 
    //         const char* numberprocs = std::getenv("SLURM_CPUS_ON_NODE");
    //         if ( nCases > 1 )
    //             std::cout << "\n Running model for " << nCases << " cases on ";
    //         else
    //             std::cout << "\n Running model for " << nCases << " case on ";
    //         std::cout << numberprocs << " processors." << std::endl;
    //     #else
    //         if ( nCases > 1 )
    //             std::cout << "\n Running model for " << nCases << " cases." << std::endl; 
    //         else
    //             std::cout << "\n Running model for " << nCases << " case." << std::endl; 
    //     #endif /* OMP */
    // }

    //PARALLEL_CASES = Input_Opt.SIMULATION_PARAMETER_SWEEP;

    /* ====================================================================== */
    /* ---- CASE LOOP STARTS HERE ------------------------------------------- */
    /* ====================================================================== */

    #pragma omp parallel for schedule(dynamic, 1) shared(Input_Opt, parameters, nCases) if( PARALLEL_CASES )
    for ( iCase = 0; iCase < nCases; iCase++ ) {

        unsigned int jCase = iOFFSET + iCase;

        std::string fullPath, fullPath_ADJ, fullPath_BOX, fullPath_micro;
        std::stringstream ss, ss_ADJ, ss_BOX, ss_micro;
        ss     << std::setw(6) << std::setfill('0') << jCase;
        std::string file     = Input_Opt.SIMULATION_FORWARD_FILENAME + ss.str();
        ss_ADJ << std::setw(6) << std::setfill('0') << jCase;
        std::string file_ADJ = Input_Opt.SIMULATION_ADJOINT_FILENAME + ss_ADJ.str();
        ss_BOX << std::setw(6) << std::setfill('0') << jCase;
        std::string file_BOX = Input_Opt.SIMULATION_BOX_FILENAME + ss_BOX.str();
        ss_micro << std::setw(6) << std::setfill('0') << jCase;
        std::string file_micro = "Micro" + ss_micro.str();

        if ( Input_Opt.SIMULATION_OUTPUT_FOLDER.back() == '/' ) {
            fullPath       = Input_Opt.SIMULATION_OUTPUT_FOLDER + file;
            fullPath_ADJ   = Input_Opt.SIMULATION_OUTPUT_FOLDER + file_ADJ;
            fullPath_BOX   = Input_Opt.SIMULATION_OUTPUT_FOLDER + file_BOX;
            fullPath_micro = Input_Opt.SIMULATION_OUTPUT_FOLDER + file_micro;
        } else {
            fullPath       = Input_Opt.SIMULATION_OUTPUT_FOLDER + '/' + file;
            fullPath_ADJ   = Input_Opt.SIMULATION_OUTPUT_FOLDER + '/' + file_ADJ;
            fullPath_BOX   = Input_Opt.SIMULATION_OUTPUT_FOLDER + '/' + file_BOX;
            fullPath_micro = Input_Opt.SIMULATION_OUTPUT_FOLDER + '/' + file_micro;
        }
        fullPath = fullPath + ".nc";
        fullPath_ADJ = fullPath_ADJ + ".nc";
        fullPath_BOX = fullPath_BOX + ".nc";
        fullPath_micro = fullPath_micro + ".out";

        bool fileExist = 0;

        if ( Input_Opt.SIMULATION_ADJOINT ) {
            #pragma omp critical
            { fileExist = exist( fullPath_ADJ ); }
        } else {
            #pragma omp critical
            { fileExist = exist( fullPath ); }
        }

        // Hardcode for now
        std::string author = "Thibaud M. Fritz (fritzt@mit.edu)";

        if ( !fileExist || Input_Opt.SIMULATION_OVERWRITE ) {

            const Input inputCase( iCase, parameters, \
                                   fullPath,          \
                                   fullPath_ADJ,      \
                                   fullPath_BOX,      \
                                   fullPath_micro,    \
                                   author );

            #pragma omp critical
            { 
                std::cout << " -> Running case " << iCase;
                #ifdef OMP
                    std::cout << " on thread " << omp_get_thread_num();
                #endif /* OMP */
                std::cout << "" << std::endl;
            }
            Input_Opt.TS_AERO_FILENAME = "ts_aerosol_case" + std::to_string(iCase) + "_hhmm.nc";

            SimStatus case_status;
            switch (model) {

                /* Box Model */
                case 0:

                    std::cout << "Not implemented yet" << std::endl;
                    break;

                /* Plume Model (APCEMM) */
                case 1: {
                    std::cout << "running epm... " << std::endl;
                    LAGRIDPlumeModel LAGRID_Model(Input_Opt, inputCase);
                    case_status = LAGRID_Model.runFullModel();
                    // iERR = PlumeModel( Input_Opt, inputCase );
                    break;
                    
                }

                /* Adjoint Model */
                case 2:

                    std::cout << "Not implemented yet" << std::endl;
                    break;

                case 3:

                    std::cout << "Not implemented yet" << std::endl;
                    break;

                default:

                    std::cout << "Wrong input for model" << std::endl;
                    std::cout << "model = " << model << "" << std::endl;
                    std::cout << "Value should be between 0 and 3" << std::endl;
                    break;
                    
            }

            #pragma omp critical 
            {
                if ( case_status == SimStatus::Failed ) {
                    std::cout.precision(3);
                    std::cout << "\n APCEMM Case: " << iCase << " failed";
                    #ifdef OMP
                        std::cout << " on thread " << omp_get_thread_num();
                    #endif /* OMP */
                    std::cout << "." << std::endl;
                    // This error reporting is not being used right now
                    // std::cout << " Error: " << iERR << "" << std::endl;
                    std::cout << std::fixed;
                    std::cout << std::setprecision(3);
                    std::cout << " T   : " << std::setw(8) << inputCase.temperature_K();
                    std::cout << " [K]" << std::endl;
                    std::cout << " P   : " << std::setw(8) << inputCase.pressure_Pa()/((double) 100.0);
                    std::cout << " [hPa]" << std::endl;
                    std::cout << " RH_w: " << std::setw(8) << inputCase.relHumidity_w();
                    std::cout << " [%]" << std::endl;
                    std::cout << " LON : " << std::setw(8) << inputCase.longitude_deg();
                    std::cout << " [deg]" << std::endl;
                    std::cout << " LAT : " << std::setw(8) << inputCase.latitude_deg();
                    std::cout << " [deg]" << std::endl;
                }
                else { std::cout << " APCEMM Case: " << iCase << " completed." << std::endl; }
                
                CreateStatusOutput(Input_Opt.SIMULATION_OUTPUT_FOLDER, iCase, case_status);
            }

        }

    }
    
    /* ====================================================================== */
    /* ---- CASE LOOP ENDS HERE --------------------------------------------- */
    /* ====================================================================== */
   
    std::cout << "\n All cases have been completed!" << std::endl;

    /* ====================================================================== */
    /* ---- END NORMALLY ---------------------------------------------------- */
    /* ====================================================================== */

    return 0;


} /* End of Main */

void CreateREADME( const std::string folder, const std::string fileName, const std::string purpose )
{

    std::ofstream README;
    const std::string fullPath = folder + "/" + fileName;
    README.open( fullPath.c_str() );

    README << "############################################################################\
             \n############################################################################\
             \n###                                                                      ###\
             \n###                             APCEMM                                   ###\
             \n###                               --                                     ###\
             \n###   A(ircraft) P(lume) C(hemistry) E(mission) M(icrophysics) M(odel)   ###\
             \n###                                                                      ###\
             \n###                                                                      ###\
             \n###   Version: 5.0                                                       ###\
             \n###   Author : Thibaud M. Fritz                                          ###\
             \n###   Contact: Thibaud M. Fritz (fritzt@mit.edu),                        ###\
             \n###            Sebastian D. Eastham (seastham@mit.edu)                   ###\
             \n###                                                                      ###\
             \n############################################################################\
             \n############################################################################\
             \n###                                                                      ###\
             \n###   This project was funded by NASA and developed at                   ###\
             \n###   the laboratory for Aviation and the Environment,                   ###\
             \n###   Massachusetts Institute of Technology,                             ###\
             \n###   Cambridge, MA, USA                                                 ###\
             \n###                                                                      ###\
             \n############################################################################\
             \n############################################################################\n\n";


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

    /* Print simulation start date */
    std::time_t rawtime;
    std::tm* timeinfo;
    std::time(&rawtime);
    timeinfo = std::localtime(&rawtime);

    char buffer[50];
    std::strftime(buffer,50,"%m/%d/%Y %H:%M",timeinfo);
    README << "\n## Simulation start date " << buffer << "\n";

    /* Print source code directory */
    char cwd[PATH_MAX];
    getcwd(cwd, sizeof(cwd));
    if ( cwd != NULL )
        README << "\n## Source files: " << cwd << "\n";
    else
        std::cout << "\n Failed to get current working directory" << std::endl;

    /* Print destination folder */
    README << "\n## Destination folder: " << folder << "\n";

    /* Getting hostname and username */
    const char* username = std::getenv("SLURM_JOB_USER");
    const char* node = std::getenv("SLURM_NODELIST");
    README << "\n## Running as " << username << " on " << node << "\n";

    /* Printing purpose */
    README << "\n## Purpose: " << purpose;

    /* Print empty lines and force flush */
    README << "\n\n\n" << std::endl;

    README.close();

} /* End of PrintMessage */

void CreateStatusOutput(const std::string folder, const int caseNumber, const SimStatus status)
{
    std::string fileName = "status_case" + std::to_string(caseNumber);
    std::ofstream statusFile;

    const std::string fullPath = folder + "/" + fileName;
    statusFile.open( fullPath.c_str() );

    switch (status)
    {
    case SimStatus::Complete:
        statusFile << "Complete" << std::endl;
        break;
    case SimStatus::Incomplete:
        statusFile << "Incomplete" << std::endl;
        break;
    case SimStatus::NoWaterSaturation:
        statusFile << "NoWaterSaturation" << std::endl;
        break;
    case SimStatus::NoPersistence:
        statusFile << "NoPersistence" << std::endl;
        break;
    case SimStatus::NoSurvivalVortex:
        statusFile << "NoSurvivalVortex" << std::endl;
        break;
    case SimStatus::Failed:
        statusFile << "Failed" << std::endl;
        break;
    default:
        statusFile << "Unknown exit code" << std::endl;
        break;
    }

    statusFile.close();

} /* End of CreateStatusOutput */

/* End of Main.cpp */
