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
#ifndef FMT_HEADER_ONLY
#define FMT_HEADER_ONLY
#endif

#include <iostream>
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
#include <fmt/core.h>
#include "YamlInputReader/YamlInputReader.hpp"
#include "Core/Parameters.hpp"
#include "Core/Input.hpp"
#include "Core/LAGRIDPlumeModel.hpp"
#include "Core/Status.hpp"
#include "Core/Diag_Mod.hpp"
#include "Util/MC_Rand.hpp"

void CreateREADME( const std::string folder, const std::string fileName, \
                   const std::string purpose );
void CreateStatusOutput(const std::string folder, const SimStatus status);
int PlumeModel( OptInput &Input_Opt, const Input &inputCase );

/* One run owns one output folder. Refuse to write into a folder that already
 * holds the output of an earlier run unless overwriting is enabled. Every run
 * writes a status file, and every run that reaches the plume model writes
 * epm-output.nc, so those two catch a finished run. The timeseries and the EPM
 * microphysics file catch a run that was interrupted before its status was
 * written. */
bool holdsAPCEMMOutput( const std::string &folder )
{

    std::error_code error;
    if ( !std::filesystem::is_directory( folder, error ) )
        return false;

    for ( const auto &entry : std::filesystem::directory_iterator( folder, error ) ) {
        if ( !entry.is_regular_file() )
            continue;
        const std::string name = entry.path().filename().string();
        if ( name == "status" || name == "epm-output.nc" || name == "Micro.out" )
            return true;
        if ( name.starts_with( "ts" ) && entry.path().extension() == ".nc" )
            return true;
    }
    return false;

} /* End of holdsAPCEMMOutput */

int main( int argc, char* argv[])
{

    const unsigned int model = 1;

    #ifdef DEBUG
        std::cout << "-------- DEBUG is enabled --------" << std::endl;
    #endif

    /* Model configuration and the scenario to simulate */
    OptInput Input_Opt;
    Input scenario;

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
        std::cout << "No Input Files Detected!" << std::endl;
        std::cout << "Exiting ... " << std::endl;
        return 1;
    }

    /* Later input files override values from the earlier ones. This layers one
     * run, it does not create several. */
    std::vector<std::string> INPUT_FILE_PATHS;
    for (int i = 1; i < argc; ++i) {
      std::string FILENAME = argv[i];
      std::filesystem::path INPUT_FILE_PATH(FILENAME);
      INPUT_FILE_PATH = std::filesystem::canonical(INPUT_FILE_PATH);
      INPUT_FILE_PATHS.push_back(INPUT_FILE_PATH);
    }

    YamlInputReader::readYamlInputFiles( Input_Opt, scenario, INPUT_FILE_PATHS );

    if (Input_Opt.ADV_SAVE_PSD_GRID){
        Diag::set_storePSD(true);
    }

    // Set the seed once at the top-level
    setSeed( Input_Opt );

    if ( holdsAPCEMMOutput( Input_Opt.SIMULATION_OUTPUT_FOLDER ) \
            && !Input_Opt.SIMULATION_OVERWRITE ) {
        std::cout << " Output folder already holds APCEMM output: ";
        std::cout << Input_Opt.SIMULATION_OUTPUT_FOLDER << std::endl;
        std::cout << " Point the run at an empty folder, or set";
        std::cout << " 'Overwrite if folder exists (T/F)' to T." << std::endl;
        std::cout << "Exiting ... " << std::endl;
        return 1;
    }

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
            return 1;
        }

        /* Create README */
        const std::string description = "";
        CreateREADME( Input_Opt.SIMULATION_OUTPUT_FOLDER, "README", description );

    }

    /* The switch is a placeholder for future models. Only the plume model is
     * implemented, so every other arm leaves the status at Failed. */
    SimStatus status = SimStatus::Failed;
    switch (model) {

        /* Box Model */
        case 0:

            std::cout << "Not implemented yet" << std::endl;
            break;

        /* Plume Model (APCEMM) */
        case 1: {
            std::cout << "running epm... " << std::endl;
            LAGRIDPlumeModel LAGRID_Model(Input_Opt, scenario);
            status = LAGRID_Model.runFullModel();
            // iERR = PlumeModel( Input_Opt, scenario );
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

    if ( status == SimStatus::Failed ) {
        std::cout << "\n APCEMM failed." << std::endl;
        // This error reporting is not being used right now
        // std::cout << " Error: " << iERR << "" << std::endl;

        // Report contrail location
        // {:>8.2f} = right align with a width of 8 with 2 decimals
        fmt::print(" LON [deg]: {:>8.2f}\n", scenario.longitude_deg());
        fmt::print(" LAT [deg]: {:>8.2f}\n", scenario.latitude_deg());
        fmt::print(" P   [hPa]: {:>8.2f}\n", scenario.pressure_Pa()/100.0);

        // Report relevant input met data when crashing
        if (Input_Opt.MET_LOADMET)
        {
            fmt::print(" Met file : {:>}\n",  Input_Opt.MET_FILENAME);
        }
    }
    else { std::cout << "\n APCEMM completed." << std::endl; }

    CreateStatusOutput(Input_Opt.SIMULATION_OUTPUT_FOLDER, status);

    /* Every physical outcome is a success. Only a failure of the model itself
     * reports a non-zero exit code, so a batch driver can rely on it.
     * Incomplete means the run reached its configured Plume Process end time
     * with the contrail still alive, which is a normal result. */
    return status == SimStatus::Failed ? 1 : 0;

} /* End of Main */

// TODO rewrite all of this with std::filesystem
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
    // Calling with args (NULL, 0) allocates a new buffer of the correct size
    char *cwd = getcwd(NULL, 0);
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

void CreateStatusOutput(const std::string folder, const SimStatus status)
{
    std::ofstream statusFile;

    const std::string fullPath = folder + "/status";
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
