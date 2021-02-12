/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* ReadInput Program File                                           */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 12/10/2018                                */
/* File                 : ReadInput.cpp                             */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "Core/ReadInput.hpp"

const unsigned int SUCCESS = 1;

const bool VERBOSE = false;
const char* EMPTY = "";
const char* SPACE = " ";
const char* TAB   = "\t";
const char* COLON = ":";
const char* COMMA = ",";
const char* FILENAME = "input.apcemm";
const char* FILESEP = "/";
std::ifstream inputFile;
const unsigned int FIRSTCOL = 25;

std::string line;


void Read_Input_File( OptInput &Input_Opt )
{

    /* Read\_Input\_File is the driver program for reading
     * APCEMM input file "input.apcemm" from disk */

    std::string TOPTITLE;
    const std::string HEADER(80, '*');
    bool RC;
    std::string fullPath = "";

    if ( const char* simDir = std::getenv("APCEMM_runDir") )
        fullPath += simDir;
    else {
        std::cout << " \n Simulation Directory is not defined!" << std::endl;
        const char* currDir = std::getenv("PWD");
        fullPath += currDir;
        std::cout << " Reading from PWD: " << currDir << std::endl;
        std::cout << " For future runs, make sure that the variable";
        std::cout << " 'APCEMM_runDir' is exported" << std::endl;
    }

    fullPath += FILESEP;
    fullPath += FILENAME;

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Read\_Input\_File begins here !
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    /* Echo output */
    std::cout << "\n" << HEADER << std::endl;
    std::cout << " A P C E M M   U S E R   I N P U T" << std::endl;
    std::cout << HEADER << std::endl;

    /* Assume success */
    RC = SUCCESS;

    /* Open file */
    std::cout << "\n Reading from: " << fullPath << "\n" << std::endl;
    inputFile.open( fullPath.c_str() );
    if ( !inputFile ) {
        /* Call error */
        std::cout << " APCEMM I/O ERROR: No such file in directory" << std::endl;
        exit(1);
    }

    /* Read TOPTITLE */
    getline( inputFile, TOPTITLE, '\n' );

    /* Print to the console? */
    if ( VERBOSE )
        std::cout << TOPTITLE << std::endl;

    /* Loop until EOF */
    while ( getline( inputFile, line, '\n' ) ) {

        /* Read a line from the file or exit if EOF */

        /* Print to the console? */
        if ( VERBOSE )
            std::cout << " " << line << std::endl;

        if ( strstr( line.c_str(), "SIMULATION MENU" ) != NULL ) {

            /* ================================ */
            /* ==== Read SIMULATION MENU ====== */
            /* ================================ */

            Read_Simulation_Menu( Input_Opt, RC );

        } else if ( strstr( line.c_str(), "PARAMETER SWEEP" ) != NULL ) {

            /* ================================ */
            /* ==== Read PARAMETER SWEEP ====== */
            /* ================================ */

            Read_Parameters( Input_Opt, RC );

        } else if ( strstr( line.c_str(), "TRANSPORT MENU" ) != NULL ) {

            /* ================================ */
            /* ==== Read TRANSPORT MENU ======= */
            /* ================================ */

            Read_Transport_Menu( Input_Opt, RC );

        } else if ( strstr( line.c_str(), "CHEMISTRY MENU" ) != NULL ) {

            /* ================================ */
            /* ==== Read CHEMISTRY MENU ======= */
            /* ================================ */

            Read_Chemistry_Menu( Input_Opt, RC );

        } else if ( strstr( line.c_str(), "AEROSOL MENU" ) != NULL ) {

            /* ================================ */
            /* ==== Read AEROSOL MENU ========= */
            /* ================================ */

            Read_Aerosol_Menu( Input_Opt, RC );

        } else if ( strstr( line.c_str(), "METEOROLOGY MENU" ) != NULL ) {

            /* ================================ */
            /* ==== Read METEOROLOGY MENU ===== */
            /* ================================ */

            Read_Meteorology_Menu( Input_Opt, RC );

        } else if ( strstr( line.c_str(), "DIAGNOSTIC MENU" ) != NULL ) {

            /* ================================ */
            /* ==== Read DIAGNOSTIC MENU ====== */
            /* ================================ */

            Read_Diagnostic_Menu( Input_Opt, RC );

        } else if ( strstr( line.c_str(), "TIMESERIES MENU" ) != NULL ) {

            /* ================================ */
            /* ==== Read TIMESERIES MENU ====== */
            /* ================================ */

            Read_Timeseries_Menu( Input_Opt, RC );

        } else if ( strstr( line.c_str(), "PROD & LOSS MENU" ) != NULL ) {

            /* ================================ */
            /* ==== Read PROD&LOSS MENU ======= */
            /* ================================ */

            Read_PL_Menu( Input_Opt, RC );

        } else if ( strstr( line.c_str(), "END OF FILE" ) != NULL ) {
            std::cout << " " << line << std::endl;
            break;
        }

    }

    getline( inputFile, line, '\n' );
    std::cout << " " << line << std::endl;

    /* Close input file */
    inputFile.close();

    /* Check that input files do not conflict.
     * If a conflict is detected, abort! */
    Are_Flags_Valid( Input_Opt );

    return;

} /* End of Read_Input_File */

std::vector<std::string> Split_Line( std::string line2split, const std::string delimiter )
{

    /* DESCRIPTION: Function Split\_Line separates a string into substrings
     * according to the delimiter */

    std::vector<std::string> substring;

    size_t pos = 0;
    std::string token;

    /* ==================================================== */
    /* Split_Line begins here!                              */
    /* ==================================================== */

    while ( (pos = line2split.find(delimiter) ) != std::string::npos ) {
        token = line2split.substr(0, pos);
        if ( token.length() > 0 )
            substring.push_back(token);
        line2split.erase(0, pos + delimiter.length());
    }
    substring.push_back(line2split);

    return substring;

} /* End of Split_Line */

void Read_Simulation_Menu( OptInput &Input_Opt, bool &RC )
{

    /* DESCRIPTION: Function Read\_Simulation\_Menu reads the SIMULATION MENU
     * section of the APCEMM input file. */

    /* INPUT/OUTPUT PARAMETERS:
     * - Input_Opt: Input options
     * - RC: Success or failure
     */

    /* Read until all lines from the menu are read */

    std::vector<std::string> tokens;

    /* ==================================================== */
    /* Parameter sweep?                                     */
    /* ==================================================== */

    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable range */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    if ( ( strcmp(tokens[0].c_str(), "T" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "t" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "TRUE" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "true" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "True" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "YES" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Y" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "y" )    == 0 ) )
        Input_Opt.SIMULATION_PARAMETER_SWEEP = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "f" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "FALSE" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "false" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "False" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "NO" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "No" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "no" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "N" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "n" )     == 0 ) )
        Input_Opt.SIMULATION_PARAMETER_SWEEP = 0;
    else {
        std::cout << " Wrong input for: " << "Parameter sweep?" << std::endl;
        exit(1);
    }

    /* ==================================================== */
    /* Monte Carlo?                                         */
    /* ==================================================== */

    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    if ( ( strcmp(tokens[0].c_str(), "T" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "t" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "TRUE" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "true" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "True" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "YES" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Y" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "y" )    == 0 ) )
        Input_Opt.SIMULATION_MONTECARLO = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "f" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "FALSE" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "false" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "False" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "NO" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "No" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "no" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "N" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "n" )     == 0 ) )
        Input_Opt.SIMULATION_MONTECARLO = 0;
    else {
        std::cout << " Wrong input for " << "Monte Carlo?" << std::endl;
        exit(1);
    }

    /* ==================================================== */
    /* Number of runs?                                      */
    /* ==================================================== */

    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    try {
        Input_Opt.SIMULATION_MCRUNS = std::stoi( tokens[0] );
        if ( Input_Opt.SIMULATION_MCRUNS < 1 ) {
            std::cout << " Wrong input for " << "MC runs" << std::endl;
            std::cout << " Number of runs needs to be strictly positive" << std::endl;
            exit(1);
        }
    } catch(std::exception& e) {
        std::cout << " Could not convert string '" << tokens[0] << "' to int for " << "MC runs" << std::endl;
        exit(1);
    }

    /* ==================================================== */
    /* Output folder                                        */
    /* ==================================================== */

    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    Input_Opt.SIMULATION_OUTPUT_FOLDER = tokens[0];

    /* ==================================================== */
    /* Overwrite?                                           */
    /* ==================================================== */

    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable range */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    if ( ( strcmp(tokens[0].c_str(), "T" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "t" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "TRUE" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "true" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "True" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "YES" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Y" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "y" )    == 0 ) )
        Input_Opt.SIMULATION_OVERWRITE = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "f" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "FALSE" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "false" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "False" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "NO" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "No" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "no" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "N" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "n" )     == 0 ) )
        Input_Opt.SIMULATION_OVERWRITE = 0;
    else {
        std::cout << " Wrong input for: " << "Overwrite?" << std::endl;
        exit(1);
    }

    /* ==================================================== */
    /* Run directory                                        */
    /* ==================================================== */

    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    Input_Opt.SIMULATION_RUN_DIRECTORY = tokens[0];

    /* ==================================================== */
    /* Use threaded FFT?                                    */
    /* ==================================================== */

    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable range */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    if ( ( strcmp(tokens[0].c_str(), "T" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "t" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "TRUE" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "true" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "True" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "YES" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Y" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "y" )    == 0 ) )
        Input_Opt.SIMULATION_THREADED_FFT = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "f" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "FALSE" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "false" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "False" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "NO" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "No" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "no" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "N" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "n" )     == 0 ) )
        Input_Opt.SIMULATION_THREADED_FFT = 0;
    else {
        std::cout << " Wrong input for: " << "Use threaded FFT?" << std::endl;
        exit(1);
    }

    /* ==================================================== */
    /* Use FFTW WISDOM?                                     */
    /* ==================================================== */

    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable range */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    if ( ( strcmp(tokens[0].c_str(), "T" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "t" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "TRUE" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "true" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "True" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "YES" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Y" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "y" )    == 0 ) )
        Input_Opt.SIMULATION_USE_FFTW_WISDOM = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "f" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "FALSE" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "false" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "False" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "NO" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "No" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "no" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "N" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "n" )     == 0 ) )
        Input_Opt.SIMULATION_USE_FFTW_WISDOM = 0;
    else {
        std::cout << " Wrong input for: " << "Use FFTW WISDOM?" << std::endl;
        exit(1);
    }

    /* ==================================================== */
    /* Directory with write permission                      */
    /* ==================================================== */

    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    Input_Opt.SIMULATION_DIRECTORY_W_WRITE_PERMISSION = tokens[0];

    /* Check that folder exists and that user has write permissions */
    if ( Input_Opt.SIMULATION_USE_FFTW_WISDOM ) {

        /* Create output directory */
        struct stat sb;
        if ( !( stat( Input_Opt.SIMULATION_DIRECTORY_W_WRITE_PERMISSION.c_str(), &sb) == 0 \
                    && S_ISDIR(sb.st_mode) ) ) {

            /* Create directory */
            const int dir_err = \
                    mkdir( Input_Opt.SIMULATION_DIRECTORY_W_WRITE_PERMISSION.c_str(), \
                            S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );

            if ( dir_err == -1 ) {
                std::cout << " Could not create directory: ";
                std::cout << Input_Opt.SIMULATION_DIRECTORY_W_WRITE_PERMISSION << std::endl;
                std::cout << " You may not have write permission" << std::endl;
                std::cout << " Turning SIMULATION_USE_FFTW_WISDOM off" << std::endl;
                Input_Opt.SIMULATION_USE_FFTW_WISDOM = 0;
            }

        }

        if ( Input_Opt.SIMULATION_DIRECTORY_W_WRITE_PERMISSION.back() != '/' )
            Input_Opt.SIMULATION_DIRECTORY_W_WRITE_PERMISSION += '/';

    } else
        Input_Opt.SIMULATION_DIRECTORY_W_WRITE_PERMISSION = "";

    /* ==================================================== */
    /* Input background condition                           */
    /* ==================================================== */

    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    Input_Opt.SIMULATION_INPUT_BACKG_COND = tokens[0];

    /* ==================================================== */
    /* Save forward results                                 */
    /* ==================================================== */

    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    if ( ( strcmp(tokens[0].c_str(), "T" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "t" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "TRUE" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "true" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "True" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "YES" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Y" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "y" )    == 0 ) )
        Input_Opt.SIMULATION_SAVE_FORWARD = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "f" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "FALSE" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "false" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "False" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "NO" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "No" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "no" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "N" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "n" )     == 0 ) )
        Input_Opt.SIMULATION_SAVE_FORWARD = 0;
    else {
        std::cout << " Wrong input for: " << "Save forward results" << std::endl;
        exit(1);
    }

    /* ==================================================== */
    /* netCDF file name                                     */
    /* ==================================================== */

    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    tokens[0].erase(std::remove(tokens[0].begin(), tokens[0].end(), '*'), tokens[0].end());
    Input_Opt.SIMULATION_FORWARD_FILENAME = tokens[0];

    /* ==================================================== */
    /* Adjoint Optimization                                 */
    /* ==================================================== */

    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    if ( ( strcmp(tokens[0].c_str(), "T" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "t" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "TRUE" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "true" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "True" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "YES" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Y" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "y" )    == 0 ) )
        Input_Opt.SIMULATION_ADJOINT = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "f" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "FALSE" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "false" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "False" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "NO" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "No" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "no" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "N" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "n" )     == 0 ) )
        Input_Opt.SIMULATION_ADJOINT = 0;
    else {
        std::cout << " Wrong input for: " << "Adjoint Optimization" << std::endl;
        exit(1);
    }

    /* ==================================================== */
    /* netCDF file name                                     */
    /* ==================================================== */

    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    tokens[0].erase(std::remove(tokens[0].begin(), tokens[0].end(), '*'), tokens[0].end());
    Input_Opt.SIMULATION_ADJOINT_FILENAME = tokens[0];

    /* ==================================================== */
    /* Run box model                                        */
    /* ==================================================== */

    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    if ( ( strcmp(tokens[0].c_str(), "T" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "t" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "TRUE" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "true" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "True" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "YES" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Y" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "y" )    == 0 ) )
        Input_Opt.SIMULATION_BOXMODEL = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "f" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "FALSE" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "false" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "False" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "NO" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "No" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "no" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "N" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "n" )     == 0 ) )
        Input_Opt.SIMULATION_BOXMODEL = 0;
    else {
        std::cout << " Wrong input for: " << "Run box model" << std::endl;
        exit(1);
    }

    /* ==================================================== */
    /* netCDF file name                                     */
    /* ==================================================== */

    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    tokens[0].erase(std::remove(tokens[0].begin(), tokens[0].end(), '*'), tokens[0].end());
    Input_Opt.SIMULATION_BOX_FILENAME = tokens[0];


    /* Return success */
    RC = SUCCESS;

    /* ==================================================== */
    /* Print to screen                                      */
    /* ==================================================== */

    std::cout << " ------------------------+------------------------------------------------------ " << std::endl;
    std::cout << " %%% SIMULATION MENU %%% :"                                                        << std::endl;
    std::cout << " ------------------------+------------------------------------------------------ " << std::endl;
    std::cout << " Parameter sweep?        : " << Input_Opt.SIMULATION_PARAMETER_SWEEP               << std::endl;
    std::cout << "  => Monte Carlo?        : " << Input_Opt.SIMULATION_MONTECARLO                    << std::endl;
    std::cout << "   => Number of runs     : " << Input_Opt.SIMULATION_MCRUNS                        << std::endl;
    std::cout << " Output folder           : " << Input_Opt.SIMULATION_OUTPUT_FOLDER                 << std::endl;
    std::cout << "  => Overwrite? if exists: " << Input_Opt.SIMULATION_OVERWRITE                     << std::endl;
    std::cout << " Run directory           : " << Input_Opt.SIMULATION_RUN_DIRECTORY                 << std::endl;
    std::cout << " Use threaded FFT?       : " << Input_Opt.SIMULATION_THREADED_FFT               << std::endl;
    std::cout << " Use FFTW WISDOM?        : " << Input_Opt.SIMULATION_USE_FFTW_WISDOM               << std::endl;
    std::cout << " => Dir w/ w permission  : " << Input_Opt.SIMULATION_DIRECTORY_W_WRITE_PERMISSION  << std::endl;
    std::cout << " Input backgrd condition : " << Input_Opt.SIMULATION_INPUT_BACKG_COND              << std::endl;
    std::cout << " Save Forward results    : " << Input_Opt.SIMULATION_SAVE_FORWARD                  << std::endl;
    std::cout << "  => netCDF file name    : " << Input_Opt.SIMULATION_FORWARD_FILENAME              << std::endl;
    std::cout << " Turn on adjoint optim.  : " << Input_Opt.SIMULATION_ADJOINT                       << std::endl;
    std::cout << "  => netCDF file name    : " << Input_Opt.SIMULATION_ADJOINT_FILENAME              << std::endl;
    std::cout << " Run box model           : " << Input_Opt.SIMULATION_BOXMODEL                      << std::endl;
    std::cout << "  => netCDF file name    : " << Input_Opt.SIMULATION_BOX_FILENAME                  << std::endl;
    std::cout << " ------------------------+------------------------------------------------------ " << std::endl;

} /* End of Read_Simulation_Menu */

void Read_Parameters( OptInput &Input_Opt, bool &RC )
{

    /* DESCRIPTION: Function Read\_Parameters reads the PARAMETER SWEEP
     * section of the APCEMM input file.
     *
     * If Parameter Sweep? is turned off in SIMULATION MENU, only the
     * first element of each parameter is considered */

    /* INPUT/OUTPUT PARAMETERS:
     * - Input_Opt: Input options
     * - RC: Success or failure
     */

    /* Read until all lines from the menu are read */

    std::vector<std::string> tokens;
    std::string subline;
    std::string variable;
    std::string unit;
    unsigned first, last;
    std::size_t found;
    double value;

    bool FIRST_PARAM_SWEEP = 1;

    /* Skip menu header lines */
    getline( inputFile, line, '\n' );

    /* Variable ranges can be defined in two ways:
     *  - min:step:max
     *  - val1 val2 val3 ... */

    /* However, for MC simulations the ranges should be defined as
     * min max or min:max */

    /* Indices for parameter input file */
    int iPlumeProcess = -1;
    int iTemperature  = -1;
    int iRHw          = -1;
    int iHDiff        = -1;
    int iVDiff        = -1;
    int iShear        = -1;
    int iLON          = -1;
    int iLAT          = -1;
    int iPressure     = -1;
    int iEDay         = -1;
    int iETime        = -1;
    int ibNOx         = -1;
    int ibHNO3        = -1;
    int ibO3          = -1;
    int ibCO          = -1;
    int ibCH4         = -1;
    int ibSO2         = -1;
    int iEINOx        = -1;
    int iNOxFlow      = -1;
    int iEICO         = -1;
    int iEIUHC        = -1;
    int iEISO2        = -1;
    int iSO2toSO4     = -1;
    int iEISoot       = -1;
    int iEISootRad    = -1;
    int itotFuelFlow  = -1;
    int iaircraftMass = -1;
    int iflightSpeed  = -1;
    int inumEngines   = -1;
    int iwingspan     = -1;
    int icoreExitTemp = -1;
    int ibypassArea   = -1;

    /* ==================================================== */
    /* Read parameters from file?                           */
    /* ==================================================== */

    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    if ( ( strcmp(tokens[0].c_str(), "T" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "t" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "TRUE" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "true" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "True" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "YES" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Y" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "y" )    == 0 ) )
        Input_Opt.PARAMETER_FILEINPUT = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "f" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "FALSE" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "false" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "False" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "NO" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "No" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "no" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "N" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "n" )     == 0 ) )
        Input_Opt.PARAMETER_FILEINPUT = 0;
    else {
        std::cout << " Wrong input for: " << "Read parameters from file?" << std::endl;
        exit(1);
    }

    /* ==================================================== */
    /* File name                                            */
    /* ==================================================== */

    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    Vector_2D fileInput( 26, Vector_1D( 0, 0.0E+00 ) );

    if ( Input_Opt.PARAMETER_FILEINPUT == 1 ) {
        std::ifstream f(tokens[0].c_str());
        if ( f.good() )
            Input_Opt.PARAMETER_FILENAME = tokens[0];
        else {
            std::cout << " Could not find parameter input file." << std::endl;
            std::string fullPath = "";
            if ( const char* simDir = std::getenv("APCEMM_runDir") )
                fullPath += simDir;
            else {
                const char* currDir = std::getenv("PWD");
                fullPath += currDir;
            }
            fullPath += tokens[0].c_str();
            std::cout << " The file does not exist. Path: " << fullPath << std::endl;
            exit(-1);
        }
        std::string currLine;
        if ( f.is_open() ) {

            /* Read header line and find fields */
            getline( f, currLine, '\n' );
            tokens = Split_Line( currLine, COMMA );

            std::cout << "\n Reading parameter input file..." << std::endl;
            for ( unsigned int i = 0; i < tokens.size(); i++ ) {
                std::transform(tokens[i].begin(), tokens[i].end(), tokens[i].begin(), ::tolower);

                std::cout << "     ";
                if ( strcmp( tokens[i].c_str(), "plume process" ) == 0 ) {
                    std::cout << std::left << std::setw(40) << " Found Plume process field at index " << std::setw(4) << i;
                    std::cout << std::right << std::setw(12) << "  (" << tokens[i].c_str() << ")" << std::endl;
                    iPlumeProcess = i;
                } else if ( ( strcmp( tokens[i].c_str(), "temperature" ) == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "temp" )        == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "t" )           == 0 ) ) {
                    std::cout << std::left << std::setw(40) << " Found Temperature field at index " << std::setw(4) << i;
                    std::cout << std::right << std::setw(12) << "  (" << tokens[i].c_str() << ")" << std::endl;
                    iTemperature = i;
                } else if ( ( strcmp( tokens[i].c_str(), "relative humidity" ) == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "rhw" )      == 0 ) ) {
                    std::cout << std::left << std::setw(40) << " Found Relative humidity field at index " << std::setw(4) << i;
                    std::cout << std::right << std::setw(12) << "  (" << tokens[i].c_str() << ")" << std::endl;
                    iRHw = i;
                } else if ( ( strcmp( tokens[i].c_str(), "horizontal diffusion" ) == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "hdiff" )      == 0 ) ) {
                    std::cout << std::left << std::setw(40) << " Found Horizontal diffusion field at index " << std::setw(4) << i;
                    std::cout << std::right << std::setw(12) << "  (" << tokens[i].c_str() << ")" << std::endl;
                    iHDiff = i;
                } else if ( ( strcmp( tokens[i].c_str(), "vertical diffusion" ) == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "vdiff" )      == 0 ) ) {
                    std::cout << std::left << std::setw(40) << " Found Vertical diffusion field at index " << std::setw(4) << i;
                    std::cout << std::right << std::setw(12) << "  (" << tokens[i].c_str() << ")" << std::endl;
                    iVDiff = i;
                } else if ( strcmp( tokens[i].c_str(), "shear" ) == 0 ) {
                    std::cout << std::left << std::setw(40) << " Found Shear field at index " << std::setw(4) << i;
                    std::cout << std::right << std::setw(12) << "  (" << tokens[i].c_str() << ")" << std::endl;
                    iShear = i;
                } else if ( ( strcmp( tokens[i].c_str(), "latitude" ) == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "lat" )      == 0 ) ) {
                    std::cout << std::left << std::setw(40) << " Found Latitude field at index " << std::setw(4) << i;
                    std::cout << std::right << std::setw(12) << "  (" << tokens[i].c_str() << ")" << std::endl;
                    iLAT = i;
                } else if ( ( strcmp( tokens[i].c_str(), "longitude" ) == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "long" )      == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "lon" )       == 0 ) ) {
                    std::cout << std::left << std::setw(40) << " Found Longitude field at index " << std::setw(4) << i;
                    std::cout << std::right << std::setw(12) << "  (" << tokens[i].c_str() << ")" << std::endl;
                    iLON = i;
                } else if ( ( strcmp( tokens[i].c_str(), "pressure" ) == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "press" )    == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "p" )        == 0 ) ) {
                    std::cout << std::left << std::setw(40) << " Found Pressure field at index " << std::setw(4) << i;
                    std::cout << std::right << std::setw(12) << "  (" << tokens[i].c_str() << ")" << std::endl;
                    iPressure = i;
                } else if ( ( strcmp( tokens[i].c_str(), "emission day" ) == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "e.day" )        == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "eday" )         == 0 ) ) {
                    std::cout << std::left << std::setw(40) << " Found Emission day field at index " << std::setw(4) << i;
                    std::cout << std::right << std::setw(12) << "  (" << tokens[i].c_str() << ")" << std::endl;
                    iEDay = i;
                } else if ( ( strcmp( tokens[i].c_str(), "emission time" ) == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "e.time" )        == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "etime" )         == 0 ) ) {
                    std::cout << std::left << std::setw(40) << " Found Emission time field at index " << std::setw(4) << i;
                    std::cout << std::right << std::setw(12) << "  (" << tokens[i].c_str() << ")" << std::endl;
                    iETime = i;
                } else if ( ( strcmp( tokens[i].c_str(), "background nox" ) == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "b.nox" )          == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "bnox" )           == 0 ) ) {
                    std::cout << std::left << std::setw(40) << " Found Background NOx field at index " << std::setw(4) << i;
                    std::cout << std::right << std::setw(12) << "  (" << tokens[i].c_str() << ")" << std::endl;
                    ibNOx = i;
                } else if ( ( strcmp( tokens[i].c_str(), "background hno3" ) == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "b.hno3" )          == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "bhno3" )           == 0 ) ) {
                    std::cout << std::left << std::setw(40) << " Found Background HNO3 field at index " << std::setw(4) << i;
                    std::cout << std::right << std::setw(12) << "  (" << tokens[i].c_str() << ")" << std::endl;
                    ibHNO3 = i;
                } else if ( ( strcmp( tokens[i].c_str(), "background o3" ) == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "b.o3" )          == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "bo3" )           == 0 ) ) {
                    std::cout << std::left << std::setw(40) << " Found Background O3 field at index " << std::setw(4) << i;
                    std::cout << std::right << std::setw(12) << "  (" << tokens[i].c_str() << ")" << std::endl;
                    ibO3 = i;
                } else if ( ( strcmp( tokens[i].c_str(), "background co" ) == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "b.co" )          == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "bco" )           == 0 ) ) {
                    std::cout << std::left << std::setw(40) << " Found Background CO field at index " << std::setw(4) << i;
                    std::cout << std::right << std::setw(12) << "  (" << tokens[i].c_str() << ")" << std::endl;
                    ibCO = i;
                } else if ( ( strcmp( tokens[i].c_str(), "background ch4" ) == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "b.ch4" )          == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "bch4" )           == 0 ) ) {
                    std::cout << std::left << std::setw(40) << " Found Background CH4 field at index " << std::setw(4) << i;
                    std::cout << std::right << std::setw(12) << "  (" << tokens[i].c_str() << ")" << std::endl;
                    ibCH4 = i;
                } else if ( ( strcmp( tokens[i].c_str(), "background so2" ) == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "b.so2" )          == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "bso2" )           == 0 ) ) {
                    std::cout << std::left << std::setw(40) << " Found Background SO2 field at index " << std::setw(4) << i;
                    std::cout << std::right << std::setw(12) << "  (" << tokens[i].c_str() << ")" << std::endl;
                    ibSO2 = i;
                } else if ( ( strcmp( tokens[i].c_str(), "emission index nox" ) == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "ei.nox" )             == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "einox" )              == 0 ) ) {
                    std::cout << std::left << std::setw(40) << " Found Emission index NOx field at index " << std::setw(4) << i;
                    std::cout << std::right << std::setw(12) << "  (" << tokens[i].c_str() << ")" << std::endl;
                    iEINOx = i;
                } else if ( ( strcmp( tokens[i].c_str(), "nox flow" ) == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "noxflow" )  == 0 ) ) {
                    std::cout << std::left << std::setw(40) << " Found NOx flow field at index " << std::setw(4) << i;
                    std::cout << std::right << std::setw(12) << "  (" << tokens[i].c_str() << ")" << std::endl;
                    iNOxFlow = i;
                } else if ( ( strcmp( tokens[i].c_str(), "emission index co" ) == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "ei.co" )             == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "eico" )              == 0 ) ) {
                    std::cout << std::left << std::setw(40) << " Found Emission index CO field at index " << std::setw(4) << i;
                    std::cout << std::right << std::setw(12) << "  (" << tokens[i].c_str() << ")" << std::endl;
                    iEICO = i;
                } else if ( ( strcmp( tokens[i].c_str(), "emission index uhc" ) == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "ei.uhc" )             == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "eiuhc" )              == 0 ) ) {
                    std::cout << std::left << std::setw(40) << " Found Emission index UHC field at index " << std::setw(4) << i;
                    std::cout << std::right << std::setw(12) << "  (" << tokens[i].c_str() << ")" << std::endl;
                    iEIUHC = i;
                } else if ( ( strcmp( tokens[i].c_str(), "emission index SO2" ) == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "ei.so2" )             == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "eiso2" )              == 0 ) ) {
                    std::cout << std::left << std::setw(40) << " Found Emission index SO2 field at index " << std::setw(4) << i;
                    std::cout << std::right << std::setw(12) << "  (" << tokens[i].c_str() << ")" << std::endl;
                    iEISO2 = i;
                } else if ( ( strcmp( tokens[i].c_str(), "so2 to so4 conversion" ) == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "so2 to so4 conv" )       == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "so2 to so4" )            == 0 ) ) {
                    std::cout << std::left << std::setw(40) << " Found SO2 to SO4 conversion field at index " << std::setw(4) << i;
                    std::cout << std::right << std::setw(12) << "  (" << tokens[i].c_str() << ")" << std::endl;
                    iSO2toSO4 = i;
                } else if ( ( strcmp( tokens[i].c_str(), "emission index soot" ) == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "ei.soot" )             == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "eisoot" )              == 0 ) ) {
                    std::cout << std::left << std::setw(40) << " Found Emission index Soot field at index " << std::setw(4) << i;
                    std::cout << std::right << std::setw(12) << "  (" << tokens[i].c_str() << ")" << std::endl;
                    iEISoot = i;
                } else if ( ( strcmp( tokens[i].c_str(), "emission index sootrad" ) == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "ei.sootrad" )             == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "eisootrad" )              == 0 ) ) {
                    std::cout << std::left << std::setw(40) << " Found Emission index SootRad field at index " << std::setw(4) << i;
                    std::cout << std::right << std::setw(12) << "  (" << tokens[i].c_str() << ")" << std::endl;
                    iEISootRad = i;
                } else if ( ( strcmp( tokens[i].c_str(), "total fuel flow" ) == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "fuel flow" )       == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "fuel" )            == 0 ) ) {
                    std::cout << std::left << std::setw(40) << " Found Total fuel flow field at index " << std::setw(4) << i;
                    std::cout << std::right << std::setw(12) << "  (" << tokens[i].c_str() << ")" << std::endl;
                    itotFuelFlow = i;
                } else if ( ( strcmp( tokens[i].c_str(), "aircraft mass" ) == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "ac mass" )       == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "acmass" )        == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "amass" )         == 0 ) ) {
                    std::cout << std::left << std::setw(40) << " Found Aircraft mass field at index " << std::setw(4) << i;
                    std::cout << std::right << std::setw(12) << "  (" << tokens[i].c_str() << ")" << std::endl;
                    iaircraftMass = i;                                          
                } else if ( ( strcmp( tokens[i].c_str(), "flight speed" ) == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "flt speed" )       == 0 ) || \              
                            ( strcmp( tokens[i].c_str(), "fltspeed" )        == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "fspeed" )         == 0 ) ) {               
                    std::cout << std::left << std::setw(40) << " Found Aircraft flight speed field at index " << std::setw(4) << i;                                   
                    std::cout << std::right << std::setw(12) << "  (" << tokens[i].c_str() << ")" << std::endl;
                    iflightSpeed = i;
                } else if ( ( strcmp( tokens[i].c_str(), "number of engines" ) == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "num. of engines" )       == 0 ) || \              
                            ( strcmp( tokens[i].c_str(), "num. eng." )        == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "numeng" )         == 0 ) ) {               
                    std::cout << std::left << std::setw(40) << " Found Aircraft number of engines field at index " << std::setw(4) << i;                                   
                    std::cout << std::right << std::setw(12) << "  (" << tokens[i].c_str() << ")" << std::endl;
                    inumEngines = i;                                          
                } else if ( ( strcmp( tokens[i].c_str(), "wing span" ) == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "wingspan" )         == 0 ) ) {               
                    std::cout << std::left << std::setw(40) << " Found Aircraft wingspan field at index " << std::setw(4) << i;                                   
                    std::cout << std::right << std::setw(12) << "  (" << tokens[i].c_str() << ")" << std::endl;
                    iwingspan = i;                                          
                } else if ( ( strcmp( tokens[i].c_str(), "core exit temperature" ) == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "core exit temp." )       == 0 ) || \              
                            ( strcmp( tokens[i].c_str(), "core exit T" )        == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "core T" )         == 0 ) ) {               
                    std::cout << std::left << std::setw(40) << " Found core exit temperature field at index " << std::setw(4) << i;                                   
                    std::cout << std::right << std::setw(12) << "  (" << tokens[i].c_str() << ")" << std::endl;
                    icoreExitTemp = i;
                } else if ( ( strcmp( tokens[i].c_str(), "exit bypass area" ) == 0 ) || \
                            ( strcmp( tokens[i].c_str(), "bypass area" )       == 0 ) || \              
                            ( strcmp( tokens[i].c_str(), "area" )         == 0 ) ) {               
                    std::cout << std::left << std::setw(40) << " Found bypass exit area field at index " << std::setw(4) << i;                                   
                    std::cout << std::right << std::setw(12) << "  (" << tokens[i].c_str() << ")" << std::endl;
                    ibypassArea = i;
                } else 
                    std::cout << " Ignoring field: " << tokens[i] << std::endl;
            }

            /* Read-in values */
            while ( 1 ) {

                /* Read in current line */
                f >> currLine;

                /* If we have reached EOF, then break; */
                if ( f.eof() ) break;

                tokens = Split_Line( currLine, COMMA );

                Input_Opt.PARAMETER_FILECASES += 1;
                for ( unsigned int i = 0; i < tokens.size(); i++ ) {
                    try {
                        if ( iPlumeProcess == i )
                            fileInput[0].push_back( std::stod( tokens[i] ) );
                        else if ( iTemperature == i )
                            fileInput[1].push_back( std::stod( tokens[i] ) );
                        else if ( iRHw == i )
                            fileInput[2].push_back( std::stod( tokens[i] ) );
                        else if ( iHDiff == i )
                            fileInput[3].push_back( std::stod( tokens[i] ) );
                        else if ( iVDiff == i )
                            fileInput[4].push_back( std::stod( tokens[i] ) );
                        else if ( iShear == i )
                            fileInput[5].push_back( std::stod( tokens[i] ) );
                        else if ( iLON == i )
                            fileInput[6].push_back( std::stod( tokens[i] ) );
                        else if ( iLAT == i )
                            fileInput[7].push_back( std::stod( tokens[i] ) );
                        else if ( iPressure == i )
                            fileInput[8].push_back( std::stod( tokens[i] ) );
                        else if ( iEDay == i )
                            fileInput[9].push_back( std::stoi( tokens[i] ) );
                        else if ( iETime == i )
                            fileInput[10].push_back( std::stod( tokens[i] ) );
                        else if ( iEINOx == i )
                            fileInput[11].push_back( std::stod( tokens[i] ) );
                        else if ( iNOxFlow == i )
                            fileInput[11].push_back( std::stod( tokens[i] ) );
                        else if ( iEICO == i )
                            fileInput[12].push_back( std::stod( tokens[i] ) );
                        else if ( iEIUHC == i )
                            fileInput[13].push_back( std::stod( tokens[i] ) );
                        else if ( iEISO2 == i )
                            fileInput[14].push_back( std::stod( tokens[i] ) );
                        else if ( iSO2toSO4 == i )
                            fileInput[15].push_back( std::stod( tokens[i] ) );
                        else if ( iEISoot == i )
                            fileInput[16].push_back( std::stod( tokens[i] ) );
                        else if ( iEISootRad == i )
                            fileInput[17].push_back( std::stod( tokens[i] ) );
                        else if ( itotFuelFlow == i )
                            fileInput[18].push_back( std::stod( tokens[i] ) );
                        else if ( iaircraftMass == i )
                            fileInput[19].push_back( std::stod( tokens[i] ) );
                        else if ( ibNOx == i )
                            fileInput[20].push_back( std::stod( tokens[i] ) );
                        else if ( ibHNO3 == i )
                            fileInput[21].push_back( std::stod( tokens[i] ) );
                        else if ( ibO3 == i )
                            fileInput[22].push_back( std::stod( tokens[i] ) );
                        else if ( ibCO == i )
                            fileInput[23].push_back( std::stod( tokens[i] ) );
                        else if ( ibCH4 == i )
                            fileInput[24].push_back( std::stod( tokens[i] ) );
                        else if ( ibSO2 == i )
                            fileInput[25].push_back( std::stod( tokens[i] ) );
                        else if ( iflightSpeed == i )
                            fileInput[26].push_back( std::stod( tokens[i] ) );
                        else if ( inumEngines == i )
                            fileInput[27].push_back( std::stod( tokens[i] ) );
                        else if ( iwingspan == i )
                            fileInput[28].push_back( std::stod( tokens[i] ) );
                        else if ( icoreExitTemp == i )
                            fileInput[29].push_back( std::stod( tokens[i] ) );
                        else if ( ibypassArea == i )
                            fileInput[30].push_back( std::stod( tokens[i] ) );
                    } catch( std::exception& e ) {
                        std::cout << " Could not convert string '" << tokens[i] << "' to double/int for index " << i << " in parameter input file." << std::endl;
                        exit(-1);
                    }
                }

//                // FOR DEBUG PURPOSES!
//                for ( unsigned int i = 0; i < fileInput.size(); i++ ) {
//                    std::cout << "Looking at line " << i << std::endl;
//                    for ( unsigned int j = 0; j < fileInput[i].size(); j++ ) {
//                        std::cout << fileInput[i][j] << ", ";
//                    }
//                    std::cout << std::endl;
//                }

            }
            f.close();
            std::cout << std::endl;
        }
    }

    getline( inputFile, line, '\n' );
    getline( inputFile, line, '\n' );
    getline( inputFile, line, '\n' );

    /* ==================================================== */
    /* Plume processing time                                */
    /* ==================================================== */

    /* Variable */
    variable = "Plume process.";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Get line past the delimiter */
    subline = line.substr(FIRSTCOL);
    /* Look for colon */
    found = subline.find( COLON );
    if ( found != std::string::npos ) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

        if ( !Input_Opt.SIMULATION_MONTECARLO ) {
            if (tokens.size() != 3) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: " << std::endl;
                std::cout << "   --> begin:step:end" << std::endl;
                std::cout << "       or" << std::endl;
                std::cout << "   --> val1 val2 val3 ..." << std::endl;
                exit(1);
            }
            if ( ( std::stod(tokens[2]) <  std::stod(tokens[0]) ) || \
                 ( std::stod(tokens[1]) <= 0.0E+00 )              || \
                 ( std::stod(tokens[1]) >  ( std::stod(tokens[2]) - std::stod(tokens[0]) ) ) ) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: begin:step:end with begin < end, 0 < step < end - begin" << std::endl;
                exit(1);
            }
        }
        Input_Opt.PARAMETER_PLUMEPROCESS_RANGE = 1;
    }
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for " << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_PLUMEPROCESS_RANGE = 0;
    }

    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) && ( FIRST_PARAM_SWEEP ) ) {
        std::cout << " Multiple cases need to be run through APCEMM while the 'Parameter sweep?' argument is turned off!" << std::endl;
        std::cout << " These cases will be run in serial!" << std::endl;
        std::cout << " To enable the processing of multiple cases simultaneously, turn on the 'Parameter sweep' input!" << std::endl;
        FIRST_PARAM_SWEEP = 0;
    }

    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( tokens.size() > 2 ) {
            std::cout << " Wrong input for " << variable << " when MC is turned on!" << std::endl;
            std::cout << " Expected format is min max or min:max representing the range of possible values" << std::endl;
            exit(1);
        } else if ( tokens.size() == 2 ) {
            Input_Opt.PARAMETER_PLUMEPROCESS_RANGE = 1;
            sort(tokens.begin(), tokens.end());
        } else if ( tokens.size() == 1 )
            Input_Opt.PARAMETER_PLUMEPROCESS_RANGE = 0;
    }

    /* Find unit in between "[" and "]" */
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_PLUMEPROCESS_UNIT.assign( unit );

    if ( ( Input_Opt.PARAMETER_FILEINPUT ) && ( iPlumeProcess != -1 ) ) {
        /* Found plume process in file */
        Input_Opt.PARAMETER_PLUMEPROCESS_RANGE = 0;
        if ( fileInput[0].size() > 0 ) {
            for ( unsigned int i = 0; i < fileInput[0].size(); i++ ) {
                if ( fileInput[0][i] > 0.0E+00 )
                    Input_Opt.PARAMETER_PLUMEPROCESS.push_back( fileInput[0][i] );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            }
        } else {
            try {
                value = std::stod( tokens[0] );
                if ( value > 0.0E+00 )
                    Input_Opt.PARAMETER_PLUMEPROCESS.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[0] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    } else {
        /* Store in values for variable */
        for ( unsigned int i = 0; i < tokens.size(); i++ ) {
            try {
                value = std::stod( tokens[i] );
                if ( value > 0.0E+00 )
                    Input_Opt.PARAMETER_PLUMEPROCESS.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    }

    /* Skip header line */
    getline( inputFile, line, '\n' );

    /* ==================================================== */
    /* Temperature                                          */
    /* ==================================================== */

    /* Variable */
    variable = "Temperature";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Get line past the delimiter */
    subline = line.substr(FIRSTCOL);
    /* Look for colon */
    found = subline.find( COLON );
    if ( found != std::string::npos ) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

        if ( !Input_Opt.SIMULATION_MONTECARLO ) {
            if (tokens.size() != 3) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: " << std::endl;
                std::cout << "   --> begin:step:end" << std::endl;
                std::cout << "       or" << std::endl;
                std::cout << "   --> val1 val2 val3 ..." << std::endl;
                exit(1);
            }
            if ( ( std::stod(tokens[2]) <  std::stod(tokens[0]) ) || \
                 ( std::stod(tokens[1]) <= 0.0E+00 )              || \
                 ( std::stod(tokens[1]) >  ( std::stod(tokens[2]) - std::stod(tokens[0]) ) ) ) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: begin:step:end with begin < end, 0 < step < end - begin" << std::endl;
                exit(1);
            }
        }
        Input_Opt.PARAMETER_TEMPERATURE_RANGE = 1;
    }
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for " << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_TEMPERATURE_RANGE = 0;
    }

    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) && ( FIRST_PARAM_SWEEP ) ) {
        std::cout << " Multiple cases need to be run through APCEMM while the 'Parameter sweep?' argument is turned off!" << std::endl;
        std::cout << " These cases will be run in serial!" << std::endl;
        std::cout << " To enable the processing of multiple cases simultaneously, turn on the 'Parameter sweep' input!" << std::endl;
        FIRST_PARAM_SWEEP = 0;
    }

    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( tokens.size() > 2 ) {
            std::cout << " Wrong input for " << variable << " when MC is turned on!" << std::endl;
            std::cout << " Expected format is min max or min:max representing the range of possible values" << std::endl;
            exit(1);
        } else if ( tokens.size() == 2 ) {
            Input_Opt.PARAMETER_TEMPERATURE_RANGE = 1;
            sort(tokens.begin(), tokens.end());
        } else if ( tokens.size() == 1 )
            Input_Opt.PARAMETER_TEMPERATURE_RANGE = 0;
    }

    /* Find unit in between "[" and "]" */
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_TEMPERATURE_UNIT.assign( unit );

    if ( ( Input_Opt.PARAMETER_FILEINPUT ) && ( iTemperature != -1 ) ) {
        /* Found temperature in parameter input file */
        Input_Opt.PARAMETER_TEMPERATURE_RANGE = 0;
        if ( fileInput[1].size() > 0 ) {
            for ( unsigned int i = 0; i < fileInput[1].size(); i++ ) {
                if ( fileInput[1][i] > 0.0E+00 )
                    Input_Opt.PARAMETER_TEMPERATURE.push_back( fileInput[1][i] );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            }
        } else {
            try {
                value = std::stod( tokens[0] );
                if ( value > 0.0E+00 )
                    Input_Opt.PARAMETER_TEMPERATURE.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[0] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    } else {
        /* Store in values for variable */
        for ( unsigned int i = 0; i < tokens.size(); i++ ) {
            try {
                value = std::stod( tokens[i] );
                if ( value > 0.0E+00 )
                    Input_Opt.PARAMETER_TEMPERATURE.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    }

    /* ==================================================== */
    /* Relative humidity                                    */
    /* ==================================================== */

    /* Variable */
    variable = "Relative humidity";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Get line past the delimiter */
    subline = line.substr(FIRSTCOL);
    /* Look for colon */
    found = subline.find( COLON );
    if ( found != std::string::npos ) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

        if ( !Input_Opt.SIMULATION_MONTECARLO ) {
            if (tokens.size() != 3) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: " << std::endl;
                std::cout << "   --> begin:step:end" << std::endl;
                std::cout << "       or" << std::endl;
                std::cout << "   --> val1 val2 val3 ..." << std::endl;
                exit(1);
            }
            if ( ( std::stod(tokens[2]) <  std::stod(tokens[0]) ) || \
                 ( std::stod(tokens[1]) <= 0.0E+00 )              || \
                 ( std::stod(tokens[1]) >  ( std::stod(tokens[2]) - std::stod(tokens[0]) ) ) ) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: begin:step:end with begin < end, 0 < step < end - begin" << std::endl;
                exit(1);
            }
        }
        Input_Opt.PARAMETER_RHW_RANGE = 1;
    }
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for " << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_RHW_RANGE = 0;
    }

    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) && ( FIRST_PARAM_SWEEP ) ) {
        std::cout << " Multiple cases need to be run through APCEMM while the 'Parameter sweep?' argument is turned off!" << std::endl;
        std::cout << " These cases will be run in serial!" << std::endl;
        std::cout << " To enable the processing of multiple cases simultaneously, turn on the 'Parameter sweep' input!" << std::endl;
        FIRST_PARAM_SWEEP = 0;
    }

    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( tokens.size() > 2 ) {
            std::cout << " Wrong input for " << variable << " when MC is turned on!" << std::endl;
            std::cout << " Expected format is min max or min:max representing the range of possible values" << std::endl;
            exit(1);
        } else if ( tokens.size() == 2 ) {
            Input_Opt.PARAMETER_RHW_RANGE = 1;
            sort(tokens.begin(), tokens.end());
        } else if ( tokens.size() == 1 )
            Input_Opt.PARAMETER_RHW_RANGE = 0;
    }

    /* Find unit in between "[" and "]" */
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_RHW_UNIT.assign( unit );

    if ( ( Input_Opt.PARAMETER_FILEINPUT ) && ( iRHw != -1 ) ) {
        /* Found relative humidity in parameter input file */
        Input_Opt.PARAMETER_RHW_RANGE = 0;
        if ( fileInput[2].size() > 0 ) {
            for ( unsigned int i = 0; i < fileInput[2].size(); i++ ) {
                if ( fileInput[2][i] >= 0.0E+00 )
                    Input_Opt.PARAMETER_RHW.push_back( fileInput[2][i] );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            }
        } else {
            try {
                value = std::stod( tokens[0] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_RHW.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[0] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    } else {
        /* Store in values for variable */
        for ( unsigned int i = 0; i < tokens.size(); i++ ) {
            try {
                value = std::stod( tokens[i] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_RHW.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    }

    /* ==================================================== */
    /* Horiz. diffusion parameter                           */
    /* ==================================================== */

    /* Variable */
    variable = "Horiz. diffusion parameter";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Get line past the delimiter */
    subline = line.substr(FIRSTCOL);
    /* Look for colon */
    found = subline.find( COLON );
    if ( found != std::string::npos ) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

        if ( !Input_Opt.SIMULATION_MONTECARLO ) {
            if (tokens.size() != 3) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: " << std::endl;
                std::cout << "   --> begin:step:end" << std::endl;
                std::cout << "       or" << std::endl;
                std::cout << "   --> val1 val2 val3 ..." << std::endl;
                exit(1);
            }
            if ( ( std::stod(tokens[2]) <  std::stod(tokens[0]) ) || \
                 ( std::stod(tokens[1]) <= 0.0E+00 )              || \
                 ( std::stod(tokens[1]) >  ( std::stod(tokens[2]) - std::stod(tokens[0]) ) ) ) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: begin:step:end with begin < end, 0 < step < end - begin" << std::endl;
                exit(1);
            }
        }
        Input_Opt.PARAMETER_DH_RANGE = 1;
    }
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for " << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_DH_RANGE = 0;
    }

    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) && ( FIRST_PARAM_SWEEP ) ) {
        std::cout << " Multiple cases need to be run through APCEMM while the 'Parameter sweep?' argument is turned off!" << std::endl;
        std::cout << " These cases will be run in serial!" << std::endl;
        std::cout << " To enable the processing of multiple cases simultaneously, turn on the 'Parameter sweep' input!" << std::endl;
        FIRST_PARAM_SWEEP = 0;
    }

    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( tokens.size() > 2 ) {
            std::cout << " Wrong input for " << variable << " when MC is turned on!" << std::endl;
            std::cout << " Expected format is min max or min:max representing the range of possible values" << std::endl;
            exit(1);
        } else if ( tokens.size() == 2 ) {
            Input_Opt.PARAMETER_DH_RANGE = 1;
            sort(tokens.begin(), tokens.end());
        } else if ( tokens.size() == 1 )
            Input_Opt.PARAMETER_DH_RANGE = 0;
    }

    /* Find unit in between "[" and "]" */
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_DH_UNIT.assign( unit );

    if ( ( Input_Opt.PARAMETER_FILEINPUT ) && ( iHDiff != -1 ) ) {
        /* Found horizontal diffusion coefficient in parameter input file */
        Input_Opt.PARAMETER_DH_RANGE = 0;
        if ( fileInput[3].size() > 0 ) {
            for ( unsigned int i = 0; i < fileInput[3].size(); i++ ) {
                if ( fileInput[3][i] >= 0.0E+00 )
                    Input_Opt.PARAMETER_DH.push_back( fileInput[3][i] );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            }
        } else {
            try {
                value = std::stod( tokens[0] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_DH.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[0] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    } else {
        /* Store in values for variable */
        for ( unsigned int i = 0; i < tokens.size(); i++ ) {
            try {
                value = std::stod( tokens[i] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_DH.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    }

    /* ==================================================== */
    /* Verti. diffusion parameter                           */
    /* ==================================================== */

    /* Variable */
    variable = "Verti. diffusion parameter";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Get line past the delimiter */
    subline = line.substr(FIRSTCOL);
    /* Look for colon */
    found = subline.find( COLON );
    if ( found != std::string::npos ) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

        if ( !Input_Opt.SIMULATION_MONTECARLO ) {
            if (tokens.size() != 3) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: " << std::endl;
                std::cout << "   --> begin:step:end" << std::endl;
                std::cout << "       or" << std::endl;
                std::cout << "   --> val1 val2 val3 ..." << std::endl;
                exit(1);
            }
            if ( ( std::stod(tokens[2]) <  std::stod(tokens[0]) ) || \
                 ( std::stod(tokens[1]) <= 0.0E+00 )              || \
                 ( std::stod(tokens[1]) >  ( std::stod(tokens[2]) - std::stod(tokens[0]) ) ) ) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: begin:step:end with begin < end, 0 < step < end - begin" << std::endl;
                exit(1);
            }
        }
        Input_Opt.PARAMETER_DV_RANGE = 1;
    }
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for " << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_DV_RANGE = 0;
    }

    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) && ( FIRST_PARAM_SWEEP ) ) {
        std::cout << " Multiple cases need to be run through APCEMM while the 'Parameter sweep?' argument is turned off!" << std::endl;
        std::cout << " These cases will be run in serial!" << std::endl;
        std::cout << " To enable the processing of multiple cases simultaneously, turn on the 'Parameter sweep' input!" << std::endl;
        FIRST_PARAM_SWEEP = 0;
    }

    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( tokens.size() > 2 ) {
            std::cout << " Wrong input for " << variable << " when MC is turned on!" << std::endl;
            std::cout << " Expected format is min max or min:max representing the range of possible values" << std::endl;
            exit(1);
        } else if ( tokens.size() == 2 ) {
            Input_Opt.PARAMETER_DV_RANGE = 1;
            sort(tokens.begin(), tokens.end());
        } else if ( tokens.size() == 1 )
            Input_Opt.PARAMETER_DV_RANGE = 0;
    }

    /* Find unit in between "[" and "]" */
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_DV_UNIT.assign( unit );

    if ( ( Input_Opt.PARAMETER_FILEINPUT ) && ( iVDiff != -1 ) ) {
        /* Found vertical diffusion coefficient in parameter input file */
        Input_Opt.PARAMETER_DV_RANGE = 0;
        if ( fileInput[4].size() > 0 ) {
            for ( unsigned int i = 0; i < fileInput[4].size(); i++ ) {
                if ( fileInput[4][i] >= 0.0E+00 )
                    Input_Opt.PARAMETER_DV.push_back( fileInput[4][i] );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            }
        } else {
            try {
                value = std::stod( tokens[0] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_DV.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[0] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    } else {
        /* Store in values for variable */
        for ( unsigned int i = 0; i < tokens.size(); i++ ) {
            try {
                value = std::stod( tokens[i] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_DV.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    }

    /* ==================================================== */
    /* Shear                                                */
    /* ==================================================== */

    /* Variable */
    variable = "Shear";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Get line past the delimiter */
    subline = line.substr(FIRSTCOL);
    /* Look for colon */
    found = subline.find( COLON );
    if ( found != std::string::npos ) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

        if ( !Input_Opt.SIMULATION_MONTECARLO ) {
            if (tokens.size() != 3) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: " << std::endl;
                std::cout << "   --> begin:step:end" << std::endl;
                std::cout << "       or" << std::endl;
                std::cout << "   --> val1 val2 val3 ..." << std::endl;
                exit(1);
            }
            if ( ( std::stod(tokens[2]) <  std::stod(tokens[0]) ) || \
                 ( std::stod(tokens[1]) <= 0.0E+00 )              || \
                 ( std::stod(tokens[1]) >  ( std::stod(tokens[2]) - std::stod(tokens[0]) ) ) ) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: begin:step:end with begin < end, 0 < step < end - begin" << std::endl;
                exit(1);
            }
        }
        Input_Opt.PARAMETER_SHEAR_RANGE = 1;
    }
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for " << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_SHEAR_RANGE = 0;
    }

    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) && ( FIRST_PARAM_SWEEP ) ) {
        std::cout << " Multiple cases need to be run through APCEMM while the 'Parameter sweep?' argument is turned off!" << std::endl;
        std::cout << " These cases will be run in serial!" << std::endl;
        std::cout << " To enable the processing of multiple cases simultaneously, turn on the 'Parameter sweep' input!" << std::endl;
        FIRST_PARAM_SWEEP = 0;
    }

    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( tokens.size() > 2 ) {
            std::cout << " Wrong input for " << variable << " when MC is turned on!" << std::endl;
            std::cout << " Expected format is min max or min:max representing the range of possible values" << std::endl;
            exit(1);
        } else if ( tokens.size() == 2 ) {
            Input_Opt.PARAMETER_SHEAR_RANGE = 1;
            sort(tokens.begin(), tokens.end());
        } else if ( tokens.size() == 1 )
            Input_Opt.PARAMETER_SHEAR_RANGE = 0;
    }

    /* Find unit in between "[" and "]" */
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_SHEAR_UNIT.assign( unit );

    if ( ( Input_Opt.PARAMETER_FILEINPUT ) && ( iShear != -1 ) ) {
        /* Found shear in parameter input file */
        Input_Opt.PARAMETER_SHEAR_RANGE = 0;
        if ( fileInput[5].size() > 0 ) {
            for ( unsigned int i = 0; i < fileInput[5].size(); i++ )
                Input_Opt.PARAMETER_SHEAR.push_back( fileInput[5][i] );
        } else {
            try {
                value = std::stod( tokens[0] );
                Input_Opt.PARAMETER_SHEAR.push_back( value );
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[0] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    } else {
        /* Store in values for variable */
        for ( unsigned int i = 0; i < tokens.size(); i++ ) {
            try {
                value = std::stod( tokens[i] );
                Input_Opt.PARAMETER_SHEAR.push_back( value );
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    }

    /* Skipping line */
    getline( inputFile, line, '\n' );

    /* ==================================================== */
    /* Longitude                                            */
    /* ==================================================== */

    /* Variable */
    variable = "Longitude";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Get line past the delimiter */
    subline = line.substr(FIRSTCOL);
    /* Look for colon */
    found = subline.find( COLON );
    if ( found != std::string::npos ) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

        if ( !Input_Opt.SIMULATION_MONTECARLO ) {
            if (tokens.size() != 3) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: " << std::endl;
                std::cout << "   --> begin:step:end" << std::endl;
                std::cout << "       or" << std::endl;
                std::cout << "   --> val1 val2 val3 ..." << std::endl;
                exit(1);
            }
            if ( ( std::stod(tokens[2]) <  std::stod(tokens[0]) ) || \
                 ( std::stod(tokens[1]) <= 0.0E+00 )              || \
                 ( std::stod(tokens[1]) >  ( std::stod(tokens[2]) - std::stod(tokens[0]) ) ) ) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: begin:step:end with begin < end, 0 < step < end - begin" << std::endl;
                exit(1);
            }
        }
        Input_Opt.PARAMETER_LONGITUDE_RANGE = 1;
    }
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for " << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_LONGITUDE_RANGE = 0;
    }

    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) && ( FIRST_PARAM_SWEEP ) ) {
        std::cout << " Multiple cases need to be run through APCEMM while the 'Parameter sweep?' argument is turned off!" << std::endl;
        std::cout << " These cases will be run in serial!" << std::endl;
        std::cout << " To enable the processing of multiple cases simultaneously, turn on the 'Parameter sweep' input!" << std::endl;
        FIRST_PARAM_SWEEP = 0;
    }

    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( tokens.size() > 2 ) {
            std::cout << " Wrong input for " << variable << " when MC is turned on!" << std::endl;
            std::cout << " Expected format is min max or min:max representing the range of possible values" << std::endl;
            exit(1);
        } else if ( tokens.size() == 2 ) {
            Input_Opt.PARAMETER_LONGITUDE_RANGE = 1;
            sort(tokens.begin(), tokens.end());
        } else if ( tokens.size() == 1 )
            Input_Opt.PARAMETER_LONGITUDE_RANGE = 0;
    }

    /* Find unit in between "[" and "]" */
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_LONGITUDE_UNIT.assign( unit );

    if ( ( Input_Opt.PARAMETER_FILEINPUT ) && ( iLON != -1 ) ) {
        /* Found longitude in parameter input file */
        Input_Opt.PARAMETER_LONGITUDE_RANGE = 0;
        if ( fileInput[6].size() > 0 ) {
            for ( unsigned int i = 0; i < fileInput[6].size(); i++ )
                Input_Opt.PARAMETER_LONGITUDE.push_back( fileInput[6][i] );
        } else {
            try {
                value = std::stod( tokens[0] );
                Input_Opt.PARAMETER_LONGITUDE.push_back( value );
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[0] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    } else {
        /* Store in values for variable */
        for ( unsigned int i = 0; i < tokens.size(); i++ ) {
            try {
                Input_Opt.PARAMETER_LONGITUDE.push_back(std::stod( tokens[i] ));
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    }

    /* ==================================================== */
    /* Latitude                                             */
    /* ==================================================== */

    /* Variable */
    variable = "Latitude";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Get line past the delimiter */
    subline = line.substr(FIRSTCOL);
    /* Look for colon */
    found = subline.find( COLON );
    if ( found != std::string::npos ) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

        if ( !Input_Opt.SIMULATION_MONTECARLO ) {
            if (tokens.size() != 3) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: " << std::endl;
                std::cout << "   --> begin:step:end" << std::endl;
                std::cout << "       or" << std::endl;
                std::cout << "   --> val1 val2 val3 ..." << std::endl;
                exit(1);
            }
            if ( ( std::stod(tokens[2]) <  std::stod(tokens[0]) ) || \
                 ( std::stod(tokens[1]) <= 0.0E+00 )              || \
                 ( std::stod(tokens[1]) >  ( std::stod(tokens[2]) - std::stod(tokens[0]) ) ) ) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: begin:step:end with begin < end, 0 < step < end - begin" << std::endl;
                exit(1);
            }
        }
        Input_Opt.PARAMETER_LATITUDE_RANGE = 1;
    }
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for " << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_LATITUDE_RANGE = 0;
    }

    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) && ( FIRST_PARAM_SWEEP ) ) {
        std::cout << " Multiple cases need to be run through APCEMM while the 'Parameter sweep?' argument is turned off!" << std::endl;
        std::cout << " These cases will be run in serial!" << std::endl;
        std::cout << " To enable the processing of multiple cases simultaneously, turn on the 'Parameter sweep' input!" << std::endl;
        FIRST_PARAM_SWEEP = 0;
    }

    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( tokens.size() > 2 ) {
            std::cout << " Wrong input for " << variable << " when MC is turned on!" << std::endl;
            std::cout << " Expected format is min max or min:max representing the range of possible values" << std::endl;
            exit(1);
        } else if ( tokens.size() == 2 ) {
            Input_Opt.PARAMETER_LATITUDE_RANGE = 1;
            sort(tokens.begin(), tokens.end());
        } else if ( tokens.size() == 1 )
            Input_Opt.PARAMETER_LATITUDE_RANGE = 0;
    }

    /* Find unit in between "[" and "]" */
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_LATITUDE_UNIT.assign( unit );

    if ( ( Input_Opt.PARAMETER_FILEINPUT ) && ( iLAT != -1 ) ) {
        /* Found latitude in parameter input file */
        Input_Opt.PARAMETER_LATITUDE_RANGE = 0;
        if ( fileInput[7].size() > 0 ) {
            for ( unsigned int i = 0; i < fileInput[7].size(); i++ )
                Input_Opt.PARAMETER_LATITUDE.push_back( fileInput[7][i] );
        } else {
            try {
                value = std::stod( tokens[0] );
                Input_Opt.PARAMETER_LATITUDE.push_back( value );
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[0] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    } else {
        /* Store in values for variable */
        for ( unsigned int i = 0; i < tokens.size(); i++ ) {
            try {
                Input_Opt.PARAMETER_LATITUDE.push_back(std::stod( tokens[i] ));
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    }

    /* ==================================================== */
    /* Pressure                                             */
    /* ==================================================== */

    /* Variable */
    variable = "Pressure";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Get line past the delimiter */
    subline = line.substr(FIRSTCOL);
    /* Look for colon */
    found = subline.find( COLON );
    if ( found != std::string::npos ) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

        if ( !Input_Opt.SIMULATION_MONTECARLO ) {
            if (tokens.size() != 3) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: " << std::endl;
                std::cout << "   --> begin:step:end" << std::endl;
                std::cout << "       or" << std::endl;
                std::cout << "   --> val1 val2 val3 ..." << std::endl;
                exit(1);
            }
            if ( ( std::stod(tokens[2]) <  std::stod(tokens[0]) ) || \
                 ( std::stod(tokens[1]) <= 0.0E+00 )              || \
                 ( std::stod(tokens[1]) >  ( std::stod(tokens[2]) - std::stod(tokens[0]) ) ) ) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: begin:step:end with begin < end, 0 < step < end - begin" << std::endl;
                exit(1);
            }
        }
        Input_Opt.PARAMETER_PRESSURE_RANGE = 1;
    }
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for " << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_PRESSURE_RANGE = 0;
    }

    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) && ( FIRST_PARAM_SWEEP ) ) {
        std::cout << " Multiple cases need to be run through APCEMM while the 'Parameter sweep?' argument is turned off!" << std::endl;
        std::cout << " These cases will be run in serial!" << std::endl;
        std::cout << " To enable the processing of multiple cases simultaneously, turn on the 'Parameter sweep' input!" << std::endl;
        FIRST_PARAM_SWEEP = 0;
    }

    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( tokens.size() > 2 ) {
            std::cout << " Wrong input for " << variable << " when MC is turned on!" << std::endl;
            std::cout << " Expected format is min max or min:max representing the range of possible values" << std::endl;
            exit(1);
        } else if ( tokens.size() == 2 ) {
            Input_Opt.PARAMETER_PRESSURE_RANGE = 1;
            sort(tokens.begin(), tokens.end());
        } else if ( tokens.size() == 1 )
            Input_Opt.PARAMETER_PRESSURE_RANGE = 0;
    }

    /* Find unit in between "[" and "]" */
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_PRESSURE_UNIT.assign( unit );

    if ( ( Input_Opt.PARAMETER_FILEINPUT ) && ( iPressure != -1 ) ) {
        /* Found pressure in parameter input file */
        Input_Opt.PARAMETER_PRESSURE_RANGE = 0;
        if ( fileInput[8].size() > 0 ) {
            for ( unsigned int i = 0; i < fileInput[8].size(); i++ ) {
                if ( fileInput[8][i] >= 0.0E+00 )
                    Input_Opt.PARAMETER_PRESSURE.push_back( fileInput[8][i] );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            }
        } else {
            try {
                value = std::stod( tokens[0] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_PRESSURE.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[0] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    } else {
        /* Store in values for variable */
        for ( unsigned int i = 0; i < tokens.size(); i++ ) {
            try {
                value = std::stod( tokens[i] );
                if ( value > 0.0E+00 )
                    Input_Opt.PARAMETER_PRESSURE.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    }


    /* Skipping line */
    getline( inputFile, line, '\n' );

    /* ==================================================== */
    /* Emission day                                         */
    /* ==================================================== */

    /* Variable */
    variable = "Emission day";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Get line past the delimiter */
    subline = line.substr(FIRSTCOL);
    /* Look for colon */
    found = subline.find( COLON );
    if ( found != std::string::npos ) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

        if ( !Input_Opt.SIMULATION_MONTECARLO ) {
            if (tokens.size() != 3) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: " << std::endl;
                std::cout << "   --> begin:step:end" << std::endl;
                std::cout << "       or" << std::endl;
                std::cout << "   --> val1 val2 val3 ..." << std::endl;
                exit(1);
            }
            if ( ( std::stod(tokens[2]) <  std::stod(tokens[0]) ) || \
                 ( std::stod(tokens[1]) <= 0.0E+00 )              || \
                 ( std::stod(tokens[1]) >  ( std::stod(tokens[2]) - std::stod(tokens[0]) ) ) ) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: begin:step:end with begin < end, 0 < step < end - begin" << std::endl;
                exit(1);
            }
        }
        Input_Opt.PARAMETER_EDAY_RANGE = 1;
    }
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for " << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_EDAY_RANGE = 0;
    }

    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) && ( FIRST_PARAM_SWEEP ) ) {
        std::cout << " Multiple cases need to be run through APCEMM while the 'Parameter sweep?' argument is turned off!" << std::endl;
        std::cout << " These cases will be run in serial!" << std::endl;
        std::cout << " To enable the processing of multiple cases simultaneously, turn on the 'Parameter sweep' input!" << std::endl;
        FIRST_PARAM_SWEEP = 0;
    }

    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( tokens.size() > 2 ) {
            std::cout << " Wrong input for " << variable << " when MC is turned on!" << std::endl;
            std::cout << " Expected format is min max or min:max representing the range of possible values" << std::endl;
            exit(1);
        } else if ( tokens.size() == 2 ) {
            Input_Opt.PARAMETER_EDAY_RANGE = 1;
            sort(tokens.begin(), tokens.end());
        } else if ( tokens.size() == 1 )
            Input_Opt.PARAMETER_EDAY_RANGE = 0;
    }

    /* Find unit in between "[" and "]" */
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_EDAY_UNIT.assign( unit );

    if ( ( Input_Opt.PARAMETER_FILEINPUT ) && ( iEDay != -1 ) ) {
        /* Found emission day in parameter input file */
        Input_Opt.PARAMETER_EDAY_RANGE = 0;
        if ( fileInput[9].size() > 0 ) {
            for ( unsigned int i = 0; i < fileInput[9].size(); i++ ) {
                if ( fileInput[9][i] >= 0.0E+00 )
                    Input_Opt.PARAMETER_EDAY.push_back( fileInput[9][i] );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            }
        } else {
            try {
                value = std::stod( tokens[0] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_EDAY.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[0] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    } else {
        /* Store in values for variable */
        for ( unsigned int i = 0; i < tokens.size(); i++ ) {
            try {
                value = ( std::stoi( tokens[i] ) - 1 % 365 ) + 1;
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_EDAY.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[i] << "' to int for " << variable << std::endl;
                exit(1);
            }
        }
    }

    /* ==================================================== */
    /* Emission time                                         */
    /* ==================================================== */

    /* Variable */
    variable = "Emission time";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Get line past the delimiter */
    subline = line.substr(FIRSTCOL);
    /* Look for colon */
    found = subline.find( COLON );
    if ( found != std::string::npos ) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

        if ( !Input_Opt.SIMULATION_MONTECARLO ) {
            if (tokens.size() != 3) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: " << std::endl;
                std::cout << "   --> begin:step:end" << std::endl;
                std::cout << "       or" << std::endl;
                std::cout << "   --> val1 val2 val3 ..." << std::endl;
                exit(1);
            }
            if ( ( std::stod(tokens[2]) <  std::stod(tokens[0]) ) || \
                 ( std::stod(tokens[1]) <= 0.0E+00 )              || \
                 ( std::stod(tokens[1]) >  ( std::stod(tokens[2]) - std::stod(tokens[0]) ) ) ) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: begin:step:end with begin < end, 0 < step < end - begin" << std::endl;
                exit(1);
            }
        }
        Input_Opt.PARAMETER_ETIME_RANGE = 1;
    }
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for " << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_ETIME_RANGE = 0;
    }

    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) && ( FIRST_PARAM_SWEEP ) ) {
        std::cout << " Multiple cases need to be run through APCEMM while the 'Parameter sweep?' argument is turned off!" << std::endl;
        std::cout << " These cases will be run in serial!" << std::endl;
        std::cout << " To enable the processing of multiple cases simultaneously, turn on the 'Parameter sweep' input!" << std::endl;
        FIRST_PARAM_SWEEP = 0;
    }

    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( tokens.size() > 2 ) {
            std::cout << " Wrong input for " << variable << " when MC is turned on!" << std::endl;
            std::cout << " Expected format is min max or min:max representing the range of possible values" << std::endl;
            exit(1);
        } else if ( tokens.size() == 2 ) {
            Input_Opt.PARAMETER_ETIME_RANGE = 1;
            sort(tokens.begin(), tokens.end());
        } else if ( tokens.size() == 1 )
            Input_Opt.PARAMETER_ETIME_RANGE = 0;
    }

    /* Find unit in between "[" and "]" */
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_ETIME_UNIT.assign( unit );

    if ( ( Input_Opt.PARAMETER_FILEINPUT ) && ( iETime != -1 ) ) {
        /* Found emission time in parameter input file */
        Input_Opt.PARAMETER_ETIME_RANGE = 0;
        if ( fileInput[10].size() > 0 ) {
            for ( unsigned int i = 0; i < fileInput[10].size(); i++ ) {
                if ( fileInput[10][i] >= 0.0E+00 )
                    Input_Opt.PARAMETER_ETIME.push_back( fileInput[10][i] );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            }
        } else {
            try {
                value = std::stod( tokens[0] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_ETIME.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[0] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    } else {
        /* Store in values for variable */
        for ( unsigned int i = 0; i < tokens.size(); i++ ) {
            try {
                value = std::stod(tokens[i]);
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_ETIME.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    }

    /* Skipping line */
    getline( inputFile, line, '\n' );

    /* ==================================================== */
    /* NOx mixing ratio                                     */
    /* ==================================================== */

    /* Variable */
    variable = "NOx mixing ratio";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Get line past the delimiter */
    subline = line.substr(FIRSTCOL);
    /* Look for colon */
    found = subline.find( COLON );
    if ( found != std::string::npos ) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

        if ( !Input_Opt.SIMULATION_MONTECARLO ) {
            if (tokens.size() != 3) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: " << std::endl;
                std::cout << "   --> begin:step:end" << std::endl;
                std::cout << "       or" << std::endl;
                std::cout << "   --> val1 val2 val3 ..." << std::endl;
                exit(1);
            }
            if ( ( std::stod(tokens[2]) <  std::stod(tokens[0]) ) || \
                 ( std::stod(tokens[1]) <= 0.0E+00 )              || \
                 ( std::stod(tokens[1]) >  ( std::stod(tokens[2]) - std::stod(tokens[0]) ) ) ) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: begin:step:end with begin < end, 0 < step < end - begin" << std::endl;
                exit(1);
            }
        }
        Input_Opt.PARAMETER_BACKG_NOX_RANGE = 1;
    }
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for " << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_BACKG_NOX_RANGE = 0;
    }

    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) && ( FIRST_PARAM_SWEEP ) ) {
        std::cout << " Multiple cases need to be run through APCEMM while the 'Parameter sweep?' argument is turned off!" << std::endl;
        std::cout << " These cases will be run in serial!" << std::endl;
        std::cout << " To enable the processing of multiple cases simultaneously, turn on the 'Parameter sweep' input!" << std::endl;
        FIRST_PARAM_SWEEP = 0;
    }

    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( tokens.size() > 2 ) {
            std::cout << " Wrong input for " << variable << " when MC is turned on!" << std::endl;
            std::cout << " Expected format is min max or min:max representing the range of possible values" << std::endl;
            exit(1);
        } else if ( tokens.size() == 2 ) {
            Input_Opt.PARAMETER_BACKG_NOX_RANGE = 1;
            sort(tokens.begin(), tokens.end());
        } else if ( tokens.size() == 1 )
            Input_Opt.PARAMETER_BACKG_NOX_RANGE = 0;
    }

    /* Find unit in between "[" and "]" */
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_BACKG_NOX_UNIT.assign( unit );

    if ( ( Input_Opt.PARAMETER_FILEINPUT ) && ( ibNOx != -1 ) ) {
        /* Found background NOx in parameter input file */
        Input_Opt.PARAMETER_BACKG_NOX_RANGE = 0;
        if ( fileInput[20].size() > 0 ) {
            for ( unsigned int i = 0; i < fileInput[20].size(); i++ ) {
                if ( fileInput[20][i] >= 0.0E+00 )
                    Input_Opt.PARAMETER_BACKG_NOX.push_back( fileInput[20][i] );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            }
        } else {
            try {
                value = std::stod( tokens[0] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_BACKG_NOX.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[0] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    } else {
        /* Store in values for variable */
        for ( unsigned int i = 0; i < tokens.size(); i++ ) {
            try {
                value = std::stod( tokens[i] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_BACKG_NOX.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    }

    /* ==================================================== */
    /* HNO3 mixing ratio                                      */
    /* ==================================================== */

    /* Variable */
    variable = "HNO3 mixing ratio";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Get line past the delimiter */
    subline = line.substr(FIRSTCOL);
    /* Look for colon */
    found = subline.find( COLON );
    if ( found != std::string::npos ) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

        if ( !Input_Opt.SIMULATION_MONTECARLO ) {
            if (tokens.size() != 3) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: " << std::endl;
                std::cout << "   --> begin:step:end" << std::endl;
                std::cout << "       or" << std::endl;
                std::cout << "   --> val1 val2 val3 ..." << std::endl;
                exit(1);
            }
            if ( ( std::stod(tokens[2]) <  std::stod(tokens[0]) ) || \
                 ( std::stod(tokens[1]) <= 0.0E+00 )              || \
                 ( std::stod(tokens[1]) >  ( std::stod(tokens[2]) - std::stod(tokens[0]) ) ) ) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: begin:step:end with begin < end, 0 < step < end - begin" << std::endl;
                exit(1);
            }
        }
        Input_Opt.PARAMETER_BACKG_HNO3_RANGE = 1;
    }
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for " << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_BACKG_HNO3_RANGE = 0;
    }

    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) && ( FIRST_PARAM_SWEEP ) ) {
        std::cout << " Multiple cases need to be run through APCEMM while the 'Parameter sweep?' argument is turned off!" << std::endl;
        std::cout << " These cases will be run in serial!" << std::endl;
        std::cout << " To enable the processing of multiple cases simultaneously, turn on the 'Parameter sweep' input!" << std::endl;
        FIRST_PARAM_SWEEP = 0;
    }

    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( tokens.size() > 2 ) {
            std::cout << " Wrong input for " << variable << " when MC is turned on!" << std::endl;
            std::cout << " Expected format is min max or min:max representing the range of possible values" << std::endl;
            exit(1);
        } else if ( tokens.size() == 2 ) {
            Input_Opt.PARAMETER_BACKG_HNO3_RANGE = 1;
            sort(tokens.begin(), tokens.end());
        } else if ( tokens.size() == 1 )
            Input_Opt.PARAMETER_BACKG_HNO3_RANGE = 0;
    }

    /* Find unit in between "[" and "]" */
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_BACKG_HNO3_UNIT.assign( unit );

    if ( ( Input_Opt.PARAMETER_FILEINPUT ) && ( ibHNO3 != -1 ) ) {
        /* Found background HNO3 in parameter input file */
        Input_Opt.PARAMETER_BACKG_HNO3_RANGE = 0;
        if ( fileInput[21].size() > 0 ) {
            for ( unsigned int i = 0; i < fileInput[21].size(); i++ ) {
                if ( fileInput[21][i] >= 0.0E+00 )
                    Input_Opt.PARAMETER_BACKG_HNO3.push_back( fileInput[21][i] );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            }
        } else {
            try {
                value = std::stod( tokens[0] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_BACKG_HNO3.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[0] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    } else {
        /* Store in values for variable */
        for ( unsigned int i = 0; i < tokens.size(); i++ ) {
            try {
                value = std::stod( tokens[i] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_BACKG_HNO3.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    }

    /* ==================================================== */
    /* O3 mixing ratio                                      */
    /* ==================================================== */

    /* Variable */
    variable = "O3 mixing ratio";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Get line past the delimiter */
    subline = line.substr(FIRSTCOL);
    /* Look for colon */
    found = subline.find( COLON );
    if ( found != std::string::npos ) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

        if ( !Input_Opt.SIMULATION_MONTECARLO ) {
            if (tokens.size() != 3) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: " << std::endl;
                std::cout << "   --> begin:step:end" << std::endl;
                std::cout << "       or" << std::endl;
                std::cout << "   --> val1 val2 val3 ..." << std::endl;
                exit(1);
            }
            if ( ( std::stod(tokens[2]) <  std::stod(tokens[0]) ) || \
                 ( std::stod(tokens[1]) <= 0.0E+00 )              || \
                 ( std::stod(tokens[1]) >  ( std::stod(tokens[2]) - std::stod(tokens[0]) ) ) ) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: begin:step:end with begin < end, 0 < step < end - begin" << std::endl;
                exit(1);
            }
        }
        Input_Opt.PARAMETER_BACKG_O3_RANGE = 1;
    }
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for " << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_BACKG_O3_RANGE = 0;
    }

    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) && ( FIRST_PARAM_SWEEP ) ) {
        std::cout << " Multiple cases need to be run through APCEMM while the 'Parameter sweep?' argument is turned off!" << std::endl;
        std::cout << " These cases will be run in serial!" << std::endl;
        std::cout << " To enable the processing of multiple cases simultaneously, turn on the 'Parameter sweep' input!" << std::endl;
        FIRST_PARAM_SWEEP = 0;
    }

    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( tokens.size() > 2 ) {
            std::cout << " Wrong input for " << variable << " when MC is turned on!" << std::endl;
            std::cout << " Expected format is min max or min:max representing the range of possible values" << std::endl;
            exit(1);
        } else if ( tokens.size() == 2 ) {
            Input_Opt.PARAMETER_BACKG_O3_RANGE = 1;
            sort(tokens.begin(), tokens.end());
        } else if ( tokens.size() == 1 )
            Input_Opt.PARAMETER_BACKG_O3_RANGE = 0;
    }


    /* Find unit in between "[" and "]" */
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_BACKG_O3_UNIT.assign( unit );

    if ( ( Input_Opt.PARAMETER_FILEINPUT ) && ( ibO3 != -1 ) ) {
        /* Found background O3 in parameter input file */
        Input_Opt.PARAMETER_BACKG_O3_RANGE = 0;
        if ( fileInput[22].size() > 0 ) {
            for ( unsigned int i = 0; i < fileInput[22].size(); i++ ) {
                if ( fileInput[22][i] >= 0.0E+00 )
                    Input_Opt.PARAMETER_BACKG_O3.push_back( fileInput[22][i] );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            }
        } else {
            try {
                value = std::stod( tokens[0] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_BACKG_O3.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[0] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    } else {
        /* Store in values for variable */
        for ( unsigned int i = 0; i < tokens.size(); i++ ) {
            try {
                value = std::stod( tokens[i] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_BACKG_O3.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    }

    /* ==================================================== */
    /* CO mixing ratio                                      */
    /* ==================================================== */

    /* Variable */
    variable = "CO mixing ratio";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Get line past the delimiter */
    subline = line.substr(FIRSTCOL);
    /* Look for colon */
    found = subline.find( COLON );
    if ( found != std::string::npos ) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

        if ( !Input_Opt.SIMULATION_MONTECARLO ) {
            if (tokens.size() != 3) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: " << std::endl;
                std::cout << "   --> begin:step:end" << std::endl;
                std::cout << "       or" << std::endl;
                std::cout << "   --> val1 val2 val3 ..." << std::endl;
                exit(1);
            }
            if ( ( std::stod(tokens[2]) <  std::stod(tokens[0]) ) || \
                 ( std::stod(tokens[1]) <= 0.0E+00 )              || \
                 ( std::stod(tokens[1]) >  ( std::stod(tokens[2]) - std::stod(tokens[0]) ) ) ) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: begin:step:end with begin < end, 0 < step < end - begin" << std::endl;
                exit(1);
            }
        }
        Input_Opt.PARAMETER_BACKG_CO_RANGE = 1;
    }
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for " << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_BACKG_CO_RANGE = 0;
    }

    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) && ( FIRST_PARAM_SWEEP ) ) {
        std::cout << " Multiple cases need to be run through APCEMM while the 'Parameter sweep?' argument is turned off!" << std::endl;
        std::cout << " These cases will be run in serial!" << std::endl;
        std::cout << " To enable the processing of multiple cases simultaneously, turn on the 'Parameter sweep' input!" << std::endl;
        FIRST_PARAM_SWEEP = 0;
    }

    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( tokens.size() > 2 ) {
            std::cout << " Wrong input for " << variable << " when MC is turned on!" << std::endl;
            std::cout << " Expected format is min max or min:max representing the range of possible values" << std::endl;
            exit(1);
        } else if ( tokens.size() == 2 ) {
            Input_Opt.PARAMETER_BACKG_CO_RANGE = 1;
            sort(tokens.begin(), tokens.end());
        } else if ( tokens.size() == 1 )
            Input_Opt.PARAMETER_BACKG_CO_RANGE = 0;
    }

    /* Find unit in between "[" and "]" */
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_BACKG_CO_UNIT.assign( unit );

    if ( ( Input_Opt.PARAMETER_FILEINPUT ) && ( ibCO != -1 ) ) {
        /* Found background CO in parameter input file */
        Input_Opt.PARAMETER_BACKG_CO_RANGE = 0;
        if ( fileInput[23].size() > 0 ) {
            for ( unsigned int i = 0; i < fileInput[23].size(); i++ ) {
                if ( fileInput[23][i] >= 0.0E+00 )
                    Input_Opt.PARAMETER_BACKG_CO.push_back( fileInput[23][i] );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            }
        } else {
            try {
                value = std::stod( tokens[0] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_BACKG_CO.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[0] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    } else {
        /* Store in values for variable */
        for ( unsigned int i = 0; i < tokens.size(); i++ ) {
            try {
                value = std::stod( tokens[i] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_BACKG_CO.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    }

    /* ==================================================== */
    /* CH4 mixing ratio                                      */
    /* ==================================================== */

    /* Variable */
    variable = "CH4 mixing ratio";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Get line past the delimiter */
    subline = line.substr(FIRSTCOL);
    /* Look for colon */
    found = subline.find( COLON );
    if ( found != std::string::npos ) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

        if ( !Input_Opt.SIMULATION_MONTECARLO ) {
            if (tokens.size() != 3) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: " << std::endl;
                std::cout << "   --> begin:step:end" << std::endl;
                std::cout << "       or" << std::endl;
                std::cout << "   --> val1 val2 val3 ..." << std::endl;
                exit(1);
            }
            if ( ( std::stod(tokens[2]) <  std::stod(tokens[0]) ) || \
                 ( std::stod(tokens[1]) <= 0.0E+00 )              || \
                 ( std::stod(tokens[1]) >  ( std::stod(tokens[2]) - std::stod(tokens[0]) ) ) ) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: begin:step:end with begin < end, 0 < step < end - begin" << std::endl;
                exit(1);
            }
        }
        Input_Opt.PARAMETER_BACKG_CH4_RANGE = 1;
    }
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for " << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_BACKG_CH4_RANGE = 0;
    }

    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) && ( FIRST_PARAM_SWEEP ) ) {
        std::cout << " Multiple cases need to be run through APCEMM while the 'Parameter sweep?' argument is turned off!" << std::endl;
        std::cout << " These cases will be run in serial!" << std::endl;
        std::cout << " To enable the processing of multiple cases simultaneously, turn on the 'Parameter sweep' input!" << std::endl;
        FIRST_PARAM_SWEEP = 0;
    }

    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( tokens.size() > 2 ) {
            std::cout << " Wrong input for " << variable << " when MC is turned on!" << std::endl;
            std::cout << " Expected format is min max or min:max representing the range of possible values" << std::endl;
            exit(1);
        } else if ( tokens.size() == 2 ) {
            Input_Opt.PARAMETER_BACKG_CH4_RANGE = 1;
            sort(tokens.begin(), tokens.end());
        } else if ( tokens.size() == 1 )
            Input_Opt.PARAMETER_BACKG_CH4_RANGE = 0;
    }

    /* Find unit in between "[" and "]" */
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_BACKG_CH4_UNIT.assign( unit );

    if ( ( Input_Opt.PARAMETER_FILEINPUT ) && ( ibCH4 != -1 ) ) {
        /* Found background CH4 in parameter input file */
        Input_Opt.PARAMETER_BACKG_CH4_RANGE = 0;
        if ( fileInput[24].size() > 0 ) {
            for ( unsigned int i = 0; i < fileInput[24].size(); i++ ) {
                if ( fileInput[24][i] >= 0.0E+00 )
                    Input_Opt.PARAMETER_BACKG_CH4.push_back( fileInput[24][i] );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            }
        } else {
            try {
                value = std::stod( tokens[0] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_BACKG_CH4.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[0] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    } else {
        /* Store in values for variable */
        for ( unsigned int i = 0; i < tokens.size(); i++ ) {
            try {
                value = std::stod( tokens[i] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_BACKG_CH4.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    }

    /* ==================================================== */
    /* SO2 mixing ratio                                      */
    /* ==================================================== */

    /* Variable */
    variable = "SO2 mixing ratio";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Get line past the delimiter */
    subline = line.substr(FIRSTCOL);
    /* Look for colon */
    found = subline.find( COLON );
    if ( found != std::string::npos ) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

        if ( !Input_Opt.SIMULATION_MONTECARLO ) {
            if (tokens.size() != 3) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: " << std::endl;
                std::cout << "   --> begin:step:end" << std::endl;
                std::cout << "       or" << std::endl;
                std::cout << "   --> val1 val2 val3 ..." << std::endl;
                exit(1);
            }
            if ( ( std::stod(tokens[2]) <  std::stod(tokens[0]) ) || \
                 ( std::stod(tokens[1]) <= 0.0E+00 )              || \
                 ( std::stod(tokens[1]) >  ( std::stod(tokens[2]) - std::stod(tokens[0]) ) ) ) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: begin:step:end with begin < end, 0 < step < end - begin" << std::endl;
                exit(1);
            }
        }
        Input_Opt.PARAMETER_BACKG_SO2_RANGE = 1;
    }
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for " << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_BACKG_SO2_RANGE = 0;
    }

    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) && ( FIRST_PARAM_SWEEP ) ) {
        std::cout << " Multiple cases need to be run through APCEMM while the 'Parameter sweep?' argument is turned off!" << std::endl;
        std::cout << " These cases will be run in serial!" << std::endl;
        std::cout << " To enable the processing of multiple cases simultaneously, turn on the 'Parameter sweep' input!" << std::endl;
        FIRST_PARAM_SWEEP = 0;
    }

    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( tokens.size() > 2 ) {
            std::cout << " Wrong input for " << variable << " when MC is turned on!" << std::endl;
            std::cout << " Expected format is min max or min:max representing the range of possible values" << std::endl;
            exit(1);
        } else if ( tokens.size() == 2 ) {
            Input_Opt.PARAMETER_BACKG_SO2_RANGE = 1;
            sort(tokens.begin(), tokens.end());
        } else if ( tokens.size() == 1 )
            Input_Opt.PARAMETER_BACKG_SO2_RANGE = 0;
    }

    /* Find unit in between "[" and "]" */
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_BACKG_SO2_UNIT.assign( unit );

    if ( ( Input_Opt.PARAMETER_FILEINPUT ) && ( ibSO2 != -1 ) ) {
        /* Found background SO2 in parameter input file */
        Input_Opt.PARAMETER_BACKG_SO2_RANGE = 0;
        if ( fileInput[25].size() > 0 ) {
            for ( unsigned int i = 0; i < fileInput[25].size(); i++ ) {
                if ( fileInput[25][i] >= 0.0E+00 )
                    Input_Opt.PARAMETER_BACKG_SO2.push_back( fileInput[25][i] );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            }
        } else {
            try {
                value = std::stod( tokens[0] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_BACKG_SO2.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[0] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    } else {
        /* Store in values for variable */
        for ( unsigned int i = 0; i < tokens.size(); i++ ) {
            try {
                value = std::stod( tokens[i] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_BACKG_SO2.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    }


    /* Skipping line */
    getline( inputFile, line, '\n' );

    /* ==================================================== */
    /* NOx emission index                                   */
    /* ==================================================== */

    /* Variable */
    variable = "NOx emission index";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Get line past the delimiter */
    subline = line.substr(FIRSTCOL);
    /* Look for colon */
    found = subline.find( COLON );
    if ( found != std::string::npos ) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

        if ( !Input_Opt.SIMULATION_MONTECARLO ) {
            if (tokens.size() != 3) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: " << std::endl;
                std::cout << "   --> begin:step:end" << std::endl;
                std::cout << "       or" << std::endl;
                std::cout << "   --> val1 val2 val3 ..." << std::endl;
                exit(1);
            }
            if ( ( std::stod(tokens[2]) <  std::stod(tokens[0]) ) || \
                 ( std::stod(tokens[1]) <= 0.0E+00 )              || \
                 ( std::stod(tokens[1]) >  ( std::stod(tokens[2]) - std::stod(tokens[0]) ) ) ) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: begin:step:end with begin < end, 0 < step < end - begin" << std::endl;
                exit(1);
            }
        }
        Input_Opt.PARAMETER_EI_NOX_RANGE = 1;
    }
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for " << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_EI_NOX_RANGE = 0;
    }

    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) && ( FIRST_PARAM_SWEEP ) ) {
        std::cout << " Multiple cases need to be run through APCEMM while the 'Parameter sweep?' argument is turned off!" << std::endl;
        std::cout << " These cases will be run in serial!" << std::endl;
        std::cout << " To enable the processing of multiple cases simultaneously, turn on the 'Parameter sweep' input!" << std::endl;
        FIRST_PARAM_SWEEP = 0;
    }

    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( tokens.size() > 2 ) {
            std::cout << " Wrong input for " << variable << " when MC is turned on!" << std::endl;
            std::cout << " Expected format is min max or min:max representing the range of possible values" << std::endl;
            exit(1);
        } else if ( tokens.size() == 2 ) {
            Input_Opt.PARAMETER_EI_NOX_RANGE = 1;
            sort(tokens.begin(), tokens.end());
        } else if ( tokens.size() == 1 )
            Input_Opt.PARAMETER_EI_NOX_RANGE = 0;
    }

    /* Find unit in between "[" and "]" */
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_EI_NOX_UNIT.assign( unit );

    if ( ( Input_Opt.PARAMETER_FILEINPUT ) && ( ( iEINOx != -1 ) || ( iNOxFlow != -1 ) ) ) {
        /* Found EI NOx in parameter input file
         * If we found NOx flow, we are storing the NOx flow in
         * PARAMETER_EI_NOX for now but convert to EI once we have data for
         * fuel flow */
        Input_Opt.PARAMETER_EI_NOX_RANGE = 0;
        if ( fileInput[11].size() > 0 ) {
            for ( unsigned int i = 0; i < fileInput[11].size(); i++ ) {
                if ( fileInput[11][i] >= 0.0E+00 )
                    Input_Opt.PARAMETER_EI_NOX.push_back( fileInput[11][i] );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            }
        } else {
            try {
                value = std::stod( tokens[0] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_EI_NOX.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[0] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    } else {
        /* Store in values for variable */
        for ( unsigned int i = 0; i < tokens.size(); i++ ) {
            try {
                value = std::stod( tokens[i] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_EI_NOX.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    }

    /* ==================================================== */
    /* CO emission index                                    */
    /* ==================================================== */

    /* Variable */
    variable = "CO emission index";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Get line past the delimiter */
    subline = line.substr(FIRSTCOL);
    /* Look for colon */
    found = subline.find( COLON );
    if ( found != std::string::npos ) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

        if ( !Input_Opt.SIMULATION_MONTECARLO ) {
            if (tokens.size() != 3) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: " << std::endl;
                std::cout << "   --> begin:step:end" << std::endl;
                std::cout << "       or" << std::endl;
                std::cout << "   --> val1 val2 val3 ..." << std::endl;
                exit(1);
            }
            if ( ( std::stod(tokens[2]) <  std::stod(tokens[0]) ) || \
                 ( std::stod(tokens[1]) <= 0.0E+00 )              || \
                 ( std::stod(tokens[1]) >  ( std::stod(tokens[2]) - std::stod(tokens[0]) ) ) ) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: begin:step:end with begin < end, 0 < step < end - begin" << std::endl;
                exit(1);
            }
        }
        Input_Opt.PARAMETER_EI_CO_RANGE = 1;
    }
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for " << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_EI_CO_RANGE = 0;
    }

    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) && ( FIRST_PARAM_SWEEP ) ) {
        std::cout << " Multiple cases need to be run through APCEMM while the 'Parameter sweep?' argument is turned off!" << std::endl;
        std::cout << " These cases will be run in serial!" << std::endl;
        std::cout << " To enable the processing of multiple cases simultaneously, turn on the 'Parameter sweep' input!" << std::endl;
        FIRST_PARAM_SWEEP = 0;
    }

    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( tokens.size() > 2 ) {
            std::cout << " Wrong input for " << variable << " when MC is turned on!" << std::endl;
            std::cout << " Expected format is min max or min:max representing the range of possible values" << std::endl;
            exit(1);
        } else if ( tokens.size() == 2 ) {
            Input_Opt.PARAMETER_EI_CO_RANGE = 1;
            sort(tokens.begin(), tokens.end());
        } else if ( tokens.size() == 1 )
            Input_Opt.PARAMETER_EI_CO_RANGE = 0;
    }

    /* Find unit in between "[" and "]" */
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_EI_CO_UNIT.assign( unit );

    if ( ( Input_Opt.PARAMETER_FILEINPUT ) && ( iEICO != -1 ) ) {
        /* Found EI CO in parameter input file */
        Input_Opt.PARAMETER_EI_CO_RANGE = 0;
        if ( fileInput[12].size() > 0 ) {
            for ( unsigned int i = 0; i < fileInput[12].size(); i++ ) {
                if ( fileInput[12][i] >= 0.0E+00 )
                    Input_Opt.PARAMETER_EI_CO.push_back( fileInput[12][i] );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            }
        } else {
            try {
                value = std::stod( tokens[0] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_EI_CO.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[0] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    } else {
        /* Store in values for variable */
        for ( unsigned int i = 0; i < tokens.size(); i++ ) {
            try {
                value = std::stod( tokens[i] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_EI_CO.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    }

    /* ==================================================== */
    /* UHC emission index                                   */
    /* ==================================================== */

    /* Variable */
    variable = "UHC emission index";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Get line past the delimiter */
    subline = line.substr(FIRSTCOL);
    /* Look for colon */
    found = subline.find( COLON );
    if ( found != std::string::npos ) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

        if ( !Input_Opt.SIMULATION_MONTECARLO ) {
            if (tokens.size() != 3) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: " << std::endl;
                std::cout << "   --> begin:step:end" << std::endl;
                std::cout << "       or" << std::endl;
                std::cout << "   --> val1 val2 val3 ..." << std::endl;
                exit(1);
            }
            if ( ( std::stod(tokens[2]) <  std::stod(tokens[0]) ) || \
                 ( std::stod(tokens[1]) <= 0.0E+00 )              || \
                 ( std::stod(tokens[1]) >  ( std::stod(tokens[2]) - std::stod(tokens[0]) ) ) ) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: begin:step:end with begin < end, 0 < step < end - begin" << std::endl;
                exit(1);
            }
        }
        Input_Opt.PARAMETER_EI_UHC_RANGE = 1;
    }
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for " << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_EI_UHC_RANGE = 0;
    }

    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) && ( FIRST_PARAM_SWEEP ) ) {
        std::cout << " Multiple cases need to be run through APCEMM while the 'Parameter sweep?' argument is turned off!" << std::endl;
        std::cout << " These cases will be run in serial!" << std::endl;
        std::cout << " To enable the processing of multiple cases simultaneously, turn on the 'Parameter sweep' input!" << std::endl;
        FIRST_PARAM_SWEEP = 0;
    }

    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( tokens.size() > 2 ) {
            std::cout << " Wrong input for " << variable << " when MC is turned on!" << std::endl;
            std::cout << " Expected format is min max or min:max representing the range of possible values" << std::endl;
            exit(1);
        } else if ( tokens.size() == 2 ) {
            Input_Opt.PARAMETER_EI_UHC_RANGE = 1;
            sort(tokens.begin(), tokens.end());
        } else if ( tokens.size() == 1 )
            Input_Opt.PARAMETER_EI_UHC_RANGE = 0;
    }

    /* Find unit in between "[" and "]" */
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_EI_UHC_UNIT.assign( unit );

    if ( ( Input_Opt.PARAMETER_FILEINPUT ) && ( iEIUHC != -1 ) ) {
        /* Found EI UHC in parameter input file */
        Input_Opt.PARAMETER_EI_UHC_RANGE = 0;
        if ( fileInput[13].size() > 0 ) {
            for ( unsigned int i = 0; i < fileInput[13].size(); i++ ) {
                if ( fileInput[13][i] >= 0.0E+00 )
                    Input_Opt.PARAMETER_EI_UHC.push_back( fileInput[13][i] );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            }
        } else {
            try {
                value = std::stod( tokens[0] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_EI_UHC.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[0] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    } else {
        /* Store in values for variable */
        for ( unsigned int i = 0; i < tokens.size(); i++ ) {
            try {
                value = std::stod( tokens[i] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_EI_UHC.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    }

    /* ==================================================== */
    /* SO2 emission index                                   */
    /* ==================================================== */

    /* Variable */
    variable = "SO2 emission index";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Get line past the delimiter */
    subline = line.substr(FIRSTCOL);
    /* Look for colon */
    found = subline.find( COLON );
    if ( found != std::string::npos ) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

        if ( !Input_Opt.SIMULATION_MONTECARLO ) {
            if (tokens.size() != 3) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: " << std::endl;
                std::cout << "   --> begin:step:end" << std::endl;
                std::cout << "       or" << std::endl;
                std::cout << "   --> val1 val2 val3 ..." << std::endl;
                exit(1);
            }
            if ( ( std::stod(tokens[2]) <  std::stod(tokens[0]) ) || \
                 ( std::stod(tokens[1]) <= 0.0E+00 )              || \
                 ( std::stod(tokens[1]) >  ( std::stod(tokens[2]) - std::stod(tokens[0]) ) ) ) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: begin:step:end with begin < end, 0 < step < end - begin" << std::endl;
                exit(1);
            }
        }
        Input_Opt.PARAMETER_EI_SO2_RANGE = 1;
    }
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for " << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_EI_SO2_RANGE = 0;
    }

    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) && ( FIRST_PARAM_SWEEP ) ) {
        std::cout << " Multiple cases need to be run through APCEMM while the 'Parameter sweep?' argument is turned off!" << std::endl;
        std::cout << " These cases will be run in serial!" << std::endl;
        std::cout << " To enable the processing of multiple cases simultaneously, turn on the 'Parameter sweep' input!" << std::endl;
        FIRST_PARAM_SWEEP = 0;
    }

    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( tokens.size() > 2 ) {
            std::cout << " Wrong input for " << variable << " when MC is turned on!" << std::endl;
            std::cout << " Expected format is min max or min:max representing the range of possible values" << std::endl;
            exit(1);
        } else if ( tokens.size() == 2 ) {
            Input_Opt.PARAMETER_EI_SO2_RANGE = 1;
            sort(tokens.begin(), tokens.end());
        } else if ( tokens.size() == 1 )
            Input_Opt.PARAMETER_EI_SO2_RANGE = 0;
    }

    /* Find unit in between "[" and "]" */
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_EI_SO2_UNIT.assign( unit );

    if ( ( Input_Opt.PARAMETER_FILEINPUT ) && ( iEISO2 != -1 ) ) {
        /* Found EI SO2 in parameter input file */
        Input_Opt.PARAMETER_EI_SO2_RANGE = 0;
        if ( fileInput[14].size() > 0 ) {
            for ( unsigned int i = 0; i < fileInput[14].size(); i++ ) {
                if ( fileInput[14][i] >= 0.0E+00 )
                    Input_Opt.PARAMETER_EI_SO2.push_back( fileInput[14][i] );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            }
        } else {
            try {
                value = std::stod( tokens[0] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_EI_SO2.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[0] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    } else {
        /* Store in values for variable */
        for ( unsigned int i = 0; i < tokens.size(); i++ ) {
            try {
                value = std::stod( tokens[i] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_EI_SO2.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    }

    /* ==================================================== */
    /* SO2 to SO4 conversion                                */
    /* ==================================================== */

    /* Variable */
    variable = "SO2 to SO4 conversion";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Get line past the delimiter */
    subline = line.substr(FIRSTCOL);
    /* Look for colon */
    found = subline.find( COLON );
    if ( found != std::string::npos ) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

        if ( !Input_Opt.SIMULATION_MONTECARLO ) {
            if (tokens.size() != 3) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: " << std::endl;
                std::cout << "   --> begin:step:end" << std::endl;
                std::cout << "       or" << std::endl;
                std::cout << "   --> val1 val2 val3 ..." << std::endl;
                exit(1);
            }
            if ( ( std::stod(tokens[2]) <  std::stod(tokens[0]) ) || \
                 ( std::stod(tokens[1]) <= 0.0E+00 )              || \
                 ( std::stod(tokens[1]) >  ( std::stod(tokens[2]) - std::stod(tokens[0]) ) ) ) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: begin:step:end with begin < end, 0 < step < end - begin" << std::endl;
                exit(1);
            }
        }
        Input_Opt.PARAMETER_EI_SO2TOSO4_RANGE = 1;
    }
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for " << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_EI_SO2TOSO4_RANGE = 0;
    }

    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) && ( FIRST_PARAM_SWEEP ) ) {
        std::cout << " Multiple cases need to be run through APCEMM while the 'Parameter sweep?' argument is turned off!" << std::endl;
        std::cout << " These cases will be run in serial!" << std::endl;
        std::cout << " To enable the processing of multiple cases simultaneously, turn on the 'Parameter sweep' input!" << std::endl;
        FIRST_PARAM_SWEEP = 0;
    }

    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( tokens.size() > 2 ) {
            std::cout << " Wrong input for " << variable << " when MC is turned on!" << std::endl;
            std::cout << " Expected format is min max or min:max representing the range of possible values" << std::endl;
            exit(1);
        } else if ( tokens.size() == 2 ) {
            Input_Opt.PARAMETER_EI_SO2TOSO4_RANGE = 1;
            sort(tokens.begin(), tokens.end());
        } else if ( tokens.size() == 1 )
            Input_Opt.PARAMETER_EI_SO2TOSO4_RANGE = 0;
    }

    /* Find unit in between "[" and "]" */
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_EI_SO2TOSO4_UNIT.assign( unit );

    if ( ( Input_Opt.PARAMETER_FILEINPUT ) && ( iSO2toSO4 != -1 ) ) {
        /* Found EI SO2 to SO4 in parameter input file */
        Input_Opt.PARAMETER_EI_SO2TOSO4_RANGE = 0;
        if ( fileInput[15].size() > 0 ) {
            for ( unsigned int i = 0; i < fileInput[15].size(); i++ ) {
                if ( fileInput[15][i] >= 0.0E+00 )
                    Input_Opt.PARAMETER_EI_SO2TOSO4.push_back( fileInput[15][i] );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            }
        } else {
            try {
                value = std::stod( tokens[0] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_EI_SO2TOSO4.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[0] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    } else {
        /* Store in values for variable */
        for ( unsigned int i = 0; i < tokens.size(); i++ ) {
            try {
                value = std::stod( tokens[i] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_EI_SO2TOSO4.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    }

    /* ==================================================== */
    /* Soot emission index                                  */
    /* ==================================================== */

    /* Variable */
    variable = "Soot emission index";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Get line past the delimiter */
    subline = line.substr(FIRSTCOL);
    /* Look for colon */
    found = subline.find( COLON );
    if ( found != std::string::npos ) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

        if ( !Input_Opt.SIMULATION_MONTECARLO ) {
            if (tokens.size() != 3) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: " << std::endl;
                std::cout << "   --> begin:step:end" << std::endl;
                std::cout << "       or" << std::endl;
                std::cout << "   --> val1 val2 val3 ..." << std::endl;
                exit(1);
            }
            if ( ( std::stod(tokens[2]) <  std::stod(tokens[0]) ) || \
                 ( std::stod(tokens[1]) <= 0.0E+00 )              || \
                 ( std::stod(tokens[1]) >  ( std::stod(tokens[2]) - std::stod(tokens[0]) ) ) ) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: begin:step:end with begin < end, 0 < step < end - begin" << std::endl;
                exit(1);
            }
        }
        Input_Opt.PARAMETER_EI_SOOT_RANGE = 1;
    }
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for " << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_EI_SOOT_RANGE = 0;
    }

    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) && ( FIRST_PARAM_SWEEP ) ) {
        std::cout << " Multiple cases need to be run through APCEMM while the 'Parameter sweep?' argument is turned off!" << std::endl;
        std::cout << " These cases will be run in serial!" << std::endl;
        std::cout << " To enable the processing of multiple cases simultaneously, turn on the 'Parameter sweep' input!" << std::endl;
        FIRST_PARAM_SWEEP = 0;
    }

    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( tokens.size() > 2 ) {
            std::cout << " Wrong input for " << variable << " when MC is turned on!" << std::endl;
            std::cout << " Expected format is min max or min:max representing the range of possible values" << std::endl;
            exit(1);
        } else if ( tokens.size() == 2 ) {
            Input_Opt.PARAMETER_EI_SOOT_RANGE = 1;
            sort(tokens.begin(), tokens.end());
        } else if ( tokens.size() == 1 )
            Input_Opt.PARAMETER_EI_SOOT_RANGE = 0;
    }

    /* Find unit in between "[" and "]" */
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_EI_SOOT_UNIT.assign( unit );

    if ( ( Input_Opt.PARAMETER_FILEINPUT ) && ( iEISoot != -1 ) ) {
        /* Found EI Soot in parameter input file */
        Input_Opt.PARAMETER_EI_SOOT_RANGE = 0;
        if ( fileInput[16].size() > 0 ) {
            for ( unsigned int i = 0; i < fileInput[16].size(); i++ ) {
                if ( fileInput[16][i] >= 0.0E+00 )
                    Input_Opt.PARAMETER_EI_SOOT.push_back( fileInput[16][i] );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            }
        } else {
            try {
                value = std::stod( tokens[0] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_EI_SOOT.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[0] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    } else {
        /* Store in values for variable */
        for ( unsigned int i = 0; i < tokens.size(); i++ ) {
            try {
                value = std::stod( tokens[i] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_EI_SOOT.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    }

    /* ==================================================== */
    /* Soot radius                                          */
    /* ==================================================== */

    /* Variable */
    variable = "Soot Radius";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Get line past the delimiter */
    subline = line.substr(FIRSTCOL);
    /* Look for colon */
    found = subline.find( COLON );
    if ( found != std::string::npos ) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

        if ( !Input_Opt.SIMULATION_MONTECARLO ) {
            if (tokens.size() != 3) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: " << std::endl;
                std::cout << "   --> begin:step:end" << std::endl;
                std::cout << "       or" << std::endl;
                std::cout << "   --> val1 val2 val3 ..." << std::endl;
                exit(1);
            }
            if ( ( std::stod(tokens[2]) <  std::stod(tokens[0]) ) || \
                 ( std::stod(tokens[1]) <= 0.0E+00 )              || \
                 ( std::stod(tokens[1]) >  ( std::stod(tokens[2]) - std::stod(tokens[0]) ) ) ) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: begin:step:end with begin < end, 0 < step < end - begin" << std::endl;
                exit(1);
            }
        }
        Input_Opt.PARAMETER_EI_SOOTRAD_RANGE = 1;
    }
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for " << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_EI_SOOTRAD_RANGE = 0;
    }

    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) && ( FIRST_PARAM_SWEEP ) ) {
        std::cout << " Multiple cases need to be run through APCEMM while the 'Parameter sweep?' argument is turned off!" << std::endl;
        std::cout << " These cases will be run in serial!" << std::endl;
        std::cout << " To enable the processing of multiple cases simultaneously, turn on the 'Parameter sweep' input!" << std::endl;
        FIRST_PARAM_SWEEP = 0;
    }

    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( tokens.size() > 2 ) {
            std::cout << " Wrong input for " << variable << " when MC is turned on!" << std::endl;
            std::cout << " Expected format is min max or min:max representing the range of possible values" << std::endl;
            exit(1);
        } else if ( tokens.size() == 2 ) {
            Input_Opt.PARAMETER_EI_SOOTRAD_RANGE = 1;
            sort(tokens.begin(), tokens.end());
        } else if ( tokens.size() == 1 )
            Input_Opt.PARAMETER_EI_SOOTRAD_RANGE = 0;
    }

    /* Find unit in between "[" and "]" */
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_EI_SOOTRAD_UNIT.assign( unit );

    if ( ( Input_Opt.PARAMETER_FILEINPUT ) && ( iEISootRad != -1 ) ) {
        /* Found EI SootRad in parameter input file */
        Input_Opt.PARAMETER_EI_SOOTRAD_RANGE = 0;
        if ( fileInput[17].size() > 0 ) {
            for ( unsigned int i = 0; i < fileInput[17].size(); i++ ) {
                if ( fileInput[17][i] >= 0.0E+00 )
                    Input_Opt.PARAMETER_EI_SOOTRAD.push_back( fileInput[17][i] );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            }
        } else {
            try {
                value = std::stod( tokens[0] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_EI_SOOTRAD.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[0] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    } else {
        /* Store in values for variable */
        for ( unsigned int i = 0; i < tokens.size(); i++ ) {
            try {
                value = std::stod( tokens[i] );
                if ( value > 0.0E+00 )
                    Input_Opt.PARAMETER_EI_SOOTRAD.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    }

    /* ==================================================== */
    /* Fuel flow                                            */
    /* ==================================================== */

    /* Variable */
    variable = "Fuel flow";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Get line past the delimiter */
    subline = line.substr(FIRSTCOL);
    /* Look for colon */
    found = subline.find( COLON );
    if ( found != std::string::npos ) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

        if ( !Input_Opt.SIMULATION_MONTECARLO ) {
            if (tokens.size() != 3) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: " << std::endl;
                std::cout << "   --> begin:step:end" << std::endl;
                std::cout << "       or" << std::endl;
                std::cout << "   --> val1 val2 val3 ..." << std::endl;
                exit(1);
            }
            if ( ( std::stod(tokens[2]) <  std::stod(tokens[0]) ) || \
                 ( std::stod(tokens[1]) <= 0.0E+00 )              || \
                 ( std::stod(tokens[1]) >  ( std::stod(tokens[2]) - std::stod(tokens[0]) ) ) ) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: begin:step:end with begin < end, 0 < step < end - begin" << std::endl;
                exit(1);
            }
        }
        Input_Opt.PARAMETER_FF_RANGE = 1;
    }
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for " << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_FF_RANGE = 0;
    }

    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) && ( FIRST_PARAM_SWEEP ) ) {
        std::cout << " Multiple cases need to be run through APCEMM while the 'Parameter sweep?' argument is turned off!" << std::endl;
        std::cout << " These cases will be run in serial!" << std::endl;
        std::cout << " To enable the processing of multiple cases simultaneously, turn on the 'Parameter sweep' input!" << std::endl;
        FIRST_PARAM_SWEEP = 0;
    }

    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( tokens.size() > 2 ) {
            std::cout << " Wrong input for " << variable << " when MC is turned on!" << std::endl;
            std::cout << " Expected format is min max or min:max representing the range of possible values" << std::endl;
            exit(1);
        } else if ( tokens.size() == 2 ) {
            Input_Opt.PARAMETER_FF_RANGE = 1;
            sort(tokens.begin(), tokens.end());
        } else if ( tokens.size() == 1 )
            Input_Opt.PARAMETER_FF_RANGE = 0;
    }

    /* Find unit in between "[" and "]" */
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_FF_UNIT.assign( unit );

    if ( ( Input_Opt.PARAMETER_FILEINPUT ) && ( itotFuelFlow != -1 ) ) {
        /* Found Fuel flow in parameter input file */
        Input_Opt.PARAMETER_FF_RANGE = 0;
        if ( fileInput[18].size() > 0 ) {
            for ( unsigned int i = 0; i < fileInput[18].size(); i++ ) {
                if ( fileInput[18][i] >= 0.0E+00 )
                    Input_Opt.PARAMETER_FF.push_back( fileInput[18][i] );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            }
        } else {
            try {
                value = std::stod( tokens[0] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_FF.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[0] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    } else {
        /* Store in values for variable */
        for ( unsigned int i = 0; i < tokens.size(); i++ ) {
            try {
                value = std::stod( tokens[i] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_FF.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    }
    
    /* ==================================================== */
    /* Aircraft mass                                        */
    /* ==================================================== */

    /* Variable */
    variable = "Aircraft mass";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Get line past the delimiter */
    subline = line.substr(FIRSTCOL);
    /* Look for colon */
    found = subline.find( COLON );
    if ( found != std::string::npos ) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

        if ( !Input_Opt.SIMULATION_MONTECARLO ) {
            if (tokens.size() != 3) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: " << std::endl;
                std::cout << "   --> begin:step:end" << std::endl;
                std::cout << "       or" << std::endl;
                std::cout << "   --> val1 val2 val3 ..." << std::endl;
                exit(1);
            }
            if ( ( std::stod(tokens[2]) <  std::stod(tokens[0]) ) || \
                 ( std::stod(tokens[1]) <= 0.0E+00 )              || \
                 ( std::stod(tokens[1]) >  ( std::stod(tokens[2]) - std::stod(tokens[0]) ) ) ) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: begin:step:end with begin < end, 0 < step < end - begin" << std::endl;
                exit(1);
            }
        }
        Input_Opt.PARAMETER_AMASS_RANGE = 1;
    } 
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for " << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_AMASS_RANGE = 0;
    }
   
    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) ) {
        std::cout << " APCEMM cannot accept multiple cases when the 'parameter sweep?' argument is turned off! Aborting." << std::endl;
        exit(1);
    }
    
    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( tokens.size() > 2 ) {
            std::cout << " Wrong input for " << variable << " when MC is turned on!" << std::endl;
            std::cout << " Expected format is min max or min:max representing the range of possible values" << std::endl;
            exit(1);
        } else if ( tokens.size() == 2 ) {
            Input_Opt.PARAMETER_AMASS_RANGE = 1;
            sort(tokens.begin(), tokens.end());
        } else if ( tokens.size() == 1 )
            Input_Opt.PARAMETER_AMASS_RANGE = 0;
    }
    
    /* Find unit in between "[" and "]" */ 
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_AMASS_UNIT.assign( unit );

    if ( ( Input_Opt.PARAMETER_FILEINPUT ) && ( iaircraftMass != -1 ) ) {
        /* Found Fuel flow in parameter input file */
        Input_Opt.PARAMETER_AMASS_RANGE = 0;
        if ( fileInput[19].size() > 0 ) {
            for ( unsigned int i = 0; i < fileInput[19].size(); i++ ) {
                if ( fileInput[19][i] >= 0.0E+00 )
                    Input_Opt.PARAMETER_AMASS.push_back( fileInput[19][i] );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            }
        } else {
            try {
                value = std::stod( tokens[0] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_AMASS.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[0] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    } else {
        /* Store in values for variable */
        for ( unsigned int i = 0; i < tokens.size(); i++ ) {
            try {
                value = std::stod( tokens[i] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_AMASS.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    }

    /* ==================================================== */
    /* Flight speed                                         */
    /* ==================================================== */

    /* Variable */
    variable = "Flight speed";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Get line past the delimiter */
    subline = line.substr(FIRSTCOL);
    /* Look for colon */
    found = subline.find( COLON );
    if ( found != std::string::npos ) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

        if ( !Input_Opt.SIMULATION_MONTECARLO ) {
            if (tokens.size() != 3) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: " << std::endl;
                std::cout << "   --> begin:step:end" << std::endl;
                std::cout << "       or" << std::endl;
                std::cout << "   --> val1 val2 val3 ..." << std::endl;
                exit(1);
            }
            if ( ( std::stod(tokens[2]) <  std::stod(tokens[0]) ) || \
                 ( std::stod(tokens[1]) <= 0.0E+00 )              || \
                 ( std::stod(tokens[1]) >  ( std::stod(tokens[2]) - std::stod(tokens[0]) ) ) ) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: begin:step:end with begin < end, 0 < step < end - begin" << std::endl;
                exit(1);
            }
        }
        Input_Opt.PARAMETER_FSPEED_RANGE = 1;
    } 
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for " << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_FSPEED_RANGE = 0;
    }
   
    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) ) {
        std::cout << " APCEMM cannot accept multiple cases when the 'parameter sweep?' argument is turned off! Aborting." << std::endl;
        exit(1);
    }
    
    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( tokens.size() > 2 ) {
            std::cout << " Wrong input for " << variable << " when MC is turned on!" << std::endl;
            std::cout << " Expected format is min max or min:max representing the range of possible values" << std::endl;
            exit(1);
        } else if ( tokens.size() == 2 ) {
            Input_Opt.PARAMETER_FSPEED_RANGE = 1;
            sort(tokens.begin(), tokens.end());
        } else if ( tokens.size() == 1 )
            Input_Opt.PARAMETER_FSPEED_RANGE = 0;
    }
    
    /* Find unit in between "[" and "]" */ 
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_FSPEED_UNIT.assign( unit );

    if ( ( Input_Opt.PARAMETER_FILEINPUT ) && ( iflightSpeed != -1 ) ) {
        /* Found flight speed in parameter input file */
        Input_Opt.PARAMETER_FSPEED_RANGE = 0;
        if ( fileInput[26].size() > 0 ) {
            for ( unsigned int i = 0; i < fileInput[26].size(); i++ ) {
                if ( fileInput[26][i] >= 0.0E+00 )
                    Input_Opt.PARAMETER_FSPEED.push_back( fileInput[26][i] );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            }
        } else {
            try {
                value = std::stod( tokens[0] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_FSPEED.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[0] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    } else {
        /* Store in values for variable */
        for ( unsigned int i = 0; i < tokens.size(); i++ ) {
            try {
                value = std::stod( tokens[i] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_FSPEED.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    }

    /* ==================================================== */
    /* Number of engines                                    */
    /* ==================================================== */

    /* Variable */
    variable = "Number of engines";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Get line past the delimiter */
    subline = line.substr(FIRSTCOL);
    /* Look for colon */
    found = subline.find( COLON );
    if ( found != std::string::npos ) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

        if ( !Input_Opt.SIMULATION_MONTECARLO ) {
            if (tokens.size() != 3) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: " << std::endl;
                std::cout << "   --> begin:step:end" << std::endl;
                std::cout << "       or" << std::endl;
                std::cout << "   --> val1 val2 val3 ..." << std::endl;
                exit(1);
            }
            if ( ( std::stod(tokens[2]) <  std::stod(tokens[0]) ) || \
                 ( std::stod(tokens[1]) <= 0.0E+00 )              || \
                 ( std::stod(tokens[1]) >  ( std::stod(tokens[2]) - std::stod(tokens[0]) ) ) ) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: begin:step:end with begin < end, 0 < step < end - begin" << std::endl;
                exit(1);
            }
        }
        Input_Opt.PARAMETER_NUMENG_RANGE = 1;
    } 
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for " << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_NUMENG_RANGE = 0;
    }
   
    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) ) {
        std::cout << " APCEMM cannot accept multiple cases when the 'parameter sweep?' argument is turned off! Aborting." << std::endl;
        exit(1);
    }
    
    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( tokens.size() > 2 ) {
            std::cout << " Wrong input for " << variable << " when MC is turned on!" << std::endl;
            std::cout << " Expected format is min max or min:max representing the range of possible values" << std::endl;
            exit(1);
        } else if ( tokens.size() == 2 ) {
            Input_Opt.PARAMETER_NUMENG_RANGE = 1;
            sort(tokens.begin(), tokens.end());
        } else if ( tokens.size() == 1 )
            Input_Opt.PARAMETER_NUMENG_RANGE = 0;
    }
    
    /* Find unit in between "[" and "]" */ 
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_NUMENG_UNIT.assign( unit );

    if ( ( Input_Opt.PARAMETER_FILEINPUT ) && ( inumEngines != -1 ) ) {
        /* Found number of engines in parameter input file */
        Input_Opt.PARAMETER_NUMENG_RANGE = 0;
        if ( fileInput[27].size() > 0 ) {
            for ( unsigned int i = 0; i < fileInput[27].size(); i++ ) {
                if ( ( fileInput[27][i] == 2 ) || ( fileInput[27][i] == 4 ) )
                    Input_Opt.PARAMETER_NUMENG.push_back( fileInput[27][i] );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            }
        } else {
            try {
                value = std::stod( tokens[0] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_NUMENG.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[0] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    } else {
        /* Store in values for variable */
        for ( unsigned int i = 0; i < tokens.size(); i++ ) {
            try {
                value = std::stod( tokens[i] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_NUMENG.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    }

    /* ==================================================== */
    /* Wingspan                                             */
    /* ==================================================== */

    /* Variable */
    variable = "Wingspan";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Get line past the delimiter */
    subline = line.substr(FIRSTCOL);
    /* Look for colon */
    found = subline.find( COLON );
    if ( found != std::string::npos ) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

        if ( !Input_Opt.SIMULATION_MONTECARLO ) {
            if (tokens.size() != 3) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: " << std::endl;
                std::cout << "   --> begin:step:end" << std::endl;
                std::cout << "       or" << std::endl;
                std::cout << "   --> val1 val2 val3 ..." << std::endl;
                exit(1);
            }
            if ( ( std::stod(tokens[2]) <  std::stod(tokens[0]) ) || \
                 ( std::stod(tokens[1]) <= 0.0E+00 )              || \
                 ( std::stod(tokens[1]) >  ( std::stod(tokens[2]) - std::stod(tokens[0]) ) ) ) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: begin:step:end with begin < end, 0 < step < end - begin" << std::endl;
                exit(1);
            }
        }
        Input_Opt.PARAMETER_WINGSPAN_RANGE = 1;
    } 
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for " << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_WINGSPAN_RANGE = 0;
    }
   
    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) ) {
        std::cout << " APCEMM cannot accept multiple cases when the 'parameter sweep?' argument is turned off! Aborting." << std::endl;
        exit(1);
    }
    
    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( tokens.size() > 2 ) {
            std::cout << " Wrong input for " << variable << " when MC is turned on!" << std::endl;
            std::cout << " Expected format is min max or min:max representing the range of possible values" << std::endl;
            exit(1);
        } else if ( tokens.size() == 2 ) {
            Input_Opt.PARAMETER_WINGSPAN_RANGE = 1;
            sort(tokens.begin(), tokens.end());
        } else if ( tokens.size() == 1 )
            Input_Opt.PARAMETER_WINGSPAN_RANGE = 0;
    }
    
    /* Find unit in between "[" and "]" */ 
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_WINGSPAN_UNIT.assign( unit );

    if ( ( Input_Opt.PARAMETER_FILEINPUT ) && ( iwingspan != -1 ) ) {
        /* Found wingspan in parameter input file */
        Input_Opt.PARAMETER_WINGSPAN_RANGE = 0;
        if ( fileInput[28].size() > 0 ) {
            for ( unsigned int i = 0; i < fileInput[28].size(); i++ ) {
                if ( fileInput[28][i] > 0.0E+00 )
                    Input_Opt.PARAMETER_WINGSPAN.push_back( fileInput[28][i] );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            }
        } else {
            try {
                value = std::stod( tokens[0] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_WINGSPAN.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[0] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    } else {
        /* Store in values for variable */
        for ( unsigned int i = 0; i < tokens.size(); i++ ) {
            try {
                value = std::stod( tokens[i] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_WINGSPAN.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    }

    /* ==================================================== */
    /* Core exit temperature                                */
    /* ==================================================== */

    /* Variable */
    variable = "Core exit temperature";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Get line past the delimiter */
    subline = line.substr(FIRSTCOL);
    /* Look for colon */
    found = subline.find( COLON );
    if ( found != std::string::npos ) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

        if ( !Input_Opt.SIMULATION_MONTECARLO ) {
            if (tokens.size() != 3) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: " << std::endl;
                std::cout << "   --> begin:step:end" << std::endl;
                std::cout << "       or" << std::endl;
                std::cout << "   --> val1 val2 val3 ..." << std::endl;
                exit(1);
            }
            if ( ( std::stod(tokens[2]) <  std::stod(tokens[0]) ) || \
                 ( std::stod(tokens[1]) <= 0.0E+00 )              || \
                 ( std::stod(tokens[1]) >  ( std::stod(tokens[2]) - std::stod(tokens[0]) ) ) ) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: begin:step:end with begin < end, 0 < step < end - begin" << std::endl;
                exit(1);
            }
        }
        Input_Opt.PARAMETER_COREEXITTEMP_RANGE = 1;
    } 
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for " << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_COREEXITTEMP_RANGE = 0;
    }
   
    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) ) {
        std::cout << " APCEMM cannot accept multiple cases when the 'parameter sweep?' argument is turned off! Aborting." << std::endl;
        exit(1);
    }
    
    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( tokens.size() > 2 ) {
            std::cout << " Wrong input for " << variable << " when MC is turned on!" << std::endl;
            std::cout << " Expected format is min max or min:max representing the range of possible values" << std::endl;
            exit(1);
        } else if ( tokens.size() == 2 ) {
            Input_Opt.PARAMETER_COREEXITTEMP_RANGE = 1;
            sort(tokens.begin(), tokens.end());
        } else if ( tokens.size() == 1 )
            Input_Opt.PARAMETER_COREEXITTEMP_RANGE = 0;
    }
    
    /* Find unit in between "[" and "]" */ 
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_COREEXITTEMP_UNIT.assign( unit );

    if ( ( Input_Opt.PARAMETER_FILEINPUT ) && ( icoreExitTemp != -1 ) ) {
        /* Found core exit temperature in parameter input file */
        Input_Opt.PARAMETER_COREEXITTEMP_RANGE = 0;
        if ( fileInput[29].size() > 0 ) {
            for ( unsigned int i = 0; i < fileInput[29].size(); i++ ) {
                if ( fileInput[29][i] > 0.0E+00 )
                    Input_Opt.PARAMETER_COREEXITTEMP.push_back( fileInput[29][i] );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            }
        } else {
            try {
                value = std::stod( tokens[0] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_COREEXITTEMP.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[0] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    } else {
        /* Store in values for variable */
        for ( unsigned int i = 0; i < tokens.size(); i++ ) {
            try {
                value = std::stod( tokens[i] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_COREEXITTEMP.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    }

    /* ==================================================== */
    /* Exit bypass area                                     */
    /* ==================================================== */

    /* Variable */
    variable = "Exit bypass area";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Get line past the delimiter */
    subline = line.substr(FIRSTCOL);
    /* Look for colon */
    found = subline.find( COLON );
    if ( found != std::string::npos ) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

        if ( !Input_Opt.SIMULATION_MONTECARLO ) {
            if (tokens.size() != 3) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: " << std::endl;
                std::cout << "   --> begin:step:end" << std::endl;
                std::cout << "       or" << std::endl;
                std::cout << "   --> val1 val2 val3 ..." << std::endl;
                exit(1);
            }
            if ( ( std::stod(tokens[2]) <  std::stod(tokens[0]) ) || \
                 ( std::stod(tokens[1]) <= 0.0E+00 )              || \
                 ( std::stod(tokens[1]) >  ( std::stod(tokens[2]) - std::stod(tokens[0]) ) ) ) {
                std::cout << " Wrong input for " << variable << std::endl;
                std::cout << " Expected format is: begin:step:end with begin < end, 0 < step < end - begin" << std::endl;
                exit(1);
            }
        }
        Input_Opt.PARAMETER_BYPASSAREA_RANGE = 1;
    } 
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for " << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_BYPASSAREA_RANGE = 0;
    }
   
    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) ) {
        std::cout << " APCEMM cannot accept multiple cases when the 'parameter sweep?' argument is turned off! Aborting." << std::endl;
        exit(1);
    }
    
    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( tokens.size() > 2 ) {
            std::cout << " Wrong input for " << variable << " when MC is turned on!" << std::endl;
            std::cout << " Expected format is min max or min:max representing the range of possible values" << std::endl;
            exit(1);
        } else if ( tokens.size() == 2 ) {
            Input_Opt.PARAMETER_BYPASSAREA_RANGE = 1;
            sort(tokens.begin(), tokens.end());
        } else if ( tokens.size() == 1 )
            Input_Opt.PARAMETER_BYPASSAREA_RANGE = 0;
    }
    
    /* Find unit in between "[" and "]" */ 
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_BYPASSAREA_UNIT.assign( unit );

    if ( ( Input_Opt.PARAMETER_FILEINPUT ) && ( ibypassArea != -1 ) ) {
        /* Found bypass area in parameter input file */
        Input_Opt.PARAMETER_BYPASSAREA_RANGE = 0;
        if ( fileInput[29].size() > 0 ) {
            for ( unsigned int i = 0; i < fileInput[29].size(); i++ ) {
                if ( fileInput[29][i] > 0.0E+00 )
                    Input_Opt.PARAMETER_BYPASSAREA.push_back( fileInput[29][i] );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            }
        } else {
            try {
                value = std::stod( tokens[0] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_BYPASSAREA.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[0] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    } else {
        /* Store in values for variable */
        for ( unsigned int i = 0; i < tokens.size(); i++ ) {
            try {
                value = std::stod( tokens[i] );
                if ( value >= 0.0E+00 )
                    Input_Opt.PARAMETER_BYPASSAREA.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << variable << " needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
                exit(1);
            }
        }
    }

    /* If we found NOx flow, convert to EI */
    if ( ( Input_Opt.PARAMETER_FILEINPUT ) && ( iNOxFlow != -1 ) ) {
        if ( Input_Opt.PARAMETER_FF.size() == Input_Opt.PARAMETER_EI_NOX.size() ) {
            for ( unsigned int i = 0; i < Input_Opt.PARAMETER_EI_NOX.size(); i++ )
                Input_Opt.PARAMETER_EI_NOX[i] /= Input_Opt.PARAMETER_FF[i];
        } else {
            for ( unsigned int i = 0; i < Input_Opt.PARAMETER_EI_NOX.size(); i++ )
                Input_Opt.PARAMETER_EI_NOX[i] /= Input_Opt.PARAMETER_FF[0];
        }
        /* Do we have NOx flow in [g/s] or [g/m]?
         * If [g/m], then convert to [g/s] assuming average flight speed */
        for ( unsigned int i = 0; i < Input_Opt.PARAMETER_EI_NOX.size(); i++ )
            Input_Opt.PARAMETER_EI_NOX[i] *= 250.0;
    }


    /* Return success */
    RC = SUCCESS;

    /* ==================================================== */
    /* Print to screen                                      */
    /* ==================================================== */

    std::cout << " %%% PARAMETER SWEEP %%% :"                                                        << std::endl;
    std::cout << " ------------------------+------------------------------------------------------ " << std::endl;

    std::cout << " Input file?             : " << Input_Opt.PARAMETER_FILEINPUT << std::endl;
    std::cout << "  => File name           : " << Input_Opt.PARAMETER_FILENAME << std::endl;

    std::cout << " Parameter entries --->  : " << std::endl;

    std::cout << "  => Plume process. [" << Input_Opt.PARAMETER_PLUMEPROCESS_UNIT << "] : ";
    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( Input_Opt.PARAMETER_PLUMEPROCESS_RANGE )
            std::cout << "[" << Input_Opt.PARAMETER_PLUMEPROCESS[0] << "," << Input_Opt.PARAMETER_PLUMEPROCESS[1] << "]" << std::endl;
        else
            std::cout << Input_Opt.PARAMETER_PLUMEPROCESS[0] << std::endl;
    } else {
        if ( Input_Opt.PARAMETER_PLUMEPROCESS_RANGE )
            std::cout << Input_Opt.PARAMETER_PLUMEPROCESS[0] << ":" << Input_Opt.PARAMETER_PLUMEPROCESS[1] << ":" << Input_Opt.PARAMETER_PLUMEPROCESS[2] << std::endl;
        else {
            for ( unsigned int i = 0; i < Input_Opt.PARAMETER_PLUMEPROCESS.size(); i++ )
                std::cout << Input_Opt.PARAMETER_PLUMEPROCESS[i] << " ";
            std::cout << std::endl;
        }
    }
    /* ---- Meteorological parameters --------------------- */
    std::cout << " Meteorological parameter:" << std::endl;
    std::cout << "  => Temperature [" << Input_Opt.PARAMETER_TEMPERATURE_UNIT << "]     : ";
    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( Input_Opt.PARAMETER_TEMPERATURE_RANGE )
            std::cout << "[" << Input_Opt.PARAMETER_TEMPERATURE[0] << "," << Input_Opt.PARAMETER_TEMPERATURE[1] << "]" << std::endl;
        else
            std::cout << Input_Opt.PARAMETER_TEMPERATURE[0] << std::endl;
    } else {
        if ( Input_Opt.PARAMETER_TEMPERATURE_RANGE )
            std::cout << Input_Opt.PARAMETER_TEMPERATURE[0] << ":" << Input_Opt.PARAMETER_TEMPERATURE[1] << ":" << Input_Opt.PARAMETER_TEMPERATURE[2] << std::endl;
        else {
            for ( unsigned int i = 0; i < Input_Opt.PARAMETER_TEMPERATURE.size(); i++ )
                std::cout << Input_Opt.PARAMETER_TEMPERATURE[i] << " ";
            std::cout << std::endl;
        }
    }
    std::cout << "  => R.Hum wrt water [" << Input_Opt.PARAMETER_RHW_UNIT << "] : ";
    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( Input_Opt.PARAMETER_RHW_RANGE )
            std::cout << "[" << Input_Opt.PARAMETER_RHW[0] << "," << Input_Opt.PARAMETER_RHW[1] << "]" << std::endl;
        else
            std::cout << Input_Opt.PARAMETER_RHW[0] << std::endl;
    } else {
        if ( Input_Opt.PARAMETER_RHW_RANGE )
            std::cout << Input_Opt.PARAMETER_RHW[0] << ":" << Input_Opt.PARAMETER_RHW[1] << ":" << Input_Opt.PARAMETER_RHW[2] << std::endl;
        else {
            for ( unsigned int i = 0; i < Input_Opt.PARAMETER_RHW.size(); i++ )
                std::cout << Input_Opt.PARAMETER_RHW[i] << " ";
            std::cout << std::endl;
        }
    }
    std::cout << "  => Horiz. diff. [" << Input_Opt.PARAMETER_DH_UNIT << "]: ";
    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( Input_Opt.PARAMETER_DH_RANGE )
            std::cout << "[" << Input_Opt.PARAMETER_DH[0] << "," << Input_Opt.PARAMETER_DH[1] << "]" << std::endl;
        else
            std::cout << Input_Opt.PARAMETER_DH[0] << std::endl;
    } else {
        if ( Input_Opt.PARAMETER_DH_RANGE )
            std::cout << Input_Opt.PARAMETER_DH[0] << ":" << Input_Opt.PARAMETER_DH[1] << ":" << Input_Opt.PARAMETER_DH[2] << std::endl;
        else {
            for ( unsigned int i = 0; i < Input_Opt.PARAMETER_DH.size(); i++ )
                std::cout << Input_Opt.PARAMETER_DH[i] << " ";
            std::cout << std::endl;
        }
    }
    std::cout << "  => Verti. diff. [" << Input_Opt.PARAMETER_DV_UNIT << "]: ";
    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( Input_Opt.PARAMETER_DV_RANGE )
            std::cout << "[" << Input_Opt.PARAMETER_DV[0] << "," << Input_Opt.PARAMETER_DV[1] << "]" << std::endl;
        else
            std::cout << Input_Opt.PARAMETER_DV[0] << std::endl;
    } else {
        if ( Input_Opt.PARAMETER_DV_RANGE )
            std::cout << Input_Opt.PARAMETER_DV[0] << ":" << Input_Opt.PARAMETER_DV[1] << ":" << Input_Opt.PARAMETER_DV[2] << std::endl;
        else {
            for ( unsigned int i = 0; i < Input_Opt.PARAMETER_DV.size(); i++ )
                std::cout << Input_Opt.PARAMETER_DV[i] << " ";
            std::cout << std::endl;
        }
    }
    std::cout << "  => Shear [" << Input_Opt.PARAMETER_SHEAR_UNIT << "]         : ";
    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( Input_Opt.PARAMETER_SHEAR_RANGE )
            std::cout << "[" << Input_Opt.PARAMETER_SHEAR[0] << "," << Input_Opt.PARAMETER_SHEAR[1] << "]" << std::endl;
        else
            std::cout << Input_Opt.PARAMETER_SHEAR[0] << std::endl;
    } else {
        if ( Input_Opt.PARAMETER_SHEAR_RANGE )
            std::cout << Input_Opt.PARAMETER_SHEAR[0] << ":" << Input_Opt.PARAMETER_SHEAR[1] << ":" << Input_Opt.PARAMETER_SHEAR[2] << std::endl;
        else {
            for ( unsigned int i = 0; i < Input_Opt.PARAMETER_SHEAR.size(); i++ )
                std::cout << Input_Opt.PARAMETER_SHEAR[i] << " ";
            std::cout << std::endl;
        }
    }

    /* ---- Geographical location ------------------------- */
    std::cout << " Geographical parameters :" << std::endl;
    std::cout << "  => LON [" << Input_Opt.PARAMETER_LONGITUDE_UNIT << "]           : ";
    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( Input_Opt.PARAMETER_LONGITUDE_RANGE )
            std::cout << "[" << Input_Opt.PARAMETER_LONGITUDE[0] << "," << Input_Opt.PARAMETER_LONGITUDE[1] << "]" << std::endl;
        else
            std::cout << Input_Opt.PARAMETER_LONGITUDE[0] << std::endl;
    } else {
        if ( Input_Opt.PARAMETER_LONGITUDE_RANGE )
            std::cout << Input_Opt.PARAMETER_LONGITUDE[0] << ":" << Input_Opt.PARAMETER_LONGITUDE[1] << ":" << Input_Opt.PARAMETER_LONGITUDE[2] << std::endl;
        else {
            for ( unsigned int i = 0; i < Input_Opt.PARAMETER_LONGITUDE.size(); i++ )
                std::cout << Input_Opt.PARAMETER_LONGITUDE[i] << " ";
            std::cout << std::endl;
        }
    }
    std::cout << "  => LAT [" << Input_Opt.PARAMETER_LATITUDE_UNIT << "]           : ";
    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( Input_Opt.PARAMETER_LATITUDE_RANGE )
            std::cout << "[" << Input_Opt.PARAMETER_LATITUDE[0] << "," << Input_Opt.PARAMETER_LATITUDE[1] << "]" << std::endl;
        else
            std::cout << Input_Opt.PARAMETER_LATITUDE[0] << std::endl;
    } else {
        if ( Input_Opt.PARAMETER_LATITUDE_RANGE )
            std::cout << Input_Opt.PARAMETER_LATITUDE[0] << ":" << Input_Opt.PARAMETER_LATITUDE[1] << ":" << Input_Opt.PARAMETER_LATITUDE[2] << std::endl;
        else {
            for ( unsigned int i = 0; i < Input_Opt.PARAMETER_LATITUDE.size(); i++ )
                std::cout << Input_Opt.PARAMETER_LATITUDE[i] << " ";
            std::cout << std::endl;
        }
    }
    std::cout << "  => Pressure [" << Input_Opt.PARAMETER_PRESSURE_UNIT << "]      : ";
    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( Input_Opt.PARAMETER_PRESSURE_RANGE )
            std::cout << "[" << Input_Opt.PARAMETER_PRESSURE[0] << "," << Input_Opt.PARAMETER_PRESSURE[1] << "]" << std::endl;
        else
            std::cout << Input_Opt.PARAMETER_PRESSURE[0] << std::endl;
    } else {
        if ( Input_Opt.PARAMETER_PRESSURE_RANGE )
            std::cout << Input_Opt.PARAMETER_PRESSURE[0] << ":" << Input_Opt.PARAMETER_PRESSURE[1] << ":" << Input_Opt.PARAMETER_PRESSURE[2] << std::endl;
        else {
            for ( unsigned int i = 0; i < Input_Opt.PARAMETER_PRESSURE.size(); i++ )
                std::cout << Input_Opt.PARAMETER_PRESSURE[i] << " ";
            std::cout << std::endl;
        }
    }

    /* ---- Time of emission ------------------------------ */
    std::cout << " Time of emission        :" << std::endl;
    std::cout << "  => Emission day [" << Input_Opt.PARAMETER_EDAY_UNIT << "]: ";
    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( Input_Opt.PARAMETER_EDAY_RANGE )
            std::cout << "[" << Input_Opt.PARAMETER_EDAY[0] << "," << Input_Opt.PARAMETER_EDAY[1] << "]" << std::endl;
        else
            std::cout << Input_Opt.PARAMETER_EDAY[0] << std::endl;
    } else {
        if ( Input_Opt.PARAMETER_EDAY_RANGE )
            std::cout << Input_Opt.PARAMETER_EDAY[0] << ":" << Input_Opt.PARAMETER_EDAY[1] << ":" << Input_Opt.PARAMETER_EDAY[2] << std::endl;
        else {
            for ( unsigned int i = 0; i < Input_Opt.PARAMETER_EDAY.size(); i++ )
                std::cout << Input_Opt.PARAMETER_EDAY[i] << " ";
            std::cout << std::endl;
        }
    }
    std::cout << "  => Emission time [" << Input_Opt.PARAMETER_ETIME_UNIT << "]  : ";
    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( Input_Opt.PARAMETER_ETIME_RANGE )
            std::cout << "[" << Input_Opt.PARAMETER_ETIME[0] << "," << Input_Opt.PARAMETER_ETIME[1] << "]" << std::endl;
        else
            std::cout << Input_Opt.PARAMETER_ETIME[0] << std::endl;
    } else {
        if ( Input_Opt.PARAMETER_ETIME_RANGE )
            std::cout << Input_Opt.PARAMETER_ETIME[0] << ":" << Input_Opt.PARAMETER_ETIME[1] << ":" << Input_Opt.PARAMETER_ETIME[2] << std::endl;
        else {
            for ( unsigned int i = 0; i < Input_Opt.PARAMETER_ETIME.size(); i++ )
                std::cout << Input_Opt.PARAMETER_ETIME[i] << " ";
            std::cout << std::endl;
        }
    }

    /* ---- Background mixing ratios ---------------------- */
    std::cout << " Background mix. ratios  :" << std::endl;
    std::cout << "  => NOx  [" << Input_Opt.PARAMETER_BACKG_NOX_UNIT << "]          : ";
    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( Input_Opt.PARAMETER_BACKG_NOX_RANGE )
            std::cout << "[" << Input_Opt.PARAMETER_BACKG_NOX[0] << "," << Input_Opt.PARAMETER_BACKG_NOX[1] << "]" << std::endl;
        else
            std::cout << Input_Opt.PARAMETER_BACKG_NOX[0] << std::endl;
    } else {
        if ( Input_Opt.PARAMETER_BACKG_NOX_RANGE )
            std::cout << Input_Opt.PARAMETER_BACKG_NOX[0] << ":" << Input_Opt.PARAMETER_BACKG_NOX[1] << ":" << Input_Opt.PARAMETER_BACKG_NOX[2] << std::endl;
        else {
            for ( unsigned int i = 0; i < Input_Opt.PARAMETER_BACKG_NOX.size(); i++ )
                std::cout << Input_Opt.PARAMETER_BACKG_NOX[i] << " ";
            std::cout << std::endl;
        }
    }
    std::cout << "  => HNO3 [" << Input_Opt.PARAMETER_BACKG_HNO3_UNIT << "]          : ";
    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( Input_Opt.PARAMETER_BACKG_HNO3_RANGE )
            std::cout << "[" << Input_Opt.PARAMETER_BACKG_HNO3[0] << "," << Input_Opt.PARAMETER_BACKG_HNO3[1] << "]" << std::endl;
        else
            std::cout << Input_Opt.PARAMETER_BACKG_HNO3[0] << std::endl;
    } else {
        if ( Input_Opt.PARAMETER_BACKG_HNO3_RANGE )
            std::cout << Input_Opt.PARAMETER_BACKG_HNO3[0] << ":" << Input_Opt.PARAMETER_BACKG_HNO3[1] << ":" << Input_Opt.PARAMETER_BACKG_HNO3[2] << std::endl;
        else {
            for ( unsigned int i = 0; i < Input_Opt.PARAMETER_BACKG_HNO3.size(); i++ )
                std::cout << Input_Opt.PARAMETER_BACKG_HNO3[i] << " ";
            std::cout << std::endl;
        }
    }
    std::cout << "  => O3   [" << Input_Opt.PARAMETER_BACKG_O3_UNIT << "]          : ";
    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( Input_Opt.PARAMETER_BACKG_O3_RANGE )
            std::cout << "[" << Input_Opt.PARAMETER_BACKG_O3[0] << "," << Input_Opt.PARAMETER_BACKG_O3[1] << "]" << std::endl;
        else
            std::cout << Input_Opt.PARAMETER_BACKG_O3[0] << std::endl;
    } else {
        if ( Input_Opt.PARAMETER_BACKG_O3_RANGE )
            std::cout << Input_Opt.PARAMETER_BACKG_O3[0] << ":" << Input_Opt.PARAMETER_BACKG_O3[1] << ":" << Input_Opt.PARAMETER_BACKG_O3[2] << std::endl;
        else {
            for ( unsigned int i = 0; i < Input_Opt.PARAMETER_BACKG_O3.size(); i++ )
                std::cout << Input_Opt.PARAMETER_BACKG_O3[i] << " ";
            std::cout << std::endl;
        }
    }
    std::cout << "  => CO   [" << Input_Opt.PARAMETER_BACKG_CO_UNIT << "]          : ";
    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( Input_Opt.PARAMETER_BACKG_CO_RANGE )
            std::cout << "[" << Input_Opt.PARAMETER_BACKG_CO[0] << "," << Input_Opt.PARAMETER_BACKG_CO[1] << "]" << std::endl;
        else
            std::cout << Input_Opt.PARAMETER_BACKG_CO[0] << std::endl;
    } else {
        if ( Input_Opt.PARAMETER_BACKG_CO_RANGE )
            std::cout << Input_Opt.PARAMETER_BACKG_CO[0] << ":" << Input_Opt.PARAMETER_BACKG_CO[1] << ":" << Input_Opt.PARAMETER_BACKG_CO[2] << std::endl;
        else {
            for ( unsigned int i = 0; i < Input_Opt.PARAMETER_BACKG_CO.size(); i++ )
                std::cout << Input_Opt.PARAMETER_BACKG_CO[i] << " ";
            std::cout << std::endl;
        }
    }
    std::cout << "  => CH4  [" << Input_Opt.PARAMETER_BACKG_CH4_UNIT << "]          : ";
    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( Input_Opt.PARAMETER_BACKG_CH4_RANGE )
            std::cout << "[" << Input_Opt.PARAMETER_BACKG_CH4[0] << "," << Input_Opt.PARAMETER_BACKG_CH4[1] << "]" << std::endl;
        else
            std::cout << Input_Opt.PARAMETER_BACKG_CH4[0] << std::endl;
    } else {
        if ( Input_Opt.PARAMETER_BACKG_CH4_RANGE )
            std::cout << Input_Opt.PARAMETER_BACKG_CH4[0] << ":" << Input_Opt.PARAMETER_BACKG_CH4[1] << ":" << Input_Opt.PARAMETER_BACKG_CH4[2] << std::endl;
        else {
            for ( unsigned int i = 0; i < Input_Opt.PARAMETER_BACKG_CH4.size(); i++ )
                std::cout << Input_Opt.PARAMETER_BACKG_CH4[i] << " ";
            std::cout << std::endl;
        }
    }
    std::cout << "  => SO2  [" << Input_Opt.PARAMETER_BACKG_SO2_UNIT << "]          : ";
    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( Input_Opt.PARAMETER_BACKG_SO2_RANGE )
            std::cout << "[" << Input_Opt.PARAMETER_BACKG_SO2[0] << "," << Input_Opt.PARAMETER_BACKG_SO2[1] << "]" << std::endl;
        else
            std::cout << Input_Opt.PARAMETER_BACKG_SO2[0] << std::endl;
    } else {
        if ( Input_Opt.PARAMETER_BACKG_SO2_RANGE )
            std::cout << Input_Opt.PARAMETER_BACKG_SO2[0] << ":" << Input_Opt.PARAMETER_BACKG_SO2[1] << ":" << Input_Opt.PARAMETER_BACKG_SO2[2] << std::endl;
        else {
            for ( unsigned int i = 0; i < Input_Opt.PARAMETER_BACKG_SO2.size(); i++ )
                std::cout << Input_Opt.PARAMETER_BACKG_SO2[i] << " ";
            std::cout << std::endl;
        }
    }

    /* ---- Emission indices ------------------------------ */
    std::cout << " Emission indices        :" << std::endl;
    std::cout << "  => NOx  [" << Input_Opt.PARAMETER_EI_NOX_UNIT << "]    : ";
    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( Input_Opt.PARAMETER_EI_NOX_RANGE )
            std::cout << "[" << Input_Opt.PARAMETER_EI_NOX[0] << "," << Input_Opt.PARAMETER_EI_NOX[1] << "]" << std::endl;
        else
            std::cout << Input_Opt.PARAMETER_EI_NOX[0] << std::endl;
    } else {
        if ( Input_Opt.PARAMETER_EI_NOX_RANGE )
            std::cout << Input_Opt.PARAMETER_EI_NOX[0] << ":" << Input_Opt.PARAMETER_EI_NOX[1] << ":" << Input_Opt.PARAMETER_EI_NOX[2] << std::endl;
        else {
            for ( unsigned int i = 0; i < Input_Opt.PARAMETER_EI_NOX.size(); i++ )
                std::cout << Input_Opt.PARAMETER_EI_NOX[i] << " ";
            std::cout << std::endl;
        }
    }
    std::cout << "  => CO   [" << Input_Opt.PARAMETER_EI_CO_UNIT << "]    : ";
    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( Input_Opt.PARAMETER_EI_CO_RANGE )
            std::cout << "[" << Input_Opt.PARAMETER_EI_CO[0] << "," << Input_Opt.PARAMETER_EI_CO[1] << "]" << std::endl;
        else
            std::cout << Input_Opt.PARAMETER_EI_CO[0] << std::endl;
    } else {
        if ( Input_Opt.PARAMETER_EI_CO_RANGE )
            std::cout << Input_Opt.PARAMETER_EI_CO[0] << ":" << Input_Opt.PARAMETER_EI_CO[1] << ":" << Input_Opt.PARAMETER_EI_CO[2] << std::endl;
        else {
            for ( unsigned int i = 0; i < Input_Opt.PARAMETER_EI_CO.size(); i++ )
                std::cout << Input_Opt.PARAMETER_EI_CO[i] << " ";
            std::cout << std::endl;
        }
    }
    std::cout << "  => UHC  [" << Input_Opt.PARAMETER_EI_UHC_UNIT << "]    : ";
    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( Input_Opt.PARAMETER_EI_UHC_RANGE )
            std::cout << "[" << Input_Opt.PARAMETER_EI_UHC[0] << "," << Input_Opt.PARAMETER_EI_UHC[1] << "]" << std::endl;
        else
            std::cout << Input_Opt.PARAMETER_EI_UHC[0] << std::endl;
    } else {
        if ( Input_Opt.PARAMETER_EI_UHC_RANGE )
            std::cout << Input_Opt.PARAMETER_EI_UHC[0] << ":" << Input_Opt.PARAMETER_EI_UHC[1] << ":" << Input_Opt.PARAMETER_EI_UHC[2] << std::endl;
        else {
            for ( unsigned int i = 0; i < Input_Opt.PARAMETER_EI_UHC.size(); i++ )
                std::cout << Input_Opt.PARAMETER_EI_UHC[i] << " ";
            std::cout << std::endl;
        }
    }
    std::cout << "  => SO2  [" << Input_Opt.PARAMETER_EI_SO2_UNIT << "]    : ";
    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( Input_Opt.PARAMETER_EI_SO2_RANGE )
            std::cout << "[" << Input_Opt.PARAMETER_EI_SO2[0] << "," << Input_Opt.PARAMETER_EI_SO2[1] << "]" << std::endl;
        else
            std::cout << Input_Opt.PARAMETER_EI_SO2[0] << std::endl;
    } else {
        if ( Input_Opt.PARAMETER_EI_SO2_RANGE )
            std::cout << Input_Opt.PARAMETER_EI_SO2[0] << ":" << Input_Opt.PARAMETER_EI_SO2[1] << ":" << Input_Opt.PARAMETER_EI_SO2[2] << std::endl;
        else {
            for ( unsigned int i = 0; i < Input_Opt.PARAMETER_EI_SO2.size(); i++ )
                std::cout << Input_Opt.PARAMETER_EI_SO2[i] << " ";
            std::cout << std::endl;
        }
    }
    std::cout << "  => SO2 to SO4 conv [" << Input_Opt.PARAMETER_EI_SO2TOSO4_UNIT << "] : ";
    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( Input_Opt.PARAMETER_EI_SO2TOSO4_RANGE )
            std::cout << "[" << Input_Opt.PARAMETER_EI_SO2TOSO4[0] << "," << Input_Opt.PARAMETER_EI_SO2TOSO4[1] << "]" << std::endl;
        else
            std::cout << Input_Opt.PARAMETER_EI_SO2TOSO4[0] << std::endl;
    } else {
        if ( Input_Opt.PARAMETER_EI_SO2TOSO4_RANGE )
            std::cout << Input_Opt.PARAMETER_EI_SO2TOSO4[0] << ":" << Input_Opt.PARAMETER_EI_SO2TOSO4[1] << ":" << Input_Opt.PARAMETER_EI_SO2TOSO4[2] << std::endl;
        else {
            for ( unsigned int i = 0; i < Input_Opt.PARAMETER_EI_SO2TOSO4.size(); i++ )
                std::cout << Input_Opt.PARAMETER_EI_SO2TOSO4[i] << " ";
            std::cout << std::endl;
        }
    }
    std::cout << "  => Soot [" << Input_Opt.PARAMETER_EI_SOOT_UNIT << "]    : ";
    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( Input_Opt.PARAMETER_EI_SOOT_RANGE )
            std::cout << "[" << Input_Opt.PARAMETER_EI_SOOT[0] << "," << Input_Opt.PARAMETER_EI_SOOT[1] << "]" << std::endl;
        else
            std::cout << Input_Opt.PARAMETER_EI_SOOT[0] << std::endl;
    } else {
        if ( Input_Opt.PARAMETER_EI_SOOT_RANGE )
            std::cout << Input_Opt.PARAMETER_EI_SOOT[0] << ":" << Input_Opt.PARAMETER_EI_SOOT[1] << ":" << Input_Opt.PARAMETER_EI_SOOT[2] << std::endl;
        else {
            for ( unsigned int i = 0; i < Input_Opt.PARAMETER_EI_SOOT.size(); i++ )
                std::cout << Input_Opt.PARAMETER_EI_SOOT[i] << " ";
            std::cout << std::endl;
        }
    }
    std::cout << "  => Soot Radius [" << Input_Opt.PARAMETER_EI_SOOTRAD_UNIT << "]     : ";
    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( Input_Opt.PARAMETER_EI_SOOTRAD_RANGE )
            std::cout << "[" << Input_Opt.PARAMETER_EI_SOOTRAD[0] << "," << Input_Opt.PARAMETER_EI_SOOTRAD[1] << "]" << std::endl;
        else
            std::cout << Input_Opt.PARAMETER_EI_SOOTRAD[0] << std::endl;
    } else {
        if ( Input_Opt.PARAMETER_EI_SOOTRAD_RANGE )
            std::cout << Input_Opt.PARAMETER_EI_SOOTRAD[0] << ":" << Input_Opt.PARAMETER_EI_SOOTRAD[1] << ":" << Input_Opt.PARAMETER_EI_SOOTRAD[2] << std::endl;
        else {
            for ( unsigned int i = 0; i < Input_Opt.PARAMETER_EI_SOOTRAD.size(); i++ )
                std::cout << Input_Opt.PARAMETER_EI_SOOTRAD[i] << " ";
            std::cout << std::endl;
        }
    }

    /* ---- Total fuel flow ------------------------------- */
    std::cout << "  Total fuel flow [" << Input_Opt.PARAMETER_FF_UNIT << "] : ";
    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( Input_Opt.PARAMETER_FF_RANGE )
            std::cout << "[" << Input_Opt.PARAMETER_FF[0] << "," << Input_Opt.PARAMETER_FF[1] << "]" << std::endl;
        else
            std::cout << Input_Opt.PARAMETER_FF[0] << std::endl;
    } else {
        if ( Input_Opt.PARAMETER_FF_RANGE )
            std::cout << Input_Opt.PARAMETER_FF[0] << ":" << Input_Opt.PARAMETER_FF[1] << ":" << Input_Opt.PARAMETER_FF[2] << std::endl;
        else {
            for ( unsigned int i = 0; i < Input_Opt.PARAMETER_FF.size(); i++ )
                std::cout << Input_Opt.PARAMETER_FF[i] << " ";
            std::cout << std::endl;
        }
    }

    /* ---- Aircraft mass ------------------------------- */
    std::cout << "  Aircraft mass [" << Input_Opt.PARAMETER_AMASS_UNIT << "] : ";
    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( Input_Opt.PARAMETER_AMASS_RANGE )
            std::cout << "[" << Input_Opt.PARAMETER_AMASS[0] << "," << Input_Opt.PARAMETER_AMASS[1] << "]" << std::endl;
        else
            std::cout << Input_Opt.PARAMETER_AMASS[0] << std::endl;
    } else {
        if ( Input_Opt.PARAMETER_AMASS_RANGE )
            std::cout << Input_Opt.PARAMETER_AMASS[0] << ":" << Input_Opt.PARAMETER_AMASS[1] << ":" << Input_Opt.PARAMETER_AMASS[2] << std::endl;
        else {
            for ( unsigned int i = 0; i < Input_Opt.PARAMETER_AMASS.size(); i++ )
                std::cout << Input_Opt.PARAMETER_AMASS[i] << " ";
            std::cout << std::endl;
        }
    }

    /* ---- Flight speed ------------------------------- */
    std::cout << "  Flight speed [" << Input_Opt.PARAMETER_AMASS_UNIT << "] : ";
    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( Input_Opt.PARAMETER_FSPEED_RANGE )
            std::cout << "[" << Input_Opt.PARAMETER_FSPEED[0] << "," << Input_Opt.PARAMETER_FSPEED[1] << "]" << std::endl;
        else
            std::cout << Input_Opt.PARAMETER_FSPEED[0] << std::endl;
    } else {
        if ( Input_Opt.PARAMETER_FSPEED_RANGE )
            std::cout << Input_Opt.PARAMETER_FSPEED[0] << ":" << Input_Opt.PARAMETER_FSPEED[1] << ":" << Input_Opt.PARAMETER_FSPEED[2] << std::endl;
        else {
            for ( unsigned int i = 0; i < Input_Opt.PARAMETER_FSPEED.size(); i++ )
                std::cout << Input_Opt.PARAMETER_FSPEED[i] << " ";
            std::cout << std::endl;
        }
    }

    /* ---- Number of engines ------------------------------- */
    std::cout << "  Num. of engines [" << Input_Opt.PARAMETER_NUMENG_UNIT << "] : ";
    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( Input_Opt.PARAMETER_NUMENG_RANGE )
            std::cout << "[" << Input_Opt.PARAMETER_NUMENG[0] << "," << Input_Opt.PARAMETER_NUMENG[1] << "]" << std::endl;
        else
            std::cout << Input_Opt.PARAMETER_NUMENG[0] << std::endl;
    } else {
        if ( Input_Opt.PARAMETER_NUMENG_RANGE )
            std::cout << Input_Opt.PARAMETER_NUMENG[0] << ":" << Input_Opt.PARAMETER_NUMENG[1] << ":" << Input_Opt.PARAMETER_NUMENG[2] << std::endl;
        else {
            for ( unsigned int i = 0; i < Input_Opt.PARAMETER_NUMENG.size(); i++ )
                std::cout << Input_Opt.PARAMETER_NUMENG[i] << " ";
            std::cout << std::endl;
        }
    }

    /* ---- Wingspan ------------------------------- */
    std::cout << "  Wingspan [" << Input_Opt.PARAMETER_WINGSPAN_UNIT << "] : ";
    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( Input_Opt.PARAMETER_WINGSPAN_RANGE )
            std::cout << "[" << Input_Opt.PARAMETER_WINGSPAN[0] << "," << Input_Opt.PARAMETER_WINGSPAN[1] << "]" << std::endl;
        else
            std::cout << Input_Opt.PARAMETER_WINGSPAN[0] << std::endl;
    } else {
        if ( Input_Opt.PARAMETER_WINGSPAN_RANGE )
            std::cout << Input_Opt.PARAMETER_WINGSPAN[0] << ":" << Input_Opt.PARAMETER_WINGSPAN[1] << ":" << Input_Opt.PARAMETER_WINGSPAN[2] << std::endl;
        else {
            for ( unsigned int i = 0; i < Input_Opt.PARAMETER_WINGSPAN.size(); i++ )
                std::cout << Input_Opt.PARAMETER_WINGSPAN[i] << " ";
            std::cout << std::endl;
        }
    }

    /* ---- Core exit temperature ------------------------------- */
    std::cout << "  Core exit temp. [" << Input_Opt.PARAMETER_COREEXITTEMP_UNIT << "] : ";
    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( Input_Opt.PARAMETER_COREEXITTEMP_RANGE )
            std::cout << "[" << Input_Opt.PARAMETER_COREEXITTEMP[0] << "," << Input_Opt.PARAMETER_COREEXITTEMP[1] << "]" << std::endl;
        else
            std::cout << Input_Opt.PARAMETER_COREEXITTEMP[0] << std::endl;
    } else {
        if ( Input_Opt.PARAMETER_COREEXITTEMP_RANGE )
            std::cout << Input_Opt.PARAMETER_COREEXITTEMP[0] << ":" << Input_Opt.PARAMETER_COREEXITTEMP[1] << ":" << Input_Opt.PARAMETER_COREEXITTEMP[2] << std::endl;
        else {
            for ( unsigned int i = 0; i < Input_Opt.PARAMETER_COREEXITTEMP.size(); i++ )
                std::cout << Input_Opt.PARAMETER_COREEXITTEMP[i] << " ";
            std::cout << std::endl;
        }
    }

    /* ---- Bypass area ------------------------------- */
    std::cout << "  Exit bypass area [" << Input_Opt.PARAMETER_BYPASSAREA_UNIT << "] : ";
    if ( Input_Opt.SIMULATION_MONTECARLO ) {
        if ( Input_Opt.PARAMETER_BYPASSAREA_RANGE )
            std::cout << "[" << Input_Opt.PARAMETER_BYPASSAREA[0] << "," << Input_Opt.PARAMETER_BYPASSAREA[1] << "]" << std::endl;
        else
            std::cout << Input_Opt.PARAMETER_BYPASSAREA[0] << std::endl;
    } else {
        if ( Input_Opt.PARAMETER_BYPASSAREA_RANGE )
            std::cout << Input_Opt.PARAMETER_BYPASSAREA[0] << ":" << Input_Opt.PARAMETER_BYPASSAREA[1] << ":" << Input_Opt.PARAMETER_BYPASSAREA[2] << std::endl;
        else {
            for ( unsigned int i = 0; i < Input_Opt.PARAMETER_BYPASSAREA.size(); i++ )
                std::cout << Input_Opt.PARAMETER_BYPASSAREA[i] << " ";
            std::cout << std::endl;
        }
    }


} /* End of Read_Parameters */

void Read_Transport_Menu( OptInput &Input_Opt, bool &RC )
{

    /* DESCRIPTION: Function Read\_Transport\_Menu reads the TRANSPORT MENU
     * section of the APCEMM input file. */

    /* INPUT/OUTPUT PARAMETERS:
     * - Input_Opt: Input options
     * - RC: Success or failure
     */

    /* Read until all lines from the menu are read */

    std::vector<std::string> tokens;
    std::string variable;
    double value;
    std::string unit;
    unsigned first, last;

    /* ==================================================== */
    /* Turn on Transport?                                   */
    /* ==================================================== */

    variable = "Turn on Transport?";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable range */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    if ( ( strcmp(tokens[0].c_str(), "T" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "t" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "TRUE" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "true" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "True" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "YES" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Y" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "y" )    == 0 ) )
        Input_Opt.TRANSPORT_TRANSPORT = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "f" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "FALSE" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "false" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "False" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "NO" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "No" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "no" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "N" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "n" )     == 0 ) )
        Input_Opt.TRANSPORT_TRANSPORT = 0;
    else {
        std::cout << " Wrong input for: " << variable << std::endl;
        exit(1);
    }

    /* ==================================================== */
    /* Fill negative values?                                */
    /* ==================================================== */

    variable = "Fill negative values?";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable range */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    if ( ( strcmp(tokens[0].c_str(), "T" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "t" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "TRUE" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "true" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "True" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "YES" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Y" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "y" )    == 0 ) )
        Input_Opt.TRANSPORT_FILL = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "f" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "FALSE" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "false" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "False" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "NO" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "No" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "no" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "N" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "n" )     == 0 ) )
        Input_Opt.TRANSPORT_FILL = 0;
    else {
        std::cout << " Wrong input for: " << variable << std::endl;
        exit(1);
    }

    /* ==================================================== */
    /* Transport timestep                                   */
    /* ==================================================== */

    variable = "Transport timestep";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    try {
        value = std::stod( tokens[0] );
        if ( Input_Opt.TRANSPORT_TRANSPORT ) {
            if ( value > 0.0E+00 )
                Input_Opt.TRANSPORT_TIMESTEP = value;
            else {
                std::cout << " Wrong input for: " << variable << std::endl;
                std::cout << " Timestep needs to be positive" << std::endl;
                exit(1);
            }
        }
    } catch(std::exception& e) {
        std::cout << " Could not convert string to double for " << variable << std::endl;
        exit(1);
    }

    /* ==================================================== */
    /* Particle flux correction                             */
    /* ==================================================== */

    variable = "Particle flux correction?";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    if ( ( strcmp(tokens[0].c_str(), "T" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "t" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "TRUE" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "true" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "True" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "YES" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Y" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "y" )    == 0 ) )
        Input_Opt.TRANSPORT_PART_FLUX = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "f" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "FALSE" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "false" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "False" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "NO" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "No" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "no" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "N" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "n" )     == 0 ) )
        Input_Opt.TRANSPORT_PART_FLUX = 0;
    else {
        std::cout << " Wrong input for: " << variable << std::endl;
        exit(1);
    }

    /* ==================================================== */
    /* Turn on plume updraft?                               */
    /* ==================================================== */

    variable = "Turn on plume updraft?";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    if ( ( strcmp(tokens[0].c_str(), "T" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "t" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "TRUE" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "true" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "True" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "YES" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Y" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "y" )    == 0 ) )
        Input_Opt.TRANSPORT_UPDRAFT = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "f" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "FALSE" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "false" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "False" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "NO" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "No" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "no" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "N" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "n" )     == 0 ) )
        Input_Opt.TRANSPORT_UPDRAFT = 0;
    else {
        std::cout << " Wrong input for: " << variable << std::endl;
        exit(1);
    }

    /* ==================================================== */
    /* Updraft timescale                                    */
    /* ==================================================== */

    variable = "Updraft timescale";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    try {
        value = std::stod( tokens[0] );
        if ( value > 0.0E+00 )
            Input_Opt.TRANSPORT_UPDRAFT_TIMESCALE = value;
        else {
            std::cout << " Wrong input for: " << variable << std::endl;
            std::cout << " Timescale needs to be positive" << std::endl;
            exit(1);
        }
    } catch(std::exception& e) {
        std::cout << " Could not convert string to double for " << variable << std::endl;
        exit(1);
    }

    /* Find unit in between "[" and "]" */
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );

    if ( unit.compare( "s" ) == 0 ) {
        /* Do nothing. Default unit */
    } else if ( unit.compare( "min" ) == 0 ) {
        /* Convert min to s */
        Input_Opt.TRANSPORT_UPDRAFT_TIMESCALE *= 60.0;
    } else if ( unit.compare( "hr" ) == 0 ) {
        /* Convert hr to s */
        Input_Opt.TRANSPORT_UPDRAFT_TIMESCALE *= 3600.0;
    } else {
        std::cout << " Unknown unit for variable 'Updraft timescale': ";
        std::cout << unit << std::endl;
        exit(1);
    }

    /* ==================================================== */
    /* Updraft velocity                                     */
    /* ==================================================== */

    variable = "Updraft velocity";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    try {
        value = std::stod( tokens[0] );
        Input_Opt.TRANSPORT_UPDRAFT_VELOCITY = value;
    } catch(std::exception& e) {
        std::cout << " Could not convert string to double for " << variable << std::endl;
        exit(1);
    }

    /* Find unit in between "[" and "]" */
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );

    if ( unit.compare( "m/s" ) == 0 ) {
        /* Do nothing. Default unit */
    } else if ( unit.compare( "cm/s" ) == 0 ) {
        /* Convert cm/s to m/s */
        Input_Opt.TRANSPORT_UPDRAFT_VELOCITY *= 1.0E-02;
    } else if ( unit.compare( "mm/s" ) == 0 ) {
        /* Convert mm/s to m/s */
        Input_Opt.TRANSPORT_UPDRAFT_VELOCITY *= 1.00E-03;
    } else {
        std::cout << " Unknown unit for variable 'Updraft velocity': ";
        std::cout << unit << std::endl;
        exit(1);
    }

    /* Return success */
    RC = SUCCESS;

    /* ==================================================== */
    /* Print to screen                                      */
    /* ==================================================== */

    std::cout << " ------------------------+------------------------------------------------------ " << std::endl;
    std::cout << " %%% TRANSPORT MENU %%%  :"                                                        << std::endl;
    std::cout << " ------------------------+------------------------------------------------------ " << std::endl;
    std::cout << " Turn on Transport?      : " << Input_Opt.TRANSPORT_TRANSPORT                      << std::endl;
    std::cout << "  => Fill Negative Values: " << Input_Opt.TRANSPORT_FILL                           << std::endl;
    std::cout << " Transport Timestep [min]: " << Input_Opt.TRANSPORT_TIMESTEP                       << std::endl;
    std::cout << " Particle flux correction: " << Input_Opt.TRANSPORT_PART_FLUX                      << std::endl;
    std::cout << " Turn on plume updraft?  : " << Input_Opt.TRANSPORT_UPDRAFT                        << std::endl;
    std::cout << "  => Updraft timescale[s]: " << Input_Opt.TRANSPORT_UPDRAFT_TIMESCALE              << std::endl;
    std::cout << "  => Updraft vel.   [m/s]: " << Input_Opt.TRANSPORT_UPDRAFT_VELOCITY               << std::endl;

} /* End of Read_Transport_Menu */

void Read_Chemistry_Menu( OptInput &Input_Opt, bool &RC )
{

    /* DESCRIPTION: Function Read\_Chemistry\_Menu reads the CHEMISTRY MENU
     * section of the APCEMM input file. */

    /* INPUT/OUTPUT PARAMETERS:
     * - Input_Opt: Input options
     * - RC: Success or failure
     */

    /* Read until all lines from the menu are read */

    std::vector<std::string> tokens;
    std::string variable;
    double value;

    /* ==================================================== */
    /* Turn on Chemistry?                                   */
    /* ==================================================== */

    variable = "Turn on Chemistry?";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable range */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    if ( ( strcmp(tokens[0].c_str(), "T" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "t" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "TRUE" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "true" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "True" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "YES" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Y" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "y" )    == 0 ) )
        Input_Opt.CHEMISTRY_CHEMISTRY = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "f" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "FALSE" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "false" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "False" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "NO" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "No" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "no" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "N" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "n" )     == 0 ) )
        Input_Opt.CHEMISTRY_CHEMISTRY = 0;
    else {
        std::cout << " Wrong input for: " << variable << std::endl;
        exit(1);
    }

    if ( ( Input_Opt.SIMULATION_ADJOINT == 1 ) && \
         ( Input_Opt.CHEMISTRY_CHEMISTRY == 0 ) ) {
        std::cout << " Cannot perform an adjoint simulation with chemistry turned off. Abort!" << std::endl;
        exit(1);
    }

    /* ==================================================== */
    /* Perform het. chem.?                                  */
    /* ==================================================== */

    variable = "Perform het. chem.?";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable range */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    if ( ( strcmp(tokens[0].c_str(), "T" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "t" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "TRUE" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "true" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "True" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "YES" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Y" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "y" )    == 0 ) )
        Input_Opt.CHEMISTRY_HETCHEM = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "f" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "FALSE" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "false" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "False" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "NO" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "No" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "no" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "N" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "n" )     == 0 ) )
        Input_Opt.CHEMISTRY_HETCHEM = 0;
    else {
        std::cout << " Wrong input for: " << variable << std::endl;
        exit(1);
    }

    /* ==================================================== */
    /* Chemistry Timestep                                   */
    /* ==================================================== */

    variable = "Chemistry Timestep";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    try {
        value = std::stod( tokens[0] );
        if ( Input_Opt.CHEMISTRY_CHEMISTRY ) {
            if ( value > 0.0E+00 )
                Input_Opt.CHEMISTRY_TIMESTEP = value;
            else {
                std::cout << " Wrong input for: " << variable << std::endl;
                std::cout << " Timestep needs to be positive" << std::endl;
                exit(1);
            }
        }
    } catch(std::exception& e) {
        std::cout << " Could not convert string to double for " << variable << std::endl;
        exit(1);
    }

    /* ==================================================== */
    /* Photolysis rates folder                              */
    /* ==================================================== */

    variable = "Photolysis rates folder";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    Input_Opt.CHEMISTRY_JRATE_FOLDER = tokens[0];

    /* Return success */
    RC = SUCCESS;

    /* ==================================================== */
    /* Print to screen                                      */
    /* ==================================================== */

    std::cout << " ------------------------+------------------------------------------------------ " << std::endl;
    std::cout << " %%% CHEMISTRY MENU %%%  :"                                                        << std::endl;
    std::cout << " ------------------------+------------------------------------------------------ " << std::endl;
    std::cout << " Turn on Chemistry?      : " << Input_Opt.CHEMISTRY_CHEMISTRY                      << std::endl;
    std::cout << " Perform het. chem.?     : " << Input_Opt.CHEMISTRY_HETCHEM                        << std::endl;
    std::cout << " Chemistry Timestep [min]: " << Input_Opt.CHEMISTRY_TIMESTEP                       << std::endl;
    std::cout << " Photolysis rates folder : " << Input_Opt.CHEMISTRY_JRATE_FOLDER                   << std::endl;

} /* End of Read_Chemistry_Menu */

void Read_Aerosol_Menu( OptInput &Input_Opt, bool &RC )
{

    /* DESCRIPTION: Function Read\_Aerosol\_Menu reads the AEROSOL MENU
     * section of the APCEMM input file. */

    /* INPUT/OUTPUT PARAMETERS:
     * - Input_Opt: Input options
     * - RC: Success or failure
     */

    /* Read until all lines from the menu are read */

    std::vector<std::string> tokens;
    std::string variable;
    double value;

    /* ==================================================== */
    /* Turn on grav. settling?                              */
    /* ==================================================== */

    variable = "Turn on grav. settling?";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable range */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    if ( ( strcmp(tokens[0].c_str(), "T" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "t" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "TRUE" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "true" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "True" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "YES" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Y" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "y" )    == 0 ) )
        Input_Opt.AEROSOL_GRAVSETTLING = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "f" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "FALSE" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "false" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "False" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "NO" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "No" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "no" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "N" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "n" )     == 0 ) )
        Input_Opt.AEROSOL_GRAVSETTLING = 0;
    else {
        std::cout << " Wrong input for: " << variable << std::endl;
        exit(1);
    }

    /* ==================================================== */
    /* Turn on solid coag.?                                 */
    /* ==================================================== */

    variable = "Turn on solid coag.?";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable range */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    if ( ( strcmp(tokens[0].c_str(), "T" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "t" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "TRUE" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "true" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "True" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "YES" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Y" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "y" )    == 0 ) )
        Input_Opt.AEROSOL_COAGULATION_SOLID = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "f" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "FALSE" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "false" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "False" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "NO" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "No" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "no" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "N" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "n" )     == 0 ) )
        Input_Opt.AEROSOL_COAGULATION_SOLID = 0;
    else {
        std::cout << " Wrong input for: " << variable << std::endl;
        exit(1);
    }

    /* ==================================================== */
    /* Turn on liquid coag.?                                */
    /* ==================================================== */

    variable = "Turn on liquid coag.?";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable range */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    if ( ( strcmp(tokens[0].c_str(), "T" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "t" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "TRUE" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "true" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "True" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "YES" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Y" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "y" )    == 0 ) )
        Input_Opt.AEROSOL_COAGULATION_LIQUID = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "f" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "FALSE" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "false" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "False" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "NO" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "No" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "no" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "N" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "n" )     == 0 ) )
        Input_Opt.AEROSOL_COAGULATION_LIQUID = 0;
    else {
        std::cout << " Wrong input for: " << variable << std::endl;
        exit(1);
    }

    /* ==================================================== */
    /* Coagulation Timestep                                 */
    /* ==================================================== */

    variable = "Coagulation Timestep";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    try {
        value = std::stod( tokens[0] );
        if ( value > 0.0E+00 )
            Input_Opt.AEROSOL_COAGULATION_TIMESTEP = value;
        else {
            std::cout << " Wrong input for: " << variable << std::endl;
            std::cout << " Timestep needs to be positive" << std::endl;
            exit(1);
        }
    } catch(std::exception& e) {
        std::cout << " Could not convert string to double for " << variable << std::endl;
        exit(1);
    }

    /* ==================================================== */
    /* Turn on ice growth?                                  */
    /* ==================================================== */

    variable = "Turn on ice growth?";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    if ( ( strcmp(tokens[0].c_str(), "T" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "t" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "TRUE" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "true" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "True" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "YES" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Y" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "y" )    == 0 ) )
        Input_Opt.AEROSOL_ICE_GROWTH = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "f" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "FALSE" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "false" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "False" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "NO" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "No" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "no" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "N" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "n" )     == 0 ) )
        Input_Opt.AEROSOL_ICE_GROWTH = 0;
    else {
        std::cout << " Wrong input for: " << variable << std::endl;
        exit(1);
    }

    /* Return success */
    RC = SUCCESS;

    /* ==================================================== */
    /* Print to screen                                      */
    /* ==================================================== */

    std::cout << " ------------------------+------------------------------------------------------ " << std::endl;
    std::cout << " %%% AEROSOL MENU %%%    :"                                                        << std::endl;
    std::cout << " ------------------------+------------------------------------------------------ " << std::endl;
    std::cout << " Turn on grav. settling? : " << Input_Opt.AEROSOL_GRAVSETTLING                     << std::endl;
    std::cout << " Turn on solid coag.?    : " << Input_Opt.AEROSOL_COAGULATION_SOLID                << std::endl;
    std::cout << " Turn on liquid coag.?   : " << Input_Opt.AEROSOL_COAGULATION_LIQUID               << std::endl;
    std::cout << "  => Coag. timestep [min]: " << Input_Opt.AEROSOL_COAGULATION_TIMESTEP             << std::endl;
    std::cout << " Turn on ice growth?     : " << Input_Opt.AEROSOL_ICE_GROWTH                       << std::endl;

} /* End of Read_Aerosol_Menu */

void Read_Meteorology_Menu( OptInput &Input_Opt, bool &RC )
{

    /* DESCRIPTION: Function Read\_Meteorology\_Menu reads the METEOROLOGY MENU
     * section of the APCEMM input file. */

    /* INPUT/OUTPUT PARAMETERS:
     * - Input_Opt: Input options
     * - RC: Success or failure
     */

    /* Read until all lines from the menu are read */

    std::vector<std::string> tokens;
    std::string variable;
    double value;

    /* ==================================================== */
    /* Do we have MET input?                                */
    /* ==================================================== */

    variable = "Do we have MET input?";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    if ( ( strcmp(tokens[0].c_str(), "T" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "t" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "TRUE" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "true" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "True" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "YES" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Y" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "y" )    == 0 ) )
        Input_Opt.MET_LOADMET = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "f" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "FALSE" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "false" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "False" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "NO" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "No" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "no" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "N" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "n" )     == 0 ) )
        Input_Opt.MET_LOADMET = 0;
    else {
        std::cout << " Wrong input for: " << variable << std::endl;
        exit(1);
    }

    /* ==================================================== */
    /* Met file                                             */
    /* ==================================================== */

    variable = "Met file";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    tokens[0].erase(std::remove(tokens[0].begin(), tokens[0].end(), '*'), tokens[0].end());
    Input_Opt.MET_FILENAME = tokens[0];

    /* ==================================================== */
    /* Initialize T from MET?                               */
    /* ==================================================== */

    variable = "Initialize T from MET?";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    if ( ( strcmp(tokens[0].c_str(), "T" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "t" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "TRUE" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "true" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "True" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "YES" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Y" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "y" )    == 0 ) )
        Input_Opt.MET_LOADTEMP = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "f" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "FALSE" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "false" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "False" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "NO" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "No" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "no" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "N" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "n" )     == 0 ) )
        Input_Opt.MET_LOADTEMP = 0;
    else {
        std::cout << " Wrong input for: " << variable << std::endl;
        exit(1);
    }

    /* ==================================================== */
    /* Interpolate temperature in time?                     */
    /* ==================================================== */

    variable = "Interpolate T in time?";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    if ( ( strcmp(tokens[0].c_str(), "T" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "t" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "TRUE" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "true" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "True" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "YES" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Y" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "y" )    == 0 ) )
        Input_Opt.MET_INTERPTEMP = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "f" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "FALSE" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "false" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "False" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "NO" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "No" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "no" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "N" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "n" )     == 0 ) )
        Input_Opt.MET_INTERPTEMP = 0;
    else {
        std::cout << " Wrong input for: " << variable << std::endl;
        exit(1);
    }

    /* ==================================================== */
    /* Initialize H2O from MET?                             */
    /* ==================================================== */

    variable = "Initialize H2O from MET?";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    if ( ( strcmp(tokens[0].c_str(), "T" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "t" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "TRUE" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "true" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "True" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "YES" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Y" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "y" )    == 0 ) )
        Input_Opt.MET_LOADH2O = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "f" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "FALSE" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "false" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "False" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "NO" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "No" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "no" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "N" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "n" )     == 0 ) )
        Input_Opt.MET_LOADH2O = 0;
    else {
        std::cout << " Wrong input for: " << variable << std::endl;
        exit(1);
    }

    /* ==================================================== */
    /* Interpolate H2O in time?                             */
    /* ==================================================== */

    variable = "Interpolate H2O in time?";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    if ( ( strcmp(tokens[0].c_str(), "T" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "t" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "TRUE" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "true" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "True" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "YES" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Y" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "y" )    == 0 ) )
        Input_Opt.MET_INTERPH2O = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "f" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "FALSE" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "false" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "False" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "NO" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "No" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "no" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "N" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "n" )     == 0 ) )
        Input_Opt.MET_INTERPH2O = 0;
    else {
        std::cout << " Wrong input for: " << variable << std::endl;
        exit(1);
    }

    /* ==================================================== */
    /* Initialize shear from MET?                           */
    /* ==================================================== */

    variable = "Initialize shear from MET?";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    if ( ( strcmp(tokens[0].c_str(), "T" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "t" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "TRUE" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "true" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "True" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "YES" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Y" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "y" )    == 0 ) )
        Input_Opt.MET_LOADSHEAR = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "f" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "FALSE" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "false" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "False" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "NO" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "No" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "no" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "N" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "n" )     == 0 ) )
        Input_Opt.MET_LOADSHEAR = 0;
    else {
        std::cout << " Wrong input for: " << variable << std::endl;
        exit(1);
    }

    /* ==================================================== */
    /* Interpolate shear in time?                           */
    /* ==================================================== */

    variable = "Interpolate shear in time?";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    if ( ( strcmp(tokens[0].c_str(), "T" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "t" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "TRUE" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "true" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "True" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "YES" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Y" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "y" )    == 0 ) )
        Input_Opt.MET_INTERPSHEAR = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "f" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "FALSE" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "false" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "False" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "NO" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "No" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "no" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "N" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "n" )     == 0 ) )
        Input_Opt.MET_INTERPSHEAR = 0;
    else {
        std::cout << " Wrong input for: " << variable << std::endl;
        exit(1);
    }

    /* Skip line */
    getline( inputFile, line, '\n' );

    /* ==================================================== */
    /* Impose moist layer depth                             */
    /* ==================================================== */

    variable = "Impose moist layer depth?";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    if ( ( strcmp(tokens[0].c_str(), "T" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "t" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "TRUE" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "true" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "True" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "YES" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Y" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "y" )    == 0 ) )
        Input_Opt.MET_FIXDEPTH = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "f" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "FALSE" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "false" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "False" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "NO" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "No" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "no" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "N" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "n" )     == 0 ) )
        Input_Opt.MET_FIXDEPTH = 0;
    else {
        std::cout << " Wrong input for: " << variable << std::endl;
        exit(1);
    }

    /* ==================================================== */
    /* Moist layer depth                                    */
    /* ==================================================== */

    variable = "Moist layer depth [m]";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    try {
        value = std::stod( tokens[0] );
        if ( value > 0.0E+00 )
            Input_Opt.MET_DEPTH = value;
        else {
            std::cout << " Wrong input for: " << variable << std::endl;
            std::cout << " Timestep needs to be positive" << std::endl;
            exit(1);
        }
    } catch(std::exception& e) {
        std::cout << " Could not convert string to double for " << variable << std::endl;
        exit(1);
    }

    /* Skip line */
    getline( inputFile, line, '\n' );

    /* ==================================================== */
    /* Impose temperature lapse rate                        */
    /* ==================================================== */

    variable = "Impose temperature lapse rate?";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    if ( ( strcmp(tokens[0].c_str(), "T" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "t" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "TRUE" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "true" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "True" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "YES" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Y" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "y" )    == 0 ) )
        Input_Opt.MET_FIXLAPSERATE = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "f" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "FALSE" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "false" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "False" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "NO" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "No" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "no" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "N" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "n" )     == 0 ) )
        Input_Opt.MET_FIXLAPSERATE = 0;
    else {
        std::cout << " Wrong input for: " << variable << std::endl;
        exit(1);
    }

    /* ==================================================== */
    /* Temperature lapse rate                               */
    /* ==================================================== */

    variable = "Temperature lapse rate [K/m]";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    try {
        value = std::stod( tokens[0] );
        Input_Opt.MET_LAPSERATE = value;
    } catch(std::exception& e) {
        std::cout << " Could not convert string to double for " << variable << std::endl;
        exit(1);
    }

    /* Skip line */
    getline( inputFile, line, '\n' );

    /* ==================================================== */
    /* Add diurnal variations?                              */
    /* ==================================================== */

    variable = "Diurnal variations?";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    if ( ( strcmp(tokens[0].c_str(), "T" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "t" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "TRUE" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "true" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "True" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "YES" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Y" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "y" )    == 0 ) )
        Input_Opt.MET_DIURNAL = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "f" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "FALSE" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "false" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "False" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "NO" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "No" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "no" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "N" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "n" )     == 0 ) )
        Input_Opt.MET_DIURNAL = 0;
    else {
        std::cout << " Wrong input for: " << variable << std::endl;
        exit(1);
    }


    /* Return success */
    RC = SUCCESS;

    /* ==================================================== */
    /* Print to screen                                      */
    /* ==================================================== */

    std::cout << " ------------------------+------------------------------------------------------ " << std::endl;
    std::cout << " %%% METEOROLOGY MENU %%%:"                                                        << std::endl;
    std::cout << " ------------------------+------------------------------------------------------ " << std::endl;
    std::cout << " Do we have MET input?   : " << Input_Opt.MET_LOADMET                              << std::endl;
    std::cout << "  => Met file            : " << Input_Opt.MET_FILENAME                             << std::endl;
    std::cout << "  => Init T from MET?    : " << Input_Opt.MET_LOADTEMP                             << std::endl;
    std::cout << "    => Interpolate T?    : " << Input_Opt.MET_INTERPTEMP                           << std::endl;
    std::cout << "  => Init H2O from MET?  : " << Input_Opt.MET_LOADH2O                              << std::endl;
    std::cout << "    => Interpolate H2O?  : " << Input_Opt.MET_INTERPH2O                            << std::endl;
    std::cout << "  => Init shear from MET?: " << Input_Opt.MET_LOADSHEAR                            << std::endl;
    std::cout << "    => Interpolate shear?: " << Input_Opt.MET_INTERPSHEAR                          << std::endl;
    std::cout << " - OR -------------------: " << std::endl;
    std::cout << " Impose moist layer depth: " << Input_Opt.MET_FIXDEPTH                             << std::endl;
    std::cout << "  => Moist layer depth[m]: " << Input_Opt.MET_DEPTH                                << std::endl;
    std::cout << " ---- OR ----------------: " << std::endl;
    std::cout << " Impose temp. lapse rate : " << Input_Opt.MET_FIXLAPSERATE                         << std::endl;
    std::cout << "  => Lapse rate [K/m]    : " << Input_Opt.MET_LAPSERATE                            << std::endl;
    std::cout << " ------------------------: " << std::endl;
    std::cout << " Add diurnal variations  : " << Input_Opt.MET_DIURNAL                              << std::endl;

} /* End of Read_Meteorology_Menu */

void Read_Diagnostic_Menu( OptInput &Input_Opt, bool &RC )
{

    /* DESCRIPTION: Function Read\_Diagnostic\_Menu reads the DIAGNOSTIC MENU
     * section of the APCEMM input file. */

    /* INPUT/OUTPUT PARAMETERS:
     * - Input_Opt: Input options
     * - RC: Success or failure
     */

    /* Read until all lines from the menu are read */

    std::vector<std::string> tokens;
    std::string variable;

    /* ==================================================== */
    /* netCDF file name                                     */
    /* ==================================================== */

    variable = "netCDF file name";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    tokens[0].erase(std::remove(tokens[0].begin(), tokens[0].end(), '*'), tokens[0].end());
    Input_Opt.DIAG_FILENAME = tokens[0];

    /* Skip line */
    getline( inputFile, line, '\n' );

    /* Return success */
    RC = SUCCESS;

    /* ==================================================== */
    /* Print to screen                                      */
    /* ==================================================== */

    std::cout << " ------------------------+------------------------------------------------------ " << std::endl;
    std::cout << " %%% DIAGNOSTIC MENU %%% :"                                                        << std::endl;
    std::cout << " ------------------------+------------------------------------------------------ " << std::endl;
    std::cout << " netCDF file name        : " << Input_Opt.DIAG_FILENAME                            << std::endl;
    std::cout << " Diagnostic Entries ---> : L"                                                      << std::endl;

} /* End of Read_Diagnostic_Menu */

void Read_Timeseries_Menu( OptInput &Input_Opt, bool &RC )
{

    /* DESCRIPTION: Function Read\_Timeseries\_Menu reads the TIMESERIES MENU
     * section of the APCEMM input file. */

    /* INPUT/OUTPUT PARAMETERS:
     * - Input_Opt: Input options
     * - RC: Success or failure
     */

    /* Read until all lines from the menu are read */

    std::vector<std::string> tokens;
    std::string variable;
    int value;

    /* ==================================================== */
    /* Save species timeseries?                             */
    /* ==================================================== */

    variable = "Save species timeseries?";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    if ( ( strcmp(tokens[0].c_str(), "T" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "t" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "TRUE" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "true" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "True" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "YES" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Y" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "y" )    == 0 ) )
        Input_Opt.TS_SPEC = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "f" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "FALSE" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "false" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "False" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "NO" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "No" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "no" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "N" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "n" )     == 0 ) )
        Input_Opt.TS_SPEC = 0;
    else {
        std::cout << " Wrong input for: " << variable << std::endl;
        exit(1);
    }

    /* ==================================================== */
    /* Inst timeseries file                                 */
    /* ==================================================== */

    variable = "Inst timeseries file";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    tokens[0].erase(std::remove(tokens[0].begin(), tokens[0].end(), '*'), tokens[0].end());
    Input_Opt.TS_FILENAME = tokens[0];

    /* ==================================================== */
    /* Species to include                                   */
    /* ==================================================== */

    variable = "Species to include";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    /* Store in values for variable */
    for ( unsigned int i = 0; i < tokens.size(); i++ ) {
        if ( strcmp(tokens[i].c_str(), EMPTY) != 0 ) {
            try {
                value = std::stoi( tokens[i] );
                if ( value > 0 )
                    Input_Opt.TS_SPECIES.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << " Index needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[i] << "' to int for " << variable << std::endl;
                exit(1);
            }
        }
    }

    /* ==================================================== */
    /* Frequency [min]                                      */
    /* ==================================================== */

    variable = "Frequency [min]";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    if ( ( strcmp(tokens[0].c_str(), "-") == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Inf") == 0 ) )
        Input_Opt.TS_FREQ = 0.0E+00;
    else {
        try {
            value = std::stod( tokens[0] );
            if ( value >= 0.0E+00 )
                Input_Opt.TS_FREQ = value;
            else {
                std::cout << " Wrong input for: " << variable << std::endl;
                std::cout << " Frequency needs to be positive" << std::endl;
                exit(1);
            }
        } catch(std::exception& e) {
            std::cout << " Could not convert string to double for " << variable << std::endl;
            exit(1);
        }
    }

    /* ==================================================== */
    /* Save aerosol timeseries?                             */
    /* ==================================================== */

    variable = "Save aerosol timeseries?";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    if ( ( strcmp(tokens[0].c_str(), "T" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "t" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "TRUE" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "true" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "True" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "YES" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Y" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "y" )    == 0 ) )
        Input_Opt.TS_AERO = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "f" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "FALSE" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "false" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "False" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "NO" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "No" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "no" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "N" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "n" )     == 0 ) )
        Input_Opt.TS_AERO = 0;
    else {
        std::cout << " Wrong input for: " << variable << std::endl;
        exit(1);
    }

    /* ==================================================== */
    /* Inst timeseries file                                 */
    /* ==================================================== */

    variable = "Inst timeseries file";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    tokens[0].erase(std::remove(tokens[0].begin(), tokens[0].end(), '*'), tokens[0].end());
    Input_Opt.TS_AERO_FILENAME = tokens[0];

    /* ==================================================== */
    /* Aerosol to include                                   */
    /* ==================================================== */

    variable = "Aerosol to include";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    /* Store in values for variable */
    for ( unsigned int i = 0; i < tokens.size(); i++ ) {
        if ( strcmp(tokens[i].c_str(), EMPTY) != 0 ) {
            try {
                value = std::stoi( tokens[i] );
                if ( value > 0 )
                    Input_Opt.TS_AEROSOL.push_back( value );
                else {
                    std::cout << " Wrong input for: " << variable << std::endl;
                    std::cout << " Index needs to be positive" << std::endl;
                    exit(1);
                }
            } catch(std::exception& e) {
                std::cout << " Could not convert string '" << tokens[i] << "' to int for " << variable << std::endl;
                exit(1);
            }
        }
    }

    /* ==================================================== */
    /* Frequency [min]                                      */
    /* ==================================================== */

    variable = "Frequency [min]";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    if ( ( strcmp(tokens[0].c_str(), "-") == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Inf") == 0 ) )
        Input_Opt.TS_AERO_FREQ = 0.0E+00;
    else {
        try {
            value = std::stod( tokens[0] );
            if ( value >= 0.0E+00 )
                Input_Opt.TS_AERO_FREQ = value;
            else {
                std::cout << " Wrong input for: " << variable << std::endl;
                std::cout << " Frequency needs to be positive" << std::endl;
                exit(1);
            }
        } catch(std::exception& e) {
            std::cout << " Could not convert string to double for " << variable << std::endl;
            exit(1);
        }
    }

    /* Return success */
    RC = SUCCESS;

    /* ==================================================== */
    /* Print to screen                                      */
    /* ==================================================== */

    std::cout << " ------------------------+------------------------------------------------------ " << std::endl;
    std::cout << " %%% TIMESERIES MENU %%% :"                                                        << std::endl;
    std::cout << " ------------------------+------------------------------------------------------ " << std::endl;
    std::cout << " Save species timeseries?: " << Input_Opt.TS_SPEC                                  << std::endl;
    std::cout << "  => Inst timeseries file: " << Input_Opt.TS_FILENAME                              << std::endl;
    std::cout << "  => Species to include  : ";
    for ( unsigned int i = 0; i < Input_Opt.TS_SPECIES.size(); i++ )
        std::cout << Input_Opt.TS_SPECIES[i] << " ";
    std::cout                                                                << std::endl;
    std::cout << "  => Frequency [min]     : " << Input_Opt.TS_FREQ          << std::endl;
    std::cout << " Save aerosol timeseries?: " << Input_Opt.TS_AERO          << std::endl;
    std::cout << "  => Inst timeseries file: " << Input_Opt.TS_AERO_FILENAME << std::endl;
    std::cout << "  => Aerosol to include  : ";
    for ( unsigned int i = 0; i < Input_Opt.TS_AEROSOL.size(); i++ )
        std::cout << Input_Opt.TS_AEROSOL[i] << " ";
    std::cout                                                            << std::endl;
    std::cout << "  => Frequency [min]     : " << Input_Opt.TS_AERO_FREQ << std::endl;

} /* End of Read_Timeseries_Menu */

void Read_PL_Menu( OptInput &Input_Opt, bool &RC )
{

    /* DESCRIPTION: Function Read\_PL\_Menu reads the PROD & LOSS MENU
     * section of the APCEMM input file. */

    /* INPUT/OUTPUT PARAMETERS:
     * - Input_Opt: Input options
     * - RC: Success or failure
     */

    /* Read until all lines from the menu are read */

    std::vector<std::string> tokens;
    std::string variable;

    /* ==================================================== */
    /* Turn on P/L diag?                                    */
    /* ==================================================== */

    variable = "Turn on P/L diag?";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable range */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    if ( ( strcmp(tokens[0].c_str(), "T" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "t" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "TRUE" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "true" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "True" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "YES" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Y" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "y" )    == 0 ) )
        Input_Opt.PL_PL = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "f" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "FALSE" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "false" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "False" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "NO" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "No" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "no" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "N" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "n" )     == 0 ) )
        Input_Opt.PL_PL = 0;
    else {
        std::cout << " Wrong input for: " << variable << std::endl;
        exit(1);
    }

    /* ==================================================== */
    /* Save O3 P/L?                                         */
    /* ==================================================== */

    variable = "Save O3 P/L?";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;

    /* Extract variable range */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    if ( ( strcmp(tokens[0].c_str(), "T" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "t" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "TRUE" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "true" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "True" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "YES" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "yes" )  == 0 ) || \
         ( strcmp(tokens[0].c_str(), "Y" )    == 0 ) || \
         ( strcmp(tokens[0].c_str(), "y" )    == 0 ) )
        Input_Opt.PL_O3 = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "f" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "FALSE" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "false" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "False" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "NO" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "No" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "no" )    == 0 ) || \
              ( strcmp(tokens[0].c_str(), "N" )     == 0 ) || \
              ( strcmp(tokens[0].c_str(), "n" )     == 0 ) )
        Input_Opt.PL_O3 = 0;
    else {
        std::cout << " Wrong input for: " << variable << std::endl;
        exit(1);
    }

    if ( ( Input_Opt.CHEMISTRY_CHEMISTRY == 0 ) && \
         ( ( Input_Opt.PL_PL == 1 ) || ( Input_Opt.PL_O3 == 1 ) ) ) {
         std::cout << " Chemistry if turned off but rate output is on!" << std::endl;
         std::cout << " Turning off rate output!" << std::endl;
         Input_Opt.PL_PL = 0;
         Input_Opt.PL_O3 = 0;
    }

    /* Return success */
    RC = SUCCESS;

    /* ==================================================== */
    /* Print to screen                                      */
    /* ==================================================== */

    std::cout << " ------------------------+------------------------------------------------------ " << std::endl;
    std::cout << " %%% PROD & LOSS MENU %%%:"                                                        << std::endl;
    std::cout << " ------------------------+------------------------------------------------------ " << std::endl;
    std::cout << " Turn on P/L diag?       : " << Input_Opt.PL_PL                                    << std::endl;
    std::cout << " Save O3 P/L?            : " << Input_Opt.PL_O3                                    << std::endl;
    std::cout << " ------------------------+------------------------------------------------------ " << std::endl;

} /* End of Read_PL_Menu */

Vector_2D CombVec( OptInput &Input_Opt )
{

    const bool print = 0;

    unsigned int counter = 1;
    unsigned int nCases;

    unsigned int i, j;

    double currVal = 0.0E+00;

    Vector_1D cases;
    Vector_2D y, z;
    Vector_2D u, v;

    if ( Input_Opt.SIMULATION_MONTECARLO ) {

        /* Initialize seed */
        setSeed();

        /* ======================================================================= */
        /* ---- PLUME PROCESSING TIME ( SIMULATION TIME ) ------------------------ */
        /* ---- Accepted units are: hr (default)                                   */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_PLUMEPROCESS_RANGE ) {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( fRand(Input_Opt.PARAMETER_PLUMEPROCESS[0], \
                                       Input_Opt.PARAMETER_PLUMEPROCESS[1]) );
        } else {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( Input_Opt.PARAMETER_PLUMEPROCESS[0] );
        }

        /* Do unit conversion to default unit */
        if ( ( Input_Opt.PARAMETER_PLUMEPROCESS_UNIT.compare( "0-24" ) == 0 ) || \
             ( Input_Opt.PARAMETER_PLUMEPROCESS_UNIT.compare( "hr" ) == 0 ) ) {
            /* Do nothing. Default unit */
        } else if ( Input_Opt.PARAMETER_PLUMEPROCESS_UNIT.compare( "min" ) == 0 ) {
            /* Convert from min to hr */
            for ( i = 0; i < cases.size(); i++ )
                cases[i] /= 60.0;
        } else {
            std::cout << " Unknown unit for variable 'Plume Processing Time': ";
            std::cout << Input_Opt.PARAMETER_PLUMEPROCESS_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_PLUMEPROCESS_UNIT = "hr";

        y.push_back(Vector_1D(Input_Opt.SIMULATION_MCRUNS));
        for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
            y[counter-1][i] = cases[i];
        cases.clear();

        /* ======================================================================= */
        /* ---- TEMPERATURE ------------------------------------------------------ */
        /* ---- Accepted units are: Kelvin (default), Celsius, Fahrenheit, Rankine */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_TEMPERATURE_RANGE ) {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( fRand(Input_Opt.PARAMETER_TEMPERATURE[0], \
                                       Input_Opt.PARAMETER_TEMPERATURE[1]) );
        } else {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( Input_Opt.PARAMETER_TEMPERATURE[0] );
        }

        /* Do unit conversion to default unit */
        if ( Input_Opt.PARAMETER_TEMPERATURE_UNIT.compare( "K" ) == 0 ) {
            /* Do nothing. Default unit */
        } else if ( Input_Opt.PARAMETER_TEMPERATURE_UNIT.compare( "C" ) == 0 ) {
            /* Convert C to K */
            for ( i = 0; i < cases.size(); i++ )
                cases[i] += 273.15;
        } else if ( Input_Opt.PARAMETER_TEMPERATURE_UNIT.compare( "F" ) == 0 ) {
            /* Convert F to K */
            for ( i = 0; i < cases.size(); i++ )
                cases[i] = ( cases[i] + 459.67 ) * (double) 5.0/9;
        } else if ( Input_Opt.PARAMETER_TEMPERATURE_UNIT.compare( "R" ) == 0 ) {
            /* Convert R to K */
            for ( i = 0; i < cases.size(); i++ )
                cases[i] *= (double) 5.0/9;
        } else {
            std::cout << " Unknown unit for variable 'Temperature': ";
            std::cout << Input_Opt.PARAMETER_TEMPERATURE_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_TEMPERATURE_UNIT = "K";

        counter += 1;
        y.push_back(Vector_1D(Input_Opt.SIMULATION_MCRUNS));
        for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
            y[counter-1][i] = cases[i];
        cases.clear();

        /* ======================================================================= */
        /* ---- RELATIVE HUMIDITY ------------------------------------------------ */
        /* ---- Accepted units are: % (0-100, default), - (0-1)                    */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_RHW_RANGE ) {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( fRand(Input_Opt.PARAMETER_RHW[0], \
                                       Input_Opt.PARAMETER_RHW[1]) );
        } else {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( Input_Opt.PARAMETER_RHW[0] );
        }

        if ( Input_Opt.PARAMETER_RHW_UNIT.compare( "%" ) == 0 ) {
            /* Do nothing. Default unit */
        } else if ( Input_Opt.PARAMETER_RHW_UNIT.compare( "-" ) == 0 ) {
            /* Convert - to % */
            for ( i = 0; i < cases.size(); i++ )
                cases[i] *= 100.0;
        } else {
            std::cout << " Unknown unit for variable 'Relative humidity': ";
            std::cout << Input_Opt.PARAMETER_RHW_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_RHW_UNIT = "%";

        counter += 1;
        y.push_back(Vector_1D(Input_Opt.SIMULATION_MCRUNS));
        for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
            y[counter-1][i] = cases[i];
        cases.clear();

        /* ======================================================================= */
        /* ---- HORIZONTAL DIFFUSION PARAMETER ----------------------------------- */
        /* ---- Accepted units are: m^2/s (default)                                */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_DH_RANGE ) {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( fRand(Input_Opt.PARAMETER_DH[0], \
                                       Input_Opt.PARAMETER_DH[1]) );
        } else {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( Input_Opt.PARAMETER_DH[0] );
        }

        if ( Input_Opt.PARAMETER_DH_UNIT.compare( "m^2/s" ) == 0 ) {
            /* Do nothing. Default unit */
        } else {
            std::cout << " Unknown unit for variable 'Horiz. diffusion parameter': ";
            std::cout << Input_Opt.PARAMETER_DH_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_DH_UNIT = "m^2/s";

        counter += 1;
        y.push_back(Vector_1D(Input_Opt.SIMULATION_MCRUNS));
        for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
            y[counter-1][i] = cases[i];
        cases.clear();

        /* ======================================================================= */
        /* ---- VERTICAL DIFFUSION PARAMETER ------------------------------------- */
        /* ---- Accepted units are: m^2/s (default)                                */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_DV_RANGE ) {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( fRand(Input_Opt.PARAMETER_DV[0], \
                                       Input_Opt.PARAMETER_DV[1]) );
        } else {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( Input_Opt.PARAMETER_DV[0] );
        }

        if ( Input_Opt.PARAMETER_DV_UNIT.compare( "m^2/s" ) == 0 ) {
            /* Do nothing. Default unit */
        } else {
            std::cout << " Unknown unit for variable 'Vert. diffusion parameter': ";
            std::cout << Input_Opt.PARAMETER_DV_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_DV_UNIT = "m^2/s";

        counter += 1;
        y.push_back(Vector_1D(Input_Opt.SIMULATION_MCRUNS));
        for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
            y[counter-1][i] = cases[i];
        cases.clear();

        /* ======================================================================= */
        /* ---- SHEAR ------------------------------------------------------------ */
        /* ---- Accepted units are: 1/s (default)                                  */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_SHEAR_RANGE ) {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( fRand(Input_Opt.PARAMETER_SHEAR[0], \
                                       Input_Opt.PARAMETER_SHEAR[1]) );
        } else {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( Input_Opt.PARAMETER_SHEAR[0] );
        }

        if ( Input_Opt.PARAMETER_SHEAR_UNIT.compare( "1/s" ) == 0 ) {
            /* Do nothing. Default unit */
        } else {
            std::cout << " Unknown unit for variable 'Shear': ";
            std::cout << Input_Opt.PARAMETER_SHEAR_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_SHEAR_UNIT = "1/s";

        counter += 1;
        y.push_back(Vector_1D(Input_Opt.SIMULATION_MCRUNS));
        for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
            y[counter-1][i] = cases[i];
        cases.clear();

        /* ======================================================================= */
        /* ---- LONGITUDE -------------------------------------------------------- */
        /* ---- Accepted units are: degree (default)                               */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_LONGITUDE_RANGE ) {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( fRand(Input_Opt.PARAMETER_LONGITUDE[0], \
                                       Input_Opt.PARAMETER_LONGITUDE[1]) );
        } else {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( Input_Opt.PARAMETER_LONGITUDE[0] );
        }

        if ( Input_Opt.PARAMETER_LONGITUDE_UNIT.compare( "deg" ) == 0 ) {
            /* Do nothing. Default unit */
        } else {
            std::cout << " Unknown unit for variable 'Longitude': ";
            std::cout << Input_Opt.PARAMETER_LONGITUDE_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_LONGITUDE_UNIT = "deg";

        counter += 1;
        y.push_back(Vector_1D(Input_Opt.SIMULATION_MCRUNS));
        for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
            y[counter-1][i] = cases[i];
        cases.clear();

        /* ======================================================================= */
        /* ---- LATITUDE --------------------------------------------------------- */
        /* ---- Accepted units are: degree (default)                               */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_LATITUDE_RANGE ) {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( fRand(Input_Opt.PARAMETER_LATITUDE[0], \
                                       Input_Opt.PARAMETER_LATITUDE[1]) );
        } else {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( Input_Opt.PARAMETER_LATITUDE[0] );
        }

        if ( Input_Opt.PARAMETER_LATITUDE_UNIT.compare( "deg" ) == 0 ) {
            /* Do nothing. Default unit */
        } else {
            std::cout << " Unknown unit for variable 'Latitude': ";
            std::cout << Input_Opt.PARAMETER_LATITUDE_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_LATITUDE_UNIT = "deg";

        counter += 1;
        y.push_back(Vector_1D(Input_Opt.SIMULATION_MCRUNS));
        for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
            y[counter-1][i] = cases[i];
        cases.clear();

        /* ======================================================================= */
        /* ---- PRESSURE --------------------------------------------------------- */
        /* ---- Accepted units are: Pa (default), hPa                              */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_PRESSURE_RANGE ) {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( fRand(Input_Opt.PARAMETER_PRESSURE[0], \
                                       Input_Opt.PARAMETER_PRESSURE[1]) );
        } else {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( Input_Opt.PARAMETER_PRESSURE[0] );
        }

        if ( Input_Opt.PARAMETER_PRESSURE_UNIT.compare( "Pa" ) == 0 ) {
            /* Do nothing. Default unit */
        } else if ( Input_Opt.PARAMETER_PRESSURE_UNIT.compare( "hPa" ) == 0 ) {
            for ( i = 0; i < cases.size(); i++ )
                cases[i] *= 100.0;
        } else {
            std::cout << " Unknown unit for variable 'Pressure': ";
            std::cout << Input_Opt.PARAMETER_PRESSURE_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_PRESSURE_UNIT = "Pa";

        counter += 1;
        y.push_back(Vector_1D(Input_Opt.SIMULATION_MCRUNS));
        for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
            y[counter-1][i] = cases[i];
        cases.clear();

        /* ======================================================================= */
        /* ---- EMISSION DAY ----------------------------------------------------- */
        /* ---- Accepted units are: 1-365 (default)                                */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_EDAY_RANGE ) {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( fRand((int) Input_Opt.PARAMETER_EDAY[0], \
                                       (int) Input_Opt.PARAMETER_EDAY[1]) );
        } else {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( Input_Opt.PARAMETER_EDAY[0] );
        }

        if ( ( Input_Opt.PARAMETER_EDAY_UNIT.compare( "1-365" ) == 0 ) || \
             ( Input_Opt.PARAMETER_EDAY_UNIT.compare( "-" ) == 0 ) ) {
            /* Do nothing. Default unit */
        } else {
            std::cout << " Unknown unit for variable 'Emission Day': ";
            std::cout << Input_Opt.PARAMETER_EDAY_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_EDAY_UNIT = "1-365";

        counter += 1;
        y.push_back(Vector_1D(Input_Opt.SIMULATION_MCRUNS));
        for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
            y[counter-1][i] = cases[i];
        cases.clear();

        /* ======================================================================= */
        /* ---- EMISSION TIME ---------------------------------------------------- */
        /* ---- Accepted units are: 0-24 (default)                                 */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_ETIME_RANGE ) {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( fRand(Input_Opt.PARAMETER_ETIME[0], \
                                       Input_Opt.PARAMETER_ETIME[1]) );
        } else {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( Input_Opt.PARAMETER_ETIME[0] );
        }

        if ( ( Input_Opt.PARAMETER_ETIME_UNIT.compare( "0-24" ) == 0 ) || \
             ( Input_Opt.PARAMETER_ETIME_UNIT.compare( "hr" ) == 0 ) ) {
            /* Do nothing. Default unit */
        } else {
            std::cout << " Unknown unit for variable 'Emission Time': ";
            std::cout << Input_Opt.PARAMETER_ETIME_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_ETIME_UNIT = "0-24";

        counter += 1;
        y.push_back(Vector_1D(Input_Opt.SIMULATION_MCRUNS));
        for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
            y[counter-1][i] = cases[i];
        cases.clear();

        /* ======================================================================= */
        /* ---- EI_NOX ----------------------------------------------------------- */
        /* ---- Accepted units are: g/kg_fuel (default)                            */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_EI_NOX_RANGE ) {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( fRand(Input_Opt.PARAMETER_EI_NOX[0], \
                                       Input_Opt.PARAMETER_EI_NOX[1]) );
        } else {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( Input_Opt.PARAMETER_EI_NOX[0] );
        }

        if ( ( Input_Opt.PARAMETER_EI_NOX_UNIT.compare( "g/kg_fuel" )      == 0 ) || \
             ( Input_Opt.PARAMETER_EI_NOX_UNIT.compare( "g/kg_f" )         == 0 ) || \
             ( Input_Opt.PARAMETER_EI_NOX_UNIT.compare( "g/kg" )           == 0 ) || \
             ( Input_Opt.PARAMETER_EI_NOX_UNIT.compare( "gNO2/kg_fuel" )   == 0 ) || \
             ( Input_Opt.PARAMETER_EI_NOX_UNIT.compare( "gNO2/kg_f" )      == 0 ) || \
             ( Input_Opt.PARAMETER_EI_NOX_UNIT.compare( "gNO2/kg" )        == 0 ) || \
             ( Input_Opt.PARAMETER_EI_NOX_UNIT.compare( "g(NO2)/kg_fuel" ) == 0 ) || \
             ( Input_Opt.PARAMETER_EI_NOX_UNIT.compare( "g(NO2)/kg_f" )    == 0 ) || \
             ( Input_Opt.PARAMETER_EI_NOX_UNIT.compare( "g(NO2)/kg" )      == 0 ) ) {
            /* Do nothing. Default unit */
        } else {
            std::cout << " Unknown unit for variable 'EI_NOX': ";
            std::cout << Input_Opt.PARAMETER_EI_NOX_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_EI_NOX_UNIT = "g/kg_fuel";

        counter += 1;
        y.push_back(Vector_1D(Input_Opt.SIMULATION_MCRUNS));
        for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
            y[counter-1][i] = cases[i];
        cases.clear();

        /* ======================================================================= */
        /* ---- EI_CO ------------------------------------------------------------ */
        /* ---- Accepted units are: g/kg_fuel (default)                            */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_EI_CO_RANGE ) {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( fRand(Input_Opt.PARAMETER_EI_CO[0], \
                                       Input_Opt.PARAMETER_EI_CO[1]) );
        } else {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( Input_Opt.PARAMETER_EI_CO[0] );
        }

        if ( ( Input_Opt.PARAMETER_EI_CO_UNIT.compare( "g/kg_fuel" ) == 0 ) || \
             ( Input_Opt.PARAMETER_EI_CO_UNIT.compare( "g/kg_f" )    == 0 ) || \
             ( Input_Opt.PARAMETER_EI_CO_UNIT.compare( "g/kg" )      == 0 ) ) {
            /* Do nothing. Default unit */
        } else {
            std::cout << " Unknown unit for variable 'EI_CO': ";
            std::cout << Input_Opt.PARAMETER_EI_CO_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_EI_CO_UNIT = "g/kg_fuel";

        counter += 1;
        y.push_back(Vector_1D(Input_Opt.SIMULATION_MCRUNS));
        for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
            y[counter-1][i] = cases[i];
        cases.clear();

        /* ======================================================================= */
        /* ---- EI_UHC ----------------------------------------------------------- */
        /* ---- Accepted units are: g/kg_fuel (default)                            */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_EI_UHC_RANGE ) {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( fRand(Input_Opt.PARAMETER_EI_UHC[0], \
                                       Input_Opt.PARAMETER_EI_UHC[1]) );
        } else {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( Input_Opt.PARAMETER_EI_UHC[0] );
        }

        if ( ( Input_Opt.PARAMETER_EI_UHC_UNIT.compare( "g/kg_fuel" ) == 0 ) || \
             ( Input_Opt.PARAMETER_EI_UHC_UNIT.compare( "g/kg_f" )    == 0 ) || \
             ( Input_Opt.PARAMETER_EI_UHC_UNIT.compare( "g/kg" )      == 0 ) ) {
            /* Do nothing. Default unit */
        } else {
            std::cout << " Unknown unit for variable 'EI_UHC': ";
            std::cout << Input_Opt.PARAMETER_EI_UHC_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_EI_UHC_UNIT = "g/kg_fuel";

        counter += 1;
        y.push_back(Vector_1D(Input_Opt.SIMULATION_MCRUNS));
        for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
            y[counter-1][i] = cases[i];
        cases.clear();

        /* ======================================================================= */
        /* ---- EI_SO2 ----------------------------------------------------------- */
        /* ---- Accepted units are: g/kg_fuel (default)                            */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_EI_SO2_RANGE ) {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( fRand(Input_Opt.PARAMETER_EI_SO2[0], \
                                       Input_Opt.PARAMETER_EI_SO2[1]) );
        } else {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( Input_Opt.PARAMETER_EI_SO2[0] );
        }

        if ( ( Input_Opt.PARAMETER_EI_SO2_UNIT.compare( "g/kg_fuel" ) == 0 ) || \
             ( Input_Opt.PARAMETER_EI_SO2_UNIT.compare( "g/kg_f" )    == 0 ) || \
             ( Input_Opt.PARAMETER_EI_SO2_UNIT.compare( "g/kg" )      == 0 ) ) {
            /* Do nothing. Default unit */
        } else {
            std::cout << " Unknown unit for variable 'EI_SO2': ";
            std::cout << Input_Opt.PARAMETER_EI_SO2_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_EI_SO2_UNIT = "g/kg_fuel";

        counter += 1;
        y.push_back(Vector_1D(Input_Opt.SIMULATION_MCRUNS));
        for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
            y[counter-1][i] = cases[i];
        cases.clear();

        /* ======================================================================= */
        /* ---- EI_SO2TOSO4 ------------------------------------------------------ */
        /* ---- Accepted units are: - (default), %                                 */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_EI_SO2TOSO4_RANGE ) {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( fRand(Input_Opt.PARAMETER_EI_SO2TOSO4[0], \
                                       Input_Opt.PARAMETER_EI_SO2TOSO4[1]) );
        } else {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( Input_Opt.PARAMETER_EI_SO2TOSO4[0] );
        }

        if ( Input_Opt.PARAMETER_EI_SO2TOSO4_UNIT.compare( "-" ) == 0 ) {
            /* Do nothing. Default unit */
        } else if ( Input_Opt.PARAMETER_EI_SO2TOSO4_UNIT.compare( "%" ) == 0 ) {
            /* Convert % to - */
            for ( i = 0; i < cases.size(); i++ )
                cases[i] /= 100.0;
        } else {
            std::cout << " Unknown unit for variable 'EI_SO2TOSO4': ";
            std::cout << Input_Opt.PARAMETER_EI_SO2TOSO4_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_EI_SO2TOSO4_UNIT = "-";

        counter += 1;
        y.push_back(Vector_1D(Input_Opt.SIMULATION_MCRUNS));
        for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
            y[counter-1][i] = cases[i];
        cases.clear();

        /* ======================================================================= */
        /* ---- EI_SOOT ---------------------------------------------------------- */
        /* ---- Accepted units are: g/kg_fuel (default)                            */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_EI_SOOT_RANGE ) {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( fRand(Input_Opt.PARAMETER_EI_SOOT[0], \
                                       Input_Opt.PARAMETER_EI_SOOT[1]) );
        } else {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( Input_Opt.PARAMETER_EI_SOOT[0] );
        }

        if ( ( Input_Opt.PARAMETER_EI_SOOT_UNIT.compare( "g/kg_fuel" ) == 0 ) || \
             ( Input_Opt.PARAMETER_EI_SOOT_UNIT.compare( "g/kg_f" )    == 0 ) || \
             ( Input_Opt.PARAMETER_EI_SOOT_UNIT.compare( "g/kg" )      == 0 ) ) {
            /* Do nothing. Default unit */
        } else if ( ( Input_Opt.PARAMETER_EI_SOOT_UNIT.compare( "mg/kg_fuel" ) == 0 ) || \
             ( Input_Opt.PARAMETER_EI_SOOT_UNIT.compare( "mg/kg_f" )    == 0 ) || \
             ( Input_Opt.PARAMETER_EI_SOOT_UNIT.compare( "mg/kg" )      == 0 ) ) {
            /* Convert mg to g */
            for ( i = 0; i < cases.size(); i++ )
                cases[i] /= 1.00E+03;
        } else {
            std::cout << " Unknown unit for variable 'EI_SOOT': ";
            std::cout << Input_Opt.PARAMETER_EI_SOOT_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_EI_SOOT_UNIT = "g/kg_fuel";

        counter += 1;
        y.push_back(Vector_1D(Input_Opt.SIMULATION_MCRUNS));
        for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
            y[counter-1][i] = cases[i];
        cases.clear();

        /* ======================================================================= */
        /* ---- EI_SOOTRAD ------------------------------------------------------- */
        /* ---- Accepted units are: m (default), nm                                */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_EI_SOOTRAD_RANGE ) {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( fRand(Input_Opt.PARAMETER_EI_SOOTRAD[0], \
                                       Input_Opt.PARAMETER_EI_SOOTRAD[1]) );
        } else {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( Input_Opt.PARAMETER_EI_SOOTRAD[0] );
        }

        if ( Input_Opt.PARAMETER_EI_SOOTRAD_UNIT.compare( "m" ) == 0 ) {
            /* Do nothing. Default unit */
        } else if ( Input_Opt.PARAMETER_EI_SOOTRAD_UNIT.compare( "nm" ) == 0 ) {
            /* Convert nm to m */
            for ( i = 0; i < cases.size(); i++ )
                cases[i] /= 1.00E+09;
        } else {
            std::cout << " Unknown unit for variable 'EI_SOOTRAD': ";
            std::cout << Input_Opt.PARAMETER_EI_SOOTRAD_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_EI_SOOTRAD_UNIT = "m";

        counter += 1;
        y.push_back(Vector_1D(Input_Opt.SIMULATION_MCRUNS));
        for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
            y[counter-1][i] = cases[i];
        cases.clear();

        /* ======================================================================= */
        /* ---- TOTAL FUEL FLOW -------------------------------------------------- */
        /* ---- Accepted units are: kg/s (default)                                 */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_FF_RANGE ) {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( fRand(Input_Opt.PARAMETER_FF[0], \
                                       Input_Opt.PARAMETER_FF[1]) );
        } else {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( Input_Opt.PARAMETER_FF[0] );
        }

        if ( ( Input_Opt.PARAMETER_FF_UNIT.compare( "kg_fuel/s" ) == 0 ) || \
             ( Input_Opt.PARAMETER_FF_UNIT.compare( "kg_f/s" )    == 0 ) || \
             ( Input_Opt.PARAMETER_FF_UNIT.compare( "kg/s" )      == 0 ) ) {
            /* Do nothing. Default unit */
        } else {
            std::cout << " Unknown unit for variable 'FF': ";
            std::cout << Input_Opt.PARAMETER_FF_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_FF_UNIT = "kg/s";

        counter += 1;
        y.push_back(Vector_1D(Input_Opt.SIMULATION_MCRUNS));
        for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
            y[counter-1][i] = cases[i];
        cases.clear();

        /* ======================================================================= */
        /* ---- AIRCRAFT MASS ---------------------------------------------------- */
        /* ---- Accepted units are: kg (default), tonnes                           */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_AMASS_RANGE ) {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( fRand(Input_Opt.PARAMETER_AMASS[0], \
                                       Input_Opt.PARAMETER_AMASS[1]) );
        } else {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( Input_Opt.PARAMETER_AMASS[0] );
        }

        if ( Input_Opt.PARAMETER_AMASS_UNIT.compare( "kg" ) == 0 ) {
            /* Do nothing. Default unit */
        } else if ( ( Input_Opt.PARAMETER_AMASS_UNIT.compare( "to" )      == 0 ) || \
                    ( Input_Opt.PARAMETER_AMASS_UNIT.compare( "tonne " )  == 0 ) || \
                    ( Input_Opt.PARAMETER_AMASS_UNIT.compare( "tonnes " ) == 0 ) ) {
            /* Convert tonnes to kg */
            for ( i = 0; i < cases.size(); i++ )
                cases[i] *= 1.00E+03;
        } else {
            std::cout << " Unknown unit for variable 'AMASS': ";
            std::cout << Input_Opt.PARAMETER_AMASS_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_AMASS_UNIT = "kg";

        counter += 1;
        y.push_back(Vector_1D(Input_Opt.SIMULATION_MCRUNS));
        for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
            y[counter-1][i] = cases[i];
        cases.clear();

        /* ======================================================================= */
        /* ---- FLIGHT SPEED ----------------------------------------------------- */
        /* ---- Accepted units are: m/s (default)                           */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_FSPEED_RANGE ) {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( fRand(Input_Opt.PARAMETER_FSPEED[0], \
                                       Input_Opt.PARAMETER_FSPEED[1]) );
        } else {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( Input_Opt.PARAMETER_FSPEED[0] );
        }

        if ( Input_Opt.PARAMETER_FSPEED_UNIT.compare( "m/s" ) == 0 ) {
            /* Do nothing. Default unit */
        } else {
            std::cout << " Unknown unit for variable 'FSPEED': ";
            std::cout << Input_Opt.PARAMETER_FSPEED_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_FSPEED_UNIT = "m/s";

        counter += 1;
        y.push_back(Vector_1D(Input_Opt.SIMULATION_MCRUNS));
        for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
            y[counter-1][i] = cases[i];
        cases.clear();

        /* ======================================================================= */
        /* ---- NUMBER OF ENGINES ------------------------------------------------ */
        /* ---- Accepted units are: 2/4 (default)                           */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_NUMENG_RANGE ) {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( fRand(Input_Opt.PARAMETER_NUMENG[0], \
                                       Input_Opt.PARAMETER_NUMENG[1]) );
        } else {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( Input_Opt.PARAMETER_NUMENG[0] );
        }

        if ( Input_Opt.PARAMETER_NUMENG_UNIT.compare( "2/4" ) == 0 ) {
            /* Do nothing. Default unit */
        } else {
            std::cout << " Unknown unit for variable 'NUMENG': ";
            std::cout << Input_Opt.PARAMETER_NUMENG_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_NUMENG_UNIT = "2/4";

        counter += 1;
        y.push_back(Vector_1D(Input_Opt.SIMULATION_MCRUNS));
        for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
            y[counter-1][i] = cases[i];
        cases.clear();

        /* ======================================================================= */
        /* ---- WINGSPAN --------------------------------------------------------- */
        /* ---- Accepted units are: m (default), ft                                */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_WINGSPAN_RANGE ) {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( fRand(Input_Opt.PARAMETER_WINGSPAN[0], \
                                       Input_Opt.PARAMETER_WINGSPAN[1]) );
        } else {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( Input_Opt.PARAMETER_WINGSPAN[0] );
        }

        if ( Input_Opt.PARAMETER_WINGSPAN_UNIT.compare( "m" ) == 0 ) {
            /* Do nothing. Default unit */
        } else if ( Input_Opt.PARAMETER_WINGSPAN_UNIT.compare( "ft" )      == 0 ) {
            /* Convert tonnes to kg */
            for ( i = 0; i < cases.size(); i++ )
                cases[i] /= 3.281;
        } else {
            std::cout << " Unknown unit for variable 'WINGSPAN': ";
            std::cout << Input_Opt.PARAMETER_WINGSPAN_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_AMASS_UNIT = "m";

        counter += 1;
        y.push_back(Vector_1D(Input_Opt.SIMULATION_MCRUNS));
        for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
            y[counter-1][i] = cases[i];
        cases.clear();

        /* ======================================================================= */
        /* ---- CORE EXIT TEMPERATURE -------------------------------------------- */
        /* ---- Accepted units are: K (default), C                                 */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_COREEXITTEMP_RANGE ) {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( fRand(Input_Opt.PARAMETER_COREEXITTEMP[0], \
                                       Input_Opt.PARAMETER_COREEXITTEMP[1]) );
        } else {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( Input_Opt.PARAMETER_COREEXITTEMP[0] );
        }

        if ( Input_Opt.PARAMETER_COREEXITTEMP_UNIT.compare( "K" ) == 0 ) {
            /* Do nothing. Default unit */
        } else if ( Input_Opt.PARAMETER_COREEXITTEMP_UNIT.compare( "C" )      == 0 ) {
            /* Convert tonnes to kg */
            for ( i = 0; i < cases.size(); i++ )
                cases[i] += 273.15;
        } else {
            std::cout << " Unknown unit for variable 'COREEXITTEMP': ";
            std::cout << Input_Opt.PARAMETER_COREEXITTEMP_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_COREEXITTEMP_UNIT = "K";

        counter += 1;
        y.push_back(Vector_1D(Input_Opt.SIMULATION_MCRUNS));
        for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
            y[counter-1][i] = cases[i];
        cases.clear();

        /* ======================================================================= */
        /* ---- EXIT BYPASS AREA ------------------------------------------------- */
        /* ---- Accepted units are: m^2 (default)                                 */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_BYPASSAREA_RANGE ) {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( fRand(Input_Opt.PARAMETER_BYPASSAREA[0], \
                                       Input_Opt.PARAMETER_BYPASSAREA[1]) );
        } else {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( Input_Opt.PARAMETER_BYPASSAREA[0] );
        }

        if ( Input_Opt.PARAMETER_BYPASSAREA_UNIT.compare( "m^2" ) == 0 ) {
            /* Do nothing. Default unit */
        } else {
            std::cout << " Unknown unit for variable 'BYPASSAREA': ";
            std::cout << Input_Opt.PARAMETER_BYPASSAREA_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_BYPASSAREA_UNIT = "m^2";

        counter += 1;
        y.push_back(Vector_1D(Input_Opt.SIMULATION_MCRUNS));
        for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
            y[counter-1][i] = cases[i];
        cases.clear();


        /* ======================================================================= */
        /* ---- BACKGROUND NOX MIXING RATIO -------------------------------------- */
        /* ---- Accepted units are: ppb (default), ppt, ppm                        */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_BACKG_NOX_RANGE ) {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( fRand(Input_Opt.PARAMETER_BACKG_NOX[0], \
                                       Input_Opt.PARAMETER_BACKG_NOX[1]) );
        } else {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( Input_Opt.PARAMETER_BACKG_NOX[0] );
        }

        if ( Input_Opt.PARAMETER_BACKG_NOX_UNIT.compare( "ppb" ) == 0 ) {
            /* Do nothing. Default unit */
        } else if ( Input_Opt.PARAMETER_BACKG_NOX_UNIT.compare( "ppt" ) == 0 ) {
            /* Convert ppt to ppb */
            for ( i = 0; i < cases.size(); i++ )
                cases[i] *= 1.00E-03;
        } else if ( Input_Opt.PARAMETER_BACKG_NOX_UNIT.compare( "ppm" ) == 0 ) {
            /* Convert ppm to ppb */
            for ( i = 0; i < cases.size(); i++ )
                cases[i] *= 1.00E+03;
        } else {
            std::cout << " Unknown unit for variable 'BACKG_NOX': ";
            std::cout << Input_Opt.PARAMETER_BACKG_NOX_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_BACKG_NOX_UNIT = "ppb";

        counter += 1;
        y.push_back(Vector_1D(Input_Opt.SIMULATION_MCRUNS));
        for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
            y[counter-1][i] = cases[i];
        cases.clear();

        /* ======================================================================= */
        /* ---- BACKGROUND HNO3 MIXING RATIO ------------------------------------- */
        /* ---- Accepted units are: ppb (default), ppt, ppm                        */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_BACKG_HNO3_RANGE ) {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( fRand(Input_Opt.PARAMETER_BACKG_HNO3[0], \
                                       Input_Opt.PARAMETER_BACKG_HNO3[1]) );
        } else {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( Input_Opt.PARAMETER_BACKG_HNO3[0] );
        }

        if ( Input_Opt.PARAMETER_BACKG_HNO3_UNIT.compare( "ppb" ) == 0 ) {
            /* Do nothing. Default unit */
        } else if ( Input_Opt.PARAMETER_BACKG_HNO3_UNIT.compare( "ppt" ) == 0 ) {
            /* Convert ppt to ppb */
            for ( i = 0; i < cases.size(); i++ )
                cases[i] *= 1.00E-03;
        } else if ( Input_Opt.PARAMETER_BACKG_HNO3_UNIT.compare( "ppm" ) == 0 ) {
            /* Convert ppm to ppb */
            for ( i = 0; i < cases.size(); i++ )
                cases[i] *= 1.00E+03;
        } else {
            std::cout << " Unknown unit for variable 'BACKG_HNO3': ";
            std::cout << Input_Opt.PARAMETER_BACKG_HNO3_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_BACKG_HNO3_UNIT = "HNO3";

        counter += 1;
        y.push_back(Vector_1D(Input_Opt.SIMULATION_MCRUNS));
        for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
            y[counter-1][i] = cases[i];
        cases.clear();

        /* ======================================================================= */
        /* ---- BACKGROUND O3 MIXING RATIO --------------------------------------- */
        /* ---- Accepted units are: ppb (default), ppt, ppm                        */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_BACKG_O3_RANGE ) {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( fRand(Input_Opt.PARAMETER_BACKG_O3[0], \
                                       Input_Opt.PARAMETER_BACKG_O3[1]) );
        } else {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( Input_Opt.PARAMETER_BACKG_O3[0] );
        }

        if ( Input_Opt.PARAMETER_BACKG_O3_UNIT.compare( "ppb" ) == 0 ) {
            /* Do nothing. Default unit */
        } else if ( Input_Opt.PARAMETER_BACKG_O3_UNIT.compare( "ppt" ) == 0 ) {
            /* Convert ppt to ppb */
            for ( i = 0; i < cases.size(); i++ )
                cases[i] *= 1.00E-03;
        } else if ( Input_Opt.PARAMETER_BACKG_O3_UNIT.compare( "ppm" ) == 0 ) {
            /* Convert ppm to ppb */
            for ( i = 0; i < cases.size(); i++ )
                cases[i] *= 1.00E+03;
        } else {
            std::cout << " Unknown unit for variable 'BACKG_O3': ";
            std::cout << Input_Opt.PARAMETER_BACKG_O3_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_BACKG_O3_UNIT = "ppb";

        counter += 1;
        y.push_back(Vector_1D(Input_Opt.SIMULATION_MCRUNS));
        for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
            y[counter-1][i] = cases[i];
        cases.clear();

        /* ======================================================================= */
        /* ---- BACKGROUND CO MIXING RATIO --------------------------------------- */
        /* ---- Accepted units are: ppb (default), ppt, ppm                        */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_BACKG_CO_RANGE ) {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( fRand(Input_Opt.PARAMETER_BACKG_CO[0], \
                                       Input_Opt.PARAMETER_BACKG_CO[1]) );
        } else {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( Input_Opt.PARAMETER_BACKG_CO[0] );
        }

        if ( Input_Opt.PARAMETER_BACKG_CO_UNIT.compare( "ppb" ) == 0 ) {
            /* Do nothing. Default unit */
        } else if ( Input_Opt.PARAMETER_BACKG_CO_UNIT.compare( "ppt" ) == 0 ) {
            /* Convert ppt to ppb */
            for ( i = 0; i < cases.size(); i++ )
                cases[i] *= 1.00E-03;
        } else if ( Input_Opt.PARAMETER_BACKG_CO_UNIT.compare( "ppm" ) == 0 ) {
            /* Convert ppm to ppb */
            for ( i = 0; i < cases.size(); i++ )
                cases[i] *= 1.00E+03;
        } else {
            std::cout << " Unknown unit for variable 'BACKG_CO': ";
            std::cout << Input_Opt.PARAMETER_BACKG_CO_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_BACKG_CO_UNIT = "ppb";

        counter += 1;
        y.push_back(Vector_1D(Input_Opt.SIMULATION_MCRUNS));
        for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
            y[counter-1][i] = cases[i];
        cases.clear();

        /* ======================================================================= */
        /* ---- BACKGROUND CH4 MIXING RATIO -------------------------------------- */
        /* ---- Accepted units are: ppb (default), ppt, ppm                        */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_BACKG_CH4_RANGE ) {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( fRand(Input_Opt.PARAMETER_BACKG_CH4[0], \
                                       Input_Opt.PARAMETER_BACKG_CH4[1]) );
        } else {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( Input_Opt.PARAMETER_BACKG_CH4[0] );
        }

        if ( Input_Opt.PARAMETER_BACKG_CH4_UNIT.compare( "ppb" ) == 0 ) {
            /* Do nothing. Default unit */
        } else if ( Input_Opt.PARAMETER_BACKG_CH4_UNIT.compare( "ppt" ) == 0 ) {
            /* Convert ppt to ppb */
            for ( i = 0; i < cases.size(); i++ )
                cases[i] *= 1.00E-03;
        } else if ( Input_Opt.PARAMETER_BACKG_CH4_UNIT.compare( "ppm" ) == 0 ) {
            /* Convert ppm to ppb */
            for ( i = 0; i < cases.size(); i++ )
                cases[i] *= 1.00E+03;
        } else {
            std::cout << " Unknown unit for variable 'BACKG_CH4': ";
            std::cout << Input_Opt.PARAMETER_BACKG_CH4_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_BACKG_CH4_UNIT = "ppb";

        counter += 1;
        y.push_back(Vector_1D(Input_Opt.SIMULATION_MCRUNS));
        for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
            y[counter-1][i] = cases[i];
        cases.clear();

        /* ======================================================================= */
        /* ---- BACKGROUND SO2 MIXING RATIO -------------------------------------- */
        /* ---- Assumed units are: ppb (default), ppt, ppm                        */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_BACKG_SO2_RANGE ) {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( fRand(Input_Opt.PARAMETER_BACKG_SO2[0], \
                                       Input_Opt.PARAMETER_BACKG_SO2[1]) );
        } else {
            for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
                cases.push_back( Input_Opt.PARAMETER_BACKG_SO2[0] );
        }

        if ( Input_Opt.PARAMETER_BACKG_SO2_UNIT.compare( "ppb" ) == 0 ) {
            /* Do nothing. Default unit */
        } else if ( Input_Opt.PARAMETER_BACKG_SO2_UNIT.compare( "ppt" ) == 0 ) {
            /* Convert ppt to ppb */
            for ( i = 0; i < cases.size(); i++ )
                cases[i] *= 1.00E-03;
        } else if ( Input_Opt.PARAMETER_BACKG_SO2_UNIT.compare( "ppm" ) == 0 ) {
            /* Convert ppm to ppb */
            for ( i = 0; i < cases.size(); i++ )
                cases[i] *= 1.00E+03;
        } else {
            std::cout << " Unknown unit for variable 'BACKG_SO2': ";
            std::cout << Input_Opt.PARAMETER_BACKG_SO2_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_BACKG_SO2_UNIT = "ppb";

        counter += 1;
        y.push_back(Vector_1D(Input_Opt.SIMULATION_MCRUNS));
        for ( i = 0; i < Input_Opt.SIMULATION_MCRUNS; i++ )
            y[counter-1][i] = cases[i];

        return y;

    } else if ( Input_Opt.PARAMETER_FILEINPUT ) {

        /* ======================================================================= */
        /* ---- PLUME PROCESSING TIME ( SIMULATION TIME ) ------------------------ */
        /* ---- Assumed units are: hr                                              */
        /* ======================================================================= */

        y.push_back( Vector_1D( Input_Opt.PARAMETER_FILECASES ) );
        if ( Input_Opt.PARAMETER_PLUMEPROCESS.size() == Input_Opt.PARAMETER_FILECASES ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_PLUMEPROCESS[i];
        } else if ( Input_Opt.PARAMETER_PLUMEPROCESS.size() == 1 ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_PLUMEPROCESS[0];
        } else {
            std::cout << " In CombVec:";
            std::cout << " PARAMETER_PLUMEPROCESS has the wrong shape: ";
            std::cout << Input_Opt.PARAMETER_PLUMEPROCESS.size() << std::endl;
        }

        counter += 1;

        /* ======================================================================= */
        /* ---- TEMPERATURE ------------------------------------------------------ */
        /* ---- Assumed units are: Kelvin                                          */
        /* ======================================================================= */

        y.push_back( Vector_1D( Input_Opt.PARAMETER_FILECASES ) );
        if ( Input_Opt.PARAMETER_TEMPERATURE.size() == Input_Opt.PARAMETER_FILECASES ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_TEMPERATURE[i];
        } else if ( Input_Opt.PARAMETER_TEMPERATURE.size() == 1 ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_TEMPERATURE[0];
        } else {
            std::cout << " In CombVec:";
            std::cout << " PARAMETER_TEMPERATURE has the wrong shape: ";
            std::cout << Input_Opt.PARAMETER_TEMPERATURE.size() << std::endl;
        }

        counter += 1;

        /* ======================================================================= */
        /* ---- RELATIVE HUMIDITY ------------------------------------------------ */
        /* ---- Assumed units are: %                                               */
        /* ======================================================================= */

        y.push_back( Vector_1D( Input_Opt.PARAMETER_FILECASES ) );
        if ( Input_Opt.PARAMETER_RHW.size() == Input_Opt.PARAMETER_FILECASES ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_RHW[i];
        } else if ( Input_Opt.PARAMETER_RHW.size() == 1 ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_RHW[0];
        } else {
            std::cout << " In CombVec:";
            std::cout << " PARAMETER_RHW has the wrong shape: ";
            std::cout << Input_Opt.PARAMETER_RHW.size() << std::endl;
        }

        counter += 1;

        /* ======================================================================= */
        /* ---- HORIZONTAL DIFFUSION PARAMETER ----------------------------------- */
        /* ---- Assumed units are: m^2/s                                           */
        /* ======================================================================= */

        y.push_back( Vector_1D( Input_Opt.PARAMETER_FILECASES ) );
        if ( Input_Opt.PARAMETER_DH.size() == Input_Opt.PARAMETER_FILECASES ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_DH[i];
        } else if ( Input_Opt.PARAMETER_DH.size() == 1 ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_DH[0];
        } else {
            std::cout << " In CombVec:";
            std::cout << " PARAMETER_DH has the wrong shape: ";
            std::cout << Input_Opt.PARAMETER_DH.size() << std::endl;
        }

        counter += 1;

        /* ======================================================================= */
        /* ---- VERTICAL DIFFUSION PARAMETER ------------------------------------- */
        /* ---- Assumed units are: m^2/s                                           */
        /* ======================================================================= */

        y.push_back( Vector_1D( Input_Opt.PARAMETER_FILECASES ) );
        if ( Input_Opt.PARAMETER_DV.size() == Input_Opt.PARAMETER_FILECASES ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_DV[i];
        } else if ( Input_Opt.PARAMETER_DV.size() == 1 ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_DV[0];
        } else {
            std::cout << " In CombVec:";
            std::cout << " PARAMETER_DV has the wrong shape: ";
            std::cout << Input_Opt.PARAMETER_DV.size() << std::endl;
        }

        counter += 1;

        /* ======================================================================= */
        /* ---- SHEAR ------------------------------------------------------------ */
        /* ---- Assumed units are: 1/s                                             */
        /* ======================================================================= */

        y.push_back( Vector_1D( Input_Opt.PARAMETER_FILECASES ) );
        if ( Input_Opt.PARAMETER_SHEAR.size() == Input_Opt.PARAMETER_FILECASES ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_SHEAR[i];
        } else if ( Input_Opt.PARAMETER_SHEAR.size() == 1 ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_SHEAR[0];
        } else {
            std::cout << " In CombVec:";
            std::cout << " PARAMETER_SHEAR has the wrong shape: ";
            std::cout << Input_Opt.PARAMETER_SHEAR.size() << std::endl;
        }

        counter += 1;

        /* ======================================================================= */
        /* ---- LONGITUDE -------------------------------------------------------- */
        /* ---- Assumed units are: degree                                          */
        /* ======================================================================= */

        y.push_back( Vector_1D( Input_Opt.PARAMETER_FILECASES ) );
        if ( Input_Opt.PARAMETER_LONGITUDE.size() == Input_Opt.PARAMETER_FILECASES ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_LONGITUDE[i];
        } else if ( Input_Opt.PARAMETER_LONGITUDE.size() == 1 ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_LONGITUDE[0];
        } else {
            std::cout << " In CombVec:";
            std::cout << " PARAMETER_LONGITUDE has the wrong shape: ";
            std::cout << Input_Opt.PARAMETER_LONGITUDE.size() << std::endl;
        }

        counter += 1;

        /* ======================================================================= */
        /* ---- LATITUDE --------------------------------------------------------- */
        /* ---- Assumed units are: degree                                          */
        /* ======================================================================= */

        y.push_back( Vector_1D( Input_Opt.PARAMETER_FILECASES ) );
        if ( Input_Opt.PARAMETER_LATITUDE.size() == Input_Opt.PARAMETER_FILECASES ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_LATITUDE[i];
        } else if ( Input_Opt.PARAMETER_LATITUDE.size() == 1 ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_LATITUDE[0];
        } else {
            std::cout << " In CombVec:";
            std::cout << " PARAMETER_LATITUDE has the wrong shape: ";
            std::cout << Input_Opt.PARAMETER_LATITUDE.size() << std::endl;
        }

        counter += 1;

        /* ======================================================================= */
        /* ---- PRESSURE --------------------------------------------------------- */
        /* ---- Assumed units are: hPa                                             */
        /* ======================================================================= */

        y.push_back( Vector_1D( Input_Opt.PARAMETER_FILECASES ) );
        if ( Input_Opt.PARAMETER_PRESSURE.size() == Input_Opt.PARAMETER_FILECASES ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_PRESSURE[i] * 100.0;
        } else if ( Input_Opt.PARAMETER_PRESSURE.size() == 1 ) {
            if ( Input_Opt.PARAMETER_PRESSURE_UNIT.compare( "Pa" ) == 0 ) {
                /* Do nothing. Default unit */
                for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                    y[counter-1][i] = Input_Opt.PARAMETER_PRESSURE[0];
            } else if ( Input_Opt.PARAMETER_PRESSURE_UNIT.compare( "hPa" ) == 0 ) {
                for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                    y[counter-1][i] = Input_Opt.PARAMETER_PRESSURE[0] * 100.0;
            } else {
                std::cout << " Unknown unit for variable 'Pressure': ";
                std::cout << Input_Opt.PARAMETER_PRESSURE_UNIT << std::endl;
                exit(1);
            }
        } else {
            std::cout << " In CombVec:";
            std::cout << " PARAMETER_PRESSURE has the wrong shape: ";
            std::cout << Input_Opt.PARAMETER_PRESSURE.size() << std::endl;
        }

        counter += 1;

        /* ======================================================================= */
        /* ---- EMISSION DAY ----------------------------------------------------- */
        /* ---- Assumed units are: 1-365                                           */
        /* ======================================================================= */

        y.push_back( Vector_1D( Input_Opt.PARAMETER_FILECASES ) );
        if ( Input_Opt.PARAMETER_EDAY.size() == Input_Opt.PARAMETER_FILECASES ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_EDAY[i];
        } else if ( Input_Opt.PARAMETER_EDAY.size() == 1 ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_EDAY[0];
        } else {
            std::cout << " In CombVec:";
            std::cout << " PARAMETER_EDAY has the wrong shape: ";
            std::cout << Input_Opt.PARAMETER_EDAY.size() << std::endl;
        }

        counter += 1;

        /* ======================================================================= */
        /* ---- EMISSION TIME ---------------------------------------------------- */
        /* ---- Assumed units are: 0-24                                            */
        /* ======================================================================= */

        y.push_back( Vector_1D( Input_Opt.PARAMETER_FILECASES ) );
        if ( Input_Opt.PARAMETER_ETIME.size() == Input_Opt.PARAMETER_FILECASES ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_ETIME[i];
        } else if ( Input_Opt.PARAMETER_ETIME.size() == 1 ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_ETIME[0];
        } else {
            std::cout << " In CombVec:";
            std::cout << " PARAMETER_ETIME has the wrong shape: ";
            std::cout << Input_Opt.PARAMETER_ETIME.size() << std::endl;
        }

        counter += 1;

        /* ======================================================================= */
        /* ---- EI_NOX ----------------------------------------------------------- */
        /* ---- Assumed units are: g/kg_fuel                                       */
        /* ======================================================================= */

        y.push_back( Vector_1D( Input_Opt.PARAMETER_FILECASES ) );
        if ( Input_Opt.PARAMETER_EI_NOX.size() == Input_Opt.PARAMETER_FILECASES ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_EI_NOX[i];
        } else if ( Input_Opt.PARAMETER_EI_NOX.size() == 1 ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_EI_NOX[0];
        } else {
            std::cout << " In CombVec:";
            std::cout << " PARAMETER_EI_NOX has the wrong shape: ";
            std::cout << Input_Opt.PARAMETER_EI_NOX.size() << std::endl;
        }

        counter += 1;

        /* ======================================================================= */
        /* ---- EI_CO ------------------------------------------------------------ */
        /* ---- Assumed units are: g/kg_fuel                                       */
        /* ======================================================================= */

        y.push_back( Vector_1D( Input_Opt.PARAMETER_FILECASES ) );
        if ( Input_Opt.PARAMETER_EI_CO.size() == Input_Opt.PARAMETER_FILECASES ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_EI_CO[i];
        } else if ( Input_Opt.PARAMETER_EI_CO.size() == 1 ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_EI_CO[0];
        } else {
            std::cout << " In CombVec:";
            std::cout << " PARAMETER_EI_CO has the wrong shape: ";
            std::cout << Input_Opt.PARAMETER_EI_CO.size() << std::endl;
        }

        counter += 1;

        /* ======================================================================= */
        /* ---- EI_UHC ----------------------------------------------------------- */
        /* ---- Assumed units are: g/kg_fuel                                       */
        /* ======================================================================= */

        y.push_back( Vector_1D( Input_Opt.PARAMETER_FILECASES ) );
        if ( Input_Opt.PARAMETER_EI_UHC.size() == Input_Opt.PARAMETER_FILECASES ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_EI_UHC[i];
        } else if ( Input_Opt.PARAMETER_EI_UHC.size() == 1 ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_EI_UHC[0];
        } else {
            std::cout << " In CombVec:";
            std::cout << " PARAMETER_EI_UHC has the wrong shape: ";
            std::cout << Input_Opt.PARAMETER_EI_UHC.size() << std::endl;
        }

        counter += 1;

        /* ======================================================================= */
        /* ---- EI_SO2 ----------------------------------------------------------- */
        /* ---- Assumed units are: g/kg_fuel                                       */
        /* ======================================================================= */

        y.push_back( Vector_1D( Input_Opt.PARAMETER_FILECASES ) );
        if ( Input_Opt.PARAMETER_EI_SO2.size() == Input_Opt.PARAMETER_FILECASES ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_EI_SO2[i];
        } else if ( Input_Opt.PARAMETER_EI_SO2.size() == 1 ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_EI_SO2[0];
        } else {
            std::cout << " In CombVec:";
            std::cout << " PARAMETER_EI_SO2 has the wrong shape: ";
            std::cout << Input_Opt.PARAMETER_EI_SO2.size() << std::endl;
        }

        counter += 1;

        /* ======================================================================= */
        /* ---- EI_SO2TOSO4 ------------------------------------------------------ */
        /* ---- Assumed units are: -                                               */
        /* ======================================================================= */

        y.push_back( Vector_1D( Input_Opt.PARAMETER_FILECASES ) );
        if ( Input_Opt.PARAMETER_EI_SO2TOSO4.size() == Input_Opt.PARAMETER_FILECASES ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_EI_SO2TOSO4[i];
        } else if ( Input_Opt.PARAMETER_EI_SO2TOSO4.size() == 1 ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_EI_SO2TOSO4[0];
        } else {
            std::cout << " In CombVec:";
            std::cout << " PARAMETER_EI_SO2TOSO4 has the wrong shape: ";
            std::cout << Input_Opt.PARAMETER_EI_SO2TOSO4.size() << std::endl;
        }

        counter += 1;

        /* ======================================================================= */
        /* ---- EI_SOOT ---------------------------------------------------------- */
        /* ---- Assumed units are: g/kg_fuel                                       */
        /* ======================================================================= */

        y.push_back( Vector_1D( Input_Opt.PARAMETER_FILECASES ) );
        if ( Input_Opt.PARAMETER_EI_SOOT.size() == Input_Opt.PARAMETER_FILECASES ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_EI_SOOT[i];
        } else if ( Input_Opt.PARAMETER_EI_SOOT.size() == 1 ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_EI_SOOT[0];
        } else {
            std::cout << " In CombVec:";
            std::cout << " PARAMETER_EI_SOOT has the wrong shape: ";
            std::cout << Input_Opt.PARAMETER_EI_SOOT.size() << std::endl;
        }

        counter += 1;

        /* ======================================================================= */
        /* ---- EI_SOOTRAD ------------------------------------------------------- */
        /* ---- Assumed units are: m                                               */
        /* ======================================================================= */

        y.push_back( Vector_1D( Input_Opt.PARAMETER_FILECASES ) );
        if ( Input_Opt.PARAMETER_EI_SOOTRAD.size() == Input_Opt.PARAMETER_FILECASES ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_EI_SOOTRAD[i];
        } else if ( Input_Opt.PARAMETER_EI_SOOTRAD.size() == 1 ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_EI_SOOTRAD[0];
        } else {
            std::cout << " In CombVec:";
            std::cout << " PARAMETER_EI_SOOTRAD has the wrong shape: ";
            std::cout << Input_Opt.PARAMETER_EI_SOOTRAD.size() << std::endl;
        }

        counter += 1;

        /* ======================================================================= */
        /* ---- TOTAL FUEL FLOW -------------------------------------------------- */
        /* ---- Assumed units are: kg/s                                            */
        /* ======================================================================= */

        y.push_back( Vector_1D( Input_Opt.PARAMETER_FILECASES ) );
        if ( Input_Opt.PARAMETER_FF.size() == Input_Opt.PARAMETER_FILECASES ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_FF[i];
        } else if ( Input_Opt.PARAMETER_FF.size() == 1 ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_FF[0];
        } else {
            std::cout << " In CombVec:";
            std::cout << " PARAMETER_FF has the wrong shape: ";
            std::cout << Input_Opt.PARAMETER_FF.size() << std::endl;
        }

        counter += 1;

        /* ======================================================================= */
        /* ---- FLIGHT SPEED ----------------------------------------------------- */
        /* ---- Assumed units are: m/s                                             */
        /* ======================================================================= */

        y.push_back( Vector_1D( Input_Opt.PARAMETER_FILECASES ) );
        if ( Input_Opt.PARAMETER_FSPEED.size() == Input_Opt.PARAMETER_FILECASES ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_FSPEED[i];
        } else if ( Input_Opt.PARAMETER_FSPEED.size() == 1 ) {
            if ( Input_Opt.PARAMETER_FSPEED_UNIT.compare( "m/s" ) == 0 ) {
                /* Do nothing. Default unit */
                for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                    y[counter-1][i] = Input_Opt.PARAMETER_FSPEED[0];
            } else {
                std::cout << " Unknown unit for variable 'FSPEED': ";
                std::cout << Input_Opt.PARAMETER_FSPEED_UNIT << std::endl;
                exit(1);
            }

            /* Updating unit now that conversion has been taken care of */
            Input_Opt.PARAMETER_FSPEED_UNIT = "m/s";
        } else {
            std::cout << " In CombVec:";
            std::cout << " PARAMETER_FSPEED has the wrong shape: ";
            std::cout << Input_Opt.PARAMETER_FSPEED.size() << std::endl;
        }

        counter += 1;

        /* ======================================================================= */
        /* ---- NUMBER OF ENGINES ------------------------------------------------ */
        /* ---- Assumed units are: 2/4                                             */
        /* ======================================================================= */

        y.push_back( Vector_1D( Input_Opt.PARAMETER_FILECASES ) );
        if ( Input_Opt.PARAMETER_NUMENG.size() == Input_Opt.PARAMETER_FILECASES ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_NUMENG[i];
        } else if ( Input_Opt.PARAMETER_NUMENG.size() == 1 ) {
            if ( Input_Opt.PARAMETER_NUMENG_UNIT.compare( "2/4" ) == 0 ) {
                /* Do nothing. Default unit */
                for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                    y[counter-1][i] = Input_Opt.PARAMETER_FSPEED[0];
            } else {
                std::cout << " Unknown unit for variable 'NUMENG': ";
                std::cout << Input_Opt.PARAMETER_NUMENG_UNIT << std::endl;
                exit(1);
            }

            /* Updating unit now that conversion has been taken care of */
            Input_Opt.PARAMETER_NUMENG_UNIT = "2/4";
        } else {
            std::cout << " In CombVec:";
            std::cout << " PARAMETER_NUMENG has the wrong shape: ";
            std::cout << Input_Opt.PARAMETER_NUMENG.size() << std::endl;
        }

        counter += 1;

        /* ======================================================================= */
        /* ---- WINGSPAN --------------------------------------------------------- */
        /* ---- Assumed units are: m                                               */
        /* ======================================================================= */

        y.push_back( Vector_1D( Input_Opt.PARAMETER_FILECASES ) );
        if ( Input_Opt.PARAMETER_WINGSPAN.size() == Input_Opt.PARAMETER_FILECASES ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_WINGSPAN[i];
        } else if ( Input_Opt.PARAMETER_WINGSPAN.size() == 1 ) {
            if ( Input_Opt.PARAMETER_WINGSPAN_UNIT.compare( "m" ) == 0 ) {
                /* Do nothing. Default unit */
                for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                    y[counter-1][i] = Input_Opt.PARAMETER_FSPEED[0];
            } else if ( Input_Opt.PARAMETER_WINGSPAN_UNIT.compare( "ft" ) == 0 ) {
                /* Do nothing. Default unit */
                for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                    y[counter-1][i] = Input_Opt.PARAMETER_FSPEED[0] / 3.281;
            } else {
                std::cout << " Unknown unit for variable 'WINGSPAN': ";
                std::cout << Input_Opt.PARAMETER_WINGSPAN_UNIT << std::endl;
                exit(1);
            }

            /* Updating unit now that conversion has been taken care of */
            Input_Opt.PARAMETER_WINGSPAN_UNIT = "m";
        } else {
            std::cout << " In CombVec:";
            std::cout << " PARAMETER_WINGSPAN has the wrong shape: ";
            std::cout << Input_Opt.PARAMETER_WINGSPAN.size() << std::endl;
        }

        counter += 1;

        /* ======================================================================= */
        /* ---- CORE EXIT TEMPERATURE -------------------------------------------- */
        /* ---- Assumed units are: K                                               */
        /* ======================================================================= */

        y.push_back( Vector_1D( Input_Opt.PARAMETER_FILECASES ) );
        if ( Input_Opt.PARAMETER_COREEXITTEMP.size() == Input_Opt.PARAMETER_FILECASES ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_COREEXITTEMP[i];
        } else if ( Input_Opt.PARAMETER_COREEXITTEMP.size() == 1 ) {
            if ( Input_Opt.PARAMETER_COREEXITTEMP_UNIT.compare( "K" ) == 0 ) {
                /* Do nothing. Default unit */
                for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                    y[counter-1][i] = Input_Opt.PARAMETER_COREEXITTEMP[0];
            } else if ( Input_Opt.PARAMETER_COREEXITTEMP_UNIT.compare( "C" ) == 0 ) {
                /* Do nothing. Default unit */
                for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                    y[counter-1][i] = Input_Opt.PARAMETER_COREEXITTEMP[0] + 273.15;
            } else {
                std::cout << " Unknown unit for variable 'COREEXITTEMP': ";
                std::cout << Input_Opt.PARAMETER_COREEXITTEMP_UNIT << std::endl;
                exit(1);
            }

            /* Updating unit now that conversion has been taken care of */
            Input_Opt.PARAMETER_COREEXITTEMP_UNIT = "K";
        } else {
            std::cout << " In CombVec:";
            std::cout << " PARAMETER_COREEXITTEMP has the wrong shape: ";
            std::cout << Input_Opt.PARAMETER_COREEXITTEMP.size() << std::endl;
        }

        counter += 1;

        /* ======================================================================= */
        /* ---- EXIT BYPASS AREA -------------------------------------------- */
        /* ---- Assumed units are: m^2                                               */
        /* ======================================================================= */

        y.push_back( Vector_1D( Input_Opt.PARAMETER_FILECASES ) );
        if ( Input_Opt.PARAMETER_BYPASSAREA.size() == Input_Opt.PARAMETER_FILECASES ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_BYPASSAREA[i];
        } else if ( Input_Opt.PARAMETER_BYPASSAREA.size() == 1 ) {
            if ( Input_Opt.PARAMETER_BYPASSAREA_UNIT.compare( "m^2" ) == 0 ) {
                /* Do nothing. Default unit */
                for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                    y[counter-1][i] = Input_Opt.PARAMETER_BYPASSAREA[0];
            } else {
                std::cout << " Unknown unit for variable 'BYPASSAREA': ";
                std::cout << Input_Opt.PARAMETER_BYPASSAREA_UNIT << std::endl;
                exit(1);
            }

            /* Updating unit now that conversion has been taken care of */
            Input_Opt.PARAMETER_BYPASSAREA_UNIT = "m^2";
        } else {
            std::cout << " In CombVec:";
            std::cout << " PARAMETER_BYPASSAREA has the wrong shape: ";
            std::cout << Input_Opt.PARAMETER_BYPASSAREA.size() << std::endl;
        }

        counter += 1;


        /* ======================================================================= */
        /* ---- AIRCRAFT MASS ---------------------------------------------------- */
        /* ---- Assumed units are: kg                                              */
        /* ======================================================================= */

        y.push_back( Vector_1D( Input_Opt.PARAMETER_FILECASES ) );
        if ( Input_Opt.PARAMETER_AMASS.size() == Input_Opt.PARAMETER_FILECASES ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_AMASS[i];
        } else if ( Input_Opt.PARAMETER_AMASS.size() == 1 ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_AMASS[0];
        } else {
            std::cout << " In CombVec:";
            std::cout << " PARAMETER_AMASS has the wrong shape: ";
            std::cout << Input_Opt.PARAMETER_AMASS.size() << std::endl;
        }

        counter += 1;

        /* ======================================================================= */
        /* ---- BACKGROUND NOX MIXING RATIO -------------------------------------- */
        /* ---- Assumed units are: ppt                                             */
        /* ======================================================================= */

        y.push_back( Vector_1D( Input_Opt.PARAMETER_FILECASES ) );
        if ( Input_Opt.PARAMETER_BACKG_NOX.size() == Input_Opt.PARAMETER_FILECASES ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_BACKG_NOX[i] * 1.00E-03;
        } else if ( Input_Opt.PARAMETER_BACKG_NOX.size() == 1 ) {
            if ( Input_Opt.PARAMETER_BACKG_NOX_UNIT.compare( "ppb" ) == 0 ) {
                /* Do nothing. Default unit */
                for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                    y[counter-1][i] = Input_Opt.PARAMETER_BACKG_NOX[0];
            } else if ( Input_Opt.PARAMETER_BACKG_NOX_UNIT.compare( "ppt" ) == 0 ) {
                /* Convert ppt to ppb */
                for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                    y[counter-1][i] = Input_Opt.PARAMETER_BACKG_NOX[0] * 1.00E-03;
            } else if ( Input_Opt.PARAMETER_BACKG_NOX_UNIT.compare( "ppm" ) == 0 ) {
                /* Convert ppm to ppb */
                for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                    y[counter-1][i] = Input_Opt.PARAMETER_BACKG_NOX[0] * 1.00E+03;
            } else {
                std::cout << " Unknown unit for variable 'BACKG_NOX': ";
                std::cout << Input_Opt.PARAMETER_BACKG_NOX_UNIT << std::endl;
                exit(1);
            }

            /* Updating unit now that conversion has been taken care of */
            Input_Opt.PARAMETER_BACKG_NOX_UNIT = "ppb";
        } else {
            std::cout << " In CombVec:";
            std::cout << " PARAMETER_BACKG_NOX has the wrong shape: ";
            std::cout << Input_Opt.PARAMETER_BACKG_NOX.size() << std::endl;
        }

        counter += 1;

        /* ======================================================================= */
        /* ---- BACKGROUND HNO3 MIXING RATIO ------------------------------------- */
        /* ---- Assumed units are: ppt                                             */
        /* ======================================================================= */

        y.push_back( Vector_1D( Input_Opt.PARAMETER_FILECASES ) );
        if ( Input_Opt.PARAMETER_BACKG_HNO3.size() == Input_Opt.PARAMETER_FILECASES ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_BACKG_HNO3[i] * 1.00E-03;
        } else if ( Input_Opt.PARAMETER_BACKG_HNO3.size() == 1 ) {
            if ( Input_Opt.PARAMETER_BACKG_HNO3_UNIT.compare( "ppb" ) == 0 ) {
                /* Do nothing. Default unit */
                for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                    y[counter-1][i] = Input_Opt.PARAMETER_BACKG_HNO3[0];
            } else if ( Input_Opt.PARAMETER_BACKG_HNO3_UNIT.compare( "ppt" ) == 0 ) {
                /* Convert ppt to ppb */
                for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                    y[counter-1][i] = Input_Opt.PARAMETER_BACKG_HNO3[0] * 1.00E-03;
            } else if ( Input_Opt.PARAMETER_BACKG_HNO3_UNIT.compare( "ppm" ) == 0 ) {
                /* Convert ppm to ppb */
                for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                    y[counter-1][i] = Input_Opt.PARAMETER_BACKG_HNO3[0] * 1.00E+03;
            } else {
                std::cout << " Unknown unit for variable 'BACKG_HNO3': ";
                std::cout << Input_Opt.PARAMETER_BACKG_HNO3_UNIT << std::endl;
                exit(1);
            }

            /* Updating unit now that conversion has been taken care of */
            Input_Opt.PARAMETER_BACKG_HNO3_UNIT = "ppb";
        } else {
            std::cout << " In CombVec:";
            std::cout << " PARAMETER_BACKG_HNO3 has the wrong shape: ";
            std::cout << Input_Opt.PARAMETER_BACKG_HNO3.size() << std::endl;
        }

        counter += 1;

        /* ======================================================================= */
        /* ---- BACKGROUND O3 MIXING RATIO --------------------------------------- */
        /* ---- Assumed units are: ppb                                             */
        /* ======================================================================= */

        y.push_back( Vector_1D( Input_Opt.PARAMETER_FILECASES ) );
        if ( Input_Opt.PARAMETER_BACKG_O3.size() == Input_Opt.PARAMETER_FILECASES ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_BACKG_O3[i] * 1.00E+00;
        } else if ( Input_Opt.PARAMETER_BACKG_O3.size() == 1 ) {
            if ( Input_Opt.PARAMETER_BACKG_O3_UNIT.compare( "ppb" ) == 0 ) {
                /* Do nothing. Default unit */
                for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                    y[counter-1][i] = Input_Opt.PARAMETER_BACKG_O3[0];
            } else if ( Input_Opt.PARAMETER_BACKG_O3_UNIT.compare( "ppt" ) == 0 ) {
                /* Convert ppt to ppb */
                for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                    y[counter-1][i] = Input_Opt.PARAMETER_BACKG_O3[0] * 1.00E-03;
            } else if ( Input_Opt.PARAMETER_BACKG_O3_UNIT.compare( "ppm" ) == 0 ) {
                /* Convert ppm to ppb */
                for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                    y[counter-1][i] = Input_Opt.PARAMETER_BACKG_O3[0] * 1.00E+03;
            } else {
                std::cout << " Unknown unit for variable 'BACKG_O3': ";
                std::cout << Input_Opt.PARAMETER_BACKG_O3_UNIT << std::endl;
                exit(1);
            }

            /* Updating unit now that conversion has been taken care of */
            Input_Opt.PARAMETER_BACKG_O3_UNIT = "ppb";
        } else {
            std::cout << " In CombVec:";
            std::cout << " PARAMETER_BACKG_O3 has the wrong shape: ";
            std::cout << Input_Opt.PARAMETER_BACKG_O3.size() << std::endl;
        }

        counter += 1;

        /* ======================================================================= */
        /* ---- BACKGROUND CO MIXING RATIO --------------------------------------- */
        /* ---- Assumed units are: ppb                                             */
        /* ======================================================================= */

        y.push_back( Vector_1D( Input_Opt.PARAMETER_FILECASES ) );
        if ( Input_Opt.PARAMETER_BACKG_CO.size() == Input_Opt.PARAMETER_FILECASES ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_BACKG_CO[i] * 1.00E+00;
        } else if ( Input_Opt.PARAMETER_BACKG_CO.size() == 1 ) {
            if ( Input_Opt.PARAMETER_BACKG_CO_UNIT.compare( "ppb" ) == 0 ) {
                /* Do nothing. Default unit */
                for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                    y[counter-1][i] = Input_Opt.PARAMETER_BACKG_CO[0];
            } else if ( Input_Opt.PARAMETER_BACKG_CO_UNIT.compare( "ppt" ) == 0 ) {
                /* Convert ppt to ppb */
                for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                    y[counter-1][i] = Input_Opt.PARAMETER_BACKG_CO[0] * 1.00E-03;
            } else if ( Input_Opt.PARAMETER_BACKG_CO_UNIT.compare( "ppm" ) == 0 ) {
                /* Convert ppm to ppb */
                for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                    y[counter-1][i] = Input_Opt.PARAMETER_BACKG_CO[0] * 1.00E+03;
            } else {
                std::cout << " Unknown unit for variable 'BACKG_CO': ";
                std::cout << Input_Opt.PARAMETER_BACKG_CO_UNIT << std::endl;
                exit(1);
            }

            /* Updating unit now that conversion has been taken care of */
            Input_Opt.PARAMETER_BACKG_CO_UNIT = "ppb";
        } else {
            std::cout << " In CombVec:";
            std::cout << " PARAMETER_BACKG_CO has the wrong shape: ";
            std::cout << Input_Opt.PARAMETER_BACKG_CO.size() << std::endl;
        }

        counter += 1;

        /* ======================================================================= */
        /* ---- BACKGROUND CH4 MIXING RATIO -------------------------------------- */
        /* ---- Assumed units are: ppm                                             */
        /* ======================================================================= */

        y.push_back( Vector_1D( Input_Opt.PARAMETER_FILECASES ) );
        if ( Input_Opt.PARAMETER_BACKG_CH4.size() == Input_Opt.PARAMETER_FILECASES ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_BACKG_CH4[i] * 1.00E+03;
        } else if ( Input_Opt.PARAMETER_BACKG_CH4.size() == 1 ) {
            if ( Input_Opt.PARAMETER_BACKG_CH4_UNIT.compare( "ppb" ) == 0 ) {
                /* Do nothing. Default unit */
                for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                    y[counter-1][i] = Input_Opt.PARAMETER_BACKG_CH4[0];
            } else if ( Input_Opt.PARAMETER_BACKG_CH4_UNIT.compare( "ppt" ) == 0 ) {
                /* Convert ppt to ppb */
                for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                    y[counter-1][i] = Input_Opt.PARAMETER_BACKG_CH4[0] * 1.00E-03;
            } else if ( Input_Opt.PARAMETER_BACKG_CH4_UNIT.compare( "ppm" ) == 0 ) {
                /* Convert ppm to ppb */
                for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                    y[counter-1][i] = Input_Opt.PARAMETER_BACKG_CH4[0] * 1.00E+03;
            } else {
                std::cout << " Unknown unit for variable 'BACKG_CH4': ";
                std::cout << Input_Opt.PARAMETER_BACKG_CH4_UNIT << std::endl;
                exit(1);
            }

            /* Updating unit now that conversion has been taken care of */
            Input_Opt.PARAMETER_BACKG_CH4_UNIT = "ppb";
        } else {
            std::cout << " In CombVec:";
            std::cout << " PARAMETER_BACKG_CH4 has the wrong shape: ";
            std::cout << Input_Opt.PARAMETER_BACKG_CH4.size() << std::endl;
        }

        counter += 1;

        /* ======================================================================= */
        /* ---- BACKGROUND SO2 MIXING RATIO -------------------------------------- */
        /* ---- Assumed units are: ppt                                             */
        /* ======================================================================= */

        y.push_back( Vector_1D( Input_Opt.PARAMETER_FILECASES ) );
        if ( Input_Opt.PARAMETER_BACKG_SO2.size() == Input_Opt.PARAMETER_FILECASES ) {
            for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                y[counter-1][i] = Input_Opt.PARAMETER_BACKG_SO2[i] * 1.00E-03;
        } else if ( Input_Opt.PARAMETER_BACKG_SO2.size() == 1 ) {
            if ( Input_Opt.PARAMETER_BACKG_SO2_UNIT.compare( "ppb" ) == 0 ) {
                /* Do nothing. Default unit */
                for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                    y[counter-1][i] = Input_Opt.PARAMETER_BACKG_SO2[0];
            } else if ( Input_Opt.PARAMETER_BACKG_SO2_UNIT.compare( "ppt" ) == 0 ) {
                /* Convert ppt to ppb */
                for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                    y[counter-1][i] = Input_Opt.PARAMETER_BACKG_SO2[0] * 1.00E-03;
            } else if ( Input_Opt.PARAMETER_BACKG_SO2_UNIT.compare( "ppm" ) == 0 ) {
                /* Convert ppm to ppb */
                for ( i = 0; i < Input_Opt.PARAMETER_FILECASES; i++ )
                    y[counter-1][i] = Input_Opt.PARAMETER_BACKG_SO2[0] * 1.00E+03;
            } else {
                std::cout << " Unknown unit for variable 'BACKG_SO2': ";
                std::cout << Input_Opt.PARAMETER_BACKG_SO2_UNIT << std::endl;
                exit(1);
            }

            /* Updating unit now that conversion has been taken care of */
            Input_Opt.PARAMETER_BACKG_SO2_UNIT = "ppb";
        } else {
            std::cout << " In CombVec:";
            std::cout << " PARAMETER_BACKG_SO2 has the wrong shape: ";
            std::cout << Input_Opt.PARAMETER_BACKG_SO2.size() << std::endl;
        }

        counter += 1;

        if ( VERBOSE ) {
            for ( i = 0; i < y.size(); i++ ) {
                for ( j = 0; j < y[i].size(); j++ )
                    std::cout << y[i][j] << ", ";
                std::cout << std::endl;
            }
        }

        return y;

    } else {
        std::cout << "combvec" << std::endl;
        /* ======================================================================= */
        /* ---- PLUME PROCESSING TIME ( SIMULATION TIME ) ------------------------ */
        /* ---- Accepted units are: hr (default)                                   */
        /* ======================================================================= */
        
	if ( Input_Opt.PARAMETER_PLUMEPROCESS_RANGE ) {
            currVal = Input_Opt.PARAMETER_PLUMEPROCESS[0];
            while ( currVal <= Input_Opt.PARAMETER_PLUMEPROCESS[2] ) {
                cases.push_back( currVal );
                currVal += Input_Opt.PARAMETER_PLUMEPROCESS[1];
            }
        } else {
            for ( i = 0; i < Input_Opt.PARAMETER_PLUMEPROCESS.size(); i++ )
                cases.push_back(Input_Opt.PARAMETER_PLUMEPROCESS[i]);
        }
        nCases = cases.size();

        /* Do unit conversion to default unit */
        if ( ( Input_Opt.PARAMETER_PLUMEPROCESS_UNIT.compare( "0-24" ) == 0 ) || \
             ( Input_Opt.PARAMETER_PLUMEPROCESS_UNIT.compare( "hr" ) == 0 ) ) {
            /* Do nothing. Default unit */
        } else if ( Input_Opt.PARAMETER_PLUMEPROCESS_UNIT.compare( "min" ) == 0 ) {
            /* Convert from min to hr */
            for ( i = 0; i < cases.size(); i++ )
                cases[i] /= 60.0;
        } else {
            std::cout << " Unknown unit for variable 'Plume Processing Time': ";
            std::cout << Input_Opt.PARAMETER_PLUMEPROCESS_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_PLUMEPROCESS_UNIT = "hr";

        y.push_back(Vector_1D(cases.size()));
        for ( i = 0; i < cases.size(); i++ )
            y[counter-1][i] = cases[i];
        cases.clear();

        /* ======================================================================= */
        /* ---- TEMPERATURE ------------------------------------------------------ */
        /* ---- Accepted units are: Kelvin (default), Celsius, Fahrenheit          */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_TEMPERATURE_RANGE ) {
            currVal = Input_Opt.PARAMETER_TEMPERATURE[0];
            while ( currVal <= Input_Opt.PARAMETER_TEMPERATURE[2] ) {
                cases.push_back( currVal );
                currVal += Input_Opt.PARAMETER_TEMPERATURE[1];
            }
        } else {
            for ( i = 0; i < Input_Opt.PARAMETER_TEMPERATURE.size(); i++ )
                cases.push_back(Input_Opt.PARAMETER_TEMPERATURE[i]);
        }

        /* Do unit conversion to default unit */
        if ( Input_Opt.PARAMETER_TEMPERATURE_UNIT.compare( "K" ) == 0 ) {
            /* Do nothing. Default unit */
        } else if ( Input_Opt.PARAMETER_TEMPERATURE_UNIT.compare( "C" ) == 0 ) {
            /* Convert C to K */
            for ( i = 0; i < cases.size(); i++ )
                cases[i] += 273.15;
        } else if ( Input_Opt.PARAMETER_TEMPERATURE_UNIT.compare( "F" ) == 0 ) {
            /* Convert F to K */
            for ( i = 0; i < cases.size(); i++ )
                cases[i] = ( cases[i] + 459.67 ) * (double) 5.0/9;
        } else {
            std::cout << " Unknown unit for variable 'Temperature': ";
            std::cout << Input_Opt.PARAMETER_TEMPERATURE_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_TEMPERATURE_UNIT = "K";

	std::cout << cases.size() << std::endl;
        z.push_back( Vector_1D(cases.size() ) );
        for ( i = 0; i < cases.size(); i++ ) {
            z[0][i] = cases[i];
            std::cout << i << ", " << cases[i] << std::endl;
        }
        nCases *= cases.size();

        u = Copy_blocked(y,z[0].size());
        v = Copy_interleaved(z,y[0].size());

        for ( i = 0; i < counter; i++ )
            y[i].clear();
        z[0].clear();
        y.clear(); z.clear();

        counter += 1;
        for ( i = 0; i < counter; i++ )
            y.push_back(Vector_1D( nCases ));

        for ( i = 0; i < nCases; i++ ) {
            for ( j = 0; j < counter - 1; j++ )
                y[j][i] = u[j][i];
            y[counter-1][i] = v[0][i];
        }
        cases.clear();

        /* ======================================================================= */
        /* ---- RELATIVE HUMIDITY ------------------------------------------------ */
        /* ---- Accepted units are: % (0-100, default), - (0-1)                    */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_RHW_RANGE ) {
            currVal = Input_Opt.PARAMETER_RHW[0];
            while ( currVal <= Input_Opt.PARAMETER_RHW[2] ) {
                cases.push_back( currVal );
                currVal += Input_Opt.PARAMETER_RHW[1];
            }
        } else {
            for ( i = 0; i < Input_Opt.PARAMETER_RHW.size(); i++ )
                cases.push_back(Input_Opt.PARAMETER_RHW[i]);
        }

        if ( Input_Opt.PARAMETER_RHW_UNIT.compare( "%" ) == 0 ) {
            /* Do nothing. Default unit */
        } else if ( Input_Opt.PARAMETER_RHW_UNIT.compare( "-" ) == 0 ) {
            /* Convert - to % */
            for ( i = 0; i < cases.size(); i++ )
                cases[i] *= 100.0;
        } else {
            std::cout << " Unknown unit for variable 'Relative humidity': ";
            std::cout << Input_Opt.PARAMETER_RHW_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_RHW_UNIT = "%";

        z.push_back( Vector_1D(cases.size() ) );
        for ( i = 0; i < cases.size(); i++ )
            z[0][i] = cases[i];
        nCases *= cases.size();

        u = Copy_blocked(y,z[0].size());
        v = Copy_interleaved(z,y[0].size());

        for ( i = 0; i < counter; i++ )
            y[i].clear();
        z[0].clear();
        y.clear(); z.clear();

        counter += 1;
        for ( i = 0; i < counter; i++ )
            y.push_back(Vector_1D( nCases ));

        for ( i = 0; i < nCases; i++ ) {
            for ( j = 0; j < counter - 1; j++ )
                y[j][i] = u[j][i];
            y[counter-1][i] = v[0][i];
        }
        cases.clear();

        /* ======================================================================= */
        /* ---- HORIZONTAL DIFFUSION PARAMETER ----------------------------------- */
        /* ---- Accepted units are: m^2/s (default)                                */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_DH_RANGE ) {
            currVal = Input_Opt.PARAMETER_DH[0];
            while ( currVal <= Input_Opt.PARAMETER_DH[2] ) {
                cases.push_back( currVal );
                currVal += Input_Opt.PARAMETER_DH[1];
            }
        } else {
            for ( i = 0; i < Input_Opt.PARAMETER_DH.size(); i++ )
                cases.push_back(Input_Opt.PARAMETER_DH[i]);
        }

        if ( Input_Opt.PARAMETER_DH_UNIT.compare( "m^2/s" ) == 0 ) {
            /* Do nothing. Default unit */
        } else {
            std::cout << " Unknown unit for variable 'Horiz. diffusion parameter': ";
            std::cout << Input_Opt.PARAMETER_DH_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_DH_UNIT = "m^2/s";

        z.push_back( Vector_1D(cases.size() ) );
        for ( i = 0; i < cases.size(); i++ )
            z[0][i] = cases[i];
        nCases *= cases.size();

        u = Copy_blocked(y,z[0].size());
        v = Copy_interleaved(z,y[0].size());

        for ( i = 0; i < counter; i++ )
            y[i].clear();
        z[0].clear();
        y.clear(); z.clear();

        counter += 1;
        for ( i = 0; i < counter; i++ )
            y.push_back(Vector_1D( nCases ));

        for ( i = 0; i < nCases; i++ ) {
            for ( j = 0; j < counter - 1; j++ )
                y[j][i] = u[j][i];
            y[counter-1][i] = v[0][i];
        }
        cases.clear();

        /* ======================================================================= */
        /* ---- VERTICAL DIFFUSION PARAMETER ------------------------------------- */
        /* ---- Accepted units are: m^2/s (default)                                */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_DV_RANGE ) {
            currVal = Input_Opt.PARAMETER_DV[0];
            while ( currVal <= Input_Opt.PARAMETER_DV[2] ) {
                cases.push_back( currVal );
                currVal += Input_Opt.PARAMETER_DV[1];
            }
        } else {
            for ( i = 0; i < Input_Opt.PARAMETER_DV.size(); i++ )
                cases.push_back(Input_Opt.PARAMETER_DV[i]);
        }

        if ( Input_Opt.PARAMETER_DV_UNIT.compare( "m^2/s" ) == 0 ) {
            /* Do nothing. Default unit */
        } else {
            std::cout << " Unknown unit for variable 'Vert. diffusion parameter': ";
            std::cout << Input_Opt.PARAMETER_DV_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_DV_UNIT = "m^2/s";

        z.push_back( Vector_1D(cases.size() ) );
        for ( i = 0; i < cases.size(); i++ )
            z[0][i] = cases[i];
        nCases *= cases.size();

        u = Copy_blocked(y,z[0].size());
        v = Copy_interleaved(z,y[0].size());

        for ( i = 0; i < counter; i++ )
            y[i].clear();
        z[0].clear();
        y.clear(); z.clear();

        counter += 1;
        for ( i = 0; i < counter; i++ )
            y.push_back(Vector_1D( nCases ));

        for ( i = 0; i < nCases; i++ ) {
            for ( j = 0; j < counter - 1; j++ )
                y[j][i] = u[j][i];
            y[counter-1][i] = v[0][i];
        }
        cases.clear();

        /* ======================================================================= */
        /* ---- SHEAR ------------------------------------------------------------ */
        /* ---- Accepted units are: 1/s (default)                                  */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_SHEAR_RANGE ) {
            currVal = Input_Opt.PARAMETER_SHEAR[0];
            while ( currVal <= Input_Opt.PARAMETER_SHEAR[2] ) {
                cases.push_back( currVal );
                currVal += Input_Opt.PARAMETER_SHEAR[1];
            }
        } else {
            for ( i = 0; i < Input_Opt.PARAMETER_SHEAR.size(); i++ )
                cases.push_back(Input_Opt.PARAMETER_SHEAR[i]);
        }

        if ( Input_Opt.PARAMETER_SHEAR_UNIT.compare( "1/s" ) == 0 ) {
            /* Do nothing. Default unit */
        } else {
            std::cout << " Unknown unit for variable 'Shear': ";
            std::cout << Input_Opt.PARAMETER_SHEAR_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_SHEAR_UNIT = "1/s";

        z.push_back( Vector_1D(cases.size() ) );
        for ( i = 0; i < cases.size(); i++ )
            z[0][i] = cases[i];
        nCases *= cases.size();

        u = Copy_blocked(y,z[0].size());
        v = Copy_interleaved(z,y[0].size());

        for ( i = 0; i < counter; i++ )
            y[i].clear();
        z[0].clear();
        y.clear(); z.clear();

        counter += 1;
        for ( i = 0; i < counter; i++ )
            y.push_back(Vector_1D( nCases ));

        for ( i = 0; i < nCases; i++ ) {
            for ( j = 0; j < counter - 1; j++ )
                y[j][i] = u[j][i];
            y[counter-1][i] = v[0][i];
        }
        cases.clear();

        /* ======================================================================= */
        /* ---- LONGITUDE -------------------------------------------------------- */
        /* ---- Accepted units are: degree (default)                               */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_LONGITUDE_RANGE ) {
            currVal = Input_Opt.PARAMETER_LONGITUDE[0];
            while ( currVal <= Input_Opt.PARAMETER_LONGITUDE[2] ) {
                cases.push_back( currVal );
                currVal += Input_Opt.PARAMETER_LONGITUDE[1];
            }
        } else {
            for ( i = 0; i < Input_Opt.PARAMETER_LONGITUDE.size(); i++ )
                cases.push_back(Input_Opt.PARAMETER_LONGITUDE[i]);
        }

        if ( Input_Opt.PARAMETER_LONGITUDE_UNIT.compare( "deg" ) == 0 ) {
            /* Do nothing. Default unit */
        } else {
            std::cout << " Unknown unit for variable 'Longitude': ";
            std::cout << Input_Opt.PARAMETER_LONGITUDE_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_LONGITUDE_UNIT = "deg";

        z.push_back( Vector_1D(cases.size() ) );
        for ( i = 0; i < cases.size(); i++ )
            z[0][i] = cases[i];
        nCases *= cases.size();

        u = Copy_blocked(y,z[0].size());
        v = Copy_interleaved(z,y[0].size());

        for ( i = 0; i < counter; i++ )
            y[i].clear();
        z[0].clear();
        y.clear(); z.clear();

        counter += 1;
        for ( i = 0; i < counter; i++ )
            y.push_back(Vector_1D( nCases ));

        for ( i = 0; i < nCases; i++ ) {
            for ( j = 0; j < counter - 1; j++ )
                y[j][i] = u[j][i];
            y[counter-1][i] = v[0][i];
        }
        cases.clear();

        /* ======================================================================= */
        /* ---- LATITUDE --------------------------------------------------------- */
        /* ---- Accepted units are: degree (default)                               */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_LATITUDE_RANGE ) {
            currVal = Input_Opt.PARAMETER_LATITUDE[0];
            while ( currVal <= Input_Opt.PARAMETER_LATITUDE[2] ) {
                cases.push_back( currVal );
                currVal += Input_Opt.PARAMETER_LATITUDE[1];
            }
        } else {
            for ( i = 0; i < Input_Opt.PARAMETER_LATITUDE.size(); i++ )
                cases.push_back(Input_Opt.PARAMETER_LATITUDE[i]);
        }

        if ( Input_Opt.PARAMETER_LATITUDE_UNIT.compare( "deg" ) == 0 ) {
            /* Do nothing. Default unit */
        } else {
            std::cout << " Unknown unit for variable 'Latitude': ";
            std::cout << Input_Opt.PARAMETER_LATITUDE_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_LATITUDE_UNIT = "deg";

        z.push_back( Vector_1D(cases.size() ) );
        for ( i = 0; i < cases.size(); i++ )
            z[0][i] = cases[i];
        nCases *= cases.size();

        u = Copy_blocked(y,z[0].size());
        v = Copy_interleaved(z,y[0].size());

        for ( i = 0; i < counter; i++ )
            y[i].clear();
        z[0].clear();
        y.clear(); z.clear();

        counter += 1;
        for ( i = 0; i < counter; i++ )
            y.push_back(Vector_1D( nCases ));

        for ( i = 0; i < nCases; i++ ) {
            for ( j = 0; j < counter - 1; j++ )
                y[j][i] = u[j][i];
            y[counter-1][i] = v[0][i];
        }
        cases.clear();

        /* ======================================================================= */
        /* ---- PRESSURE --------------------------------------------------------- */
        /* ---- Accepted units are: Pa (default), hPa                              */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_PRESSURE_RANGE ) {
            currVal = Input_Opt.PARAMETER_PRESSURE[0];
            while ( currVal <= Input_Opt.PARAMETER_PRESSURE[2] ) {
                cases.push_back( currVal );
                currVal += Input_Opt.PARAMETER_PRESSURE[1];
            }
        } else {
            for ( i = 0; i < Input_Opt.PARAMETER_PRESSURE.size(); i++ )
                cases.push_back(Input_Opt.PARAMETER_PRESSURE[i]);
        }

        if ( Input_Opt.PARAMETER_PRESSURE_UNIT.compare( "Pa" ) == 0 ) {
            /* Do nothing. Default unit */
        } else if ( Input_Opt.PARAMETER_PRESSURE_UNIT.compare( "hPa" ) == 0 ) {
            for ( i = 0; i < cases.size(); i++ )
                cases[i] *= 100.0;
        } else {
            std::cout << " Unknown unit for variable 'Pressure': ";
            std::cout << Input_Opt.PARAMETER_PRESSURE_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_PRESSURE_UNIT = "Pa";

        z.push_back( Vector_1D(cases.size() ) );
        for ( i = 0; i < cases.size(); i++ )
            z[0][i] = cases[i];
        nCases *= cases.size();

        u = Copy_blocked(y,z[0].size());
        v = Copy_interleaved(z,y[0].size());

        for ( i = 0; i < counter; i++ )
            y[i].clear();
        z[0].clear();
        y.clear(); z.clear();

        counter += 1;
        for ( i = 0; i < counter; i++ )
            y.push_back(Vector_1D( nCases ));

        for ( i = 0; i < nCases; i++ ) {
            for ( j = 0; j < counter - 1; j++ )
                y[j][i] = u[j][i];
            y[counter-1][i] = v[0][i];
        }
        cases.clear();

        /* ======================================================================= */
        /* ---- EMISSION DAY ----------------------------------------------------- */
        /* ---- Accepted units are: 1-365 (default)                                */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_EDAY_RANGE ) {
            currVal = Input_Opt.PARAMETER_EDAY[0];
            while ( currVal <= Input_Opt.PARAMETER_EDAY[2] ) {
                cases.push_back( currVal );
                currVal += Input_Opt.PARAMETER_EDAY[1];
            }
        } else {
            for ( i = 0; i < Input_Opt.PARAMETER_EDAY.size(); i++ )
                cases.push_back(Input_Opt.PARAMETER_EDAY[i]);
        }

        if ( ( Input_Opt.PARAMETER_EDAY_UNIT.compare( "1-365" ) == 0 ) || \
             ( Input_Opt.PARAMETER_EDAY_UNIT.compare( "-" ) == 0 ) ) {
            /* Do nothing. Default unit */
        } else {
            std::cout << " Unknown unit for variable 'Emission Day': ";
            std::cout << Input_Opt.PARAMETER_EDAY_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_EDAY_UNIT = "1-365";

        z.push_back( Vector_1D(cases.size() ) );
        for ( i = 0; i < cases.size(); i++ )
            z[0][i] = cases[i];
        nCases *= cases.size();

        u = Copy_blocked(y,z[0].size());
        v = Copy_interleaved(z,y[0].size());

        for ( i = 0; i < counter; i++ )
            y[i].clear();
        z[0].clear();
        y.clear(); z.clear();

        counter += 1;
        for ( i = 0; i < counter; i++ )
            y.push_back(Vector_1D( nCases ));

        for ( i = 0; i < nCases; i++ ) {
            for ( j = 0; j < counter - 1; j++ )
                y[j][i] = u[j][i];
            y[counter-1][i] = v[0][i];
        }
        cases.clear();

        /* ======================================================================= */
        /* ---- EMISSION TIME ---------------------------------------------------- */
        /* ---- Accepted units are: 0-24 (default)                                 */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_ETIME_RANGE ) {
            currVal = std::fmod( Input_Opt.PARAMETER_ETIME[0], 24.0 );
            while ( currVal <= Input_Opt.PARAMETER_ETIME[2] ) {
                cases.push_back( currVal );
                currVal += std::fmod( Input_Opt.PARAMETER_ETIME[1], 24.0 );
            }
        } else {
            for ( i = 0; i < Input_Opt.PARAMETER_ETIME.size(); i++ )
                cases.push_back( std::fmod( Input_Opt.PARAMETER_ETIME[i], 24.0 ) );
        }

        if ( ( Input_Opt.PARAMETER_ETIME_UNIT.compare( "0-24" ) == 0 ) || \
             ( Input_Opt.PARAMETER_ETIME_UNIT.compare( "hr" ) == 0 ) ) {
            /* Do nothing. Default unit */
        } else {
            std::cout << " Unknown unit for variable 'Emission Time': ";
            std::cout << Input_Opt.PARAMETER_ETIME_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_ETIME_UNIT = "0-24";

        z.push_back( Vector_1D(cases.size() ) );
        for ( i = 0; i < cases.size(); i++ )
            z[0][i] = cases[i];
        nCases *= cases.size();

        u = Copy_blocked(y,z[0].size());
        v = Copy_interleaved(z,y[0].size());

        for ( i = 0; i < counter; i++ )
            y[i].clear();
        z[0].clear();
        y.clear(); z.clear();

        counter += 1;
        for ( i = 0; i < counter; i++ )
            y.push_back(Vector_1D( nCases ));

        for ( i = 0; i < nCases; i++ ) {
            for ( j = 0; j < counter - 1; j++ )
                y[j][i] = u[j][i];
            y[counter-1][i] = v[0][i];
        }
        cases.clear();

        /* ======================================================================= */
        /* ---- EI_NOX ----------------------------------------------------------- */
        /* ---- Accepted units are: g/kg_fuel (default)                            */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_EI_NOX_RANGE ) {
            currVal = Input_Opt.PARAMETER_EI_NOX[0];
            while ( currVal <= Input_Opt.PARAMETER_EI_NOX[2] ) {
                cases.push_back( currVal );
                currVal += Input_Opt.PARAMETER_EI_NOX[1];
            }
        } else {
            for ( i = 0; i < Input_Opt.PARAMETER_EI_NOX.size(); i++ )
                cases.push_back(Input_Opt.PARAMETER_EI_NOX[i]);
        }

        if ( ( Input_Opt.PARAMETER_EI_NOX_UNIT.compare( "g/kg_fuel" )      == 0 ) || \
             ( Input_Opt.PARAMETER_EI_NOX_UNIT.compare( "g/kg_f" )         == 0 ) || \
             ( Input_Opt.PARAMETER_EI_NOX_UNIT.compare( "g/kg" )           == 0 ) || \
             ( Input_Opt.PARAMETER_EI_NOX_UNIT.compare( "gNO2/kg_fuel" )   == 0 ) || \
             ( Input_Opt.PARAMETER_EI_NOX_UNIT.compare( "gNO2/kg_f" )      == 0 ) || \
             ( Input_Opt.PARAMETER_EI_NOX_UNIT.compare( "gNO2/kg" )        == 0 ) || \
             ( Input_Opt.PARAMETER_EI_NOX_UNIT.compare( "g(NO2)/kg_fuel" ) == 0 ) || \
             ( Input_Opt.PARAMETER_EI_NOX_UNIT.compare( "g(NO2)/kg_f" )    == 0 ) || \
             ( Input_Opt.PARAMETER_EI_NOX_UNIT.compare( "g(NO2)/kg" )      == 0 ) ) {
            /* Do nothing. Default unit */
        } else {
            std::cout << " Unknown unit for variable 'EI_NOX': ";
            std::cout << Input_Opt.PARAMETER_EI_NOX_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_EI_NOX_UNIT = "g/kg_fuel";

        z.push_back( Vector_1D(cases.size() ) );
        for ( i = 0; i < cases.size(); i++ )
            z[0][i] = cases[i];
        nCases *= cases.size();

        u = Copy_blocked(y,z[0].size());
        v = Copy_interleaved(z,y[0].size());

        for ( i = 0; i < counter; i++ )
            y[i].clear();
        z[0].clear();
        y.clear(); z.clear();

        counter += 1;
        for ( i = 0; i < counter; i++ )
            y.push_back(Vector_1D( nCases ));

        for ( i = 0; i < nCases; i++ ) {
            for ( j = 0; j < counter - 1; j++ )
                y[j][i] = u[j][i];
            y[counter-1][i] = v[0][i];
        }
        cases.clear();

        /* ======================================================================= */
        /* ---- EI_CO ------------------------------------------------------------ */
        /* ---- Accepted units are: g/kg_fuel (default)                            */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_EI_CO_RANGE ) {
            currVal = Input_Opt.PARAMETER_EI_CO[0];
            while ( currVal <= Input_Opt.PARAMETER_EI_CO[2] ) {
                cases.push_back( currVal );
                currVal += Input_Opt.PARAMETER_EI_CO[1];
            }
        } else {
            for ( i = 0; i < Input_Opt.PARAMETER_EI_CO.size(); i++ )
                cases.push_back(Input_Opt.PARAMETER_EI_CO[i]);
        }

        if ( ( Input_Opt.PARAMETER_EI_CO_UNIT.compare( "g/kg_fuel" ) == 0 ) || \
             ( Input_Opt.PARAMETER_EI_CO_UNIT.compare( "g/kg_f" )    == 0 ) || \
             ( Input_Opt.PARAMETER_EI_CO_UNIT.compare( "g/kg" )      == 0 ) ) {
            /* Do nothing. Default unit */
        } else {
            std::cout << " Unknown unit for variable 'EI_CO': ";
            std::cout << Input_Opt.PARAMETER_EI_CO_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_EI_CO_UNIT = "g/kg_fuel";

        z.push_back( Vector_1D(cases.size() ) );
        for ( i = 0; i < cases.size(); i++ )
            z[0][i] = cases[i];
        nCases *= cases.size();

        u = Copy_blocked(y,z[0].size());
        v = Copy_interleaved(z,y[0].size());

        for ( i = 0; i < counter; i++ )
            y[i].clear();
        z[0].clear();
        y.clear(); z.clear();

        counter += 1;
        for ( i = 0; i < counter; i++ )
            y.push_back(Vector_1D( nCases ));

        for ( i = 0; i < nCases; i++ ) {
            for ( j = 0; j < counter - 1; j++ )
                y[j][i] = u[j][i];
            y[counter-1][i] = v[0][i];
        }
        cases.clear();

        /* ======================================================================= */
        /* ---- EI_UHC ----------------------------------------------------------- */
        /* ---- Accepted units are: g/kg_fuel (default)                            */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_EI_UHC_RANGE ) {
            currVal = Input_Opt.PARAMETER_EI_UHC[0];
            while ( currVal <= Input_Opt.PARAMETER_EI_UHC[2] ) {
                cases.push_back( currVal );
                currVal += Input_Opt.PARAMETER_EI_UHC[1];
            }
        } else {
            for ( i = 0; i < Input_Opt.PARAMETER_EI_UHC.size(); i++ )
                cases.push_back(Input_Opt.PARAMETER_EI_UHC[i]);
        }

        if ( ( Input_Opt.PARAMETER_EI_UHC_UNIT.compare( "g/kg_fuel" ) == 0 ) || \
             ( Input_Opt.PARAMETER_EI_UHC_UNIT.compare( "g/kg_f" )    == 0 ) || \
             ( Input_Opt.PARAMETER_EI_UHC_UNIT.compare( "g/kg" )      == 0 ) ) {
            /* Do nothing. Default unit */
        } else {
            std::cout << " Unknown unit for variable 'EI_UHC': ";
            std::cout << Input_Opt.PARAMETER_EI_UHC_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_EI_UHC_UNIT = "g/kg_fuel";

        z.push_back( Vector_1D(cases.size() ) );
        for ( i = 0; i < cases.size(); i++ )
            z[0][i] = cases[i];
        nCases *= cases.size();

        u = Copy_blocked(y,z[0].size());
        v = Copy_interleaved(z,y[0].size());

        for ( i = 0; i < counter; i++ )
            y[i].clear();
        z[0].clear();
        y.clear(); z.clear();

        counter += 1;
        for ( i = 0; i < counter; i++ )
            y.push_back(Vector_1D( nCases ));

        for ( i = 0; i < nCases; i++ ) {
            for ( j = 0; j < counter - 1; j++ )
                y[j][i] = u[j][i];
            y[counter-1][i] = v[0][i];
        }
        cases.clear();

        /* ======================================================================= */
        /* ---- EI_SO2 ----------------------------------------------------------- */
        /* ---- Accepted units are: g/kg_fuel (default)                            */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_EI_SO2_RANGE ) {
            currVal = Input_Opt.PARAMETER_EI_SO2[0];
            while ( currVal <= Input_Opt.PARAMETER_EI_SO2[2] ) {
                cases.push_back( currVal );
                currVal += Input_Opt.PARAMETER_EI_SO2[1];
            }
        } else {
            for ( i = 0; i < Input_Opt.PARAMETER_EI_SO2.size(); i++ )
                cases.push_back(Input_Opt.PARAMETER_EI_SO2[i]);
        }

        if ( ( Input_Opt.PARAMETER_EI_SO2_UNIT.compare( "g/kg_fuel" ) == 0 ) || \
             ( Input_Opt.PARAMETER_EI_SO2_UNIT.compare( "g/kg_f" )    == 0 ) || \
             ( Input_Opt.PARAMETER_EI_SO2_UNIT.compare( "g/kg" )      == 0 ) ) {
            /* Do nothing. Default unit */
        } else {
            std::cout << " Unknown unit for variable 'EI_SO2': ";
            std::cout << Input_Opt.PARAMETER_EI_SO2_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_EI_SO2_UNIT = "g/kg_fuel";

        z.push_back( Vector_1D(cases.size() ) );
        for ( i = 0; i < cases.size(); i++ )
            z[0][i] = cases[i];
        nCases *= cases.size();

        u = Copy_blocked(y,z[0].size());
        v = Copy_interleaved(z,y[0].size());

        for ( i = 0; i < counter; i++ )
            y[i].clear();
        z[0].clear();
        y.clear(); z.clear();

        counter += 1;
        for ( i = 0; i < counter; i++ )
            y.push_back(Vector_1D( nCases ));

        for ( i = 0; i < nCases; i++ ) {
            for ( j = 0; j < counter - 1; j++ )
                y[j][i] = u[j][i];
            y[counter-1][i] = v[0][i];
        }
        cases.clear();

        /* ======================================================================= */
        /* ---- EI_SO2TOSO4 ------------------------------------------------------ */
        /* ---- Accepted units are: - (default), %                                 */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_EI_SO2TOSO4_RANGE ) {
            currVal = Input_Opt.PARAMETER_EI_SO2TOSO4[0];
            while ( currVal <= Input_Opt.PARAMETER_EI_SO2TOSO4[2] ) {
                cases.push_back( currVal );
                currVal += Input_Opt.PARAMETER_EI_SO2TOSO4[1];
            }
        } else {
            for ( i = 0; i < Input_Opt.PARAMETER_EI_SO2TOSO4.size(); i++ )
                cases.push_back(Input_Opt.PARAMETER_EI_SO2TOSO4[i]);
        }

        if ( Input_Opt.PARAMETER_EI_SO2TOSO4_UNIT.compare( "-" ) == 0 ) {
            /* Do nothing. Default unit */
        } else if ( Input_Opt.PARAMETER_EI_SO2TOSO4_UNIT.compare( "%" ) == 0 ) {
            /* Convert % to - */
            for ( i = 0; i < cases.size(); i++ )
                cases[i] /= 100.0;
        } else {
            std::cout << " Unknown unit for variable 'EI_SO2TOSO4': ";
            std::cout << Input_Opt.PARAMETER_EI_SO2TOSO4_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_EI_SO2TOSO4_UNIT = "-";

        z.push_back( Vector_1D(cases.size() ) );
        for ( i = 0; i < cases.size(); i++ )
            z[0][i] = cases[i];
        nCases *= cases.size();

        u = Copy_blocked(y,z[0].size());
        v = Copy_interleaved(z,y[0].size());

        for ( i = 0; i < counter; i++ )
            y[i].clear();
        z[0].clear();
        y.clear(); z.clear();

        counter += 1;
        for ( i = 0; i < counter; i++ )
            y.push_back(Vector_1D( nCases ));

        for ( i = 0; i < nCases; i++ ) {
            for ( j = 0; j < counter - 1; j++ )
                y[j][i] = u[j][i];
            y[counter-1][i] = v[0][i];
        }
        cases.clear();

        /* ======================================================================= */
        /* ---- EI_SOOT ---------------------------------------------------------- */
        /* ---- Accepted units are: g/kg_fuel (default)                            */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_EI_SOOT_RANGE ) {
            currVal = Input_Opt.PARAMETER_EI_SOOT[0];
            while ( currVal <= Input_Opt.PARAMETER_EI_SOOT[2] ) {
                cases.push_back( currVal );
                currVal += Input_Opt.PARAMETER_EI_SOOT[1];
            }
        } else {
            for ( i = 0; i < Input_Opt.PARAMETER_EI_SOOT.size(); i++ )
                cases.push_back(Input_Opt.PARAMETER_EI_SOOT[i]);
        }

        if ( ( Input_Opt.PARAMETER_EI_SOOT_UNIT.compare( "g/kg_fuel" ) == 0 ) || \
             ( Input_Opt.PARAMETER_EI_SOOT_UNIT.compare( "g/kg_f" )    == 0 ) || \
             ( Input_Opt.PARAMETER_EI_SOOT_UNIT.compare( "g/kg" )      == 0 ) ) {
            /* Do nothing. Default unit */
        } else if ( ( Input_Opt.PARAMETER_EI_SOOT_UNIT.compare( "mg/kg_fuel" ) == 0 ) || \
             ( Input_Opt.PARAMETER_EI_SOOT_UNIT.compare( "mg/kg_f" )    == 0 ) || \
             ( Input_Opt.PARAMETER_EI_SOOT_UNIT.compare( "mg/kg" )      == 0 ) ) {
            /* Convert mg to g */
            for ( i = 0; i < cases.size(); i++ )
                cases[i] /= 1.00E+03;
        } else {
            std::cout << " Unknown unit for variable 'EI_SOOT': ";
            std::cout << Input_Opt.PARAMETER_EI_SOOT_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_EI_SOOT_UNIT = "g/kg_fuel";

        z.push_back( Vector_1D(cases.size() ) );
        for ( i = 0; i < cases.size(); i++ )
            z[0][i] = cases[i];
        nCases *= cases.size();

        u = Copy_blocked(y,z[0].size());
        v = Copy_interleaved(z,y[0].size());

        for ( i = 0; i < counter; i++ )
            y[i].clear();
        z[0].clear();
        y.clear(); z.clear();

        counter += 1;
        for ( i = 0; i < counter; i++ )
            y.push_back(Vector_1D( nCases ));

        for ( i = 0; i < nCases; i++ ) {
            for ( j = 0; j < counter - 1; j++ )
                y[j][i] = u[j][i];
            y[counter-1][i] = v[0][i];
        }
        cases.clear();

        /* ======================================================================= */
        /* ---- EI_SOOTRAD ------------------------------------------------------- */
        /* ---- Accepted units are: m (default), nm                                */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_EI_SOOTRAD_RANGE ) {
            currVal = Input_Opt.PARAMETER_EI_SOOTRAD[0];
            while ( currVal <= Input_Opt.PARAMETER_EI_SOOTRAD[2] ) {
                cases.push_back( currVal );
                currVal += Input_Opt.PARAMETER_EI_SOOTRAD[1];
            }
        } else {
            for ( i = 0; i < Input_Opt.PARAMETER_EI_SOOTRAD.size(); i++ )
                cases.push_back(Input_Opt.PARAMETER_EI_SOOTRAD[i]);
        }

        if ( Input_Opt.PARAMETER_EI_SOOTRAD_UNIT.compare( "m" ) == 0 ) {
            /* Do nothing. Default unit */
        } else if ( Input_Opt.PARAMETER_EI_SOOTRAD_UNIT.compare( "nm" ) == 0 ) {
            /* Convert nm to m */
            for ( i = 0; i < cases.size(); i++ )
                cases[i] /= 1.00E+09;
        } else {
            std::cout << " Unknown unit for variable 'EI_SOOTRAD': ";
            std::cout << Input_Opt.PARAMETER_EI_SOOTRAD_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_EI_SOOTRAD_UNIT = "m";

        z.push_back( Vector_1D(cases.size() ) );
        for ( i = 0; i < cases.size(); i++ )
            z[0][i] = cases[i];
        nCases *= cases.size();

        u = Copy_blocked(y,z[0].size());
        v = Copy_interleaved(z,y[0].size());

        for ( i = 0; i < counter; i++ )
            y[i].clear();
        z[0].clear();
        y.clear(); z.clear();

        counter += 1;
        for ( i = 0; i < counter; i++ )
            y.push_back(Vector_1D( nCases ));

        for ( i = 0; i < nCases; i++ ) {
            for ( j = 0; j < counter - 1; j++ )
                y[j][i] = u[j][i];
            y[counter-1][i] = v[0][i];
        }
        cases.clear();

        /* ======================================================================= */
        /* ---- TOTAL FUEL FLOW -------------------------------------------------- */
        /* ---- Accepted units are: kg/s (default)                                 */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_FF_RANGE ) {
            currVal = Input_Opt.PARAMETER_FF[0];
            while ( currVal <= Input_Opt.PARAMETER_FF[2] ) {
                cases.push_back( currVal );
                currVal += Input_Opt.PARAMETER_FF[1];
            }
        } else {
            for ( i = 0; i < Input_Opt.PARAMETER_FF.size(); i++ )
                cases.push_back(Input_Opt.PARAMETER_FF[i]);
        }

        if ( ( Input_Opt.PARAMETER_FF_UNIT.compare( "kg_fuel/s" ) == 0 ) || \
             ( Input_Opt.PARAMETER_FF_UNIT.compare( "kg_f/s" )    == 0 ) || \
             ( Input_Opt.PARAMETER_FF_UNIT.compare( "kg/s" )      == 0 ) ) {
            /* Do nothing. Default unit */
        } else {
            std::cout << " Unknown unit for variable 'FF': ";
            std::cout << Input_Opt.PARAMETER_FF_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_FF_UNIT = "kg/s";

        z.push_back( Vector_1D(cases.size() ) );
        for ( i = 0; i < cases.size(); i++ )
            z[0][i] = cases[i];
        nCases *= cases.size();

        u = Copy_blocked(y,z[0].size());
        v = Copy_interleaved(z,y[0].size());

        for ( i = 0; i < counter; i++ )
            y[i].clear();
        z[0].clear();
        y.clear(); z.clear();

        counter += 1;
        for ( i = 0; i < counter; i++ )
            y.push_back(Vector_1D( nCases ));

        for ( i = 0; i < nCases; i++ ) {
            for ( j = 0; j < counter - 1; j++ )
                y[j][i] = u[j][i];
            y[counter-1][i] = v[0][i];
        }
        cases.clear();

        /* ======================================================================= */
        /* ---- AIRCRAFT MASS ---------------------------------------------------- */
        /* ---- Accepted units are: kg (default), tonnes                           */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_AMASS_RANGE ) {
            currVal = Input_Opt.PARAMETER_AMASS[0];
            while ( currVal <= Input_Opt.PARAMETER_AMASS[2] ) {
                cases.push_back( currVal );
                currVal += Input_Opt.PARAMETER_AMASS[1];
            }
        } else {
            for ( i = 0; i < Input_Opt.PARAMETER_AMASS.size(); i++ )
                cases.push_back(Input_Opt.PARAMETER_AMASS[i]);
        }

        if ( Input_Opt.PARAMETER_AMASS_UNIT.compare( "kg" ) == 0 ) {
            /* Do nothing. Default unit */
        } else if ( ( Input_Opt.PARAMETER_AMASS_UNIT.compare( "to" )      == 0 ) || \
                    ( Input_Opt.PARAMETER_AMASS_UNIT.compare( "tonne " )  == 0 ) || \
                    ( Input_Opt.PARAMETER_AMASS_UNIT.compare( "tonnes " ) == 0 ) ) {
            /* Convert tonnes to kg */
            for ( i = 0; i < cases.size(); i++ )
                cases[i] *= 1.00E+03;
        } else {
            std::cout << " Unknown unit for variable 'AMASS': ";
            std::cout << Input_Opt.PARAMETER_AMASS_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_AMASS_UNIT = "kg";

        z.push_back( Vector_1D(cases.size() ) );
        for ( i = 0; i < cases.size(); i++ )
            z[0][i] = cases[i];
        nCases *= cases.size();

        u = Copy_blocked(y,z[0].size());
        v = Copy_interleaved(z,y[0].size());

        for ( i = 0; i < counter; i++ )
            y[i].clear();
        z[0].clear();
        y.clear(); z.clear();

        counter += 1;
        for ( i = 0; i < counter; i++ )
            y.push_back(Vector_1D( nCases ));

        for ( i = 0; i < nCases; i++ ) {
            for ( j = 0; j < counter - 1; j++ )
                y[j][i] = u[j][i];
            y[counter-1][i] = v[0][i];
        }
        cases.clear();

        /* ======================================================================= */
        /* ---- BACKGROUND NOX MIXING RATIO -------------------------------------- */
        /* ---- Accepted units are: ppb (default), ppt, ppm                        */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_BACKG_NOX_RANGE ) {
            currVal = Input_Opt.PARAMETER_BACKG_NOX[0];
            while ( currVal <= Input_Opt.PARAMETER_BACKG_NOX[2] ) {
                cases.push_back( currVal );
                currVal += Input_Opt.PARAMETER_BACKG_NOX[1];
            }
        } else {
            for ( i = 0; i < Input_Opt.PARAMETER_BACKG_NOX.size(); i++ )
                cases.push_back(Input_Opt.PARAMETER_BACKG_NOX[i]);
        }

        if ( Input_Opt.PARAMETER_BACKG_NOX_UNIT.compare( "ppb" ) == 0 ) {
            /* Do nothing. Default unit */
        } else if ( Input_Opt.PARAMETER_BACKG_NOX_UNIT.compare( "ppt" ) == 0 ) {
            /* Convert ppt to ppb */
            for ( i = 0; i < cases.size(); i++ )
                cases[i] *= 1.00E-03;
        } else if ( Input_Opt.PARAMETER_BACKG_NOX_UNIT.compare( "ppm" ) == 0 ) {
            /* Convert ppm to ppb */
            for ( i = 0; i < cases.size(); i++ )
                cases[i] *= 1.00E+03;
        } else {
            std::cout << " Unknown unit for variable 'BACKG_NOX': ";
            std::cout << Input_Opt.PARAMETER_BACKG_NOX_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_BACKG_NOX_UNIT = "ppb";

        z.push_back( Vector_1D(cases.size() ) );
        for ( i = 0; i < cases.size(); i++ )
            z[0][i] = cases[i];
        nCases *= cases.size();

        u = Copy_blocked(y,z[0].size());
        v = Copy_interleaved(z,y[0].size());

        for ( i = 0; i < counter; i++ )
            y[i].clear();
        z[0].clear();
        y.clear(); z.clear();

        counter += 1;
        for ( i = 0; i < counter; i++ )
            y.push_back(Vector_1D( nCases ));

        for ( i = 0; i < nCases; i++ ) {
            for ( j = 0; j < counter - 1; j++ )
                y[j][i] = u[j][i];
            y[counter-1][i] = v[0][i];
        }
        cases.clear();

        /* ======================================================================= */
        /* ---- BACKGROUND HNO3 MIXING RATIO ------------------------------------- */
        /* ---- Accepted units are: ppb (default), ppt, ppm                        */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_BACKG_HNO3_RANGE ) {
            currVal = Input_Opt.PARAMETER_BACKG_HNO3[0];
            while ( currVal <= Input_Opt.PARAMETER_BACKG_HNO3[2] ) {
                cases.push_back( currVal );
                currVal += Input_Opt.PARAMETER_BACKG_HNO3[1];
            }
        } else {
            for ( i = 0; i < Input_Opt.PARAMETER_BACKG_HNO3.size(); i++ )
                cases.push_back(Input_Opt.PARAMETER_BACKG_HNO3[i]);
        }

        if ( Input_Opt.PARAMETER_BACKG_HNO3_UNIT.compare( "ppb" ) == 0 ) {
            /* Do nothing. Default unit */
        } else if ( Input_Opt.PARAMETER_BACKG_HNO3_UNIT.compare( "ppt" ) == 0 ) {
            /* Convert ppt to ppb */
            for ( i = 0; i < cases.size(); i++ )
                cases[i] *= 1.00E-03;
        } else if ( Input_Opt.PARAMETER_BACKG_HNO3_UNIT.compare( "ppm" ) == 0 ) {
            /* Convert ppm to ppb */
            for ( i = 0; i < cases.size(); i++ )
                cases[i] *= 1.00E+03;
        } else {
            std::cout << " Unknown unit for variable 'BACKG_HNO3': ";
            std::cout << Input_Opt.PARAMETER_BACKG_HNO3_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_BACKG_HNO3_UNIT = "HNO3";

        z.push_back( Vector_1D(cases.size() ) );
        for ( i = 0; i < cases.size(); i++ )
            z[0][i] = cases[i];
        nCases *= cases.size();

        u = Copy_blocked(y,z[0].size());
        v = Copy_interleaved(z,y[0].size());

        for ( i = 0; i < counter; i++ )
            y[i].clear();
        z[0].clear();
        y.clear(); z.clear();

        counter += 1;
        for ( i = 0; i < counter; i++ )
            y.push_back(Vector_1D( nCases ));

        for ( i = 0; i < nCases; i++ ) {
            for ( j = 0; j < counter - 1; j++ )
                y[j][i] = u[j][i];
            y[counter-1][i] = v[0][i];
        }
        cases.clear();

        /* ======================================================================= */
        /* ---- BACKGROUND O3 MIXING RATIO --------------------------------------- */
        /* ---- Accepted units are: ppb (default), ppt, ppm                        */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_BACKG_O3_RANGE ) {
            currVal = Input_Opt.PARAMETER_BACKG_O3[0];
            while ( currVal <= Input_Opt.PARAMETER_BACKG_O3[2] ) {
                cases.push_back( currVal );
                currVal += Input_Opt.PARAMETER_BACKG_O3[1];
            }
        } else {
            for ( i = 0; i < Input_Opt.PARAMETER_BACKG_O3.size(); i++ )
                cases.push_back(Input_Opt.PARAMETER_BACKG_O3[i]);
        }

        if ( Input_Opt.PARAMETER_BACKG_O3_UNIT.compare( "ppb" ) == 0 ) {
            /* Do nothing. Default unit */
        } else if ( Input_Opt.PARAMETER_BACKG_O3_UNIT.compare( "ppt" ) == 0 ) {
            /* Convert ppt to ppb */
            for ( i = 0; i < cases.size(); i++ )
                cases[i] *= 1.00E-03;
        } else if ( Input_Opt.PARAMETER_BACKG_O3_UNIT.compare( "ppm" ) == 0 ) {
            /* Convert ppm to ppb */
            for ( i = 0; i < cases.size(); i++ )
                cases[i] *= 1.00E+03;
        } else {
            std::cout << " Unknown unit for variable 'BACKG_O3': ";
            std::cout << Input_Opt.PARAMETER_BACKG_O3_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_BACKG_O3_UNIT = "ppb";

        z.push_back( Vector_1D(cases.size() ) );
        for ( i = 0; i < cases.size(); i++ )
            z[0][i] = cases[i];
        nCases *= cases.size();

        u = Copy_blocked(y,z[0].size());
        v = Copy_interleaved(z,y[0].size());

        for ( i = 0; i < counter; i++ )
            y[i].clear();
        z[0].clear();
        y.clear(); z.clear();

        counter += 1;
        for ( i = 0; i < counter; i++ )
            y.push_back(Vector_1D( nCases ));

        for ( i = 0; i < nCases; i++ ) {
            for ( j = 0; j < counter - 1; j++ )
                y[j][i] = u[j][i];
            y[counter-1][i] = v[0][i];
        }
        cases.clear();

        /* ======================================================================= */
        /* ---- BACKGROUND CO MIXING RATIO --------------------------------------- */
        /* ---- Accepted units are: ppb (default), ppt, ppm                        */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_BACKG_CO_RANGE ) {
            currVal = Input_Opt.PARAMETER_BACKG_CO[0];
            while ( currVal <= Input_Opt.PARAMETER_BACKG_CO[2] ) {
                cases.push_back( currVal );
                currVal += Input_Opt.PARAMETER_BACKG_CO[1];
            }
        } else {
            for ( i = 0; i < Input_Opt.PARAMETER_BACKG_CO.size(); i++ )
                cases.push_back(Input_Opt.PARAMETER_BACKG_CO[i]);
        }

        if ( Input_Opt.PARAMETER_BACKG_CO_UNIT.compare( "ppb" ) == 0 ) {
            /* Do nothing. Default unit */
        } else if ( Input_Opt.PARAMETER_BACKG_CO_UNIT.compare( "ppt" ) == 0 ) {
            /* Convert ppt to ppb */
            for ( i = 0; i < cases.size(); i++ )
                cases[i] *= 1.00E-03;
        } else if ( Input_Opt.PARAMETER_BACKG_CO_UNIT.compare( "ppm" ) == 0 ) {
            /* Convert ppm to ppb */
            for ( i = 0; i < cases.size(); i++ )
                cases[i] *= 1.00E+03;
        } else {
            std::cout << " Unknown unit for variable 'BACKG_CO': ";
            std::cout << Input_Opt.PARAMETER_BACKG_CO_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_BACKG_CO_UNIT = "ppb";

        z.push_back( Vector_1D(cases.size() ) );
        for ( i = 0; i < cases.size(); i++ )
            z[0][i] = cases[i];
        nCases *= cases.size();

        u = Copy_blocked(y,z[0].size());
        v = Copy_interleaved(z,y[0].size());

        for ( i = 0; i < counter; i++ )
            y[i].clear();
        z[0].clear();
        y.clear(); z.clear();

        counter += 1;
        for ( i = 0; i < counter; i++ )
            y.push_back(Vector_1D( nCases ));

        for ( i = 0; i < nCases; i++ ) {
            for ( j = 0; j < counter - 1; j++ )
                y[j][i] = u[j][i];
            y[counter-1][i] = v[0][i];
        }
        cases.clear();

        /* ======================================================================= */
        /* ---- BACKGROUND CH4 MIXING RATIO -------------------------------------- */
        /* ---- Accepted units are: ppb (default), ppt, ppm                        */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_BACKG_CH4_RANGE ) {
            currVal = Input_Opt.PARAMETER_BACKG_CH4[0];
            while ( currVal <= Input_Opt.PARAMETER_BACKG_CH4[2] ) {
                cases.push_back( currVal );
                currVal += Input_Opt.PARAMETER_BACKG_CH4[1];
            }
        } else {
            for ( i = 0; i < Input_Opt.PARAMETER_BACKG_CH4.size(); i++ )
                cases.push_back(Input_Opt.PARAMETER_BACKG_CH4[i]);
        }

        if ( Input_Opt.PARAMETER_BACKG_CH4_UNIT.compare( "ppb" ) == 0 ) {
            /* Do nothing. Default unit */
        } else if ( Input_Opt.PARAMETER_BACKG_CH4_UNIT.compare( "ppt" ) == 0 ) {
            /* Convert ppt to ppb */
            for ( i = 0; i < cases.size(); i++ )
                cases[i] *= 1.00E-03;
        } else if ( Input_Opt.PARAMETER_BACKG_CH4_UNIT.compare( "ppm" ) == 0 ) {
            /* Convert ppm to ppb */
            for ( i = 0; i < cases.size(); i++ )
                cases[i] *= 1.00E+03;
        } else {
            std::cout << " Unknown unit for variable 'BACKG_CH4': ";
            std::cout << Input_Opt.PARAMETER_BACKG_CH4_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_BACKG_CH4_UNIT = "ppb";

        z.push_back( Vector_1D(cases.size() ) );
        for ( i = 0; i < cases.size(); i++ )
            z[0][i] = cases[i];
        nCases *= cases.size();

        u = Copy_blocked(y,z[0].size());
        v = Copy_interleaved(z,y[0].size());

        for ( i = 0; i < counter; i++ )
            y[i].clear();
        z[0].clear();
        y.clear(); z.clear();

        counter += 1;
        for ( i = 0; i < counter; i++ )
            y.push_back(Vector_1D( nCases ));

        for ( i = 0; i < nCases; i++ ) {
            for ( j = 0; j < counter - 1; j++ )
                y[j][i] = u[j][i];
            y[counter-1][i] = v[0][i];
        }
        cases.clear();

        /* ======================================================================= */
        /* ---- BACKGROUND SO2 MIXING RATIO -------------------------------------- */
        /* ---- Accepted units are: ppb (default), ppt, ppm                        */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_BACKG_SO2_RANGE ) {
            currVal = Input_Opt.PARAMETER_BACKG_SO2[0];
            while ( currVal <= Input_Opt.PARAMETER_BACKG_SO2[2] ) {
                cases.push_back( currVal );
                currVal += Input_Opt.PARAMETER_BACKG_SO2[1];
            }
        } else {
            for ( i = 0; i < Input_Opt.PARAMETER_BACKG_SO2.size(); i++ )
                cases.push_back(Input_Opt.PARAMETER_BACKG_SO2[i]);
        }

        if ( Input_Opt.PARAMETER_BACKG_SO2_UNIT.compare( "ppb" ) == 0 ) {
            /* Do nothing. Default unit */
        } else if ( Input_Opt.PARAMETER_BACKG_SO2_UNIT.compare( "ppt" ) == 0 ) {
            /* Convert ppt to ppb */
            for ( i = 0; i < cases.size(); i++ )
                cases[i] *= 1.00E-03;
        } else if ( Input_Opt.PARAMETER_BACKG_SO2_UNIT.compare( "ppm" ) == 0 ) {
            /* Convert ppm to ppb */
            for ( i = 0; i < cases.size(); i++ )
                cases[i] *= 1.00E+03;
        } else {
            std::cout << " Unknown unit for variable 'BACKG_SO2': ";
            std::cout << Input_Opt.PARAMETER_BACKG_SO2_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_BACKG_SO2_UNIT = "ppb";

        z.push_back( Vector_1D(cases.size() ) );
        for ( i = 0; i < cases.size(); i++ )
            z[0][i] = cases[i];
        nCases *= cases.size();

        u = Copy_blocked(y,z[0].size());
        v = Copy_interleaved(z,y[0].size());

        for ( i = 0; i < counter; i++ )
            y[i].clear();
        z[0].clear();
        y.clear(); z.clear();

        counter += 1;
        for ( i = 0; i < counter; i++ )
            y.push_back(Vector_1D( nCases ));

        for ( i = 0; i < nCases; i++ ) {
            for ( j = 0; j < counter - 1; j++ )
                y[j][i] = u[j][i];
            y[counter-1][i] = v[0][i];
        }

        if ( VERBOSE ) {
            for ( i = 0; i < y.size(); i++ ) {
                for ( j = 0; j < y[0].size(); j++ )
                    std::cout << y[i][j] << ", ";
                std::cout << std::endl;
            }
        }
        cases.clear();

        /* ======================================================================= */
        /* ---- FLIGHT SPEED  ---------------------------------------------------- */
        /* ---- Accepted units are: m/s (default)                                  */
        /* ======================================================================= */
	
        if ( Input_Opt.PARAMETER_FSPEED_RANGE ) {
            currVal = Input_Opt.PARAMETER_FSPEED[0];
            while ( currVal <= Input_Opt.PARAMETER_FSPEED[2] ) {
                cases.push_back( currVal );
                currVal += Input_Opt.PARAMETER_FSPEED[1];
            }
        } else {
            for ( i = 0; i < Input_Opt.PARAMETER_FSPEED.size(); i++ )
                cases.push_back(Input_Opt.PARAMETER_FSPEED[i]);
        }

        if ( Input_Opt.PARAMETER_FSPEED_UNIT.compare( "m/s" ) == 0 ) {
            /* Do nothing. Default unit */
        } else {
            std::cout << " Unknown unit for variable 'FSPEED': ";
            std::cout << Input_Opt.PARAMETER_FSPEED_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_FSPEED_UNIT = "m/s";

        z.push_back( Vector_1D(cases.size() ) );
        nCases *= cases.size();

        u = Copy_blocked(y,z[0].size());
        v = Copy_interleaved(z,y[0].size());

        for ( i = 0; i < counter; i++ )
            y[i].clear();
        z[0].clear();
        y.clear(); z.clear();

        counter += 1;
        for ( i = 0; i < counter; i++ )
            y.push_back(Vector_1D( nCases ));

        for ( i = 0; i < nCases; i++ ) {
            for ( j = 0; j < counter - 1; j++ )
                y[j][i] = u[j][i];
            y[counter-1][i] = v[0][i];
        }
        cases.clear();

        /* ======================================================================= */
        /* ---- NUMBER OF ENGINES  ---------------------------------------------------- */
        /* ---- Accepted units are: 2/4 (default)                                  */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_NUMENG_RANGE ) {
            currVal = Input_Opt.PARAMETER_NUMENG[0];
            while ( currVal <= Input_Opt.PARAMETER_NUMENG[2] ) {
                cases.push_back( currVal );
                currVal += Input_Opt.PARAMETER_NUMENG[1];
            }
        } else {
            for ( i = 0; i < Input_Opt.PARAMETER_NUMENG.size(); i++ )
                cases.push_back(Input_Opt.PARAMETER_NUMENG[i]);
        }

        if ( Input_Opt.PARAMETER_NUMENG_UNIT.compare( "2/4" ) == 0 ) {
            /* Do nothing. Default unit */
        } else {
            std::cout << " Unknown unit for variable 'NUMENG': ";
            std::cout << Input_Opt.PARAMETER_NUMENG_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_NUMENG_UNIT = "2/4";

        z.push_back( Vector_1D(cases.size() ) );
        for ( i = 0; i < cases.size(); i++ )
            z[0][i] = cases[i];
        nCases *= cases.size();

        u = Copy_blocked(y,z[0].size());
        v = Copy_interleaved(z,y[0].size());

        for ( i = 0; i < counter; i++ )
            y[i].clear();
        z[0].clear();
        y.clear(); z.clear();

        counter += 1;
        for ( i = 0; i < counter; i++ )
            y.push_back(Vector_1D( nCases ));

        for ( i = 0; i < nCases; i++ ) {
            for ( j = 0; j < counter - 1; j++ )
                y[j][i] = u[j][i];
            y[counter-1][i] = v[0][i];
        }
        cases.clear();

        /* ======================================================================= */
        /* ---- WINGSPAN      ---------------------------------------------------- */
        /* ---- Accepted units are: m (default), ft                                */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_WINGSPAN_RANGE ) {
            currVal = Input_Opt.PARAMETER_WINGSPAN[0];
            while ( currVal <= Input_Opt.PARAMETER_WINGSPAN[2] ) {
                cases.push_back( currVal );
                currVal += Input_Opt.PARAMETER_WINGSPAN[1];
            }
        } else {
            for ( i = 0; i < Input_Opt.PARAMETER_WINGSPAN.size(); i++ )
                cases.push_back(Input_Opt.PARAMETER_WINGSPAN[i]);
        }

        if ( Input_Opt.PARAMETER_WINGSPAN_UNIT.compare( "m" ) == 0 ) {
            /* Do nothing. Default unit */
        } else if ( Input_Opt.PARAMETER_WINGSPAN_UNIT.compare( "ft" ) == 0 ) {
            /* Convert m to ft */
            for ( i = 0; i < cases.size(); i++ )
                cases[i] /= 3.281;
        } else {
            std::cout << " Unknown unit for variable 'WINGSPAN': ";
            std::cout << Input_Opt.PARAMETER_WINGSPAN_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_WINGSPAN_UNIT = "m";

        z.push_back( Vector_1D(cases.size() ) );
        for ( i = 0; i < cases.size(); i++ )
            z[0][i] = cases[i];
        nCases *= cases.size();

        u = Copy_blocked(y,z[0].size());
        v = Copy_interleaved(z,y[0].size());

        for ( i = 0; i < counter; i++ )
            y[i].clear();
        z[0].clear();
        y.clear(); z.clear();

        counter += 1;
        for ( i = 0; i < counter; i++ )
            y.push_back(Vector_1D( nCases ));

        for ( i = 0; i < nCases; i++ ) {
            for ( j = 0; j < counter - 1; j++ )
                y[j][i] = u[j][i];
            y[counter-1][i] = v[0][i];
        }
        cases.clear();

        /* ======================================================================= */
        /* ---- CORE EXIT TEMPERATURE -------------------------------------------- */
        /* ---- Accepted units are: K (default), C                                */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_COREEXITTEMP_RANGE ) {
            currVal = Input_Opt.PARAMETER_COREEXITTEMP[0];
            while ( currVal <= Input_Opt.PARAMETER_COREEXITTEMP[2] ) {
                cases.push_back( currVal );
                currVal += Input_Opt.PARAMETER_COREEXITTEMP[1];
            }
        } else {
            for ( i = 0; i < Input_Opt.PARAMETER_COREEXITTEMP.size(); i++ )
                cases.push_back(Input_Opt.PARAMETER_COREEXITTEMP[i]);
        }

        if ( Input_Opt.PARAMETER_COREEXITTEMP_UNIT.compare( "K" ) == 0 ) {
            /* Do nothing. Default unit */
        } else if ( Input_Opt.PARAMETER_COREEXITTEMP_UNIT.compare( "C" ) == 0 ) {
            /* Convert C to K */
            for ( i = 0; i < cases.size(); i++ )
                cases[i] += 273.15;
        } else {
            std::cout << " Unknown unit for variable 'COREEXITTEMP': ";
            std::cout << Input_Opt.PARAMETER_COREEXITTEMP_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_COREEXITTEMP_UNIT = "K";

        z.push_back( Vector_1D(cases.size() ) );
        for ( i = 0; i < cases.size(); i++ )
            z[0][i] = cases[i];
        nCases *= cases.size();

        u = Copy_blocked(y,z[0].size());
        v = Copy_interleaved(z,y[0].size());

        for ( i = 0; i < counter; i++ )
            y[i].clear();
        z[0].clear();
        y.clear(); z.clear();

        counter += 1;
        for ( i = 0; i < counter; i++ )
            y.push_back(Vector_1D( nCases ));

        for ( i = 0; i < nCases; i++ ) {
            for ( j = 0; j < counter - 1; j++ )
                y[j][i] = u[j][i];
            y[counter-1][i] = v[0][i];
        }
        cases.clear();

        /* ======================================================================= */
        /* ---- Exit bypass area ---- -------------------------------------------- */
        /* ---- Accepted units are: m^2 (default)                                */
        /* ======================================================================= */

        if ( Input_Opt.PARAMETER_BYPASSAREA_RANGE ) {
            currVal = Input_Opt.PARAMETER_BYPASSAREA[0];
            while ( currVal <= Input_Opt.PARAMETER_BYPASSAREA[2] ) {
                cases.push_back( currVal );
                currVal += Input_Opt.PARAMETER_BYPASSAREA[1];
            }
        } else {
            for ( i = 0; i < Input_Opt.PARAMETER_BYPASSAREA.size(); i++ )
                cases.push_back(Input_Opt.PARAMETER_BYPASSAREA[i]);
        }

        if ( Input_Opt.PARAMETER_BYPASSAREA_UNIT.compare( "m^2" ) == 0 ) {
            /* Do nothing. Default unit */
        } else {
            std::cout << " Unknown unit for variable 'BYPASSAREA': ";
            std::cout << Input_Opt.PARAMETER_BYPASSAREA_UNIT << std::endl;
            exit(1);
        }

        /* Updating unit now that conversion has been taken care of */
        Input_Opt.PARAMETER_BYPASSAREA_UNIT = "m^2";

        z.push_back( Vector_1D(cases.size() ) );
        for ( i = 0; i < cases.size(); i++ )
            z[0][i] = cases[i];
        nCases *= cases.size();

        u = Copy_blocked(y,z[0].size());
        v = Copy_interleaved(z,y[0].size());

        for ( i = 0; i < counter; i++ )
            y[i].clear();
        z[0].clear();
        y.clear(); z.clear();

        counter += 1;
        for ( i = 0; i < counter; i++ )
            y.push_back(Vector_1D( nCases ));

        for ( i = 0; i < nCases; i++ ) {
            for ( j = 0; j < counter - 1; j++ )
                y[j][i] = u[j][i];
            y[counter-1][i] = v[0][i];
        }
        cases.clear();

        return y;

    }

} /* End of CombVec */

Vector_2D Copy_blocked( Vector_2D& m, int n )
{

    Vector_2D b;
    const unsigned int mr = m.size();
    const unsigned int mc = m[0].size();
    unsigned int i, j, k;
    Vector_1D tmp( mc*n, 0.0E+00 );

    for ( i = 0; i < mr; i++ )
        b.push_back(tmp);

//    unsigned int ind[mc];
    unsigned int *ind = new unsigned int[mc];

    for ( i = 0; i < mc; i++ )
        ind[i] = i;

    for ( j = 0; j < (n-1) * mc + 1; j+=mc ) {
        for ( i = 0; i < mr; i++ ) {
            for ( k = 0; k < mc; k++ ) {
                b[i][ind[k] + j] = m[i][k];
            }
        }
    }

    delete[] ind;

    return b;

} /* End of Copy_blocked */

Vector_2D Copy_interleaved( Vector_2D& m, int n )
{

    Vector_2D b;
    const unsigned int mr = m.size();
    const unsigned int mc = m[0].size();
    unsigned int i, j, k;
    Vector_1D tmp( mc, 0.0E+00 );

    for ( i = 0; i < mr*n; i++ )
        b.push_back(tmp);

//    unsigned int ind[mr];
    unsigned int *ind = new unsigned int[mr];

    for ( i = 0; i < mr; i++ )
        ind[i] = i;

    for ( i = 0; i < (n-1) * mr + 1; i+=mr ) {
        for ( j = 0; j < mc; j++ ) {
            for ( k = 0; k < mr; k++ ) {
                b[ind[k] + i][j] = m[k][j]; //m[k][j]
            }
        }
    }

    delete[] ind;

    return Reshape_Vector( b, mr, n*mc );

} /* End of Copy_interleaved */

Vector_2D Reshape_Vector( Vector_2D& vector_2D, int n_x, int n_y )
{

    Vector_2D output;
    int size_x = vector_2D.size();
    int size_y = vector_2D[0].size();

    if ( n_x * n_y != size_x * size_y ) {
        std::cout << "Invalid dimensions specified" << std::endl;
        return output;
    }
    else {
        int counter;
        for ( counter = 0; counter < n_x; counter++ )
            output.push_back(Vector_1D(n_y));

        int counter_x = 0;
        int counter_y = 0;
        int orig_counter_x = 0;
        int orig_counter_y = 0;
        for ( counter = 0; counter < n_x * n_y; counter++ ) {
            counter_x = counter%n_x;
            orig_counter_x = counter%size_x;
            output[counter_x][counter_y] = vector_2D[orig_counter_x][orig_counter_y];
            if ( counter_x == n_x - 1 )
                counter_y += 1;
            if ( orig_counter_x == size_x - 1 )
                orig_counter_y += 1;
        }
    }

    return output;

} /* End of Reshape_Vector */

void Are_Flags_Valid( const OptInput &Input_Opt )
{

    /* Are_Flags_Valid checks to make sure that flags are valid and
     * do not conflict */

    /* If MONTECARLO flag is on and parameters are read from file, then print
     * statement and exit... */
    if ( Input_Opt.PARAMETER_FILEINPUT && Input_Opt.SIMULATION_MONTECARLO ) {
        std::cout << " In Are_Flags_Valid:";
        std::cout << " PARAMETER_FILEINPUT and SIMULATION_MONTECARLO are turned on!" << std::endl;
        std::cout << " Aborting!" << std::endl;
    }

    /* If chemistry is turned off and P&L chemistry rate output is turned on,
     * print statement and exit... */
    if ( !Input_Opt.CHEMISTRY_CHEMISTRY && ( Input_Opt.PL_PL || \
                                             Input_Opt.PL_O3 ) ) {
        std::cout << " In Are_Flags_Valid:";
        std::cout << " CHEMISTRY is turned off while P&L rates output is turned on!" << std::endl;
        std::cout << " Aborting!" << std::endl;
        exit(-1);
    }

    /* The depth of the moist layer and the temperature lapse rate cannot be
     * set simultaneously as they depend on each other. Print statement and
     * exit... */
    if ( !Input_Opt.MET_LOADMET && Input_Opt.MET_FIXDEPTH && Input_Opt.MET_FIXLAPSERATE ) {
        std::cout << " In Are_Flags_Valid:";
        std::cout << " LOADMET is turned off and both FIXDEPTH and FIXLAPSERATE are turned on!" << std::endl;
        std::cout << " Aborting!" << std::endl;
        exit(-1);
    }

    /* Treatment of other conflicts goes here ... */


} /* End of AreFlagsValid */

/* End of Input_Reader.cpp */
