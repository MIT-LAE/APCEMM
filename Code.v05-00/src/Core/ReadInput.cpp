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
const char* SPACE = " ";
const char* TAB   = "\t";
const char* COLON = ":";
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
        std::cout << " Simulation Directory is not defined!" << std::endl;
        std::cout << " Make sure that the variable 'APCEMM_runDir' is exported" << std::endl;
        exit(1);
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
    std::cout << " \nReading from: " << fullPath << "\n" << std::endl;
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
            std::cout << line << std::endl;

        if ( strstr( line.c_str(), "SIMULATION MENU" ) != NULL ) {
            Read_Simulation_Menu( Input_Opt, RC ); 

        } else if ( strstr( line.c_str(), "PARAMETER SWEEP" ) != NULL ) {
            Read_Parameters( Input_Opt, RC );

        } else if ( strstr( line.c_str(), "TRANSPORT MENU" ) != NULL ) {
            Read_Transport_Menu( Input_Opt, RC );

        } else if ( strstr( line.c_str(), "CHEMISTRY MENU" ) != NULL ) {
            Read_Chemistry_Menu( Input_Opt, RC );

        } else if ( strstr( line.c_str(), "AEROSOL MENU" ) != NULL ) {
            Read_Aerosol_Menu( Input_Opt, RC );
        
        } else if ( strstr( line.c_str(), "METEOROLOGY MENU" ) != NULL ) {
            Read_Meteorology_Menu( Input_Opt, RC );
        
        } else if ( strstr( line.c_str(), "DIAGNOSTIC MENU" ) != NULL ) {
            Read_Diagnostic_Menu( Input_Opt, RC );
        
        } else if ( strstr( line.c_str(), "TIMESERIES MENU" ) != NULL ) {
            Read_Timeseries_Menu( Input_Opt, RC );

        } else if ( strstr( line.c_str(), "PROD & LOSS MENU" ) != NULL ) {
            Read_PL_Menu( Input_Opt, RC );

        } else if ( strstr( line.c_str(), "END OF FILE" ) != NULL )
            break;
                
    }

    inputFile.close();

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

    while ( (pos = line2split.find(delimiter) ) != std::string::npos) {
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

    if ( ( strcmp(tokens[0].c_str(), "T" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" ) == 0 ) )
        Input_Opt.SIMULATION_PARAMETER_SWEEP = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" ) == 0 ) )
        Input_Opt.SIMULATION_PARAMETER_SWEEP = 0;
    else {
        std::cout << " Wrong input for: " << "Parameter sweep?" << std::endl;
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

    if ( ( strcmp(tokens[0].c_str(), "T" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" ) == 0 ) )
        Input_Opt.SIMULATION_OVERWRITE = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" ) == 0 ) )
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

    if ( ( strcmp(tokens[0].c_str(), "T" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" ) == 0 ) )
        Input_Opt.SIMULATION_SAVE_FORWARD = 0;
    else if ( ( strcmp(tokens[0].c_str(), "F" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" ) == 0 ) )
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

    if ( ( strcmp(tokens[0].c_str(), "T" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" ) == 0 ) )
        Input_Opt.SIMULATION_ADJOINT = 0;
    else if ( ( strcmp(tokens[0].c_str(), "F" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" ) == 0 ) )
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

    /* Return success */
    RC = SUCCESS; 

    /* ==================================================== */
    /* Print to screen                                      */
    /* ==================================================== */

    std::cout << " ------------------------+---------------------------------------- "  << std::endl;
    std::cout << " %%% SIMULATION MENU %%% :"                                           << std::endl;
    std::cout << " ------------------------+---------------------------------------- "  << std::endl;
    std::cout << " Parameter sweep?        : " << Input_Opt.SIMULATION_PARAMETER_SWEEP  << std::endl;
    std::cout << " Output folder           : " << Input_Opt.SIMULATION_OUTPUT_FOLDER    << std::endl;
    std::cout << "  => Overwrite? if exists: " << Input_Opt.SIMULATION_OVERWRITE        << std::endl;
    std::cout << " Run directory           : " << Input_Opt.SIMULATION_RUN_DIRECTORY    << std::endl;
    std::cout << " Input backgrd condition : " << Input_Opt.SIMULATION_INPUT_BACKG_COND << std::endl;
    std::cout << " Save Forward results    : " << Input_Opt.SIMULATION_SAVE_FORWARD     << std::endl;
    std::cout << "  => netCDF file name    : " << Input_Opt.SIMULATION_FORWARD_FILENAME << std::endl;
    std::cout << " Turn on adjoint optim.  : " << Input_Opt.SIMULATION_ADJOINT          << std::endl;
    std::cout << "  => netCDF file name    : " << Input_Opt.SIMULATION_ADJOINT_FILENAME << std::endl;
    std::cout << " ------------------------+---------------------------------------- "  << std::endl;
    
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

    /* Skip menu header lines */
    getline( inputFile, line, '\n' );
    getline( inputFile, line, '\n' );

    /* Variable ranges can be defined in two ways:
     *  - min:step:max
     *  - val1 val2 val3 ... */

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
    if ( found != std::string::npos) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

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
        Input_Opt.PARAMETER_TEMPERATURE_RANGE = 1;
    } 
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for" << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_TEMPERATURE_RANGE = 0;
    }
    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) ) {
        std::cout << " APCEMM can not accept multiple cases when the 'parameter sweep?' argument is turned off! Aborting.";
        exit(1);
    }
  
    /* Find unit in between "[" and "]" */ 
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_TEMPERATURE_UNIT.assign( unit );

    /* Store in values for variable */
    for ( unsigned int i = 0; i < tokens.size(); i++ ) {
        try {
            value = std::stod( tokens[i] );
            if ( value > 0.0E+00 )
                Input_Opt.PARAMETER_TEMPERATURE.push_back( value );
            else {
                std::cout << " Wrong input for: " << variable << std::endl;
                std::cout << " Index needs to be positive" << std::endl;
                exit(1);
            }
        } catch(std::exception& e) {
            std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
            exit(1);
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
    if ( found != std::string::npos) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

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
        Input_Opt.PARAMETER_RHW_RANGE = 1;
    } 
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for" << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_RHW_RANGE = 0;
    }
    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) ) {
        std::cout << " APCEMM can not accept multiple cases when the 'parameter sweep?' argument is turned off! Aborting.";
        exit(1);
    }
    
    /* Find unit in between "[" and "]" */ 
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_RHW_UNIT.assign( unit );

    /* Store in values for variable */
    for ( unsigned int i = 0; i < tokens.size(); i++ ) {
        try {
            value = std::stod( tokens[i] );
            if ( value >= 0.0E+00 )
                Input_Opt.PARAMETER_RHW.push_back( value );
            else {
                std::cout << " Wrong input for: " << variable << std::endl;
                std::cout << " Index needs to be positive" << std::endl;
                exit(1);
            }
        } catch(std::exception& e) {
            std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
            exit(1);
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
    if ( found != std::string::npos) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

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
        Input_Opt.PARAMETER_LONGITUDE_RANGE = 1;
    } 
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for" << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_LONGITUDE_RANGE = 0;
    }
    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) ) {
        std::cout << " APCEMM can not accept multiple cases when the 'parameter sweep?' argument is turned off! Aborting.";
        exit(1);
    }
   
    /* Find unit in between "[" and "]" */ 
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_LONGITUDE_UNIT.assign( unit );

    /* Store in values for variable */
    for ( unsigned int i = 0; i < tokens.size(); i++ ) {
        try {
            Input_Opt.PARAMETER_LONGITUDE.push_back(std::stod( tokens[i] ));
        } catch(std::exception& e) {
            std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
            exit(1);
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
    if ( found != std::string::npos) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

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
        Input_Opt.PARAMETER_LATITUDE_RANGE = 1;
    } 
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for" << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_LATITUDE_RANGE = 0;
    }
    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) ) {
        std::cout << " APCEMM can not accept multiple cases when the 'parameter sweep?' argument is turned off! Aborting.";
        exit(1);
    }
    
    /* Find unit in between "[" and "]" */ 
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_LATITUDE_UNIT.assign( unit );

    /* Store in values for variable */
    for ( unsigned int i = 0; i < tokens.size(); i++ ) {
        try {
            Input_Opt.PARAMETER_LATITUDE.push_back(std::stod( tokens[i] ));
        } catch(std::exception& e) {
            std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
            exit(1);
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
    if ( found != std::string::npos) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

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
        Input_Opt.PARAMETER_PRESSURE_RANGE = 1;
    } 
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for" << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_PRESSURE_RANGE = 0;
    }
    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) ) {
        std::cout << " APCEMM can not accept multiple cases when the 'parameter sweep?' argument is turned off! Aborting.";
        exit(1);
    }
    
    /* Find unit in between "[" and "]" */ 
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_PRESSURE_UNIT.assign( unit );

    /* Store in values for variable */
    for ( unsigned int i = 0; i < tokens.size(); i++ ) {
        try {
            value = std::stod( tokens[i] );
            if ( value > 0.0E+00 )
                Input_Opt.PARAMETER_PRESSURE.push_back( value );
            else {
                std::cout << " Wrong input for: " << variable << std::endl;
                std::cout << " Index needs to be positive" << std::endl;
                exit(1);
            }
        } catch(std::exception& e) {
            std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
            exit(1);
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
    if ( found != std::string::npos) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

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
        Input_Opt.PARAMETER_EDAY_RANGE = 1;
    } 
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for" << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_EDAY_RANGE = 0;
    }
    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) ) {
        std::cout << " APCEMM can not accept multiple cases when the 'parameter sweep?' argument is turned off! Aborting.";
        exit(1);
    }
    
    /* Find unit in between "[" and "]" */ 
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_EDAY_UNIT.assign( unit );

    /* Store in values for variable */
    for ( unsigned int i = 0; i < tokens.size(); i++ ) {
        try {
            value = ( std::stoi( tokens[i] ) % 365 );
            if ( value >= 0.0E+00 )
                Input_Opt.PARAMETER_EDAY.push_back( value );
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
    if ( found != std::string::npos) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

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
        Input_Opt.PARAMETER_ETIME_RANGE = 1;
    } 
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for" << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_ETIME_RANGE = 0;
    }
    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) ) {
        std::cout << " APCEMM can not accept multiple cases when the 'parameter sweep?' argument is turned off! Aborting.";
        exit(1);
    }
    
    /* Find unit in between "[" and "]" */ 
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_ETIME_UNIT.assign( unit );

    /* Store in values for variable */
    for ( unsigned int i = 0; i < tokens.size(); i++ ) {
        try {
            value = std::fmod(std::stod(tokens[i]), 2.40E+01);
            if ( value >= 0.0E+00 )
                Input_Opt.PARAMETER_ETIME.push_back( value );
            else {
                std::cout << " Wrong input for: " << variable << std::endl;
                std::cout << " Index needs to be positive" << std::endl;
                exit(1);
            }
        } catch(std::exception& e) {
            std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
            exit(1);
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
    if ( found != std::string::npos) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

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
        Input_Opt.PARAMETER_BACKG_NOX_RANGE = 1;
    } 
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for" << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_BACKG_NOX_RANGE = 0;
    }
    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) ) {
        std::cout << " APCEMM can not accept multiple cases when the 'parameter sweep?' argument is turned off! Aborting.";
        exit(1);
    }
    
    /* Find unit in between "[" and "]" */ 
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_BACKG_NOX_UNIT.assign( unit );

    /* Store in values for variable */
    for ( unsigned int i = 0; i < tokens.size(); i++ ) {
        try {
            value = std::stod( tokens[i] );
            if ( value >= 0.0E+00 )
                Input_Opt.PARAMETER_BACKG_NOX.push_back( value );
            else {
                std::cout << " Wrong input for: " << variable << std::endl;
                std::cout << " Index needs to be positive" << std::endl;
                exit(1);
            }
        } catch(std::exception& e) {
            std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
            exit(1);
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
    if ( found != std::string::npos) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

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
        Input_Opt.PARAMETER_BACKG_HNO3_RANGE = 1;
    } 
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for" << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_BACKG_HNO3_RANGE = 0;
    }
    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) ) {
        std::cout << " APCEMM can not accept multiple cases when the 'parameter sweep?' argument is turned off! Aborting.";
        exit(1);
    }
    
    /* Find unit in between "[" and "]" */ 
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_BACKG_HNO3_UNIT.assign( unit );

    /* Store in values for variable */
    for ( unsigned int i = 0; i < tokens.size(); i++ ) {
        try {
            value = std::stod( tokens[i] );
            if ( value >= 0.0E+00 )
                Input_Opt.PARAMETER_BACKG_HNO3.push_back( value );
            else {
                std::cout << " Wrong input for: " << variable << std::endl;
                std::cout << " Index needs to be positive" << std::endl;
                exit(1);
            }
        } catch(std::exception& e) {
            std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
            exit(1);
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
    if ( found != std::string::npos) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

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
        Input_Opt.PARAMETER_BACKG_O3_RANGE = 1;
    } 
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for" << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_BACKG_O3_RANGE = 0;
    }
    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) ) {
        std::cout << " APCEMM can not accept multiple cases when the 'parameter sweep?' argument is turned off! Aborting.";
        exit(1);
    }
    
    /* Find unit in between "[" and "]" */ 
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_BACKG_O3_UNIT.assign( unit );

    /* Store in values for variable */
    for ( unsigned int i = 0; i < tokens.size(); i++ ) {
        try {
            value = std::stod( tokens[i] );
            if ( value >= 0.0E+00 )
                Input_Opt.PARAMETER_BACKG_O3.push_back( value );
            else {
                std::cout << " Wrong input for: " << variable << std::endl;
                std::cout << " Index needs to be positive" << std::endl;
                exit(1);
            }
        } catch(std::exception& e) {
            std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
            exit(1);
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
    if ( found != std::string::npos) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

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
        Input_Opt.PARAMETER_BACKG_CO_RANGE = 1;
    } 
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for" << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_BACKG_CO_RANGE = 0;
    }
    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) ) {
        std::cout << " APCEMM can not accept multiple cases when the 'parameter sweep?' argument is turned off! Aborting.";
        exit(1);
    }
    
    /* Find unit in between "[" and "]" */ 
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_BACKG_CO_UNIT.assign( unit );

    /* Store in values for variable */
    for ( unsigned int i = 0; i < tokens.size(); i++ ) {
        try {
            value = std::stod( tokens[i] );
            if ( value >= 0.0E+00 )
                Input_Opt.PARAMETER_BACKG_CO.push_back( value );
            else {
                std::cout << " Wrong input for: " << variable << std::endl;
                std::cout << " Index needs to be positive" << std::endl;
                exit(1);
            }
        } catch(std::exception& e) {
            std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
            exit(1);
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
    if ( found != std::string::npos) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

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
        Input_Opt.PARAMETER_BACKG_CH4_RANGE = 1;
    } 
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for" << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_BACKG_CH4_RANGE = 0;
    }
    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) ) {
        std::cout << " APCEMM can not accept multiple cases when the 'parameter sweep?' argument is turned off! Aborting.";
        exit(1);
    }
    
    /* Find unit in between "[" and "]" */ 
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_BACKG_CH4_UNIT.assign( unit );

    /* Store in values for variable */
    for ( unsigned int i = 0; i < tokens.size(); i++ ) {
        try {
            value = std::stod( tokens[i] );
            if ( value >= 0.0E+00 )
                Input_Opt.PARAMETER_BACKG_CH4.push_back( value );
            else {
                std::cout << " Wrong input for: " << variable << std::endl;
                std::cout << " Index needs to be positive" << std::endl;
                exit(1);
            }
        } catch(std::exception& e) {
            std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
            exit(1);
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
    if ( found != std::string::npos) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

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
        Input_Opt.PARAMETER_BACKG_SO2_RANGE = 1;
    } 
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for" << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_BACKG_SO2_RANGE = 0;
    }
    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) ) {
        std::cout << " APCEMM can not accept multiple cases when the 'parameter sweep?' argument is turned off! Aborting.";
        exit(1);
    }
    
    /* Find unit in between "[" and "]" */ 
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_BACKG_SO2_UNIT.assign( unit );

    /* Store in values for variable */
    for ( unsigned int i = 0; i < tokens.size(); i++ ) {
        try {
            value = std::stod( tokens[i] );
            if ( value >= 0.0E+00 )
                Input_Opt.PARAMETER_BACKG_SO2.push_back( value );
            else {
                std::cout << " Wrong input for: " << variable << std::endl;
                std::cout << " Index needs to be positive" << std::endl;
                exit(1);
            }
        } catch(std::exception& e) {
            std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
            exit(1);
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
    if ( found != std::string::npos) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

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
        Input_Opt.PARAMETER_EI_NOX_RANGE = 1;
    } 
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for" << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_EI_NOX_RANGE = 0;
    }
    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) ) {
        std::cout << " APCEMM can not accept multiple cases when the 'parameter sweep?' argument is turned off! Aborting.";
        exit(1);
    }
    
    /* Find unit in between "[" and "]" */ 
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_EI_NOX_UNIT.assign( unit );

    /* Store in values for variable */
    for ( unsigned int i = 0; i < tokens.size(); i++ ) {
        try {
            value = std::stod( tokens[i] );
            if ( value >= 0.0E+00 )
                Input_Opt.PARAMETER_EI_NOX.push_back( value );
            else {
                std::cout << " Wrong input for: " << variable << std::endl;
                std::cout << " Index needs to be positive" << std::endl;
                exit(1);
            }
        } catch(std::exception& e) {
            std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
            exit(1);
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
    if ( found != std::string::npos) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

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
        Input_Opt.PARAMETER_EI_CO_RANGE = 1;
    } 
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for" << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_EI_CO_RANGE = 0;
    }
    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) ) {
        std::cout << " APCEMM can not accept multiple cases when the 'parameter sweep?' argument is turned off! Aborting.";
        exit(1);
    }
    
    /* Find unit in between "[" and "]" */ 
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_EI_CO_UNIT.assign( unit );

    /* Store in values for variable */
    for ( unsigned int i = 0; i < tokens.size(); i++ ) {
        try {
            value = std::stod( tokens[i] );
            if ( value >= 0.0E+00 )
                Input_Opt.PARAMETER_EI_CO.push_back( value );
            else {
                std::cout << " Wrong input for: " << variable << std::endl;
                std::cout << " Index needs to be positive" << std::endl;
                exit(1);
            }
        } catch(std::exception& e) {
            std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
            exit(1);
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
    if ( found != std::string::npos) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

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
        Input_Opt.PARAMETER_EI_UHC_RANGE = 1;
    } 
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for" << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_EI_UHC_RANGE = 0;
    }
    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) ) {
        std::cout << " APCEMM can not accept multiple cases when the 'parameter sweep?' argument is turned off! Aborting.";
        exit(1);
    }
    
    /* Find unit in between "[" and "]" */ 
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_EI_UHC_UNIT.assign( unit );

    /* Store in values for variable */
    for ( unsigned int i = 0; i < tokens.size(); i++ ) {
        try {
            value = std::stod( tokens[i] );
            if ( value >= 0.0E+00 )
                Input_Opt.PARAMETER_EI_UHC.push_back( value );
            else {
                std::cout << " Wrong input for: " << variable << std::endl;
                std::cout << " Index needs to be positive" << std::endl;
                exit(1);
            }
        } catch(std::exception& e) {
            std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
            exit(1);
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
    if ( found != std::string::npos) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

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
        Input_Opt.PARAMETER_EI_SO2_RANGE = 1;
    } 
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for" << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_EI_SO2_RANGE = 0;
    }
    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) ) {
        std::cout << " APCEMM can not accept multiple cases when the 'parameter sweep?' argument is turned off! Aborting.";
        exit(1);
    }
    
    /* Find unit in between "[" and "]" */ 
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_EI_SO2_UNIT.assign( unit );

    /* Store in values for variable */
    for ( unsigned int i = 0; i < tokens.size(); i++ ) {
        try {
            value = std::stod( tokens[i] );
            if ( value >= 0.0E+00 )
                Input_Opt.PARAMETER_EI_SO2.push_back( value );
            else {
                std::cout << " Wrong input for: " << variable << std::endl;
                std::cout << " Index needs to be positive" << std::endl;
                exit(1);
            }
        } catch(std::exception& e) {
            std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
            exit(1);
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
    if ( found != std::string::npos) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

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
        Input_Opt.PARAMETER_EI_SO2TOSO4_RANGE = 1;
    } 
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for" << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_EI_SO2TOSO4_RANGE = 0;
    }
    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) ) {
        std::cout << " APCEMM can not accept multiple cases when the 'parameter sweep?' argument is turned off! Aborting.";
        exit(1);
    }
    
    /* Find unit in between "[" and "]" */ 
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_EI_SO2TOSO4_UNIT.assign( unit );

    /* Store in values for variable */
    for ( unsigned int i = 0; i < tokens.size(); i++ ) {
        try {
            value = std::stod( tokens[i] );
            if ( value >= 0.0E+00 )
                Input_Opt.PARAMETER_EI_SO2TOSO4.push_back( value );
            else {
                std::cout << " Wrong input for: " << variable << std::endl;
                std::cout << " Index needs to be positive" << std::endl;
                exit(1);
            }
        } catch(std::exception& e) {
            std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
            exit(1);
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
    if ( found != std::string::npos) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

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
        Input_Opt.PARAMETER_EI_SOOT_RANGE = 1;
    } 
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for" << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_EI_SOOT_RANGE = 0;
    }
    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) ) {
        std::cout << " APCEMM can not accept multiple cases when the 'parameter sweep?' argument is turned off! Aborting.";
        exit(1);
    }
    
    /* Find unit in between "[" and "]" */ 
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_EI_SOOT_UNIT.assign( unit );

    /* Store in values for variable */
    for ( unsigned int i = 0; i < tokens.size(); i++ ) {
        try {
            value = std::stod( tokens[i] );
            if ( value >= 0.0E+00 )
                Input_Opt.PARAMETER_EI_SOOT.push_back( value );
            else {
                std::cout << " Wrong input for: " << variable << std::endl;
                std::cout << " Index needs to be positive" << std::endl;
                exit(1);
            }
        } catch(std::exception& e) {
            std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
            exit(1);
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
    if ( found != std::string::npos) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

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
        Input_Opt.PARAMETER_EI_SOOTRAD_RANGE = 1;
    } 
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for" << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_EI_SOOTRAD_RANGE = 0;
    }
    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) ) {
        std::cout << " APCEMM can not accept multiple cases when the 'parameter sweep?' argument is turned off! Aborting.";
        exit(1);
    }
    
    /* Find unit in between "[" and "]" */ 
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_EI_SOOTRAD_UNIT.assign( unit );

    /* Store in values for variable */
    for ( unsigned int i = 0; i < tokens.size(); i++ ) {
        try {
            value = std::stod( tokens[i] );
            if ( value > 0.0E+00 )
                Input_Opt.PARAMETER_EI_SOOTRAD.push_back( value );
            else {
                std::cout << " Wrong input for: " << variable << std::endl;
                std::cout << " Index needs to be positive" << std::endl;
                exit(1);
            }
        } catch(std::exception& e) {
            std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
            exit(1);
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
    if ( found != std::string::npos) {
        subline.erase(std::remove(subline.begin(), subline.end(), ' '), subline.end());
        /* Extract variable range */
        tokens = Split_Line( subline, COLON );

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
        Input_Opt.PARAMETER_FF_RANGE = 1;
    } 
    else {
        /* Extract variable range */
        tokens = Split_Line( subline, SPACE );

        if ( tokens.size() < 1 ) {
            std::cout << " Expected at least one value for" << variable << std::endl;
            exit(1);
        }
        Input_Opt.PARAMETER_FF_RANGE = 0;
    }
    if ( ( tokens.size() > 1 ) && ( !Input_Opt.SIMULATION_PARAMETER_SWEEP ) ) {
        std::cout << " APCEMM can not accept multiple cases when the 'parameter sweep?' argument is turned off! Aborting.";
        exit(1);
    }
    
    /* Find unit in between "[" and "]" */ 
    first = line.find("[");
    last  = line.find("]");
    unit = line.substr( first+1, last-first-1 );
    Input_Opt.PARAMETER_FF_UNIT.assign( unit );

    /* Store in values for variable */
    for ( unsigned int i = 0; i < tokens.size(); i++ ) {
        try {
            value = std::stod( tokens[i] );
            if ( value >= 0.0E+00 )
                Input_Opt.PARAMETER_FF.push_back( value );
            else {
                std::cout << " Wrong input for: " << variable << std::endl;
                std::cout << " Index needs to be positive" << std::endl;
                exit(1);
            }
        } catch(std::exception& e) {
            std::cout << " Could not convert string '" << tokens[i] << "' to double for " << variable << std::endl;
            exit(1);
        }
    }


    /* Return success */
    RC = SUCCESS; 

    /* ==================================================== */
    /* Print to screen                                      */
    /* ==================================================== */

    std::cout << " %%% PARAMETER SWEEP %%% :"                                           << std::endl;
    std::cout << " ------------------------+---------------------------------------- "  << std::endl;

    /* ---- Meteorological parameters --------------------- */
    std::cout << " Meteorological parameter:" << std::endl;
    std::cout << "  => Temperature [" << Input_Opt.PARAMETER_TEMPERATURE_UNIT << "]     : ";
    if ( Input_Opt.PARAMETER_TEMPERATURE_RANGE )
        std::cout << Input_Opt.PARAMETER_TEMPERATURE[0] << ":" << Input_Opt.PARAMETER_TEMPERATURE[1] << ":" << Input_Opt.PARAMETER_TEMPERATURE[2] << std::endl;
    else {
        for ( unsigned int i = 0; i < Input_Opt.PARAMETER_TEMPERATURE.size(); i++ )
            std::cout << Input_Opt.PARAMETER_TEMPERATURE[i] << " ";
        std::cout << std::endl;
    }
    std::cout << "  => R.Hum wrt water [" << Input_Opt.PARAMETER_RHW_UNIT << "] : ";
    if ( Input_Opt.PARAMETER_RHW_RANGE )
        std::cout << Input_Opt.PARAMETER_RHW[0] << ":" << Input_Opt.PARAMETER_RHW[1] << ":" << Input_Opt.PARAMETER_RHW[2] << std::endl;
    else {
        for ( unsigned int i = 0; i < Input_Opt.PARAMETER_RHW.size(); i++ )
            std::cout << Input_Opt.PARAMETER_RHW[i] << " ";
        std::cout << std::endl;
    }
    
    /* ---- Geographical location ------------------------- */
    std::cout << " Geographical parameters :" << std::endl;
    std::cout << "  => LON [" << Input_Opt.PARAMETER_LONGITUDE_UNIT << "]           : ";
    if ( Input_Opt.PARAMETER_LONGITUDE_RANGE )
        std::cout << Input_Opt.PARAMETER_LONGITUDE[0] << ":" << Input_Opt.PARAMETER_LONGITUDE[1] << ":" << Input_Opt.PARAMETER_LONGITUDE[2] << std::endl;
    else {
        for ( unsigned int i = 0; i < Input_Opt.PARAMETER_LONGITUDE.size(); i++ )
            std::cout << Input_Opt.PARAMETER_LONGITUDE[i] << " ";
        std::cout << std::endl;
    }
    std::cout << "  => LAT [" << Input_Opt.PARAMETER_LATITUDE_UNIT << "]           : ";
    if ( Input_Opt.PARAMETER_LATITUDE_RANGE )
        std::cout << Input_Opt.PARAMETER_LATITUDE[0] << ":" << Input_Opt.PARAMETER_LATITUDE[1] << ":" << Input_Opt.PARAMETER_LATITUDE[2] << std::endl;
    else {
        for ( unsigned int i = 0; i < Input_Opt.PARAMETER_LATITUDE.size(); i++ )
            std::cout << Input_Opt.PARAMETER_LATITUDE[i] << " ";
        std::cout << std::endl;
    }
    std::cout << "  => Pressure [" << Input_Opt.PARAMETER_PRESSURE_UNIT << "]      : ";
    if ( Input_Opt.PARAMETER_PRESSURE_RANGE )
        std::cout << Input_Opt.PARAMETER_PRESSURE[0] << ":" << Input_Opt.PARAMETER_PRESSURE[1] << ":" << Input_Opt.PARAMETER_PRESSURE[2] << std::endl;
    else {
        for ( unsigned int i = 0; i < Input_Opt.PARAMETER_PRESSURE.size(); i++ )
            std::cout << Input_Opt.PARAMETER_PRESSURE[i] << " ";
        std::cout << std::endl;
    }

    /* ---- Time of emission ------------------------------ */
    std::cout << " Time of emission        :" << std::endl;
    std::cout << "  => Emission day [" << Input_Opt.PARAMETER_EDAY_UNIT << "]: ";
    if ( Input_Opt.PARAMETER_EDAY_RANGE )
        std::cout << Input_Opt.PARAMETER_EDAY[0] << ":" << Input_Opt.PARAMETER_EDAY[1] << ":" << Input_Opt.PARAMETER_EDAY[2] << std::endl;
    else {
        for ( unsigned int i = 0; i < Input_Opt.PARAMETER_EDAY.size(); i++ )
            std::cout << Input_Opt.PARAMETER_EDAY[i] << " ";
        std::cout << std::endl;
    }
    std::cout << "  => Emission time [" << Input_Opt.PARAMETER_ETIME_UNIT << "]: ";
    if ( Input_Opt.PARAMETER_ETIME_RANGE )
        std::cout << Input_Opt.PARAMETER_ETIME[0] << ":" << Input_Opt.PARAMETER_ETIME[1] << ":" << Input_Opt.PARAMETER_ETIME[2] << std::endl;
    else {
        for ( unsigned int i = 0; i < Input_Opt.PARAMETER_ETIME.size(); i++ )
            std::cout << Input_Opt.PARAMETER_ETIME[i] << " ";
        std::cout << std::endl;
    }

    /* ---- Background mixing ratios ---------------------- */
    std::cout << " Background mix. ratios  :" << std::endl;
    std::cout << "  => NOx  [" << Input_Opt.PARAMETER_BACKG_NOX_UNIT << "]          : ";
    if ( Input_Opt.PARAMETER_BACKG_NOX_RANGE )
        std::cout << Input_Opt.PARAMETER_BACKG_NOX[0] << ":" << Input_Opt.PARAMETER_BACKG_NOX[1] << ":" << Input_Opt.PARAMETER_BACKG_NOX[2] << std::endl;
    else {
        for ( unsigned int i = 0; i < Input_Opt.PARAMETER_BACKG_NOX.size(); i++ )
            std::cout << Input_Opt.PARAMETER_BACKG_NOX[i] << " ";
        std::cout << std::endl;
    }
    std::cout << "  => HNO3 [" << Input_Opt.PARAMETER_BACKG_HNO3_UNIT << "]          : ";
    if ( Input_Opt.PARAMETER_BACKG_HNO3_RANGE )
        std::cout << Input_Opt.PARAMETER_BACKG_HNO3[0] << ":" << Input_Opt.PARAMETER_BACKG_HNO3[1] << ":" << Input_Opt.PARAMETER_BACKG_HNO3[2] << std::endl;
    else {
        for ( unsigned int i = 0; i < Input_Opt.PARAMETER_BACKG_HNO3.size(); i++ )
            std::cout << Input_Opt.PARAMETER_BACKG_HNO3[i] << " ";
        std::cout << std::endl;
    }
    std::cout << "  => O3   [" << Input_Opt.PARAMETER_BACKG_O3_UNIT << "]          : ";
    if ( Input_Opt.PARAMETER_BACKG_O3_RANGE )
        std::cout << Input_Opt.PARAMETER_BACKG_O3[0] << ":" << Input_Opt.PARAMETER_BACKG_O3[1] << ":" << Input_Opt.PARAMETER_BACKG_O3[2] << std::endl;
    else {
        for ( unsigned int i = 0; i < Input_Opt.PARAMETER_BACKG_O3.size(); i++ )
            std::cout << Input_Opt.PARAMETER_BACKG_O3[i] << " ";
        std::cout << std::endl;
    }
    std::cout << "  => CO   [" << Input_Opt.PARAMETER_BACKG_CO_UNIT << "]          : ";
    if ( Input_Opt.PARAMETER_BACKG_CO_RANGE )
        std::cout << Input_Opt.PARAMETER_BACKG_CO[0] << ":" << Input_Opt.PARAMETER_BACKG_CO[1] << ":" << Input_Opt.PARAMETER_BACKG_CO[2] << std::endl;
    else {
        for ( unsigned int i = 0; i < Input_Opt.PARAMETER_BACKG_CO.size(); i++ )
            std::cout << Input_Opt.PARAMETER_BACKG_CO[i] << " ";
        std::cout << std::endl;
    }
    std::cout << "  => CH4  [" << Input_Opt.PARAMETER_BACKG_CH4_UNIT << "]          : ";
    if ( Input_Opt.PARAMETER_BACKG_CH4_RANGE )
        std::cout << Input_Opt.PARAMETER_BACKG_CH4[0] << ":" << Input_Opt.PARAMETER_BACKG_CH4[1] << ":" << Input_Opt.PARAMETER_BACKG_CH4[2] << std::endl;
    else {
        for ( unsigned int i = 0; i < Input_Opt.PARAMETER_BACKG_CH4.size(); i++ )
            std::cout << Input_Opt.PARAMETER_BACKG_CH4[i] << " ";
        std::cout << std::endl;
    }
    std::cout << "  => SO2  [" << Input_Opt.PARAMETER_BACKG_SO2_UNIT << "]          : ";
    if ( Input_Opt.PARAMETER_BACKG_SO2_RANGE )
        std::cout << Input_Opt.PARAMETER_BACKG_SO2[0] << ":" << Input_Opt.PARAMETER_BACKG_SO2[1] << ":" << Input_Opt.PARAMETER_BACKG_SO2[2] << std::endl;
    else {
        for ( unsigned int i = 0; i < Input_Opt.PARAMETER_BACKG_SO2.size(); i++ )
            std::cout << Input_Opt.PARAMETER_BACKG_SO2[i] << " ";
        std::cout << std::endl;
    }

    /* ---- Emission indices ------------------------------ */
    std::cout << " Emission indices        :" << std::endl;
    std::cout << "  => NOx  [" << Input_Opt.PARAMETER_EI_NOX_UNIT << "]    : ";
    if ( Input_Opt.PARAMETER_EI_NOX_RANGE )
        std::cout << Input_Opt.PARAMETER_EI_NOX[0] << ":" << Input_Opt.PARAMETER_EI_NOX[1] << ":" << Input_Opt.PARAMETER_EI_NOX[2] << std::endl;
    else {
        for ( unsigned int i = 0; i < Input_Opt.PARAMETER_EI_NOX.size(); i++ )
            std::cout << Input_Opt.PARAMETER_EI_NOX[i] << " ";
        std::cout << std::endl;
    }
    std::cout << "  => CO   [" << Input_Opt.PARAMETER_EI_CO_UNIT << "]    : ";
    if ( Input_Opt.PARAMETER_EI_CO_RANGE )
        std::cout << Input_Opt.PARAMETER_EI_CO[0] << ":" << Input_Opt.PARAMETER_EI_CO[1] << ":" << Input_Opt.PARAMETER_EI_CO[2] << std::endl;
    else {
        for ( unsigned int i = 0; i < Input_Opt.PARAMETER_EI_CO.size(); i++ )
            std::cout << Input_Opt.PARAMETER_EI_CO[i] << " ";
        std::cout << std::endl;
    }
    std::cout << "  => UHC  [" << Input_Opt.PARAMETER_EI_UHC_UNIT << "]    : ";
    if ( Input_Opt.PARAMETER_EI_UHC_RANGE )
        std::cout << Input_Opt.PARAMETER_EI_UHC[0] << ":" << Input_Opt.PARAMETER_EI_UHC[1] << ":" << Input_Opt.PARAMETER_EI_UHC[2] << std::endl;
    else {
        for ( unsigned int i = 0; i < Input_Opt.PARAMETER_EI_UHC.size(); i++ )
            std::cout << Input_Opt.PARAMETER_EI_UHC[i] << " ";
        std::cout << std::endl;
    }
    std::cout << "  => SO2  [" << Input_Opt.PARAMETER_EI_SO2_UNIT << "]    : ";
    if ( Input_Opt.PARAMETER_EI_SO2_RANGE )
        std::cout << Input_Opt.PARAMETER_EI_SO2[0] << ":" << Input_Opt.PARAMETER_EI_SO2[1] << ":" << Input_Opt.PARAMETER_EI_SO2[2] << std::endl;
    else {
        for ( unsigned int i = 0; i < Input_Opt.PARAMETER_EI_SO2.size(); i++ )
            std::cout << Input_Opt.PARAMETER_EI_SO2[i] << " ";
        std::cout << std::endl;
    }
    std::cout << "  => SO2 to SO4 conv [" << Input_Opt.PARAMETER_EI_SO2TOSO4_UNIT << "] : ";
    if ( Input_Opt.PARAMETER_EI_SO2TOSO4_RANGE )
        std::cout << Input_Opt.PARAMETER_EI_SO2TOSO4[0] << ":" << Input_Opt.PARAMETER_EI_SO2TOSO4[1] << ":" << Input_Opt.PARAMETER_EI_SO2TOSO4[2] << std::endl;
    else {
        for ( unsigned int i = 0; i < Input_Opt.PARAMETER_EI_SO2TOSO4.size(); i++ )
            std::cout << Input_Opt.PARAMETER_EI_SO2TOSO4[i] << " ";
        std::cout << std::endl;
    }
    std::cout << "  => Soot [" << Input_Opt.PARAMETER_EI_SOOT_UNIT << "]    : ";
    if ( Input_Opt.PARAMETER_EI_SOOT_RANGE )
        std::cout << Input_Opt.PARAMETER_EI_SOOT[0] << ":" << Input_Opt.PARAMETER_EI_SOOT[1] << ":" << Input_Opt.PARAMETER_EI_SOOT[2] << std::endl;
    else {
        for ( unsigned int i = 0; i < Input_Opt.PARAMETER_EI_SOOT.size(); i++ )
            std::cout << Input_Opt.PARAMETER_EI_SOOT[i] << " ";
        std::cout << std::endl;
    }
    std::cout << "  => Soot Radius [" << Input_Opt.PARAMETER_EI_SOOTRAD_UNIT << "]     : ";
    if ( Input_Opt.PARAMETER_EI_SOOTRAD_RANGE )
        std::cout << Input_Opt.PARAMETER_EI_SOOTRAD[0] << ":" << Input_Opt.PARAMETER_EI_SOOTRAD[1] << ":" << Input_Opt.PARAMETER_EI_SOOTRAD[2] << std::endl;
    else {
        for ( unsigned int i = 0; i < Input_Opt.PARAMETER_EI_SOOTRAD.size(); i++ )
            std::cout << Input_Opt.PARAMETER_EI_SOOTRAD[i] << " ";
        std::cout << std::endl;
    }

    /* ---- Total fuel flow ------------------------------- */
    std::cout << "  Total fuel flow [" << Input_Opt.PARAMETER_FF_UNIT << "] : ";
    if ( Input_Opt.PARAMETER_FF_RANGE )
        std::cout << Input_Opt.PARAMETER_FF[0] << ":" << Input_Opt.PARAMETER_FF[1] << ":" << Input_Opt.PARAMETER_FF[2] << std::endl;
    else {
        for ( unsigned int i = 0; i < Input_Opt.PARAMETER_FF.size(); i++ )
            std::cout << Input_Opt.PARAMETER_FF[i] << " ";
        std::cout << std::endl;
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

    /* ==================================================== */
    /* Turn on Transport?                                   */
    /* ==================================================== */

    variable = "Turn on Transport?";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;
   
    /* Extract variable range */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    if ( ( strcmp(tokens[0].c_str(), "T" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" ) == 0 ) )
        Input_Opt.TRANSPORT_TRANSPORT = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" ) == 0 ) )
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

    if ( ( strcmp(tokens[0].c_str(), "T" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" ) == 0 ) )
        Input_Opt.TRANSPORT_FILL = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" ) == 0 ) )
        Input_Opt.TRANSPORT_FILL = 0;
    else {
        std::cout << " Wrong input for: " << variable << std::endl;
        exit(1);
    }
    
    /* ==================================================== */
    /* Convect timestep                                     */
    /* ==================================================== */

    variable = "Convect timestep";
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

    /* Return success */
    RC = SUCCESS;

    /* ==================================================== */
    /* Print to screen                                      */
    /* ==================================================== */

    std::cout << " ------------------------+---------------------------------------- "  << std::endl;
    std::cout << " %%% TRANSPORT MENU %%%  :"                                           << std::endl;
    std::cout << " ------------------------+---------------------------------------- "  << std::endl;
    std::cout << " Turn on Transport?      : " << Input_Opt.TRANSPORT_TRANSPORT         << std::endl;
    std::cout << "  => Fill Negative Values: " << Input_Opt.TRANSPORT_FILL              << std::endl;
    std::cout << " Convect Timestep [min]  : " << Input_Opt.TRANSPORT_TIMESTEP          << std::endl;

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

    if ( ( strcmp(tokens[0].c_str(), "T" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" ) == 0 ) )
        Input_Opt.CHEMISTRY_CHEMISTRY = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" ) == 0 ) )
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
    /* Use ring structure?                                  */
    /* ==================================================== */

    variable = "Use ring structure?";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;
   
    /* Extract variable range */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    if ( ( strcmp(tokens[0].c_str(), "T" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" ) == 0 ) )
        Input_Opt.CHEMISTRY_RINGS = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" ) == 0 ) )
        Input_Opt.CHEMISTRY_RINGS = 0;
    else {
        std::cout << " Wrong input for: " << variable << std::endl;
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

    if ( ( strcmp(tokens[0].c_str(), "T" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" ) == 0 ) )
        Input_Opt.CHEMISTRY_HETCHEM = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" ) == 0 ) )
        Input_Opt.CHEMISTRY_HETCHEM = 0;
    else {
        std::cout << " Wrong input for: " << variable << std::endl;
        exit(1);
    }
    
    /* ==================================================== */
    /* Read in J-Rates?                                     */
    /* ==================================================== */

    variable = "Read in J-Rates?";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;
    
    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );
    
    if ( ( strcmp(tokens[0].c_str(), "T" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" ) == 0 ) )
        Input_Opt.CHEMISTRY_READ_JRATES = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" ) == 0 ) )
        Input_Opt.CHEMISTRY_READ_JRATES = 0;
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

    /* Return success */
    RC = SUCCESS;

    /* ==================================================== */
    /* Print to screen                                      */
    /* ==================================================== */

    std::cout << " ------------------------+---------------------------------------- "  << std::endl;
    std::cout << " %%% CHEMISTRY MENU %%%  :"                                           << std::endl;
    std::cout << " ------------------------+---------------------------------------- "  << std::endl;
    std::cout << " Turn on Chemistry?      : " << Input_Opt.CHEMISTRY_CHEMISTRY         << std::endl;
    std::cout << "  => Use ring structure? : " << Input_Opt.CHEMISTRY_RINGS             << std::endl;
    std::cout << " Perform het. chem.?     : " << Input_Opt.CHEMISTRY_HETCHEM           << std::endl;
    std::cout << " Read in J-Rates?        : " << Input_Opt.CHEMISTRY_READ_JRATES       << std::endl;
    std::cout << " Chemistry Timestep [min]: " << Input_Opt.CHEMISTRY_TIMESTEP          << std::endl;

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

    if ( ( strcmp(tokens[0].c_str(), "T" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" ) == 0 ) )
        Input_Opt.AEROSOL_GRAVSETTLING = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" ) == 0 ) )
        Input_Opt.AEROSOL_GRAVSETTLING = 0;
    else {
        std::cout << " Wrong input for: " << variable << std::endl;
        exit(1);
    }

    /* ==================================================== */
    /* Turn on coagulation?                                 */
    /* ==================================================== */

    variable = "Turn on coagulation?";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;
   
    /* Extract variable range */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );

    if ( ( strcmp(tokens[0].c_str(), "T" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" ) == 0 ) )
        Input_Opt.AEROSOL_COAGULATION = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" ) == 0 ) )
        Input_Opt.AEROSOL_COAGULATION = 0;
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
    
    if ( ( strcmp(tokens[0].c_str(), "T" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" ) == 0 ) )
        Input_Opt.AEROSOL_ICE_GROWTH = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" ) == 0 ) )
        Input_Opt.AEROSOL_ICE_GROWTH = 0;
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
    
    if ( ( strcmp(tokens[0].c_str(), "T" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" ) == 0 ) )
        Input_Opt.AEROSOL_PLUME_UPDRAFT = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" ) == 0 ) )
        Input_Opt.AEROSOL_PLUME_UPDRAFT = 0;
    else {
        std::cout << " Wrong input for: " << variable << std::endl;
        exit(1);
    }
    
    /* Return success */
    RC = SUCCESS;

    /* ==================================================== */
    /* Print to screen                                      */
    /* ==================================================== */

    std::cout << " ------------------------+---------------------------------------- "  << std::endl;
    std::cout << " %%% AEROSOL MENU %%%    :"                                           << std::endl;
    std::cout << " ------------------------+---------------------------------------- "  << std::endl;
    std::cout << " Turn on grav. settling? : " << Input_Opt.AEROSOL_GRAVSETTLING        << std::endl;
    std::cout << " Turn on coagulation?    : " << Input_Opt.AEROSOL_COAGULATION         << std::endl;
    std::cout << "  => Coag. timestep [min]: " << Input_Opt.AEROSOL_COAGULATION_TIMESTEP<< std::endl;
    std::cout << " Turn on ice growth?     : " << Input_Opt.AEROSOL_ICE_GROWTH          << std::endl;
    std::cout << " Turn on plume updraft?  : " << Input_Opt.AEROSOL_PLUME_UPDRAFT       << std::endl;

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

    /* ==================================================== */
    /* Initialize T from MET?                               */
    /* ==================================================== */
    
    variable = "Initialize T from MET?";
    getline( inputFile, line, '\n' );
    if ( VERBOSE )
        std::cout << line << std::endl;
    
    /* Extract variable */
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );
    
    if ( ( strcmp(tokens[0].c_str(), "T" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" ) == 0 ) )
        Input_Opt.MET_TEMP_INIT = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" ) == 0 ) )
        Input_Opt.MET_TEMP_INIT = 0;
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
    
    if ( ( strcmp(tokens[0].c_str(), "T" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" ) == 0 ) )
        Input_Opt.MET_H2O_INIT = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" ) == 0 ) )
        Input_Opt.MET_H2O_INIT = 0;
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

    /* Return success */
    RC = SUCCESS;

    /* ==================================================== */
    /* Print to screen                                      */
    /* ==================================================== */

    std::cout << " ------------------------+---------------------------------------- "  << std::endl;
    std::cout << " %%% METEOROLOGY MENU %%%:"                                           << std::endl;
    std::cout << " ------------------------+---------------------------------------- "  << std::endl;
    std::cout << " Initialize T from MET?  : " << Input_Opt.MET_TEMP_INIT               << std::endl;
    std::cout << " Initialize H2O from MET?: " << Input_Opt.MET_H2O_INIT                << std::endl;
    std::cout << " Met file                : " << Input_Opt.MET_FILENAME                << std::endl;

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

    std::cout << " ------------------------+---------------------------------------- "  << std::endl;
    std::cout << " %%% DIAGNOSTIC MENU %%% :"                                           << std::endl;
    std::cout << " ------------------------+---------------------------------------- "  << std::endl;
    std::cout << " netCDF file name        : " << Input_Opt.DIAG_FILENAME               << std::endl;
    std::cout << " Diagnostic Entries ---> : L" << std::endl;

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
    
    if ( ( strcmp(tokens[0].c_str(), "T" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" ) == 0 ) )
        Input_Opt.TS_SPEC = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" ) == 0 ) )
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
            if ( value > 0.0E+00 )
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
    
    if ( ( strcmp(tokens[0].c_str(), "T" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" ) == 0 ) )
        Input_Opt.TS_AERO = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" ) == 0 ) )
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
            if ( value > 0.0E+00 )
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

    std::cout << " ------------------------+---------------------------------------- "  << std::endl;
    std::cout << " %%% TIMESERIES MENU %%% :"                                           << std::endl;
    std::cout << " ------------------------+---------------------------------------- "  << std::endl;
    std::cout << " Save species timeseries?: " << Input_Opt.TS_SPEC                     << std::endl;
    std::cout << "  => Inst timeseries file: " << Input_Opt.TS_FILENAME                 << std::endl;
    std::cout << "  => Species to include  : ";
    for ( unsigned int i = 0; i < Input_Opt.TS_SPECIES.size(); i++ )
        std::cout << Input_Opt.TS_SPECIES[i] << " ";
    std::cout << std::endl;
    std::cout << "  => Frequency [min]     : " << Input_Opt.TS_FREQ                     << std::endl;
    std::cout << " Save aerosol timeseries?: " << Input_Opt.TS_AERO                     << std::endl;
    std::cout << "  => Inst timeseries file: " << Input_Opt.TS_AERO_FILENAME                 << std::endl;
    std::cout << "  => Aerosol to include  : ";
    for ( unsigned int i = 0; i < Input_Opt.TS_AEROSOL.size(); i++ )
        std::cout << Input_Opt.TS_AEROSOL[i] << " ";
    std::cout << std::endl;
    std::cout << "  => Frequency [min]     : " << Input_Opt.TS_AERO_FREQ                     << std::endl;

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

    if ( ( strcmp(tokens[0].c_str(), "T" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" ) == 0 ) )
        Input_Opt.PL_PL = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" ) == 0 ) )
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

    if ( ( strcmp(tokens[0].c_str(), "T" ) == 0 ) || \
         ( strcmp(tokens[0].c_str(), "1" ) == 0 ) )
        Input_Opt.PL_O3 = 1;
    else if ( ( strcmp(tokens[0].c_str(), "F" ) == 0 ) || \
              ( strcmp(tokens[0].c_str(), "0" ) == 0 ) )
        Input_Opt.PL_O3 = 0;
    else {
        std::cout << " Wrong input for: " << variable << std::endl;
        exit(1);
    }

    /* Return success */
    RC = SUCCESS;
    
    /* ==================================================== */
    /* Print to screen                                      */
    /* ==================================================== */

    std::cout << " ------------------------+---------------------------------------- "  << std::endl;
    std::cout << " %%% PROD & LOSS MENU %%%:"                                           << std::endl;
    std::cout << " ------------------------+---------------------------------------- "  << std::endl;
    std::cout << " Turn on P/L diag?       : " << Input_Opt.PL_PL                       << std::endl;
    std::cout << " Save O3 P/L?            : " << Input_Opt.PL_O3                       << std::endl;

} /* End of Read_PL_Menu */

Vector_2D CombVec( OptInput &Input_Opt )
{

    unsigned int counter = 1;
    unsigned int nCases;  

    unsigned int i, j;

    double currVal = 0.0E+00;

    Vector_1D cases;
    Vector_2D y, z;
    Vector_2D u, v;

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
    nCases = cases.size();

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

    y.push_back(Vector_1D(cases.size()));
    for ( i = 0; i < cases.size(); i++ )
        y[0][i] = cases[i];
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
    Input_Opt.PARAMETER_TEMPERATURE_UNIT = "K";

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
    Input_Opt.PARAMETER_TEMPERATURE_UNIT = "K";

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
    Input_Opt.PARAMETER_TEMPERATURE_UNIT = "K";
    
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
        for ( i = 0; i < nCases; i++ )
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
    
    if ( Input_Opt.PARAMETER_EDAY_UNIT.compare( "1-365" ) == 0 ) {
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
        currVal = Input_Opt.PARAMETER_ETIME[0];
        while ( currVal <= Input_Opt.PARAMETER_ETIME[2] ) {
            cases.push_back( currVal );
            currVal += Input_Opt.PARAMETER_ETIME[1];
        }
    } else {
        for ( i = 0; i < Input_Opt.PARAMETER_ETIME.size(); i++ )
            cases.push_back(Input_Opt.PARAMETER_ETIME[i]);
    }
    
    if ( Input_Opt.PARAMETER_ETIME_UNIT.compare( "0-24" ) == 0 ) {
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
    
    if ( Input_Opt.PARAMETER_EI_NOX_UNIT.compare( "g/kg_fuel" ) == 0 ) {
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
    
    if ( Input_Opt.PARAMETER_EI_CO_UNIT.compare( "g/kg_fuel" ) == 0 ) {
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
    
    if ( Input_Opt.PARAMETER_EI_UHC_UNIT.compare( "g/kg_fuel" ) == 0 ) {
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
    
    if ( Input_Opt.PARAMETER_EI_SO2_UNIT.compare( "g/kg_fuel" ) == 0 ) {
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
    
    if ( Input_Opt.PARAMETER_EI_SOOT_UNIT.compare( "g/kg_fuel" ) == 0 ) {
        /* Do nothing. Default unit */
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
    
    if ( Input_Opt.PARAMETER_FF_UNIT.compare( "kg/s" ) == 0 ) {
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
    
    return y;

} /* End of CombVec */

Vector_2D Copy_blocked( Vector_2D& m, int n )
{
    
    Vector_2D b;
    const unsigned int mr = m.size();
    const unsigned int mc = m[0].size();
    unsigned int i, j, k;

    for ( i = 0; i < mr; i++ )
        b.push_back(Vector_1D(mc*n));

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

    for ( i = 0; i < mr*n; i++ )
        b.push_back(Vector_1D(mc));

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
            //counter_y = counter%n_y;
            orig_counter_x = counter%size_x;
            //orig_counter_y = counter%size_y;
            output[counter_x][counter_y] = vector_2D[orig_counter_x][orig_counter_y];
            //std::cout << orig_counter_x << ", " << orig_counter_y << ": " << vector_2D[orig_counter_x][orig_counter_y] << " -> " << counter_x << ", " << counter_y << std::endl;
            if ( counter_x == n_x - 1 )
                counter_y += 1;
            if ( orig_counter_x == size_x - 1 )
                orig_counter_y += 1;
//            if ( counter_y == n_y - 1 )
//                counter_x += 1;
//            if ( orig_counter_y == size_y - 1 )
//                orig_counter_x += 1;
        }
//        int k, l;
//        for ( i = 0; i < n_x; i++ ) {
//            for ( j = 0; j < n_y; j++ ) {
//                k = i%size_x;
//                l = j%size_y;
//                output[i][j] = vector_2D[k][l];
//            }
//        }
    }

    return output;

} /* End of Reshape_Vector */


/* End of Input_Reader.cpp */
