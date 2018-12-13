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
        std::cout << " Make sure that the variable 'currDir' is exported" << std::endl;
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

    if ( strcmp(tokens[0].c_str(), "T" ) == 0 )
        Input_Opt.SIMULATION_PARAMETER_SWEEP = 1;
    else if ( strcmp(tokens[0].c_str(), "F" ) == 0 )
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

    if ( strcmp(tokens[0].c_str(), "T" ) == 0 )
        Input_Opt.SIMULATION_OVERWRITE = 1;
    else if ( strcmp(tokens[0].c_str(), "F" ) == 0 )
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

    if ( strcmp(tokens[0].c_str(), "T" ) == 0 )
        Input_Opt.SIMULATION_SAVE_FORWARD = 0;
    else if ( strcmp(tokens[0].c_str(), "F" ) == 0 )
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

    if ( strcmp(tokens[0].c_str(), "T" ) == 0 )
        Input_Opt.SIMULATION_ADJOINT = 0;
    else if ( strcmp(tokens[0].c_str(), "F" ) == 0 )
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
    for ( unsigned int i = 0; i < tokens.size(); i++ )
        Input_Opt.PARAMETER_TEMPERATURE.push_back(std::stod(tokens[i]));
    
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
    for ( unsigned int i = 0; i < tokens.size(); i++ )
        Input_Opt.PARAMETER_RHW.push_back(std::stod(tokens[i]));


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
    for ( unsigned int i = 0; i < tokens.size(); i++ )
        Input_Opt.PARAMETER_LONGITUDE.push_back(std::stod(tokens[i]));

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
    for ( unsigned int i = 0; i < tokens.size(); i++ )
        Input_Opt.PARAMETER_LATITUDE.push_back(std::stod(tokens[i]));
    
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
    for ( unsigned int i = 0; i < tokens.size(); i++ )
        Input_Opt.PARAMETER_PRESSURE.push_back(std::stod(tokens[i]));
    
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
    for ( unsigned int i = 0; i < tokens.size(); i++ )
        Input_Opt.PARAMETER_EDAY.push_back(std::stoi(tokens[i]) % 365);

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
    for ( unsigned int i = 0; i < tokens.size(); i++ )
        Input_Opt.PARAMETER_ETIME.push_back(std::fmod(std::stod(tokens[i]),2.40E+01));
    
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
    for ( unsigned int i = 0; i < tokens.size(); i++ )
        Input_Opt.PARAMETER_BACKG_NOX.push_back(std::stod(tokens[i]));

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
    for ( unsigned int i = 0; i < tokens.size(); i++ )
        Input_Opt.PARAMETER_BACKG_HNO3.push_back(std::stod(tokens[i]));
    
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
    for ( unsigned int i = 0; i < tokens.size(); i++ )
        Input_Opt.PARAMETER_BACKG_O3.push_back(std::stod(tokens[i]));

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
    for ( unsigned int i = 0; i < tokens.size(); i++ )
        Input_Opt.PARAMETER_BACKG_CO.push_back(std::stod(tokens[i]));

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
    for ( unsigned int i = 0; i < tokens.size(); i++ )
        Input_Opt.PARAMETER_BACKG_CH4.push_back(std::stod(tokens[i]));

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
    for ( unsigned int i = 0; i < tokens.size(); i++ )
        Input_Opt.PARAMETER_BACKG_SO2.push_back(std::stod(tokens[i]));

    
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
    for ( unsigned int i = 0; i < tokens.size(); i++ )
        Input_Opt.PARAMETER_EI_NOX.push_back(std::stod(tokens[i]));

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
    for ( unsigned int i = 0; i < tokens.size(); i++ )
        Input_Opt.PARAMETER_EI_CO.push_back(std::stod(tokens[i]));

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
    for ( unsigned int i = 0; i < tokens.size(); i++ )
        Input_Opt.PARAMETER_EI_UHC.push_back(std::stod(tokens[i]));

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
    for ( unsigned int i = 0; i < tokens.size(); i++ )
        Input_Opt.PARAMETER_EI_SO2.push_back(std::stod(tokens[i]));

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
    for ( unsigned int i = 0; i < tokens.size(); i++ )
        Input_Opt.PARAMETER_EI_SO2TOSO4.push_back(std::stod(tokens[i]));

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
    for ( unsigned int i = 0; i < tokens.size(); i++ )
        Input_Opt.PARAMETER_EI_SOOT.push_back(std::stod(tokens[i]));

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
    for ( unsigned int i = 0; i < tokens.size(); i++ )
        Input_Opt.PARAMETER_EI_SOOTRAD.push_back(std::stod(tokens[i]));

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
    for ( unsigned int i = 0; i < tokens.size(); i++ )
        Input_Opt.PARAMETER_FF.push_back(std::stod(tokens[i]));


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

    if ( strcmp(tokens[0].c_str(), "T" ) == 0 )
        Input_Opt.TRANSPORT_TRANSPORT = 1;
    else if ( strcmp(tokens[0].c_str(), "F" ) == 0 )
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

    if ( strcmp(tokens[0].c_str(), "T" ) == 0 )
        Input_Opt.TRANSPORT_FILL = 1;
    else if ( strcmp(tokens[0].c_str(), "F" ) == 0 )
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
        if ( value > 0.0E+00 )
            Input_Opt.TRANSPORT_TIMESTEP = value;
        else {
            std::cout << " Wrong input for: " << variable << std::endl;
            std::cout << " Timestep needs to be positive" << std::endl;
            exit(1);
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

    if ( strcmp(tokens[0].c_str(), "T" ) == 0 )
        Input_Opt.CHEMISTRY_CHEMISTRY = 1;
    else if ( strcmp(tokens[0].c_str(), "F" ) == 0 )
        Input_Opt.CHEMISTRY_CHEMISTRY = 0;
    else {
        std::cout << " Wrong input for: " << variable << std::endl;
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

    if ( strcmp(tokens[0].c_str(), "T" ) == 0 )
        Input_Opt.CHEMISTRY_RINGS = 1;
    else if ( strcmp(tokens[0].c_str(), "F" ) == 0 )
        Input_Opt.CHEMISTRY_RINGS = 0;
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
    
    if ( strcmp(tokens[0].c_str(), "T" ) == 0 )
        Input_Opt.CHEMISTRY_READ_JRATES = 1;
    else if ( strcmp(tokens[0].c_str(), "F" ) == 0 )
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
        if ( value > 0.0E+00 )
            Input_Opt.CHEMISTRY_TIMESTEP = value;
        else {
            std::cout << " Wrong input for: " << variable << std::endl;
            std::cout << " Timestep needs to be positive" << std::endl;
            exit(1);
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

    if ( strcmp(tokens[0].c_str(), "T" ) == 0 )
        Input_Opt.AEROSOL_GRAVSETTLING = 1;
    else if ( strcmp(tokens[0].c_str(), "F" ) == 0 )
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

    if ( strcmp(tokens[0].c_str(), "T" ) == 0 )
        Input_Opt.AEROSOL_COAGULATION = 1;
    else if ( strcmp(tokens[0].c_str(), "F" ) == 0 )
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
    
    if ( strcmp(tokens[0].c_str(), "T" ) == 0 )
        Input_Opt.AEROSOL_ICE_GROWTH = 1;
    else if ( strcmp(tokens[0].c_str(), "F" ) == 0 )
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
    
    if ( strcmp(tokens[0].c_str(), "T" ) == 0 )
        Input_Opt.AEROSOL_PLUME_UPDRAFT = 1;
    else if ( strcmp(tokens[0].c_str(), "F" ) == 0 )
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

Vector_2D ReadParameters( const OptInput Input_Opt )
{

    Vector_1D temperature_K, pressure_Pa, relHumidity_w, longitude_deg, latitude_deg;
    Vector_1D EI_NOx, EI_CO, EI_HC, EI_SO2, EI_Soot, SootRad, ff;
    Vector_1D dayGMT;
    Vector_1D emissionTime;
    Vector_1D backgNOx, backgHNO3, backgO3, backgCO, backgCH4, backgSO2;

    Vector_2D parameters;

    /* #######################################
     * ## Ambient meteorological conditions ## 
     * ####################################### */

    /* Temperature array in [K] */

    temperature_K.push_back( 220.0 );
//    for ( unsigned int i = 0; i < 31; i++ ) {
//        temperature_K.push_back( 200.0 + 2.0 * i );
//    }

    /* Pressure array in [Pa] */

    pressure_Pa.push_back( 220.0E2 );

    /* Relative humidity w.r.t liquid water array in [\%] */

    relHumidity_w.push_back( 30.0 );

    /* #######################
     * ## Coord. (lon-lat)  ## 
     * ####################### */

    /* Longitude array expressed in deg */

    longitude_deg.push_back( -15.0 );

    /* Latitude array expressed in deg */

    latitude_deg.push_back( 60.0 );

    /* #######################
     * ## Emiss day (0-365) ## 
     * ####################### */

    dayGMT.push_back( 81.0 );

    /* #######################
     * ## Emiss time  [hrs] ## 
     * ####################### */

    emissionTime.push_back( 8.0 );
    emissionTime.push_back( 20.0 );

    /* ####################### 
     * ## EMISSION INDICES  ##
     * ####################### */

    EI_NOx.push_back( 1.0E-05 );
    EI_NOx.push_back( 1.0E+01 );

//    EI_CO.push_back( 0.5E+00 );
    EI_CO.push_back( 1.0E+00 );
//    EI_CO.push_back( 1.5E+00 );
//    EI_CO.push_back( 2.0E+00 );

    EI_HC.push_back( 0.6E+00 );
    
    EI_SO2.push_back( 0.8E+00 );
    EI_SO2.push_back( 1.2E+00 );
    EI_SO2.push_back( 2.0E+00 );
    EI_SO2.push_back( 1.0E+01 );
    EI_SO2.push_back( 2.0E+01 );
    EI_SO2.push_back( 1.0E+02 );
    
    EI_Soot.push_back( 0.0E+00 );
    
    SootRad.push_back( 0.0E+00 );
    
    ff.push_back( 0.0E+00 );

    /* #######################
     * ## BACKGRD MIX RATIO ##
     * ####################### */

    backgNOx.push_back( 50.0E-03 );
//    backgNOx.push_back( 100.0E-03 );
//    backgNOx.push_back( 150.0E-03 );
//    backgNOx.push_back( 200.0E-03 );
    
    backgHNO3.push_back( 81.5E-03 );

    backgO3.push_back( 50.0E+00 );
//    backgO3.push_back( 75.0E+00 );
//    backgO3.push_back( 100.0E+00 );

//    backgCO.push_back( 20.0E+00 );
    backgCO.push_back( 40.0E+00 );
//    backgCO.push_back( 60.0E+00 );
//    backgCO.push_back( 80.0E+00 );
    
    backgCH4.push_back( 1.76E+03 );
    backgSO2.push_back( 7.25E-03 );

    parameters = CombVec( temperature_K, \
                          pressure_Pa,   \
                          relHumidity_w, \
                          longitude_deg, \
                          latitude_deg,  \
                          dayGMT,        \
                          emissionTime,  \
                          EI_NOx,        \
                          EI_CO,         \
                          EI_HC,         \
                          EI_SO2,        \
                          EI_Soot,       \
                          SootRad,       \
                          ff,            \
                          backgNOx,      \
                          backgHNO3,     \
                          backgO3,       \
                          backgCO,       \
                          backgCH4,      \
                          backgSO2 );

    /* For debug */
    if ( false ) {
        unsigned int i, j;
        std::cout.precision(2);
        for ( i = 0; i < parameters.size(); i++ ) {
            for ( j = 0; j < parameters[0].size(); j++) {
                std::cout << parameters[i][j] << " ";
            }
            std::cout << "" << std::endl;
        }
    }

    return parameters;

} /* End of ReadParameters */

Vector_2D CombVec( const Vector_1D& temperature_K, \
                   const Vector_1D& pressure_Pa,   \
                   const Vector_1D& relHumidity_w, \
                   const Vector_1D& longitude_deg, \
                   const Vector_1D& latitude_deg,  \
                   const Vector_1D& dayGMT,        \
                   const Vector_1D& emissionTime,  \
                   const Vector_1D& EI_NOx,        \
                   const Vector_1D& EI_CO,         \
                   const Vector_1D& EI_HC,         \
                   const Vector_1D& EI_SO2,        \
                   const Vector_1D& EI_Soot,       \
                   const Vector_1D& SootRad,       \
                   const Vector_1D& ff,            \
                   const Vector_1D& backgNOx,      \
                   const Vector_1D& backgHNO3,     \
                   const Vector_1D& backgO3,       \
                   const Vector_1D& backgCO,       \
                   const Vector_1D& backgCH4,      \
                   const Vector_1D& backgSO2 )
{
    Vector_2D combinations;

    unsigned int counter = 1;
    unsigned int nCases = temperature_K.size();

    unsigned int i, j;

    Vector_2D y, z;
    Vector_2D u, v;

    y.push_back(Vector_1D(temperature_K.size()));
    for ( i = 0; i < temperature_K.size(); i++ )
        y[0][i] = temperature_K[i];

    /* z = Pressure_Pa */
    z.push_back(Vector_1D(pressure_Pa.size()));
    for ( i = 0; i < pressure_Pa.size(); i++ )
        z[0][i] = pressure_Pa[i];

    u = Copy_blocked(y,z[0].size());
    v = Copy_interleaved(z,y[0].size());

    for ( i = 0; i < counter; i++ ) {
        y[i].clear();
    }
    z[0].clear();
    y.clear(); z.clear();

    nCases *= pressure_Pa.size();
    counter += 1;
    for ( i = 0; i < counter; i++ )
        y.push_back(Vector_1D( nCases ));

    for ( i = 0; i < nCases; i++ ) {
        for ( j = 0; j < counter - 1; j++ ) {
            y[j][i] = u[j][i];
        }
        y[counter-1][i] = v[0][i];
    }

    for ( i = 0; i < counter - 1; i++ ) {
        u[i].clear();
    }
    u.clear();
    v[0].clear(); v.clear();

    /* z = relHumidity_w */
    z.push_back(Vector_1D(relHumidity_w.size()));
    for ( i = 0; i < relHumidity_w.size(); i++ )
        z[0][i] = relHumidity_w[i];

    u = Copy_blocked(y,z[0].size());
    v = Copy_interleaved(z,y[0].size());

    for ( i = 0; i < counter; i++ ) {
        y[i].clear();
    }
    z[0].clear();
    y.clear(); z.clear();

    nCases *= relHumidity_w.size();
    counter += 1;
    for ( i = 0; i < counter; i++ )
        y.push_back(Vector_1D( nCases ));

    for ( i = 0; i < nCases; i++ ) {
        for ( j = 0; j < counter - 1; j++ ) {
            y[j][i] = u[j][i];
        }
        y[counter-1][i] = v[0][i];
    }

    for ( i = 0; i < counter - 1; i++ ) {
        u[i].clear();
    }
    u.clear();
    v[0].clear(); v.clear();

    /* z = longitude_deg */
    z.push_back(Vector_1D(longitude_deg.size()));
    for ( i = 0; i < longitude_deg.size(); i++ )
        z[0][i] = longitude_deg[i];

    u = Copy_blocked(y,z[0].size());
    v = Copy_interleaved(z,y[0].size());

    for ( i = 0; i < counter; i++ ) {
        y[i].clear();
    }
    z[0].clear();
    y.clear(); z.clear();

    nCases *= longitude_deg.size();
    counter += 1;
    for ( i = 0; i < counter; i++ )
        y.push_back(Vector_1D( nCases ));

    for ( i = 0; i < nCases; i++ ) {
        for ( j = 0; j < counter - 1; j++ ) {
            y[j][i] = u[j][i];
        }
        y[counter-1][i] = v[0][i];
    }

    for ( i = 0; i < counter - 1; i++ ) {
        u[i].clear();
    }
    u.clear();
    v[0].clear(); v.clear();

    /* z = latitude_deg */
    z.push_back(Vector_1D(latitude_deg.size()));
    for ( i = 0; i < latitude_deg.size(); i++ )
        z[0][i] = latitude_deg[i];

    u = Copy_blocked(y,z[0].size());
    v = Copy_interleaved(z,y[0].size());

    for ( i = 0; i < counter; i++ ) {
        y[i].clear();
    }
    z[0].clear();
    y.clear(); z.clear();

    nCases *= latitude_deg.size();
    counter += 1;
    for ( i = 0; i < counter; i++ )
        y.push_back(Vector_1D( nCases ));

    for ( i = 0; i < nCases; i++ ) {
        for ( j = 0; j < counter - 1; j++ ) {
            y[j][i] = u[j][i];
        }
        y[counter-1][i] = v[0][i];
    }
    
    /* z = dayGMT */
    z.push_back(Vector_1D(dayGMT.size()));
    for ( i = 0; i < dayGMT.size(); i++ )
        z[0][i] = dayGMT[i];

    u = Copy_blocked(y,z[0].size());
    v = Copy_interleaved(z,y[0].size());

    for ( i = 0; i < counter; i++ ) {
        y[i].clear();
    }
    z[0].clear();
    y.clear(); z.clear();

    nCases *= dayGMT.size();
    counter += 1;
    for ( i = 0; i < counter; i++ )
        y.push_back(Vector_1D( nCases ));

    for ( i = 0; i < nCases; i++ ) {
        for ( j = 0; j < counter - 1; j++ ) {
            y[j][i] = u[j][i];
        }
        y[counter-1][i] = v[0][i];
    }
    
    /* z = emissionTime */
    z.push_back(Vector_1D(emissionTime.size()));
    for ( i = 0; i < emissionTime.size(); i++ )
        z[0][i] = emissionTime[i];

    u = Copy_blocked(y,z[0].size());
    v = Copy_interleaved(z,y[0].size());

    for ( i = 0; i < counter; i++ ) {
        y[i].clear();
    }
    z[0].clear();
    y.clear(); z.clear();

    nCases *= emissionTime.size();
    counter += 1;
    for ( i = 0; i < counter; i++ )
        y.push_back(Vector_1D( nCases ));

    for ( i = 0; i < nCases; i++ ) {
        for ( j = 0; j < counter - 1; j++ ) {
            y[j][i] = u[j][i];
        }
        y[counter-1][i] = v[0][i];
    }

    /* z = EI_NOx */
    z.push_back(Vector_1D(EI_NOx.size()));
    for ( i = 0; i < EI_NOx.size(); i++ )
        z[0][i] = EI_NOx[i];

    u = Copy_blocked(y,z[0].size());
    v = Copy_interleaved(z,y[0].size());

    for ( i = 0; i < counter; i++ ) {
        y[i].clear();
    }
    z[0].clear();
    y.clear(); z.clear();

    nCases *= EI_NOx.size();
    counter += 1;
    for ( i = 0; i < counter; i++ )
        y.push_back(Vector_1D( nCases ));

    for ( i = 0; i < nCases; i++ ) {
        for ( j = 0; j < counter - 1; j++ ) {
            y[j][i] = u[j][i];
        }
        y[counter-1][i] = v[0][i];
    }

    /* z = EI_CO */
    z.push_back(Vector_1D(EI_CO.size()));
    for ( i = 0; i < EI_CO.size(); i++ )
        z[0][i] = EI_CO[i];

    u = Copy_blocked(y,z[0].size());
    v = Copy_interleaved(z,y[0].size());

    for ( i = 0; i < counter; i++ ) {
        y[i].clear();
    }
    z[0].clear();
    y.clear(); z.clear();

    nCases *= EI_CO.size();
    counter += 1;
    for ( i = 0; i < counter; i++ )
        y.push_back(Vector_1D( nCases ));

    for ( i = 0; i < nCases; i++ ) {
        for ( j = 0; j < counter - 1; j++ ) {
            y[j][i] = u[j][i];
        }
        y[counter-1][i] = v[0][i];
    }

    /* z = EI_HC */
    z.push_back(Vector_1D(EI_HC.size()));
    for ( i = 0; i < EI_HC.size(); i++ )
        z[0][i] = EI_HC[i];

    u = Copy_blocked(y,z[0].size());
    v = Copy_interleaved(z,y[0].size());

    for ( i = 0; i < counter; i++ ) {
        y[i].clear();
    }
    z[0].clear();
    y.clear(); z.clear();

    nCases *= EI_HC.size();
    counter += 1;
    for ( i = 0; i < counter; i++ )
        y.push_back(Vector_1D( nCases ));

    for ( i = 0; i < nCases; i++ ) {
        for ( j = 0; j < counter - 1; j++ ) {
            y[j][i] = u[j][i];
        }
        y[counter-1][i] = v[0][i];
    }

    /* z = EI_SO2 */
    z.push_back(Vector_1D(EI_SO2.size()));
    for ( i = 0; i < EI_SO2.size(); i++ )
        z[0][i] = EI_SO2[i];

    u = Copy_blocked(y,z[0].size());
    v = Copy_interleaved(z,y[0].size());

    for ( i = 0; i < counter; i++ ) {
        y[i].clear();
    }
    z[0].clear();
    y.clear(); z.clear();

    nCases *= EI_SO2.size();
    counter += 1;
    for ( i = 0; i < counter; i++ )
        y.push_back(Vector_1D( nCases ));

    for ( i = 0; i < nCases; i++ ) {
        for ( j = 0; j < counter - 1; j++ ) {
            y[j][i] = u[j][i];
        }
        y[counter-1][i] = v[0][i];
    }

    /* z = EI_Soot */
    z.push_back(Vector_1D(EI_Soot.size()));
    for ( i = 0; i < EI_Soot.size(); i++ )
        z[0][i] = EI_Soot[i];

    u = Copy_blocked(y,z[0].size());
    v = Copy_interleaved(z,y[0].size());

    for ( i = 0; i < counter; i++ ) {
        y[i].clear();
    }
    z[0].clear();
    y.clear(); z.clear();

    nCases *= EI_Soot.size();
    counter += 1;
    for ( i = 0; i < counter; i++ )
        y.push_back(Vector_1D( nCases ));

    for ( i = 0; i < nCases; i++ ) {
        for ( j = 0; j < counter - 1; j++ ) {
            y[j][i] = u[j][i];
        }
        y[counter-1][i] = v[0][i];
    }

    /* z = SootRad */
    z.push_back(Vector_1D(SootRad.size()));
    for ( i = 0; i < SootRad.size(); i++ )
        z[0][i] = SootRad[i];

    u = Copy_blocked(y,z[0].size());
    v = Copy_interleaved(z,y[0].size());

    for ( i = 0; i < counter; i++ ) {
        y[i].clear();
    }
    z[0].clear();
    y.clear(); z.clear();

    nCases *= SootRad.size();
    counter += 1;
    for ( i = 0; i < counter; i++ )
        y.push_back(Vector_1D( nCases ));

    for ( i = 0; i < nCases; i++ ) {
        for ( j = 0; j < counter - 1; j++ ) {
            y[j][i] = u[j][i];
        }
        y[counter-1][i] = v[0][i];
    }

    /* z = ff */
    z.push_back(Vector_1D(ff.size()));
    for ( i = 0; i < ff.size(); i++ )
        z[0][i] = ff[i];

    u = Copy_blocked(y,z[0].size());
    v = Copy_interleaved(z,y[0].size());

    for ( i = 0; i < counter; i++ ) {
        y[i].clear();
    }
    z[0].clear();
    y.clear(); z.clear();

    nCases *= ff.size();
    counter += 1;
    for ( i = 0; i < counter; i++ )
        y.push_back(Vector_1D( nCases ));

    for ( i = 0; i < nCases; i++ ) {
        for ( j = 0; j < counter - 1; j++ ) {
            y[j][i] = u[j][i];
        }
        y[counter-1][i] = v[0][i];
    }

    /* z = backgNOx */
    z.push_back(Vector_1D(backgNOx.size()));
    for ( i = 0; i < backgNOx.size(); i++ )
        z[0][i] = backgNOx[i];

    u = Copy_blocked(y,z[0].size());
    v = Copy_interleaved(z,y[0].size());

    for ( i = 0; i < counter; i++ ) {
        y[i].clear();
    }
    z[0].clear();
    y.clear(); z.clear();

    nCases *= backgNOx.size();
    counter += 1;
    for ( i = 0; i < counter; i++ )
        y.push_back(Vector_1D( nCases ));

    for ( i = 0; i < nCases; i++ ) {
        for ( j = 0; j < counter - 1; j++ ) {
            y[j][i] = u[j][i];
        }
        y[counter-1][i] = v[0][i];
    }
    
    /* z = backgHNO3 */
    z.push_back(Vector_1D(backgHNO3.size()));
    for ( i = 0; i < backgHNO3.size(); i++ )
        z[0][i] = backgHNO3[i];

    u = Copy_blocked(y,z[0].size());
    v = Copy_interleaved(z,y[0].size());

    for ( i = 0; i < counter; i++ ) {
        y[i].clear();
    }
    z[0].clear();
    y.clear(); z.clear();

    nCases *= backgHNO3.size();
    counter += 1;
    for ( i = 0; i < counter; i++ )
        y.push_back(Vector_1D( nCases ));

    for ( i = 0; i < nCases; i++ ) {
        for ( j = 0; j < counter - 1; j++ ) {
            y[j][i] = u[j][i];
        }
        y[counter-1][i] = v[0][i];
    }
    
    /* z = backgO3 */
    z.push_back(Vector_1D(backgO3.size()));
    for ( i = 0; i < backgO3.size(); i++ )
        z[0][i] = backgO3[i];

    u = Copy_blocked(y,z[0].size());
    v = Copy_interleaved(z,y[0].size());

    for ( i = 0; i < counter; i++ ) {
        y[i].clear();
    }
    z[0].clear();
    y.clear(); z.clear();

    nCases *= backgO3.size();
    counter += 1;
    for ( i = 0; i < counter; i++ )
        y.push_back(Vector_1D( nCases ));

    for ( i = 0; i < nCases; i++ ) {
        for ( j = 0; j < counter - 1; j++ ) {
            y[j][i] = u[j][i];
        }
        y[counter-1][i] = v[0][i];
    }
    
    /* z = backgCO */
    z.push_back(Vector_1D(backgCO.size()));
    for ( i = 0; i < backgCO.size(); i++ )
        z[0][i] = backgCO[i];

    u = Copy_blocked(y,z[0].size());
    v = Copy_interleaved(z,y[0].size());

    for ( i = 0; i < counter; i++ ) {
        y[i].clear();
    }
    z[0].clear();
    y.clear(); z.clear();

    nCases *= backgCO.size();
    counter += 1;
    for ( i = 0; i < counter; i++ )
        y.push_back(Vector_1D( nCases ));

    for ( i = 0; i < nCases; i++ ) {
        for ( j = 0; j < counter - 1; j++ ) {
            y[j][i] = u[j][i];
        }
        y[counter-1][i] = v[0][i];
    }
    
    /* z = backgCH4 */
    z.push_back(Vector_1D(backgCH4.size()));
    for ( i = 0; i < backgCH4.size(); i++ )
        z[0][i] = backgCH4[i];

    u = Copy_blocked(y,z[0].size());
    v = Copy_interleaved(z,y[0].size());

    for ( i = 0; i < counter; i++ ) {
        y[i].clear();
    }
    z[0].clear();
    y.clear(); z.clear();

    nCases *= backgCH4.size();
    counter += 1;
    for ( i = 0; i < counter; i++ )
        y.push_back(Vector_1D( nCases ));

    for ( i = 0; i < nCases; i++ ) {
        for ( j = 0; j < counter - 1; j++ ) {
            y[j][i] = u[j][i];
        }
        y[counter-1][i] = v[0][i];
    }
    
    /* z = backgSO2 */
    z.push_back(Vector_1D(backgSO2.size()));
    for ( i = 0; i < backgSO2.size(); i++ )
        z[0][i] = backgSO2[i];

    u = Copy_blocked(y,z[0].size());
    v = Copy_interleaved(z,y[0].size());

    for ( i = 0; i < counter; i++ ) {
        y[i].clear();
    }
    z[0].clear();
    y.clear(); z.clear();

    nCases *= backgSO2.size();
    counter += 1;
    for ( i = 0; i < counter; i++ )
        y.push_back(Vector_1D( nCases ));

    for ( i = 0; i < nCases; i++ ) {
        for ( j = 0; j < counter - 1; j++ ) {
            y[j][i] = u[j][i];
        }
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
