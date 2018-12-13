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

    if ( const char* simDir = std::getenv("currDir") )
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



/* End of Input_Reader.cpp */
