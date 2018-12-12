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

const bool VERBOSE = true;
const char* SPACE = " ";
const char* TAB   = "\t";
const char* FILENAME = "input.apcemm";
const char* FILESEP = "/";
std::ifstream inputFile;
const unsigned int FIRSTCOL = 26;

std::string line;


void Read_Input_File( const std::string folderDir )
{

    /* Read\_Input\_File is the driver program for reading 
     * APCEMM input file "input.apcemm" from disk */

    char* TOPTITLE;
    const std::string HEADER(80, '*');

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Read\_Input\_File begins here ! 
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    /* Echo output */
    std::cout << HEADER << std::endl;
    std::cout << " A P C E M M   U S E R   I N P U T" << std::endl;
    std::cout << fileName << std::endl;

    /* Assume success */
    RC = SUCCESS; 

    /* Open file */
    std::string fullPath = folderDir + FILESEP + FILENAME;
    inputFile.open( fullPath.c_str() );
    if ( !inputFile ) {
        /* Call error */
        std::cout << " APCEMM I/O ERROR: No such file in directory" << std::endl;
        exit(1);
    }

    /* Read TOPTITLE */
    inputFile >> TOPTITLE;

    /* Loop until EOF */
    while ( inputFile.good() ) {

        /* Read a line from the file, print it to the console, exit if EOF */
        inputFile >> line;
        line.replace(line.begin(), line.end(), TAB, SPACE);
        if ( VERBOSE )
            std::cout << line << std::endl;


        if ( strstr( line, "SIMULATION MENU" ) != NULL ) {
            Read_Simulation_Menu( Input_Opt, RC ); 

        } else if ( strstr( line, "PARAMETER SWEEP" ) != NULL ) {
//            Read_Parameters( Input_Opt );

        } else if ( strstr( line, "END OF FILE" ) != NULL )
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
        substring.push_back(token);
        line2split.erase(0, pos + delimiter.length());
    }
    
    return substring;

} /* End of Split_Line */

void Read_Simulation_Menu( Option &Input_Opt, bool RC )
{
    
    /* DESCRIPTION: Function Read\_Simulation\_Menu reads the SIMULATION MENU
     * section of the APCEMM input file. */

    /* INPUT/OUTPUT PARAMETERS:
     * - Input_Opt: Input options
     *
     * OUTPUT PARAMETERS:
     * - RC: Success or failure
     */

    /* Read until all lines from the menu are read */

    std::vector<std::string> tokens;

    /* ==================================================== */
    /* Parameter sweep?                                     */
    /* ==================================================== */

    inputFile >> line; line.replace(line.begin(), line.end(), TAB, SPACE);
    if ( VERBOSE )
        std::cout << line << std::endl;
    
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );
    Input_Opt.SIMULATION_PARAMETER_SWEEP = tokens[0].c_str();

    /* ==================================================== */
    /* Output folder                                        */
    /* ==================================================== */

    inputFile >> line; line.replace(line.begin(), line.end(), TAB, SPACE);
    if ( VERBOSE )
        std::cout << line << std::endl;
    
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );
    Input_Opt.SIMULATION_OUTPUT_FOLDER = tokens[0].c_str();

    /* ==================================================== */
    /* Run directory                                        */
    /* ==================================================== */

    inputFile >> line; line.replace(line.begin(), line.end(), TAB, SPACE);
    if ( VERBOSE )
        std::cout << line << std::endl;
    
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );
    Input_Opt.SIMULATION_RUN_DIRECTORY = tokens[0].c_str();
    
    /* ==================================================== */
    /* Input background condition                           */
    /* ==================================================== */

    inputFile >> line; line.replace(line.begin(), line.end(), TAB, SPACE);
    if ( VERBOSE )
        std::cout << line << std::endl;
    
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );
    Input_Opt.SIMULATION_INPUT_BACKG_COND = tokens[0].c_str();
    
    /* ==================================================== */
    /* Save forward results                                 */
    /* ==================================================== */

    inputFile >> line; line.replace(line.begin(), line.end(), TAB, SPACE);
    if ( VERBOSE )
        std::cout << line << std::endl;
    
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );
    Input_Opt.SIMULATION_SAVE_FORWARD = tokens[0].c_str();
    
    /* ==================================================== */
    /* netCDF file name                                     */
    /* ==================================================== */

    inputFile >> line; line.replace(line.begin(), line.end(), TAB, SPACE);
    if ( VERBOSE )
        std::cout << line << std::endl;
    
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );
    Input_Opt.SIMULATION_FORWARD_FILENAME = \
            tokens[0].replace(tokens[0].begin(), tokens[0].end(), "*", "");
    
    /* ==================================================== */
    /* Adjoint Optimization                                 */
    /* ==================================================== */

    inputFile >> line; line.replace(line.begin(), line.end(), TAB, SPACE);
    if ( VERBOSE )
        std::cout << line << std::endl;
    
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );
    Input_Opt.SIMULATION_ADJOINT = tokens[0].c_str();
    
    /* ==================================================== */
    /* netCDF file name                                     */
    /* ==================================================== */

    inputFile >> line; line.replace(line.begin(), line.end(), TAB, SPACE);
    if ( VERBOSE )
        std::cout << line << std::endl;
    
    tokens = Split_Line( line.substr(FIRSTCOL), SPACE );
    Input_Opt.SIMULATION_ADJOINT_FILENAME = tokens[0].c_str();

    /* Return success */
    RC = SUCCESS; 

    /* ==================================================== */
    /* Print to screen                                      */
    /* ==================================================== */

    std::cout << " %%% SIMULATION MENU %%% :" << std::endl;
    std::cout << " ------------------------+---------------------------------------- " << std::endl; 
    std::cout << " Parameter sweep?        : " << Input_Opt.SIMULATION_PARAMETER_SWEEP << std::endl;
    std::cout << " Output folder           : " << Input_Opt.SIMULATION_OUTPUT_FOLDER << std::endl; 
    std::cout << " Run directory           : " << Input_Opt.SIMULATION_RUN_DIRECTORY << std::endl; 
    std::cout << " Input backgrd condition : " << Input_Opt.SIMULATION_INPUT_BACKG_COND << std::endl; 
    std::cout << " Save Forward results    : " << Input_Opt.SIMULATION_SAVE_FORWARD << std::endl; 
    std::cout << "  => netCDF file name    : " << Input_Opt.SIMULATION_FORWARD_FILENAME << std::endl; 
    std::cout << " Turn on adjoint optim.  : " << Input_Opt.SIMULATION_ADJOINT << std::endl; 
    std::cout << "  => netCDF file name    : " << Input_OPT.SIMULATION_ADJOINT_FILENAME << std::endl; 
    
} /* End of Read_Simulation_Menu */


/* End of Input_Reader.cpp */
