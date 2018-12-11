/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Input_Reader Program File                                        */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 12/10/2018                                */
/* File                 : Input_Reader.cpp                          */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "Core/Input_Reader.hpp"

void Read_Input_File( const char* fileName )
{

    /* Read\_Input\_File is the driver program for reading 
     * APCEMM input file "input.apcemm" from disk */

    const char* SPACE = " ";
    char* TOPTITLE;
    char* line;
    char*

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Read\_Input\_File begins here ! 
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    /* Echo output */
    std::cout << std::string(100, '*') << std::endl;
    std::cout << " A P C E M M   U S E R   I N P U T" << std::endl;
    std::cout << fileName << std::endl;

    ifstream inputFile;

    /* Open file */
    inputFile.open( fileName );
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
        std::cout << line << std::endl;

        if ( strstr( line, "SIMULATION MENU" ) != NULL ) {


        } else if ( strstr( line, "PARAMETER SWEEP" ) != NULL ) {


        } else if ( strstr( line, "END OF FILE" ) != NULL )
            break;
                

    }

    inputFile.close();

    return;

} /* End of Read_Input_File */

/* End of Input_Reader.cpp */
