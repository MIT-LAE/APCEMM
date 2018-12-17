/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Diag_Mod Program File                                            */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 12/17/2018                                */
/* File                 : Diag_Mod.cpp                              */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "Core/Diag_Mod.hpp"

bool Diag_TS( const char* rootName, const int hh, const int mm, \
              const std::vector<int> speciesIndices, \
              const Solution& Data, const Mesh& m )
{

    const bool doWrite   = 1;
    const bool doRead    = 1;
    const bool overWrite = 1;

    std::string fileName( rootName );
    size_t start_pos;
    bool found;

    char hh_string[10];
    sprintf(hh_string, "%02d", hh );
    char mm_string[10];
    sprintf(mm_string, "%02d", mm );

    /* Replacing "hh" with hour number since start */
    start_pos = 0; 
    found = 0;
    while ( (start_pos = fileName.find("hh", start_pos)) != std::string::npos ) {
        fileName.replace(start_pos, 2, hh_string);
        start_pos += 2;
        found = 1;
    }

    if ( found = 0 ) {
        std::cout << " Diagnostic filename must be of the form *hhmm.nc. Aborting!" << std::endl;
        std::cout << " filename: " << rootName << std::endl;
        exit(1);
    }

    /* Replacing "mm" with minute number since start */
    start_pos = 0; 
    found = 0;
    while ( (start_pos = fileName.find("mm", start_pos)) != std::string::npos ) {
        fileName.replace(start_pos, 2, mm_string);
        start_pos += 2;
    }
    
    if ( found = 0 ) {
        std::cout << " Diagnostic filename must be of the form *hhmm.nc. Aborting!" << std::endl;
        std::cout << " filename: " << rootName << std::endl;
        exit(1);
    }

    const char* outFile = fileName.c_str();

    FileHandler fileHandler( outFile, doWrite, doRead, overWrite );
    NcFile currFile = fileHandler.openFile();
    if ( !fileHandler.isFileOpen() ) {
        std::cout << " File " << outFile << " didn't open!" << "\n";
        std::cout << " Do you have write permission?" << std::endl;                                 
        return 0;
    } else {
        
        int didSaveSucceed = 1;
        time_t rawtime;
        char buffer[80];
        time( &rawtime );
        strftime(buffer, sizeof(buffer),"%d-%m-%Y %H:%M:%S", localtime(&rawtime));

        const NcDim *xDim = fileHandler.addDim( currFile, "X Center", long(m.Nx()) );
        didSaveSucceed *= fileHandler.addVar( currFile, &(m.x())[0], "Grid cell horizontal center values", xDim, "float", "m", "X Center values");

        const NcDim *yDim = fileHandler.addDim( currFile, "X Center", long(m.Nx()) );
        didSaveSucceed *= fileHandler.addVar( currFile, &(m.y())[0], "Grid cell vertical center values", yDim, "float", "m", "Y Center values");
            
        didSaveSucceed *= fileHandler.addAtt( currFile, "FileName", fileHandler.getFileName() );
        didSaveSucceed *= fileHandler.addAtt( currFile, "Author", "Thibaud M. Fritz (fritzt@mit.edu)" );
        didSaveSucceed *= fileHandler.addAtt( currFile, "Contact", "Thibaud M. Fritz (fritzt@mit.edu)" );
        didSaveSucceed *= fileHandler.addAtt( currFile, "Generation Date", buffer );
        didSaveSucceed *= fileHandler.addAtt( currFile, "Format", "NetCDF-4" );

#if ( SAVE_TO_DOUBLE )
        double* array;
        const char* outputType = "double";
#else
        float* array;
        const char* outputType = "float";
#endif /* SAVE_TO_DOUBLE */

        /* Start saving species ... */

        /* Define conversion factor */
        double scalingFactor;
        char charSpc[30];
        char charName[50];
        char charUnit[20];

        strncpy( charSpc, "O3", sizeof(charSpc) );
        strncpy( charName, "O3 molecular concentration", sizeof(charName) );
        scalingFactor = 1.0E+00;
        strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );
        
#if ( SAVE_TO_DOUBLE )
        array = util::vect2double( Data.O3, m.Nx(), m.Ny(), scalingFactor );
#else
        array = util::vect2float ( Data.O3, m.Nx(), m.Ny(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

        didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0], (const char*)charSpc, xDim, yDim, outputType, (const char*)charUnit, (const char*)charName );

    }

    return 1;

} /* End of Diag_TS */

/* End of Diag_Mod.cpp */
