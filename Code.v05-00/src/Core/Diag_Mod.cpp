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

bool Diag_TS_Chem( const char* rootName,                     \
                   const std::vector<int> speciesIndices,    \
                   const int hh, const int mm, const int ss, \
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
    char ss_string[10];
    sprintf(ss_string, "%02d", ss );

    /* Replacing "hh" with hour number since start */
    start_pos = 0; 
    found = 0;
    while ( (start_pos = fileName.find("hh", start_pos)) != std::string::npos ) {
        fileName.replace(start_pos, 2, hh_string);
        start_pos += 2;
        found = 1;
    }

    if ( found = 0 ) {
        std::cout << " Diagnostic filename must be of the form *hhmmss.nc or *hhmm.nc. Aborting!" << std::endl;
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
        std::cout << " Diagnostic filename must be of the form *hhmmss.nc or *hhmm.nc. Aborting!" << std::endl;
        std::cout << " filename: " << rootName << std::endl;
        exit(1);
    }

    /* Replacing "ss" with minute number since start if possible */
    start_pos = 0;
    found = 0;
    while ( (start_pos = fileName.find("ss", start_pos)) != std::string::npos ) {
        fileName.replace(start_pos, 2, ss_string);
        start_pos += 2;
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
        didSaveSucceed *= fileHandler.addVar( currFile, &(m.x())[0], "X Centers", xDim, "float", "m", "Grid cell horizontal centers");

        const NcDim *yDim = fileHandler.addDim( currFile, "Y Center", long(m.Ny()) );
        didSaveSucceed *= fileHandler.addVar( currFile, &(m.y())[0], "Y Centers", yDim, "float", "m", "Grid cell vertical centers");
            
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

        for ( unsigned int i = 0; i < speciesIndices.size(); i++ ) {
            
            if ( speciesIndices[i] - 1 == ind_CO2 ) {

                strncpy( charSpc, "CO2", sizeof(charSpc) );
                strncpy( charName, "CO2 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );
            
#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.CO2, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.CO2, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_PPN ) {

                strncpy( charSpc, "PPN", sizeof(charSpc) );
                strncpy( charName, "PPN molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.PPN, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.PPN, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_BrNO2 ) {

                strncpy( charSpc, "BrNO2", sizeof(charSpc) );
                strncpy( charName, "BrNO2 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.BrNO2, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.BrNO2, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_IEPOX ) {

                strncpy( charSpc, "IEPOX", sizeof(charSpc) );
                strncpy( charName, "IEPOX molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.IEPOX, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.IEPOX, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_PMNN ) {

                strncpy( charSpc, "PMNN", sizeof(charSpc) );
                strncpy( charName, "PMNN molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.PMNN, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.PMNN, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_N2O ) {

                strncpy( charSpc, "N2O", sizeof(charSpc) );
                strncpy( charName, "N2O molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.N2O, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.N2O, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_N ) {

                strncpy( charSpc, "N", sizeof(charSpc) );
                strncpy( charName, "N molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.N, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.N, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_PAN ) {

                strncpy( charSpc, "PAN", sizeof(charSpc) );
                strncpy( charName, "PAN molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.PAN, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.PAN, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_ALK4 ) {

                strncpy( charSpc, "ALK4", sizeof(charSpc) );
                strncpy( charName, "ALK4 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.ALK4, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.ALK4, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_MAP ) {

                strncpy( charSpc, "MAP", sizeof(charSpc) );
                strncpy( charName, "MAP molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.MAP, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.MAP, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_MPN ) {

                strncpy( charSpc, "MPN", sizeof(charSpc) );
                strncpy( charName, "MPN molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.MPN, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.MPN, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_Cl2O2 ) {

                strncpy( charSpc, "Cl2O2", sizeof(charSpc) );
                strncpy( charName, "Cl2O2 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.Cl2O2, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.Cl2O2, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_ETP ) {

                strncpy( charSpc, "ETP", sizeof(charSpc) );
                strncpy( charName, "ETP molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.ETP, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.ETP, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_HNO2 ) {

                strncpy( charSpc, "HNO2", sizeof(charSpc) );
                strncpy( charName, "HNO2 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.HNO2, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.HNO2, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_C3H8 ) {

                strncpy( charSpc, "C3H8", sizeof(charSpc) );
                strncpy( charName, "C3H8 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.C3H8, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.C3H8, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_RA3P ) {

                strncpy( charSpc, "RA3P", sizeof(charSpc) );
                strncpy( charName, "RA3P molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.RA3P, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.RA3P, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_RB3P ) {

                strncpy( charSpc, "RB3P", sizeof(charSpc) );
                strncpy( charName, "RB3P molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.RB3P, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.RB3P, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_OClO ) {

                strncpy( charSpc, "OClO", sizeof(charSpc) );
                strncpy( charName, "OClO molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.OClO, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.OClO, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_ClNO2 ) {

                strncpy( charSpc, "ClNO2", sizeof(charSpc) );
                strncpy( charName, "ClNO2 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.ClNO2, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.ClNO2, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_ISOP ) {

                strncpy( charSpc, "ISOP", sizeof(charSpc) );
                strncpy( charName, "ISOP molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.ISOP, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.ISOP, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_HNO4 ) {

                strncpy( charSpc, "HNO4", sizeof(charSpc) );
                strncpy( charName, "HNO4 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.HNO4, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.HNO4, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_MAOP ) {

                strncpy( charSpc, "MAOP", sizeof(charSpc) );
                strncpy( charName, "MAOP molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.MAOP, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.MAOP, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_MP ) {

                strncpy( charSpc, "MP", sizeof(charSpc) );
                strncpy( charName, "MP molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.MP, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.MP, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_ClOO ) {

                strncpy( charSpc, "ClOO", sizeof(charSpc) );
                strncpy( charName, "ClOO molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.ClOO, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.ClOO, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_RP ) {

                strncpy( charSpc, "RP", sizeof(charSpc) );
                strncpy( charName, "RP molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.RP, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.RP, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_BrCl ) {

                strncpy( charSpc, "BrCl", sizeof(charSpc) );
                strncpy( charName, "BrCl molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.BrCl, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.BrCl, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_PP ) {

                strncpy( charSpc, "PP", sizeof(charSpc) );
                strncpy( charName, "PP molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.PP, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.PP, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_PRPN ) {

                strncpy( charSpc, "PRPN", sizeof(charSpc) );
                strncpy( charName, "PRPN molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.PRPN, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.PRPN, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_SO4 ) {

                strncpy( charSpc, "SO4", sizeof(charSpc) );
                strncpy( charName, "SO4 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.SO4, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.SO4, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_Br2 ) {

                strncpy( charSpc, "Br2", sizeof(charSpc) );
                strncpy( charName, "Br2 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.Br2, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.Br2, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_ETHLN ) {

                strncpy( charSpc, "ETHLN", sizeof(charSpc) );
                strncpy( charName, "ETHLN molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.ETHLN, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.ETHLN, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_MVKN ) {

                strncpy( charSpc, "MVKN", sizeof(charSpc) );
                strncpy( charName, "MVKN molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.MVKN, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.MVKN, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_R4P ) {

                strncpy( charSpc, "R4P", sizeof(charSpc) );
                strncpy( charName, "R4P molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.R4P, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.R4P, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_C2H6 ) {

                strncpy( charSpc, "C2H6", sizeof(charSpc) );
                strncpy( charName, "C2H6 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.C2H6, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.C2H6, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_RIP ) {

                strncpy( charSpc, "RIP", sizeof(charSpc) );
                strncpy( charName, "RIP molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.RIP, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.RIP, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_VRP ) {

                strncpy( charSpc, "VRP", sizeof(charSpc) );
                strncpy( charName, "VRP molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.VRP, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.VRP, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_ATOOH ) {

                strncpy( charSpc, "ATOOH", sizeof(charSpc) );
                strncpy( charName, "ATOOH molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.ATOOH, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.ATOOH, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_IAP ) {

                strncpy( charSpc, "IAP", sizeof(charSpc) );
                strncpy( charName, "IAP molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.IAP, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.IAP, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_DHMOB ) {

                strncpy( charSpc, "DHMOB", sizeof(charSpc) );
                strncpy( charName, "DHMOB molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.DHMOB, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.DHMOB, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_MOBA ) {

                strncpy( charSpc, "MOBA", sizeof(charSpc) );
                strncpy( charName, "MOBA molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.MOBA, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.MOBA, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_MRP ) {

                strncpy( charSpc, "MRP", sizeof(charSpc) );
                strncpy( charName, "MRP molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.MRP, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.MRP, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_N2O5 ) {

                strncpy( charSpc, "N2O5", sizeof(charSpc) );
                strncpy( charName, "N2O5 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.N2O5, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.N2O5, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_ISNOHOO ) {

                strncpy( charSpc, "ISNOHOO", sizeof(charSpc) );
                strncpy( charName, "ISNOHOO molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.ISNOHOO, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.ISNOHOO, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_ISNP ) {

                strncpy( charSpc, "ISNP", sizeof(charSpc) );
                strncpy( charName, "ISNP molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.ISNP, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.ISNP, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_ISOPNB ) {

                strncpy( charSpc, "ISOPNB", sizeof(charSpc) );
                strncpy( charName, "ISOPNB molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.ISOPNB, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.ISOPNB, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_IEPOXOO ) {

                strncpy( charSpc, "IEPOXOO", sizeof(charSpc) );
                strncpy( charName, "IEPOXOO molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.IEPOXOO, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.IEPOXOO, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_MACRNO2 ) {

                strncpy( charSpc, "MACRNO2", sizeof(charSpc) );
                strncpy( charName, "MACRNO2 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.MACRNO2, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.MACRNO2, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_ROH ) {

                strncpy( charSpc, "ROH", sizeof(charSpc) );
                strncpy( charName, "ROH molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.ROH, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.ROH, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_MOBAOO ) {

                strncpy( charSpc, "MOBAOO", sizeof(charSpc) );
                strncpy( charName, "MOBAOO molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.MOBAOO, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.MOBAOO, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_DIBOO ) {

                strncpy( charSpc, "DIBOO", sizeof(charSpc) );
                strncpy( charName, "DIBOO molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.DIBOO, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.DIBOO, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_PMN ) {

                strncpy( charSpc, "PMN", sizeof(charSpc) );
                strncpy( charName, "PMN molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.PMN, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.PMN, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_ISNOOB ) {

                strncpy( charSpc, "ISNOOB", sizeof(charSpc) );
                strncpy( charName, "ISNOOB molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.ISNOOB, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.ISNOOB, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_INPN ) {

                strncpy( charSpc, "INPN", sizeof(charSpc) );
                strncpy( charName, "INPN molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.INPN, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.INPN, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_H ) {

                strncpy( charSpc, "H", sizeof(charSpc) );
                strncpy( charName, "H molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.H, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.H, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_BrNO3 ) {

                strncpy( charSpc, "BrNO3", sizeof(charSpc) );
                strncpy( charName, "BrNO3 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.BrNO3, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.BrNO3, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_PRPE ) {

                strncpy( charSpc, "PRPE", sizeof(charSpc) );
                strncpy( charName, "PRPE molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.PRPE, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.PRPE, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_MVKOO ) {

                strncpy( charSpc, "MVKOO", sizeof(charSpc) );
                strncpy( charName, "MVKOO molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.MVKOO, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.MVKOO, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_Cl2 ) {

                strncpy( charSpc, "Cl2", sizeof(charSpc) );
                strncpy( charName, "Cl2 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.Cl2, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.Cl2, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_ISOPND ) {

                strncpy( charSpc, "ISOPND", sizeof(charSpc) );
                strncpy( charName, "ISOPND molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.ISOPND, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.ISOPND, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_HOBr ) {

                strncpy( charSpc, "HOBr", sizeof(charSpc) );
                strncpy( charName, "HOBr molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.HOBr, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.HOBr, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_A3O2 ) {

                strncpy( charSpc, "A3O2", sizeof(charSpc) );
                strncpy( charName, "A3O2 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.A3O2, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.A3O2, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_PROPNN ) {

                strncpy( charSpc, "PROPNN", sizeof(charSpc) );
                strncpy( charName, "PROPNN molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.PROPNN, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.PROPNN, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_GLYX ) {

                strncpy( charSpc, "GLYX", sizeof(charSpc) );
                strncpy( charName, "GLYX molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.GLYX, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.GLYX, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_MAOPO2 ) {

                strncpy( charSpc, "MAOPO2", sizeof(charSpc) );
                strncpy( charName, "MAOPO2 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.MAOPO2, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.MAOPO2, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_CH4 ) {

                strncpy( charSpc, "CH4", sizeof(charSpc) );
                strncpy( charName, "CH4 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.CH4, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.CH4, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_GAOO ) {

                strncpy( charSpc, "GAOO", sizeof(charSpc) );
                strncpy( charName, "GAOO molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.GAOO, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.GAOO, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_B3O2 ) {

                strncpy( charSpc, "B3O2", sizeof(charSpc) );
                strncpy( charName, "B3O2 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.B3O2, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.B3O2, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_ACET ) {

                strncpy( charSpc, "ACET", sizeof(charSpc) );
                strncpy( charName, "ACET molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.ACET, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.ACET, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_MACRN ) {

                strncpy( charSpc, "MACRN", sizeof(charSpc) );
                strncpy( charName, "MACRN molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.MACRN, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.MACRN, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_CH2OO ) {

                strncpy( charSpc, "CH2OO", sizeof(charSpc) );
                strncpy( charName, "CH2OO molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.CH2OO, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.CH2OO, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_MGLYOO ) {

                strncpy( charSpc, "MGLYOO", sizeof(charSpc) );
                strncpy( charName, "MGLYOO molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.MGLYOO, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.MGLYOO, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_VRO2 ) {

                strncpy( charSpc, "VRO2", sizeof(charSpc) );
                strncpy( charName, "VRO2 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.VRO2, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.VRO2, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_MGLOO ) {

                strncpy( charSpc, "MGLOO", sizeof(charSpc) );
                strncpy( charName, "MGLOO molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.MGLOO, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.MGLOO, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_MACROO ) {

                strncpy( charSpc, "MACROO", sizeof(charSpc) );
                strncpy( charName, "MACROO molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.MACROO, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.MACROO, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_PO2 ) {

                strncpy( charSpc, "PO2", sizeof(charSpc) );
                strncpy( charName, "PO2 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.PO2, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.PO2, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_CH3CHOO ) {

                strncpy( charSpc, "CH3CHOO", sizeof(charSpc) );
                strncpy( charName, "CH3CHOO molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.CH3CHOO, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.CH3CHOO, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_MAN2 ) {

                strncpy( charSpc, "MAN2", sizeof(charSpc) );
                strncpy( charName, "MAN2 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.MAN2, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.MAN2, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_ISNOOA ) {

                strncpy( charSpc, "ISNOOA", sizeof(charSpc) );
                strncpy( charName, "ISNOOA molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.ISNOOA, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.ISNOOA, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_H2O2 ) {

                strncpy( charSpc, "H2O2", sizeof(charSpc) );
                strncpy( charName, "H2O2 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.H2O2, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.H2O2, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_PRN1 ) {

                strncpy( charSpc, "PRN1", sizeof(charSpc) );
                strncpy( charName, "PRN1 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.PRN1, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.PRN1, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_ETO2 ) {

                strncpy( charSpc, "ETO2", sizeof(charSpc) );
                strncpy( charName, "ETO2 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.ETO2, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.ETO2, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_KO2 ) {

                strncpy( charSpc, "KO2", sizeof(charSpc) );
                strncpy( charName, "KO2 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.KO2, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.KO2, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_RCO3 ) {

                strncpy( charSpc, "RCO3", sizeof(charSpc) );
                strncpy( charName, "RCO3 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.RCO3, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.RCO3, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_HC5OO ) {

                strncpy( charSpc, "HC5OO", sizeof(charSpc) );
                strncpy( charName, "HC5OO molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.HC5OO, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.HC5OO, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_GLYC ) {

                strncpy( charSpc, "GLYC", sizeof(charSpc) );
                strncpy( charName, "GLYC molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.GLYC, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.GLYC, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_ClNO3 ) {

                strncpy( charSpc, "ClNO3", sizeof(charSpc) );
                strncpy( charName, "ClNO3 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.ClNO3, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.ClNO3, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_RIO2 ) {

                strncpy( charSpc, "RIO2", sizeof(charSpc) );
                strncpy( charName, "RIO2 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.RIO2, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.RIO2, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_R4N1 ) {

                strncpy( charSpc, "R4N1", sizeof(charSpc) );
                strncpy( charName, "R4N1 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.R4N1, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.R4N1, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_HOCl ) {

                strncpy( charSpc, "HOCl", sizeof(charSpc) );
                strncpy( charName, "HOCl molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.HOCl, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.HOCl, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_ATO2 ) {

                strncpy( charSpc, "ATO2", sizeof(charSpc) );
                strncpy( charName, "ATO2 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.ATO2, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.ATO2, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_HNO3 ) {

                strncpy( charSpc, "HNO3", sizeof(charSpc) );
                strncpy( charName, "HNO3 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.HNO3, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.HNO3, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_ISN1 ) {

                strncpy( charSpc, "ISN1", sizeof(charSpc) );
                strncpy( charName, "ISN1 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.ISN1, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.ISN1, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_MAO3 ) {

                strncpy( charSpc, "MAO3", sizeof(charSpc) );
                strncpy( charName, "MAO3 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.MAO3, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.MAO3, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_MRO2 ) {

                strncpy( charSpc, "MRO2", sizeof(charSpc) );
                strncpy( charName, "MRO2 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.MRO2, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.MRO2, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_INO2 ) {

                strncpy( charSpc, "INO2", sizeof(charSpc) );
                strncpy( charName, "INO2 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.INO2, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.INO2, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_HAC ) {

                strncpy( charSpc, "HAC", sizeof(charSpc) );
                strncpy( charName, "HAC molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.HAC, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.HAC, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_HC5 ) {

                strncpy( charSpc, "HC5", sizeof(charSpc) );
                strncpy( charName, "HC5 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.HC5, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.HC5, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_MGLY ) {

                strncpy( charSpc, "MGLY", sizeof(charSpc) );
                strncpy( charName, "MGLY molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.MGLY, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.MGLY, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_ISOPNBO2 ) {

                strncpy( charSpc, "ISOPNBO2", sizeof(charSpc) );
                strncpy( charName, "ISOPNBO2 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.ISOPNBO2, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.ISOPNBO2, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_ISOPNDO2 ) {

                strncpy( charSpc, "ISOPNDO2", sizeof(charSpc) );
                strncpy( charName, "ISOPNDO2 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.ISOPNDO2, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.ISOPNDO2, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_R4O2 ) {

                strncpy( charSpc, "R4O2", sizeof(charSpc) );
                strncpy( charName, "R4O2 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.R4O2, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.R4O2, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_R4N2 ) {

                strncpy( charSpc, "R4N2", sizeof(charSpc) );
                strncpy( charName, "R4N2 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.R4N2, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.R4N2, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_BrO ) {

                strncpy( charSpc, "BrO", sizeof(charSpc) );
                strncpy( charName, "BrO molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.BrO, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.BrO, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_RCHO ) {

                strncpy( charSpc, "RCHO", sizeof(charSpc) );
                strncpy( charName, "RCHO molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.RCHO, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.RCHO, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_MEK ) {

                strncpy( charSpc, "MEK", sizeof(charSpc) );
                strncpy( charName, "MEK molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.MEK, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.MEK, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_ClO ) {

                strncpy( charSpc, "ClO", sizeof(charSpc) );
                strncpy( charName, "ClO molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.ClO, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.ClO, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_MACR ) {

                strncpy( charSpc, "MACR", sizeof(charSpc) );
                strncpy( charName, "MACR molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.MACR, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.MACR, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_SO2 ) {

                strncpy( charSpc, "SO2", sizeof(charSpc) );
                strncpy( charName, "SO2 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.SO2, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.SO2, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_MVK ) {

                strncpy( charSpc, "MVK", sizeof(charSpc) );
                strncpy( charName, "MVK molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.MVK, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.MVK, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_ALD2 ) {

                strncpy( charSpc, "ALD2", sizeof(charSpc) );
                strncpy( charName, "ALD2 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.ALD2, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.ALD2, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_MCO3 ) {

                strncpy( charSpc, "MCO3", sizeof(charSpc) );
                strncpy( charName, "MCO3 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.MCO3, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.MCO3, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_CH2O ) {

                strncpy( charSpc, "CH2O", sizeof(charSpc) );
                strncpy( charName, "CH2O molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.CH2O, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.CH2O, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_H2O ) {

                strncpy( charSpc, "H2O", sizeof(charSpc) );
                strncpy( charName, "H2O molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.H2O, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.H2O, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_Br ) {

                strncpy( charSpc, "Br", sizeof(charSpc) );
                strncpy( charName, "Br molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.Br, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.Br, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_NO ) {

                strncpy( charSpc, "NO", sizeof(charSpc) );
                strncpy( charName, "NO molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.NO, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.NO, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_NO3 ) {

                strncpy( charSpc, "NO3", sizeof(charSpc) );
                strncpy( charName, "NO3 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.NO3, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.NO3, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_Cl ) {

                strncpy( charSpc, "Cl", sizeof(charSpc) );
                strncpy( charName, "Cl molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.Cl, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.Cl, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_O ) {

                strncpy( charSpc, "O", sizeof(charSpc) );
                strncpy( charName, "O molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.O, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.O, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_O1D ) {

                strncpy( charSpc, "O1D", sizeof(charSpc) );
                strncpy( charName, "O1D molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.O1D, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.O1D, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_O3 ) {

                strncpy( charSpc, "O3", sizeof(charSpc) );
                strncpy( charName, "O3 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.O3, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.O3, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_HO2 ) {

                strncpy( charSpc, "HO2", sizeof(charSpc) );
                strncpy( charName, "HO2 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.HO2, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.HO2, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_NO2 ) {

                strncpy( charSpc, "NO2", sizeof(charSpc) );
                strncpy( charName, "NO2 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.NO2, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.NO2, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_OH ) {

                strncpy( charSpc, "OH", sizeof(charSpc) );
                strncpy( charName, "OH molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.OH, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.OH, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_HBr ) {

                strncpy( charSpc, "HBr", sizeof(charSpc) );
                strncpy( charName, "HBr molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.HBr, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.HBr, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_HCl ) {

                strncpy( charSpc, "HCl", sizeof(charSpc) );
                strncpy( charName, "HCl molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.HCl, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.HCl, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_CO ) {

                strncpy( charSpc, "CO", sizeof(charSpc) );
                strncpy( charName, "CO molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.CO, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.CO, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else if ( speciesIndices[i] - 1 == ind_MO2 ) {

                strncpy( charSpc, "MO2", sizeof(charSpc) );
                strncpy( charName, "MO2 molecular concentration", sizeof(charName) );
                scalingFactor = 1.0E+00;
                strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                array = util::vect2double( Data.MO2, m.Ny(), m.Nx(), scalingFactor );
#else
                array = util::vect2float ( Data.MO2, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                             (const char*)charSpc, yDim, xDim,  \
                                             outputType, (const char*)charUnit, \
                                             (const char*)charName );

            } else {
                std::cout << " In Diag_Mod for timeseries: Unexpected index: " << speciesIndices[i] << std::endl;
                std::cout << " Ignoring that index..." << std::endl;
            }

        }

        delete[] array; array = NULL;

        if ( didSaveSucceed == NC_SUCCESS ) {
//            std::cout << " Done saving to netCDF!" << "\n";
        } else {
            std::cout << " Error occured in save data: didSaveSucceed: " << didSaveSucceed << std::endl;
            return SAVE_FAILURE;
        }
            
        fileHandler.closeFile( currFile );
        if ( fileHandler.isFileOpen() ) {
            std::cout << "File " << outFile << " did not close properly!" << "\n";
            return SAVE_FAILURE;
        }

    }

    return SAVE_SUCCESS;

} /* End of Diag_TS_Chem */

bool Diag_TS_Phys( const char* rootName,                     \
                   const std::vector<int> aerosolIndices,    \
                   const int hh, const int mm, const int ss, \
                   const Solution& Data, const Mesh& m,      \
                   const Meteorology &met,                   \
                   const int outputPDF )
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
    char ss_string[10];
    sprintf(ss_string, "%02d", ss );

    /* Replacing "hh" with hour number since start */
    start_pos = 0; 
    found = 0;
    while ( (start_pos = fileName.find("hh", start_pos)) != std::string::npos ) {
        fileName.replace(start_pos, 2, hh_string);
        start_pos += 2;
        found = 1;
    }

    if ( found = 0 ) {
        std::cout << " Diagnostic filename must be of the form *hhmmss.nc or *hhmm.nc. Aborting!" << std::endl;
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
        std::cout << " Diagnostic filename must be of the form *hhmmss.nc or *hhmm.nc. Aborting!" << std::endl;
        std::cout << " filename: " << rootName << std::endl;
        exit(1);
    }

    /* Replacing "ss" with minute number since start if possible */
    start_pos = 0;
    found = 0;
    while ( (start_pos = fileName.find("ss", start_pos)) != std::string::npos ) {
        fileName.replace(start_pos, 2, ss_string);
        start_pos += 2;
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
        didSaveSucceed *= fileHandler.addVar( currFile, &(m.x())[0], "X Centers", xDim, "float", "m", "Grid cell horizontal centers");

        const NcDim *yDim = fileHandler.addDim( currFile, "Y Center", long(m.Ny()) );
        didSaveSucceed *= fileHandler.addVar( currFile, &(m.y())[0], "Y Centers", yDim, "float", "m", "Grid cell vertical centers");

        const NcDim *binRad = fileHandler.addDim( currFile, "Bin mid-radius", Data.nBin_PA );
        didSaveSucceed *= fileHandler.addVar( currFile, &(Data.solidAerosol.binCenters())[0], "Bin mid-radius", binRad, "float", "m", "Ice bin center radius" );
        const NcDim *binEdge = fileHandler.addDim( currFile, "Bin edge radius", Data.nBin_PA+1 );
        didSaveSucceed *= fileHandler.addVar( currFile, &(Data.solidAerosol.binEdges())[0], "Bin edge radius", binEdge, "float", "m", "Ice bin edge radius" );

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


        /* Define conversion factor */
        double scalingFactor = 1.0;
        char charSpc[30];
        char charName[50];
        char charUnit[20];

        /* Output met */

        /* Check if evolving met. If not only save met in one file */

        /* Saving meteorological pressure */
        /* TODO: Implement 2D met pressure or at least check if press is a 2D
         * vector... */

        strncpy( charSpc, "Pressure", sizeof(charSpc) );
        strncpy( charName, "Pressure", sizeof(charName) );
        strncpy( charUnit, "Pa", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
        array = util::vect2double( met.Press(), m.Ny(), scalingFactor );
#else
        array = util::vect2float ( met.Press(), m.Ny(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

        didSaveSucceed *= fileHandler.addVar( currFile, &(array)[0],    \
                                     (const char*)charSpc,              \
                                     yDim, outputType,                  \
                                     (const char*)charUnit,             \
                                     (const char*)charName );

        /* Saving meteorological temperature */

        strncpy( charSpc, "Temperature", sizeof(charSpc) );
        strncpy( charName, "Temperature", sizeof(charName) );
        strncpy( charUnit, "K", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
        array = util::vect2double( met.Temp(), m.Ny(), m.Nx(), scalingFactor );
#else
        array = util::vect2float ( met.Temp(), m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

        didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                     (const char*)charSpc,              \
                                     yDim, xDim, outputType,            \
                                     (const char*)charUnit,             \
                                     (const char*)charName );

        /* Saving H2O gaseous concentration */

        strncpy( charSpc, "H2O", sizeof(charSpc) );
        strncpy( charName, "H2O molecular concentration", sizeof(charName) );
        strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
        array = util::vect2double( Data.H2O, m.Ny(), m.Nx(), scalingFactor );
#else
        array = util::vect2float ( Data.H2O, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

        didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                     (const char*)charSpc,              \
                                     yDim, xDim, outputType,            \
                                     (const char*)charUnit,             \
                                     (const char*)charName );


        if ( outputPDF ) {

            /* This might require a lot of disk space. Instead outputting 
             * particle number and volume might be better */

            /* Saving ice aerosol probability density function */

            strncpy( charSpc, "Ice aerosol PDF", sizeof(charSpc) );
            strncpy( charName, "Ice aerosol probability density function (dN/dlogr)", sizeof(charName) );
            strncpy( charUnit, "part/cm^3/log(r)", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
            array = util::vect2double( Data.solidAerosol.pdf, Data.nBin_PA, m.Ny(), m.Nx(), scalingFactor );
#else
            array = util::vect2float ( Data.solidAerosol.pdf, Data.nBin_PA, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

            didSaveSucceed *= fileHandler.addVar3D( currFile, &(array)[0],  \
                                         (const char*)charSpc,              \
                                         binRad, yDim, xDim, outputType,    \
                                         (const char*)charUnit,             \
                                         (const char*)charName );

            /* Saving ice aerosol bin centers.
             * A moving bin structure is adopted in APCEMM. Each grid-cell thus 
             * has a moving bin center.
             * Array: nBin x NY x NX */

            strncpy( charSpc, "Ice aerosol bin centers", sizeof(charSpc) );
            strncpy( charName, "Ice aerosol bin centers", sizeof(charName) );
            strncpy( charUnit, "m^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
            array = util::vect2double( Data.solidAerosol.bin_VCenters, Data.nBin_PA, m.Ny(), m.Nx(), scalingFactor );
#else
            array = util::vect2float ( Data.solidAerosol.bin_VCenters, Data.nBin_PA, m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

            didSaveSucceed *= fileHandler.addVar3D( currFile, &(array)[0],  \
                                         (const char*)charSpc,              \
                                         binRad, yDim, xDim, outputType,    \
                                         (const char*)charUnit,             \
                                         (const char*)charName );

        } else {

            /* This might be a better approach as this requires less disk 
             * space */

            /* Saving ice aerosol particle number 
             * Size: NY x NX */

            strncpy( charSpc, "Ice aerosol particle number", sizeof(charSpc) );
            strncpy( charName, "Ice aerosol particle number", sizeof(charName) );
            strncpy( charUnit, "part/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
            array = util::vect2double( Data.solidAerosol.TotalNumber(), m.Ny(), m.Nx(), scalingFactor );
#else
            array = util::vect2float ( Data.solidAerosol.TotalNumber(), m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

            didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0], \
                                         (const char*)charSpc,             \
                                         yDim, xDim, outputType,           \
                                         (const char*)charUnit,            \
                                         (const char*)charName );

            /* Saving ice aerosol volume
             * Size: NY x NX */

            strncpy( charSpc, "Ice aerosol volume", sizeof(charSpc) );
            strncpy( charName, "Ice aerosol volume", sizeof(charName) );
            strncpy( charUnit, "m^3/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
            array = util::vect2double( Data.solidAerosol.TotalVolume(), m.Ny(), m.Nx(), scalingFactor );
#else
            array = util::vect2float ( Data.solidAerosol.TotalVolume(), m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

            didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0], \
                                         (const char*)charSpc,             \
                                         yDim, xDim, outputType,           \
                                         (const char*)charUnit,            \
                                         (const char*)charName );

            /* Saving horizontal optical depth
             * Size: NY */

            strncpy( charSpc, "Horizontal optical depth", sizeof(charSpc) );
            strncpy( charName, "Horizontally-integrated optical depth", sizeof(charName) );
            strncpy( charUnit, "-", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
            array = util::vect2double( Data.solidAerosol.xOD( m.x_edge() ), m.Ny(), scalingFactor );
#else
            array = util::vect2float ( Data.solidAerosol.xOD( m.x_edge() ), m.Ny(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

            didSaveSucceed *= fileHandler.addVar( currFile, &(array)[0],   \
                                         (const char*)charSpc,             \
                                         yDim, outputType,                 \
                                         (const char*)charUnit,            \
                                         (const char*)charName );

            /* Saving vertical optical depth
             * Size: NX */

            strncpy( charSpc, "Vertical optical depth", sizeof(charSpc) );
            strncpy( charName, "Vertically-integrated optical depth", sizeof(charName) );
            strncpy( charUnit, "-", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
            array = util::vect2double( Data.solidAerosol.yOD( m.y_edge() ), m.Nx(), scalingFactor );
#else
            array = util::vect2float ( Data.solidAerosol.yOD( m.y_edge() ), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

            didSaveSucceed *= fileHandler.addVar( currFile, &(array)[0],   \
                                         (const char*)charSpc,             \
                                         xDim, outputType,                 \
                                         (const char*)charUnit,            \
                                         (const char*)charName );


        }

        delete[] array; array = NULL;

    }

    return SAVE_SUCCESS;

} /* End of Diag_TS_Phys */

/* End of Diag_Mod.cpp */
