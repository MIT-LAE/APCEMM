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

        for ( UInt N = 0; N < NSPECALL; N++ ) {
            for ( UInt i = 0; i < speciesIndices.size(); i++ ) {
            
                if ( speciesIndices[i] - 1 == N ) {

                    strncpy( charSpc, SPC_NAMES[N], sizeof(charSpc) );
                    strncpy( charName, SPC_NAMES[N], sizeof(charName) );
                    strcat(  charName, " molecular concentration" );
                    scalingFactor = 1.0E+00;
                    strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );
            
#if ( SAVE_TO_DOUBLE )
                    array = util::vect2double( Data.Species[N], m.Ny(), m.Nx(), scalingFactor );
#else
                    array = util::vect2float ( Data.Species[N], m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                    #pragma omp critical
                    {
                    didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                                 (const char*)charSpc, yDim, xDim,  \
                                                 outputType, (const char*)charUnit, \
                                                 (const char*)charName );
                    }

                } else {
                    std::cout << " In Diag_Mod for timeseries: Unexpected index: " << speciesIndices[i] << std::endl;
                    std::cout << " Ignoring that index..." << std::endl;
                }
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

        #pragma omp critical
        {
        didSaveSucceed *= fileHandler.addVar( currFile, &(array)[0],    \
                                     (const char*)charSpc,              \
                                     yDim, outputType,                  \
                                     (const char*)charUnit,             \
                                     (const char*)charName );
        }

        /* Saving meteorological temperature */

        strncpy( charSpc, "Temperature", sizeof(charSpc) );
        strncpy( charName, "Temperature", sizeof(charName) );
        strncpy( charUnit, "K", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
        array = util::vect2double( met.Temp(), m.Ny(), m.Nx(), scalingFactor );
#else
        array = util::vect2float ( met.Temp(), m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

        #pragma omp critical
        {
        didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                     (const char*)charSpc,              \
                                     yDim, xDim, outputType,            \
                                     (const char*)charUnit,             \
                                     (const char*)charName );
        }

        /* Saving H2O gaseous concentration */

        strncpy( charSpc, "H2O", sizeof(charSpc) );
        strncpy( charName, "H2O molecular concentration", sizeof(charName) );
        strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
        array = util::vect2double( Data.Species[ind_H2O], m.Ny(), m.Nx(), scalingFactor );
#else
        array = util::vect2float ( Data.Species[ind_H2O], m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

        #pragma omp critical
        {
        didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0],  \
                                     (const char*)charSpc,              \
                                     yDim, xDim, outputType,            \
                                     (const char*)charUnit,             \
                                     (const char*)charName );
        }


        if ( outputPDF == 2 ) {

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

            #pragma omp critical
            {
            didSaveSucceed *= fileHandler.addVar3D( currFile, &(array)[0],  \
                                         (const char*)charSpc,              \
                                         binRad, yDim, xDim, outputType,    \
                                         (const char*)charUnit,             \
                                         (const char*)charName );
            }

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

            #pragma omp critical
            {
            didSaveSucceed *= fileHandler.addVar3D( currFile, &(array)[0],  \
                                         (const char*)charSpc,              \
                                         binRad, yDim, xDim, outputType,    \
                                         (const char*)charUnit,             \
                                         (const char*)charName );
            }

        } else if ( outputPDF == 1 ) {

            /* Saving ice aerosol probability density function */

            strncpy( charSpc, "Aggregated ice aerosol PDF", sizeof(charSpc) );
            strncpy( charName, "Ice aerosol probability density function (dN/dlogr)", sizeof(charName) );
            strncpy( charUnit, "part/m/log(r)", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
            array = util::vect2double( Data.solidAerosol.PDF_Total( m ), Data.nBin_PA, scalingFactor );
#else
            array = util::vect2float ( Data.solidAerosol.PDF_Total( m ), Data.nBin_PA, scalingFactor );
#endif /* SAVE_TO_DOUBLE */

            #pragma omp critical
            {
            didSaveSucceed *= fileHandler.addVar( currFile, &(array)[0],    \
                                         (const char*)charSpc,              \ 
                                         binRad, outputType,                \
                                         (const char*)charUnit,             \
                                         (const char*)charName );
            }

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

            #pragma omp critical
            {
            didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0], \
                                         (const char*)charSpc,             \
                                         yDim, xDim, outputType,           \
                                         (const char*)charUnit,            \
                                         (const char*)charName );
            }

            /* Saving ice aerosol surface area 
             * Size: NY x NX */

            strncpy( charSpc, "Ice aerosol surface area", sizeof(charSpc) );
            strncpy( charName, "Ice aerosol surface area", sizeof(charName) );
            strncpy( charUnit, "m^2/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
            array = util::vect2double( Data.solidAerosol.TotalArea(), m.Ny(), m.Nx(), scalingFactor );
#else
            array = util::vect2float ( Data.solidAerosol.TotalArea(), m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

            #pragma omp critical
            {
            didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0], \
                                         (const char*)charSpc,             \
                                         yDim, xDim, outputType,           \
                                         (const char*)charUnit,            \
                                         (const char*)charName );
            }

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

            #pragma omp critical
            {
            didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0], \
                                         (const char*)charSpc,             \
                                         yDim, xDim, outputType,           \
                                         (const char*)charUnit,            \
                                         (const char*)charName );
            }

            /* Saving ice aerosol effective radius
             * Size: NY x NX */

            strncpy( charSpc, "Effective radius", sizeof(charSpc) );
            strncpy( charName, "Aerosol effective radius", sizeof(charName) );
            strncpy( charUnit, "m", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
            array = util::vect2double( Data.solidAerosol.EffRadius(), m.Ny(), m.Nx(), scalingFactor );
#else
            array = util::vect2float ( Data.solidAerosol.EffRadius(), m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

            #pragma omp critical
            {
            didSaveSucceed *= fileHandler.addVar2D( currFile, &(array)[0], \
                                         (const char*)charSpc,             \
                                         yDim, xDim, outputType,           \
                                         (const char*)charUnit,            \
                                         (const char*)charName );
            }

            /* Saving horizontal optical depth
             * Size: NY */

            strncpy( charSpc, "Horizontal optical depth", sizeof(charSpc) );
            strncpy( charName, "Horizontally-integrated optical depth", sizeof(charName) );
            strncpy( charUnit, "-", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
            array = util::vect2double( Data.solidAerosol.xOD( m.dx() ), m.Ny(), scalingFactor );
#else
            array = util::vect2float ( Data.solidAerosol.xOD( m.dx() ), m.Ny(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

            #pragma omp critical
            {
            didSaveSucceed *= fileHandler.addVar( currFile, &(array)[0],   \
                                         (const char*)charSpc,             \
                                         yDim, outputType,                 \
                                         (const char*)charUnit,            \
                                         (const char*)charName );
            }

            /* Saving vertical optical depth
             * Size: NX */

            strncpy( charSpc, "Vertical optical depth", sizeof(charSpc) );
            strncpy( charName, "Vertically-integrated optical depth", sizeof(charName) );
            strncpy( charUnit, "-", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
            array = util::vect2double( Data.solidAerosol.yOD( m.dy() ), m.Nx(), scalingFactor );
#else
            array = util::vect2float ( Data.solidAerosol.yOD( m.dy() ), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

            #pragma omp critical
            {
            didSaveSucceed *= fileHandler.addVar( currFile, &(array)[0],   \
                                         (const char*)charSpc,             \
                                         xDim, outputType,                 \
                                         (const char*)charUnit,            \
                                         (const char*)charName );
            }

            /* Saving overall size distributionh
             * Size: NX */

            strncpy( charSpc, "Overall size distribution", sizeof(charSpc) );
            strncpy( charName, "Overall size distribution of ice particles", sizeof(charName) );
            strncpy( charUnit, "particles/m", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
            array = util::vect2double( Data.solidAerosol.Overall_Size_Dist( m.areas() ), Data.nBin_PA, scalingFactor );
#else
            array = util::vect2float ( Data.solidAerosol.Overall_Size_Dist( m.areas() ), Data.nBin_PA, scalingFactor );
#endif /* SAVE_TO_DOUBLE */

            #pragma omp critical
            {
            didSaveSucceed *= fileHandler.addVar( currFile, &(array)[0],   \
                                         (const char*)charSpc,             \
                                         binRad, outputType,                 \
                                         (const char*)charUnit,            \
                                         (const char*)charName );
            }

        }

        delete[] array; array = NULL;

    }

    return SAVE_SUCCESS;

} /* End of Diag_TS_Phys */

/* End of Diag_Mod.cpp */
