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

#if ( SAVE_TO_DOUBLE )
    static const NcType& varDataType = ncDouble;
#else
    static const NcType& varDataType = ncFloat;
#endif

bool Diag_TS_Chem( const char* rootName,                     \
                   const std::vector<int> speciesIndices,    \
                   const int hh, const int mm, const int ss, \
                   const Solution& Data, const Mesh& m )
{

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

    // Open the file for writing - replacing anything already there
    NcFile currFile(outFile,NcFile::replace);

    time_t rawtime;
    char buffer[80];
    time( &rawtime );
    strftime(buffer, sizeof(buffer),"%d-%m-%Y %H:%M:%S", localtime(&rawtime));

    // Create time as seconds since start
    const float cur_time = hh+mm/60.0+ss/3600.0;

    // Create dimensions - make time unlimited (record dimension)
    const NcDim xDim       = currFile.addDim( "x", long(m.Nx()) );
    const NcDim yDim       = currFile.addDim( "y", long(m.Ny()) );
    //const NcDim tDim       = currFile.addDim( "Time" );

    // Add the variables corresponding to the dimensions
    NcVar xVar       = currFile.addVar( "x", ncFloat, xDim );
    NcVar yVar       = currFile.addVar( "y", ncFloat, yDim );
    //NcVar tVar       = currFile.addVar( "Time", ncFloat, tDim );

    // Put the data values and attributes into the dimension variables
    xVar.putAtt("units", "m");
    xVar.putAtt("long_name", "Grid cell horizontal centers");
    xVar.putVar(&(m.x())); // Had [0] in the original (i.e. &(m.x())[0])?
    yVar.putAtt("units", "m");
    yVar.putAtt("long_name", "Grid cell vertical centers");
    yVar.putVar(&(m.y()));
    //tVar.putAtt("units", "seconds since simulation start");
    //tVar.putAtt("long_name", "time");
    //tVar.putVar(&(cur_time));

    std::string author = "Thibaud M. Fritz (fritzt@mit.edu)";
    currFile.putAtt( "FileName", outFile );
    currFile.putAtt( "Author", author );
    currFile.putAtt( "Contact", author );
    currFile.putAtt( "Generation Date", buffer );
    currFile.putAtt( "Format", "NetCDF-4" );


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

    // Define vector of dimensions
    std::vector<NcDim> xyDims{ yDim, xDim };
    // Start point and counts for the data being added
    std::vector<size_t> startp2( 2, 0 );
    std::vector<size_t> countp2{ m.Ny(), m.Nx() };

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
                //This call causes a segfault if Chemistry is not turned on but "save species" is because
                //the number of species is set to 1, not NSPECALL, if Chemistry is not turned on.
                array = util::vect2float ( Data.Species[N], m.Ny(), m.Nx(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                #pragma omp critical
                {
                //NcVar var = currFile.addVar( &(array)[0],  \
                //                             (const char*)charSpc, yDim, xDim,  \
                //                             outputType, (const char*)charUnit, \
                //                             (const char*)charName );
                NcVar var = currFile.addVar( charSpc, varDataType, xyDims );
                var.putAtt("units", charUnit );
                var.putAtt("long_name", charName );
                var.putVar( startp2, countp2, &(array)[0] );
                }
                delete[] array;

            } else {
                std::cout << " In Diag_Mod for timeseries: Unexpected index: " << speciesIndices[i] << std::endl;
                std::cout << " Ignoring that index..." << std::endl;
            }
        }

    }
    
    currFile.close();

    return SAVE_SUCCESS;

} /* End of Diag_TS_Chem */

bool Diag_TS_Phys( const char* rootName,                     \
                   const std::vector<int> aerosolIndices,    \
                   const int hh, const int mm, const int ss, \
                   const Solution& Data, const Mesh& m,      \
                   const Meteorology &met,                   \
                   const int outputPDF, \
		   const float partNum_lost, const float iceMass_lost )
{

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

    // Open file and don't worry about overwrite
    NcFile currFile(outFile,NcFile::replace);

    time_t rawtime;
    char buffer[80];
    time( &rawtime );
    strftime(buffer, sizeof(buffer),"%d-%m-%Y %H:%M:%S", localtime(&rawtime));

    // Create time as seconds since start
    const float cur_time = hh+mm/60.0+ss/3600.0;

    // Create dimensions - make time unlimited (record dimension)
    const NcDim xDim       = currFile.addDim( "x", long(m.Nx()) );
    const NcDim yDim       = currFile.addDim( "y", long(m.Ny()) );
    const NcDim binEdgeDim = currFile.addDim( "r_b", Data.nBin_PA+1 );
    const NcDim binRadDim  = currFile.addDim( "r", Data.nBin_PA );
    const NcDim tDim       = currFile.addDim( "t" , 1);

    // Add the variables corresponding to the dimensions
    NcVar xVar       = currFile.addVar( "x", ncFloat, xDim );
    NcVar yVar       = currFile.addVar( "y", ncFloat, yDim );
    NcVar binEdgeVar = currFile.addVar( "r_b", ncFloat, binEdgeDim );
    NcVar binRadVar  = currFile.addVar( "r", ncFloat, binRadDim );
    NcVar tVar       = currFile.addVar( "t", ncFloat, tDim );

    // Put the data values and attributes into the dimension variables
    xVar.putAtt("units", "m");
    xVar.putAtt("long_name", "Grid cell horizontal centers");
    xVar.putVar(&(m.x())[0]); // &(m.x())[0] in original code
    yVar.putAtt("units", "m");
    yVar.putAtt("long_name", "Grid cell vertical centers");
    yVar.putVar(&(m.y())[0]);
    binEdgeVar.putAtt("units", "m");
    binEdgeVar.putAtt("long_name", "ice bin edge radius");
    binEdgeVar.putVar(&(Data.solidAerosol.binEdges())[0]);
    binRadVar.putAtt("units", "m");
    binRadVar.putAtt("long_name", "Ice bin center radius");
    binRadVar.putVar(&(Data.solidAerosol.binCenters())[0]);
    tVar.putAtt("units", "seconds since simulation start");
    tVar.putAtt("long_name", "time");
    tVar.putVar(&(cur_time));

    std::string author = "Thibaud M. Fritz (fritzt@mit.edu)";
    currFile.putAtt( "FileName", outFile );
    currFile.putAtt( "Author", author );
    currFile.putAtt( "Contact", author );
    currFile.putAtt( "Generation Date", buffer );
    currFile.putAtt( "Format", "NetCDF-4" );

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
        NcVar var = currFile.addVar( charSpc, varDataType, yDim );
        var.putAtt("units", charUnit );
        var.putAtt("long_name", charName );
        var.putVar( &(array)[0] );
    }
    delete[] array;

    /* Saving H2O gaseous concentration */

    strncpy( charSpc, "H2O", sizeof(charSpc) );
    strncpy( charName, "H2O molecular concentration", sizeof(charName) );
    strncpy( charUnit, "molec/cm^3", sizeof(charUnit) );

    // For storing 2D data
    std::vector<NcDim> xyDims{ yDim, xDim };
    std::vector<size_t> startp2( 2, 0 );
    std::vector<size_t> countp2{ m.Ny(), m.Nx() };

    // For storing bin radius data
    std::vector<NcDim> xycDims{ binRadDim, yDim, xDim };
    std::vector<size_t> startp3( 3, 0 );
    std::vector<size_t> countp3{ Data.nBin_PA, m.Ny(), m.Nx() };

    #if ( SAVE_TO_DOUBLE )
        array = util::vect2double( Data.Species[ind_H2O], m.Ny(), m.Nx(), scalingFactor );
    #else
        array = util::vect2float ( Data.Species[ind_H2O], m.Ny(), m.Nx(), scalingFactor );
    #endif /* SAVE_TO_DOUBLE */

    #pragma omp critical
    {
        NcVar var = currFile.addVar( charSpc, varDataType, xyDims );
        var.putAtt("units", charUnit );
        var.putAtt("long_name", charName );
        var.putVar( startp2, countp2, &(array)[0] );
    }
    delete[] array;

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
            NcVar var = currFile.addVar( charSpc, varDataType, xyDims );
            var.putAtt("units", charUnit );
            var.putAtt("long_name", charName );
            var.putVar( startp2, countp2, &(array)[0] );
        }
        delete[] array;

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

            NcVar var = currFile.addVar( charSpc, varDataType, xycDims );
            var.putAtt("units", charUnit );
            var.putAtt("long_name", charName );
            var.putVar( startp3, countp3, &(array)[0] );
        }
        delete[] array;

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
            NcVar var = currFile.addVar( charSpc, varDataType, binRadDim );
            var.putAtt("units", charUnit );
            var.putAtt("long_name", charName );
            var.putVar( &(array)[0] );
        }
        delete[] array;

    } else {

        /* This might be a better approach as this requires less disk 
         * space */

        /* Saving particle number that is flux corrected 
         * Size: 1 */

        strncpy( charSpc, "Particle number lost", sizeof(charSpc) );
        strncpy( charName, "Particle number lost", sizeof(charName) );
        strncpy( charUnit, "part/m", sizeof(charUnit) );

        #pragma omp critical
        {
            NcVar var = currFile.addVar( charSpc, varDataType, tDim );
            var.putAtt("units", charUnit );
            var.putAtt("long_name", charName );
            var.putVar( &(partNum_lost) );
        }

        /* Saving ice mass that is flux corrected 
         * Size: 1 */

        strncpy( charSpc, "Ice mass lost", sizeof(charSpc) );
        strncpy( charName, "Ice mass lost", sizeof(charName) );
        strncpy( charUnit, "kg/m", sizeof(charUnit) );

        #pragma omp critical
        {
            NcVar var = currFile.addVar( charSpc, varDataType, tDim );
            var.putAtt("units", charUnit );
            var.putAtt("long_name", charName );
            var.putVar( &(iceMass_lost) );
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
        NcVar var = currFile.addVar( charSpc, varDataType, xyDims );
        var.putAtt("units", charUnit );
        var.putAtt("long_name", charName );
        var.putVar( &(array)[0] );
    }
    delete[] array;


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
            NcVar var = currFile.addVar( charSpc, varDataType, xyDims );
            var.putAtt("units", charUnit );
            var.putAtt("long_name", charName );
            var.putVar( startp2, countp2, &(array)[0] );
        }
        delete[] array;

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
            NcVar var = currFile.addVar( charSpc, varDataType, xyDims );
            var.putAtt("units", charUnit );
            var.putAtt("long_name", charName );
            var.putVar( startp2, countp2, &(array)[0] );
        }
        delete[] array;

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
            NcVar var = currFile.addVar( charSpc, varDataType, xyDims );
            var.putAtt("units", charUnit );
            var.putAtt("long_name", charName );
            var.putVar( startp2, countp2, &(array)[0] );
        }
        delete[] array;

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
            NcVar var = currFile.addVar( charSpc, varDataType, xyDims );
            var.putAtt("units", charUnit );
            var.putAtt("long_name", charName );
            var.putVar( startp2, countp2, &(array)[0] );
        }
        delete[] array;

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
            NcVar var = currFile.addVar( charSpc, varDataType, yDim );
            var.putAtt("units", charUnit );
            var.putAtt("long_name", charName );
            var.putVar( &(array)[0] );
        }
        delete[] array;

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
            NcVar var = currFile.addVar( charSpc, varDataType, xDim );
            var.putAtt("units", charUnit );
            var.putAtt("long_name", charName );
            var.putVar( &(array)[0] );
        }
        delete[] array;

        /* Saving overall size distributionh
         * Size: n Bins */

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
            NcVar var = currFile.addVar( charSpc, varDataType, binRadDim );
            var.putAtt("units", charUnit );
            var.putAtt("long_name", charName );
            var.putVar( &(array)[0] );
        }
        delete[] array;

        /* Saving Total Ice Mass [kg/m]
         * Size: 1 */

        strncpy( charSpc, "Ice Mass", sizeof(charSpc) );
        strncpy( charName, "Total Mass of Ice Crystals", sizeof(charName) );
        strncpy( charUnit, "kg/m", sizeof(charUnit) );

        #pragma omp critical
        {
        NcVar var = currFile.addVar( charSpc, varDataType, tDim );
        var.putAtt("units", charUnit );
        var.putAtt("long_name", charName );
        const float iceMassSum = (float)(Data.solidAerosol.TotalIceMass_sum(m.areas()));
        var.putVar(&iceMassSum);
        }

        /* Saving Num Ice Particles [#/m]
         * Size: 1 */

        strncpy( charSpc, "Number Ice Particles", sizeof(charSpc) );
        strncpy( charName, "Total Number of Ice Particles", sizeof(charName) );
        strncpy( charUnit, "#/m", sizeof(charUnit) );

        #pragma omp critical
        {
        NcVar var = currFile.addVar( charSpc, varDataType, tDim );
        var.putAtt("units", charUnit );
        var.putAtt("long_name", charName );
        const float iceNumPart = (float)(Data.solidAerosol.TotalNumber_sum(m.areas()) );
        var.putVar(&iceNumPart);
        }

        /* Saving Extinction
        Size: NY x NX */

        strncpy( charSpc, "Extinction", sizeof(charSpc) );
        strncpy( charName, "Extinction", sizeof(charName) );
        strncpy( charUnit, "m^-1", sizeof(charUnit) );

    #if ( SAVE_TO_DOUBLE )
            array = util::vect2double( Data.solidAerosol.Extinction(), m.Ny(), m.Nx(), scalingFactor );
    #else
            array = util::vect2float ( Data.solidAerosol.Extinction(), m.Ny(), m.Nx(), scalingFactor );
    #endif /* SAVE_TO_DOUBLE */

        #pragma omp critical
        {
        NcVar var = currFile.addVar( charSpc, varDataType, xyDims );
        var.putAtt("units", charUnit );
        var.putAtt("long_name", charName );
        var.putVar( startp2, countp2, &(array)[0] );
        }
        delete[] array;

        /* Saving IWC
        Size: NY x NX */

        strncpy( charSpc, "IWC", sizeof(charSpc) );
        strncpy( charName, "IWC", sizeof(charName) );
        strncpy( charUnit, "kg/m^3", sizeof(charUnit) );

    #if ( SAVE_TO_DOUBLE )
            array = util::vect2double( Data.solidAerosol.IWC(), m.Ny(), m.Nx(), scalingFactor );
    #else
            array = util::vect2float ( Data.solidAerosol.IWC(), m.Ny(), m.Nx(), scalingFactor );
    #endif /* SAVE_TO_DOUBLE */

        #pragma omp critical
        {
        NcVar var = currFile.addVar( charSpc, varDataType, xyDims );
        var.putAtt("units", charUnit );
        var.putAtt("long_name", charName );
        var.putVar( startp2, countp2, &(array)[0] );
        }
        delete[] array;

    }

        /* Saving RH
        Size: NY x NX */

        strncpy( charSpc, "RHi", sizeof(charSpc) );
        strncpy( charName, "RHi", sizeof(charName) );
        strncpy( charUnit, "%", sizeof(charUnit) );

    #if ( SAVE_TO_DOUBLE )
            array = util::vect2double(physFunc::RHi_Field(Data.Species[ind_H2O], met.Temp(), met.Press()), m.Ny(), m.Nx(), scalingFactor );
    #else
            array = util::vect2float (physFunc::RHi_Field(Data.Species[ind_H2O], met.Temp(), met.Press()), m.Ny(), m.Nx(), scalingFactor );
    #endif /* SAVE_TO_DOUBLE */

        #pragma omp critical
        {
        NcVar var = currFile.addVar( charSpc, varDataType, xyDims );
        var.putAtt("units", charUnit );
        var.putAtt("long_name", charName );
        var.putVar( startp2, countp2, &(array)[0] );
        }
        delete[] array;


    return SAVE_SUCCESS;

} /* End of Diag_TS_Phys */

/* End of Diag_Mod.cpp */
