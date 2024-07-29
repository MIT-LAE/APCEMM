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
#include <format>

#include "Core/Diag_Mod.hpp"
namespace Diag {

    static const NcType& varDataType = ncFloat;

    bool Diag_TS_Chem( const char* rootName,                     \
                    const std::vector<int> speciesIndices,    \
                    const int hh, const int mm, const int ss, \
                    const Solution& Data, const Mesh& m )
    {

        std::filesystem::path rootPath( rootName );
        std::string fileName = rootPath.filename().generic_string();

        replace_hhmmss(fileName, hh, mm, ss);
        
        fileName = std::filesystem::path(rootPath.parent_path() / fileName).generic_string();
        const char* outFile = fileName.c_str();

        // Open the file for writing - replacing anything already there
        NcFile currFile(outFile,NcFile::replace);

        time_t rawtime;
        char buffer[80];
        time( &rawtime );
        strftime(buffer, sizeof(buffer),"%d-%m-%Y %H:%M:%S", localtime(&rawtime));

        // Create time as seconds since start
        // const float cur_time = hh+mm/60.0+ss/3600.0;

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
        // const char* outputType = "float";
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
                    /*
                    NcVar var = currFile.addVar( &(array)[0],  \
                                                (const char*)charSpc, yDim, xDim,  \
                                                outputType, (const char*)charUnit, \
                                                (const char*)charName );
                    */
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

    void replace_hhmmss(string& fileName, int hh, int mm, int ss) {
        size_t start_pos;
        bool found;

        /* Replacing "hh" with hour number since start */
        start_pos = 0; 
        found = 0;

        // Make zero padded string representations with C++20 format
        auto hh_string = std::format("{:02}", hh);
        auto mm_string = std::format("{:02}", mm);
        auto ss_string = std::format("{:02}", ss);

        while ( (start_pos = fileName.find("hh", start_pos)) != std::string::npos ) {
            fileName.replace(start_pos, 2, hh_string);
            start_pos += 2;
            found = 1;
        }

        if ( found == 0 ) {
            std::cout << " Diagnostic filename must be of the form *hhmmss.nc or *hhmm.nc. Aborting!" << std::endl;
            std::cout << " filename: " << fileName << std::endl;
            throw std::runtime_error("Did not find pattern hh");
        }

        /* Replacing "mm" with minute number since start */
        start_pos = 0; 
        found = 0;
        while ( (start_pos = fileName.find("mm", start_pos)) != std::string::npos ) {
            fileName.replace(start_pos, 2, mm_string);
            start_pos += 2;
            found = 1;
        }
        
        if ( found == 0 ) {
            std::cout << " Diagnostic filename must be of the form *hhmmss.nc or *hhmm.nc. Aborting!" << std::endl;
            std::cout << " filename: " << fileName << std::endl;
            throw std::runtime_error("Did not find pattern mm");
        }

        /* Replacing "ss" with minute number since start if possible */
        start_pos = 0;
        found = 0;
        while ( (start_pos = fileName.find("ss", start_pos)) != std::string::npos ) {
            fileName.replace(start_pos, 2, ss_string);
            start_pos += 2;
        }
    }

    void add0DVar(NcFile& currFile, const float toSave, const NcDim& dim, const string& name, const string& desc, const string& units) {

        NcVar var = currFile.addVar( name, varDataType, dim );
        var.putAtt("units", units );
        var.putAtt("long_name", desc );
        var.putVar(&toSave);
    }

    void add1DVar(NcFile& currFile, const Vector_1D& toSave, const NcDim& dim, const string& name, const string& desc, const string& units) {
        if(toSave.size() != dim.getSize()) {
            throw std::runtime_error("Save failed! NcDim dimension size and array sizes don't match! Variable name: " + name );
        }
        const double scalingFactor = 1;
        float* array = util::vect2float ( toSave, toSave.size(), scalingFactor );
        NcVar var = currFile.addVar( name, varDataType, dim );
        var.putAtt("units", units );
        var.putAtt("long_name", desc );
        var.putVar(array);
        delete[] array;

    }
    void add2DVar(NcFile& currFile, const Vector_2D& toSave, const vector<NcDim> dims, const string& name, const string& desc, const string& units) {
        if((toSave.size() != dims[0].getSize()) || (toSave[0].size() != dims[1].getSize())) {
            throw std::runtime_error("Save failed! NcDim dimension size and array sizes don't match! Variable name: " + name );
        }
        const double scalingFactor = 1;
        float* array = util::vect2float (toSave, toSave.size(), toSave[0].size(), scalingFactor );
        NcVar var = currFile.addVar( name, varDataType, dims );
        var.putAtt("units", units );
        var.putAtt("long_name", desc );
        var.putVar( array );

        delete[] array;
    }

    void Diag_TS_Phys( const char* rootName,
                    const int hh, const int mm, const int ss,
                    const AIM::Grid_Aerosol& iceAer, const Vector_2D& H2O,
                    const Vector_1D& xCoord, const Vector_1D& yCoord,
                    const Vector_1D& xEdges, const Vector_1D& yEdges,
                    const Meteorology &met)
    {   
        long unsigned int nBin = iceAer.getNBin();
        long unsigned int nx = xCoord.size();
        long unsigned int ny = yCoord.size();

        Vector_2D areas = VectorUtils::cellAreas(xEdges, yEdges);
        Vector_1D dx_vec(nx, xCoord[1] - xCoord[0]);
        Vector_1D dy_vec(ny, yCoord[1] - yCoord[0]);
        std::filesystem::path rootPath( rootName );
        std::string fileName = rootPath.filename().generic_string();

        replace_hhmmss(fileName, hh, mm, ss);

        fileName = std::filesystem::path(rootPath.parent_path() / fileName).generic_string();
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
        const NcDim xDim       = currFile.addDim( "x", nx );
        const NcDim yDim       = currFile.addDim( "y", ny );
        const NcDim binEdgeDim = currFile.addDim( "r_b", nBin +1 );
        const NcDim binRadDim  = currFile.addDim( "r", nBin );
        const NcDim tDim       = currFile.addDim( "t" , 1);

        // For storing 2D data
        std::vector<NcDim> xyDims{ yDim, xDim };

        // For storing bin radius data
        std::vector<NcDim> xycDims{ binRadDim, yDim, xDim };

        // Add the variables corresponding to the dimensions
        NcVar xVar       = currFile.addVar( "x", ncFloat, xDim );
        NcVar yVar       = currFile.addVar( "y", ncFloat, yDim );
        NcVar binEdgeVar = currFile.addVar( "r_e", ncFloat, binEdgeDim );
        NcVar binRadVar  = currFile.addVar( "r", ncFloat, binRadDim );
        NcVar tVar       = currFile.addVar( "t", ncFloat, tDim );

        // Put the data values and attributes into the dimension variables
        xVar.putAtt("units", "m");
        xVar.putAtt("long_name", "Grid cell horizontal centers");
        xVar.putVar(xCoord.data());
        yVar.putAtt("units", "m");
        yVar.putAtt("long_name", "Grid cell vertical centers");
        yVar.putVar(&(yCoord)[0]);
        binEdgeVar.putAtt("units", "m");
        binEdgeVar.putAtt("long_name", "ice bin edge radius");
        binEdgeVar.putVar(&(iceAer.getBinEdges())[0]);
        binRadVar.putAtt("units", "m");
        binRadVar.putAtt("long_name", "Ice bin center radius");
        binRadVar.putVar(&(iceAer.getBinCenters())[0]);
        tVar.putAtt("units", "seconds since simulation start");
        tVar.putAtt("long_name", "time");
        tVar.putVar(&(cur_time));

        std::string author = "Thibaud M. Fritz (fritzt@mit.edu)";
        currFile.putAtt( "FileName", outFile );
        currFile.putAtt( "Author", author );
        currFile.putAtt( "Contact", author );
        currFile.putAtt( "Generation Date", buffer );
        currFile.putAtt( "Format", "NetCDF-4" );

        /* Output met */

        /* Check if evolving met. If not only save met in one file */

        /* Saving meteorological pressure */
        add1DVar(currFile, met.Press(), yDim, "Pressure", "Pressure", "Pa");

        //Save Altitude
        add1DVar(currFile, met.Altitude(), yDim, "Altitude", "Altitude", "m");

        /* Saving H2O gaseous concentration */
        add2DVar(currFile, H2O, xyDims, "H2O", "H2O molecular concentration", "molec / cm^3");

        /* Saving meteorological temperature */
        add2DVar(currFile, met.Temp(), xyDims, "Temperature", "Temperature", "K");

        /* Saving ice aerosol particle number */
        add2DVar(currFile, iceAer.TotalNumber(), xyDims, "Ice aerosol particle number", "Ice aerosol particle number concentration", "# / cm^3");

        // /* Saving ice aerosol surface area 
        add2DVar(currFile, iceAer.TotalArea(), xyDims, "Ice aerosol surface area", "Ice aerosol surface area", "m^2 / cm^3");
    
        /* Saving ice aerosol volume */
        add2DVar(currFile, iceAer.TotalVolume(), xyDims, "Ice aerosol volume", "Ice aerosol volume", "m^3 / cm^3");

        /* Saving ice aerosol effective radius */
        add2DVar(currFile, iceAer.EffRadius(), xyDims, "Effective radius", "Ice aerosol effective radius", "m");

        /* Saving horizontal optical depth */
        add1DVar(currFile, iceAer.xOD(dx_vec), yDim, "Horizontal optical depth", "Horizontally-integrated optical depth", "-");

        /* Saving vertical optical depth */
        add1DVar(currFile, iceAer.yOD(dy_vec), xDim, "Vertical optical depth", "Vertically-integrated optical depth", "-");
    
        /* Saving overall size distribution */ 
        add1DVar(currFile, iceAer.Overall_Size_Dist(areas), binRadDim, "Overall size distribution", "Overall size distribution of ice particles", "part / m");

        /* Saving Total Ice Mass [kg/m] */
        add0DVar(currFile, iceAer.TotalIceMass_sum(areas), tDim, "Ice Mass", "Total Mass of Ice Crystals of Cross Section", "kg / m");

        /* Saving Num Ice Particles [#/m] */
        add0DVar(currFile, iceAer.TotalNumber_sum(areas), tDim, "Number Ice Particles", "Total Number of Ice Particles of Cross Section", "# / m");

        /* Saving Extinction [-/m]*/
        add2DVar(currFile, iceAer.Extinction(), xyDims, "Extinction", "Extinction", "m^-1");

        /* Saving IWC */
        add2DVar(currFile, iceAer.IWC(), xyDims, "IWC", "Ice Water Content", "kg / m^3");

        /* Saving RHi */
        add2DVar(currFile, physFunc::RHi_Field(H2O, met.Temp(), met.Press()), xyDims, "RHi", "Relative Humidity w.r.t. Ice", "%");

        //Contrail width, depth, and integrated OD
        add0DVar(currFile, iceAer.extinctionWidth(xCoord), tDim, "width", "Contrail Extinction-Defined Width", "m");
        add0DVar(currFile, iceAer.extinctionDepth(yCoord), tDim, "depth", "Contrail Extinction-Defined Depth", "m");
        add0DVar(currFile, iceAer.intYOD(dx_vec, dy_vec), tDim, "intOD", "Integrated Vertical Optical Depth", "m");
    } /* End of Diag_TS_Phys */

}

/* End of Diag_Mod.cpp */
