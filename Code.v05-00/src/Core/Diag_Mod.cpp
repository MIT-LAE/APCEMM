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
#ifndef FMT_HEADER_ONLY
#define FMT_HEADER_ONLY
#endif

#include <filesystem>
#include <fmt/core.h>
#include "KPP/KPP_Parameters.h"
#include "Util/PhysFunction.hpp"
#include "Core/Util.hpp"
#include "Core/Diag_Mod.hpp"

namespace Diag {

    namespace {
        bool storePSD = false;
    }

    void set_storePSD(bool val){
        storePSD = val;
    }

    static const NcType& varDataType = ncFloat;

    void replace_hhmmss(string& fileName, int hh, int mm, int ss) {
        size_t start_pos;
        bool found;

        /* Replacing "hh" with hour number since start */
        start_pos = 0; 
        found = 0;

        // Make zero padded string representations with fmt
        auto hh_string = fmt::format("{:02}", hh);
        auto mm_string = fmt::format("{:02}", mm);
        auto ss_string = fmt::format("{:02}", ss);

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
        float* array = util::vect2float ( toSave, toSave.size() );
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
        float* array = util::vect2float (toSave, toSave.size(), toSave[0].size() );
        NcVar var = currFile.addVar( name, varDataType, dims );
        var.putAtt("units", units );
        var.putAtt("long_name", desc );
        var.putVar( array );

        delete[] array;
    }

    void add3DVar(NcFile& currFile, const Vector_3D& toSave, const vector<NcDim> dims, const string& name, const string& desc, const string& units) {
        if((toSave.size() != dims[0].getSize()) || (toSave[0].size() != dims[1].getSize()) || (toSave[0][0].size() != dims[2].getSize())) {
            throw std::runtime_error("Save failed! NcDim dimension size and array sizes don't match! Variable name: " + name );
        }
        float* array = util::vect2float (toSave, toSave.size(), toSave[0].size(), toSave[0][0].size() );
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

        if (storePSD)
        {
            add3DVar(currFile, iceAer.Number(), xycDims, "n_aer", "Ice aerosol particle number concentration by radius", "# / cm^3");
        }
    } /* End of Diag_TS_Phys */

}

/* End of Diag_Mod.cpp */
