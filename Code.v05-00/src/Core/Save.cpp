/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Save Program File                                                */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : Save.cpp                                  */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "Core/Save.hpp"

#if ( SAVE_TO_DOUBLE )
    static const NcType& varDataType = ncDouble;
#else
    static const NcType& varDataType = ncFloat;
#endif

namespace output 
{
    void add_ncvar( const NcFile& currFile,           \
                     const std::string varName,       \
                     const std::vector<float>& value, \
                     const NcDim dim,                 \
                     const std::string units,         \
                     const std::string long_name=""   )
    {
        const NcVar var = currFile.addVar( varName, ncFloat, dim );
        var.putVar( &value );
        var.putAtt("units", units);
        var.putAtt("long_name", long_name);
        return;
    }

    void add_ncvar( const NcFile& currFile,            \
                     const std::string varName,        \
                     const std::vector<double>& value, \
                     const NcDim dim,                  \
                     const std::string units,          \
                     const std::string long_name=""    )
    {
        const NcVar var = currFile.addVar( varName, ncDouble, dim );
        var.putVar( &value );
        var.putAtt("units", units);
        var.putAtt("long_name", long_name);
        return;
    }

    void add_ncvar( const NcFile& currFile,         \
                     const std::string varName,      \
                     const float& value,             \
                     const std::string units,        \
                     const std::string long_name=""  )
    {
        const NcVar var = currFile.addVar( varName, ncFloat );
        var.putVar( &value );
        var.putAtt("units", units);
        var.putAtt("long_name", long_name);
        return;
    }

    void add_ncvar( const NcFile& currFile,        \
                     const std::string varName,     \
                     const double& value,           \
                     const std::string units,       \
                     const std::string long_name="" )
    {
        const NcVar var = currFile.addVar( varName, ncDouble );
        var.putVar( &value );
        var.putAtt("units", units);
        var.putAtt("long_name", long_name);
        return;
    }

    void add_ncvar( const NcFile& currFile,         \
                     const std::string varName,      \
                     const int& value,             \
                     const std::string units,        \
                     const std::string long_name=""  )
    {
        const NcVar var = currFile.addVar( varName, ncInt );
        var.putVar( &value );
        var.putAtt("units", units);
        var.putAtt("long_name", long_name);
        return;
    }

    void add_ncvar( const NcFile& currFile,         \
                     const std::string varName,      \
                     const UInt& value,             \
                     const std::string units,        \
                     const std::string long_name=""  )
    {
        const NcVar var = currFile.addVar( varName, ncUint );
        var.putVar( &value );
        var.putAtt("units", units);
        var.putAtt("long_name", long_name);
        return;
    }

    int Write( const char* outFile,                                              \
               const OptInput &Input_Opt,                                        \
               const std::vector<int> speciesIndices,                            \
               const SpeciesArray &ringData, const Ambient &ambientData,         \
               const Cluster &ringCluster, const std::vector<double> &timeArray, \
               const Input &input,                                               \
               const double &airDens, const double &relHumidity_i,               \
               const double &sunRise, const double &sunSet,                      \
               const Vector_3D &plumeRates, const Vector_2D &ambientRates )
    {

        time_t rawtime;
        char buffer[80];

        NcFile currFile(outFile, NcFile::replace);
        time( &rawtime );
        strftime(buffer, sizeof(buffer),"%d-%m-%Y %H:%M:%S", localtime(&rawtime));

        // Define the time mid points
        Vector_1D time_midStep( timeArray.size()-1, 0.0E+00 );

        for ( unsigned int iTime = 0; iTime < timeArray.size() - 1; iTime++ )
            time_midStep[iTime] = 0.5 * (timeArray[iTime] + timeArray[iTime+1]);

        // Two dimensions defined: time at edges of time steps and at mid points
        const NcDim tDim       = currFile.addDim( "t_b", timeArray.size() );
        add_ncvar(currFile, "t_b", timeArray, tDim, "s", "Time at time step ends");
        const NcDim tMidDim    = currFile.addDim( "t", timeArray.size() - 1 );
        add_ncvar(currFile, "t", time_midStep, tMidDim, "s", "Time at time step midpoints");
        //NcVar tVar    = currFile.addVar( "t_b", ncFloat, tDim );
        //NcVar tMidVar = currFile.addVar( "t", ncFloat, tMidDim );
        //tVar.putAtt("units", "s");
        //tVar.putAtt("long_name", "Time at time step ends");
        //tMidVar.putAtt("units", "s");
        //tMidVar.putAtt("long_name", "Time at time step midpoints");

        //tVar.putVar( &timeArray[0] );
        //tMidVar.putVar( &time_midStep[0] );
           
#ifdef RINGS

        const NcDim ringDim = currFile.addDim( "ring", int(ringCluster.getnRing()) );
        NcVar ringVar = currFile.addVar( "ring", ncInt, ringDim );
        ringVar.putAtt("units", "-");
        ringVar.putVar(&((ringCluster.getRingIndex()))[0]);

#endif /* RINGS */

        std::string author = "Thibaud M. Fritz (fritzt@mit.edu)";
        currFile.putAtt( "FileName", outFile );
        currFile.putAtt( "Author", author );
        currFile.putAtt( "Contact", author );
        currFile.putAtt( "Generation Date", buffer );
        currFile.putAtt( "Format", "NetCDF-4" );

        add_ncvar( currFile, "temperature", input.temperature_K(),       "K",          "Ambient temperature" );
        add_ncvar( currFile, "pressure",    input.pressure_Pa() / 100.0, "hPa",        "Ambient pressure");
        add_ncvar( currFile, "airden",      airDens,                     "molec cm-3", "Number density of air");
        add_ncvar( currFile, "rh_w",        input.relHumidity_w(),       "-",          "Ambient relative humidity w.r.t. water");
        add_ncvar( currFile, "rh_i",        relHumidity_i,               "-",          "Ambient relative humidity w.r.t. ice");
        add_ncvar( currFile, "lon",         input.longitude_deg(),       "deg East",   "Longitude");
        add_ncvar( currFile, "lat",         input.latitude_deg(),        "deg North",  "Latitude");
        add_ncvar( currFile, "sunrise",     sunRise,                     "hours",      "Time of local sunrise");
        add_ncvar( currFile, "sunset",      sunSet,                      "hours",      "Time of local sunset");
        add_ncvar( currFile, "shear",       input.shear(),               "s-1",        "Vertical wind shear");
        add_ncvar( currFile, "doy",         input.emissionDOY(),         "-",          "Day of the year");
        add_ncvar( currFile, "time_init",   input.emissionTime(),        "hours",      "Time of emission");
        add_ncvar( currFile, "eiNOx",       input.EI_NOx(),              "g kg-1",     "NOx emissions index on an NO2 mass basis");
        add_ncvar( currFile, "eiCO",        input.EI_CO(),               "g kg-1",     "CO emissions index");
        add_ncvar( currFile, "eiHC",        input.EI_HC(),               "g kg-1",     "HC emissions index");
        add_ncvar( currFile, "eiSO2",       input.EI_SO2(),              "g kg-1",     "SO2 emissions index");
        add_ncvar( currFile, "eiBC",        input.EI_Soot(),             "g kg-1",     "Soot emissions index");
        add_ncvar( currFile, "rBC",         input.sootRad(),             "m",          "Soot radius");
        add_ncvar( currFile, "fuelflow",    input.fuelFlow(),            "kg s-1",     "Fuel flow rate");
        add_ncvar( currFile, "bgNOx",       input.backgNOx(),            "ppbv",       "Background NOx mixing ratio");
        add_ncvar( currFile, "bgHNO3",      input.backgHNO3(),           "ppbv",       "Background HNO3 mixing ratio");
        add_ncvar( currFile, "bgO3",        input.backgO3(),             "ppbv",       "Background O3 mixing ratio");
        add_ncvar( currFile, "bgCO",        input.backgCO(),             "ppbv",       "Background CO mixing ratio");
        add_ncvar( currFile, "bgCH4",       input.backgCH4(),            "ppbv",       "Background CH4 mixing ratio");
        add_ncvar( currFile, "bgSO2",       input.backgSO2(),            "ppbv",       "Background SO2 mixing ratio");

        // Time-varying
        add_ncvar( currFile, "cossza",      ambientData.cosSZA, tMidDim, "-",          "Cosine of the solar zenith angle");

#ifdef RINGS

        //currFile.addVar( &(ringCluster.getRingArea())[0], "Ring Area", ringDim, "float", "m^2", "Ring Area" );
        add_ncvar( currFile, "ring_area", ringCluster.getRingArea(), ringDim, "m2", "Ring area" );

#if ( SAVE_TO_DOUBLE )
        double* spcArray;
        double* ambArray;
        const char* outputType = "double";
#else
        float* spcArray;
        float* ambArray;
        const char* outputType = "float";
#endif /* SAVE_TO_DOUBLE */

        if ( Input_Opt.PL_PL ) {
       
            //TODO: FIX 
            const NcDim *famDim = fileHandler.addDim( currFile, "Family", NFAM );
            std::vector<int> family( NFAM, 0 );

            for ( unsigned int iFam = 0; iFam < NFAM; iFam++ )
                family[iFam] = iFam;

            currFile.addVar( currFile, &family[0], "Family", famDim, "int", "-", "Family");

#if ( SAVE_TO_DOUBLE )
            double* ratesArray;
            ratesArray = util::vect2double( plumeRates, time_midStep.size(), ringCluster.getnRing(), NFAM );
#else
            float* ratesArray;
            ratesArray = util::vect2float ( plumeRates, time_midStep.size(), ringCluster.getnRing(), NFAM );
#endif /* SAVE_TO_DOUBLE */
                
            currFile.addVar( currFile, &(ratesArray)[0], "Rates", tDim_midStep, ringDim, famDim, outputType, "molec/cm^3/s" );

            util::delete1D( ratesArray );
                
                
#if ( SAVE_TO_DOUBLE )
            double* ratesArray_;
            ratesArray_ = util::vect2double( ambientRates, time_midStep.size(), NFAM );
#else
            float* ratesArray_;
            ratesArray_ = util::vect2float ( ambientRates, time_midStep.size(), NFAM );
#endif /* SAVE_TO_DOUBLE */
            
            currFile.addVar( currFile, &(ratesArray_)[0], "Ambient Rates", tDim_midStep, famDim, outputType, "molec/cm^3/s" );

            util::delete1D( ratesArray_ );

        } else {

            if ( Input_Opt.PL_O3 ) {

                const NcDim *famDim = fileHandler.addDim( currFile, "Family", 2 );
                std::vector<int> family( 2, 0 );

                for ( unsigned int iFam = 0; iFam < 2; iFam++ )
                    family[iFam] = iFam;

                currFile.addVar( currFile, &family[0], "Family", famDim, "int", "-", "Family");

#if ( SAVE_TO_DOUBLE )
                double* ratesArray;
                ratesArray = util::vect2double( plumeRates, time_midStep.size(), ringCluster.getnRing(), 2 );
#else
                float* ratesArray;
                ratesArray = util::vect2float ( plumeRates, time_midStep.size(), ringCluster.getnRing(), 2 );
#endif /* SAVE_TO_DOUBLE */
                    
                currFile.addVar( currFile, &(ratesArray)[0], "Rates", tDim, ringDim, famDim, outputType, "molec/cm^3/s" );

                util::delete1D( ratesArray );
                 
#if ( SAVE_TO_DOUBLE )
                double* ratesArray_;
                ratesArray_ = util::vect2double( ambientRates, time_midStep.size(), 2 );
#else
                float* ratesArray_;
                ratesArray_ = util::vect2float ( ambientRates, time_midStep.size(), 2 );
#endif /* SAVE_TO_DOUBLE */
                    
                currFile.addVar( currFile, &(ratesArray_)[0], "Ambient Rates", tDim_midStep, famDim, outputType, "molec/cm^3/s" );

                util::delete1D( ratesArray_ );

            }
        }

#endif /* RINGS */

/* Define conversion factors */
#define TO_PPTH         1.0 / airDens * 1.0E+03 /* Conversion factor from molecule/cm^3 to PPTH */
#define TO_PPM          1.0 / airDens * 1.0E+06 /* Conversion factor from molecule/cm^3 to PPM  */
#define TO_PPB          1.0 / airDens * 1.0E+09 /* Conversion factor from molecule/cm^3 to PPB  */
#define TO_PPT          1.0 / airDens * 1.0E+12 /* Conversion factor from molecule/cm^3 to PPT  */

#ifdef RINGS


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
                    strcat(  charName, " mixing ratio" );
                    scalingFactor = TO_PPB;
                    strncpy( charUnit, "ppb", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                    spcArray = util::vect2double( ringData.Species[N],    timeArray.size(), ringCluster.getnRing(), scalingFactor );
                    ambArray = util::vect2double( ambientData.Species[N], timeArray.size(), scalingFactor );
#else
                    spcArray = util::vect2float ( ringData.Species[N],    timeArray.size(), ringCluster.getnRing(), scalingFactor );
                    ambArray = util::vect2float ( ambientData.Species[N], timeArray.size(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */


                    currFile.addVar( currFile, &(spcArray)[0], (const char*)charSpc, tDim, ringDim, outputType, (const char*)charUnit, (const char*)charName ); 

                    strcat(  charSpc, "_a" );
                    strncpy( charName, SPC_NAMES[N], sizeof(charName) );
                    strcat(  charName, " ambient mixing ratio" );
                    currFile.addVar( currFile, &(ambArray)[0], (const char*)charSpc, tDim, outputType, (const char*)charUnit, (const char*)charName );
                }
            }
        }

        strncpy( charSpc, "NOx", sizeof(charSpc) );
        strncpy( charName, "NOx mixing ratio", sizeof(charName) );
        scalingFactor = TO_PPB;
        strncpy( charUnit, "ppb", sizeof(charUnit) );

        std::vector<std::vector<double>> NOx = util::add2D( ringData.Species[ind_NO], ringData.Species[ind_NO2] );
        std::vector<double> NOx_a = util::add1D( ambientData.Species[ind_NO], ambientData.Species[ind_NO2] );
#if ( SAVE_TO_DOUBLE )
        spcArray = util::vect2double( NOx,   timeArray.size(), ringCluster.getnRing(), scalingFactor );
        ambArray = util::vect2double( NOx_a, timeArray.size(), scalingFactor );
#else
        spcArray = util::vect2float( NOx,   timeArray.size(), ringCluster.getnRing(), scalingFactor );
        ambArray = util::vect2float( NOx_a, timeArray.size(), scalingFactor );
#endif /* SAVE_TO_DOUBLE */
            
        currFile.addVar( currFile, &(spcArray)[0], (const char*)charSpc, tDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

        // This is wrong. NOy shuold include (eg) 2xN2O5
        /*
        strncpy( charSpc, "NOx_a", sizeof(charSpc) );
        strncpy( charName, "NOx ambient mixing ratio", sizeof(charName) );
        currFile.addVar( currFile, &(ambArray)[0], (const char*)charSpc, tDim, outputType, (const char*)charUnit, (const char*)charName );
        
        strncpy( charSpc, "NOy", sizeof(charSpc) );
        strncpy( charName, "NOy mixing ratio", sizeof(charName) );
        scalingFactor = TO_PPB;
        strncpy( charUnit, "ppb", sizeof(charUnit) );

        std::vector<std::vector<double>> NOy = util::add2D( util::add2D( util::add2D( util::add2D( util::add2D( util::add2D( util::add2D( util::add2D( util::add2D( util::add2D( util::add2D( util::add2D( util::add2D( util::add2D( util::add2D( util::add2D(util::add2D( util::add2D( util::add2D( util::add2D( ringData.Species[ind_NO], ringData.Species[ind_NO2] ), ringData.Species[ind_NO3] ), ringData.Species[ind_HNO2] ), ringData.Species[ind_HNO3] ), ringData.Species[ind_HNO4] ), ringData.Species[ind_N2O5] ), ringData.Species[ind_N2O5] ), ringData.Species[ind_PAN] ), ringData.Species[ind_BrNO2] ), ringData.Species[ind_BrNO3] ), ringData.Species[ind_ClNO2] ), ringData.Species[ind_ClNO3] ), ringData.Species[ind_PPN] ), ringData.Species[ind_N] ), ringData.Species[ind_MPN] ), ringData.Species[ind_PROPNN] ), ringData.Species[ind_PRPN] ), ringData.Species[ind_R4N1] ), ringData.Species[ind_PRN1] ), ringData.Species[ind_R4N2] );
        std::vector<double> NOy_a = util::add1D( util::add1D( util::add1D( util::add1D( util::add1D( util::add1D( util::add1D( util::add1D( util::add1D( util::add1D( util::add1D( util::add1D( util::add1D( util::add1D( util::add1D( util::add1D( util::add1D( util::add1D( util::add1D( util::add1D( ambientData.Species[ind_NO], ambientData.Species[ind_NO2] ), ambientData.Species[ind_NO3] ), ambientData.Species[ind_HNO2] ), ambientData.Species[ind_HNO3] ), ambientData.Species[ind_HNO4] ), ambientData.Species[ind_N2O5] ), ambientData.Species[ind_N2O5] ), ambientData.Species[ind_PAN] ), ambientData.Species[ind_BrNO2] ), ambientData.Species[ind_BrNO3] ), ambientData.Species[ind_ClNO2] ), ambientData.Species[ind_ClNO3] ), ambientData.Species[ind_PPN] ), ambientData.Species[ind_N] ), ambientData.Species[ind_MPN] ), ambientData.Species[ind_PROPNN] ), ambientData.Species[ind_PRPN] ), ambientData.Species[ind_R4N1] ), ambientData.Species[ind_PRN1] ), ambientData.Species[ind_R4N2] );
#if ( SAVE_TO_DOUBLE )
        spcArray = util::vect2double( NOy,   timeArray.size(), ringCluster.getnRing(), scalingFactor );
        ambArray = util::vect2double( NOy_a, timeArray.size(), scalingFactor );
#else
        spcArray = util::vect2float( NOy,   timeArray.size(), ringCluster.getnRing(), scalingFactor );
        ambArray = util::vect2float( NOy_a, timeArray.size(), scalingFactor );
#endif // SAVE_TO_DOUBLE
            
        currFile.addVar( currFile, &(spcArray)[0], (const char*)charSpc, tDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

        strncpy( charSpc, "NOy_a", sizeof(charSpc) );
        strncpy( charName, "NOy ambient mixing ratio", sizeof(charName) );
        currFile.addVar( currFile, &(ambArray)[0], (const char*)charSpc, tDim, outputType, (const char*)charUnit, (const char*)charName );

        strncpy( charSpc, "NOy_N2O", sizeof(charSpc) );
        strncpy( charName, "NOy + N2O mixing ratio", sizeof(charName) );
        scalingFactor = TO_PPB;
        strncpy( charUnit, "ppb", sizeof(charUnit) );

        std::vector<std::vector<double>> NOy_N2O = util::add2D( NOy, ringData.Species[ind_N2O] );
        std::vector<double> NOy_N2O_a = util::add1D( NOy_a, ambientData.Species[ind_N2O] );
#if ( SAVE_TO_DOUBLE )
        spcArray = util::vect2double( NOy_N2O, timeArray.size(), ringCluster.getnRing(), scalingFactor );
        ambArray = util::vect2double( NOy_N2O_a, timeArray.size(), scalingFactor );
#else
        spcArray = util::vect2float( NOy_N2O, timeArray.size(), ringCluster.getnRing(), scalingFactor );
        ambArray = util::vect2float( NOy_N2O_a, timeArray.size(), scalingFactor );
#endif // SAVE_TO_DOUBLE
            
        currFile.addVar( currFile, &(spcArray)[0], (const char*)charSpc, tDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

        strncpy( charSpc, "NOy_N2Oa", sizeof(charSpc) );
        strncpy( charName, "NOy + N2O ambient mixing ratio", sizeof(charName) );
        currFile.addVar( currFile, &(ambArray)[0], (const char*)charSpc, tDim, outputType, (const char*)charUnit, (const char*)charName );
            */
        
        util::delete1D( spcArray );
        util::delete1D( ambArray );

#endif /* RINGS */

        currFile.close();

        return SAVE_SUCCESS;

    } /* End of Write */

    int Write_MicroPhys( const char* outputFile, \
                         const std::vector<std::vector<std::vector<std::vector<double>>>> &output_MicroPhys, \
                         const std::vector<double> &timeArray, const std::vector<double> &binCenters, \
                         const std::vector<double> &horizDim, const std::vector<double> &verticDim, \
                         const double temperature_K, const double pressure_Pa, const double lapseRate,
                         const double relHumidity_w, const double relHumidity_i )
    {
        const char* currFileName( outputFile );

        int didSaveSucceed = 1;
        time_t rawtime;
        char buffer[80];

        NcFile currFile(outputFile,NcFile::replace);
//        std::cout << "\n Starting saving to netCDF (file name: " << fileHandler.getFileName() <<  ") \n";
        time( &rawtime );
        strftime(buffer, sizeof(buffer),"%d-%m-%Y %H:%M:%S", localtime(&rawtime));

        const NcDim xDim       = currFile.addDim( "x", long(horizDim.size()) );
        const NcDim yDim       = currFile.addDim( "y", long(verticDim.size()) );
        const NcDim binRadDim  = currFile.addDim( "r", binCenters.size() );
        const NcDim tDim       = currFile.addDim( "t", timeArray.size() );

        add_ncvar(currFile, "x", horizDim,   xDim,      "m", "Grid cell x centers");
        add_ncvar(currFile, "y", verticDim,  yDim,      "m", "Grid cell y centers");
        add_ncvar(currFile, "r", binCenters, binRadDim, "m", "Grid cell y centers");
        add_ncvar(currFile, "t", timeArray,  tDim,      "s", "Time since simulation start");
        //const NcVar xVar       = currFile.addVar( "x", ncFloat, xDim );
        //const NcVar yVar       = currFile.addVar( "y", ncFloat, yDim );
        //const NcVar binRadVar  = currFile.addVar( "r", ncFloat, binRadDim );
        //const NcVar tVar       = currFile.addVar( "t", ncFloat, tDim );

        //// Put the data values and attributes into the dimension variables
        //xVar.putAtt("units", "m");
        //xVar.putAtt("long_name", "Grid cell horizontal centers");
        //xVar.putVar(&horizDim[0]);
        //yVar.putAtt("units", "m");
        //yVar.putAtt("long_name", "Grid cell vertical centers");
        //yVar.putVar(&verticDim[0]);
        //binRadVar.putAtt("units", "m");
        //binRadVar.putAtt("long_name", "Ice bin center radius");
        //binRadVar.putVar(&binCenters[0]);
        //tVar.putAtt("units", "seconds since simulation start");
        //tVar.putAtt("long_name", "time");
        //tVar.putVar(&timeArray[0]);
        //const NcDim *tDim = fileHandler.addDim( currFile, "Time", long(timeArray.size()) );
        //currFile.addVar( currFile, &timeArray[0], "Time", tDim, "float", "s", "Time");
        //const NcDim *binDim = fileHandler.addDim( currFile, "Bin Centers", long(binCenters.size()) );
        //currFile.addVar( currFile, &binCenters[0], "bin Centers", binDim, "float", "m", "Bin Centers");
        //const NcDim *XDim = fileHandler.addDim( currFile, "Horizontal dimension", long(horizDim.size()) );
        //currFile.addVar( currFile, &horizDim[0], "X coordinate", XDim, "float", "m", "X coordinate of the grid cell centers");
        //const NcDim *YDim = fileHandler.addDim( currFile, "Vertical dimension", long(verticDim.size()) );
        //currFile.addVar( currFile, &verticDim[0], "Y coordinate", YDim, "float", "m", "Y coordinate of the grid cell centers");
        
        std::string author = "Thibaud M. Fritz (fritzt@mit.edu)";
        currFile.putAtt( "FileName", outputFile );
        currFile.putAtt( "Author", author );
        currFile.putAtt( "Contact", author );
        currFile.putAtt( "Generation Date", buffer );
        currFile.putAtt( "Format", "NetCDF-4" );

        // Add fixed values
        add_ncvar( currFile, "temperature", temperature_K,              "K",          "Ambient temperature" );
        add_ncvar( currFile, "pressure",    pressure_Pa / 100.0,        "hPa",        "Ambient pressure");
        add_ncvar( currFile, "lapse",       lapseRate,                  "K/km",       "Vertical lapse rate");
        add_ncvar( currFile, "rh_w",        relHumidity_w,              "-",          "Ambient relative humidity w.r.t. water");
        add_ncvar( currFile, "rh_i",        relHumidity_i,              "-",          "Ambient relative humidity w.r.t. water");
        //currFile.addVar( currFile, &temperature_K, "Temperature", 1, "float", "K"   , "Ambient Temperature" );
        //currFile.addVar( currFile, &pressure_Pa  , "Pressure"   , 1, "float", "Pa"  , "Ambient Pressure" );
        //currFile.addVar( currFile, &lapseRate    , "Lapse Rate" , 1, "float", "K/km", "Ambient temperature lapse rate" );
        //currFile.addVar( currFile, &relHumidity_w, "RHW"        , 1, "float", "-"   , "Ambient Rel. Humidity w.r.t water" );
        //currFile.addVar( currFile, &relHumidity_i, "RHI"        , 1, "float", "-"   , "Ambient Rel. Humidity w.r.t ice" );

#if ( SAVE_TO_DOUBLE )
        double* aerArray;
        const char* outputType = "double";
#else
        float* aerArray;
        const char* outputType = "float";
#endif /* SAVE_TO_DOUBLE */
    
        char charSpc[30];
        char charName[50];
        char charUnit[20];
            
        strncpy( charSpc, "Aerosol", sizeof(charSpc) );
        strncpy( charName, "Aerosol number concentration", sizeof(charName) );
        strncpy( charUnit, "#/cm^3", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
        aerArray = util::vect2double( output_MicroPhys , timeArray.size(), binCenters.size(), verticDim.size(), horizDim.size(), 1.0 );
#else
        aerArray = util::vect2float( output_MicroPhys, timeArray.size(), binCenters.size(), verticDim.size(), horizDim.size(), 1.0 );
#endif /* SAVE_TO_DOUBLE */
         
        std::vector<NcDim> dimVector{ tDim, binRadDim, yDim, xDim };
        std::vector<size_t> startp4( 4, 0 );
        std::vector<size_t> countp4{ timeArray.size(), binCenters.size(), verticDim.size(), horizDim.size() };
        #pragma omp critical
        {
        NcVar var = currFile.addVar( charSpc, varDataType, dimVector );
        var.putAtt("units", charUnit);
        var.putAtt("long_name", charName);
        var.putVar( startp4, countp4, &(aerArray)[0] );
        }

        util::delete1D( aerArray );

        currFile.close();

        return SAVE_SUCCESS;

    } /* End of Write_MicroPhys */

    int Write_Adjoint( const char* outputFile,                                    \
                       const std::vector<int> speciesIndices,                     \
                       const SpeciesArray &ringData, const Ambient &ambientData,  \
                       const Ambient &adjointData,                                \
                       const std::vector<double> &ringArea, const double totArea, \
                       const std::vector<double> &timeArray,                      \
                       const Input &input,                                        \
                       const double &airDens, const double &relHumidity_i )
    {

        std::cout << " Starting to save adjoint output!" << std::endl;

        time_t rawtime;
        char buffer[80];

        NcFile currFile( outputFile, NcFile::replace );
//        std::cout << "\n Starting saving to netCDF (file name: " << fileHandler.getFileName() <<  ") \n";
        time( &rawtime );
        strftime(buffer, sizeof(buffer),"%d-%m-%Y %H:%M:%S", localtime(&rawtime));
       
        // Calculate time step midpoints 
        Vector_1D time_midStep( timeArray.size()-1, 0.0E+00 );
        for ( unsigned int iTime = 0; iTime < timeArray.size() - 1; iTime++ )
            time_midStep[iTime] = 0.5 * (timeArray[iTime] + timeArray[iTime+1]);

        const NcDim tDim       = currFile.addDim( "t_b", timeArray.size() );
        add_ncvar(currFile, "t_b", timeArray, tDim, "s", "Time at time step ends");
        const NcDim tMidDim    = currFile.addDim( "t", timeArray.size() - 1 );
        add_ncvar(currFile, "t", time_midStep, tMidDim, "s", "Time at time step midpoints");
        
        currFile.putAtt( "FileName", outputFile );
        currFile.putAtt( "Author", input.author() );
        currFile.putAtt( "Contact", input.author() );
        currFile.putAtt( "Generation Date", buffer );
        currFile.putAtt( "Format", "NetCDF-4" );

        add_ncvar( currFile, "temperature", input.temperature_K(),       "K",          "Ambient temperature" );
        add_ncvar( currFile, "pressure",    input.pressure_Pa() / 100.0, "hPa",        "Ambient pressure");
        add_ncvar( currFile, "airden",      airDens,                     "molec cm-3", "Number density of air");
        add_ncvar( currFile, "rh_w",        input.relHumidity_w(),       "-",          "Ambient relative humidity w.r.t. water");
        add_ncvar( currFile, "rh_i",        relHumidity_i,               "-",          "Ambient relative humidity w.r.t. ice");
        add_ncvar( currFile, "lon",         input.longitude_deg(),       "deg East",   "Longitude");
        add_ncvar( currFile, "lat",         input.latitude_deg(),        "deg North",  "Latitude");
        //add_ncvar( currFile, "sunrise",     sunRise,                     "hours",      "Time of local sunrise");
        //add_ncvar( currFile, "sunset",      sunSet,                      "hours",      "Time of local sunset");
        add_ncvar( currFile, "shear",       input.shear(),               "s-1",        "Vertical wind shear");
        add_ncvar( currFile, "doy",         input.emissionDOY(),         "-",          "Day of the year");
        add_ncvar( currFile, "time_init",   input.emissionTime(),        "hours",      "Time of emission");
        add_ncvar( currFile, "eiNOx",       input.EI_NOx(),              "g kg-1",     "NOx emissions index on an NO2 mass basis");
        add_ncvar( currFile, "eiCO",        input.EI_CO(),               "g kg-1",     "CO emissions index");
        add_ncvar( currFile, "eiHC",        input.EI_HC(),               "g kg-1",     "HC emissions index");
        add_ncvar( currFile, "eiSO2",       input.EI_SO2(),              "g kg-1",     "SO2 emissions index");
        add_ncvar( currFile, "eiBC",        input.EI_Soot(),             "g kg-1",     "Soot emissions index");
        add_ncvar( currFile, "rBC",         input.sootRad(),             "m",          "Soot radius");
        add_ncvar( currFile, "fuelflow",    input.fuelFlow(),            "kg s-1",     "Fuel flow rate");
        add_ncvar( currFile, "bgNOx",       input.backgNOx(),            "ppbv",       "Background NOx mixing ratio");
        add_ncvar( currFile, "bgHNO3",      input.backgHNO3(),           "ppbv",       "Background HNO3 mixing ratio");
        add_ncvar( currFile, "bgO3",        input.backgO3(),             "ppbv",       "Background O3 mixing ratio");
        add_ncvar( currFile, "bgCO",        input.backgCO(),             "ppbv",       "Background CO mixing ratio");
        add_ncvar( currFile, "bgCH4",       input.backgCH4(),            "ppbv",       "Background CH4 mixing ratio");
        add_ncvar( currFile, "bgSO2",       input.backgSO2(),            "ppbv",       "Background SO2 mixing ratio");

        // Time-varying
        add_ncvar( currFile, "cossza",      ambientData.cosSZA, tMidDim, "-",          "Cosine of the solar zenith angle");

        /* Define conversion factors */
        #define TO_PPTH         1.0 / airDens * 1.0E+03 /* Conversion factor from molecule/cm^3 to PPTH */
        #define TO_PPM          1.0 / airDens * 1.0E+06 /* Conversion factor from molecule/cm^3 to PPM  */
        #define TO_PPB          1.0 / airDens * 1.0E+09 /* Conversion factor from molecule/cm^3 to PPB  */
        #define TO_PPT          1.0 / airDens * 1.0E+12 /* Conversion factor from molecule/cm^3 to PPT  */

        #if ( SAVE_TO_DOUBLE )
            double* spcArray;
            double* ambArray;
            double* adjArray;
            const char* outputType = "double";
        #else
            float* spcArray;
            float* ambArray;
            float* adjArray;
            const char* outputType = "float";
        #endif /* SAVE_TO_DOUBLE */

        /* ringAverage contains all ring-averaged concentrations in molec/cm^3 and is indexed as:
         * ringAverage[iTime][iSpecies] */
        std::vector<std::vector<double>> ringAverage = ringData.RingAverage( ringArea, totArea );

        const unsigned int NT = timeArray.size();
        std::vector<double> plumeData( NT, 0.0E+00 );

        /* Start saving species ... */

        /* Define conversion factor */
        double scalingFactor;
        char charSpc[30];
        char charName[50];
        char charUnit[20];

        for ( UInt N = 0; N < NSPECALL; N++ ) {
            for ( UInt i = 0; i < speciesIndices.size(); i++ ) {
                if ( speciesIndices[i] - 1 == N ) {
                    scalingFactor = TO_PPB;
                    strncpy( charUnit, "ppb", sizeof(charUnit) );
           
                    for ( unsigned int iNt = 0; iNt < NT; iNt++ )
                        plumeData[iNt] = ringAverage[iNt][N];
          
                    // Takes a vector and returns a pointer to the first element after
                    // multiplying by a scaling factor. Now recasting to a vector
                    // afterwards. 
#if ( SAVE_TO_DOUBLE )
                    spcArray = util::vect2double( plumeData             , NT, scalingFactor );
                    adjArray = util::vect2double( adjointData.Species[N], NT, scalingFactor );
                    ambArray = util::vect2double( ambientData.Species[N], NT, scalingFactor );
                    std::vector<double> spcVec{ spcArray, spcArray + NT };
                    std::vector<double> adjVec{ adjArray, adjArray + NT };
                    std::vector<double> ambVec{ ambArray, ambArray + NT };
#else
                    spcArray = util::vect2float ( plumeData             , NT, scalingFactor );
                    adjArray = util::vect2float ( adjointData.Species[N], NT, scalingFactor );
                    ambArray = util::vect2float ( ambientData.Species[N], NT, scalingFactor );
                    std::vector<float> spcVec{ spcArray, spcArray + NT };
                    std::vector<float> adjVec{ adjArray, adjArray + NT };
                    std::vector<float> ambVec{ ambArray, ambArray + NT };
#endif /* SAVE_TO_DOUBLE */


                    strncpy( charSpc, SPC_NAMES[N], sizeof(charSpc) );
                    strcat(  charSpc, "_Plume" );
                    strncpy( charName, SPC_NAMES[N], sizeof(charName) );
                    strcat(  charSpc, " plume-averaged mixing ratio" );

                    //currFile.addVar( currFile, &(spcArray)[0], (const char*)charSpc, tDim, outputType, (const char*)charUnit, (const char*)charName );
                    //{
                    //    NcVar var = currFile.addVar( charSpc, varDataType, tDim );
                    //    var.putAtt("units", charUnit);
                    //    var.putAtt("long_name", charName);
                    //    var.putVar(&(spcArray)[0]);
                    //}
                    //add_ncvar(currFile, charSpc, spcArray, tDim, charUnit, charName);
                    add_ncvar(currFile, charSpc, spcVec, tDim, charUnit, charName);

                    strncpy( charSpc, SPC_NAMES[N], sizeof(charSpc) );
                    strcat(  charSpc, "_Adjoint" );
                    strncpy( charName, SPC_NAMES[N], sizeof(charName) );
                    strcat(  charSpc, " optimized mixing ratio" );

                    //currFile.addVar( currFile, &(adjArray)[0], (const char*)charSpc, tDim, outputType, (const char*)charUnit, (const char*)charName ); 
                    //{
                    //    NcVar var = currFile.addVar( charSpc, varDataType, tDim );
                    //    var.putAtt("units", charUnit);
                    //    var.putAtt("long_name", charName);
                    //    var.putVar(&(spcArray)[0]);
                    //}
                    add_ncvar(currFile, charSpc, adjVec, tDim, charUnit, charName);

                    strncpy( charSpc, SPC_NAMES[N], sizeof(charSpc) );
                    strcat(  charSpc, "_Ambient" );
                    strncpy( charName, SPC_NAMES[N], sizeof(charName) );
                    strcat(  charSpc, " ambient mixing ratio" );

                    //currFile.addVar( currFile, &(ambArray)[0], (const char*)charSpc, tDim, outputType, (const char*)charUnit, (const char*)charName ); 
                    add_ncvar(currFile, charSpc, ambVec, tDim, charUnit, charName);
                }
            }
        }

        delete[] spcArray; spcArray = NULL;
        delete[] ambArray; ambArray = NULL;
        delete[] adjArray; adjArray = NULL;

        scalingFactor = TO_PPB;
        strncpy( charUnit, "ppb", sizeof(charUnit) );

        for ( unsigned int iNt = 0; iNt < NT; iNt++ )
            plumeData[iNt] = ringAverage[iNt][ind_NO] + ringAverage[iNt][ind_NO2];

        std::vector<double> ambientNOx = util::add1D( ambientData.Species[ind_NO], ambientData.Species[ind_NO2] );
        std::vector<double> adjointNOx = util::add1D( adjointData.Species[ind_NO], adjointData.Species[ind_NO2] );

        {
#if ( SAVE_TO_DOUBLE )
        spcArray = util::vect2double( plumeData , NT, scalingFactor );
        adjArray = util::vect2double( adjointNOx, NT, scalingFactor );
        ambArray = util::vect2double( ambientNOx, NT, scalingFactor );
        std::vector<double> spcVec{ spcArray, spcArray + NT };
        std::vector<double> adjVec{ adjArray, adjArray + NT };
        std::vector<double> ambVec{ ambArray, ambArray + NT };
#else
        spcArray = util::vect2float ( plumeData , NT, scalingFactor );
        adjArray = util::vect2float ( adjointNOx, NT, scalingFactor );
        ambArray = util::vect2float ( ambientNOx, NT, scalingFactor );
        std::vector<float> spcVec{ spcArray, spcArray + NT };
        std::vector<float> adjVec{ adjArray, adjArray + NT };
        std::vector<float> ambVec{ ambArray, ambArray + NT };
#endif /* SAVE_TO_DOUBLE */


        strncpy( charSpc, "NOx_Plume", sizeof(charSpc) );
        strncpy( charName, "NOx plume-averaged mixing ratio", sizeof(charName) );
        //currFile.addVar( currFile, &(spcArray)[0], (const char*)charSpc, tDim, outputType, (const char*)charUnit, (const char*)charName ); 
        add_ncvar(currFile, charSpc, spcVec, tDim, charUnit, charName);

        strncpy( charSpc, "NOx_Adjoint", sizeof(charSpc) );
        strncpy( charName, "NOx optimized mixing ratio", sizeof(charName) );
        //currFile.addVar( currFile, &(adjArray)[0], (const char*)charSpc, tDim, outputType, (const char*)charUnit, (const char*)charName ); 
        add_ncvar(currFile, charSpc, adjVec, tDim, charUnit, charName);

        strncpy( charSpc, "NOx_Ambient", sizeof(charSpc) );
        strncpy( charName, "NOx ambient mixing ratio", sizeof(charName) );
        //currFile.addVar( currFile, &(ambArray)[0], (const char*)charSpc, tDim, outputType, (const char*)charUnit, (const char*)charName ); 
        add_ncvar(currFile, charSpc, ambVec, tDim, charUnit, charName);
        }

        delete[] spcArray; spcArray = NULL;
        delete[] ambArray; ambArray = NULL;
        delete[] adjArray; adjArray = NULL;

        //TODO: Fix this - NOy should have 2xN2O5 for example
        /*
        scalingFactor = TO_PPB;
        strncpy( charUnit, "ppb", sizeof(charUnit) );

        for ( unsigned int iNt = 0; iNt < NT; iNt++ )
            plumeData[iNt] = ringAverage[iNt][ind_NO] + ringAverage[iNt][ind_NO2] + ringAverage[iNt][ind_NO3] + ringAverage[iNt][ind_HNO2] + ringAverage[iNt][ind_HNO3] + ringAverage[iNt][ind_HNO4] + 2 * ringAverage[iNt][ind_N2O5] + ringAverage[iNt][ind_PAN] + ringAverage[iNt][ind_BrNO2] + ringAverage[iNt][ind_BrNO3] + ringAverage[iNt][ind_ClNO2] + ringAverage[iNt][ind_ClNO3] + ringAverage[iNt][ind_PPN] + ringAverage[iNt][ind_N] + ringAverage[iNt][ind_MPN] + ringAverage[iNt][ind_PROPNN] + ringAverage[iNt][ind_PRPN] + ringAverage[iNt][ind_R4N1] + ringAverage[iNt][ind_PRN1] + ringAverage[iNt][ind_R4N2]; 

        std::vector<double> ambientNOy = util::add1D( util::add1D( util::add1D( util::add1D( util::add1D( util::add1D( util::add1D( util::add1D( util::add1D( util::add1D( util::add1D( util::add1D( util::add1D( util::add1D( util::add1D( util::add1D(util::add1D( util::add1D( util::add1D( util::add1D( ambientData.Species[ind_NO], ambientData.Species[ind_NO2] ), ambientData.Species[ind_NO3] ), ambientData.Species[ind_HNO2] ), ambientData.Species[ind_HNO3] ), ambientData.Species[ind_HNO4] ), ambientData.Species[ind_N2O5] ), ambientData.Species[ind_N2O5] ), ambientData.Species[ind_PAN] ), ambientData.Species[ind_BrNO2] ), ambientData.Species[ind_BrNO3] ), ambientData.Species[ind_ClNO2] ), ambientData.Species[ind_ClNO3] ), ambientData.Species[ind_PPN] ), ambientData.Species[ind_N] ), ambientData.Species[ind_MPN] ), ambientData.Species[ind_PROPNN] ), ambientData.Species[ind_PRPN] ), ambientData.Species[ind_R4N1] ), ambientData.Species[ind_PRN1] ), ambientData.Species[ind_R4N2] );

#if ( SAVE_TO_DOUBLE )
            spcArray = util::vect2double( plumeData , NT, scalingFactor );
            ambArray = util::vect2double( ambientNOy, NT, scalingFactor );
#else
            spcArray = util::vect2float ( plumeData , NT, scalingFactor );
            ambArray = util::vect2float ( ambientNOy, NT, scalingFactor );
#endif // SAVE_TO_DOUBLE

        strncpy( charSpc, "NOy_Plume", sizeof(charSpc) );
        strncpy( charName, "NOy plume-averaged mixing ratio", sizeof(charName) );
        currFile.addVar( currFile, &(spcArray)[0], (const char*)charSpc, tDim, outputType, (const char*)charUnit, (const char*)charName ); 

        strncpy( charSpc, "NOy_Ambient", sizeof(charSpc) );
        strncpy( charName, "NOy ambient mixing ratio", sizeof(charName) );
        currFile.addVar( currFile, &(ambArray)[0], (const char*)charSpc, tDim, outputType, (const char*)charUnit, (const char*)charName ); 

        delete[] spcArray; spcArray = NULL;
        delete[] ambArray; ambArray = NULL;
        delete[] adjArray; adjArray = NULL;
        */

        return SAVE_SUCCESS;

    } /* End of Write_Adjoint */
    
    int Write_Box( const char* outputFile,               \
                   const std::vector<int> speciesIndices,\
                   const Ambient &boxData,               \
                   const std::vector<double> &timeArray, \
                   const Input &input,                   \
                   const double &airDens,                \
                   const double &relHumidity_i )
    {

        std::cout << " Starting to save box output!" << std::endl;

        time_t rawtime;
        char buffer[80];

        NcFile currFile(outputFile,NcFile::replace);
//            std::cout << "\n Starting saving to netCDF (file name: " << fileHandler.getFileName() <<  ") \n";
        time( &rawtime );
        strftime(buffer, sizeof(buffer),"%d-%m-%Y %H:%M:%S", localtime(&rawtime));

        Vector_1D time_midStep( timeArray.size()-1, 0.0E+00 );
        for ( unsigned int iTime = 0; iTime < timeArray.size() - 1; iTime++ )
            time_midStep[iTime] = 0.5 * (timeArray[iTime] + timeArray[iTime+1]);

        const NcDim tDim       = currFile.addDim( "t_b", timeArray.size() );
        add_ncvar(currFile, "t_b", timeArray, tDim, "s", "Time at time step ends");
        const NcDim tMidDim    = currFile.addDim( "t", timeArray.size() - 1 );
        add_ncvar(currFile, "t", time_midStep, tMidDim, "s", "Time at time step midpoints");
        
        std::string author = "Thibaud M. Fritz (fritzt@mit.edu)";
        currFile.putAtt( "FileName", outputFile );
        currFile.putAtt( "Author", author );
        currFile.putAtt( "Contact", author );
        currFile.putAtt( "Generation Date", buffer );
        currFile.putAtt( "Format", "NetCDF-4" );

        add_ncvar( currFile, "temperature", input.temperature_K(),       "K",          "Ambient temperature" );
        add_ncvar( currFile, "pressure",    input.pressure_Pa() / 100.0, "hPa",        "Ambient pressure");
        add_ncvar( currFile, "airden",      airDens,                     "molec cm-3", "Number density of air");
        add_ncvar( currFile, "rh_w",        input.relHumidity_w(),       "-",          "Ambient relative humidity w.r.t. water");
        add_ncvar( currFile, "rh_i",        relHumidity_i,               "-",          "Ambient relative humidity w.r.t. ice");
        add_ncvar( currFile, "lon",         input.longitude_deg(),       "deg East",   "Longitude");
        add_ncvar( currFile, "lat",         input.latitude_deg(),        "deg North",  "Latitude");
        //add_ncvar( currFile, "sunrise",     sunRise,                     "hours",      "Time of local sunrise");
        //add_ncvar( currFile, "sunset",      sunSet,                      "hours",      "Time of local sunset");
        add_ncvar( currFile, "shear",       input.shear(),               "s-1",        "Vertical wind shear");
        add_ncvar( currFile, "doy",         input.emissionDOY(),         "-",          "Day of the year");
        add_ncvar( currFile, "time_init",   input.emissionTime(),        "hours",      "Time of emission");
        add_ncvar( currFile, "eiNOx",       input.EI_NOx(),              "g kg-1",     "NOx emissions index on an NO2 mass basis");
        add_ncvar( currFile, "eiCO",        input.EI_CO(),               "g kg-1",     "CO emissions index");
        add_ncvar( currFile, "eiHC",        input.EI_HC(),               "g kg-1",     "HC emissions index");
        add_ncvar( currFile, "eiSO2",       input.EI_SO2(),              "g kg-1",     "SO2 emissions index");
        add_ncvar( currFile, "eiBC",        input.EI_Soot(),             "g kg-1",     "Soot emissions index");
        add_ncvar( currFile, "rBC",         input.sootRad(),             "m",          "Soot radius");
        add_ncvar( currFile, "fuelflow",    input.fuelFlow(),            "kg s-1",     "Fuel flow rate");
        add_ncvar( currFile, "bgNOx",       input.backgNOx(),            "ppbv",       "Background NOx mixing ratio");
        add_ncvar( currFile, "bgHNO3",      input.backgHNO3(),           "ppbv",       "Background HNO3 mixing ratio");
        add_ncvar( currFile, "bgO3",        input.backgO3(),             "ppbv",       "Background O3 mixing ratio");
        add_ncvar( currFile, "bgCO",        input.backgCO(),             "ppbv",       "Background CO mixing ratio");
        add_ncvar( currFile, "bgCH4",       input.backgCH4(),            "ppbv",       "Background CH4 mixing ratio");
        add_ncvar( currFile, "bgSO2",       input.backgSO2(),            "ppbv",       "Background SO2 mixing ratio");

        // Time-varying
        add_ncvar( currFile, "cossza",      boxData.cosSZA, tMidDim,     "-",          "Cosine of the solar zenith angle");

        /* Define conversion factors */
        #define TO_PPTH         1.0 / airDens * 1.0E+03 /* Conversion factor from molecule/cm^3 to PPTH */
        #define TO_PPM          1.0 / airDens * 1.0E+06 /* Conversion factor from molecule/cm^3 to PPM  */
        #define TO_PPB          1.0 / airDens * 1.0E+09 /* Conversion factor from molecule/cm^3 to PPB  */
        #define TO_PPT          1.0 / airDens * 1.0E+12 /* Conversion factor from molecule/cm^3 to PPT  */

        #if ( SAVE_TO_DOUBLE )
        double* ambArray;
        const char* outputType = "double";
        #else
        float* ambArray;
        const char* outputType = "float";
        #endif /* SAVE_TO_DOUBLE */

        const unsigned int NT = timeArray.size();

        /* Start saving species ... */

        /* Define conversion factor */
        double scalingFactor;
        char charSpc[30];
        char charName[50];
        char charUnit[20];
           
        for ( UInt N = 0; N < NSPECALL; N++ ) {
            for ( UInt i = 0; i < speciesIndices.size(); i++ ) {
                if ( speciesIndices[i] - 1 == N ) {
                    scalingFactor = TO_PPB;
                    strncpy( charUnit, "ppb", sizeof(charUnit) );

#if ( SAVE_TO_DOUBLE )
                    ambArray = util::vect2double( boxData.Species[N] , NT, scalingFactor );
                    std::vector<double> ambVec{ ambArray, ambArray + NT };
#else
                    ambArray = util::vect2float ( boxData.Species[N] , NT, scalingFactor );
                    std::vector<float> ambVec{ ambArray, ambArray + NT };
#endif /* SAVE_TO_DOUBLE */

                    strncpy( charSpc, SPC_NAMES[N], sizeof(charSpc) );
                    strcat(  charSpc, "_Ambient" );
                    strncpy( charName, SPC_NAMES[N], sizeof(charName) );
                    strcat(  charName, " ambient mixing ratio" );
        
                    add_ncvar( currFile, charSpc, ambVec, tDim, charUnit, charName );  
                    //currFile.addVar( currFile, &(ambArray)[0], (const char*)charSpc, tDim, outputType, (const char*)charUnit, (const char*)charName ); 
                }
            }
        }

        scalingFactor = TO_PPB;
        strncpy( charUnit, "ppb", sizeof(charUnit) );
        
        std::vector<double> ambientNOx = util::add1D( boxData.Species[ind_NO], boxData.Species[ind_NO2] );
#if ( SAVE_TO_DOUBLE )
        ambArray = util::vect2double( ambientNOx, NT, scalingFactor );
        std::vector<double> ambVec{ ambArray, ambArray + NT };
#else
        ambArray = util::vect2float ( ambientNOx, NT, scalingFactor );
        std::vector<float> ambVec{ ambArray, ambArray + NT };
#endif /* SAVE_TO_DOUBLE */


        strncpy( charSpc, "NOx_Ambient", sizeof(charSpc) );
        strncpy( charName, "NOx ambient mixing ratio", sizeof(charName) );
        //currFile.addVar( currFile, &(ambArray)[0], (const char*)charSpc, tDim, outputType, (const char*)charUnit, (const char*)charName ); 
        add_ncvar( currFile, charSpc, ambVec, tDim, charUnit, charName );  
         

        delete[] ambArray; ambArray = NULL;

        currFile.close();

        return SAVE_SUCCESS;

    } /* End of Write_Box */


}

/* End of Save.cpp */
