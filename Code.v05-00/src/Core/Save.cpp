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

namespace output 
{

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

        const bool doWrite = 1;
        const bool doRead = 1;
        const bool overWrite = 1;

        int didSaveSucceed = 1;
        time_t rawtime;
        char buffer[80];

        FileHandler fileHandler( outFile, doWrite, doRead, overWrite );
        NcFile currFile = fileHandler.openFile();
        if ( !fileHandler.isFileOpen() ) {
            std::cout << " File " << outFile << " didn't open!" << "\n";
            return SAVE_FAILURE;
        } else {
//            std::cout << "\n Starting saving to netCDF (file name: " << fileHandler.getFileName() <<  ") \n";
            time( &rawtime );
            strftime(buffer, sizeof(buffer),"%d-%m-%Y %H:%M:%S", localtime(&rawtime));

            const NcDim *timeDim = fileHandler.addDim( currFile, "Time", timeArray.size() );
            didSaveSucceed *= fileHandler.addVar( currFile, &timeArray[0], "Time", timeDim, "float", "s", "Time");
            
            const NcDim *timeDim_midStep = fileHandler.addDim( currFile, "Time_mid", timeArray.size() - 1 );
            Vector_1D time_midStep( timeArray.size()-1, 0.0E+00 );

            for ( unsigned int iTime = 0; iTime < timeArray.size() - 1; iTime++ )
                time_midStep[iTime] = 0.5 * (timeArray[iTime] + timeArray[iTime+1]);

            didSaveSucceed *= fileHandler.addVar( currFile, &time_midStep[0], "Time_mid", timeDim_midStep, "float", "s", "Time at mid time-step");
           
#ifdef RINGS

            const NcDim *ringDim = fileHandler.addDim( currFile, "ring", long(ringCluster.getnRing()) );
            didSaveSucceed *= fileHandler.addVar( currFile, &((ringCluster.getRingIndex()))[0], "ring index", ringDim, "short", "-", "Ring Indices");

#endif /* RINGS */

            didSaveSucceed *= fileHandler.addAtt( currFile, "FileName", fileHandler.getFileName() );
            didSaveSucceed *= fileHandler.addAtt( currFile, "Author", "Thibaud M. Fritz (fritzt@mit.edu)" );
            didSaveSucceed *= fileHandler.addAtt( currFile, "Contact", "Thibaud M. Fritz (fritzt@mit.edu)" );
            didSaveSucceed *= fileHandler.addAtt( currFile, "GenerationDate", buffer );
            didSaveSucceed *= fileHandler.addAtt( currFile, "Format", "NetCDF-4" );

            double value = 0.0E+00;
            value = input.temperature_K();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Temperature", 1, "float", "K"  , "Ambient Temperature" );
            value = input.pressure_Pa() / 100.0;
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Pressure"   , 1, "float", "hPa", "Ambient Pressure" );
            value = airDens;
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Air Density", 1, "float", "molecule / cm ^ 3", "Molecular density" );
            value = input.relHumidity_w();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "RHW"        , 1, "float", "-"  , "Ambient Rel. Humidity w.r.t water" );
            value = relHumidity_i;
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "RHI"        , 1, "float", "-"  , "Ambient Rel. Humidity w.r.t ice" );
            value = input.longitude_deg();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Longitude"  , 1, "float", "deg", "Longitude" );
            value = input.latitude_deg();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Latitude"   , 1, "float", "deg", "Latitude" );
            value = sunRise;
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Sun Rise"   , 1, "float", "hrs", "Local sun rise" );
            value = sunSet;
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Sun Set"    , 1, "float", "hrs", "Local sun set" );
            value = input.shear();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Shear"          , 1, "float", "1/s", "Ambient wind shear" );
            value = input.emissionDOY();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Emission Day"   , 1, "int"  , "-"  , "Emission day" );
            value = input.emissionTime();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Emission Time"  , 1, "float", "hr"  , "Emission time" );
            value = input.EI_NOx();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "NOx EI"         , 1, "float", "g/kg_fuel"  , "NOx Emission index" );
            value = input.EI_CO();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "CO EI"          , 1, "float", "g/kg_fuel"  , "CO Emission index" );
            value = input.EI_HC();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "HC EI"          , 1, "float", "g/kg_fuel"  , "HC Emission index" );
            value = input.EI_SO2();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "SO2 EI"          , 1, "float", "g/kg_fuel"  , "SO2 Emission index" );
            value = input.EI_Soot();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Soot EI"        , 1, "float", "g/kg_fuel"  , "Soot Emission index" );
            value = input.sootRad();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Soot Radius"    , 1, "float", "g/kg_fuel"  , "Soot radius" );
            value = input.fuelFlow();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Fuel flow"      , 1, "float", "kg/s"  , "Engine fuel flow" );
            value = input.backgNOx();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Background NOx" , 1, "float", "ppb"  , "Background NOx mixing ratio" );
            value = input.backgHNO3();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Background HNO3", 1, "float", "ppb"  , "Background HNO3 mixing ratio" );
            value = input.backgO3();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Background O3"  , 1, "float", "ppb"  , "Background O3 mixing ratio" );
            value = input.backgCO();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Background CO"  , 1, "float", "ppb"  , "Background CO mixing ratio" );
            value = input.backgCH4();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Background CH4" , 1, "float", "ppb"  , "Background CH4 mixing ratio" );
            value = input.backgSO2();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Background SO2" , 1, "float", "ppb"  , "Background SO2 mixing ratio" );

            didSaveSucceed *= fileHandler.addVar( currFile, &(ambientData.cosSZA)[0], "CSZA", timeDim_midStep, "float", "-", "Cosine of the solar zenith angle" );

#ifdef RINGS

            didSaveSucceed *= fileHandler.addVar( currFile, &(ringCluster.getRingArea())[0], "Ring Area", ringDim, "float", "m^2", "Ring Area" );

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
            
                const NcDim *famDim = fileHandler.addDim( currFile, "Family", NFAM );
                std::vector<int> family( NFAM, 0 );

                for ( unsigned int iFam = 0; iFam < NFAM; iFam++ )
                    family[iFam] = iFam;

                didSaveSucceed *= fileHandler.addVar( currFile, &family[0], "Family", famDim, "int", "-", "Family");

#if ( SAVE_TO_DOUBLE )
                double* ratesArray;
                ratesArray = util::vect2double( plumeRates, time_midStep.size(), ringCluster.getnRing(), NFAM );
#else
                float* ratesArray;
                ratesArray = util::vect2float ( plumeRates, time_midStep.size(), ringCluster.getnRing(), NFAM );
#endif /* SAVE_TO_DOUBLE */
                
                didSaveSucceed *= fileHandler.addVar3D( currFile, &(ratesArray)[0], "Rates", timeDim_midStep, ringDim, famDim, outputType, "molec/cm^3/s" );

                util::delete1D( ratesArray );
                
                
#if ( SAVE_TO_DOUBLE )
                double* ratesArray_;
                ratesArray_ = util::vect2double( ambientRates, time_midStep.size(), NFAM );
#else
                float* ratesArray_;
                ratesArray_ = util::vect2float ( ambientRates, time_midStep.size(), NFAM );
#endif /* SAVE_TO_DOUBLE */
                
                didSaveSucceed *= fileHandler.addVar2D( currFile, &(ratesArray_)[0], "Ambient Rates", timeDim_midStep, famDim, outputType, "molec/cm^3/s" );

                util::delete1D( ratesArray_ );

            } else {

                if ( Input_Opt.PL_O3 ) {

                    const NcDim *famDim = fileHandler.addDim( currFile, "Family", 2 );
                    std::vector<int> family( 2, 0 );

                    for ( unsigned int iFam = 0; iFam < 2; iFam++ )
                        family[iFam] = iFam;

                    didSaveSucceed *= fileHandler.addVar( currFile, &family[0], "Family", famDim, "int", "-", "Family");

#if ( SAVE_TO_DOUBLE )
                    double* ratesArray;
                    ratesArray = util::vect2double( plumeRates, time_midStep.size(), ringCluster.getnRing(), 2 );
#else
                    float* ratesArray;
                    ratesArray = util::vect2float ( plumeRates, time_midStep.size(), ringCluster.getnRing(), 2 );
#endif /* SAVE_TO_DOUBLE */
                    
                    didSaveSucceed *= fileHandler.addVar3D( currFile, &(ratesArray)[0], "Rates", timeDim, ringDim, famDim, outputType, "molec/cm^3/s" );

                    util::delete1D( ratesArray );
                     
#if ( SAVE_TO_DOUBLE )
                    double* ratesArray_;
                    ratesArray_ = util::vect2double( ambientRates, time_midStep.size(), 2 );
#else
                    float* ratesArray_;
                    ratesArray_ = util::vect2float ( ambientRates, time_midStep.size(), 2 );
#endif /* SAVE_TO_DOUBLE */
                    
                    didSaveSucceed *= fileHandler.addVar2D( currFile, &(ratesArray_)[0], "Ambient Rates", timeDim_midStep, famDim, outputType, "molec/cm^3/s" );

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


                            didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName ); 

                            strcat(  charSpc, "_a" );
                            strncpy( charName, SPC_NAMES[N], sizeof(charName) );
                            strcat(  charName, " ambient mixing ratio" );
                            didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );
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
                
                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "NOx_a", sizeof(charSpc) );
                strncpy( charName, "NOx ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );
            
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
#endif /* SAVE_TO_DOUBLE */
                
                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "NOy_a", sizeof(charSpc) );
                strncpy( charName, "NOy ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

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
#endif /* SAVE_TO_DOUBLE */
                
                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "NOy_N2Oa", sizeof(charSpc) );
                strncpy( charName, "NOy + N2O ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );
            
            util::delete1D( spcArray );
            util::delete1D( ambArray );

#endif /* RINGS */

            if ( didSaveSucceed == NC_SUCCESS ) {
//                std::cout << " Done saving to netCDF!" << "\n";
            } else if ( didSaveSucceed != NC_SUCCESS ) {
                std::cout << "Error occured in saving data: didSaveSucceed: " << didSaveSucceed << "\n";
                return SAVE_FAILURE;
            }

            fileHandler.closeFile( currFile );
            if ( fileHandler.isFileOpen() ) {
                std::cout << "File " << outFile << " didn't close properly!" << "\n";
                return SAVE_FAILURE;
            }
        }

        return SAVE_SUCCESS;

    } /* End of Write */

    int Write_MicroPhys( const char* outputFile, \
                         const std::vector<std::vector<std::vector<std::vector<double>>>> &output_MicroPhys, \
                         const std::vector<double> &timeArray, const std::vector<double> &binCenters, \
                         const std::vector<double> &horizDim, const std::vector<double> &verticDim, \
                         const double temperature_K, const double pressure_Pa, const double lapseRate,
                         const double relHumidity_w, const double relHumidity_i )
    {
        const bool doWrite = 1;
        const bool doRead = 1;
        const bool overWrite = 1;
        const char* currFileName( outputFile );

        int didSaveSucceed = 1;
        time_t rawtime;
        char buffer[80];

        FileHandler fileHandler( currFileName, doWrite, doRead, overWrite );
        NcFile currFile = fileHandler.openFile();
        if ( !fileHandler.isFileOpen() ) {
            std::cout << " File " << currFileName << " didn't open!" << "\n";
            return SAVE_FAILURE;
        } else {
//            std::cout << "\n Starting saving to netCDF (file name: " << fileHandler.getFileName() <<  ") \n";
            time( &rawtime );
            strftime(buffer, sizeof(buffer),"%d-%m-%Y %H:%M:%S", localtime(&rawtime));

            const NcDim *timeDim = fileHandler.addDim( currFile, "Time", long(timeArray.size()) );
            didSaveSucceed *= fileHandler.addVar( currFile, &timeArray[0], "Time", timeDim, "float", "s", "Time");
            const NcDim *binDim = fileHandler.addDim( currFile, "Bin Centers", long(binCenters.size()) );
            didSaveSucceed *= fileHandler.addVar( currFile, &binCenters[0], "bin Centers", binDim, "float", "m", "Bin Centers");
            const NcDim *XDim = fileHandler.addDim( currFile, "Horizontal dimension", long(horizDim.size()) );
            didSaveSucceed *= fileHandler.addVar( currFile, &horizDim[0], "X coordinate", XDim, "float", "m", "X coordinate of the grid cell centers");
            const NcDim *YDim = fileHandler.addDim( currFile, "Vertical dimension", long(verticDim.size()) );
            didSaveSucceed *= fileHandler.addVar( currFile, &verticDim[0], "Y coordinate", YDim, "float", "m", "Y coordinate of the grid cell centers");
            
            didSaveSucceed *= fileHandler.addAtt( currFile, "FileName", fileHandler.getFileName() );
            didSaveSucceed *= fileHandler.addAtt( currFile, "Author", "Thibaud M. Fritz (fritzt@mit.edu)" );
            didSaveSucceed *= fileHandler.addAtt( currFile, "Contact", "Thibaud M. Fritz (fritzt@mit.edu)" );
            didSaveSucceed *= fileHandler.addAtt( currFile, "GenerationDate", buffer );
            didSaveSucceed *= fileHandler.addAtt( currFile, "Format", "NetCDF-4" );

            didSaveSucceed *= fileHandler.addConst( currFile, &temperature_K, "Temperature", 1, "float", "K"   , "Ambient Temperature" );
            didSaveSucceed *= fileHandler.addConst( currFile, &pressure_Pa  , "Pressure"   , 1, "float", "Pa"  , "Ambient Pressure" );
            didSaveSucceed *= fileHandler.addConst( currFile, &lapseRate    , "Lapse Rate" , 1, "float", "K/km", "Ambient temperature lapse rate" );
            didSaveSucceed *= fileHandler.addConst( currFile, &relHumidity_w, "RHW"        , 1, "float", "-"   , "Ambient Rel. Humidity w.r.t water" );
            didSaveSucceed *= fileHandler.addConst( currFile, &relHumidity_i, "RHI"        , 1, "float", "-"   , "Ambient Rel. Humidity w.r.t ice" );

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
             
            didSaveSucceed *= fileHandler.addVar4D( currFile, &(aerArray)[0], (const char*)charSpc, timeDim, binDim, YDim, XDim, outputType, (const char*)charUnit, (const char*)charName );

            util::delete1D( aerArray );

            if ( didSaveSucceed == NC_SUCCESS ) {
//                std::cout << " Done saving to netCDF!" << "\n";
            } else if ( didSaveSucceed != NC_SUCCESS ) {
                std::cout << "Error occured in saving data: didSaveSucceed: " << didSaveSucceed << "\n";
                return SAVE_FAILURE;
            }

            fileHandler.closeFile( currFile );
            if ( fileHandler.isFileOpen() ) {
                std::cout << "File " << currFileName << " didn't close properly!" << "\n";
                return SAVE_FAILURE;
            }

        }
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
        const bool doWrite = 1;
        const bool doRead = 1;
        const bool overWrite = 1;

        int didSaveSucceed = 1;
        time_t rawtime;
        char buffer[80];

        FileHandler fileHandler( outputFile, doWrite, doRead, overWrite );
        NcFile currFile = fileHandler.openFile();
        if ( !fileHandler.isFileOpen() ) {
            std::cout << " File " << outputFile << " didn't open!" << "\n";
            return SAVE_FAILURE;
        } else {
//            std::cout << "\n Starting saving to netCDF (file name: " << fileHandler.getFileName() <<  ") \n";
            time( &rawtime );
            strftime(buffer, sizeof(buffer),"%d-%m-%Y %H:%M:%S", localtime(&rawtime));

            const NcDim *timeDim = fileHandler.addDim( currFile, "Time", long(timeArray.size()) );
            didSaveSucceed *= fileHandler.addVar( currFile, &timeArray[0], "Time", timeDim, "float", "s", "Time");
            
            const NcDim *timeDim_midStep = fileHandler.addDim( currFile, "Time_mid", timeArray.size() - 1 );
            Vector_1D time_midStep( timeArray.size()-1, 0.0E+00 );

            for ( unsigned int iTime = 0; iTime < timeArray.size() - 1; iTime++ )
                time_midStep[iTime] = 0.5 * (timeArray[iTime] + timeArray[iTime+1]);

            didSaveSucceed *= fileHandler.addVar( currFile, &time_midStep[0], "Time_mid", timeDim_midStep, "float", "s", "Time at mid time-step");
            
            didSaveSucceed *= fileHandler.addAtt( currFile, "FileName", fileHandler.getFileName() );
            didSaveSucceed *= fileHandler.addAtt( currFile, "Author", "Thibaud M. Fritz (fritzt@mit.edu)" );
            didSaveSucceed *= fileHandler.addAtt( currFile, "Contact", "Thibaud M. Fritz (fritzt@mit.edu)" );
            didSaveSucceed *= fileHandler.addAtt( currFile, "Generation Date", buffer );
            didSaveSucceed *= fileHandler.addAtt( currFile, "Format", "NetCDF-4" );

            double value = 0.0E+00;
            value = input.temperature_K();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Temperature"    , 1, "float", "K"  , "Ambient Temperature" );
            value = input.pressure_Pa() / 100.0;
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Pressure"       , 1, "float", "hPa", "Ambient Pressure" );
            value = airDens;
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Air Density"    , 1, "float", "molecule / cm ^ 3", "Molecular density" );
            value = input.relHumidity_w();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "RHW"            , 1, "float", "-"  , "Ambient Rel. Humidity w.r.t water" );
            value = relHumidity_i;
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "RHI"            , 1, "float", "-"  , "Ambient Rel. Humidity w.r.t ice" );
            value = input.shear();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Shear"          , 1, "float", "1/s", "Ambient wind shear" );
            value = input.emissionDOY();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Emission Day"   , 1, "int"  , "-"  , "Emission day" );
            value = input.emissionTime();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Emission Time"  , 1, "float", "hr"  , "Emission time" );
            value = input.longitude_deg();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Longitude"      , 1, "float", "deg", "Longitude" );
            value = input.latitude_deg();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Latitude"       , 1, "float", "deg", "Latitude" );
            value = input.EI_NOx();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "NOx EI"         , 1, "float", "g/kg_fuel"  , "NOx Emission index" );
            value = input.EI_CO();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "CO EI"          , 1, "float", "g/kg_fuel"  , "CO Emission index" );
            value = input.EI_HC();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "HC EI"          , 1, "float", "g/kg_fuel"  , "HC Emission index" );
            value = input.EI_SO2();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "SO2 EI"         , 1, "float", "g/kg_fuel"  , "SO2 Emission index" );
            value = input.EI_Soot();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Soot EI"        , 1, "float", "g/kg_fuel"  , "Soot Emission index" );
            value = input.sootRad();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Soot Radius"    , 1, "float", "g/kg_fuel"  , "Soot radius" );
            value = input.fuelFlow();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Fuel flow"      , 1, "float", "kg/s"  , "Engine fuel flow" );
            value = input.backgNOx();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Background NOx" , 1, "float", "ppb"  , "Background NOx mixing ratio" );
            value = input.backgHNO3();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Background HNO3", 1, "float", "ppb"  , "Background HNO3 mixing ratio" );
            value = input.backgO3();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Background O3"  , 1, "float", "ppb"  , "Background O3 mixing ratio" );
            value = input.backgCO();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Background CO"  , 1, "float", "ppb"  , "Background CO mixing ratio" );
            value = input.backgCH4();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Background CH4" , 1, "float", "ppb"  , "Background CH4 mixing ratio" );
            value = input.backgSO2();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Background SO2" , 1, "float", "ppb"  , "Background SO2 mixing ratio" );

            didSaveSucceed *= fileHandler.addVar( currFile, &(ambientData.cosSZA)[0], "CSZA", timeDim_midStep, "float", "-", "Cosine of the solar zenith angle" );

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
               
#if ( SAVE_TO_DOUBLE )
                        spcArray = util::vect2double( plumeData             , NT, scalingFactor );
                        adjArray = util::vect2double( adjointData.Species[N], NT, scalingFactor );
                        ambArray = util::vect2double( ambientData.Species[N], NT, scalingFactor );
#else
                        spcArray = util::vect2float ( plumeData             , NT, scalingFactor );
                        adjArray = util::vect2float ( adjointData.Species[N], NT, scalingFactor );
                        ambArray = util::vect2float ( ambientData.Species[N], NT, scalingFactor );
#endif /* SAVE_TO_DOUBLE */


                        strncpy( charSpc, SPC_NAMES[N], sizeof(charSpc) );
                        strcat(  charSpc, "_Plume" );
                        strncpy( charName, SPC_NAMES[N], sizeof(charName) );
                        strcat(  charSpc, " plume-averaged mixing ratio" );

                        didSaveSucceed *= fileHandler.addVar( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName ); 

                        strncpy( charSpc, SPC_NAMES[N], sizeof(charSpc) );
                        strcat(  charSpc, "_Adjoint" );
                        strncpy( charName, SPC_NAMES[N], sizeof(charName) );
                        strcat(  charSpc, " optimized mixing ratio" );

                        didSaveSucceed *= fileHandler.addVar( currFile, &(adjArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName ); 

                        strncpy( charSpc, SPC_NAMES[N], sizeof(charSpc) );
                        strcat(  charSpc, "_Ambient" );
                        strncpy( charName, SPC_NAMES[N], sizeof(charName) );
                        strcat(  charSpc, " ambient mixing ratio" );

                        didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName ); 
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
#if ( SAVE_TO_DOUBLE )
            spcArray = util::vect2double( plumeData , NT, scalingFactor );
            adjArray = util::vect2double( adjointNOx, NT, scalingFactor );
            ambArray = util::vect2double( ambientNOx, NT, scalingFactor );
#else
            spcArray = util::vect2float ( plumeData , NT, scalingFactor );
            adjArray = util::vect2float ( adjointNOx, NT, scalingFactor );
            ambArray = util::vect2float ( ambientNOx, NT, scalingFactor );
#endif /* SAVE_TO_DOUBLE */


            strncpy( charSpc, "NOx_Plume", sizeof(charSpc) );
            strncpy( charName, "NOx plume-averaged mixing ratio", sizeof(charName) );
            didSaveSucceed *= fileHandler.addVar( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName ); 

            strncpy( charSpc, "NOx_Adjoint", sizeof(charSpc) );
            strncpy( charName, "NOx optimized mixing ratio", sizeof(charName) );
            didSaveSucceed *= fileHandler.addVar( currFile, &(adjArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName ); 

            strncpy( charSpc, "NOx_Ambient", sizeof(charSpc) );
            strncpy( charName, "NOx ambient mixing ratio", sizeof(charName) );
            didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName ); 

            delete[] spcArray; spcArray = NULL;
            delete[] ambArray; ambArray = NULL;
            delete[] adjArray; adjArray = NULL;

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
#endif /* SAVE_TO_DOUBLE */

            strncpy( charSpc, "NOy_Plume", sizeof(charSpc) );
            strncpy( charName, "NOy plume-averaged mixing ratio", sizeof(charName) );
            didSaveSucceed *= fileHandler.addVar( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName ); 

            strncpy( charSpc, "NOy_Ambient", sizeof(charSpc) );
            strncpy( charName, "NOy ambient mixing ratio", sizeof(charName) );
            didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName ); 

            delete[] spcArray; spcArray = NULL;
            delete[] ambArray; ambArray = NULL;
            delete[] adjArray; adjArray = NULL;

            if ( didSaveSucceed == NC_SUCCESS ) {
//                std::cout << " Done saving to netCDF!" << "\n";
            } else if ( didSaveSucceed != NC_SUCCESS ) {
                std::cout << "Error occured in saving data: didSaveSucceed: " << didSaveSucceed << "\n";
                return SAVE_FAILURE;
            }

            fileHandler.closeFile( currFile );
            if ( fileHandler.isFileOpen() ) {
                std::cout << "File " << outputFile << " didn't close properly!" << "\n";
                return SAVE_FAILURE;
            }

        }

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
        const bool doWrite = 1;
        const bool doRead = 1;
        const bool overWrite = 1;

        int didSaveSucceed = 1;
        time_t rawtime;
        char buffer[80];

        FileHandler fileHandler( outputFile, doWrite, doRead, overWrite );
        NcFile currFile = fileHandler.openFile();
        if ( !fileHandler.isFileOpen() ) {
            std::cout << " File " << outputFile << " didn't open!" << "\n";
            return SAVE_FAILURE;
        } else {
//            std::cout << "\n Starting saving to netCDF (file name: " << fileHandler.getFileName() <<  ") \n";
            time( &rawtime );
            strftime(buffer, sizeof(buffer),"%d-%m-%Y %H:%M:%S", localtime(&rawtime));

            const NcDim *timeDim = fileHandler.addDim( currFile, "Time", long(timeArray.size()) );
            didSaveSucceed *= fileHandler.addVar( currFile, &timeArray[0], "Time", timeDim, "float", "s", "Time");
            
            const NcDim *timeDim_midStep = fileHandler.addDim( currFile, "Time_mid", timeArray.size() - 1 );
            Vector_1D time_midStep( timeArray.size()-1, 0.0E+00 );

            for ( unsigned int iTime = 0; iTime < timeArray.size() - 1; iTime++ )
                time_midStep[iTime] = 0.5 * (timeArray[iTime] + timeArray[iTime+1]);

            didSaveSucceed *= fileHandler.addVar( currFile, &time_midStep[0], "Time_mid", timeDim_midStep, "float", "s", "Time at mid time-step");
            
            didSaveSucceed *= fileHandler.addAtt( currFile, "FileName", fileHandler.getFileName() );
            didSaveSucceed *= fileHandler.addAtt( currFile, "Author", "Thibaud M. Fritz (fritzt@mit.edu)" );
            didSaveSucceed *= fileHandler.addAtt( currFile, "Contact", "Thibaud M. Fritz (fritzt@mit.edu)" );
            didSaveSucceed *= fileHandler.addAtt( currFile, "Generation Date", buffer );
            didSaveSucceed *= fileHandler.addAtt( currFile, "Format", "NetCDF-4" );

            double value = 0.0E+00;
            value = input.temperature_K();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Temperature"    , 1, "float", "K"  , "Ambient Temperature" );
            value = input.pressure_Pa() / 100.0;
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Pressure"       , 1, "float", "hPa", "Ambient Pressure" );
            value = airDens;
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Air Density"    , 1, "float", "molecule / cm ^ 3", "Molecular density" );
            value = input.relHumidity_w();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "RHW"            , 1, "float", "-"  , "Ambient Rel. Humidity w.r.t water" );
            value = relHumidity_i;
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "RHI"            , 1, "float", "-"  , "Ambient Rel. Humidity w.r.t ice" );
            value = input.shear();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Shear"          , 1, "float", "1/s", "Ambient wind shear" );
            value = input.emissionDOY();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Emission Day"   , 1, "int"  , "-"  , "Emission day" );
            value = input.emissionTime();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Emission Time"  , 1, "float", "hr"  , "Emission time" );
            value = input.longitude_deg();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Longitude"      , 1, "float", "deg", "Longitude" );
            value = input.latitude_deg();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Latitude"       , 1, "float", "deg", "Latitude" );
            value = input.EI_NOx();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "NOx EI"         , 1, "float", "g/kg_fuel"  , "NOx Emission index" );
            value = input.EI_CO();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "CO EI"          , 1, "float", "g/kg_fuel"  , "CO Emission index" );
            value = input.EI_HC();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "HC EI"          , 1, "float", "g/kg_fuel"  , "HC Emission index" );
            value = input.EI_SO2();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "SO2 EI"         , 1, "float", "g/kg_fuel"  , "SO2 Emission index" );
            value = input.EI_Soot();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Soot EI"        , 1, "float", "g/kg_fuel"  , "Soot Emission index" );
            value = input.sootRad();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Soot Radius"    , 1, "float", "g/kg_fuel"  , "Soot radius" );
            value = input.fuelFlow();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Fuel flow"      , 1, "float", "kg/s"  , "Engine fuel flow" );
            value = input.backgNOx();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Background NOx" , 1, "float", "ppb"  , "Background NOx mixing ratio" );
            value = input.backgHNO3();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Background HNO3", 1, "float", "ppb"  , "Background HNO3 mixing ratio" );
            value = input.backgO3();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Background O3"  , 1, "float", "ppb"  , "Background O3 mixing ratio" );
            value = input.backgCO();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Background CO"  , 1, "float", "ppb"  , "Background CO mixing ratio" );
            value = input.backgCH4();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Background CH4" , 1, "float", "ppb"  , "Background CH4 mixing ratio" );
            value = input.backgSO2();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Background SO2" , 1, "float", "ppb"  , "Background SO2 mixing ratio" );

            didSaveSucceed *= fileHandler.addVar( currFile, &(boxData.cosSZA)[0], "CSZA", timeDim_midStep, "float", "-", "Cosine of the solar zenith angle" );

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
#else
                        ambArray = util::vect2float ( boxData.Species[N] , NT, scalingFactor );
#endif /* SAVE_TO_DOUBLE */

                        strncpy( charSpc, SPC_NAMES[N], sizeof(charSpc) );
                        strcat(  charSpc, "_Ambient" );
                        strncpy( charName, SPC_NAMES[N], sizeof(charName) );
                        strcat(  charName, " ambient mixing ratio" );
                        
                        didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName ); 
                    }
                }
            }

            scalingFactor = TO_PPB;
            strncpy( charUnit, "ppb", sizeof(charUnit) );
            
            std::vector<double> ambientNOx = util::add1D( boxData.Species[ind_NO], boxData.Species[ind_NO2] );
#if ( SAVE_TO_DOUBLE )
                ambArray = util::vect2double( ambientNOx, NT, scalingFactor );
#else
                ambArray = util::vect2float ( ambientNOx, NT, scalingFactor );
#endif /* SAVE_TO_DOUBLE */


            strncpy( charSpc, "NOx_Ambient", sizeof(charSpc) );
            strncpy( charName, "NOx ambient mixing ratio", sizeof(charName) );
            didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName ); 

            delete[] ambArray; ambArray = NULL;

            if ( didSaveSucceed == NC_SUCCESS ) {
//                std::cout << " Done saving to netCDF!" << "\n";
            } else if ( didSaveSucceed != NC_SUCCESS ) {
                std::cout << "Error occured in saving data: didSaveSucceed: " << didSaveSucceed << "\n";
                return SAVE_FAILURE;
            }

            fileHandler.closeFile( currFile );
            if ( fileHandler.isFileOpen() ) {
                std::cout << "File " << outputFile << " didn't close properly!" << "\n";
                return SAVE_FAILURE;
            }

        }

        return SAVE_SUCCESS;

    } /* End of Write_Box */


}

/* End of Save.cpp */
