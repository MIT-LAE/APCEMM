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
               const SpeciesArray &ringSpecies, const Ambient ambientData,       \
               const Cluster &ringCluster, const std::vector<double> &timeArray, \
               const Input &input,                                               \
               const double &airDens, const double &relHumidity_i,               \
               const double &sunRise, const double &sunSet ) 
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

            const NcDim *timeDim = fileHandler.addDim( currFile, "Time", long(timeArray.size()) );
            didSaveSucceed *= fileHandler.addVar( currFile, &timeArray[0], "Time", timeDim, "float", "s", "Time");
           
#if ( RINGS )

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
            value = input.pressure_Pa();
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

            didSaveSucceed *= fileHandler.addVar( currFile, &(ambientData.cosSZA)[0], "CSZA", timeDim, "float", "-", "Cosine of the solar zenith angle" );

#if ( RINGS )

            didSaveSucceed *= fileHandler.addVar( currFile, &(ringCluster.getRingArea())[0], "Ring Area", ringDim, "float", "m^2", "Ring Area" );

#endif /* RINGS */

/* Define conversion factors */
#define TO_PPTH         1.0 / airDens * 1.0E+03 /* Conversion factor from molecule/cm^3 to PPTH */
#define TO_PPM          1.0 / airDens * 1.0E+06 /* Conversion factor from molecule/cm^3 to PPM  */
#define TO_PPB          1.0 / airDens * 1.0E+09 /* Conversion factor from molecule/cm^3 to PPB  */
#define TO_PPT          1.0 / airDens * 1.0E+12 /* Conversion factor from molecule/cm^3 to PPT  */

#if ( RINGS )

#if ( SAVE_TO_DOUBLE )
                double* spcArray;
                double* ambArray;
                const char* outputType = "double";
#else
                float* spcArray;
                float* ambArray;
                const char* outputType = "float";
#endif /* SAVE_TO_DOUBLE */

                /* Start saving species ... */

                /* Define conversion factor */
                double scalingFactor;
                char charSpc[30];
                char charName[50];
                char charUnit[20];

#if ( DO_SAVE_CO2 )
                strncpy( charSpc, "CO2", sizeof(charSpc) );
                strncpy( charName, "CO2 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPM;
                strncpy( charUnit, "ppm", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.CO2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.CO2, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.CO2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.CO2, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName ); 

                strncpy( charSpc, "CO2_a", sizeof(charSpc) );
                strncpy( charName, "CO2 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName ); 

#endif /* DO_SAVE_CO2 */
            
#if ( DO_SAVE_PPN )
                strncpy( charSpc, "PPN", sizeof(charSpc) );
                strncpy( charName, "PPN mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.PPN, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.PPN, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.PPN, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.PPN, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "PPN_a", sizeof(charSpc) );
                strncpy( charName, "PPN ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_PPN */

#if ( DO_SAVE_BrNO2 ) 
                strncpy( charSpc, "BrNO2", sizeof(charSpc) );
                strncpy( charName, "BrNO2 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.BrNO2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.BrNO2, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.BrNO2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.BrNO2, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "BrNO2_a", sizeof(charSpc) );
                strncpy( charName, "BrNO2 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_BrNO2 */

#if ( DO_SAVE_IEPOX )
                strncpy( charSpc, "IEPOX", sizeof(charSpc) );
                strncpy( charName, "IEPOX mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.IEPOX, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.IEPOX, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.IEPOX, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.IEPOX, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "IEPOX_a", sizeof(charSpc) );
                strncpy( charName, "IEPOX ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_IEPOX */

#if ( DO_SAVE_PMNN ) 
                strncpy( charSpc, "PMNN", sizeof(charSpc) );
                strncpy( charName, "PMNN mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.PMNN, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.PMNN, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.PMNN, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.PMNN, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "PMNN_a", sizeof(charSpc) );
                strncpy( charName, "PMNN ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_PMNN */

#if ( DO_SAVE_N2O )
                strncpy( charSpc, "N2O", sizeof(charSpc) );
                strncpy( charName, "N2O mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.N2O, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.N2O, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.N2O, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.N2O, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "N2O_a", sizeof(charSpc) );
                strncpy( charName, "N2O ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_N2O */

#if ( DO_SAVE_N )
                strncpy( charSpc, "N", sizeof(charSpc) );
                strncpy( charName, "N mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.N, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.N, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.N, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.N, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "N_a", sizeof(charSpc) );
                strncpy( charName, "N ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_N */

#if ( DO_SAVE_PAN )
                strncpy( charSpc, "PAN", sizeof(charSpc) );
                strncpy( charName, "PAN mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.PAN, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.PAN, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.PAN, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.PAN, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "PAN_a", sizeof(charSpc) );
                strncpy( charName, "PAN ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_PAN */

#if ( DO_SAVE_ALK4 )
                strncpy( charSpc, "ALK4", sizeof(charSpc) );
                strncpy( charName, "ALK4 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.ALK4, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.ALK4, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.ALK4, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.ALK4, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "ALK4_a", sizeof(charSpc) );
                strncpy( charName, "ALK4 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_ALK4 */

#if ( DO_SAVE_MAP )
                strncpy( charSpc, "MAP", sizeof(charSpc) );
                strncpy( charName, "MAP mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.MAP, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.MAP, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.MAP, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.MAP, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "MAP_a", sizeof(charSpc) );
                strncpy( charName, "MAP ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_MAP */

#if ( DO_SAVE_MPN )
                strncpy( charSpc, "MPN", sizeof(charSpc) );
                strncpy( charName, "MPN mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.MPN, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.MPN, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.MPN, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.MPN, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "MPN_a", sizeof(charSpc) );
                strncpy( charName, "MPN ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_MPN */

#if ( DO_SAVE_Cl2O2 )
                strncpy( charSpc, "Cl2O2", sizeof(charSpc) );
                strncpy( charName, "Cl2O2 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.Cl2O2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.Cl2O2, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.Cl2O2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.Cl2O2, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "Cl2O2_a", sizeof(charSpc) );
                strncpy( charName, "Cl2O2 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_Cl2O2 */

#if ( DO_SAVE_ETP )
                strncpy( charSpc, "ETP", sizeof(charSpc) );
                strncpy( charName, "ETP mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.ETP, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.ETP, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.ETP, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.ETP, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "ETP_a", sizeof(charSpc) );
                strncpy( charName, "ETP ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_ETP */

#if ( DO_SAVE_HNO2 )
                strncpy( charSpc, "HNO2", sizeof(charSpc) );
                strncpy( charName, "HNO2 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.HNO2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.HNO2, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.HNO2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.HNO2, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "HNO2_a", sizeof(charSpc) );
                strncpy( charName, "HNO2 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_HNO2 */

#if ( DO_SAVE_C3H8 )
                strncpy( charSpc, "C3H8", sizeof(charSpc) );
                strncpy( charName, "C3H8 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.C3H8, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.C3H8, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.C3H8, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.C3H8, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "C3H8_a", sizeof(charSpc) );
                strncpy( charName, "C3H8 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_C3H8 */

#if ( DO_SAVE_RA3P )
                strncpy( charSpc, "RA3P", sizeof(charSpc) );
                strncpy( charName, "RA3P mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.RA3P, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.RA3P, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.RA3P, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.RA3P, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "RA3P_a", sizeof(charSpc) );
                strncpy( charName, "RA3P ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_RA3P */

#if ( DO_SAVE_RB3P )
                strncpy( charSpc, "RB3P", sizeof(charSpc) );
                strncpy( charName, "RB3P mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.RB3P, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.RB3P, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.RB3P, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.RB3P, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "RB3P_a", sizeof(charSpc) );
                strncpy( charName, "RB3P ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_RB3P */

#if ( DO_SAVE_OClO )
                strncpy( charSpc, "OClO", sizeof(charSpc) );
                strncpy( charName, "OClO mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.OClO, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.OClO, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.OClO, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.OClO, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "OClO_a", sizeof(charSpc) );
                strncpy( charName, "OClO ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_OClO */

#if ( DO_SAVE_ClNO2 )
                strncpy( charSpc, "ClNO2", sizeof(charSpc) );
                strncpy( charName, "ClNO2 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.ClNO2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.ClNO2, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.ClNO2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.ClNO2, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "ClNO2_a", sizeof(charSpc) );
                strncpy( charName, "ClNO2 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_ClNO2 */

#if ( DO_SAVE_ISOP )
                strncpy( charSpc, "ISOP", sizeof(charSpc) );
                strncpy( charName, "ISOP mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.ISOP, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.ISOP, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.ISOP, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.ISOP, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "ISOP_a", sizeof(charSpc) );
                strncpy( charName, "ISOP ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_ISOP */

#if ( DO_SAVE_HNO4 )
                strncpy( charSpc, "HNO4", sizeof(charSpc) );
                strncpy( charName, "HNO4 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.HNO4, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.HNO4, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.HNO4, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.HNO4, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "HNO4_a", sizeof(charSpc) );
                strncpy( charName, "HNO4 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_HNO4 */

#if ( DO_SAVE_MAOP )
                strncpy( charSpc, "MAOP", sizeof(charSpc) );
                strncpy( charName, "MAOP mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.MAOP, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.MAOP, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.MAOP, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.MAOP, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "MAOP_a", sizeof(charSpc) );
                strncpy( charName, "MAOP ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_MAOP */

#if ( DO_SAVE_MP )
                strncpy( charSpc, "MP", sizeof(charSpc) );
                strncpy( charName, "MP mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.MP, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.MP, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.MP, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.MP, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "MP_a", sizeof(charSpc) );
                strncpy( charName, "MP ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_MP */

#if ( DO_SAVE_ClOO )
                strncpy( charSpc, "ClOO", sizeof(charSpc) );
                strncpy( charName, "ClOO mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.ClOO, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.ClOO, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.ClOO, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.ClOO, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "ClOO_a", sizeof(charSpc) );
                strncpy( charName, "ClOO ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_ClOO */

#if ( DO_SAVE_RP )
                strncpy( charSpc, "RP", sizeof(charSpc) );
                strncpy( charName, "RP mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.RP, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.RP, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.RP, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.RP, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "RP_a", sizeof(charSpc) );
                strncpy( charName, "RP ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_RP */

#if ( DO_SAVE_BrCl )
                strncpy( charSpc, "BrCl", sizeof(charSpc) );
                strncpy( charName, "BrCl mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.BrCl, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.BrCl, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.BrCl, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.BrCl, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "BrCl_a", sizeof(charSpc) );
                strncpy( charName, "BrCl ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_BrCl */

#if ( DO_SAVE_PP )
                strncpy( charSpc, "PP", sizeof(charSpc) );
                strncpy( charName, "PP mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.PP, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.PP, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.PP, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.PP, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "PP_a", sizeof(charSpc) );
                strncpy( charName, "PP ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_PP */

#if ( DO_SAVE_PRPN )
                strncpy( charSpc, "PRPN", sizeof(charSpc) );
                strncpy( charName, "PRPN mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.PRPN, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.PRPN, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.PRPN, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.PRPN, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "PRPN_a", sizeof(charSpc) );
                strncpy( charName, "PRPN ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_PRPN */

#if ( DO_SAVE_SO4 )
                strncpy( charSpc, "SO4", sizeof(charSpc) );
                strncpy( charName, "SO4 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.SO4, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.SO4, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.SO4, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.SO4, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "SO4_a", sizeof(charSpc) );
                strncpy( charName, "SO4 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_SO4 */

#if ( DO_SAVE_SO4_L )
                strncpy( charSpc, "SO4_L", sizeof(charSpc) );
                strncpy( charName, "SO4_L mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.SO4L, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.SO4L, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.SO4L, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.SO4L, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "SO4_L_a", sizeof(charSpc) );
                strncpy( charName, "SO4_L ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_SO4_L */

#if ( DO_SAVE_Br2 )
                strncpy( charSpc, "Br2", sizeof(charSpc) );
                strncpy( charName, "Br2 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.Br2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.Br2, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.Br2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.Br2, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "Br2_a", sizeof(charSpc) );
                strncpy( charName, "Br2 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_Br2 */

#if ( DO_SAVE_ETHLN )
                strncpy( charSpc, "ETHLN", sizeof(charSpc) );
                strncpy( charName, "ETHLN mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.ETHLN, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.ETHLN, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.ETHLN, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.ETHLN, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "ETHLN_a", sizeof(charSpc) );
                strncpy( charName, "ETHLN ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_ETHLN */

#if ( DO_SAVE_MVKN )
                strncpy( charSpc, "MVKN", sizeof(charSpc) );
                strncpy( charName, "MVKN mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.MVKN, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.MVKN, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.MVKN, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.MVKN, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "MVKN_a", sizeof(charSpc) );
                strncpy( charName, "MVKN ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_MVKN */

#if ( DO_SAVE_R4P )
                strncpy( charSpc, "R4P", sizeof(charSpc) );
                strncpy( charName, "R4P mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.R4P, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.R4P, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.R4P, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.R4P, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "R4P_a", sizeof(charSpc) );
                strncpy( charName, "R4P ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_R4P */

#if ( DO_SAVE_C2H6 )
                strncpy( charSpc, "C2H6", sizeof(charSpc) );
                strncpy( charName, "C2H6 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.C2H6, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.C2H6, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.C2H6, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.C2H6, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "C2H6_a", sizeof(charSpc) );
                strncpy( charName, "C2H6 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_C2H6 */

#if ( DO_SAVE_RIP )
                strncpy( charSpc, "RIP", sizeof(charSpc) );
                strncpy( charName, "RIP mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.RIP, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.RIP, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.RIP, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.RIP, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "RIP_a", sizeof(charSpc) );
                strncpy( charName, "RIP ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_RIP */

#if ( DO_SAVE_VRP )
                strncpy( charSpc, "VRP", sizeof(charSpc) );
                strncpy( charName, "VRP mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.VRP, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.VRP, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.VRP, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.VRP, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "VRP_a", sizeof(charSpc) );
                strncpy( charName, "VRP ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_VRP */

#if ( DO_SAVE_ATOOH )
                strncpy( charSpc, "ATOOH", sizeof(charSpc) );
                strncpy( charName, "ATOOH mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.ATOOH, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.ATOOH, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.ATOOH, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.ATOOH, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "ATOOH_a", sizeof(charSpc) );
                strncpy( charName, "ATOOH ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_ATOOH */

#if ( DO_SAVE_IAP )
                strncpy( charSpc, "IAP", sizeof(charSpc) );
                strncpy( charName, "IAP mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.IAP, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.IAP, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.IAP, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.IAP, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "IAP_a", sizeof(charSpc) );
                strncpy( charName, "IAP ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_IAP */

#if ( DO_SAVE_DHMOB )
                strncpy( charSpc, "DHMOB", sizeof(charSpc) );
                strncpy( charName, "DHMOB mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.DHMOB, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.DHMOB, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.DHMOB, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.DHMOB, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "DHMOB_a", sizeof(charSpc) );
                strncpy( charName, "DHMOB ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_DHMOB */

#if ( DO_SAVE_MOBA )
                strncpy( charSpc, "MOBA", sizeof(charSpc) );
                strncpy( charName, "MOBA mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.MOBA, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.MOBA, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.MOBA, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.MOBA, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "MOBA_a", sizeof(charSpc) );
                strncpy( charName, "MOBA ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_MOBA */

#if ( DO_SAVE_MRP )
                strncpy( charSpc, "MRP", sizeof(charSpc) );
                strncpy( charName, "MRP mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.MRP, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.MRP, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.MRP, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.MRP, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "MRP_a", sizeof(charSpc) );
                strncpy( charName, "MRP ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_MRP */

#if ( DO_SAVE_N2O5 )
                strncpy( charSpc, "N2O5", sizeof(charSpc) );
                strncpy( charName, "N2O5 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.N2O5, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.N2O5, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.N2O5, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.N2O5, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "N2O5_a", sizeof(charSpc) );
                strncpy( charName, "N2O5 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_N2O5 */

#if ( DO_SAVE_ISNOHOO )
                strncpy( charSpc, "ISNOHOO", sizeof(charSpc) );
                strncpy( charName, "ISNOHOO mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.ISNOHOO, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.ISNOHOO, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.ISNOHOO, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.ISNOHOO, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "ISNOHOO_a", sizeof(charSpc) );
                strncpy( charName, "ISNOHOO ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_ISNOHOO */

#if ( DO_SAVE_ISNP )
                strncpy( charSpc, "ISNP", sizeof(charSpc) );
                strncpy( charName, "ISNP mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.ISNP, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.ISNP, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.ISNP, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.ISNP, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "ISNP_a", sizeof(charSpc) );
                strncpy( charName, "ISNP ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_ISNP */

#if ( DO_SAVE_ISOPNB )
                strncpy( charSpc, "ISOPNB", sizeof(charSpc) );
                strncpy( charName, "ISOPNB mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.ISOPNB, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.ISOPNB, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.ISOPNB, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.ISOPNB, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "ISOPNB_a", sizeof(charSpc) );
                strncpy( charName, "ISOPNB ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_ISOPNB */

#if ( DO_SAVE_IEPOXOO )
                strncpy( charSpc, "IEPOXOO", sizeof(charSpc) );
                strncpy( charName, "IEPOXOO mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.IEPOXOO, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.IEPOXOO, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.IEPOXOO, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.IEPOXOO, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "IEPOXOO_a", sizeof(charSpc) );
                strncpy( charName, "IEPOXOO ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_IEPOXOO */

#if ( DO_SAVE_MACRNO2 )
                strncpy( charSpc, "MACRNO2", sizeof(charSpc) );
                strncpy( charName, "MACRNO2 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.MACRNO2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.MACRNO2, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.MACRNO2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.MACRNO2, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "MACRNO2_a", sizeof(charSpc) );
                strncpy( charName, "MACRNO2 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_MACRNO2 */

#if ( DO_SAVE_ROH )
                strncpy( charSpc, "ROH", sizeof(charSpc) );
                strncpy( charName, "ROH mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.ROH, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.ROH, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.ROH, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.ROH, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "ROH_a", sizeof(charSpc) );
                strncpy( charName, "ROH ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_ROH */

#if ( DO_SAVE_MOBAOO )
                strncpy( charSpc, "MOBAOO", sizeof(charSpc) );
                strncpy( charName, "MOBAOO mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.MOBAOO, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.MOBAOO, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.MOBAOO, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.MOBAOO, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "MOBAOO_a", sizeof(charSpc) );
                strncpy( charName, "MOBAOO ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_MOBAOO */

#if ( DO_SAVE_DIBOO )
                strncpy( charSpc, "DIBOO", sizeof(charSpc) );
                strncpy( charName, "DIBOO mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.DIBOO, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.DIBOO, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.DIBOO, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.DIBOO, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "DIBOO_a", sizeof(charSpc) );
                strncpy( charName, "DIBOO ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_DIBOO */

#if ( DO_SAVE_PMN )
                strncpy( charSpc, "PMN", sizeof(charSpc) );
                strncpy( charName, "PMN mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.PMN, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.PMN, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.PMN, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.PMN, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "PMN_a", sizeof(charSpc) );
                strncpy( charName, "PMN ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_PMN */

#if ( DO_SAVE_ISNOOB )
                strncpy( charSpc, "ISNOOB", sizeof(charSpc) );
                strncpy( charName, "ISNOOB mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.ISNOOB, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.ISNOOB, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.ISNOOB, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.ISNOOB, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "ISNOOB_a", sizeof(charSpc) );
                strncpy( charName, "ISNOOB ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_ISNOOB */

#if ( DO_SAVE_INPN )
                strncpy( charSpc, "INPN", sizeof(charSpc) );
                strncpy( charName, "INPN mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.INPN, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.INPN, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.INPN, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.INPN, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "INPN_a", sizeof(charSpc) );
                strncpy( charName, "INPN ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_INPN */

#if ( DO_SAVE_H )
                strncpy( charSpc, "H", sizeof(charSpc) );
                strncpy( charName, "H mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.H, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.H, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.H, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.H, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "H_a", sizeof(charSpc) );
                strncpy( charName, "H ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_H */

#if ( DO_SAVE_BrNO3 )
                strncpy( charSpc, "BrNO3", sizeof(charSpc) );
                strncpy( charName, "BrNO3 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.BrNO3, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.BrNO3, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.BrNO3, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.BrNO3, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "BrNO3_a", sizeof(charSpc) );
                strncpy( charName, "BrNO3 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_BrNO3 */

#if ( DO_SAVE_PRPE )
                strncpy( charSpc, "PRPE", sizeof(charSpc) );
                strncpy( charName, "PRPE mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.PRPE, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.PRPE, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.PRPE, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.PRPE, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "PRPE_a", sizeof(charSpc) );
                strncpy( charName, "PRPE ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_PRPE */

#if ( DO_SAVE_MVKOO )
                strncpy( charSpc, "MVKOO", sizeof(charSpc) );
                strncpy( charName, "MVKOO mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.MVKOO, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.MVKOO, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.MVKOO, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.MVKOO, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "MVKOO_a", sizeof(charSpc) );
                strncpy( charName, "MVKOO ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_MVKOO */

#if ( DO_SAVE_Cl2 )
                strncpy( charSpc, "Cl2", sizeof(charSpc) );
                strncpy( charName, "Cl2 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.Cl2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.Cl2, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.Cl2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.Cl2, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "Cl2_a", sizeof(charSpc) );
                strncpy( charName, "Cl2 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_Cl2 */

#if ( DO_SAVE_ISOPND )
                strncpy( charSpc, "ISOPND", sizeof(charSpc) );
                strncpy( charName, "ISOPND mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.ISOPND, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.ISOPND, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.ISOPND, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.ISOPND, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "ISOPND_a", sizeof(charSpc) );
                strncpy( charName, "ISOPND ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_ISOPND */

#if ( DO_SAVE_HOBr )
                strncpy( charSpc, "HOBr", sizeof(charSpc) );
                strncpy( charName, "HOBr mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.HOBr, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.HOBr, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.HOBr, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.HOBr, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "HOBr_a", sizeof(charSpc) );
                strncpy( charName, "HOBr ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_HOBr */

#if ( DO_SAVE_HOBr_L )
                strncpy( charSpc, "HOBr_L", sizeof(charSpc) );
                strncpy( charName, "HOBr_L mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.HOBrL, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.HOBrL, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.HOBrL, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.HOBrL, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "HOBr_L_a", sizeof(charSpc) );
                strncpy( charName, "HOBr_L ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_HOBr_L */

#if ( DO_SAVE_A3O2 )
                strncpy( charSpc, "A3O2", sizeof(charSpc) );
                strncpy( charName, "A3O2 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.A3O2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.A3O2, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.A3O2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.A3O2, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "A3O2_a", sizeof(charSpc) );
                strncpy( charName, "A3O2 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_A3O2 */

#if ( DO_SAVE_PROPNN )
                strncpy( charSpc, "PROPNN", sizeof(charSpc) );
                strncpy( charName, "PROPNN mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.PROPNN, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.PROPNN, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.PROPNN, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.PROPNN, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "PROPNN_a", sizeof(charSpc) );
                strncpy( charName, "PROPNN ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_PROPNN */

#if ( DO_SAVE_GLYX )
                strncpy( charSpc, "GLYX", sizeof(charSpc) );
                strncpy( charName, "GLYX mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.GLYX, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.GLYX, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.GLYX, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.GLYX, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "GLYX_a", sizeof(charSpc) );
                strncpy( charName, "GLYX ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_GLYX */

#if ( DO_SAVE_MAOPO2 )
                strncpy( charSpc, "MAOPO2", sizeof(charSpc) );
                strncpy( charName, "MAOPO2 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.MAOPO2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.MAOPO2, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.MAOPO2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.MAOPO2, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "MAOPO2_a", sizeof(charSpc) );
                strncpy( charName, "MAOPO2 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_MAOPO2 */

#if ( DO_SAVE_CH4 )
                strncpy( charSpc, "CH4", sizeof(charSpc) );
                strncpy( charName, "CH4 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.CH4, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.CH4, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.CH4, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.CH4, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "CH4_a", sizeof(charSpc) );
                strncpy( charName, "CH4 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_CH4 */

#if ( DO_SAVE_GAOO )
                strncpy( charSpc, "GAOO", sizeof(charSpc) );
                strncpy( charName, "GAOO mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.GAOO, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.GAOO, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.GAOO, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.GAOO, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "GAOO_a", sizeof(charSpc) );
                strncpy( charName, "GAOO ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_GAOO */

#if ( DO_SAVE_B3O2 )
                strncpy( charSpc, "B3O2", sizeof(charSpc) );
                strncpy( charName, "B3O2 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.B3O2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.B3O2, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.B3O2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.B3O2, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "B3O2_a", sizeof(charSpc) );
                strncpy( charName, "B3O2 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_B3O2 */

#if ( DO_SAVE_ACET )
                strncpy( charSpc, "ACET", sizeof(charSpc) );
                strncpy( charName, "ACET mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.ACET, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.ACET, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.ACET, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.ACET, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "ACET_a", sizeof(charSpc) );
                strncpy( charName, "ACET ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_ACET */

#if ( DO_SAVE_MACRN )
                strncpy( charSpc, "MACRN", sizeof(charSpc) );
                strncpy( charName, "MACRN mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.MACRN, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.MACRN, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.MACRN, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.MACRN, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "MACRN_a", sizeof(charSpc) );
                strncpy( charName, "MACRN ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_MACRN */

#if ( DO_SAVE_CH2OO )
                strncpy( charSpc, "CH2OO", sizeof(charSpc) );
                strncpy( charName, "CH2OO mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.CH2OO, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.CH2OO, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.CH2OO, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.CH2OO, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "CH2OO_a", sizeof(charSpc) );
                strncpy( charName, "CH2OO ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_CH2OO */

#if ( DO_SAVE_MGLYOO )
                strncpy( charSpc, "MGLYOO", sizeof(charSpc) );
                strncpy( charName, "MGLYOO mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.MGLYOO, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.MGLYOO, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.MGLYOO, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.MGLYOO, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "MGLYOO_a", sizeof(charSpc) );
                strncpy( charName, "MGLYOO ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_MGLYOO */

#if ( DO_SAVE_VRO2 )
                strncpy( charSpc, "VRO2", sizeof(charSpc) );
                strncpy( charName, "VRO2 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.VRO2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.VRO2, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.VRO2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.VRO2, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "VRO2_a", sizeof(charSpc) );
                strncpy( charName, "VRO2 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_VRO2 */

#if ( DO_SAVE_MGLOO )
                strncpy( charSpc, "MGLOO", sizeof(charSpc) );
                strncpy( charName, "MGLOO mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.MGLOO, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.MGLOO, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.MGLOO, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.MGLOO, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "MGLOO_a", sizeof(charSpc) );
                strncpy( charName, "MGLOO ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_MGLOO */

#if ( DO_SAVE_MACROO )
                strncpy( charSpc, "MACROO", sizeof(charSpc) );
                strncpy( charName, "MACROO mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.MACROO, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.MACROO, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.MACROO, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.MACROO, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "MACROO_a", sizeof(charSpc) );
                strncpy( charName, "MACROO ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_MACROO */

#if ( DO_SAVE_PO2 )
                strncpy( charSpc, "PO2", sizeof(charSpc) );
                strncpy( charName, "PO2 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.PO2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.PO2, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.PO2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.PO2, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "PO2_a", sizeof(charSpc) );
                strncpy( charName, "PO2 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_PO2 */

#if ( DO_SAVE_CH3CHOO )
                strncpy( charSpc, "CH3CHOO", sizeof(charSpc) );
                strncpy( charName, "CH3CHOO mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.CH3CHOO, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.CH3CHOO, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.CH3CHOO, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.CH3CHOO, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "CH3CHOO_a", sizeof(charSpc) );
                strncpy( charName, "CH3CHOO ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_CH3CHOO */

#if ( DO_SAVE_MAN2 )
                strncpy( charSpc, "MAN2", sizeof(charSpc) );
                strncpy( charName, "MAN2 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.MAN2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.MAN2, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.MAN2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.MAN2, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "MAN2_a", sizeof(charSpc) );
                strncpy( charName, "MAN2 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_MAN2 */

#if ( DO_SAVE_ISNOOA )
                strncpy( charSpc, "ISNOOA", sizeof(charSpc) );
                strncpy( charName, "ISNOOA mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.ISNOOA, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.ISNOOA, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.ISNOOA, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.ISNOOA, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "ISNOOA_a", sizeof(charSpc) );
                strncpy( charName, "ISNOOA ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_ISNOOA */

#if ( DO_SAVE_H2O2 )
                strncpy( charSpc, "H2O2", sizeof(charSpc) );
                strncpy( charName, "H2O2 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.H2O2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.H2O2, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.H2O2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.H2O2, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "H2O2_a", sizeof(charSpc) );
                strncpy( charName, "H2O2 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_H2O2 */

#if ( DO_SAVE_PRN1 )
                strncpy( charSpc, "PRN1", sizeof(charSpc) );
                strncpy( charName, "PRN1 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.PRN1, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.PRN1, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.PRN1, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.PRN1, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "PRN1_a", sizeof(charSpc) );
                strncpy( charName, "PRN1 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_PRN1 */

#if ( DO_SAVE_ETO2 )
                strncpy( charSpc, "ETO2", sizeof(charSpc) );
                strncpy( charName, "ETO2 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.ETO2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.ETO2, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.ETO2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.ETO2, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "ETO2_a", sizeof(charSpc) );
                strncpy( charName, "ETO2 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_ETO2 */

#if ( DO_SAVE_KO2 )
                strncpy( charSpc, "KO2", sizeof(charSpc) );
                strncpy( charName, "KO2 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.KO2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.KO2, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.KO2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.KO2, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "KO2_a", sizeof(charSpc) );
                strncpy( charName, "KO2 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_KO2 */

#if ( DO_SAVE_RCO3 )
                strncpy( charSpc, "RCO3", sizeof(charSpc) );
                strncpy( charName, "RCO3 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.RCO3, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.RCO3, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.RCO3, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.RCO3, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "RCO3_a", sizeof(charSpc) );
                strncpy( charName, "RCO3 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_RCO3 */

#if ( DO_SAVE_HC5OO )
                strncpy( charSpc, "HC5OO", sizeof(charSpc) );
                strncpy( charName, "HC5OO mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.HC5OO, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.HC5OO, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.HC5OO, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.HC5OO, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "HC5OO_a", sizeof(charSpc) );
                strncpy( charName, "HC5OO ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_HC5OO */

#if ( DO_SAVE_GLYC )
                strncpy( charSpc, "GLYC", sizeof(charSpc) );
                strncpy( charName, "GLYC mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.GLYC, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.GLYC, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.GLYC, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.GLYC, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "GLYC_a", sizeof(charSpc) );
                strncpy( charName, "GLYC ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_GLYC */

#if ( DO_SAVE_ClNO3 )
                strncpy( charSpc, "ClNO3", sizeof(charSpc) );
                strncpy( charName, "ClNO3 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.ClNO3, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.ClNO3, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.ClNO3, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.ClNO3, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "ClNO3_a", sizeof(charSpc) );
                strncpy( charName, "ClNO3 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_ClNO3 */

#if ( DO_SAVE_RIO2 )
                strncpy( charSpc, "RIO2", sizeof(charSpc) );
                strncpy( charName, "RIO2 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.RIO2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.RIO2, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.RIO2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.RIO2, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "RIO2_a", sizeof(charSpc) );
                strncpy( charName, "RIO2 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_RIO2 */

#if ( DO_SAVE_R4N1 )
                strncpy( charSpc, "R4N1", sizeof(charSpc) );
                strncpy( charName, "R4N1 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.R4N1, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.R4N1, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.R4N1, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.R4N1, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "R4N1_a", sizeof(charSpc) );
                strncpy( charName, "R4N1 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_R4N1 */

#if ( DO_SAVE_HOCl )
                strncpy( charSpc, "HOCl", sizeof(charSpc) );
                strncpy( charName, "HOCl mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.HOCl, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.HOCl, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.HOCl, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.HOCl, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "HOCl_a", sizeof(charSpc) );
                strncpy( charName, "HOCl ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_HOCl */

#if ( DO_SAVE_HOCl_L )
                strncpy( charSpc, "HOCl_L", sizeof(charSpc) );
                strncpy( charName, "HOCl_L mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.HOClL, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.HOClL, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.HOClL, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.HOClL, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "HOCl_L_a", sizeof(charSpc) );
                strncpy( charName, "HOCl_L ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_HOCl_L */

#if ( DO_SAVE_ATO2 )
                strncpy( charSpc, "ATO2", sizeof(charSpc) );
                strncpy( charName, "ATO2 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.ATO2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.ATO2, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.ATO2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.ATO2, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "ATO2_a", sizeof(charSpc) );
                strncpy( charName, "ATO2 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_ATO2 */

#if ( DO_SAVE_HNO3 )
                strncpy( charSpc, "HNO3", sizeof(charSpc) );
                strncpy( charName, "HNO3 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.HNO3, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.HNO3, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.HNO3, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.HNO3, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "HNO3_a", sizeof(charSpc) );
                strncpy( charName, "HNO3 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_HNO3 */

#if ( DO_SAVE_HNO3_L )
                strncpy( charSpc, "HNO3_L", sizeof(charSpc) );
                strncpy( charName, "HNO3_L mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.HNO3L, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.HNO3L, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.HNO3L, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.HNO3L, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "HNO3_L_a", sizeof(charSpc) );
                strncpy( charName, "HNO3_L ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_HNO3_L */

#if ( DO_SAVE_HNO3_S )
                strncpy( charSpc, "HNO3_S", sizeof(charSpc) );
                strncpy( charName, "HNO3_S mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.HNO3S, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.HNO3S, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.HNO3S, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.HNO3S, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "HNO3_S_a", sizeof(charSpc) );
                strncpy( charName, "HNO3_S ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_HNO3_S */

#if ( DO_SAVE_ISN1 )
                strncpy( charSpc, "ISN1", sizeof(charSpc) );
                strncpy( charName, "ISN1 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.ISN1, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.ISN1, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.ISN1, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.ISN1, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "ISN1_a", sizeof(charSpc) );
                strncpy( charName, "ISN1 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_ISN1 */

#if ( DO_SAVE_MAO3 )
                strncpy( charSpc, "MAO3", sizeof(charSpc) );
                strncpy( charName, "MAO3 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.MAO3, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.MAO3, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.MAO3, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.MAO3, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "MAO3_a", sizeof(charSpc) );
                strncpy( charName, "MAO3 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_MAO3 */

#if ( DO_SAVE_MRO2 )
                strncpy( charSpc, "MRO2", sizeof(charSpc) );
                strncpy( charName, "MRO2 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.MRO2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.MRO2, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.MRO2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.MRO2, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "MRO2_a", sizeof(charSpc) );
                strncpy( charName, "MRO2 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_MRO2 */

#if ( DO_SAVE_INO2 )
                strncpy( charSpc, "INO2", sizeof(charSpc) );
                strncpy( charName, "INO2 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.INO2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.INO2, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.INO2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.INO2, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "INO2_a", sizeof(charSpc) );
                strncpy( charName, "INO2 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_INO2 */

#if ( DO_SAVE_HAC )
                strncpy( charSpc, "HAC", sizeof(charSpc) );
                strncpy( charName, "HAC mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.HAC, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.HAC, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.HAC, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.HAC, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "HAC_a", sizeof(charSpc) );
                strncpy( charName, "HAC ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_HAC */

#if ( DO_SAVE_HC5 )
                strncpy( charSpc, "HC5", sizeof(charSpc) );
                strncpy( charName, "HC5 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.HC5, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.HC5, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.HC5, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.HC5, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "HC5_a", sizeof(charSpc) );
                strncpy( charName, "HC5 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_HC5 */

#if ( DO_SAVE_MGLY )
                strncpy( charSpc, "MGLY", sizeof(charSpc) );
                strncpy( charName, "MGLY mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.MGLY, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.MGLY, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.MGLY, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.MGLY, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "MGLY_a", sizeof(charSpc) );
                strncpy( charName, "MGLY ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_MGLY */

#if ( DO_SAVE_ISOPNBO2 )
                strncpy( charSpc, "ISOPNBO2", sizeof(charSpc) );
                strncpy( charName, "ISOPNBO2 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.ISOPNBO2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.ISOPNBO2, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.ISOPNBO2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.ISOPNBO2, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "ISOPNBO2_a", sizeof(charSpc) );
                strncpy( charName, "ISOPNBO2 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_ISOPNBO2 */

#if ( DO_SAVE_ISOPNDO2 )
                strncpy( charSpc, "ISOPNDO2", sizeof(charSpc) );
                strncpy( charName, "ISOPNDO2 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.ISOPNDO2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.ISOPNDO2, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.ISOPNDO2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.ISOPNDO2, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "ISOPNDO2_a", sizeof(charSpc) );
                strncpy( charName, "ISOPNDO2 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_ISOPNDO2 */

#if ( DO_SAVE_R4O2 )
                strncpy( charSpc, "R4O2", sizeof(charSpc) );
                strncpy( charName, "R4O2 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.R4O2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.R4O2, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.R4O2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.R4O2, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "R4O2_a", sizeof(charSpc) );
                strncpy( charName, "R4O2 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_R4O2 */

#if ( DO_SAVE_R4N2 )
                strncpy( charSpc, "R4N2", sizeof(charSpc) );
                strncpy( charName, "R4N2 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.R4N2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.R4N2, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.R4N2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.R4N2, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "R4N2_a", sizeof(charSpc) );
                strncpy( charName, "R4N2 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_R4N2 */

#if ( DO_SAVE_BrO )
                strncpy( charSpc, "BrO", sizeof(charSpc) );
                strncpy( charName, "BrO mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.BrO, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.BrO, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.BrO, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.BrO, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "BrO_a", sizeof(charSpc) );
                strncpy( charName, "BrO ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_BrO */

#if ( DO_SAVE_RCHO )
                strncpy( charSpc, "RCHO", sizeof(charSpc) );
                strncpy( charName, "RCHO mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.RCHO, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.RCHO, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.RCHO, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.RCHO, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "RCHO_a", sizeof(charSpc) );
                strncpy( charName, "RCHO ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_RCHO */

#if ( DO_SAVE_MEK )
                strncpy( charSpc, "MEK", sizeof(charSpc) );
                strncpy( charName, "MEK mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.MEK, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.MEK, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.MEK, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.MEK, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "MEK_a", sizeof(charSpc) );
                strncpy( charName, "MEK ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_MEK */

#if ( DO_SAVE_ClO )
                strncpy( charSpc, "ClO", sizeof(charSpc) );
                strncpy( charName, "ClO mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.ClO, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.ClO, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.ClO, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.ClO, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "ClO_a", sizeof(charSpc) );
                strncpy( charName, "ClO ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_ClO */

#if ( DO_SAVE_MACR )
                strncpy( charSpc, "MACR", sizeof(charSpc) );
                strncpy( charName, "MACR mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.MACR, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.MACR, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.MACR, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.MACR, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "MACR_a", sizeof(charSpc) );
                strncpy( charName, "MACR ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_MACR */

#if ( DO_SAVE_SO2 )
                strncpy( charSpc, "SO2", sizeof(charSpc) );
                strncpy( charName, "SO2 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.SO2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.SO2, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.SO2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.SO2, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "SO2_a", sizeof(charSpc) );
                strncpy( charName, "SO2 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_SO2 */

#if ( DO_SAVE_MVK )
                strncpy( charSpc, "MVK", sizeof(charSpc) );
                strncpy( charName, "MVK mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.MVK, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.MVK, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.MVK, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.MVK, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "MVK_a", sizeof(charSpc) );
                strncpy( charName, "MVK ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_MVK */

#if ( DO_SAVE_ALD2 )
                strncpy( charSpc, "ALD2", sizeof(charSpc) );
                strncpy( charName, "ALD2 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.ALD2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.ALD2, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.ALD2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.ALD2, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "ALD2_a", sizeof(charSpc) );
                strncpy( charName, "ALD2 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_ALD2 */

#if ( DO_SAVE_MCO3 )
                strncpy( charSpc, "MCO3", sizeof(charSpc) );
                strncpy( charName, "MCO3 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.MCO3, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.MCO3, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.MCO3, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.MCO3, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "MCO3_a", sizeof(charSpc) );
                strncpy( charName, "MCO3 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_MCO3 */

#if ( DO_SAVE_CH2O )
                strncpy( charSpc, "CH2O", sizeof(charSpc) );
                strncpy( charName, "CH2O mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.CH2O, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.CH2O, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.CH2O, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.CH2O, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "CH2O_a", sizeof(charSpc) );
                strncpy( charName, "CH2O ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_CH2O */

#if ( DO_SAVE_H2O )
                strncpy( charSpc, "H2O", sizeof(charSpc) );
                strncpy( charName, "H2O mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.H2O, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.H2O, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.H2O, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.H2O, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "H2O_a", sizeof(charSpc) );
                strncpy( charName, "H2O ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_H2O */

#if ( DO_SAVE_H2O_L )
                strncpy( charSpc, "H2O_L", sizeof(charSpc) );
                strncpy( charName, "H2O_L mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.H2OL, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.H2OL, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.H2OL, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.H2OL, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "H2O_L_a", sizeof(charSpc) );
                strncpy( charName, "H2O_L ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_H2O_L */

#if ( DO_SAVE_H2O_S )
                strncpy( charSpc, "H2O_S", sizeof(charSpc) );
                strncpy( charName, "H2O_S mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.H2OS, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.H2OS, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.H2OS, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.H2OS, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "H2O_S_a", sizeof(charSpc) );
                strncpy( charName, "H2O_S ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_H2O_S */

#if ( DO_SAVE_Br )
                strncpy( charSpc, "Br", sizeof(charSpc) );
                strncpy( charName, "Br mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.Br, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.Br, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.Br, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.Br, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "Br_a", sizeof(charSpc) );
                strncpy( charName, "Br ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_Br */

#if ( DO_SAVE_NO )
                strncpy( charSpc, "NO", sizeof(charSpc) );
                strncpy( charName, "NO mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.NO, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.NO, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.NO, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.NO, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "NO_a", sizeof(charSpc) );
                strncpy( charName, "NO ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_NO */

#if ( DO_SAVE_NO3 )
                strncpy( charSpc, "NO3", sizeof(charSpc) );
                strncpy( charName, "NO3 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.NO3, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.NO3, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.NO3, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.NO3, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "NO3_a", sizeof(charSpc) );
                strncpy( charName, "NO3 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_NO3 */

#if ( DO_SAVE_Cl )
                strncpy( charSpc, "Cl", sizeof(charSpc) );
                strncpy( charName, "Cl mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.Cl, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.Cl, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.Cl, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.Cl, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "Cl_a", sizeof(charSpc) );
                strncpy( charName, "Cl ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_Cl */

#if ( DO_SAVE_O )
                strncpy( charSpc, "O", sizeof(charSpc) );
                strncpy( charName, "O mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.O, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.O, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.O, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.O, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "O_a", sizeof(charSpc) );
                strncpy( charName, "O ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_O */

#if ( DO_SAVE_O1D )
                strncpy( charSpc, "O1D", sizeof(charSpc) );
                strncpy( charName, "O1D mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.O1D, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.O1D, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.O1D, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.O1D, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "O1D_a", sizeof(charSpc) );
                strncpy( charName, "O1D ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_O1D */

#if ( DO_SAVE_O3 )
                strncpy( charSpc, "O3", sizeof(charSpc) );
                strncpy( charName, "O3 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.O3, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.O3, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.O3, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.O3, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "O3_a", sizeof(charSpc) );
                strncpy( charName, "O3 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_O3 */

#if ( DO_SAVE_HO2 )
                strncpy( charSpc, "HO2", sizeof(charSpc) );
                strncpy( charName, "HO2 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.HO2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.HO2, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.HO2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.HO2, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "HO2_a", sizeof(charSpc) );
                strncpy( charName, "HO2 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_HO2 */

#if ( DO_SAVE_NO2 )
                strncpy( charSpc, "NO2", sizeof(charSpc) );
                strncpy( charName, "NO2 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.NO2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.NO2, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.NO2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.NO2, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "NO2_a", sizeof(charSpc) );
                strncpy( charName, "NO2 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_NO2 */

#if ( DO_SAVE_OH )
                strncpy( charSpc, "OH", sizeof(charSpc) );
                strncpy( charName, "OH mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.OH, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.OH, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.OH, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.OH, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "OH_a", sizeof(charSpc) );
                strncpy( charName, "OH ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_OH */

#if ( DO_SAVE_HBr )
                strncpy( charSpc, "HBr", sizeof(charSpc) );
                strncpy( charName, "HBr mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.HBr, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.HBr, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.HBr, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.HBr, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "HBr_a", sizeof(charSpc) );
                strncpy( charName, "HBr ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_HBr */

#if ( DO_SAVE_HBr_L )
                strncpy( charSpc, "HBr_L", sizeof(charSpc) );
                strncpy( charName, "HBr_L mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.HBrL, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.HBrL, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.HBrL, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.HBrL, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "HBr_L_a", sizeof(charSpc) );
                strncpy( charName, "HBr_L ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_HBr_L */

#if ( DO_SAVE_HCl )
                strncpy( charSpc, "HCl", sizeof(charSpc) );
                strncpy( charName, "HCl mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.HCl, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.HCl, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.HCl, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.HCl, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "HCl_a", sizeof(charSpc) );
                strncpy( charName, "HCl ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_HCl */

#if ( DO_SAVE_HCl_L )
                strncpy( charSpc, "HCl_L", sizeof(charSpc) );
                strncpy( charName, "HCl_L mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.HClL, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.HClL, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.HClL, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.HClL, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "HCl_L_a", sizeof(charSpc) );
                strncpy( charName, "HCl_L ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_HCl_L */

#if ( DO_SAVE_CO )
                strncpy( charSpc, "CO", sizeof(charSpc) );
                strncpy( charName, "CO mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.CO, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.CO, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.CO, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.CO, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */


                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "CO_a", sizeof(charSpc) );
                strncpy( charName, "CO ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_CO */

#if ( DO_SAVE_MO2 )
                strncpy( charSpc, "MO2", sizeof(charSpc) );
                strncpy( charName, "MO2 mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( ringSpecies.MO2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( ambientData.MO2, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float ( ringSpecies.MO2, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float ( ambientData.MO2, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */

                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "MO2_a", sizeof(charSpc) );
                strncpy( charName, "MO2 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );

#endif /* DO_SAVE_MO2 */

#if ( DO_SAVE_NOx )

                strncpy( charSpc, "NOx", sizeof(charSpc) );
                strncpy( charName, "NOx mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

                std::vector<std::vector<double>> NOx = util::add2D( ringSpecies.NO, ringSpecies.NO2 );
                std::vector<double> NOx_a = util::add1D( ambientData.NO, ambientData.NO2 );
        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( NOx, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( NOx_a, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float( NOx, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float( NOx_a, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */
                
                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "NOx_a", sizeof(charSpc) );
                strncpy( charName, "NOx ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );
            
#endif /* DO_SAVE_NOx */

#if ( DO_SAVE_NOy )

                strncpy( charSpc, "NOy", sizeof(charSpc) );
                strncpy( charName, "NOy mixing ratio", sizeof(charName) );
                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

                std::vector<std::vector<double>> NOy = util::add2D( util::add2D( util::add2D( util::add2D( util::add2D( util::add2D( util::add2D( util::add2D( util::add2D( util::add2D( util::add2D( util::add2D( util::add2D( util::add2D( util::add2D( util::add2D(util::add2D( util::add2D( util::add2D( util::add2D( ringSpecies.NO, ringSpecies.NO2 ), ringSpecies.NO3 ), ringSpecies.HNO2 ), ringSpecies.HNO3 ), ringSpecies.HNO4 ), ringSpecies.N2O5 ), ringSpecies.N2O5 ), ringSpecies.PAN ), ringSpecies.BrNO2 ), ringSpecies.BrNO3 ), ringSpecies.ClNO2 ), ringSpecies.ClNO3 ), ringSpecies.PPN ), ringSpecies.N ), ringSpecies.MPN ), ringSpecies.PROPNN ), ringSpecies.PRPN ), ringSpecies.R4N1 ), ringSpecies.PRN1 ), ringSpecies.R4N2 );
                std::vector<double> NOy_a = util::add1D( util::add1D( util::add1D( util::add1D( util::add1D( util::add1D( util::add1D( util::add1D( util::add1D( util::add1D( util::add1D( util::add1D( util::add1D( util::add1D( util::add1D( util::add1D( util::add1D( util::add1D( util::add1D( util::add1D( ambientData.NO, ambientData.NO2 ), ambientData.NO3 ), ambientData.HNO2 ), ambientData.HNO3 ), ambientData.HNO4 ), ambientData.N2O5 ), ambientData.N2O5 ), ambientData.PAN ), ambientData. BrNO2 ), ambientData.BrNO3 ), ambientData.ClNO2 ), ambientData.ClNO3 ), ambientData.PPN ), ambientData.N ), ambientData.MPN ), ambientData.PROPNN ), ambientData.PRPN ), ambientData.R4N1 ), ambientData.PRN1 ), ambientData.R4N2 );
        #if ( SAVE_TO_DOUBLE )
                spcArray = util::vect2double( NOy, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2double( NOy_a, timeArray.size(), scalingFactor );
        #else
                spcArray = util::vect2float( NOy, timeArray.size(), ringCluster.getnRing(), scalingFactor );
                ambArray = util::vect2float( NOy_a, timeArray.size(), scalingFactor );
        #endif /* SAVE_TO_DOUBLE */
                
                didSaveSucceed *= fileHandler.addVar2D( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, ringDim, outputType, (const char*)charUnit, (const char*)charName );

                strncpy( charSpc, "NOy_a", sizeof(charSpc) );
                strncpy( charName, "NOy ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName );
            
#endif /* DO_SAVE_NOy */

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
            didSaveSucceed *= fileHandler.addConst( currFile, &pressure_Pa  , "Pressure"   , 1, "float", "hPa" , "Ambient Pressure" );
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

    int Write_Adjoint( const char* outputFile,                                     \
                       const SpeciesArray &ringSpecies, const Ambient ambientData, \
                       const Ambient adjointData,                                  \
                       const std::vector<double> &ringArea, const double totArea,  \
                       const std::vector<double> &timeArray,                       \
                       const Input &input,                                         \
                       const double &airDens, const double &relHumidity_i )
    {

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
            
            didSaveSucceed *= fileHandler.addAtt( currFile, "FileName", fileHandler.getFileName() );
            didSaveSucceed *= fileHandler.addAtt( currFile, "Author", "Thibaud M. Fritz (fritzt@mit.edu)" );
            didSaveSucceed *= fileHandler.addAtt( currFile, "Contact", "Thibaud M. Fritz (fritzt@mit.edu)" );
            didSaveSucceed *= fileHandler.addAtt( currFile, "GenerationDate", buffer );
            didSaveSucceed *= fileHandler.addAtt( currFile, "Format", "NetCDF-4" );

            double value = 0.0E+00;
            value = input.temperature_K();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Temperature"    , 1, "float", "K"  , "Ambient Temperature" );
            value = input.pressure_Pa();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Pressure"       , 1, "float", "hPa", "Ambient Pressure" );
            value = airDens;
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Air Density"    , 1, "float", "molecule / cm ^ 3", "Molecular density" );
            value = input.relHumidity_w();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "RHW"            , 1, "float", "-"  , "Ambient Rel. Humidity w.r.t water" );
            value = relHumidity_i;
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "RHI"            , 1, "float", "-"  , "Ambient Rel. Humidity w.r.t ice" );
            value = input.dayGMT();
            didSaveSucceed *= fileHandler.addConst( currFile, &value, "Day GMT"        , 1, "int"  , "-"  , "Emission day" );
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

            didSaveSucceed *= fileHandler.addVar( currFile, &(ambientData.cosSZA)[0], "CSZA", timeDim, "float", "-", "Cosine of the solar zenith angle" );

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
            std::vector<std::vector<double>> ringAverage = ringSpecies.RingAverage( ringArea, totArea );

            const unsigned int NT = timeArray.size();
            std::vector<double> plumeData( NT, 0.0E+00 );

            /* Start saving species ... */

            /* Define conversion factor */
            double scalingFactor;
            char charSpc[30];
            char charName[50];
            char charUnit[20];
               

            #if ( DO_SAVE_SO2 )

                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

                for ( unsigned int iNt = 0; iNt < NT; iNt++ )
                    plumeData[iNt] = ringAverage[iNt][ind_SO2];
               
                #if ( SAVE_TO_DOUBLE )
                    spcArray = util::vect2double( plumeData      , NT, scalingFactor );
                    ambArray = util::vect2double( ambientData.SO2 , NT, scalingFactor );
                #else
                    spcArray = util::vect2float ( plumeData      , NT, scalingFactor );
                    ambArray = util::vect2float ( ambientData.SO2 , NT, scalingFactor );
                #endif /* SAVE_TO_DOUBLE */


                strncpy( charSpc, "SO2_Plume", sizeof(charSpc) );
                strncpy( charName, "SO2 plume-averaged mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName ); 

                strncpy( charSpc, "SO2_Ambient", sizeof(charSpc) );
                strncpy( charName, "SO2 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName ); 

                delete[] spcArray; spcArray = NULL;
                delete[] ambArray; ambArray = NULL;
                delete[] adjArray; adjArray = NULL;

            #endif /* DO_SAVE_SO2 */

            #if ( DO_SAVE_O3 )

                scalingFactor = TO_PPB;
                strncpy( charUnit, "ppb", sizeof(charUnit) );

                for ( unsigned int iNt = 0; iNt < NT; iNt++ )
                    plumeData[iNt] = ringAverage[iNt][ind_O3];
               
                #if ( SAVE_TO_DOUBLE )
                    spcArray = util::vect2double( plumeData      , NT, scalingFactor );
                    adjArray = util::vect2double( adjointData.O3 , NT, scalingFactor );
                    ambArray = util::vect2double( ambientData.O3 , NT, scalingFactor );
                #else
                    spcArray = util::vect2float ( plumeData      , NT, scalingFactor );
                    adjArray = util::vect2float ( adjointData.O3 , NT, scalingFactor );
                    ambArray = util::vect2float ( ambientData.O3 , NT, scalingFactor );
                #endif /* SAVE_TO_DOUBLE */


                strncpy( charSpc, "O3_Plume", sizeof(charSpc) );
                strncpy( charName, "O3 plume-averaged mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName ); 

                strncpy( charSpc, "O3_Adjoint", sizeof(charSpc) );
                strncpy( charName, "O3 optimized mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(adjArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName ); 

                strncpy( charSpc, "O3_Ambient", sizeof(charSpc) );
                strncpy( charName, "O3 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName ); 

                delete[] spcArray; spcArray = NULL;
                delete[] ambArray; ambArray = NULL;
                delete[] adjArray; adjArray = NULL;

            #endif /* DO_SAVE_O3 */

            #if ( DO_SAVE_NO )

                scalingFactor = TO_PPT;
                strncpy( charUnit, "ppt", sizeof(charUnit) );

                for ( unsigned int iNt = 0; iNt < NT; iNt++ )
                    plumeData[iNt] = ringAverage[iNt][ind_NO];


                #if ( SAVE_TO_DOUBLE )
                    spcArray = util::vect2double( plumeData      , NT, scalingFactor );
                    adjArray = util::vect2double( adjointData.NO , NT, scalingFactor );
                    ambArray = util::vect2double( ambientData.NO , NT, scalingFactor );
                #else
                    spcArray = util::vect2float ( plumeData      , NT, scalingFactor );
                    adjArray = util::vect2float ( adjointData.NO , NT, scalingFactor );
                    ambArray = util::vect2float ( ambientData.NO , NT, scalingFactor );
                #endif /* SAVE_TO_DOUBLE */


                strncpy( charSpc, "NO_Plume", sizeof(charSpc) );
                strncpy( charName, "NO plume-averaged mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName ); 

                strncpy( charSpc, "NO_Adjoint", sizeof(charSpc) );
                strncpy( charName, "NO optimized mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(adjArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName ); 

                strncpy( charSpc, "NO_Ambient", sizeof(charSpc) );
                strncpy( charName, "NO ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName ); 

                delete[] spcArray; spcArray = NULL;
                delete[] ambArray; ambArray = NULL;
                delete[] adjArray; adjArray = NULL;

            #endif /* DO_SAVE_NO */

            #if ( DO_SAVE_NO2 )

                scalingFactor = TO_PPT;
                strncpy( charUnit, "ppt", sizeof(charUnit) );

                for ( unsigned int iNt = 0; iNt < NT; iNt++ )
                    plumeData[iNt] = ringAverage[iNt][ind_NO2];


                #if ( SAVE_TO_DOUBLE )
                    spcArray = util::vect2double( plumeData      , NT, scalingFactor );
                    adjArray = util::vect2double( adjointData.NO2, NT, scalingFactor );
                    ambArray = util::vect2double( ambientData.NO2, NT, scalingFactor );
                #else
                    spcArray = util::vect2float ( plumeData      , NT, scalingFactor );
                    adjArray = util::vect2float ( adjointData.NO2, NT, scalingFactor );
                    ambArray = util::vect2float ( ambientData.NO2, NT, scalingFactor );
                #endif /* SAVE_TO_DOUBLE */


                strncpy( charSpc, "NO2_Plume", sizeof(charSpc) );
                strncpy( charName, "NO2 plume-averaged mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName ); 

                strncpy( charSpc, "NO2_Adjoint", sizeof(charSpc) );
                strncpy( charName, "NO2 optimized mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(adjArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName ); 

                strncpy( charSpc, "NO2_Ambient", sizeof(charSpc) );
                strncpy( charName, "NO2 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName ); 

                delete[] spcArray; spcArray = NULL;
                delete[] ambArray; ambArray = NULL;
                delete[] adjArray; adjArray = NULL;

            #endif /* DO_SAVE_NO2 */

            #if ( DO_SAVE_HNO3 )

                scalingFactor = TO_PPT;
                strncpy( charUnit, "ppt", sizeof(charUnit) );

                for ( unsigned int iNt = 0; iNt < NT; iNt++ )
                    plumeData[iNt] = ringAverage[iNt][ind_HNO3];


                #if ( SAVE_TO_DOUBLE )
                    spcArray = util::vect2double( plumeData       , NT, scalingFactor );
                    adjArray = util::vect2double( adjointData.HNO3, NT, scalingFactor );
                    ambArray = util::vect2double( ambientData.HNO3, NT, scalingFactor );
                #else
                    spcArray = util::vect2float ( plumeData       , NT, scalingFactor );
                    adjArray = util::vect2float ( adjointData.HNO3, NT, scalingFactor );
                    ambArray = util::vect2float ( ambientData.HNO3, NT, scalingFactor );
                #endif /* SAVE_TO_DOUBLE */


                strncpy( charSpc, "HNO3_Plume", sizeof(charSpc) );
                strncpy( charName, "HNO3 plume-averaged mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(spcArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName ); 

                strncpy( charSpc, "HNO3_Adjoint", sizeof(charSpc) );
                strncpy( charName, "HNO3 optimized mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(adjArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName ); 

                strncpy( charSpc, "HNO3_Ambient", sizeof(charSpc) );
                strncpy( charName, "HNO3 ambient mixing ratio", sizeof(charName) );
                didSaveSucceed *= fileHandler.addVar( currFile, &(ambArray)[0], (const char*)charSpc, timeDim, outputType, (const char*)charUnit, (const char*)charName ); 

                delete[] spcArray; spcArray = NULL;
                delete[] ambArray; ambArray = NULL;
                delete[] adjArray; adjArray = NULL;

            #endif /* DO_SAVE_HNO3 */

            #if ( DO_SAVE_NOx )

                scalingFactor = TO_PPT;
                strncpy( charUnit, "ppt", sizeof(charUnit) );

                for ( unsigned int iNt = 0; iNt < NT; iNt++ )
                    plumeData[iNt] = ringAverage[iNt][ind_NO] + ringAverage[iNt][ind_NO2];

                std::vector<double> ambientNOx = util::add1D( ambientData.NO, ambientData.NO2 );
                std::vector<double> adjointNOx = util::add1D( adjointData.NO, adjointData.NO2 );
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

            #endif /* DO_SAVE_NOx */
            

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

}

/* End of Save.cpp */
