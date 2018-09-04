/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* ReadParameters Program File                                      */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : ReadParameters.cpp                        */
/* Working directory    : /home/fritzt/APCEMM-SourceCode            */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <iostream>
#include <vector>

std::vector<std::vector<double> > CombVec( std::vector<double>& temperature_K, \
                                           std::vector<double>& pressure_Pa,   \
                                           std::vector<double>& relHumidity_w, \
                                           std::vector<double>& longitude_deg, \
                                           std::vector<double>& latitude_deg );

std::vector<std::vector<double> > ReadParameters( )
{

    std::cout << "Reading in parameters..." << std::endl;

    std::vector<double> temperature_K, pressure_Pa, relHumidity_w, longitude_deg, latitude_deg;
    std::vector<std::vector<double> > parameters;

    /* Temperature array in [K] */

    temperature_K.push_back(210.01); // = {210.0, 220.0};
    //temperature_K.push_back(235.2);

    /* Pressure array in [Pa] */

    pressure_Pa.push_back(240.0E2);
    //pressure_Pa.push_back(440.0E2);

    /* Relative humidity w.r.t liquid water array in [\%] */

    relHumidity_w.push_back(30.0);

    /* Longitude array expressed in deg */

    longitude_deg.push_back(-15.0);

    /* Latitude array expressed in deg */

    latitude_deg.push_back(60.0);

    parameters = CombVec( temperature_K, \
                          pressure_Pa,   \
                          relHumidity_w, \
                          longitude_deg, \
                          latitude_deg );

    /* For debug */
    if ( false ) {
        unsigned int i, j;
        std::cout.precision(2);
        for ( i = 0; i < parameters.size(); i++ ) {
            for ( j = 0; j < parameters[0].size(); j++) {
                std::cout << parameters[i][j] << " ";
            }
            std::cout << "" << std::endl;
        }
    }

    return parameters;

} /* End of ReadParameters */
