/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Main Program File                                                */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : Main.cpp                                  */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <string>
#include <vector>

using namespace std;

void PrintMessage( bool doPrint );
vector<vector<double> > ReadParameters( );
int PlumeModel( double temperature_K, double pressure_Pa, \
                 double relHumidity_w, double longitude_deg, \
                 double latitude_deg );

int main( int , char* [] )
{

    bool doPrint = false;
    vector<vector<double> > parameters;
    unsigned int nCases;

    /* Print welcome message */
    PrintMessage( doPrint );

    /* Read in parameters */
    parameters = ReadParameters();

    nCases  = 1; //parameters[0].size();

    for ( unsigned int i = 0; i < nCases; i++ ) {
        double temperature_K = parameters[0][i];
        double pressure_Pa = parameters[1][i];
        double relHumidity_w = parameters[2][i];
        double longitude_deg = parameters[3][i];
        double latitude_deg = parameters[4][i];

        unsigned int model;
        /*
        * model = 0 -> Box     Model
        * model = 1 -> Plume   Model (APCEMM)
        * model = 2 -> Adjoint Model
        * model = 3 -> Box + Plume Model
        */
        model = 1;

        int iERR;

        if (model == 1) {
            iERR = PlumeModel( temperature_K, pressure_Pa, \
                               relHumidity_w, longitude_deg, \
                               latitude_deg );

            if ( iERR < 0 ) {
                cout.precision(3);
                cout << "APCEMM Case: " << i << " failed." << endl;
                cout << "Error: " << iERR << endl;
                cout << fixed;
                cout << setprecision(3);
                cout << "T   : " << setw(8) << temperature_K << " [K]" << endl;
                cout << "P   : " << setw(8) << pressure_Pa/((double) 100.0) << " [hPa]" << endl;
                cout << "RH_w: " << setw(8) << relHumidity_w << " [%]" << endl;
                cout << "LON : " << setw(8) << longitude_deg << " [deg]" << endl;
                cout << "LAT : " << setw(8) << latitude_deg << " [deg]" << endl;
            }
            else {
                cout << "APCEMM Case: " << i << " completed." << endl;
            }
        }

    }

    return 0;

}
