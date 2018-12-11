/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Meteorology Program File                                         */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 11/2/2018                                 */
/* File                 : Meteorology.cpp                           */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "Core/Meteorology.hpp"

Meteorology::Meteorology( ) :
    LOAD( 0 )
{

    /* Default constructor */

} /* End of Meteorology::Meteorology */

Meteorology::Meteorology( const bool loadFile, \
                          const Mesh &m, \
                          const double temperature_K, \
                          const double altitude, \
                          const double LapseRate, \
                          const bool DBG ) : 
    LOAD ( loadFile )
{

    /* Constructor */

    std::vector<double> X = m.x();
    std::vector<double> Y = m.y();
    
    alt.assign( Y.size(), 0.0E+00 );
    press.assign( Y.size(), 0.0E+00 );

    for ( unsigned int iNy = 0; iNy < Y.size(); iNy++ )
        temp.push_back( std::vector<double>( X.size(), 0.0E+00 ) );

    for ( unsigned int iNy = 0; iNy < Y.size(); iNy++ )
        H2O.push_back( std::vector<double>( X.size(), 0.0E+00 ) );

    if ( LOAD ) {

    /* !@#$ */

    } else { 

        for ( unsigned int iNy = 0; iNy < Y.size(); iNy++ )
            alt[iNy] = altitude + Y[iNy];
        met::ISA( alt, press );
        /* User defined fields can be set here ! */
        for ( unsigned int iNy = 0; iNy < Y.size(); iNy++ ) {
            temp[iNy][0] = temperature_K + Y[iNy] * LapseRate;
            /* Unit check: Y [m] * LapseRate [K/m] */
            H2O[iNy][0] = 0.0E+00;
            for ( unsigned int iNx = 0; iNx < X.size(); iNx++ ) {
                temp[iNy][iNx] = temp[iNy][iNx];
                H2O[iNy][iNx] = H2O[iNy][0];
            }
        }

    }

    if ( DBG ) {
        std::cout << "\n DEBUG : Meteorology\n";
        std::cout << "         Grid number | Altitude [km] | Pressure [hPa] | Temperature [K] | H2O [-]\n";
        for ( unsigned int iNy = Y.size() - 1; iNy --> 0; ) {
            std::cout << "         " << std::setw(11) << iNy << " | " << std::setw(13) << alt[iNy] * 1.00E-03 \
                      << " | " << std::setw(14) << press[iNy] * 1.00E-02 << " | " << std::setw(15) << temp[iNy][0] \
                      << " | " << std::setw(9) << H2O[iNy][0] << "\n";
        }
        std::cout << "\n";
    }

} /* End of Meteorology::Meteorology */

Meteorology::Meteorology( const Meteorology &met ) :
    LOAD( met.LOAD )
{

    alt = met.alt;
    press = met.press;
    temp = met.temp;
    H2O = met.H2O;

} /* End of Meteorology::Meteorology */

Meteorology::~Meteorology( )
{

    /* Destructor */

} /* End of Meteorology::~Meteorology */

/* End of Meteorology.cpp */
