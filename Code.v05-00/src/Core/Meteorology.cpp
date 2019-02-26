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
    LOAD( 0 ),
    TEMPERATURE( 220 ),
    ALTITUDE( 10.5E+03 ),
    LAPSERATE( -3.0E-03 )
{

    /* Default constructor */

} /* End of Meteorology::Meteorology */

Meteorology::Meteorology( const bool loadFile, \
                          const Mesh &m, \
                          const double temperature_K, \
                          const double altitude, \
                          const double LapseRate, \
                          const bool DBG ) : 
    LOAD ( loadFile ),
    TEMPERATURE( temperature_K ),
    ALTITUDE( altitude ),
    LAPSERATE( LapseRate )
{

    /* Constructor */

    std::vector<double> X = m.x();
    std::vector<double> Y = m.y();
    
    alt_.assign( Y.size(), 0.0E+00 );
    press_.assign( Y.size(), 0.0E+00 );

    for ( unsigned int iNy = 0; iNy < Y.size(); iNy++ )
        temp_.push_back( std::vector<double>( X.size(), 0.0E+00 ) );

    for ( unsigned int iNy = 0; iNy < Y.size(); iNy++ )
        H2O_.push_back( std::vector<double>( X.size(), 0.0E+00 ) );

    if ( LOAD ) {

    /* !@#$ */

    } else { 

        for ( unsigned int iNy = 0; iNy < Y.size(); iNy++ )
            alt_[iNy] = altitude + Y[iNy];
        met::ISA( alt_, press_ );
        /* User defined fields can be set here ! */
        for ( unsigned int iNy = 0; iNy < Y.size(); iNy++ ) {
            temp_[iNy][0] = temperature_K + Y[iNy] * LapseRate;
            /* Unit check: Y [m] * LapseRate [K/m] */
            H2O_[iNy][0] = 0.0E+00;
            for ( unsigned int iNx = 1; iNx < X.size(); iNx++ ) {
                temp_[iNy][iNx] = temp_[iNy][0];
                H2O_[iNy][iNx]  = H2O_[iNy][0];
            }
        }

    }

    if ( DBG ) {
        std::cout << "\n DEBUG : Meteorology\n";
        std::cout << "         Grid number | Altitude [km] | Pressure [hPa] | Temperature [K] | H2O [-]\n";
        for ( unsigned int iNy = Y.size() - 1; iNy --> 0; ) {
            std::cout << "         " << std::setw(11) << iNy << " | " << std::setw(13) << alt_[iNy] * 1.00E-03 \
                      << " | " << std::setw(14) << press_[iNy] * 1.00E-02 << " | " << std::setw(15) << temp_[iNy][0] \
                      << " | " << std::setw(9) << H2O_[iNy][0] << "\n";
        }
        std::cout << "\n";
    }

} /* End of Meteorology::Meteorology */

Meteorology::Meteorology( const Meteorology &met ) :
    LOAD( met.LOAD ),
    TEMPERATURE( met.TEMPERATURE ),
    ALTITUDE( met.ALTITUDE ),
    LAPSERATE( met.LAPSERATE )
{

    alt_ = met.alt_;
    press_ = met.press_;
    temp_ = met.temp_;
    H2O_ = met.H2O_;

} /* End of Meteorology::Meteorology */

Meteorology::~Meteorology( )
{

    /* Destructor */

} /* End of Meteorology::~Meteorology */

void Meteorology::Update( const Mesh &m, const double dTrav_x, const double dTrav_y )
{

    std::vector<double> X = m.x();
    std::vector<double> Y = m.y();

    for ( unsigned int iNy = 0; iNy < Y.size(); iNy++ )
        alt_[iNy] = ALTITUDE + Y[iNy] + dTrav_y;
    met::ISA( alt_, press_ );
    /* User defined fields can be set here ! */
    for ( unsigned int iNy = 0; iNy < Y.size(); iNy++ ) {
        temp_[iNy][0] = TEMPERATURE + ( Y[iNy] + dTrav_y ) * LAPSERATE;
        /* Unit check: Y [m] * LapseRate [K/m] */
        H2O_[iNy][0] = 0.0E+00;
        for ( unsigned int iNx = 1; iNx < X.size(); iNx++ ) {
            temp_[iNy][iNx] = temp_[iNy][0];
            H2O_[iNy][iNx]  = H2O_[iNy][0];
        }
    }


} /* End of Meteorology::UpdateMet */

/* End of Meteorology.cpp */
