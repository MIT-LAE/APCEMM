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
    LAPSERATE( -3.0E-03 ),
    CSTDEPTH( 0 ),
    DEPTH( 0.0E+00 ),
    DIURNAL_AMPL( 0.0E+00 ),
    DIURNAL_PHASE( 0.0E+00 )
{

    /* Default constructor */

} /* End of Meteorology::Meteorology */

Meteorology::Meteorology( const bool loadFile,        \
                          const double solarTime_h,   \
                          const Mesh &m,              \
                          const double temperature_K, \
                          const double altitude,      \
                          const double LapseRate,     \
                          const bool cstDepth,        \
                          const double depth,         \
                          const bool DBG ) : 
    LOAD ( loadFile ),
    TEMPERATURE( temperature_K ),
    ALTITUDE( altitude ),
    LAPSERATE( LapseRate ),
    CSTDEPTH( cstDepth ),
    DEPTH( depth )
{

    /* Constructor */

    std::vector<double> X = m.x();
    std::vector<double> Y = m.y();
    
    alt_.assign( Y.size(), 0.0E+00 );
    press_.assign( Y.size(), 0.0E+00 );

    for ( unsigned int jNy = 0; jNy < Y.size(); jNy++ )
        temp_.push_back( std::vector<double>( X.size(), 0.0E+00 ) );

    for ( unsigned int jNy = 0; jNy < Y.size(); jNy++ )
        H2O_.push_back( std::vector<double>( X.size(), 0.0E+00 ) );
    
    /* The diurnal temperature amplitude represents one-half of the 
     * diurnal temperature range, and the phase represents the hour
     * of maximum temperature. We define the diurnal temperature variation
     * as a cosine wave centered on the phase hour. 
     * The data on the diurnal temperature amplitude and phase is gathered
     * from: 
     * Seidel, D. J., M. Free, and J. Wang (2005), Diurnal cycle of 
     * upper-air temperature estimated from radiosondes,J. Geophys.Res.,
     * 110, D09102, doi:10.1029/2004JD005526.*/

    DIURNAL_AMPL  = 0.1 ; /* [K] */
    DIURNAL_PHASE = 14.0; /* [hrs] */

    diurnalPert = DIURNAL_AMPL * cos( 2.0E+00 * physConst::PI * ( solarTime_h - DIURNAL_PHASE ) / 24.0E+00 );

    if ( LOAD ) {

    /* !@#$ */

    } else { 

        for ( unsigned int jNy = 0; jNy < Y.size(); jNy++ )
            alt_[jNy] = altitude + Y[jNy];
        met::ISA( alt_, press_ );
        /* User defined fields can be set here ! */

//        /* X-invariant met field: */
//        for ( unsigned int jNy = 0; jNy < Y.size(); jNy++ ) {
//            temp_[jNy][X.size()/2] = temperature_K + Y[jNy] * LapseRate + diurnalPert;
//            /* Unit check: Y [m] * LapseRate [K/m] */
//            H2O_[jNy][0] = 0.0E+00;
//            for ( unsigned int iNx = 0; iNx < X.size(); iNx++ ) {
//                if ( ( X[iNx] < -30.0E+03 ) || ( X[iNx] > 30.0E+03 ) ) {
//                    /* Add a 1K depression in that region */
//                    temp_[jNy][iNx] = temp_[jNy][X.size()/2] + 1;
//                    H2O_[jNy][iNx]  = H2O_[jNy][X.size()/2];
//                } else {
//                    temp_[jNy][iNx] = temp_[jNy][X.size()/2];
//                    H2O_[jNy][iNx]  = H2O_[jNy][X.size()/2];
//                }
//            }
//        }

        /* RHi bubble centered on x = 0
         * Add a ~0.5K depression to ambient air */
        const double DELTAT = 1.0;
        const double BACKGT = temperature_K + DELTAT + diurnalPert; /* [K] */
        double TOP, BOT;
        double LEFT, RIGHT;

        TOP = 600.0;
        if ( DEPTH == 0.0E+00 )
            BOT = 200.0;
        else
            BOT = DEPTH;
        //LEFT = 7.50E+03;
        LEFT = 30.00E+03;
        RIGHT= LEFT;

        for ( unsigned int jNy = 0; jNy < Y.size(); jNy++ ) {
            H2O_[jNy][0] = 0.0E+00;
            for ( unsigned int iNx = 0; iNx < X.size(); iNx++ ) {
                temp_[jNy][iNx] = BACKGT + LapseRate * Y[jNy] \
                               + ( temperature_K - BACKGT )   \
                               * ( 1.0 - 0.5 * ( std::tanh( ( X[iNx] - LEFT ) / 1.0E+03 ) + 1.0 )) \
                               * ( 0.0 + 0.5 * ( std::tanh( ( X[iNx] + RIGHT) / 1.0E+03 ) + 1.0 ));
//                               * ( 1.0 - 0.5 * ( std::tanh( ( Y[jNy] - TOP  ) / 1.0E+02 ) + 1.0 )) \
//                               * ( 0.0 + 0.5 * ( std::tanh( ( Y[jNy] + BOT  ) / 1.0E+02 ) + 1.0 ))
                H2O_[jNy][iNx]  = H2O_[jNy][0];
            }
        }

    }

    if ( DBG ) {
        std::cout << "\n DEBUG : Meteorology\n";
        std::cout << "         Grid number | Altitude [km] | Pressure [hPa] | Temperature [K] | H2O [-]\n";
        for ( unsigned int jNy = Y.size() - 1; jNy --> 0; ) {
            std::cout << "         " << std::setw(11) << jNy << " | " << std::setw(13) << alt_[jNy] * 1.00E-03 \
                      << " | " << std::setw(14) << press_[jNy] * 1.00E-02 << " | " << std::setw(15) << temp_[jNy][0] \
                      << " | " << std::setw(9) << H2O_[jNy][0] << "\n";
        }
        std::cout << "\n";
    }

} /* End of Meteorology::Meteorology */

Meteorology::Meteorology( const Meteorology &met ) :
    LOAD( met.LOAD ),
    TEMPERATURE( met.TEMPERATURE ),
    ALTITUDE( met.ALTITUDE ),
    LAPSERATE( met.LAPSERATE ),
    CSTDEPTH( met.CSTDEPTH ),
    DEPTH( met.DEPTH )
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

void Meteorology::Update( const double solarTime_h, const Mesh &m, \ 
                          const double dTrav_x, const double dTrav_y )
{

    std::vector<double> X = m.x();
    std::vector<double> Y = m.y();

    for ( unsigned int jNy = 0; jNy < Y.size(); jNy++ )
        alt_[jNy] = ALTITUDE + Y[jNy] + dTrav_y;
    met::ISA( alt_, press_ );
    
    diurnalPert = DIURNAL_AMPL * cos( 2.0E+00 * physConst::PI * ( solarTime_h - DIURNAL_PHASE ) / 24.0E+00 );

    /* User defined fields can be set here ! */

    /* X-invariant met field: */
//    for ( unsigned int jNy = 0; jNy < Y.size(); jNy++ ) {
//        temp_[jNy][X.size()/2] = TEMPERATURE + ( Y[jNy] + dTrav_y ) * LAPSERATE + diurnalPert;
//        /* Unit check: Y [m] * LapseRate [K/m] */
//        H2O_[jNy][0] = 0.0E+00;
//        for ( unsigned int iNx = 0; iNx < X.size(); iNx++ ) {
//            if ( ( X[iNx] < -30.0E+03 ) || ( X[iNx] > 30.0E+03 ) ) {
//                /* Add a 1K depression in that region */
//                temp_[jNy][iNx] = temp_[jNy][X.size()/2] + 1;
//                H2O_[jNy][iNx]  = H2O_[jNy][X.size()/2];
//            } else {
//                temp_[jNy][iNx] = temp_[jNy][X.size()/2];
//                H2O_[jNy][iNx]  = H2O_[jNy][X.size()/2];
//            }
//        }
//    }
        
    /* RHi bubble centered on x = 0
     * Add a ~0.5K depression to ambient air */
    const double DELTAT = 1.0;
    const double BACKGT = TEMPERATURE + DELTAT + diurnalPert; /* [K] */
    double TOP, BOT;
    double LEFT, RIGHT;

    TOP = 600.0;
    if ( DEPTH == 0.0E+00 )
        BOT = 200.0;
    else
        BOT = DEPTH;
    //LEFT = 7.50E+03;
    LEFT = 30.00E+03;
    RIGHT= LEFT;

    for ( unsigned int jNy = 0; jNy < Y.size(); jNy++ ) {
        H2O_[jNy][0] = 0.0E+00;
        for ( unsigned int iNx = 0; iNx < X.size(); iNx++ ) {
            temp_[jNy][iNx] = BACKGT + LAPSERATE * ( Y[jNy] + dTrav_y ) \
                           + ( TEMPERATURE - BACKGT ) \
                           * ( 1.0 - 0.5 * ( std::tanh( ( X[iNx] + dTrav_x - LEFT ) / 1.0E+03 ) + 1.0 )) \
                           * ( 0.0 + 0.5 * ( std::tanh( ( X[iNx] + dTrav_x + RIGHT) / 1.0E+03 ) + 1.0 )); 
//                           * ( 1.0 - 0.5 * ( std::tanh( ( Y[jNy] + dTrav_y - TOP  ) / 1.0E+02 ) + 1.0 )) \
//                           * ( 0.0 + 0.5 * ( std::tanh( ( Y[jNy] + dTrav_y + BOT  ) / 1.0E+02 ) + 1.0 ))
            H2O_[jNy][iNx]  = H2O_[jNy][0];
        }
    }


} /* End of Meteorology::UpdateMet */

/* End of Meteorology.cpp */
