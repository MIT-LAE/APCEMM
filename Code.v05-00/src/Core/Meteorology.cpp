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
    TEMPERATURE( 220 ),
    PRESSURE( 220.0E+02 ),
    RHI( 0.0E+00 ),
    ALTITUDE( 10.5E+03 ),
    DIURNAL_AMPL( 0.0E+00 ),
    DIURNAL_PHASE( 0.0E+00 )
{

    /* Default constructor */

} /* End of Meteorology::Meteorology */

Meteorology::Meteorology( const OptInput &USERINPUT,      \
                          const RealDouble solarTime_h,   \
                          const Mesh &m,                  \
                          const RealDouble temperature_K, \
                          const RealDouble pressure_Pa,   \
                          const RealDouble relHumidity_i, \
                          const bool DBG ) : 
    TEMPERATURE( temperature_K ),
    PRESSURE( pressure_Pa ),
    RHI( relHumidity_i )
{

    /* Constructor */

    Vector_1D X = m.x();
    Vector_1D Y = m.y();
    
    alt_.assign( Y.size(), 0.0E+00 );
    press_.assign( Y.size(), 0.0E+00 );

    for ( UInt jNy = 0; jNy < Y.size(); jNy++ )
        temp_.push_back( Vector_1D( X.size(), 0.0E+00 ) );

    for ( UInt jNy = 0; jNy < Y.size(); jNy++ )
        H2O_.push_back( Vector_1D( X.size(), 0.0E+00 ) );

    /* If imposing a fixed moist layer depth, the temperature lapse rate gets
     * overwritten */
    LAPSERATE = USERINPUT.MET_LAPSERATE; /* [K/m] */

    /* If a fixed depth of the moist layer is imposed, compute the temperature 
     * lapse rate to have the prescribed RH at flight level and a 100% RH at 
     * MET_DEPTH (expressed in m) */
    if ( USERINPUT.MET_FIXDEPTH )
        LAPSERATE = met::ComputeLapseRate( TEMPERATURE, RHI, \
                                           USERINPUT.MET_DEPTH );

    /* The diurnal temperature amplitude represents one-half of the 
     * diurnal temperature range, and the phase represents the hour
     * of maximum temperature. We define the diurnal temperature variation
     * as a cosine wave centered on the phase hour. 
     * The data on the diurnal temperature amplitude and phase is gathered
     * from: 
     * Seidel, D. J., M. Free, and J. Wang (2005), Diurnal cycle of 
     * upper-air temperature estimated from radiosondes,J. Geophys.Res.,
     * 110, D09102, doi:10.1029/2004JD005526.*/

    if ( USERINPUT.MET_DIURNAL ) {
        /* Based on values around 220 hPa */
        DIURNAL_AMPL  = 0.1 ; /* [K] */
        DIURNAL_PHASE = 14.0; /* [hrs] */
    } else {
        DIURNAL_AMPL  = 0.0 ; /* [K] */
        DIURNAL_PHASE = 12.0; /* [hrs] */
    }

    diurnalPert = DIURNAL_AMPL * cos( 2.0E+00 * physConst::PI * ( solarTime_h - DIURNAL_PHASE ) / 24.0E+00 );

    if ( USERINPUT.MET_LOADMET ) {

        /* Met is loaded in */
        TYPE = 0;

        if ( USERINPUT.MET_LOADTEMP ) {
            /* !@#$ */
        }

        if ( USERINPUT.MET_LOADH2O ) {
            /* !@#$ */
        }

    } else { 

        met::ISA_pAlt( ALTITUDE, PRESSURE );

        for ( UInt jNy = 0; jNy < Y.size(); jNy++ )
            alt_[jNy] = ALTITUDE + Y[jNy];
        met::ISA( alt_, press_ );

        /* User defined fields can be set here ! */

        /* Set temperature input type if user-defined */
        TYPE = 2;


        if ( TYPE == 1 ) {
            /* (quasi-)X-invariant met field: */

            LEFT = 35.0E+03;
            RIGHT= LEFT;

            DELTAT = 1.0; /* Add a ~1.0K depression far away */
            for ( UInt jNy = 0; jNy < Y.size(); jNy++ ) {
                temp_[jNy][X.size()/2] = TEMPERATURE + Y[jNy] * LAPSERATE + diurnalPert;
                /* Unit check: Y [m] * LapseRate [K/m] */
                H2O_[jNy][0] = 0.0E+00;
                for ( UInt iNx = 0; iNx < X.size(); iNx++ ) {
                    if ( ( X[iNx] < -LEFT ) || ( X[iNx] > RIGHT ) ) {
                        /* Add a 1K depression in that region */
                        temp_[jNy][iNx] = temp_[jNy][X.size()/2] + DELTAT;
                        H2O_[jNy][iNx]  = H2O_[jNy][X.size()/2];
                    } else {
                        temp_[jNy][iNx] = temp_[jNy][X.size()/2];
                        H2O_[jNy][iNx]  = H2O_[jNy][X.size()/2];
                    }
                }
            }

        } else if ( TYPE == 2 ) {

            /* RHi bubble centered on x = 0
             * Add a ~3.0K depression to ambient air */
            DELTAT = 3.0;

            TOP = 600.0;
            if ( USERINPUT.MET_FIXDEPTH )
                BOT = USERINPUT.MET_DEPTH;
            else
                BOT = 200.0;
            LEFT = 30.00E+03;
            RIGHT= LEFT;

            for ( UInt jNy = 0; jNy < Y.size(); jNy++ ) {
                H2O_[jNy][0] = 0.0E+00;
                for ( UInt iNx = 0; iNx < X.size(); iNx++ ) {
                    temp_[jNy][iNx] = TEMPERATURE + DELTAT + diurnalPert + LAPSERATE * Y[jNy] \
                                   - ( DELTAT ) \
                                   * ( 1.0 - 0.5 * ( std::tanh( ( X[iNx] - LEFT ) / 1.0E+03 ) + 1.0 )) \
                                   * ( 0.0 + 0.5 * ( std::tanh( ( X[iNx] + RIGHT) / 1.0E+03 ) + 1.0 ));
    //                               * ( 1.0 - 0.5 * ( std::tanh( ( Y[jNy] - TOP  ) / 1.0E+02 ) + 1.0 )) \
    //                               * ( 0.0 + 0.5 * ( std::tanh( ( Y[jNy] + BOT  ) / 1.0E+02 ) + 1.0 ))
                    H2O_[jNy][iNx]  = H2O_[jNy][0];
                }
            }

        } else {

            std::cout << " In Meteorology::Meteorology: Undefined value for user-input type: ";
            std::cout << TYPE << std::endl;
            exit(-1);

        }
    }

    if ( DBG ) {
        std::cout << "\n DEBUG : Meteorology\n";
        std::cout << "         Grid number | Altitude [km] | Pressure [hPa] | Temperature [K] | H2O [-]\n";
        for ( UInt jNy = Y.size() - 1; jNy --> 0; ) {
            std::cout << "         " << std::setw(11) << jNy << " | " << std::setw(13) << alt_[jNy] * 1.00E-03 \
                      << " | " << std::setw(14) << press_[jNy] * 1.00E-02 << " | " << std::setw(15) << temp_[jNy][0] \
                      << " | " << std::setw(9) << H2O_[jNy][0] << "\n";
        }
        std::cout << "\n";
    }

} /* End of Meteorology::Meteorology */

Meteorology::Meteorology( const Meteorology &met ) :
    TEMPERATURE( met.TEMPERATURE ),
    PRESSURE( met.PRESSURE ),
    RHI( met.RHI ),
    ALTITUDE( met.ALTITUDE )
{

    TYPE          = met.TYPE;
    LAPSERATE     = met.LAPSERATE;
    DIURNAL_AMPL  = met.DIURNAL_AMPL;
    DIURNAL_PHASE = met.DIURNAL_PHASE;
    diurnalPert   = met.diurnalPert;
    DELTAT        = met.DELTAT;
    TOP           = met.TOP;
    BOT           = met.BOT;
    LEFT          = met.LEFT;
    RIGHT         = met.RIGHT;
    alt_          = met.alt_;
    press_        = met.press_;
    temp_         = met.temp_;
    H2O_          = met.H2O_;

} /* End of Meteorology::Meteorology */

Meteorology::~Meteorology( )
{

    /* Destructor */

} /* End of Meteorology::~Meteorology */

void Meteorology::Update( const RealDouble solarTime_h, const Mesh &m, \ 
                          const RealDouble dTrav_x, const RealDouble dTrav_y )
{

    UInt iNx = 0;
    UInt jNy = 0;

    Vector_1D X = m.x();
    Vector_1D Y = m.y();

    diurnalPert = DIURNAL_AMPL * cos( 2.0E+00 * physConst::PI * ( solarTime_h - DIURNAL_PHASE ) / 24.0E+00 );

#pragma omp parallel for                \
            if      ( !PARALLEL_CASES ) \
            default ( shared          ) \
            private ( jNy             ) \
            schedule( dynamic, 1      )
    for ( jNy = 0; jNy < Y.size(); jNy++ )
        alt_[jNy] = ALTITUDE + Y[jNy] + dTrav_y;
    met::ISA( alt_, press_ );
    
    /* User defined fields can be set here ! */

    if ( TYPE == 1 ) {
        
        /* (quasi-)X-invariant met field: */

#pragma omp parallel for                \
            if      ( !PARALLEL_CASES ) \
            default ( shared          ) \
            private ( iNx, jNy        ) \
            schedule( dynamic, 1      )
        for ( jNy = 0; jNy < Y.size(); jNy++ ) {
            temp_[jNy][X.size()/2] = TEMPERATURE \
                                   + diurnalPert + LAPSERATE * ( Y[jNy] + dTrav_y );
            /* Unit check: Y [m] * LapseRate [K/m] */
            H2O_[jNy][0] = 0.0E+00;
            for ( iNx = 0; iNx < X.size(); iNx++ ) {
                if ( ( X[iNx] < -LEFT ) || ( X[iNx] > RIGHT ) ) {
                    /* Add a 1K depression in that region */
                    temp_[jNy][iNx] = temp_[jNy][X.size()/2] + DELTAT;
                    H2O_[jNy][iNx]  = H2O_[jNy][X.size()/2];
                } else {
                    temp_[jNy][iNx] = temp_[jNy][X.size()/2];
                    H2O_[jNy][iNx]  = H2O_[jNy][X.size()/2];
                }
            }
        }

    } else if ( TYPE == 2 ) {

        /* RHi bubble centered on x = 0 */

#pragma omp parallel for                \
            if      ( !PARALLEL_CASES ) \
            default ( shared          ) \
            private ( iNx, jNy        ) \
            schedule( dynamic, 1      )
        for ( jNy = 0; jNy < Y.size(); jNy++ ) {
            H2O_[jNy][0] = 0.0E+00;
            for ( iNx = 0; iNx < X.size(); iNx++ ) {
                temp_[jNy][iNx] = TEMPERATURE \
                               +  DELTAT + diurnalPert + LAPSERATE * ( Y[jNy] + dTrav_y ) \
                               - ( DELTAT ) \
                               * ( 1.0 - 0.5 * ( std::tanh( ( X[iNx] - LEFT ) / 1.0E+03 ) + 1.0 )) \
                               * ( 0.0 + 0.5 * ( std::tanh( ( X[iNx] + RIGHT) / 1.0E+03 ) + 1.0 ));
//                               * ( 1.0 - 0.5 * ( std::tanh( ( Y[jNy] - TOP  ) / 1.0E+02 ) + 1.0 )) \
//                               * ( 0.0 + 0.5 * ( std::tanh( ( Y[jNy] + BOT  ) / 1.0E+02 ) + 1.0 ))
                H2O_[jNy][iNx]  = H2O_[jNy][0];
            }
        }

    } else {

        std::cout << " In Meteorology::Meteorology: Undefined value for user-input type: ";
        std::cout << TYPE << std::endl;
        exit(-1);

    }

} /* End of Meteorology::UpdateMet */

/* End of Meteorology.cpp */
