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
    SHEAR( 0.0E+00 ),
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
                          const RealDouble shear_pm, \
                          const bool DBG ) : 
    TEMPERATURE( temperature_K ),
    PRESSURE( pressure_Pa ),
    RHI( relHumidity_i ),
    SHEAR( shear_pm )
{

    /* Constructor */

    UInt iNx = 0;
    UInt jNy = 0;

    Vector_1D X = m.x();
    Vector_1D Y = m.y();

    alt_.assign( Y.size(), 0.0E+00 );
    press_.assign( Y.size(), 0.0E+00 );
    shear_.assign( Y.size(), 0.0E+00 );

    /* Vector of 0s */
    Vector_1D v1d( X.size(), 0.0E+00 );

    for ( jNy = 0; jNy < Y.size(); jNy++ ) {
        temp_.push_back( v1d );
        airDens_.push_back( v1d );
        H2O_.push_back( v1d );
    }

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
    std::cout << "loadmet=" << USERINPUT.MET_LOADMET << std::endl;

    diurnalPert = DIURNAL_AMPL * cos( 2.0E+00 * physConst::PI * ( solarTime_h - DIURNAL_PHASE ) / 24.0E+00 );

    if ( USERINPUT.MET_LOADMET ) {

        /* Met is loaded in */
        /* TYPE = 0; */
        TYPE = 3; /* TYPE = 0 means met for every time step */

        /* Set up for error messages */
        static const int NC_ERR = 2;
        NcError err( NcError::silent_nonfatal );


        /* Open the netcdf file for read access */
        NcFile dataFile( USERINPUT.MET_FILENAME.c_str(), NcFile::ReadOnly  );
	if ( !dataFile.is_valid() ) {
	    std::cout << "Netcdf file is not valid" << std::endl;
	}

        /* Identify the length of variables in input file */
        int var_len;
        int var_var_len = 1;
        NcVar* var_len_ncVar;
        if ( !(var_len_ncVar = dataFile.get_var( "var_length" )) ) {
            std::cout << "In Meteorology:: Meteorology: getting variable length pointer" << std::endl;
        }
        var_len_ncVar->get(&var_len, var_var_len);

        /* Extract pressure and altitude from input file */
        //int var_len = 5117;
        float altitude_user[var_len], pressure_user[var_len];
	altitude_store_.assign( var_len, 0.0E+00 );
        NcVar* altitude_ncVar;
        NcVar* pressure_ncVar;
        if ( !(altitude_ncVar = dataFile.get_var( "altitude" )) ) {
            std::cout << "In Meteorology:: Meteorology: getting altitude pointer" << std::endl;
        }
        if ( !(pressure_ncVar = dataFile.get_var( "pressure" )) ) {
            std::cout << "In Meteorology:: Meteorology: getting pressure pointer" << std::endl;
        }
        altitude_ncVar->get(altitude_user, var_len);
        pressure_ncVar->get(pressure_user, var_len);
        int i;
        for ( i = 0; i < var_len; i++ ) {
            pressure_user[i] *= 100.0;
            altitude_user[i] *= 1000.0;
	    altitude_store_[i] = altitude_user[i];
        }

        /* Identify closest pressure in input file from user-defined pressure */
        UInt i_Zp = met::nearestNeighbor( pressure_user, PRESSURE );
        pres_user = pressure_user[i_Zp];
        alt_user = met::linearInterp( pressure_user, altitude_user, PRESSURE );
        for ( UInt jNy = 0; jNy < Y.size(); jNy++ ) {
             // alt_[jNy] = altitude_user[i_Zp] + Y[jNy];
             alt_[jNy] = alt_user + Y[jNy];
        }


        float relhumid_user[var_len];
        float temperature_user[var_len];
	Vector_1D v1d_store( 8, 0.0E+00 );
        for ( UInt itime = 0; itime < var_len; itime++ ) {
            temperature_store_.push_back( v1d_store );
        }

        std::cout << "Interpolate T: " << USERINPUT.MET_INTERPTEMP << std::endl;
        std::cout << "Interpolate H2O: " << USERINPUT.MET_INTERPH2O << std::endl;

        if ( USERINPUT.MET_LOADTEMP ) {
            /* !@#$ */

            /* Define temperature input dimension */
            NcVar* temp_ncVar;

	    /* Give back a pointer to the requested NcVar */
	    if ( !(temp_ncVar = dataFile.get_var( "temperature" )) ) {
	        std::cout << "In Meteorology::Meteorology: getting temperature pointer" << std::endl;
	    }

            /* Extract the temperature data */

            if ( USERINPUT.MET_INTERPTEMP ) {
	        
		float temperature_store_temp[var_len][8];
                if ( !(temp_ncVar->get(&temperature_store_temp[0][0], var_len, 8)) ) {   
                    std::cout << "In Meteorology::Meteorology: extracting temperature" << std::endl;                                                              
                }

                /* Identify temperature at above pressure */
                temp_user = temperature_store_temp[i_Zp][0];
                
		/* Extract 2D temperature data into a 1D array */
                for ( UInt i = 0; i < var_len; i++ ) {
                    temperature_user[i] = temperature_store_temp[i][0];
		    for ( UInt itime = 0; itime < 8; itime++ ) {
                        temperature_store_[i][itime] = temperature_store_temp[i][itime];
		    }
                }

            }
            else {
                if ( !(temp_ncVar->get(temperature_user, var_len)) ) {
                    std::cout << "In Meteorology::Meteorology: extracting temperature" << std::endl;
                }
                /* Identify temperature at above pressure */
                temp_user = temperature_user[i_Zp];
            }

            /* Identify closest temperature to given pressure */
            /* Loop round each vertical layer to estimate temperature */
            for ( UInt jNy = 0; jNy < Y.size(); jNy++ ) {

                /* Find the closest values above and below the central pressure */
                UInt i_Z = met::nearestNeighbor( altitude_user, alt_[jNy] );

                /* Loop round horizontal coordinates to assign temperature */
                for ( UInt iNx = 0; iNx < X.size(); iNx++ ) {
                    temp_[jNy][iNx] = temperature_user[i_Z];
                }

            }

        }

        if ( USERINPUT.MET_LOADH2O ) {
            /* !@#$ */

            /* Define temperature input dimension */
            NcVar* relhumid_ncVar;

	    /* Give back a pointer to the requested NcVar */
	    if ( !(relhumid_ncVar = dataFile.get_var( "relative_humidity" )) ) {
	        std::cout << "In Meteorology::Meteorology: getting relative humidity pointer" << std::endl;
	    }

            /* Extract the temperature data */
            if ( !(relhumid_ncVar->get(relhumid_user, var_len)) ) {
                std::cout << "In Meteorology::Meteorology: extracting relative humidity" << std::endl;
            }

            /* Identify temperature at above pressure */
            RHw_user = relhumid_user[i_Zp];

            /* Identify closest temperature to given pressure */
            /* Loop round each vertical layer to estimate temperature */
            for ( UInt jNy = 0; jNy < Y.size(); jNy++ ) {

                /* Find the closest values above and below the central pressure */
                UInt i_Z = met::nearestNeighbor( altitude_user, alt_[jNy] );
                /* Loop round horizontal coordinates to assign temperature */
                for ( UInt iNx = 0; iNx < X.size(); iNx++ ) {
                    H2O_[jNy][iNx] = relhumid_user[i_Z]/((double) 100.00) *\
                                    physFunc::pSat_H2Ol( temp_[jNy][iNx] ) / ( physConst::kB * temp_[jNy][iNx] ) * 1.00E-06;
                }

            }

        }

        /* Identify the saturation depth */
        satdepth_user = 100; //met::satdepth_calc( relhumid_user, temperature_user, altitude_user, i_Zp, var_len );
        if ( satdepth_user != 1.0 ) {
            satdepth_user = satdepth_user - ( altitude_user[i_Zp]-alt_user );
        }

        /* Assign shear vector */
        for ( jNy = 0; jNy < Y.size(); jNy++ )
            shear_[jNy] = SHEAR;

    } else { 

        met::ISA_pAlt( ALTITUDE, PRESSURE );

        for ( jNy = 0; jNy < Y.size(); jNy++ )
            alt_[jNy] = ALTITUDE + Y[jNy];
        met::ISA( alt_, press_ );

        /* Assign shear vector */
        for ( jNy = 0; jNy < Y.size(); jNy++ )
            shear_[jNy] = SHEAR;

        /* Set temperature input type if user-defined */
        TYPE = 3;


        if ( TYPE == 1 ) {
            /* (quasi-)X-invariant met field: */

            LEFT = 35.0E+03;
            RIGHT= LEFT;

            DELTAT = 1.0; /* Add a ~1.0K depression far away */
#pragma omp parallel for                \
            if      ( !PARALLEL_CASES ) \
            default ( shared          ) \
            private ( iNx, jNy        ) \
            schedule( dynamic, 1      )
            for ( jNy = 0; jNy < Y.size(); jNy++ ) {
                temp_[jNy][X.size()/2] = TEMPERATURE + Y[jNy] * LAPSERATE + diurnalPert;
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

            /* RHi bubble centered on x = 0
             * Add a ~3.0K depression to ambient air */
            DELTAT = 3.0;

            TOP = 600.0;
            if ( USERINPUT.MET_FIXDEPTH )
                BOT = USERINPUT.MET_DEPTH;
            else
                BOT = 200.0;
#ifndef XLIM
            LEFT = XLIM_LEFT-0.50E+04; /* 30.00E+03; */
            RIGHT= XLIM_RIGHT-0.50E+04;
#else
            LEFT = XLIM-0.50E+04; /* 30.00E+03; */
            RIGHT= XLIM-0.50E+04;
#endif

#pragma omp parallel for                \
            if      ( !PARALLEL_CASES ) \
            default ( shared          ) \
            private ( iNx, jNy        ) \
            schedule( dynamic, 1      )
            for ( jNy = 0; jNy < Y.size(); jNy++ ) {
                H2O_[jNy][0] = 0.0E+00;
                for ( iNx = 0; iNx < X.size(); iNx++ ) {
                    temp_[jNy][iNx] = TEMPERATURE + DELTAT + diurnalPert + LAPSERATE * Y[jNy] \
                                   - ( DELTAT ) \
                                   * ( 1.0 - 0.5 * ( std::tanh( ( X[iNx] - LEFT ) / 1.0E+03 ) + 1.0 )) \
                                   * ( 0.0 + 0.5 * ( std::tanh( ( X[iNx] + RIGHT) / 1.0E+03 ) + 1.0 ));
    //                               * ( 1.0 - 0.5 * ( std::tanh( ( Y[jNy] - TOP  ) / 1.0E+02 ) + 1.0 )) \
    //                               * ( 0.0 + 0.5 * ( std::tanh( ( Y[jNy] + BOT  ) / 1.0E+02 ) + 1.0 ))
                    H2O_[jNy][iNx]  = H2O_[jNy][0];
                }
            }

        } else if ( TYPE == 3 ) {

            /* RHi layer centered on y = 0 */
            TOP = 200.0;
            if ( USERINPUT.MET_FIXDEPTH )
                BOT = USERINPUT.MET_DEPTH;
            else
                BOT = 200.0;

            //RH_star = 110;
            RH_star = RHI;
            RH_far  = 50;

#pragma omp parallel for                \
            if      ( !PARALLEL_CASES ) \
            default ( shared          ) \
            private ( iNx, jNy        ) \
            schedule( dynamic, 1      )
            for ( jNy = 0; jNy < Y.size(); jNy++ ) {
                temp_[jNy][0] = TEMPERATURE + Y[jNy] * LAPSERATE + diurnalPert;
                /* Unit check: Y [m] * LapseRate [K/m] */

                if ( ( Y[jNy] < TOP ) && ( Y[jNy] > -BOT ) ) {
                    RH = RH_star;
                    H2O_[jNy][0] = RH / 1.0E+02 * physFunc::pSat_H2Os( temp_[jNy][0] ) / ( physConst::kB * temp_[jNy][0] * 1.0E+06 );
                } else {
                    if ( Y[jNy] < 0 ) {
                        RH = RH_star - ( Y[jNy] + BOT ) / BOT * ( RH_far - RH_star );
                    } else {
                        RH = RH_star + ( Y[jNy] - TOP ) / TOP * ( RH_far - RH_star );
                    }
                    if ( ( Y[jNy] < -2.0*BOT ) || ( Y[jNy] > 2.0*TOP ) )
                        RH = RH_far;
                    H2O_[jNy][0] = RH / 1.0E+02 * physFunc::pSat_H2Os( temp_[jNy][0] ) / ( physConst::kB * temp_[jNy][0] * 1.0E+06 );
                }

                for ( iNx = 0; iNx < X.size(); iNx++ ) {
                    temp_[jNy][iNx] = temp_[jNy][0];
                    H2O_[jNy][iNx]  = H2O_[jNy][0];
                }
            }

        } else {

            std::cout << " (1) In Meteorology::Meteorology: Undefined value for user-input type: ";
            std::cout << TYPE << std::endl;
            exit(-1);

        }
    }

    RealDouble invkB = 1.00E-06 / physConst::kB;
#pragma omp parallel for        \
    if      ( !PARALLEL_CASES ) \
    default ( shared          ) \
    private ( iNx, jNy        ) \
    schedule( dynamic, 1      )
    for ( jNy = 0; jNy < Y.size(); jNy++ ) {
        for ( iNx = 0; iNx < X.size(); iNx++ )
            airDens_[jNy][iNx] = press_[jNy] / temp_[jNy][iNx] * invkB;
    }


    if ( DBG ) {
        std::cout << "\n DEBUG : Meteorology\n";
        std::cout << "         Grid number | Altitude [km] | Pressure [hPa] | Temperature [K] | H2O [-]\n";
        for ( jNy = Y.size() - 1; jNy --> 0; ) {
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
    SHEAR( met.SHEAR ),
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
    airDens_      = met.airDens_;
    H2O_          = met.H2O_;
    temperature_store_ = met.temperature_store_;
    altitude_store_ = met.altitude_store_;

} /* End of Meteorology::Meteorology */

Meteorology::~Meteorology( )
{

    /* Destructor */

} /* End of Meteorology::~Meteorology */

void Meteorology::Update( const OptInput &USERINPUT, const RealDouble solarTime_h, \
                          const RealDouble simTime_h, const Mesh &m, \
                          const RealDouble dTrav_x, const RealDouble dTrav_y )
{

    UInt iNx = 0;
    UInt jNy = 0;

    Vector_1D X = m.x();
    Vector_1D Y = m.y();

    diurnalPert = DIURNAL_AMPL * cos( 2.0E+00 * physConst::PI * ( solarTime_h - DIURNAL_PHASE ) / 24.0E+00 );

    
    /* User defined fields can be set here ! */
    if ( USERINPUT.MET_LOADMET ) {

        std::cout << "Umm... you need to set this up, mate" << std::endl;

        /* Define variable sizes */
	UInt var_len = temperature_store_.size();
	float temperature_before[var_len];
	float temperature_after[var_len];
	float temperature_interp[var_len];

	/* Extract temperature data before and after current time, and interpolate */
        UInt itime_extract = std::floor( simTime_h / 3 );
        for ( UInt i = 0; i < var_len; i++ ) {
	    temperature_before[i] = temperature_store_[i][itime_extract];
	    temperature_after[i] = temperature_store_[i][itime_extract+1];
	    temperature_interp[i] = temperature_before[i] + ( temperature_after[i] - temperature_before[i] ) / 3.0 * ( simTime_h - itime_extract * 3.0 );
	}

        /* Identify closest temperature to given pressure */
        /* Loop round each vertical layer to estimate temperature */
        for ( UInt jNy = 0; jNy < Y.size(); jNy++ ) {
            
	    /* Find the closest values above and below the central pressure */
            UInt i_Z = met::nearestNeighbor( altitude_store_, alt_[jNy] );
            
	    /* Loop round horizontal coordinates to assign temperature */
            for ( UInt iNx = 0; iNx < X.size(); iNx++ ) {
                temp_[jNy][iNx] = temperature_interp[i_Z];
            }

        }

    }
    else {

#pragma omp parallel for        \
        if      ( !PARALLEL_CASES ) \
        default ( shared          ) \
        private ( jNy             ) \
        schedule( dynamic, 1      )
        for ( jNy = 0; jNy < Y.size(); jNy++ )
            alt_[jNy] = ALTITUDE + Y[jNy] + dTrav_y;
        met::ISA( alt_, press_ );
	
	if ( TYPE == 1 ) {
	    
	    /* (quasi-)X-invariant met field: */

#pragma omp parallel for            \
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

#pragma omp parallel for            \
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

	} else if ( TYPE == 3 ) {

    //        /* RHi layer centered on y = 0 */

    //        RH_star = RHI;
    //        RH_far  = 50;

    //#pragma omp parallel for            \
    //        if      ( !PARALLEL_CASES ) \
    //        default ( shared          ) \
    //        private ( iNx, jNy        ) \
    //        schedule( dynamic, 1      )
    //        for ( jNy = 0; jNy < Y.size(); jNy++ ) {
    //            temp_[jNy][0] = TEMPERATURE + Y[jNy] * LAPSERATE + diurnalPert;
    //            /* Unit check: Y [m] * LapseRate [K/m] */

    //            if ( ( Y[jNy] < TOP ) && ( Y[jNy] > -BOT ) ) {
    //                RH = RH_star;
    //                H2O_[jNy][0] = RH / 1.0E+02 * physFunc::pSat_H2Os( temp_[jNy][0] ) / ( physConst::kB * temp_[jNy][0] * 1.0E+06 );
    //            } else {
    //                if ( Y[jNy] < 0 ) {
    //                    RH = RH_star - ( Y[jNy] + BOT ) / BOT * ( RH_far - RH_star );
    //                } else {
    //                    RH = RH_star + ( Y[jNy] - TOP ) / TOP * ( RH_far - RH_star );
    //                }
    //                if ( ( Y[jNy] < -2.0*BOT ) || ( Y[jNy] > 2.0*TOP ) )
    //                    RH = RH_far;
    //                H2O_[jNy][0] = RH / 1.0E+02 * physFunc::pSat_H2Os( temp_[jNy][0] ) / ( physConst::kB * temp_[jNy][0] * 1.0E+06 );
    //            }

    //            for ( iNx = 0; iNx < X.size(); iNx++ ) {
    //                temp_[jNy][iNx] = temp_[jNy][0];
    //                H2O_[jNy][iNx]  = H2O_[jNy][0];
    //            }
    //        }


	} else {

	    std::cout << " In Meteorology::Meteorology: Undefined value for user-input type: ";
	    std::cout << TYPE << std::endl;
	    exit(-1);

	}
    }

RealDouble invkB = 1.00E-06 / physConst::kB;
#pragma omp parallel for        \
    if      ( !PARALLEL_CASES ) \
    default ( shared          ) \
    private ( iNx, jNy        ) \
    schedule( dynamic, 1      )
    for ( UInt jNy = 0; jNy < Y.size(); jNy++ ) {
        for ( UInt iNx = 0; iNx < X.size(); iNx++ )
            airDens_[jNy][iNx] = press_[jNy] / temp_[jNy][iNx] * invkB;
    }

    if ( 1 ) {
        std::cout << "\n DEBUG : Meteorology\n";
        std::cout << "         Grid number | Altitude [km] | Pressure [hPa] | Temperature [K] | H2O [-] | RHi [-] |\n";
        for ( jNy = Y.size() - 1; jNy --> 0; ) {
            std::cout << "         " << std::setw(11) << jNy << " | " << std::setw(13) << alt_[jNy] * 1.00E-03 \
                      << " | " << std::setw(14) << press_[jNy] * 1.00E-02 << " | " << std::setw(15) << temp_[jNy][0] \
                      << " | " << std::setw(9) << H2O_[jNy][0] << " | " << std::setw(9) << H2O_[jNy][0]*physConst::kB*temp_[jNy][0]/physFunc::pSat_H2Os( temp_[jNy][0] )*1.0E+06 << "\n";
        }
        std::cout << "\n";
    }

} /* End of Meteorology::UpdateMet */

/* End of Meteorology.cpp */
