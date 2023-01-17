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
        loadMet(USERINPUT, solarTime_h, m, temperature_K, pressure_Pa, relHumidity_i, shear_pm, 0);

    } else { 

        met::ISA_pAlt( ALTITUDE, PRESSURE );

        for ( jNy = 0; jNy < Y.size(); jNy++ )
            alt_[jNy] = ALTITUDE + Y[jNy];
        met::ISA( alt_, press_ );

        /* Assign shear vector */
        for ( jNy = 0; jNy < Y.size(); jNy++ )
            shear_[jNy] = SHEAR;

        /* Set temperature input type if user-defined */
        //Not sure what these types are trying to do -Michael
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
    shear_store_ = met.shear_store_;

} /* End of Meteorology::Meteorology */

Meteorology::~Meteorology( )
{

    /* Destructor */

} /* End of Meteorology::~Meteorology */

void Meteorology::loadMet(const OptInput &USERINPUT,      \
                          const RealDouble solarTime_h,   \
                          const Mesh &m,                  \
                          const RealDouble temperature_K, \
                          const RealDouble pressure_Pa,   \
                          const RealDouble relHumidity_i, \
                          const RealDouble shear_pm, \
                          const bool DBG )
{
    Vector_1D X = m.x();
    Vector_1D Y = m.y();
    //TODO: Rewrite to work with native ERA5/GEOS5 output
    //      Make this work with time-dependent met. data
    //      Make emission hour and temperature/shear met updating consistent...

    /*
        Met data format:
        
        2 dimensions: altitude and time
        Temp. and shear are functions of both

        RH is for some reason independent of time
    */

    /* Met is loaded in */
    /* TYPE = 0; */
    TYPE = 3; /* TYPE = 0 means met for every time step */

    /* Open the netcdf file for read access */
    NcFile dataFile( USERINPUT.MET_FILENAME.c_str(), NcFile::read  );

    /* Identify the length of variables in input file */
    NcDim zDim = dataFile.getDim("altitude");
    int altitudeDim = zDim.getSize();

    int timeDim = dataFile.getDim("time").getSize();

    std::cout << "altitudeDim = " << altitudeDim << std::endl;
    std::cout << "timeDim = " << timeDim << std::endl;

    /* Extract pressure and altitude from input file */
    float altitude_init[altitudeDim], pressure_init[altitudeDim];
    
    altitude_store_.assign( altitudeDim, 0.0E+00 );
    NcVar altitude_ncVar = dataFile.getVar("altitude");
    NcVar pressure_ncVar = dataFile.getVar("pressure");
    altitude_ncVar.getVar(altitude_init);
    pressure_ncVar.getVar(pressure_init);
    int i;
    for ( i = 0; i < altitudeDim; i++ ) {
        pressure_init[i] *= 100.0;
        altitude_init[i] *= 1000.0;
        altitude_store_[i] = altitude_init[i];
    }

    /* Identify closest pressure in input file from user-defined pressure */
    UInt i_Zp = met::nearestNeighbor( pressure_init, PRESSURE, altitudeDim);
    pres_user = pressure_init[i_Zp];
    alt_user = met::linearInterp( pressure_init, altitude_init, PRESSURE, altitudeDim);
    for ( UInt jNy = 0; jNy < Y.size(); jNy++ ) {
            // alt_[jNy] = altitude_user[i_Zp] + Y[jNy];
            alt_[jNy] = alt_user + Y[jNy];
        UInt i_Z = met::nearestNeighbor( altitude_init, alt_[jNy], altitudeDim);
            press_[jNy] = pressure_init[i_Z];
    }


    float relhumid_init[altitudeDim];
    float temperature_init[altitudeDim];
    float shear_init[altitudeDim];
    Vector_1D v1d_store(timeDim, 0.0E+00 );
    for ( UInt itime = 0; itime < altitudeDim; itime++ ) {
        temperature_store_.push_back( v1d_store );
        shear_store_.push_back( v1d_store );
    }

    if ( USERINPUT.MET_LOADTEMP ) {

        /* !@#$ */
        // SDE 2022-05-07: Not quite sure what this all does. My guess (?) is that
        // it is designed to read in multiple time slices, but this doesn't seem
        // like it would be able to do that.

        // MX Jan 11, 2023: The time slice was hardcoded to be 8 in length before, I fixed that. Not sure what the old format was.
        // This does NOT account for situations where the first entries in the time series do not match up with the emission time.
        // Therefore, it is up to the USER to ensure that the ordering of their met data is correct and to provide enough time series data to cover the duration of their simulation.

        NcVar temp_ncVar = dataFile.getVar("temperature");

        /* Extract the temperature data, must flatten 2d data into a 1d array for netcdf api to work*/
        float tempData[altitudeDim*timeDim];
        temp_ncVar.getVar(tempData);


        /* Identify temperature at above pressure at time 0 */
        temp_user = tempData[timeDim * i_Zp];    

        if ( USERINPUT.MET_TEMPTIMESERIES ) { 

            /* Extract 2D temperature data into a 1D array */
            for ( UInt i = 0; i < altitudeDim; i++ ) {
                //set initial user temp to temp at time = 0 for the given altitude
                temperature_init[i] = tempData[i * timeDim];

                //sets up time-dependent temp 2d array
                for ( UInt itime = 0; itime < timeDim; itime++ ) {
                    temperature_store_[i][itime] = tempData[i*timeDim + itime];
                }
            }
        }
        else {
            for ( UInt i = 0; i < altitudeDim; i++ ) {
                //set this to time = 0 at that altitude
                temperature_init[i] = tempData[i*timeDim];
            }      
            /* Identify temperature at above pressure */
            temp_user = temperature_init[i_Zp];
        }
        /* Identify closest temperature to given pressure */
        /* Loop round each vertical layer to estimate temperature */
        for ( UInt jNy = 0; jNy < Y.size(); jNy++ ) {

            /* Find the closest values above and below the central pressure */
            UInt i_Z = met::nearestNeighbor( altitude_init, alt_[jNy], altitudeDim);

            /* Loop round horizontal coordinates to assign temperature */
            for ( UInt iNx = 0; iNx < X.size(); iNx++ ) {
                temp_[jNy][iNx] = temperature_init[i_Z];
            }

        }
    }

    if ( USERINPUT.MET_LOADRH ) {
        /* !@#$ */
        //Assumes constant RHw for the entire domain

        /* Define temperature input dimension */
        NcVar relhumid_ncVar = dataFile.getVar("relative_humidity");
        relhumid_ncVar.getVar(relhumid_init);

        /* Identify temperature at above pressure */
        RHw_user = relhumid_init[i_Zp];

        /* Identify closest temperature to given pressure */
        /* Loop round each vertical layer to estimate temperature */
        for ( UInt jNy = 0; jNy < Y.size(); jNy++ ) {

            /* Find the closest values above and below the central pressure */
            UInt i_Z = met::nearestNeighbor( altitude_init, alt_[jNy], altitudeDim);
            /* Loop round horizontal coordinates to assign temperature */
            for ( UInt iNx = 0; iNx < X.size(); iNx++ ) {
                H2O_[jNy][iNx] = relhumid_init[i_Z]/((double) 100.00) *\
                                physFunc::pSat_H2Ol( temp_[jNy][iNx] ) / ( physConst::kB * temp_[jNy][iNx] ) * 1.00E-06;
            }

        }

    }

    /* Identify the saturation depth */
    satdepth_user = met::satdepth_calc( relhumid_init, temperature_init, altitude_init, i_Zp, altitudeDim );
    if ( satdepth_user != 1.0 ) {
        satdepth_user = satdepth_user - ( altitude_init[i_Zp]-alt_user );
    }
    std::cout << "Saturation Depth: " << satdepth_user << "m" << std::endl;
    if ( USERINPUT.MET_LOADSHEAR ) {
        /* !?*# */
        //Exactly same logic as loading temperature from met file
         float shearData[altitudeDim*timeDim];
        NcVar temp_ncVar = dataFile.getVar("shear");
        temp_ncVar.getVar(shearData);

        S_user = shearData[timeDim * i_Zp];
            
        if ( USERINPUT.MET_SHEARTIMESERIES ) {

            /* Extract 2D temperature data into a 1D array */
            for ( UInt i = 0; i < altitudeDim; i++ ) {
                shear_init[i] = shearData[i*timeDim];
                for ( UInt itime = 0; itime < timeDim; itime++ ) {
                    shear_store_[i][itime] = shearData[i*timeDim + itime];
                }
            }

        }
        else {
            for ( UInt i = 0; i < altitudeDim; i++ ) {
                shear_init[i] = shearData[i*timeDim];
            }
            /* Identify temperature at above pressure */
            S_user = shear_init[i_Zp];

        }

        /* Identify closest temperature to given pressure */
        /* Loop round each vertical layer to estimate temperature */
        for ( UInt jNy = 0;  jNy < Y.size(); jNy++ ) {

            /* Find the closest values above and below the central pressure */
            UInt i_Z = met::nearestNeighbor( altitude_init, alt_[jNy], altitudeDim);
            /* Loop round horizontal coordinates to assign temperature */
            for ( UInt iNx = 0; iNx < X.size(); iNx++ ) {
                shear_[jNy] = shear_init[i_Z];
            }

        }

    } else {

        /* Assign shear vector */
        for (UInt jNy = 0; jNy < Y.size(); jNy++ )
            shear_[jNy] = SHEAR;

    }                          
}
void Meteorology::Update( const OptInput &USERINPUT, const RealDouble solarTime_h, \
                          const RealDouble simTime_h, const Mesh &m, \
                          const RealDouble dTrav_x, const RealDouble dTrav_y, \
	       	          const bool DBG )
{
    
    UInt iNx = 0;
    UInt jNy = 0;

    Vector_1D X = m.x();
    Vector_1D Y = m.y();

    diurnalPert = DIURNAL_AMPL * cos( 2.0E+00 * physConst::PI * ( solarTime_h - DIURNAL_PHASE ) / 24.0E+00 );

    
    /* User defined fields can be set here ! */
    if ( USERINPUT.MET_LOADMET ) {

        if ( USERINPUT.MET_TEMPTIMESERIES ) {
	    
            /* Define variable sizes */
            UInt altitudeDim = temperature_store_.size();
            UInt timeDim = temperature_store_[0].size();

            float temperature_before[altitudeDim];
            float temperature_after[altitudeDim];
            float temperature_interp[altitudeDim];

            /* Extract temperature data before and after current time, and interpolate */
            UInt itime_extract = std::min((UInt)(simTime_h / USERINPUT.MET_DT), timeDim - 1);

            for ( UInt i = 0; i < altitudeDim; i++ ) {
                temperature_before[i] = temperature_store_[i][itime_extract];
                if ( itime_extract >= timeDim - 1 ) {
                    temperature_after[i] = temperature_store_[i][itime_extract];
                }
                else {
                    temperature_after[i] = temperature_store_[i][itime_extract+1];
                }
                temperature_interp[i] = temperature_before[i] + ( temperature_after[i] - temperature_before[i] ) * ( simTime_h - itime_extract ) / USERINPUT.MET_DT;
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

        if ( USERINPUT.MET_SHEARTIMESERIES ) {
        
            /* Define variable sizes */
            UInt altitudeDim = shear_store_.size();
            UInt timeDim = temperature_store_[0].size();
            float shear_before[altitudeDim];
            float shear_after[altitudeDim];
            float shear_interp[altitudeDim];

            /* Extract shear data before and after current time, and interpolate */
            UInt itime_extract = std::min((UInt)(simTime_h / USERINPUT.MET_DT), timeDim - 1);
            for ( UInt i = 0; i < altitudeDim; i++ ) {
                shear_before[i] = shear_store_[i][itime_extract];
                if ( itime_extract >= timeDim - 1 ) {
                    shear_after[i] = shear_store_[i][itime_extract];
                }
                else {
                    shear_after[i] = shear_store_[i][itime_extract+1];
                }
                shear_interp[i] = shear_before[i] + ( shear_after[i] - shear_before[i] ) * ( simTime_h - itime_extract ) / USERINPUT.MET_DT;
            }

            /* Identify closest shear to given pressure */
            /* Loop round each vertical layer to estimate shear */
            for ( UInt jNy = 0; jNy < Y.size(); jNy++ ) {
                /* Find the closest values above and below the central pressure */
                UInt i_Z = met::nearestNeighbor( altitude_store_, alt_[jNy] );
                shear_[jNy] = shear_interp[i_Z];
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

    if ( DBG ) {
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

