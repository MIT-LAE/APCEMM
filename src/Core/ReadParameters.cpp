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
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <iostream>
#include <vector>

std::vector<std::vector<double> > Copy_blocked( std::vector<std::vector<double> >& m, \
                                                int n );
std::vector<std::vector<double> > Copy_interleaved( std::vector<std::vector<double> >& m, \
                                                    int n );
std::vector<std::vector<double> > Reshape_Vector( std::vector<std::vector<double> >& vector_2D, int n_x, int n_y );

std::vector<std::vector<double> > CombVec( const std::vector<double>& temperature_K, \
                                           const std::vector<double>& pressure_Pa,   \
                                           const std::vector<double>& relHumidity_w, \
                                           const std::vector<double>& longitude_deg, \
                                           const std::vector<double>& latitude_deg,  \
                                           const std::vector<double>& dayGMT,        \
                                           const std::vector<double>& emissionTime,  \
                                           const std::vector<double>& EI_NOx,        \
                                           const std::vector<double>& EI_CO,         \
                                           const std::vector<double>& EI_HC,         \
                                           const std::vector<double>& EI_SO2,        \
                                           const std::vector<double>& EI_Soot,       \
                                           const std::vector<double>& SootRad,       \
                                           const std::vector<double>& ff,            \
                                           const std::vector<double>& backgNOx,      \
                                           const std::vector<double>& backgHNO3,     \
                                           const std::vector<double>& backgO3,       \
                                           const std::vector<double>& backgCO,       \
                                           const std::vector<double>& backgCH4,      \
                                           const std::vector<double>& backgSO2 );

std::vector<std::vector<double> > ReadParameters( )
{

    std::vector<double> temperature_K, pressure_Pa, relHumidity_w, longitude_deg, latitude_deg;
    std::vector<double> EI_NOx, EI_CO, EI_HC, EI_SO2, EI_Soot, SootRad, ff;
    std::vector<double> dayGMT;
    std::vector<double> emissionTime;
    std::vector<double> backgNOx, backgHNO3, backgO3, backgCO, backgCH4, backgSO2;

    std::vector<std::vector<double> > parameters;

    /* #######################################
     * ## Ambient meteorological conditions ## 
     * ####################################### */

    /* Temperature array in [K] */

    temperature_K.push_back( 220.0 );
//    for ( unsigned int i = 0; i < 31; i++ ) {
//        temperature_K.push_back( 200.0 + 2.0 * i );
//    }

    /* Pressure array in [Pa] */

    pressure_Pa.push_back( 220.0E2 );

    /* Relative humidity w.r.t liquid water array in [\%] */

    relHumidity_w.push_back( 30.0 );

    /* #######################
     * ## Coord. (lon-lat)  ## 
     * ####################### */

    /* Longitude array expressed in deg */

    longitude_deg.push_back( -15.0 );

    /* Latitude array expressed in deg */

    latitude_deg.push_back( 60.0 );

    /* #######################
     * ## Emiss day (0-365) ## 
     * ####################### */

    dayGMT.push_back( 81.0 );

    /* #######################
     * ## Emiss time  [hrs] ## 
     * ####################### */

    emissionTime.push_back( 8.0 );
    emissionTime.push_back( 20.0 );

    /* ####################### 
     * ## EMISSION INDICES  ##
     * ####################### */

    EI_NOx.push_back( 1.0E-05 );
    EI_NOx.push_back( 1.0E+01 );

//    EI_CO.push_back( 0.5E+00 );
    EI_CO.push_back( 1.0E+00 );
//    EI_CO.push_back( 1.5E+00 );
//    EI_CO.push_back( 2.0E+00 );

    EI_HC.push_back( 0.6E+00 );
    
    EI_SO2.push_back( 0.8E+00 );
    EI_SO2.push_back( 1.2E+00 );
    EI_SO2.push_back( 2.0E+00 );
    EI_SO2.push_back( 1.0E+01 );
    EI_SO2.push_back( 2.0E+01 );
    EI_SO2.push_back( 1.0E+02 );
    
    EI_Soot.push_back( 0.0E+00 );
    
    SootRad.push_back( 0.0E+00 );
    
    ff.push_back( 0.0E+00 );

    /* #######################
     * ## BACKGRD MIX RATIO ##
     * ####################### */

    backgNOx.push_back( 50.0E-03 );
//    backgNOx.push_back( 100.0E-03 );
//    backgNOx.push_back( 150.0E-03 );
//    backgNOx.push_back( 200.0E-03 );
    
    backgHNO3.push_back( 81.5E-03 );

    backgO3.push_back( 50.0E+00 );
//    backgO3.push_back( 75.0E+00 );
//    backgO3.push_back( 100.0E+00 );

//    backgCO.push_back( 20.0E+00 );
    backgCO.push_back( 40.0E+00 );
//    backgCO.push_back( 60.0E+00 );
//    backgCO.push_back( 80.0E+00 );
    
    backgCH4.push_back( 1.76E+03 );
    backgSO2.push_back( 7.25E-03 );

    parameters = CombVec( temperature_K, \
                          pressure_Pa,   \
                          relHumidity_w, \
                          longitude_deg, \
                          latitude_deg,  \
                          dayGMT,        \
                          emissionTime,  \
                          EI_NOx,        \
                          EI_CO,         \
                          EI_HC,         \
                          EI_SO2,        \
                          EI_Soot,       \
                          SootRad,       \
                          ff,            \
                          backgNOx,      \
                          backgHNO3,     \
                          backgO3,       \
                          backgCO,       \
                          backgCH4,      \
                          backgSO2 );

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

std::vector<std::vector<double> > CombVec( const std::vector<double>& temperature_K, \
                                           const std::vector<double>& pressure_Pa,   \
                                           const std::vector<double>& relHumidity_w, \
                                           const std::vector<double>& longitude_deg, \
                                           const std::vector<double>& latitude_deg,  \
                                           const std::vector<double>& dayGMT,        \
                                           const std::vector<double>& emissionTime,  \
                                           const std::vector<double>& EI_NOx,        \
                                           const std::vector<double>& EI_CO,         \
                                           const std::vector<double>& EI_HC,         \
                                           const std::vector<double>& EI_SO2,         \
                                           const std::vector<double>& EI_Soot,       \
                                           const std::vector<double>& SootRad,       \
                                           const std::vector<double>& ff,            \
                                           const std::vector<double>& backgNOx,      \
                                           const std::vector<double>& backgHNO3,     \
                                           const std::vector<double>& backgO3,       \
                                           const std::vector<double>& backgCO,       \
                                           const std::vector<double>& backgCH4,      \
                                           const std::vector<double>& backgSO2 )
{
    std::vector<std::vector<double> > combinations;

    unsigned int counter = 1;
    unsigned int nCases = temperature_K.size();

    unsigned int i, j;

    std::vector<std::vector<double> > y, z;
    std::vector<std::vector<double> > u, v;

    y.push_back(std::vector<double>(temperature_K.size()));
    for ( i = 0; i < temperature_K.size(); i++ )
        y[0][i] = temperature_K[i];

    /* z = Pressure_Pa */
    z.push_back(std::vector<double>(pressure_Pa.size()));
    for ( i = 0; i < pressure_Pa.size(); i++ )
        z[0][i] = pressure_Pa[i];

    u = Copy_blocked(y,z[0].size());
    v = Copy_interleaved(z,y[0].size());

    for ( i = 0; i < counter; i++ ) {
        y[i].clear();
    }
    z[0].clear();
    y.clear(); z.clear();

    nCases *= pressure_Pa.size();
    counter += 1;
    for ( i = 0; i < counter; i++ )
        y.push_back(std::vector<double>( nCases ));

    for ( i = 0; i < nCases; i++ ) {
        for ( j = 0; j < counter - 1; j++ ) {
            y[j][i] = u[j][i];
        }
        y[counter-1][i] = v[0][i];
    }

    for ( i = 0; i < counter - 1; i++ ) {
        u[i].clear();
    }
    u.clear();
    v[0].clear(); v.clear();

    /* z = relHumidity_w */
    z.push_back(std::vector<double>(relHumidity_w.size()));
    for ( i = 0; i < relHumidity_w.size(); i++ )
        z[0][i] = relHumidity_w[i];

    u = Copy_blocked(y,z[0].size());
    v = Copy_interleaved(z,y[0].size());

    for ( i = 0; i < counter; i++ ) {
        y[i].clear();
    }
    z[0].clear();
    y.clear(); z.clear();

    nCases *= relHumidity_w.size();
    counter += 1;
    for ( i = 0; i < counter; i++ )
        y.push_back(std::vector<double>( nCases ));

    for ( i = 0; i < nCases; i++ ) {
        for ( j = 0; j < counter - 1; j++ ) {
            y[j][i] = u[j][i];
        }
        y[counter-1][i] = v[0][i];
    }

    for ( i = 0; i < counter - 1; i++ ) {
        u[i].clear();
    }
    u.clear();
    v[0].clear(); v.clear();

    /* z = longitude_deg */
    z.push_back(std::vector<double>(longitude_deg.size()));
    for ( i = 0; i < longitude_deg.size(); i++ )
        z[0][i] = longitude_deg[i];

    u = Copy_blocked(y,z[0].size());
    v = Copy_interleaved(z,y[0].size());

    for ( i = 0; i < counter; i++ ) {
        y[i].clear();
    }
    z[0].clear();
    y.clear(); z.clear();

    nCases *= longitude_deg.size();
    counter += 1;
    for ( i = 0; i < counter; i++ )
        y.push_back(std::vector<double>( nCases ));

    for ( i = 0; i < nCases; i++ ) {
        for ( j = 0; j < counter - 1; j++ ) {
            y[j][i] = u[j][i];
        }
        y[counter-1][i] = v[0][i];
    }

    for ( i = 0; i < counter - 1; i++ ) {
        u[i].clear();
    }
    u.clear();
    v[0].clear(); v.clear();

    /* z = latitude_deg */
    z.push_back(std::vector<double>(latitude_deg.size()));
    for ( i = 0; i < latitude_deg.size(); i++ )
        z[0][i] = latitude_deg[i];

    u = Copy_blocked(y,z[0].size());
    v = Copy_interleaved(z,y[0].size());

    for ( i = 0; i < counter; i++ ) {
        y[i].clear();
    }
    z[0].clear();
    y.clear(); z.clear();

    nCases *= latitude_deg.size();
    counter += 1;
    for ( i = 0; i < counter; i++ )
        y.push_back(std::vector<double>( nCases ));

    for ( i = 0; i < nCases; i++ ) {
        for ( j = 0; j < counter - 1; j++ ) {
            y[j][i] = u[j][i];
        }
        y[counter-1][i] = v[0][i];
    }
    
    /* z = dayGMT */
    z.push_back(std::vector<double>(dayGMT.size()));
    for ( i = 0; i < dayGMT.size(); i++ )
        z[0][i] = dayGMT[i];

    u = Copy_blocked(y,z[0].size());
    v = Copy_interleaved(z,y[0].size());

    for ( i = 0; i < counter; i++ ) {
        y[i].clear();
    }
    z[0].clear();
    y.clear(); z.clear();

    nCases *= dayGMT.size();
    counter += 1;
    for ( i = 0; i < counter; i++ )
        y.push_back(std::vector<double>( nCases ));

    for ( i = 0; i < nCases; i++ ) {
        for ( j = 0; j < counter - 1; j++ ) {
            y[j][i] = u[j][i];
        }
        y[counter-1][i] = v[0][i];
    }
    
    /* z = emissionTime */
    z.push_back(std::vector<double>(emissionTime.size()));
    for ( i = 0; i < emissionTime.size(); i++ )
        z[0][i] = emissionTime[i];

    u = Copy_blocked(y,z[0].size());
    v = Copy_interleaved(z,y[0].size());

    for ( i = 0; i < counter; i++ ) {
        y[i].clear();
    }
    z[0].clear();
    y.clear(); z.clear();

    nCases *= emissionTime.size();
    counter += 1;
    for ( i = 0; i < counter; i++ )
        y.push_back(std::vector<double>( nCases ));

    for ( i = 0; i < nCases; i++ ) {
        for ( j = 0; j < counter - 1; j++ ) {
            y[j][i] = u[j][i];
        }
        y[counter-1][i] = v[0][i];
    }

    /* z = EI_NOx */
    z.push_back(std::vector<double>(EI_NOx.size()));
    for ( i = 0; i < EI_NOx.size(); i++ )
        z[0][i] = EI_NOx[i];

    u = Copy_blocked(y,z[0].size());
    v = Copy_interleaved(z,y[0].size());

    for ( i = 0; i < counter; i++ ) {
        y[i].clear();
    }
    z[0].clear();
    y.clear(); z.clear();

    nCases *= EI_NOx.size();
    counter += 1;
    for ( i = 0; i < counter; i++ )
        y.push_back(std::vector<double>( nCases ));

    for ( i = 0; i < nCases; i++ ) {
        for ( j = 0; j < counter - 1; j++ ) {
            y[j][i] = u[j][i];
        }
        y[counter-1][i] = v[0][i];
    }

    /* z = EI_CO */
    z.push_back(std::vector<double>(EI_CO.size()));
    for ( i = 0; i < EI_CO.size(); i++ )
        z[0][i] = EI_CO[i];

    u = Copy_blocked(y,z[0].size());
    v = Copy_interleaved(z,y[0].size());

    for ( i = 0; i < counter; i++ ) {
        y[i].clear();
    }
    z[0].clear();
    y.clear(); z.clear();

    nCases *= EI_CO.size();
    counter += 1;
    for ( i = 0; i < counter; i++ )
        y.push_back(std::vector<double>( nCases ));

    for ( i = 0; i < nCases; i++ ) {
        for ( j = 0; j < counter - 1; j++ ) {
            y[j][i] = u[j][i];
        }
        y[counter-1][i] = v[0][i];
    }

    /* z = EI_HC */
    z.push_back(std::vector<double>(EI_HC.size()));
    for ( i = 0; i < EI_HC.size(); i++ )
        z[0][i] = EI_HC[i];

    u = Copy_blocked(y,z[0].size());
    v = Copy_interleaved(z,y[0].size());

    for ( i = 0; i < counter; i++ ) {
        y[i].clear();
    }
    z[0].clear();
    y.clear(); z.clear();

    nCases *= EI_HC.size();
    counter += 1;
    for ( i = 0; i < counter; i++ )
        y.push_back(std::vector<double>( nCases ));

    for ( i = 0; i < nCases; i++ ) {
        for ( j = 0; j < counter - 1; j++ ) {
            y[j][i] = u[j][i];
        }
        y[counter-1][i] = v[0][i];
    }

    /* z = EI_SO2 */
    z.push_back(std::vector<double>(EI_SO2.size()));
    for ( i = 0; i < EI_SO2.size(); i++ )
        z[0][i] = EI_SO2[i];

    u = Copy_blocked(y,z[0].size());
    v = Copy_interleaved(z,y[0].size());

    for ( i = 0; i < counter; i++ ) {
        y[i].clear();
    }
    z[0].clear();
    y.clear(); z.clear();

    nCases *= EI_SO2.size();
    counter += 1;
    for ( i = 0; i < counter; i++ )
        y.push_back(std::vector<double>( nCases ));

    for ( i = 0; i < nCases; i++ ) {
        for ( j = 0; j < counter - 1; j++ ) {
            y[j][i] = u[j][i];
        }
        y[counter-1][i] = v[0][i];
    }

    /* z = EI_Soot */
    z.push_back(std::vector<double>(EI_Soot.size()));
    for ( i = 0; i < EI_Soot.size(); i++ )
        z[0][i] = EI_Soot[i];

    u = Copy_blocked(y,z[0].size());
    v = Copy_interleaved(z,y[0].size());

    for ( i = 0; i < counter; i++ ) {
        y[i].clear();
    }
    z[0].clear();
    y.clear(); z.clear();

    nCases *= EI_Soot.size();
    counter += 1;
    for ( i = 0; i < counter; i++ )
        y.push_back(std::vector<double>( nCases ));

    for ( i = 0; i < nCases; i++ ) {
        for ( j = 0; j < counter - 1; j++ ) {
            y[j][i] = u[j][i];
        }
        y[counter-1][i] = v[0][i];
    }

    /* z = SootRad */
    z.push_back(std::vector<double>(SootRad.size()));
    for ( i = 0; i < SootRad.size(); i++ )
        z[0][i] = SootRad[i];

    u = Copy_blocked(y,z[0].size());
    v = Copy_interleaved(z,y[0].size());

    for ( i = 0; i < counter; i++ ) {
        y[i].clear();
    }
    z[0].clear();
    y.clear(); z.clear();

    nCases *= SootRad.size();
    counter += 1;
    for ( i = 0; i < counter; i++ )
        y.push_back(std::vector<double>( nCases ));

    for ( i = 0; i < nCases; i++ ) {
        for ( j = 0; j < counter - 1; j++ ) {
            y[j][i] = u[j][i];
        }
        y[counter-1][i] = v[0][i];
    }

    /* z = ff */
    z.push_back(std::vector<double>(ff.size()));
    for ( i = 0; i < ff.size(); i++ )
        z[0][i] = ff[i];

    u = Copy_blocked(y,z[0].size());
    v = Copy_interleaved(z,y[0].size());

    for ( i = 0; i < counter; i++ ) {
        y[i].clear();
    }
    z[0].clear();
    y.clear(); z.clear();

    nCases *= ff.size();
    counter += 1;
    for ( i = 0; i < counter; i++ )
        y.push_back(std::vector<double>( nCases ));

    for ( i = 0; i < nCases; i++ ) {
        for ( j = 0; j < counter - 1; j++ ) {
            y[j][i] = u[j][i];
        }
        y[counter-1][i] = v[0][i];
    }

    /* z = backgNOx */
    z.push_back(std::vector<double>(backgNOx.size()));
    for ( i = 0; i < backgNOx.size(); i++ )
        z[0][i] = backgNOx[i];

    u = Copy_blocked(y,z[0].size());
    v = Copy_interleaved(z,y[0].size());

    for ( i = 0; i < counter; i++ ) {
        y[i].clear();
    }
    z[0].clear();
    y.clear(); z.clear();

    nCases *= backgNOx.size();
    counter += 1;
    for ( i = 0; i < counter; i++ )
        y.push_back(std::vector<double>( nCases ));

    for ( i = 0; i < nCases; i++ ) {
        for ( j = 0; j < counter - 1; j++ ) {
            y[j][i] = u[j][i];
        }
        y[counter-1][i] = v[0][i];
    }
    
    /* z = backgHNO3 */
    z.push_back(std::vector<double>(backgHNO3.size()));
    for ( i = 0; i < backgHNO3.size(); i++ )
        z[0][i] = backgHNO3[i];

    u = Copy_blocked(y,z[0].size());
    v = Copy_interleaved(z,y[0].size());

    for ( i = 0; i < counter; i++ ) {
        y[i].clear();
    }
    z[0].clear();
    y.clear(); z.clear();

    nCases *= backgHNO3.size();
    counter += 1;
    for ( i = 0; i < counter; i++ )
        y.push_back(std::vector<double>( nCases ));

    for ( i = 0; i < nCases; i++ ) {
        for ( j = 0; j < counter - 1; j++ ) {
            y[j][i] = u[j][i];
        }
        y[counter-1][i] = v[0][i];
    }
    
    /* z = backgO3 */
    z.push_back(std::vector<double>(backgO3.size()));
    for ( i = 0; i < backgO3.size(); i++ )
        z[0][i] = backgO3[i];

    u = Copy_blocked(y,z[0].size());
    v = Copy_interleaved(z,y[0].size());

    for ( i = 0; i < counter; i++ ) {
        y[i].clear();
    }
    z[0].clear();
    y.clear(); z.clear();

    nCases *= backgO3.size();
    counter += 1;
    for ( i = 0; i < counter; i++ )
        y.push_back(std::vector<double>( nCases ));

    for ( i = 0; i < nCases; i++ ) {
        for ( j = 0; j < counter - 1; j++ ) {
            y[j][i] = u[j][i];
        }
        y[counter-1][i] = v[0][i];
    }
    
    /* z = backgCO */
    z.push_back(std::vector<double>(backgCO.size()));
    for ( i = 0; i < backgCO.size(); i++ )
        z[0][i] = backgCO[i];

    u = Copy_blocked(y,z[0].size());
    v = Copy_interleaved(z,y[0].size());

    for ( i = 0; i < counter; i++ ) {
        y[i].clear();
    }
    z[0].clear();
    y.clear(); z.clear();

    nCases *= backgCO.size();
    counter += 1;
    for ( i = 0; i < counter; i++ )
        y.push_back(std::vector<double>( nCases ));

    for ( i = 0; i < nCases; i++ ) {
        for ( j = 0; j < counter - 1; j++ ) {
            y[j][i] = u[j][i];
        }
        y[counter-1][i] = v[0][i];
    }
    
    /* z = backgCH4 */
    z.push_back(std::vector<double>(backgCH4.size()));
    for ( i = 0; i < backgCH4.size(); i++ )
        z[0][i] = backgCH4[i];

    u = Copy_blocked(y,z[0].size());
    v = Copy_interleaved(z,y[0].size());

    for ( i = 0; i < counter; i++ ) {
        y[i].clear();
    }
    z[0].clear();
    y.clear(); z.clear();

    nCases *= backgCH4.size();
    counter += 1;
    for ( i = 0; i < counter; i++ )
        y.push_back(std::vector<double>( nCases ));

    for ( i = 0; i < nCases; i++ ) {
        for ( j = 0; j < counter - 1; j++ ) {
            y[j][i] = u[j][i];
        }
        y[counter-1][i] = v[0][i];
    }
    
    /* z = backgSO2 */
    z.push_back(std::vector<double>(backgSO2.size()));
    for ( i = 0; i < backgSO2.size(); i++ )
        z[0][i] = backgSO2[i];

    u = Copy_blocked(y,z[0].size());
    v = Copy_interleaved(z,y[0].size());

    for ( i = 0; i < counter; i++ ) {
        y[i].clear();
    }
    z[0].clear();
    y.clear(); z.clear();

    nCases *= backgSO2.size();
    counter += 1;
    for ( i = 0; i < counter; i++ )
        y.push_back(std::vector<double>( nCases ));

    for ( i = 0; i < nCases; i++ ) {
        for ( j = 0; j < counter - 1; j++ ) {
            y[j][i] = u[j][i];
        }
        y[counter-1][i] = v[0][i];
    }
    
    return y;

} /* End of CombVec */

std::vector<std::vector<double> > Copy_blocked( std::vector<std::vector<double> >& m, \
                                                int n )
{
    std::vector<std::vector<double> > b;
    const unsigned int mr = m.size();
    const unsigned int mc = m[0].size();
    unsigned int i, j, k;

    for ( i = 0; i < mr; i++ )
        b.push_back(std::vector<double>(mc*n));

//    unsigned int ind[mc];
    unsigned int *ind = new unsigned int[mc];

    for ( i = 0; i < mc; i++ )
        ind[i] = i;

    for ( j = 0; j < (n-1) * mc + 1; j+=mc ) {
        for ( i = 0; i < mr; i++ ) {
            for ( k = 0; k < mc; k++ ) {
                b[i][ind[k] + j] = m[i][k];
            }
        }
    }

    delete[] ind;

    return b;
} /* End of Copy_blocked */

std::vector<std::vector<double> > Copy_interleaved( std::vector<std::vector<double> >& m, \
                                                    int n )
{
    std::vector<std::vector<double> > b;
    const unsigned int mr = m.size();
    const unsigned int mc = m[0].size();
    unsigned int i, j, k;

    for ( i = 0; i < mr*n; i++ )
        b.push_back(std::vector<double>(mc));

//    unsigned int ind[mr];
    unsigned int *ind = new unsigned int[mr];

    for ( i = 0; i < mr; i++ )
        ind[i] = i;

    for ( i = 0; i < (n-1) * mr + 1; i+=mr ) {
        for ( j = 0; j < mc; j++ ) {
            for ( k = 0; k < mr; k++ ) {
                b[ind[k] + i][j] = m[k][j]; //m[k][j]
            }
        }
    }

    delete[] ind;

    return Reshape_Vector( b, mr, n*mc );

} /* End of Copy_interleaved */

std::vector<std::vector<double> > Reshape_Vector( std::vector<std::vector<double> >& vector_2D, int n_x, int n_y )
{
    std::vector<std::vector<double> > output;
    int size_x = vector_2D.size();
    int size_y = vector_2D[0].size();

    if ( n_x * n_y != size_x * size_y ) {
        std::cout << "Invalid dimensions specified" << std::endl;
        return output;
    }
    else {
        int counter;
        for ( counter = 0; counter < n_x; counter++ )
            output.push_back(std::vector<double>(n_y));

        int counter_x = 0;
        int counter_y = 0;
        int orig_counter_x = 0;
        int orig_counter_y = 0;
        for ( counter = 0; counter < n_x * n_y; counter++ ) {
            counter_x = counter%n_x;
            //counter_y = counter%n_y;
            orig_counter_x = counter%size_x;
            //orig_counter_y = counter%size_y;
            output[counter_x][counter_y] = vector_2D[orig_counter_x][orig_counter_y];
            //std::cout << orig_counter_x << ", " << orig_counter_y << ": " << vector_2D[orig_counter_x][orig_counter_y] << " -> " << counter_x << ", " << counter_y << std::endl;
            if ( counter_x == n_x - 1 )
                counter_y += 1;
            if ( orig_counter_x == size_x - 1 )
                orig_counter_y += 1;
//            if ( counter_y == n_y - 1 )
//                counter_x += 1;
//            if ( orig_counter_y == size_y - 1 )
//                orig_counter_x += 1;
        }
//        int k, l;
//        for ( i = 0; i < n_x; i++ ) {
//            for ( j = 0; j < n_y; j++ ) {
//                k = i%size_x;
//                l = j%size_y;
//                output[i][j] = vector_2D[k][l];
//            }
//        }
    }

    return output;

} /* End of Reshape_Vector */

/* End of ReadParameters.cpp */
