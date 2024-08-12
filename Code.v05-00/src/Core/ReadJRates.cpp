/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* ReadJRates Program File                                          */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 1/16/2019                                 */
/* File                 : ReadJRates.cpp                            */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "Core/ReadJRates.hpp"

void ReadJRates( const char* ROOTDIR,                          \
                 const unsigned int MM, const unsigned int DD, \
                 const double LON, const double LAT,           \
                 const double P_hPa,                           \
                 double NOON_JRATES[] )
{

    /* Set up for error messages */
    // static const int NC_ERR = 2;

    std::string fullPath;
    std::stringstream mm, dd;
    std::string varName;
    //double lonIn[72];
    //double latIn[46];
    //double pMidIn[59];
    double lonIn[72];
    double latIn[46];
    double pMidIn[59];

    NcVar data;
   
    mm << std::setw(2) << std::setfill('0') << MM;
    dd << std::setw(2) << std::setfill('0') << DD;
    fullPath += ROOTDIR;
    fullPath += "/JData_2013-" + mm.str() + "-" + dd.str() + ".nc";

    //NcFile dataFile( fullPath.c_str(), NcFile::ReadOnly );
    //if ( !dataFile.is_valid() ) {
    //    std::cout << " Photolysis rate input file '" << fullPath << "' not found!" << std::endl;
    //    exit(-1);
    //}
    NcFile dataFile( fullPath.c_str(), NcFile::read );

    varName = "lon";
    data = dataFile.getVar(varName.c_str());
    data.getVar(lonIn);

    varName = "lat";
    data = dataFile.getVar(varName.c_str());
    data.getVar(latIn);

    varName = "pmid";
    data = dataFile.getVar(varName.c_str());
    data.getVar(pMidIn);

    Vector_1D lon(lonIn, lonIn+72);
    Vector_1D lat(latIn, latIn+46);
    Vector_1D pmid(pMidIn, pMidIn+59);

    /* Find values in table closest to local longitude, latitude and pressure */
    const auto LON_IT = std::lower_bound( lon.begin(), lon.end(), LON );
    const unsigned int LON_SIZE = lon.size();
    const unsigned int LON_INDEX = std::distance( lon.begin(), LON_IT );
    bool LON_EDGE = 0;
    double LON_LOW, LON_HIGH;

    /*
    Initialize LON_LOW to help compiler understand it's never uninitialized when used
    This value will be overwritten before LON_LOW is called
    */
    LON_LOW = 0; 

    if ( ( LON_INDEX == 0 ) || ( LON_INDEX + 1 == LON_SIZE ) )
        LON_EDGE = 1;

    if ( !LON_EDGE )
        LON_LOW = lon[LON_INDEX - 1];
    LON_HIGH = lon[LON_INDEX];
    
    const auto LAT_IT = std::lower_bound( lat.begin(), lat.end(), LAT );
    const unsigned int LAT_SIZE = lat.size();
    const unsigned int LAT_INDEX = std::distance( lat.begin(), LAT_IT );
    bool LAT_EDGE = 0;
    double LAT_LOW, LAT_HIGH;

    if ( ( LAT_INDEX == 0 ) || ( LAT_INDEX + 1 == LAT_SIZE ) )
        LAT_EDGE = 1;
    
    if ( !LAT_EDGE )
        LAT_LOW = lat[LAT_INDEX - 1];
    LAT_HIGH = lat[LAT_INDEX];
    
    /* Since PMID is sorted in decreasing order, reverse it to find the correct index, ... */
    std::reverse(pmid.begin(), pmid.end());
    const auto PRES_IT = std::lower_bound( pmid.begin(), pmid.end(), P_hPa );
    /* ... and then reverse back */
    std::reverse(pmid.begin(), pmid.end());
    const unsigned int PRES_SIZE = pmid.size();
    const unsigned int PRES_INDEX = PRES_SIZE - std::distance( pmid.begin(), PRES_IT );

    const double PRES_HIGH = pmid[PRES_INDEX - 1];
    const double PRES_LOW  = pmid[PRES_INDEX];

    //double jRate[59][46][72];
    double jRate[59][46][72];
    double jRateNEL, jRateNEH, jRateNWL, jRateNWH, jRateSEL, jRateSEH, jRateSWL, jRateSWH;

    /* Read the photolysis rates ... */

    unsigned int iPhotol = 0;

    /* O2_J1 */
    varName = "O2_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* O3_J1 */
    varName = "O3_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* O3_J2 */
    varName = "O3_J2";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* H2O_J1 */
    varName = "H2O_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* HO2_J1 */
    varName = "HO2_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* NO_J1 */
    varName = "NO_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* CH2O_J1 */
    varName = "CH2O_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* CH2O_J2 */
    varName = "CH2O_J2";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* H2O2_J1 */
    varName = "H2O2_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* MP_J1 */
    varName = "MP_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* NO2_J1 */
    varName = "NO2_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* NO3_J1 */
    varName = "NO3_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* NO3_J2 */
    varName = "NO3_J2";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* N2O5_J1 */
    varName = "N2O5_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* HNO2_J1 */
    varName = "HNO2_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* HNO3_J1 */
    varName = "HNO3_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* HNO4_J1 */
    varName = "HNO4_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* HNO4_J2 */
    varName = "HNO4_J2";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* ClNO3_J1 */
    varName = "ClNO3_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* ClNO3_J2 */
    varName = "ClNO3_J2";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* ClNO2_J1 */
    varName = "ClNO2_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* Cl2_J1 */
    varName = "Cl2_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* Br2_J1 */
    varName = "Br2_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* HOCl_J1 */
    varName = "HOCl_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* OClO_J1 */
    varName = "OClO_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* Cl2O2_J1 */
    varName = "Cl2O2_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* ClO_J1 */
    varName = "ClO_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* BrO_J1 */
    varName = "BrO_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* BrNO3_J1 */
    varName = "BrNO3_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* BrNO3_J2 */
    varName = "BrNO3_J2";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* BrNO2_J1 */
    varName = "BrNO2_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* HOBr_J1 */
    varName = "HOBr_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* BrCl_J1 */
    varName = "BrCl_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* OCS_J1 */
    varName = "OCS_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* SO2_J1 */
    varName = "SO2_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* N2O_J1 */
    varName = "N2O_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* CFC11_J1 */
    varName = "CFC11_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* CFC12_J1 */
    varName = "CFC12_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* CFC113_J1 */
    varName = "CFC113_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* CFC114_J1 */
    varName = "CFC114_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* CFC115_J1 */
    varName = "CFC115_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* CCl4_J1 */
    varName = "CCl4_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* CH3Cl_J1 */
    varName = "CH3Cl_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* CH3CCl3_J1 */
    varName = "CH3CCl3_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* CH2Cl2_J1 */
    varName = "CH2Cl2_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* HCFC22_J1 */
    varName = "HCFC22_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* HCFC123_J1 */
    varName = "HCFC123_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* HCFC141b_J1 */
    varName = "HCFC141b_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* HCFC142b_J1 */
    varName = "HCFC142b_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* CH3Br_J1 */
    varName = "CH3Br_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* H1211_J1 */
    varName = "H1211_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* H12O2_J1 */
    varName = "H12O2_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* H1301_J1 */
    varName = "H1301_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* H2402_J1 */
    varName = "H2402_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* CH2Br2_J1 */
    varName = "CH2Br2_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* CHBr3_J1 */
    varName = "CHBr3_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* CH3I_J1 */
    varName = "CH3I_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* CF3I_J1 */
    varName = "CF3I_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* PAN_J1 */
    varName = "PAN_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* R4N2_J1 */
    varName = "R4N2_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* ALD2_J1 */
    varName = "ALD2_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* ALD2_J2 */
    varName = "ALD2_J2";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* MVK_J1 */
    varName = "MVK_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* MVK_J2 */
    varName = "MVK_J2";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* MVK_J3 */
    varName = "MVK_J3";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* MACR_J1 */
    varName = "MACR_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* MACR_J2 */
    varName = "MACR_J2";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* GLYC_J1 */
    varName = "GLYC_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* MEK_J1 */
    varName = "MEK_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* RCHO_J1 */
    varName = "RCHO_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* MGLY_J1 */
    varName = "MGLY_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* GLYX_J1 */
    varName = "GLYX_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* GLYX_J2 */
    varName = "GLYX_J2";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* GLYX_J3 */
    varName = "GLYX_J3";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* HAC_J1 */
    varName = "HAC_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* ACET_J1 */
    varName = "ACET_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* ACET_J2 */
    varName = "ACET_J2";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* INPN_J1 */
    varName = "INPN_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* PRPN_J1 */
    varName = "PRPN_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* ETP_J1 */
    varName = "ETP_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* RA3P_J1 */
    varName = "RA3P_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* RB3P_J1 */
    varName = "RB3P_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* R4P_J1 */
    varName = "R4P_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* PP_J1 */
    varName = "PP_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* RP_J1 */
    varName = "RP_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* RIP_J1 */
    varName = "RIP_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* IAP_J1 */
    varName = "IAP_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* ISNP_J1 */
    varName = "ISNP_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* VRP_J1 */
    varName = "VRP_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* MRP_J1 */
    varName = "MRP_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* MAOP_J1 */
    varName = "MAOP_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* MACRN_J1 */
    varName = "MACRN_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* MVKN_J1 */
    varName = "MVKN_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* ISOPNB_J1 */
    varName = "ISOPNB_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* ISOPND_J1 */
    varName = "ISOPND_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* PROPNN_J1 */
    varName = "PROPNN_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* ATOOH_J1 */
    varName = "ATOOH_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* R4N2_J1 */
    varName = "R4N2_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* MAP_J1 */
    varName = "MAP_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* SO4_J1 */
    varName = "SO4_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* ClNO2_J2 */
    varName = "ClNO2_J2";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* ClOO_J1 */
    varName = "ClOO_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* O3_J3 */
    varName = "O3_J3";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* MPN_J1 */
    varName = "MPN_J1";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    
    /* MPN_J2 */
    varName = "MPN_J2";
    //data = dataFile.get_var(varName.c_str());
    //data.getVar(&jRate[0][0][0], 59, 46, 72);
    data = dataFile.getVar(varName.c_str());
    data.getVar(jRate);
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                 ( PRES_HIGH - P_hPa ) * jRateNEH ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateSEH ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                        ( PRES_HIGH - P_hPa ) * jRateNEH ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        jRateSWL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX-1];
        jRateSWH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX-1];
        jRateSEL = jRate[PRES_INDEX-1][LAT_INDEX-1][LON_INDEX  ];
        jRateSEH = jRate[PRES_INDEX  ][LAT_INDEX-1][LON_INDEX  ];
        jRateNWL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX-1];
        jRateNWH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX-1];
        jRateNEL = jRate[PRES_INDEX-1][LAT_INDEX  ][LON_INDEX  ];
        jRateNEH = jRate[PRES_INDEX  ][LAT_INDEX  ][LON_INDEX  ];
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateSWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateSEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateSEH ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRateNWL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNWH ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRateNEL + \
                                                                               ( PRES_HIGH - P_hPa ) * jRateNEH ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    iPhotol++;
    

} /* End of ReadJRates */

/* End of ReadJRates.cpp */
