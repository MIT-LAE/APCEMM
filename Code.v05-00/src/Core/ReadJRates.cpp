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

    std::string fullPath_Grid;
    std::string fullPath;
    std::stringstream mm, dd;
   
    mm << std::setw(2) << std::setfill('0') << MM;
    dd << std::setw(2) << std::setfill('0') << DD;
    fullPath += ROOTDIR;
    fullPath += "/JData_2013-" + mm.str() + "-" + dd.str() + ".mat";
    fullPath_Grid += ROOTDIR;
    fullPath_Grid += "/Grid.mat";

    /* First deal with grid to find indices to evaluate the rates at */

    MATFile *pmat;
    mxArray *mxLON, *mxLAT, *mxPMID;
    Vector_1D lon, lat, pmid;

    /* Open the mat file */
    pmat = matOpen( fullPath_Grid.c_str(), "r" );
    if ( pmat == NULL ) {
        std::cout << "Error opening JRate grid file: " << fullPath << std::endl;
        exit(-1);
    }

    /* Read the desired variables: LON, LAT and PMID */
    mxLON  = matGetVariable( pmat, "LON"  );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            lon.reserve( num );
            lon.assign(pr, pr+num);
        }
    }
    mxDestroyArray( mxLON );

    mxLAT  = matGetVariable( pmat, "LAT"  );
    if ( ( mxLAT != NULL ) && mxIsDouble( mxLAT ) && \
            !mxIsEmpty( mxLAT ) ) {
        mwSize num = mxGetNumberOfElements( mxLAT );
        double *pr = mxGetPr( mxLAT );
        if ( pr != NULL ) {
            lat.reserve( num );
            lat.assign(pr, pr+num);
        }
    }
    mxDestroyArray( mxLAT );

    mxPMID = matGetVariable( pmat, "PMID" );
    if ( ( mxPMID != NULL ) && mxIsDouble( mxPMID ) && \
            !mxIsEmpty( mxPMID ) ) {
        mwSize num = mxGetNumberOfElements( mxPMID );
        double *pr = mxGetPr( mxPMID );
        if ( pr != NULL ) {
            pmid.reserve( num );
            pmid.assign(pr, pr+num);
        }
    }
    mxDestroyArray( mxPMID );

    /* Close mat file */
    matClose( pmat );

    /* Find values in table closest to local longitude, latitude and pressure */
    const auto LON_IT = std::lower_bound( lon.begin(), lon.end(), LON );
    const unsigned int LON_SIZE = lon.size();
    const unsigned int LON_INDEX = std::distance( lon.begin(), LON_IT );
    bool LON_EDGE = 0;
    double LON_LOW, LON_HIGH;

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

    mxArray *mxJRate;
    Vector_1D jRate;

    /* Open the mat file */
    pmat = matOpen( fullPath.c_str(), "r" );
    if ( pmat == NULL ) {
        std::cout << "Error opening JRate file: " << fullPath << std::endl;
        exit(-1);
    }

    /* Read the photolysis rates ... */

    unsigned int iPhotol = 0;

    /* O2_J1 */
    mxJRate = matGetVariable( pmat, "O2_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    /* ... and interpolate */
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* O3_J1 */
    mxJRate = matGetVariable( pmat, "O3_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;

    /* O3_J2 */
    mxJRate = matGetVariable( pmat, "O3_J2" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;

    /* H2O_J1 */
    mxJRate = matGetVariable( pmat, "H2O_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;

    /* HO2_J1 */
    mxJRate = matGetVariable( pmat, "HO2_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;

    /* NO_J1 */
    mxJRate = matGetVariable( pmat, "NO_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;

    /* CH2O_J1 */
    mxJRate = matGetVariable( pmat, "CH2O_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;

    /* CH2O_J2 */
    mxJRate = matGetVariable( pmat, "CH2O_J2" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;

    /* H2O2_J1 */
    mxJRate = matGetVariable( pmat, "H2O2_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;

    /* MP_J1 */
    mxJRate = matGetVariable( pmat, "MP_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* NO2_J1 */
    mxJRate = matGetVariable( pmat, "NO2_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;

    /* NO3_J1 */
    mxJRate = matGetVariable( pmat, "NO3_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* NO3_J2 */
    mxJRate = matGetVariable( pmat, "NO3_J2" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* N2O5_J1 */
    mxJRate = matGetVariable( pmat, "N2O5_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* HNO2_J1 */
    mxJRate = matGetVariable( pmat, "HNO2_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* HNO3_J1 */
    mxJRate = matGetVariable( pmat, "HNO3_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;

    /* HNO4_J1 */
    mxJRate = matGetVariable( pmat, "HNO4_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;

    /* HNO4_J2 */
    mxJRate = matGetVariable( pmat, "HNO4_J2" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;

    /* ClNO3_J1 */
    mxJRate = matGetVariable( pmat, "ClNO3_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;

    /* ClNO3_J2 */
    mxJRate = matGetVariable( pmat, "ClNO3_J2" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;

    /* ClNO2_J1 */
    mxJRate = matGetVariable( pmat, "ClNO2_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;

    /* Cl2_J1 */
    mxJRate = matGetVariable( pmat, "Cl2_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;

    /* Br2_J1 */
    mxJRate = matGetVariable( pmat, "Br2_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* HOCl_J1 */
    mxJRate = matGetVariable( pmat, "HOCl_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* OClO_J1 */
    mxJRate = matGetVariable( pmat, "OClO_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* Cl2O2_J1 */
    mxJRate = matGetVariable( pmat, "Cl2O2_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* ClO_J1 */
    mxJRate = matGetVariable( pmat, "ClO_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* BrO_J1 */
    mxJRate = matGetVariable( pmat, "BrO_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* BrNO3_J1 */
    mxJRate = matGetVariable( pmat, "BrNO3_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* BrNO3_J2 */
    mxJRate = matGetVariable( pmat, "BrNO3_J2" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* BrNO2_J1 */
    mxJRate = matGetVariable( pmat, "BrNO2_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* HOBr_J1 */
    mxJRate = matGetVariable( pmat, "HOBr_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* BrCl_J1 */
    mxJRate = matGetVariable( pmat, "BrCl_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* OCS_J1 */
    mxJRate = matGetVariable( pmat, "OCS_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* SO2_J1 */
    mxJRate = matGetVariable( pmat, "SO2_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* N2O_J1 */
    mxJRate = matGetVariable( pmat, "N2O_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* CFC11_J1 */
    mxJRate = matGetVariable( pmat, "CFC11_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* CFC12_J1 */
    mxJRate = matGetVariable( pmat, "CFC12_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* CFC113_J1 */
    mxJRate = matGetVariable( pmat, "CFC113_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* CFC114_J1 */
    mxJRate = matGetVariable( pmat, "CFC114_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* CFC115_J1 */
    mxJRate = matGetVariable( pmat, "CFC115_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* CCl4_J1 */
    mxJRate = matGetVariable( pmat, "CCl4_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* CH3Cl_J1 */
    mxJRate = matGetVariable( pmat, "CH3Cl_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* CH3CCl3_J1 */
    mxJRate = matGetVariable( pmat, "CH3CCl3_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* CH2Cl2_J1 */
    mxJRate = matGetVariable( pmat, "CH2Cl2_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* HCFC22_J1 */
    mxJRate = matGetVariable( pmat, "HCFC22_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* HCFC123_J1 */
    mxJRate = matGetVariable( pmat, "HCFC123_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* HCFC141b_J1 */
    mxJRate = matGetVariable( pmat, "HCFC141b_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* HCFC142b_J1 */
    mxJRate = matGetVariable( pmat, "HCFC142b_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* CH3Br_J1 */
    mxJRate = matGetVariable( pmat, "CH3Br_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* H1211_J1 */
    mxJRate = matGetVariable( pmat, "H1211_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* H12O2_J1 */
    mxJRate = matGetVariable( pmat, "H12O2_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* H1301_J1 */
    mxJRate = matGetVariable( pmat, "H1301_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* H2402_J1 */
    mxJRate = matGetVariable( pmat, "H2402_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* CH2Br2_J1 */
    mxJRate = matGetVariable( pmat, "CH2Br2_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* CHBr3_J1 */
    mxJRate = matGetVariable( pmat, "CHBr3_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* CH3I_J1 */
    mxJRate = matGetVariable( pmat, "CH3I_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* CF3I_J1 */
    mxJRate = matGetVariable( pmat, "CF3I_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;

    /* PAN_J1 */
    mxJRate = matGetVariable( pmat, "PAN_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;

    /* R4N2_J1 */
    mxJRate = matGetVariable( pmat, "R4N2_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;

    /* ALD2_J1 */
    mxJRate = matGetVariable( pmat, "ALD2_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;

    /* ALD2_J2 */
    mxJRate = matGetVariable( pmat, "ALD2_J2" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;

    /* MVK_J1 */
    mxJRate = matGetVariable( pmat, "MVK_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;

    /* MVK_J2 */
    mxJRate = matGetVariable( pmat, "MVK_J2" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* MVK_J3 */
    mxJRate = matGetVariable( pmat, "MVK_J3" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;

    /* MACR_J1 */
    mxJRate = matGetVariable( pmat, "MACR_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;

    /* MACR_J2 */
    mxJRate = matGetVariable( pmat, "MACR_J2" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;

    /* GLYC_J1 */
    mxJRate = matGetVariable( pmat, "GLYC_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;

    /* MEK_J1 */
    mxJRate = matGetVariable( pmat, "MEK_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;

    /* RCHO_J1 */
    mxJRate = matGetVariable( pmat, "RCHO_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;

    /* MGLY_J1 */
    mxJRate = matGetVariable( pmat, "MGLY_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;

    /* GLYX_J1 */
    mxJRate = matGetVariable( pmat, "GLYX_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;

    /* GLYX_J2 */
    mxJRate = matGetVariable( pmat, "GLYX_J2" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;

    /* GLYX_J3 */
    mxJRate = matGetVariable( pmat, "GLYX_J3" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;

    /* HAC_J1 */
    mxJRate = matGetVariable( pmat, "HAC_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;

    /* ACET_J1 */
    mxJRate = matGetVariable( pmat, "ACET_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;

    /* ACET_J2 */
    mxJRate = matGetVariable( pmat, "ACET_J2" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;

    /* INPN_J1 */
    mxJRate = matGetVariable( pmat, "INPN_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;

    /* PRPN_J1 */
    mxJRate = matGetVariable( pmat, "PRPN_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;

    /* ETP_J1 */
    mxJRate = matGetVariable( pmat, "ETP_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* RA3P_J1 */
    mxJRate = matGetVariable( pmat, "RA3P_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* RB3P_J1 */
    mxJRate = matGetVariable( pmat, "RB3P_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* R4P_J1 */
    mxJRate = matGetVariable( pmat, "R4P_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* PP_J1 */
    mxJRate = matGetVariable( pmat, "PP_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* RP_J1 */
    mxJRate = matGetVariable( pmat, "RP_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* RIP_J1 */
    mxJRate = matGetVariable( pmat, "RIP_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* IAP_J1 */
    mxJRate = matGetVariable( pmat, "IAP_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* ISNP_J1 */
    mxJRate = matGetVariable( pmat, "ISNP_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* VRP_J1 */
    mxJRate = matGetVariable( pmat, "VRP_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* MRP_J1 */
    mxJRate = matGetVariable( pmat, "MRP_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* MAOP_J1 */
    mxJRate = matGetVariable( pmat, "MAOP_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* MACRN_J1 */
    mxJRate = matGetVariable( pmat, "MACRN_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* MVKN_J1 */
    mxJRate = matGetVariable( pmat, "MVKN_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* ISOPNB_J1 */
    mxJRate = matGetVariable( pmat, "ISOPNB_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* ISOPND_J1 */
    mxJRate = matGetVariable( pmat, "ISOPND_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* PROPNN_J1 */
    mxJRate = matGetVariable( pmat, "PROPNN_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* ATOOH_J1 */
    mxJRate = matGetVariable( pmat, "ATOOH_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* R4N2_J1 */
    mxJRate = matGetVariable( pmat, "R4N2_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* MAP_J1 */
    mxJRate = matGetVariable( pmat, "MAP_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* SO4_J1 */
    mxJRate = matGetVariable( pmat, "SO4_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* ClNO2_J2 */
    mxJRate = matGetVariable( pmat, "ClNO2_J2" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* ClOO_J1 */
    mxJRate = matGetVariable( pmat, "ClOO_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* O3_J3 */
    mxJRate = matGetVariable( pmat, "O3_J3" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* MPN_J1 */
    mxJRate = matGetVariable( pmat, "MPN_J1" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );
    iPhotol++;
    
    /* MPN_J2 */
    mxJRate = matGetVariable( pmat, "MPN_J2" );
    if ( ( mxLON != NULL ) && mxIsDouble( mxLON ) && \
            !mxIsEmpty( mxLON ) ) {
        mwSize num = mxGetNumberOfElements( mxLON );
        double *pr = mxGetPr( mxLON );
        if ( pr != NULL ) {
            jRate.reserve( num );
            jRate.assign(pr, pr+num);
        }
    }
    
    if ( LAT_EDGE && LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                 ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) /\
                               ( PRES_HIGH - PRES_LOW );
    } else if ( LAT_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + LAT_INDEX * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW );
    } else if ( LON_EDGE ) {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                 ( LAT - LAT_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                        ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LAT_HIGH - LAT_LOW );
    } else {
        NOON_JRATES[iPhotol] = ( ( LAT_HIGH - LAT ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX - 1 ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) + \
                                 ( LAT - LAT_LOW  ) * ( ( LON_HIGH - LON ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX - 1 + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) + \
                                                        ( LON - LON_LOW  ) * ( ( P_hPa - PRES_LOW  ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX - 1 ) * LON_SIZE * LAT_SIZE] + \
                                                                               ( PRES_HIGH - P_hPa ) * jRate[LON_INDEX     + ( LAT_INDEX     ) * LON_SIZE + ( PRES_INDEX     ) * LON_SIZE * LAT_SIZE] ) ) ) /\
                               ( PRES_HIGH - PRES_LOW ) / ( LON_HIGH - LON_LOW ) / ( LAT_HIGH - LAT_LOW );
    }

    jRate.clear();
    mxDestroyArray( mxJRate );

    /* Close mat file */
    matClose( pmat );

} /* End of ReadJRates */

/* End of ReadJRates.cpp */
