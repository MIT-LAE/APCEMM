/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* SZA Program File                                                 */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : SZA.cpp                                   */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "Core/SZA.hpp"

SZA::SZA( const double lat_, const unsigned int day_ ):
    latitude( lat_ ),
    dayGMT( day_ )
{

    /* Default Constructor */

    const double r_SZA = 2*physConst::PI*(std::floor(dayGMT) - 1)/365.0;

    DEC = A0 - A1*cos(1.0*r_SZA) + B1*sin(1.0*r_SZA)\
             - A2*cos(2.0*r_SZA) + B2*sin(2.0*r_SZA)\
             - A3*cos(3.0*r_SZA) + B3*sin(3.0*r_SZA);

    SINLAT = std::sin(lat_ * physConst::PI/180.0);
    COSLAT = std::cos(lat_ * physConst::PI/180.0);
    SINDEC = std::sin(DEC);
    COSDEC = std::cos(DEC);

    sunRise = std::max((12.0 - 180.0/(physConst::PI*15.0)*acos(-(SINLAT * SINDEC)\
                                                        / (COSLAT * COSDEC))), 0.0);
    sunSet  = std::min((12.0 + 180.0/(physConst::PI*15.0)*acos(-(SINLAT * SINDEC)\
                                              / (COSLAT * COSDEC))), 24.0);

    CSZA = 0.0;

} /* End of SZA::SZA */

SZA::~SZA( )
{

    /* Destructor */

} /* End of SZA::~SZA */

void SZA::Update( const double solarTime )
{

    CSZA = std::max( SINLAT * SINDEC + COSLAT * COSDEC * std::cos( std::abs( ( solarTime/3600.0 - 12.0 ) ) * 15.0 * physConst::PI / 180.0 ), 0.0 ); 

} /* End of SZA::Update_SUN */

/* End of SZA.cpp */
