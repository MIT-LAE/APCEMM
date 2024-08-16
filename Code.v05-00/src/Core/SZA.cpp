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

#include <cmath>
#include "Util/PhysConstant.hpp"
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

    sinLAT = std::sin(lat_ * physConst::PI/180.0);
    cosLAT = std::cos(lat_ * physConst::PI/180.0);
    sinDEC = std::sin(DEC);
    cosDEC = std::cos(DEC);

    sunRise = std::max((12.0 - 180.0/(physConst::PI*15.0)*acos(-(sinLAT * sinDEC)\
                                                        / (cosLAT * cosDEC))), 0.0);
    sunSet  = std::min((12.0 + 180.0/(physConst::PI*15.0)*acos(-(sinLAT * sinDEC)\
                                              / (cosLAT * cosDEC))), 24.0);

    CSZA_max = std::max( sinLAT * sinDEC + cosLAT * cosDEC, 0.0 );
    CSZA = 0.0;

    /* Vector of coefficient 
     * CSZA = CSZA_Vector[0] + CSZA_Vector[1] * cos( CSZA_Vector{3} * (Time/3600.0 - 12.0 ) )*/
    CSZA_Vector = { sinLAT * sinDEC, cosLAT * cosDEC, 15.0 * physConst::PI / 180.0 };

} /* End of SZA::SZA */

SZA::~SZA( )
{

    /* Destructor */

} /* End of SZA::~SZA */

void SZA::Update( const double solarTime )
{

    CSZA = std::max( sinLAT * sinDEC + cosLAT * cosDEC * std::cos( std::abs( ( solarTime/3600.0 - 12.0 ) ) * 15.0 * physConst::PI / 180.0 ), 0.0 ); 

} /* End of SZA::Update */

double SZA::getCSZA( const double solarTime )
{

    return std::max( sinLAT * sinDEC + cosLAT * cosDEC * std::cos( std::abs( ( solarTime/3600.0 - 12.0 ) ) * 15.0 * physConst::PI / 180.0 ), 0.0 ); 

} /* End of SZA::getCSZA */

std::vector<double> SZA::getSZA_Vector( ) const
{

    return CSZA_Vector;

} /* End of SZA::getSZA_Vector */

void* castSZA( SZA *sun )
{

    return( reinterpret_cast< void* >( sun ) );

} /* End of castSZA */

double castCSZA( void *sun, const double solarTime )
{
    
    return( reinterpret_cast< SZA* >( sun )->getCSZA( solarTime ) );

} /* End of castCSZA */

/* End of SZA.cpp */
