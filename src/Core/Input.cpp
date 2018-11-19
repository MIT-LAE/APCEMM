/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Input Program File                                               */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 11/19/2018                                */
/* File                 : Input.cpp                                 */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "Core/Input.hpp"

Input::Input( const RealDouble temperature_K, \
              const RealDouble pressure_Pa,   \
              const RealDouble relHumidity_w, \
              const RealDouble long_deg,      \
              const RealDouble lat_deg,       \
              const UInt dayGMT,              \
              const RealDouble emissionTime,  \
              const Vector_1D emissionInput):
    temperature_K_ ( temperature_K ),
    pressure_Pa_   ( pressure_Pa) ,
    relHumidity_w_ ( relHumidity_w ),
    longitude_deg_ ( long_deg ),
    latitude_deg_  ( lat_deg ),
    dayGMT_        ( dayGMT ),
    emissionTime_  ( emissionTime ),
    EI_NOx_        ( emissionInput[0] ),
    EI_CO_         ( emissionInput[1] ),
    EI_HC_         ( emissionInput[2] ),
    EI_Soot_       ( emissionInput[3] ),
    sootRad_       ( emissionInput[4] ),
    fuelFlow_      ( emissionInput[5] )
{

    /* Constructor */

    while ( dayGMT_ > 365 )
        dayGMT_ -= 365;

    while ( emissionTime_ > 24.0 )
        emissionTime_ -= 24.0;

} /* End of Input::Input */

Input::Input( unsigned int iCase, \
              const Vector_2D &parameters ):
    temperature_K_ ( parameters[0][iCase] ),
    pressure_Pa_   ( parameters[1][iCase] ),
    relHumidity_w_ ( parameters[2][iCase] ),
    longitude_deg_ ( parameters[3][iCase] ),
    latitude_deg_  ( parameters[4][iCase] ),
    dayGMT_        ( parameters[5][iCase] ),
    emissionTime_  ( parameters[6][iCase] ),
    EI_NOx_        ( parameters[7][iCase] ),
    EI_CO_         ( parameters[8][iCase] ),
    EI_HC_         ( parameters[9][iCase] ),
    EI_Soot_       ( parameters[10][iCase] ),
    sootRad_       ( parameters[11][iCase] ),
    fuelFlow_      ( parameters[12][iCase] )
{

    /* Constructor */
  
    while ( dayGMT_ > 365 )
        dayGMT_ -= 365;

    while ( emissionTime_ > 24.0 )
        emissionTime_ -= 24.0;
    
} /* End of Input::Input */

Input::~Input()
{

    /* Destructor */

} /* End of Input::~Input */

/* End of Input.cpp */
