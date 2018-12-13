/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Input_Mod Program File                                           */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 12/12/2018                                */
/* File                 : Input_Mod.cpp                             */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "Core/Input_Mod.hpp"

OptInput::OptInput( ):
    PARAMETER_TEMPERATURE( 0 ),
    PARAMETER_TEMPERATURE_UNIT( "" ),
    PARAMETER_RHW( 0 ),
    PARAMETER_RHW_UNIT( "" ),
    PARAMETER_LONGITUDE( 0 ),
    PARAMETER_LONGITUDE_UNIT( "" ),
    PARAMETER_LATITUDE( 0 ),
    PARAMETER_LATITUDE_UNIT( "" ),
    PARAMETER_PRESSURE( 0 ),
    PARAMETER_PRESSURE_UNIT( "" ),
    PARAMETER_EDAY( 0 ),
    PARAMETER_EDAY_UNIT( "" ),
    PARAMETER_ETIME( 0 ),
    PARAMETER_ETIME_UNIT( "" ),
    PARAMETER_BACKG_NOX( 0 ),
    PARAMETER_BACKG_NOX_UNIT( "" ),
    PARAMETER_BACKG_HNO3( 0 ),
    PARAMETER_BACKG_HNO3_UNIT( "" ),
    PARAMETER_BACKG_O3( 0 ),
    PARAMETER_BACKG_O3_UNIT( "" ),
    PARAMETER_BACKG_CO( 0 ),
    PARAMETER_BACKG_CO_UNIT( "" ),
    PARAMETER_BACKG_CH4( 0 ),
    PARAMETER_BACKG_CH4_UNIT( "" ),
    PARAMETER_BACKG_SO2( 0 ),
    PARAMETER_BACKG_SO2_UNIT( "" ),
    PARAMETER_EI_NOX( 0 ),
    PARAMETER_EI_NOX_UNIT( "" ),
    PARAMETER_EI_CO( 0 ),
    PARAMETER_EI_CO_UNIT( "" ),
    PARAMETER_EI_UHC( 0 ),
    PARAMETER_EI_UHC_UNIT( "" ),
    PARAMETER_EI_SO2( 0 ),
    PARAMETER_EI_SO2_UNIT( "" ),
    PARAMETER_EI_SO2TOSO4( 0 ),
    PARAMETER_EI_SO2TOSO4_UNIT( "" ),
    PARAMETER_EI_SOOT( 0 ),
    PARAMETER_EI_SOOT_UNIT( "" ),
    PARAMETER_EI_SOOTRAD( 0 ),
    PARAMETER_EI_SOOTRAD_UNIT( "" )
{

    /* Default constructor */

} /* End of OptInput::OptInput */

OptInput::~OptInput( )
{

    /* Destructor */

} /* End of OptInput::~OptInput */

/* End of Input_Mod.cpp */

