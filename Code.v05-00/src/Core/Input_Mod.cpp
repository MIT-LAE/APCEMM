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
    SIMULATION_PARAMETER_SWEEP( 0 ),
    SIMULATION_MONTECARLO( 0 ),
    SIMULATION_MCRUNS( 0 ),
    SIMULATION_OUTPUT_FOLDER( "" ),
    SIMULATION_OVERWRITE( 0 ),
    SIMULATION_RUN_DIRECTORY( "" ),
    SIMULATION_INPUT_BACKG_COND( "" ),
    SIMULATION_SAVE_FORWARD( 0 ),
    SIMULATION_FORWARD_FILENAME( "" ),
    SIMULATION_ADJOINT( 0 ),
    SIMULATION_ADJOINT_FILENAME( "" ),
    PARAMETER_TEMPERATURE( 0 ),
    PARAMETER_TEMPERATURE_UNIT( "" ),
    PARAMETER_TEMPERATURE_RANGE( 0 ),
    PARAMETER_RHW( 0 ),
    PARAMETER_RHW_UNIT( "" ),
    PARAMETER_RHW_RANGE( 0 ),
    PARAMETER_LONGITUDE( 0 ),
    PARAMETER_LONGITUDE_UNIT( "" ),
    PARAMETER_LONGITUDE_RANGE( 0 ),
    PARAMETER_LATITUDE( 0 ),
    PARAMETER_LATITUDE_UNIT( "" ),
    PARAMETER_LATITUDE_RANGE( 0 ),
    PARAMETER_PRESSURE( 0 ),
    PARAMETER_PRESSURE_UNIT( "" ),
    PARAMETER_PRESSURE_RANGE( 0 ),
    PARAMETER_EDAY( 0 ),
    PARAMETER_EDAY_UNIT( "" ),
    PARAMETER_EDAY_RANGE( 0 ),
    PARAMETER_ETIME( 0 ),
    PARAMETER_ETIME_UNIT( "" ),
    PARAMETER_ETIME_RANGE( 0 ),
    PARAMETER_BACKG_NOX( 0 ),
    PARAMETER_BACKG_NOX_UNIT( "" ),
    PARAMETER_BACKG_NOX_RANGE( 0 ),
    PARAMETER_BACKG_HNO3( 0 ),
    PARAMETER_BACKG_HNO3_UNIT( "" ),
    PARAMETER_BACKG_HNO3_RANGE( 0 ),
    PARAMETER_BACKG_O3( 0 ),
    PARAMETER_BACKG_O3_UNIT( "" ),
    PARAMETER_BACKG_O3_RANGE( 0 ),
    PARAMETER_BACKG_CO( 0 ),
    PARAMETER_BACKG_CO_UNIT( "" ),
    PARAMETER_BACKG_CO_RANGE( 0 ),
    PARAMETER_BACKG_CH4( 0 ),
    PARAMETER_BACKG_CH4_UNIT( "" ),
    PARAMETER_BACKG_CH4_RANGE( 0 ),
    PARAMETER_BACKG_SO2( 0 ),
    PARAMETER_BACKG_SO2_UNIT( "" ),
    PARAMETER_BACKG_SO2_RANGE( 0 ),
    PARAMETER_EI_NOX( 0 ),
    PARAMETER_EI_NOX_UNIT( "" ),
    PARAMETER_EI_NOX_RANGE( 0 ),
    PARAMETER_EI_CO( 0 ),
    PARAMETER_EI_CO_UNIT( "" ),
    PARAMETER_EI_CO_RANGE( 0 ),
    PARAMETER_EI_UHC( 0 ),
    PARAMETER_EI_UHC_UNIT( "" ),
    PARAMETER_EI_UHC_RANGE( 0 ),
    PARAMETER_EI_SO2( 0 ),
    PARAMETER_EI_SO2_UNIT( "" ),
    PARAMETER_EI_SO2_RANGE( 0 ),
    PARAMETER_EI_SO2TOSO4( 0 ),
    PARAMETER_EI_SO2TOSO4_UNIT( "" ),
    PARAMETER_EI_SO2TOSO4_RANGE( 0 ),
    PARAMETER_EI_SOOT( 0 ),
    PARAMETER_EI_SOOT_UNIT( "" ),
    PARAMETER_EI_SOOT_RANGE( 0 ),
    PARAMETER_EI_SOOTRAD( 0 ),
    PARAMETER_EI_SOOTRAD_UNIT( "" ),
    PARAMETER_EI_SOOTRAD_RANGE( 0 ),
    PARAMETER_FF( 0 ),
    PARAMETER_FF_UNIT( "" ),
    PARAMETER_FF_RANGE( 0 ),
    TRANSPORT_TRANSPORT( 0 ),
    TRANSPORT_FILL( 0 ),
    TRANSPORT_TIMESTEP( 0.0E+00 ),
    CHEMISTRY_CHEMISTRY( 0 ),
    CHEMISTRY_HETCHEM( 0 ),
    CHEMISTRY_JRATE_FOLDER( "" ),
    CHEMISTRY_TIMESTEP( 0.0E+00 ),
    AEROSOL_GRAVSETTLING( 0 ),
    AEROSOL_COAGULATION( 0 ),
    AEROSOL_COAGULATION_TIMESTEP( 0.0E+00 ),
    AEROSOL_ICE_GROWTH( 0 ),
    AEROSOL_PLUME_UPDRAFT( 0 ),
    MET_TEMP_INIT( 0 ),
    MET_H2O_INIT( 0 ),
    MET_FILENAME( "" ),
    DIAG_FILENAME( "" ),
    TS_SPEC( 0 ),
    TS_FILENAME( "" ),
    TS_SPECIES( 0 ),
    TS_FREQ( 0.0E+00 ),
    TS_AERO( 0 ),
    TS_AERO_FILENAME( "" ),
    TS_AEROSOL( 0 ),
    TS_AERO_FREQ( 0.0E+00 ),
    PL_PL( 0 ),
    PL_O3( 0 )
{

    /* Default constructor */

} /* End of OptInput::OptInput */

OptInput::~OptInput( )
{

    /* Destructor */

} /* End of OptInput::~OptInput */

/* End of Input_Mod.cpp */

