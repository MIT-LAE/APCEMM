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
    SIMULATION_THREADED_FFT( 0 ),
    SIMULATION_USE_FFTW_WISDOM( 0 ),
    SIMULATION_DIRECTORY_W_WRITE_PERMISSION( "" ),
    SIMULATION_INPUT_BACKG_COND( "" ),
    SIMULATION_SAVE_FORWARD( 0 ),
    SIMULATION_FORWARD_FILENAME( "" ),
    SIMULATION_ADJOINT( 0 ),
    SIMULATION_ADJOINT_FILENAME( "" ),
    SIMULATION_BOXMODEL( 0 ),
    SIMULATION_BOX_FILENAME( "" ),
    PARAMETER_FILEINPUT( 0 ),
    PARAMETER_FILENAME( "" ),
    PARAMETER_PARAM_MAP(std::unordered_map<std::string, Vector_1D>()),
    TRANSPORT_TRANSPORT( 0 ),
    TRANSPORT_FILL( 0 ),
    TRANSPORT_TIMESTEP( 0.0E+00 ),
    TRANSPORT_PART_FLUX( 0 ),
    TRANSPORT_UPDRAFT( 0 ),
    TRANSPORT_UPDRAFT_TIMESCALE( 0.0E+00 ),
    TRANSPORT_UPDRAFT_VELOCITY( 0.0E+00 ),
    CHEMISTRY_CHEMISTRY( 0 ),
    CHEMISTRY_HETCHEM( 0 ),
    CHEMISTRY_JRATE_FOLDER( "" ),
    CHEMISTRY_TIMESTEP( 0.0E+00 ),
    AEROSOL_GRAVSETTLING( 0 ),
    AEROSOL_COAGULATION_SOLID( 0 ),
    AEROSOL_COAGULATION_LIQUID( 0 ),
    AEROSOL_COAGULATION_TIMESTEP( 0.0E+00 ),
    AEROSOL_ICE_GROWTH( 0 ),
    MET_LOADMET( 0 ),
    MET_FILENAME( "" ),
    MET_DT (0),
    MET_LOADTEMP( 0 ),
    MET_TEMPTIMESERIES( 0 ),
    MET_LOADRH( 0 ),
    MET_RHTIMESERIES( 0 ),
    MET_LOADSHEAR( 0 ),
    MET_SHEARTIMESERIES( 0 ),
    MET_FIXDEPTH( 0 ),
    MET_DEPTH( 0.0E+00 ),
    MET_FIXLAPSERATE( 0 ),
    MET_LAPSERATE( 0.0E+00 ),
    MET_DIURNAL( 0 ),
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