#include "Util/PhysFunction.hpp"
#include "Core/MPMSimVarsWrapper.hpp"

MPMSimVarsWrapper::MPMSimVarsWrapper(const Input& input, const OptInput& Input_Opt, const double depth_estimate):
	RUN_BOXMODEL(Input_Opt.SIMULATION_BOXMODEL),
	BUILD_LUT(Input_Opt.SIMULATION_PARAMETER_SWEEP),
	SAVE_FORWARD(Input_Opt.SIMULATION_SAVE_FORWARD),
	ADJOINT(Input_Opt.SIMULATION_ADJOINT),
	BACKG_FILENAME(Input_Opt.SIMULATION_INPUT_BACKG_COND),
	THREADED_FFT(Input_Opt.SIMULATION_THREADED_FFT),
	USE_WISDOM(Input_Opt.SIMULATION_USE_FFTW_WISDOM),
	FFTW_DIR(Input_Opt.SIMULATION_DIRECTORY_W_WRITE_PERMISSION),
	TRANSPORT(Input_Opt.TRANSPORT_TRANSPORT),
	FILLNEG(Input_Opt.TRANSPORT_FILL),
	UPDRAFT(Input_Opt.TRANSPORT_UPDRAFT),
	UPDRAFT_TIME(Input_Opt.TRANSPORT_UPDRAFT_TIMESCALE),
	UPDRAFT_VEL(Input_Opt.TRANSPORT_UPDRAFT_VELOCITY),
	CHEMISTRY(Input_Opt.CHEMISTRY_CHEMISTRY),
	HETCHEM(Input_Opt.CHEMISTRY_HETCHEM),
	JRATE_FOLDER(Input_Opt.CHEMISTRY_JRATE_FOLDER),
	GRAVSETTLING(Input_Opt.AEROSOL_GRAVSETTLING),
	ICE_COAG(Input_Opt.AEROSOL_COAGULATION_SOLID),
	LIQ_COAG(Input_Opt.AEROSOL_COAGULATION_LIQUID),
	ICE_GROWTH(Input_Opt.AEROSOL_ICE_GROWTH),
	TEMP_PERTURB(Input_Opt.MET_ENABLE_TEMP_PERTURB),
	metDepth(depth_estimate),
	DIAG_FILENAME(Input_Opt.DIAG_FILENAME),
	TS_FOLDER(Input_Opt.SIMULATION_OUTPUT_FOLDER),
	TS_SPEC(Input_Opt.TS_SPEC),
	TS_SPEC_FILEPATH(TS_FOLDER + "/" + Input_Opt.TS_FILENAME),
    TS_SPEC_LIST(Input_Opt.TS_SPECIES),
	TS_FREQ(Input_Opt.TS_FREQ),
	TS_AERO(Input_Opt.TS_AERO),
	TS_AERO_FILEPATH(TS_FOLDER + "/" + Input_Opt.TS_AERO_FILENAME),
    TS_AERO_LIST(Input_Opt.TS_AEROSOL),
	TS_AERO_FREQ(Input_Opt.TS_AERO_FREQ),
	SAVE_PL(Input_Opt.PL_PL),
	SAVE_O3PL(Input_Opt.PL_O3),
	temperature_K(input.temperature_K()),
	pressure_Pa(input.pressure_Pa()),
	relHumidity_w(input.relHumidity_w())


{
    /* Compute relative humidity w.r.t ice */
    relHumidity_i = relHumidity_w * physFunc::pSat_H2Ol( temperature_K )\
                                             / physFunc::pSat_H2Os( temperature_K );

}

