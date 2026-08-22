#include "Util/PhysFunction.hpp"
#include "Core/MPMSimVarsWrapper.hpp"

MPMSimVarsWrapper::MPMSimVarsWrapper(const Input& input, const OptInput& Input_Opt, const double depth_estimate):
	// Placeholder for a future look-up table build. Nothing reads it: it used
	// to mirror the parameter sweep flag, which no longer exists.
	BUILD_LUT(false),
	BACKG_FILENAME(Input_Opt.SIMULATION_INPUT_BACKG_COND),
	TRANSPORT(Input_Opt.TRANSPORT_TRANSPORT),
	GRAVSETTLING(Input_Opt.AEROSOL_GRAVSETTLING),
	ICE_COAG(Input_Opt.AEROSOL_COAGULATION_SOLID),
	ICE_GROWTH(Input_Opt.AEROSOL_ICE_GROWTH),
	TEMP_PERTURB(Input_Opt.MET_ENABLE_TEMP_PERTURB),
	metDepth(depth_estimate),
	TS_FOLDER(Input_Opt.SIMULATION_OUTPUT_FOLDER),
	TS_AERO(Input_Opt.TS_AERO),
	TS_AERO_FILEPATH(TS_FOLDER + "/" + Input_Opt.TS_AERO_FILENAME),
	TS_AERO_FREQ(Input_Opt.TS_AERO_FREQ),
	temperature_K(input.temperature_K()),
	pressure_Pa(input.pressure_Pa()),
	relHumidity_w(input.relHumidity_w())


{
    /* Compute relative humidity w.r.t ice */
    relHumidity_i = relHumidity_w * physFunc::pSat_H2Ol( temperature_K )\
                                             / physFunc::pSat_H2Os( temperature_K );

}

