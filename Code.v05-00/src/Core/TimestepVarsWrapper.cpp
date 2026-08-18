#include "Core/TimestepVarsWrapper.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
//Declaration of vars is same order as in header, so no undefined behavior
TimestepVarsWrapper::TimestepVarsWrapper(const Input& input, const OptInput& Input_Opt):
tEmission_h(input.emissionTime()),
tInitial_h(tEmission_h),
tFinal_h(tInitial_h + input.simulationTime()),
tInitial_s(tInitial_h * 3600.0),
tFinal_s(tFinal_h   * 3600.0),
curr_Time_s(tInitial_s),
dt(0),
nTime(0),
LAST_STEP(0),
ITS_TIME_FOR_TRANSPORT(0),
ITS_TIME_FOR_LIQ_COAGULATION(0),
ITS_TIME_FOR_ICE_COAGULATION(0),
ITS_TIME_FOR_ICE_GROWTH(0),
ITS_TIME_FOR_TEMPPERTURB(0),
TRANSPORT_DT(Input_Opt.TRANSPORT_TIMESTEP * 60.0),
CHEMISTRY_DT(Input_Opt.CHEMISTRY_TIMESTEP * 60.0),
COAG_DT(Input_Opt.AEROSOL_COAGULATION_TIMESTEP * 60.0),
TEMP_PERTURB_DT(Input_Opt.MET_TEMP_PERTURB_TIMESCALE * 60.0),
ICE_GROWTH_DT(Input_Opt.AEROSOL_ICE_GROWTH_TIMESTEP * 60.0),
lastTimeTransport(curr_Time_s),
lastTimeChem(curr_Time_s),
lastTimeLiqCoag(curr_Time_s),
lastTimeIceCoag(curr_Time_s),
lastTimeIceGrowth(curr_Time_s),
lastTimeTempPerturb(curr_Time_s),
totalIceParticles_before(0),
totalIceMass_before(0),
totalIceParticles_initial(0),
totalIceMass_initial(0),
totalIceParticles_now(0),
totalIceMass_now(0),
totalIceParticles_last(0),
totalIceMass_last(0),
totalIceParticles_after(0),
totalIceMass_after(0),
totPart_lost(0),
totIce_lost(0),
ABORT_THRESHOLD(1.0e-3)
{
    /* The base (heartbeat) timestep is determined from enabled transport and ice growth processes. */
    
    Vector_1D timesteps;

    if(Input_Opt.TRANSPORT_TRANSPORT) 
        timesteps.push_back(TRANSPORT_DT);
    if(Input_Opt.AEROSOL_ICE_GROWTH)
        timesteps.push_back(ICE_GROWTH_DT);
    
    if (timesteps.empty()) {
        throw std::runtime_error("No active processes to determine timestep");
    }

    dt = *(std::min_element(timesteps.begin(), timesteps.end()));
    std::cout << "Calculated Timestep: " << dt/60.0 << "[min]" << std::endl;
    if (dt <= 0) throw std::runtime_error("Invalid Timestep"); 

    for (const double &step : timesteps) {
        double ratio = step / dt;
        if (std::abs(ratio - std::round(ratio)) > 1.0e-6) {
            throw std::runtime_error("Process timesteps must be integer multiples of the heartbeat timestep");
        }
    }
}