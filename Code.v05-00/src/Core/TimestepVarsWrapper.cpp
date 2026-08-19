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
ITS_TIME_FOR_ICE_GROWTH(0),
TRANSPORT_DT(Input_Opt.TRANSPORT_TIMESTEP * 60.0),
ICE_GROWTH_DT(Input_Opt.AEROSOL_ICE_GROWTH_TIMESTEP * 60.0),
lastTimeTransport(curr_Time_s),
lastTimeIceGrowth(curr_Time_s),
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
    dt = TRANSPORT_DT;
    if (dt <= 0) throw std::runtime_error("Invalid Timestep"); 
    std::cout << "Outer Timestep: " << dt/60.0 << " [min]" << std::endl;
    std::cout << "Inner Physics Substep: " << Input_Opt.TRANSPORT_ICE_GROWTH_SUBSTEP << " [s]" << std::endl;
}