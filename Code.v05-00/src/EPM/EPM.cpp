#include "Core/Input_Mod.hpp"
#include "EPM/EPM.hpp"

#include "EPM/Models/Original.hpp"
#include "EPM/Models/External.hpp"
#include "EPM/Models/NewPhysics.hpp"

namespace EPM {

std::unique_ptr<Models::Base> make_epm(
    const OptInput &optInput, const Input &input, const Aircraft &aircraft,
    const Emission &EI, const Meteorology &met,
    const MPMSimVarsWrapper &simVars) {
    switch (optInput.SIMULATION_EPM_TYPE) {
    case epm_type::EPM_ORIGINAL:
        return std::make_unique<Models::Original>(optInput, input, aircraft, EI, met, simVars);
    case epm_type::EPM_EXTERNAL:
        return std::make_unique<Models::External>(optInput, input, aircraft, EI, met, simVars);
    case epm_type::EPM_NEW_PHYSICS:
        return std::make_unique<Models::NewPhysics>(optInput, input, aircraft, EI, met, simVars);
    default:
        throw std::invalid_argument("Unknown EPM type specified in SIMULATION MENU.");
    }
}

}
