#ifndef EPM_MODELS_NEW_PHYSICS_HPP_INCLUDED
#define EPM_MODELS_NEW_PHYSICS_HPP_INCLUDED

#include <stdexcept>

#include "EPM/Models/Base.hpp"


namespace EPM::Models {

class NewPhysics : public Base {
public:
    NewPhysics(const OptInput &optInput, const Input &input,
               const Aircraft &aircraft, const Emission &EI,
               const Meteorology &met, const MPMSimVarsWrapper &simVars) :
        Base(optInput, input, aircraft, EI, met, simVars) {
        throw std::runtime_error("New EPM model is not implemented yet.");
    }

    std::variant<EPM::Output, SimStatus> run() override;
};

} // namespace EPM::Models

#endif
