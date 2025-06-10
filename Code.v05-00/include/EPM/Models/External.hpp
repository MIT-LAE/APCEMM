#ifndef EPM_MODELS_EXTERNAL_HPP_INCLUDED
#define EPM_MODELS_EXTERNAL_HPP_INCLUDED

#include "EPM/Models/Base.hpp"


namespace EPM::Models {

class External : public Base {
public:
    External(const OptInput &optInput, const Input &input,
             const Aircraft &aircraft, const Emission &EI,
             const Meteorology &met, const MPMSimVarsWrapper &simVars):
    Base(optInput, input, aircraft, EI, met, simVars) {
      throw std::runtime_error("External EPM model is not implemented yet.");
    }

    std::variant<EPM::Output, SimStatus> run() override;
};

}

#endif
