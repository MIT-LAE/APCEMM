#ifndef EPM_MODELS_BASE_HPP_INCLUDED
#define EPM_MODELS_BASE_HPP_INCLUDED

#include <variant>

#include "Core/Input.hpp"
#include "Core/Input_Mod.hpp"
#include "Core/Status.hpp"
#include "EPM/EPM.hpp"


namespace EPM::Models {

class Base {
public:
    Base() = delete;
    Base(const OptInput &optInput, const Input &input,
         const Aircraft &aircraft, const Emission &EI,
         const Meteorology &met, const MPMSimVarsWrapper &simVars):
        optInput_(optInput), input_(input),
        aircraft_(aircraft), EI_(EI), met_(met), simVars_(simVars) {}
    virtual ~Base() = default;

    // Disable copy and move semantics: EPM model instances will only ever be
    // accessed via std::unique_ptr values.
    Base(const Base &) = delete;
    Base &operator=(const Base &) = delete;
    Base(Base &&) = delete;
    Base &operator=(Base &&) = delete;

    // Virtual function to be overridden by derived classes
    virtual std::variant<EPM::Output, SimStatus> run() = 0;

protected:
    const OptInput &optInput_;
    const Input &input_;
    const Aircraft &aircraft_;
    const Emission &EI_;
    const Meteorology &met_;
    const MPMSimVarsWrapper &simVars_;
};

} // namespace EPM::Models

#endif
