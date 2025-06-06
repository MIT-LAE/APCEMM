#ifndef EPM_MODELS_ORIGINAL_HPP_INCLUDED
#define EPM_MODELS_ORIGINAL_HPP_INCLUDED

#include "Util/ForwardDecl.hpp"
#include "EPM/Models/Base.hpp"

namespace EPM::Models {

class Original : public Base {
public:
    Original(const OptInput &optInput, const Input &input,
             const Aircraft &aircraft, const Emission &EI,
             const Meteorology &met, const MPMSimVarsWrapper &simVars);

    std::variant<EPM::Output, SimStatus> run() override;

private:
    std::variant<EPM::Output, SimStatus>
    Integrate(const Vector_2D &aerArray);

    SimStatus RunMicrophysics(
        EPM::Output &out, const Vector_2D& aerArray, double delta_T);

    Vector_1D VAR_;  // Concentration of variable species.
    Vector_1D FIX_;  // Concentration of fixed species.
};

} // namespace EPM::Models

#endif
