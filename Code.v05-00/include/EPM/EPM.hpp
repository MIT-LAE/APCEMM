#ifndef EPM_HPP_INCLUDED
#define EPM_HPP_INCLUDED

#include <memory>

#include "AIM/Aerosol.hpp"
#include "Core/Aircraft.hpp"
#include "Core/Emission.hpp"
#include "Core/Input.hpp"
#include "Core/Input_Mod.hpp"
#include "Core/MPMSimVarsWrapper.hpp"
#include "Core/Meteorology.hpp"


namespace EPM {
    struct Output {
    double finalTemp;
    double iceRadius;
    double iceDensity;
    double sootDensity;
    double H2O_mol;
    double SO4g_mol;
    double SO4l_mol;
    AIM::Aerosol SO4Aer;
    AIM::Aerosol IceAer;
    double area;
    double bypassArea;
    double coreExitTemp;
    };

    namespace Models {
        class Base;
        class Original;
        class NewPhysics;
        class External;
    }

    std::unique_ptr<Models::Base> make_epm(
        const OptInput &optInput, const Input &input,
        const Aircraft &aircraft, const Emission &EI,
        const Meteorology &met, const MPMSimVarsWrapper &simVars);

}

#endif
