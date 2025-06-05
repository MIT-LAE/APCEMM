#ifndef EPM_HPP_INCLUDED
#define EPM_HPP_INCLUDED

#include <memory>

#include "Core/Input_Mod.hpp"
#include "AIM/Aerosol.hpp"
#include "EPM/Models/Base.hpp"

namespace EPM {
    std::unique_ptr<Models::Base> make_epm(const OptInput &optInput);

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
}

#endif
