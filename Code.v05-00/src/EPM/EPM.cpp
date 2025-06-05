#include "Core/Input_Mod.hpp"
#include "EPM/EPM.hpp"

#include "EPM/Models/Original.hpp"
#include "EPM/Models/External.hpp"
#include "EPM/Models/NewPhysics.hpp"


std::unique_ptr<EPM::Models::Base> EPM::make_epm(const OptInput &optInput) {
  switch (optInput.SIMULATION_EPM_TYPE) {
  case epm_type::EPM_ORIGINAL:
    return std::make_unique<EPM::Models::Original>(optInput);
  case epm_type::EPM_EXTERNAL:
    return std::make_unique<EPM::Models::External>(optInput);
  case epm_type::EPM_NEW:
    return std::make_unique<EPM::Models::NewPhysics>(optInput);
  }
}
