#include <cctype>
#include <netcdf>

#include "Core/Input_Mod.hpp"
#include "EPM/EPM.hpp"

#include "EPM/Models/Original.hpp"
#include "EPM/Models/External.hpp"
#include "EPM/Models/NewPhysics.hpp"

namespace EPM {

struct {
    std::string name;
    std::string units;
    double Output::*value;
} scalar_vars[] = {{"finalTemp", "K", &Output::finalTemp},
                   {"iceRadius", "m", &Output::iceRadius},
                   {"iceDensity", "kg/m^3", &Output::iceDensity},
                   {"sootDensity", "kg/m^3", &Output::sootDensity},
                   {"H2O_mol", "molec/cm^3", &Output::H2O_mol},
                   {"SO4g_mol", "molec/cm^3", &Output::SO4g_mol},
                   {"SO4l_mol", "molec/cm^3", &Output::SO4l_mol},
                   {"area", "m^2", &Output::area},
                   {"bypassArea", "m^2", &Output::bypassArea},
                   {"coreExitTemp", "K", &Output::coreExitTemp}};

struct {
    std::string name;
    AIM::Aerosol Output::*value;
} pdf_vars[] = {
    {"so4", &Output::SO4Aer},
    {"ice", &Output::IceAer}};

void Output::write(std::string filename) const {
  NcFile nc(filename, NcFile::replace);

  // Scalar variables.

  for (auto &var : scalar_vars) {
    NcVar nc_var = nc.addVar(var.name, ncDouble);
    nc_var.putAtt("units", var.units);
    nc_var.putVar(&(this->*var.value));
  }

  // PDF variables: keep the SO4 and ice bin dimensions separate for
  // generality.

  for (auto &var : pdf_vars) {
    const AIM::Aerosol &aerosol = this->*var.value;

    const NcDim r_dim = nc.addDim(var.name + "_r", aerosol.getBinCenters().size());
    NcVar r_var = nc.addVar(var.name + "_r", ncDouble, r_dim);
    r_var.putAtt("units", "m");
    r_var.putAtt("long_name", var.name + " aerosol bin center radius");
    r_var.putVar(aerosol.getBinCenters().data());

    const NcDim r_e_dim = nc.addDim(var.name + "_r_e", aerosol.getBinEdges().size());
    NcVar r_e_var = nc.addVar(var.name + "_r_e", ncDouble, r_e_dim);
    r_e_var.putAtt("units", "m");
    r_e_var.putAtt("long_name", var.name + " aerosol bin edge radius");
    r_e_var.putVar(aerosol.getBinEdges().data());

    NcVar pdf_var = nc.addVar(var.name + "_pdf", ncDouble, r_dim);
    pdf_var.putAtt("long_name", var.name + " aerosol particle size distribution");
    pdf_var.putVar(aerosol.getPDF().data());
  }
}

void Output::read(std::string filename) {
    // Implementation for reading output from a file
}

std::unique_ptr<Models::Base> make_epm(
    const OptInput &optInput, const Input &input, const Aircraft &aircraft,
    const Emission &EI, const Meteorology &met,
    const MPMSimVarsWrapper &simVars) {
    switch (optInput.SIMULATION_EPM_TYPE) {
    case epm_type::EPM_ORIGINAL:
        return std::make_unique<Models::Original>(optInput, input, aircraft, EI, met, simVars);
    case epm_type::EPM_EXTERNAL:
      throw std::invalid_argument("EXTERNAL EPM NOT YET IMPLEMENTED!");
      // return std::make_unique<Models::External>(optInput, input, aircraft, EI, met, simVars);
    case epm_type::EPM_NEW_PHYSICS:
      throw std::invalid_argument("NEW PHYSICS EPM NOT YET IMPLEMENTED!");
      // return std::make_unique<Models::NewPhysics>(optInput, input, aircraft, EI, met, simVars);
    default:
        throw std::invalid_argument("Unknown EPM type specified in SIMULATION MENU.");
    }
}

}
