#include "EPM/Models/External.hpp"

namespace EPM::Models {

External::External(const OptInput &optInput, const Input &input,
                   const Aircraft &aircraft, const Emission &EI,
                   const Meteorology &met, const MPMSimVarsWrapper &simVars) :
    Base(optInput, input, aircraft, EI, met, simVars)
{
    if (optInput.SIMULATION_EXTERNAL_EPM_NETCDF_FILENAME == "=MISSING=") {
        throw std::invalid_argument(
            "External EPM NetCDF filename is missing in the input file.");
    }
}

std::variant<EPM::Output, SimStatus> External::run() {
    EPM::Output output;
    std::cout << "Reading external EPM NetCDF file: "
              << optInput_.SIMULATION_EXTERNAL_EPM_NETCDF_FILENAME << std::endl;
    output.read(optInput_.SIMULATION_EXTERNAL_EPM_NETCDF_FILENAME);
    return output;
}

}
