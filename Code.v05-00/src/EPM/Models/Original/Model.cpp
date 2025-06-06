#include <variant>

#include "EPM/Models/Original.hpp"
#include "EPM/EPM.hpp"
#include "EPM/Solution.hpp"
#include "KPP/KPP_Parameters.h"
#include "Util/PhysConstant.hpp"

using physConst::kB;


namespace EPM::Models {

Original::Original(const OptInput &optInput, const Input &input, const Aircraft &aircraft,
                    const Emission &EI, const Meteorology &met,
                    const MPMSimVarsWrapper &simVars) :
    Base(optInput, input, aircraft, EI, met, simVars), VAR_(NVAR)
{ }

std::variant<EPM::Output, SimStatus> Original::run() {
    // This is used for initializing the EPM solution. The only thing that
    // gets taken from it in the end is the soot density.
    EPM::Solution epmSolution(optInput_);

    /* Compute airDens from pressure and temperature */
    double airDens = simVars_.pressure_Pa / (kB    * met_.tempRef()) * 1.00E-06;
    // [molec/cm3] = [Pa = J/m3]          / ([J/K] * [K])            * [m3/cm3]

    /* This sets the species array values in the Solution data structure to the
    ambient data file, EXCEPT FOR H2O which is user defined via met input or
    rhw input */
    epmSolution.Initialize(
        simVars_.BACKG_FILENAME.c_str(), input_, airDens, met_,
        optInput_, VAR_, false);

    double dy =
        (optInput_.ADV_GRID_YLIM_UP + optInput_.ADV_GRID_YLIM_DOWN) /
        optInput_.ADV_GRID_NY;
    double y0 = -optInput_.ADV_GRID_YLIM_DOWN;
    unsigned int i_0 = static_cast<unsigned int>(
        std::floor(optInput_.ADV_GRID_XLIM_LEFT /
                   optInput_.ADV_GRID_NX)); // index i where x = 0
    unsigned int j_0 = static_cast<unsigned int>(std::floor(-y0 / dy)); // index j where y = 0

    // This sets the values in VAR_ to the values in the solution data
    // structure at indices i, j
    epmSolution.getData(VAR_, i_0, j_0);

    // Run EPM.
    return Integrate(epmSolution.getSootDensity());
}

} // namespace EPM::Models
