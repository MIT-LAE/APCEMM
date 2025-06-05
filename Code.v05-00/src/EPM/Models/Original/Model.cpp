#include <variant>

#include "EPM/Models/Original.hpp"
#include "EPM/Models/Original/Integrate.hpp"
#include "EPM/EPM.hpp"
#include "EPM/Solution.hpp"
#include "KPP/KPP_Parameters.h"
#include "Util/PhysConstant.hpp"

using physConst::kB;


namespace EPM::Models {

    Original::Original(const OptInput &optInput, const Input &input, const Aircraft &aircraft,
                       const Emission &EI, const Meteorology &met,
                       const MPMSimVarsWrapper &simVars) :
        Base(optInput, input, aircraft, EI, met, simVars) { }

  std::variant<EPM::Output, SimStatus> Original::run() {
    double C[NSPEC];        /* Concentration of all species */
    double *VAR = &C[0];    /* Concentration of variable species */
    double *FIX = &C[NVAR]; /* Concentration of fixed species */

    // Still need the solution data structure...
    // TODO: WHAT DOES THAT MEAN? WHY "STILL NEED"?
    EPM::Solution epmSolution(optInput_);

    /* Compute airDens from pressure and temperature */
    double airDens = simVars_.pressure_Pa / (kB * met_.tempRef()) * 1.00E-06;
    // [molec/cm3] = [Pa = J/m3]          / ([J/K] * [K])            * [m3/cm3]

    /* This sets the species array values in the Solution data structure to the
    ambient data file, EXCEPT FOR H2O which is user defined via met input or
    rhw input */
    epmSolution.Initialize(
        simVars_.BACKG_FILENAME.c_str(), input_, airDens, met_,
        optInput_, VAR, FIX, false);

    Vector_2D aerArray = epmSolution.getAerosol();

    Vector_1D yEdges(optInput_.ADV_GRID_NY + 1);
    double dy =
        (optInput_.ADV_GRID_YLIM_UP + optInput_.ADV_GRID_YLIM_DOWN) /
        optInput_.ADV_GRID_NY;
    int i_0 = std::floor(optInput_.ADV_GRID_XLIM_LEFT /
                        optInput_.ADV_GRID_NX); // index i where x = 0
    int j_0 = std::floor(-yEdges[0] / dy);      // index j where y = 0

    // This sets the values in VAR and FIX to the values in the solution data
    // structure at indices i, j
    epmSolution.getData(VAR, FIX, i_0, j_0);

    // Run EPM.
    return OriginalImpl::Integrate(
            met_.tempRef(), simVars_.pressure_Pa, met_.rhwRef(),
            input_.bypassArea(), input_.coreExitTemp(), VAR, aerArray, aircraft_,
            EI_, simVars_.CHEMISTRY, optInput_.ADV_AMBIENT_LAPSERATE,
            input_.fileName_micro());
}

} // namespace EPM::Models
