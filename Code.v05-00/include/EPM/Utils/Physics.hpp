#ifndef EPM_UTILS_PHYSICS_H_INCLUDED
#define EPM_UTILS_PHYSICS_H_INCLUDED

namespace EPM::Utils::Physics {

/* Vortex sinking timescales, taken from Unterstrasser et al., 2008 */
const double t_Vortex_0 = 8.00E+00;
const double t_Vortex_1 = 1.10E+02;

/* Dilution timescales for a B747, taken from:
* B. Kärcher, "A trajectory box model for aircraft exhaust plumes", Journal of Geophysical Research, 1995 */
const double t_0 = 1.00E-04; /* [s],  */
const double t_j = 1.00E-02; /* [s],  */
const double t_1 = 8.00E+00; /* [s], Transition to vortex regime */
const double t_2 = 6.60E+01; /* [s], Transition to dispersion regime */

const double m = 2.0;
const double n = 50.0;
const double Cv = 3.0;

/* Engine exit plane characteristics for a B747, taken from:
* B. Kärcher, "A trajectory box model for aircraft exhaust plumes", Journal of Geophysical Research, 1995 */
/* Engine exit core area im m^2 */
const double Ac0 = 0.604;
/* Engine exit core velocity in m/s */
const double uc0 = 475.7;
/* Engine exit core temperature in K */
/* const double Tc0 = 547.3; */
/* double Tc0 */
/* Engine exit bypass area in m^2 */
/* const double Ab0 = 1.804; */
/* double Ab0 */

double entrainmentRate(const double time);
double depositionRate(const double r, const double T, const double P, const double H2O,
                      const double r_0,  const double theta );
double dT_Vortex(const double time, const double delta_T, bool deriv = 0);
bool isFreezable(const double r, const double T, const double H2O,
                 const double r0);
double condensationRate(const double r, const double T, const double P,
                        const double H2O, const double theta);

} // namespace EPM::Utils::Physics

#endif
