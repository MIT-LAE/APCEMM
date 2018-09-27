/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* PhysConstant Header File                                         */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : PhysConstant.h                            */
/* Working directory    : /home/fritzt/APCEMM-SourceCode            */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


#ifndef PHYSCONSTANT_H_INCLUDED
#define PHYSCONSTANT_H_INCLUDED

#include "MolarWeights.hpp"

/* Physical constants */
#ifndef PHYS_CONSTANTS
#define PHYS_CONSTANTS

#define Na        6.022140857E+23 /* Avogadro number   , Unit : [ 1 / mol ] */
#define kB        1.380648528E-23 /* Boltzmann constant, Unit : [ J /  K ] */
#define R         8.314459848E+00 /* Ideal gas constant, Unit : [ J / ( K.mol ) ] */
#define R_Air     R/MW_Air        /* Specific gas constant, Unit : [ J / (K.kg) ] */

#define CP_Air    1.005           /* Specific heat capacity of air at 300K, Unit : [ kJ / ( kg K ) ] */
#define GAMMA_Air 1.4             /* Heat capacity ratio of air, Unit : [ - ] */

#define g         9.80665E+00     /* Acceleration due to gravity, Unit : [ m / s^2 ] */

#define PRES_SL   101325.0        /* Pressure at sea level, Unit : [ Pa ] */
#define TEMP_SL   288.15          /* Temperature at sea level, Unit : [ K ] */

#define GAMMA_AD  -g/CP_Air       /* Adiabatic temperature lapse rate, Unit : [ K / km ] */

#define EPSILON   1.00E-05        /* Turbulent dissipation rate, 
                                   * Range : ( 1.0E-08 - 1.0E-02 ), 
                                   * Unit  : [ m ^ 2 / s ^ 3 ]
                                   * Source: U. Schumann, A contrail cirrus prediction model, Geo. Sc. Dev., 2012 */

#endif /* PHYS_CONSTANTS */

#endif /* PHYSCONSTANT_H_INCLUDED */
