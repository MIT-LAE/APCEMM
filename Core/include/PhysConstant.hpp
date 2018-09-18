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
#define R_Air     R/mW_Air        /* Specific gas constant, Unit : [ J / (K.kg) ] */

#define gamma_Air 1.4             /* Heat capacity ratio of air, Unit : [ - ] */

#define g         9.80665E+00     /* Acceleration due to gravity, Unit : [ m / s^2 ] */

#define PRES_SL   101325.0        /* Pressure at sea level, Unit : [ Pa ] */
#define TEMP_SL   288.15          /* Temperature at sea leve, Unit : [ K ] */

#endif /* PHYS_CONSTANTS */

#endif /* PHYSCONSTANT_H_INCLUDED */
