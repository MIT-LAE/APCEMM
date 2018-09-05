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


/* Mathematical constants */
#ifndef MATH_CONSTANTS
#define PI        3.141592653589793238460 /* \pi */
#endif

/* Physical constants */
#ifndef PHYS_CONSTANTS
#define Na        6.022140857E+23 /* Avogadro number   , Unit : [ 1 / mol ] */
#define kB        1.380648528E-23 /* Boltzmann constant, Unit : [ J /  K ] */
#define R         8.314459848E+00 /* Ideal gas constant, Unit : [ J / ( K.mol ) ] */

#define g         9.80665E+00 /* Acceleration due to gravity, Unit : [ m / s^2 ] */
#endif

/* Molar weights, Unit : [ kg / mol ] */
#ifndef MOLAR_WEIGHTS
#define mW_Air    28.5766E-03
#define mW_O2     31.9988E-03
#define mW_N2     28.0134E-03
#define mW_H2O    18.0153E-03
#define mW_HNO3   63.0128E-03
#define mW_SO4    98.0785E-03
#define mW_SO2    64.0638E-03
#endif
