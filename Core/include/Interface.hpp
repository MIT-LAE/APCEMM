/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Interface Header File                                            */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 8/12/2018                                 */
/* File                 : Interface.hpp                             */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef INTERFACE_H_INCLUDED
#define INTERFACE_H_INCLUDED

/* TRANSPORT */
#define DIFFUSION               1    /* Is diffusion turned on? */
#define ADVECTION               1    /* Is advection turned on? */

/* FFTW */
#define FFTW_WISDOM             0    /* Find most efficient algorithm through FFTW_wisdom. Takes ~ 10s */
const char* const WISDOMFILE = "../SANDS/data/FFTW_Wisdom.out"; 

/* CHEMISTRY */
#define CHEMISTRY               1    /* Is chemistry turned on? */

/* MICROPHYSICS */
#define ICE_MICROPHYSICS        1    /* Is ice microphysics turned on? */
#define SUL_MICROPHYSICS        1    /* Is sulfate microphysics turned on? */

/* SYMMETRIES */
#define X_SYMMETRY              1    /* Is the problem symmetric around the x-axis? */
#define Y_SYMMETRY              1    /* Is the problem symmetric around the y-axis? */

/* BACKGROUND MIX RATIO */
const char* const AMBFILE    = "data/Ambient.txt";

/* OUTPUT */
#define DOSAVEPL                1    /* Save chemical rates? */

/* TIME */
#define TIME_IT                 1    /* Time simulation? */

/* MASS CHECK */
#define NOy_MASS_CHECK          1    /* NOy mass check? */

/* DEBUG */
#define DEBUG_AC_INPUT          0    /* Debug AC Input? */
#define DEBUG_BG_INPUT          0    /* Debug Background Input? */
#define DEBUG_EI_INPUT          0    /* Debug Emission Input? */
#define DEBUG_RINGS             0    /* Debug Rings? */
#define DEBUG_MAPPING           0    /* Debug Mesh to ring mapping? */

#endif /* INTERFACE_H_INCLUDED */
