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

/* SYMMETRIES */
#define X_SYMMETRY              0    /* Is the problem symmetric
                                        around the x-axis? */
#define Y_SYMMETRY              0    /* Is the problem symmetric
                                        around the y-axis? */

/* OUTPUT CONCENTRATIONS */

/* Save output as double? otherwise save as float.
 * Saving as float will reduce the memory requirements */
#define SAVE_TO_DOUBLE          0


/* MASS CHECK */
#define NOy_MASS_CHECK          0    /* NOy mass check? */
#define CO2_MASS_CHECK          0    /* CO2 mass check? */
#define H2O_MASS_CHECK          0    /* H2O mass check? */

/* DEBUG */
/* DEBUG is now specified in Makefile header */
#define DEBUG_ALL               1    /* Debug all? */
//#define DEBUG                   0    /* Debug all except variables that have
//                                        specific debug options (see below) */

#define DEBUG_AC_INPUT          0    /* Debug AC Input? */
#define DEBUG_BG_INPUT          0    /* Debug Background Input? */
#define DEBUG_EI_INPUT          0    /* Debug Emission Input? */
#define DEBUG_RINGS             0    /* Debug Rings? */
#define DEBUG_MAPPING           0    /* Debug Mesh to ring mapping? */
#define DEBUG_COAGKERNEL        0    /* Debug Coagulation Kernel? */
#define DEBUG_ADJOINT           0    /* Debug adjoint? */

#if DEBUG_ALL

    #undef DEBUG
    #undef DEBUG_AC_INPUT
    #undef DEBUG_BG_INPUT
    #undef DEBUG_EI_INPUT
    #undef DEBUG_RINGS
    #undef DEBUG_MAPPING
    #undef DEBUG_COAGKERNEL
    #define DEBUG               1
    #define DEBUG_AC_INPUT      1
    #define DEBUG_BG_INPUT      1
    #define DEBUG_EI_INPUT      1
    #define DEBUG_RINGS         1
    #define DEBUG_MAPPING       1
    #define DEBUG_COAGKERNEL    1
    #define DEBUG_ADJOINT       1

#endif 

/* If adjoint is turned off, make sure DEBUG_ADJOINT is false */
#if ( !ADJOINT )

    #undef DEBUG_ADJOINT
    #define DEBUG_ADJOINT       0

#endif


#endif /* INTERFACE_H_INCLUDED */
