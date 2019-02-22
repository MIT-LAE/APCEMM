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

/* APCEMM Look-up table */
#define APCEMM_LUT              1    /* Build look-up table? */
#define REBUILD                 0    /* Overwrite output files */

#define PSC_SIM                 0    /* Polar Stratospheric Clouds? */

/* MICROPHYSICS */
#define ICE_MICROPHYSICS        0    /* Is ice microphysics turned 
                                        on? */
#define ICECOAG_TSTEP           3600 /* Minimal coagulation time step 
                                        for ice in s */
#define LIQ_MICROPHYSICS        0    /* Is sulfate microphysics turned
                                        on? */
#define LIQCOAG_TSTEP           3600 /* Minimal coagulation time step for 
                                        liquid aerosols in s */

/* SYMMETRIES */
#define X_SYMMETRY              1    /* Is the problem symmetric 
                                        around the x-axis? */
#define Y_SYMMETRY              1    /* Is the problem symmetric 
                                        around the y-axis? */

/* METEOROLOGICAL DATA */ 
#define LOAD_MET                0    /* Load fine plume-scale met 
                                        data (2D water and 
                                        tempereature fields) */

/* OUTPUT CONCENTRATIONS */

/* Save output as double? otherwise save as float.
 * Saving as float will reduce the memory requirements */
#define SAVE_TO_DOUBLE          1 


/* OUTPUT MICROPHYSICS */
/* Save solid aerosol gridded bins? */
/* At each time step, the array is NX*NY*nBin_PA */
// This has been removed and can now be handled using timeseries output
//#define SAVE_PA_MICROPHYS       0
/* Time interval at which solid aerosol properties are saved 
 * (in seconds) */
//#define SAVE_PA_DT              3600
//
//const char* const OUT_FILE_PA = "data/IceAerosol.nc";

/* If ice microphysics is not performed, make sure that we do not 
 * output it */
//#if ( !ICE_MICROPHYSICS )
//
//    #undef SAVE_PA_MICROPHYS 
//    #define SAVE_PA_MICROPHYS  0
//
//#endif

/* Save liquid aerosol gridded bins? */
/* At each time step, the array is NX*NY*nBin_LA */
// This has been removed and can now be handled using timeseries output
//#define SAVE_LA_MICROPHYS       0
/* Time interval at which liquid aerosol properties are saved 
 * (in seconds) */
//#define SAVE_LA_DT              3600
//
//const char* const OUT_FILE_LA = "data/LiqAerosol.nc";

/* If liquid microphysics is not performed, make sure that we do not 
 * output it */
//#if ( !LIQ_MICROPHYSICS )
//
//    #undef SAVE_LA_MICROPHYS
//    #define SAVE_LA_MICROPHYS   0
//
//#endif


/* MASS CHECK */
#define NOy_MASS_CHECK          0    /* NOy mass check? */
#define CO2_MASS_CHECK          0    /* CO2 mass check? */
#define H2O_MASS_CHECK          0    /* H2O mass check? */

/* DEBUG */
/* DEBUG is now specified in Makefile header */
#define DEBUG_ALL               0    /* Debug all? */
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
