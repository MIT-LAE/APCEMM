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

/* APCEMM LUT */
#define APCEMM_LUT              1    /* Build look-up table? */

/* TRANSPORT */
#define DIFFUSION               1    /* Is diffusion turned on? */
#define ADVECTION               1    /* Is advection turned on? */

/* FFTW */
#define FFTW_WISDOM             0    /* Find most efficient algorithm through FFTW_wisdom. Takes ~ 10s */
const char* const WISDOMFILE = "data/FFTW_Wisdom.out"; 

/* CHEMISTRY */
#define CHEMISTRY               1    /* Is chemistry turned on? */
#define HETCHEMISTRY            1    /* Is heterogeneous chemistry turned on? */
#define PSC_SIM                 0    /* Polar Stratospheric Clouds? */

/* MICROPHYSICS */
#define ICE_MICROPHYSICS        1    /* Is ice microphysics turned on? */
#define ICECOAG_TSTEP           300  /* Minimal coagulation time step for ice in s */
#define LIQ_MICROPHYSICS        1    /* Is sulfate microphysics turned on? */
#define LIQCOAG_TSTEP           3600 /* Minimal coagulation time step for liquid aerosols in s */

/* SYMMETRIES */
#define X_SYMMETRY              1    /* Is the problem symmetric around the x-axis? */
#define Y_SYMMETRY              1    /* Is the problem symmetric around the y-axis? */

/* BACKGROUND MIX RATIO */
const char* const AMBFILE     = "data/Ambient.txt";

/* OUTPUT */
#define SAVE_OUTPUT             1    /* Save output? */
const char* const OUT_FILE    = "data/output.nc";
#define SAVE_TO_DOUBLE          1    /* Save output as double? otherwise float */
#define SAVEPL                  1    /* Save chemical rates */
#define SAVE_PA_MICROPHYS       0    /* Save solid aerosol gridded bins? ( Size = NX*NY*nBin_PA ) */
#define SAVE_PA_DT              600  /* Output solid aerosol every x seconds */
#define SAVE_LA_MICROPHYS       1    /* Save liquid aerosol gridded bins? ( Size = NX*NY*nBin_LA ) */
#define SAVE_LA_DT              3600 /* Output liquid aerosol every x seconds */
const char* const OUT_FILE_MICROPHYS = "data/microphys.nc";

/* TIME */
#define TIME_IT                 0    /* Time simulation? */

/* MASS CHECK */
#define NOy_MASS_CHECK          1    /* NOy mass check? */
#define CO2_MASS_CHECK          0    /* CO2 mass check? */
#define H2O_MASS_CHECK          1    /* H2O mass check? */

/* DEBUG */
#define DEBUG_AC_INPUT          0    /* Debug AC Input? */
#define DEBUG_BG_INPUT          0    /* Debug Background Input? */
#define DEBUG_EI_INPUT          0    /* Debug Emission Input? */
#define DEBUG_RINGS             0    /* Debug Rings? */
#define DEBUG_MAPPING           0    /* Debug Mesh to ring mapping? */
#define DEBUG_COAGKERNEL        0    /* Debug Coagulation Kernel? */
#define PRINT_DEBUG             0    /* Debug? */

#if ( PRINT_DEBUG )
#define DEBUG(X)         std::cout << x
#else
#define DEBUG(X)
#endif

#endif /* INTERFACE_H_INCLUDED */
