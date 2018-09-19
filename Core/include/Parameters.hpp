/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Parameter Header File                                            */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : Parameters.h                              */
/* Working directory    : /home/fritzt/APCEMM-SourceCode            */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef PARAMETERS_H_INCLUDED
#define PARAMETERS_H_INCLUDED

/* Grid parameters */
#define NX                256         /* Number of grid cells in the x-direction */
#define NY                256         /* Number of grid cells in the y-direction */
#define NCELL             NX*NY       /* Number of grid cells */
#define XLIM              5.0E+03     /* x-limits of the domain [m] */
#define YLIM              7.0E+02     /* y-limits of the domain [m] */

/* Ring structure */
#define RINGS             1           /* 1 or 0? 1: Ring structure, 0: No rings */
#define NRING             15          /* Number of rings */

/* Time */
#define TSTART            0.0E+00     /* Initial time [hrs] (after emission) */
#define TSIMUL            2.4E+01     /* Simulation time [hrs] */
#define DT                9.0E+02     /* Default time step [s] */

/* Diffusion */
#define DH                1.5E+01     /* Steady-state horizontal diffusion parameter [m2/s] */
#define DV                1.5E-01     /* Steady-state vertical diffusion parameter [m2/s] */
#define DH0               1.7E+01     /* Initial horizontal diffusion parameter [m2/s] */
#define DV0               2.0E-01     /* Initial vertical diffusion parameter [m2/s] */
#define tH0               1.8E+02     /* Timescale of initial enhanced horizontal diffusion [s] */
#define tV0               1.8E+02     /* Timescale of initial enhanced vertical diffusion [s] */
#define DPROF             1           /* Time profile of the initial diffusion parameter [0,1,2]
                                       * 0: Step
                                       * 1: Linear
                                       * 2: Exponential */

/* Chemistry parameters */
#define N_SPC                 135         /* Number of chemical species */
#define N_VAR                 127         /* Number of Variable species */
#define N_FIX                 8           /* Number of Fixed species */
#define N_REACT               475         /* Number of reactions */
#define RTOLS                 1.00E-03    /* Relative tolerances in KPP */
#define ATOLS                 1.00E-03    /* Absolute tolerances in KPP */

#endif /* PARAMETERS_H_INCLUDED */
