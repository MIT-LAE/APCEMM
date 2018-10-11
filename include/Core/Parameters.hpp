/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Parameter Header File                                            */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : Parameters.hpp                            */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef PARAMETERS_H_INCLUDED
#define PARAMETERS_H_INCLUDED

/* Grid parameters */
#define NX                256         /* Number of grid cells in the x-direction */
#define NY                256         /* Number of grid cells in the y-direction */
#define NCELL             NX*NY       /* Number of grid cells */
#define XLIM              6.5E+03     /* x-limits of the domain [m] */
#define YLIM              6.5E+02     /* y-limits of the domain [m] */

/* Ring structure */
#define RINGS             1           /* 1 or 0? 1: Ring structure, 0: No rings */
#define NRING             15          /* Number of rings */

/* Time */
#define TSTART            0.0E+00     /* Initial time [hrs] (after emission) */
#define TSIMUL            2.4E+01     /* Simulation time [hrs] */
#define DT                9.0E+02     /* Default time step [s] */

/* Atmospheric parameters */
#define GAMMA            -3.0E+00     /* Ambient temperature lapse rate [K/km] */
#define HOMOGENEOUS_NAT   1           /* Allow homogeneous NAT? */
#define T_NAT_SUPERCOOL   3.0         /* NAT supercooling requirement [K] */

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
                                       * 2: Exponential
                                       */

/* Advection  */
#define VX                0.0E+00     /* Steady-state horizontal advection velocity [m/s] */
#define VY                0.0E+00     /* Steady-state vertical advection velocity [m/s] */
#define SYNLIFT           0           /* Initial synoptic lifting? */
#define T_SYN             3.6E+03     /* Timescale associated to synoptic lifting */
#define V_SYN             5.0E-02     /* Velocity magnitude [m/s] */
#define SYNPROF           1           /* Time profile of the initial synoptic lifting [0,1]
                                       * 0: Step
                                       * 1: Linear
                                       */

/* Chemistry parameters */
#define N_SPC                 135         /* Number of chemical species */
#define N_VAR                 127         /* Number of Variable species */
#define N_FIX                 8           /* Number of Fixed species */
#define N_REACT               475         /* Number of reactions */
#define KPP_RTOLS             1.00E-03    /* Relative tolerances in KPP */
#define KPP_ATOLS             1.00E-03    /* Absolute tolerances in KPP */


/* Microphysics parameters */
#define SO4_R_LOW             1.00E-10    /* Sulfates' lower bin radius [m] */
#define SO4_R_HIG             5.00E-07    /* Sulfates' larger bin radius [m] */
#define SO4_VRAT              1.50E+00    /* Size ratio between two consecutives bins */

/* Early plume integration */
#define VORTEX_SINKING        1           /* Consider vortex sinking? */
#define EPM_RTOLS             1.00E-05    /* Relative tolerances in EPM */
#define EPM_ATOLS             1.00E-07    /* Absolute tolerances in EPM */

#endif /* PARAMETERS_H_INCLUDED */
