/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/* Main Program File                                                */
/*                                                                  */
/* Generated by KPP-2.2.3 symbolic chemistry Kinetics PreProcessor  */
/*       (http://www.cs.vt.edu/~asandu/Software/KPP)                */
/* KPP is distributed under GPL, the general public licence         */
/*       (http://www.gnu.org/copyleft/gpl.html)                     */
/* (C) 1995-1997, V. Damian & A. Sandu, CGRER, Univ. Iowa           */
/* (C) 1997-2005, A. Sandu, Michigan Tech, Virginia Tech            */
/*     With important contributions from:                           */
/*        M. Damian, Villanova University, USA                      */
/*        R. Sander, Max-Planck Institute for Chemistry, Mainz, Germany */
/*                                                                  */
/* File                 : KPP_Main.c                                */
/* Time                 : Fri Jul 13 11:26:33 2018                  */
/* Equation file        : KPP.kpp                                   */
/* Output root filename : KPP                                       */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#ifdef OMP
    #include "omp.h"
#endif /* OMP */

#include "KPP/KPP.hpp"
#include "KPP/KPP_Parameters.h"
#include "KPP/KPP_Global.h"
#include "KPP/KPP_Sparse.h"

#define MAX(a, b) \
    ({ __typeof__ (a) _a = (a); \
        __typeof__ (b) _b = (b); \
      _a > _b ? _a : _b; })
#define ABS(x) __extension__ ({ __typeof(x) tmp = x; \
                    tmp < 0 ? -tmp : tmp; })

#define TOPPM     1.00E+06
#define TOPPB     1.00E+09
#define TOPPT     1.00E+12

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                   */
/* MAIN_ADJ - Main program - driver routine                          */
/*   Arguments :                                                     */
/*                                                                   */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   Driver for the Adjoint (ADJ) model
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

int KPP_Main_ADJ( const double finalPlume[], const double initBackg[],  \
                  const double temperature_K, const double pressure_Pa, \
                  const double airDens, const double timeArray[],       \
                  const unsigned int NT,                                \
                  const double RTOLS, const double ATOLS,               \
                  double VAR_OUTPUT[], const bool verbose )
{

    int i, j;
    int IERR = 1;

    /* Did we break? */
    bool BREAK = 0;

    /* Species considered for the optimization */
    int ind_OPT[NOPT] = {ind_NO, ind_NO2, ind_O3};

    /* ---- TIME VARIABLES ------------------ */


    const double TSTART = timeArray[0];
    const double TEND   = timeArray[NT-1];

    /* ---- INITIALIZE ARRAYS --------------- */
   
    double VAR_BACKG[NVAR];
    double VAR_RUN[NVAR];

    for ( i = 0; i < NVAR; i++ )
        VAR_BACKG[i] = initBackg[i];

    /* ---- INITIALIZE SENSITIVITIES -------- */

    /* NADJ = Number of functionals for which
     * sensitivities are computed.
     * Note: the setting below is for sensitivities
     * of all final concentrations. The setting may
     * have to be changed for other applications. */
    const int NADJ = NVAR;
    double Y_adj[NADJ][NVAR], ATOL_adj[NADJ][NVAR], RTOL_adj[NADJ][NVAR];

    for( i = 0; i < NADJ; i++ ) {
        for( j = 0; j < NVAR; j++ )
            Y_adj[i][j] = (double)0.0;
    }
    for( j = 0; j < NADJ; j++ )
        Y_adj[j][j] = (double)1.0;

    /* ---- TOLERANCES ---------------------- */

    double RTOL[NVAR];
    double ATOL[NVAR];
    for( i = 0; i < NVAR; i++ ) {
        RTOL[i] = RTOLS;
        ATOL[i] = ATOLS;
    }
    
    double STEPMIN = (double)0.0;

    /* Tolerances for calculating adjoints are 
     * used for controlling adjoint truncation 
     * error and for solving the linear adjoint
     * equations by iterations.
     * Note: Adjoints typically span many orders
     * of magnitude and a careful tuning of 
     * ATOL_adj may be necessary */

    for( i = 0; i < NADJ; i++ ) {
        for( j = 0; j < NVAR; j++ ) {
            RTOL_adj[i][j] = 1.0e-5;
            ATOL_adj[i][j] = 1.0e-10;
        }
    }

    /* ---- INTEGRATION SETTINGS/STATIS ----- */
    
    /* ICNTRL , RCNTRL  = Adjoint settings
     * ISTATUS, RSTATUS = Adjoint statistics */

    double RCNTRL[20], RSTATUS[20];
    int ICNTRL[20], ISTATUS[20];

    /* Default control options */
    for( i = 0; i < 20; i++ ) {
        ICNTRL[i] = 0;
        RCNTRL[i] = (double)0.0;
        ISTATUS[i] = 0;
        RSTATUS[i] = (double)0.0;
    }

    /*    ADJOINT INPUT PARAMETERS:

    Note: For input parameters equal to zero the default values of the
       corresponding variables are used.

    ICNTRL[0]   = 1: F = F(y)   Independent of T (AUTONOMOUS)
                = 0: F = F(t,y) Depends on T (NON-AUTONOMOUS)

    ICNTRL[1]   = 0: AbsTol, RelTol are NVAR-dimensional vectors
                = 1:  AbsTol, RelTol are scalars

    ICNTRL[2]  -> selection of a particular Rosenbrock method
        = 0 :  default method is Rodas3
        = 1 :  method is  Ros2
        = 2 :  method is  Ros3
        = 3 :  method is  Ros4
        = 4 :  method is  Rodas3
        = 5:   method is  Rodas4

    ICNTRL[3]  -> maximum number of integration steps
        For ICNTRL[4]=0) the default value of BUFSIZE is used

    ICNTRL[5]  -> selection of a particular Rosenbrock method for the
                continuous adjoint integration - for cts adjoint it
                can be different than the forward method ICNTRL(3)
         Note 1: to avoid interpolation errors (which can be huge!)
                   it is recommended to use only ICNTRL[6] = 2 or 4
         Note 2: the performance of the full continuous adjoint
                   strongly depends on the forward solution accuracy Abs/RelTol

    ICNTRL[6] -> Type of adjoint algorithm
         = 0 : default is discrete adjoint ( of method ICNTRL[3] )
         = 1 : no adjoint
         = 2 : discrete adjoint ( of method ICNTRL[3] )
         = 3 : fully adaptive continuous adjoint ( with method ICNTRL[6] )
         = 4 : simplified continuous adjoint ( with method ICNTRL[6] )

    ICNTRL[7]  -> checkpointing the LU factorization at each step:
        ICNTRL[7]=0 : do *not* save LU factorization (the default)
        ICNTRL[7]=1 : save LU factorization
        Note: if ICNTRL[7]=1 the LU factorization is *not* saved

    ~~~>  Real input parameters:

    RCNTRL[0]  -> Hmin, lower bound for the integration step size
          It is strongly recommended to keep Hmin = ZERO

    RCNTRL[1]  -> Hmax, upper bound for the integration step size

    RCNTRL[2]  -> Hstart, starting value for the integration step size

    RCNTRL[3]  -> FacMin, lower bound on step decrease factor (default=0.2)

    RCNTRL[4]  -> FacMax, upper bound on step increase factor (default=6)

    RCNTRL[5]  -> FacRej, step decrease factor after multiple rejections
            (default=0.1)

    RCNTRL[6]  -> FacSafe, by which the new step is slightly smaller
         than the predicted value  (default=0.9)

    RCNTRL[7]  -> ThetaMin. If Newton convergence rate smaller
                  than ThetaMin the Jacobian is not recomputed;
                  (default=0.001) */

    /*    ADJOINT OUTPUT PARAMETERS:

    Note: each call to RosenbrockADJ adds the corrent no. of fcn calls
      to previous value of ISTATUS(1), and similar for the other params.
      Set ISTATUS[0:9] = 0 before call to avoid this accumulation.

    ISTATUS[0] = No. of function calls
    ISTATUS[1] = No. of jacobian calls
    ISTATUS[2] = No. of steps
    ISTATUS[3] = No. of accepted steps
    ISTATUS[4] = No. of rejected steps (except at the beginning)
    ISTATUS[5] = No. of LU decompositions
    ISTATUS[6] = No. of forward/backward substitutions
    ISTATUS[7] = No. of singular matrix decompositions

    ~~~>  Real output parameters:

    RSTATUS[0]  -> Texit, the time corresponding to the
                   computed Y upon return
    RSTATUS[1]  -> Hexit, last accepted step before exit
    For multiple restarts, use Hexit as Hstart in the following run */
    
    /* Running the discrete adjoint */
    ICNTRL[6] = 2;
    
    /* ---- INITIALIZE ARRAYS FOR CURR. RUN - */

    for ( i = 0; i < NVAR; i++ )
        VAR_RUN[i] = VAR_BACKG[i];
                
    for ( i = 0; i < NREACT; i++ )
        RCONST[i] = 0.0E+00;

    for ( i = 0; i < NSPEC; i++ ) {
        HET[i][0] = 0.0E+00;
        HET[i][1] = 0.0E+00;
        HET[i][2] = 0.0E+00;
    }

    Update_RCONST( temperature_K, pressure_Pa, airDens, VAR_RUN[ind_H2O] );

    /* ---- COMPUTE SENSITIVITIES ----------- */

    TIME = TSTART;
    IERR = INTEGRATE_ADJ( NADJ, VAR_RUN, Y_adj, TSTART, TEND, ATOL_adj, RTOL_adj, ATOL, RTOL, ICNTRL,
                          RCNTRL, ISTATUS, RSTATUS, STEPMIN );

    /* If integration failed, stop here */
    if ( IERR < 0 ) {
        printf(" Initial adjoint integration failed\n");
        return IERR;
    }

    /* Compute difference between forward ambient
     * run and final plume concentrations */
    double DELTAO3f, DELTANOf, DELTANO2f;
    
    DELTAO3f  = VAR_RUN[ind_O3]  - finalPlume[ind_O3];
    DELTANOf  = VAR_RUN[ind_NO]  - finalPlume[ind_NO];
    DELTANO2f = VAR_RUN[ind_NO2] - finalPlume[ind_NO2];

    const double DELTAO3f_0  = DELTAO3f;
    const double DELTANOf_0  = DELTANOf;
    const double DELTANO2f_0 = DELTANO2f;

    if ( verbose ) {
        printf(" DELTAO3 : %f [ppt]\n", DELTAO3f  / airDens * TOPPT);
        printf(" DELTANO : %f [ppt]\n", DELTANOf  / airDens * TOPPT);
        printf(" DELTANO2: %f [ppt]\n", DELTANO2f  / airDens * TOPPT);
    }

    /* Optimization method:
     * (1) Conjugate gradient */

    static const char switchOpt[] = "Conjugate Gradient";

    if ( verbose )
        printf(" Running KPP adjoint to compute sensitivities and effective emission indices using %s.\n", switchOpt);

    if (( strcmp( switchOpt, "CG" ) == 0 ) || \
        ( strcmp( switchOpt, "Conjugate Gradient" ) == 0 ) ) {

        /* Optimization problem:
         *
         * Let K : R^NVAR -> R^NVAR
         *         x0     -> x(t_f)
         * where x(t_f) is the solution array of 
         *   | dy/dt  = F(t,y)
         *   | y(t_0) = x0
         * evaluated at t = t_f.
         *
         * The goal is to find x0* such that:
         *
         *      x0* = arg min 1/2*( K(x0) - x* )^T * W * ( K(x0) - x* )
         *          = arg min f(x0)
         *
         * where x* is given by the forward plume model. 
         * W is a diagonal matrix containing weights
         *
         * We use the following algorithm:
         *
         * (1) Initialize x0
         * (2) Compute the sensitivity matrix S = (dy_i(t_f)/dy_j(t_0)) 
         *     using KPP_Adjoint. Set the reduced sensitivity matrix
         *     J to S(ind_SPC(1:NOPT),ind_SPC(1:NOPT)).
         * (3) Compute the new direction d = - J^T * W * ( K(x0) - x* )
         * (4) Solve the optimization problem:
         *     alpha* = arg min f( x0 + alpha * d )
         * (5) Update x0 += alpha* * d
         * (6) If f(x0) >= THRESH, return to (2), otherwise break;
         *
         */

        const int N_MAX = 20;
        int N = 0;

        double *wDiff, *dir, *dirold;

        wDiff  = (double *) malloc( sizeof(double) * NOPT );
        if ( wDiff == NULL ) {
            printf(" malloc of size %d failed!\n", NOPT);
            return -1;
        }
        dir    = (double *) malloc( sizeof(double) * NOPT );
        if ( dir == NULL ) {
            printf(" malloc of size %d failed!\n", NOPT);
            return -1;
        }
        dirold = (double *) malloc( sizeof(double) * NOPT );
        if ( dirold == NULL ) {
            printf(" malloc of size %d failed!\n", NOPT);
            return -1;
        }

        double **A_mat; 

        A_mat = (double **) malloc( sizeof(double *) * NOPT );
        for ( i = 0; i < NOPT; i++ ) {
            A_mat[i] = (double *) malloc( sizeof(double) * NOPT );
            if ( A_mat[i] == NULL ) {
                printf(" malloc of size %d failed!\n", NOPT);
                return -1;
            }
        }
        
        for ( i = 0; i < NOPT; i++ ) {
            wDiff[i]  = 0.0E+00; 
            dir[i]    = 0.0E+00;
            dirold[i] = 0.0E+00;
            for ( j = 0; j < NOPT; j++ )
                A_mat[i][j] = 0.0E+00;
        }

        /* Compute weights */
        double DIAG[NOPT];

        /* Assigning arrays */
        for ( i = 0; i < NOPT; i++ ) {
            for ( j = 0; j < NOPT; j++ ) {
                /* Assign the reduced sensitivity matrix: J^T */
                A_mat[i][j] = (double)Y_adj[ind_OPT[j]][ind_OPT[i]];
            }
//            if ( i <= 1 ) // NOx
//                DIAG[i] = 1.0 / (VAR_BACKG[ind_OPT[0]] + VAR_BACKG[ind_OPT[1]] - finalPlume[ind_OPT[0]] - finalPlume[ind_OPT[1]]);
//            else if ( i > 1 ) // O3
//                DIAG[i] = 1.0 / (VAR_BACKG[ind_OPT[i]] - finalPlume[ind_OPT[i]]); 
            if ( i <= 1 ) // NOx
                DIAG[i] = 1.0 / ( VAR_BACKG[ind_OPT[0]] + VAR_BACKG[ind_OPT[1]] );
            else if ( i > 1 ) // O3
                DIAG[i] = 1.0 / ( VAR_BACKG[ind_OPT[i]] ); 

            wDiff[i] = DIAG[i] * ( VAR_RUN[ind_OPT[i]] - finalPlume[ind_OPT[i]] );
        }

        /* Print reduced sensitivity matrix and weighted difference
         * to console? */
        if ( verbose ) {
            printf(" Reduced sensitivity matrix\n");
            printf(" J^T = [");
            for ( i = 0; i < NOPT; i++ ) {
                if ( i > 0 )
                    printf("       [");
                for ( j = 0; j < NOPT; j++ )
                    printf("%+e ",A_mat[i][j]);
                printf("]\n");
            }
            printf(" Weighted difference\n");
            printf(" wDiff = [");
            for ( i = 0; i < NOPT; i++ )
                printf("%+e ",wDiff[i]);
            printf("]\n");
        }

        /* Compute the new direction */
        for ( i = 0; i < NOPT; i++ ) {
            dir[i] = 0.0E+00;
            for ( j = 0; j < NOPT; j++ )
                dir[i] -= A_mat[i][j] * wDiff[j];
        }

        double VAR_INIT[NVAR];
        double VAR_DIR[NVAR];
        double METRIC_min;
        double METRIC;
        double METRIC_OLD;
        double METRIC_prev;
        double METRIC_ABS_MIN;

        double alpha[NOPT] = { 0.0 };
        double alpha_maxstep = -1.0;
        double alpha_step;

        double beta = 1.0;
        
        METRIC_min = 1.00E+40;
        METRIC = METRIC_min;
        METRIC_OLD = METRIC_min;
        METRIC_prev = METRIC_min;
        METRIC_ABS_MIN = METRIC_min;
        
        for ( i = 0; i < NVAR; i++ ) {
            VAR_RUN[i]     = VAR_BACKG[i];
            VAR_INIT[i]    = VAR_BACKG[i];
            VAR_DIR[i]     = VAR_BACKG[i];
            VAR_OUTPUT[i]  = VAR_BACKG[i];
        }
        VAR_RUN[ind_O3]    -= DELTAO3f/Y_adj[ind_O3][ind_O3];
        VAR_INIT[ind_O3]   -= DELTAO3f/Y_adj[ind_O3][ind_O3];
        VAR_DIR[ind_O3]    -= DELTAO3f/Y_adj[ind_O3][ind_O3];
        VAR_OUTPUT[ind_O3] -= DELTAO3f/Y_adj[ind_O3][ind_O3];
              
        for ( i = 0; i < NOPT; i++ ) {
            alpha[i] = ABS( VAR_DIR[ind_OPT[i]] / dir[i] );
            if ( alpha_maxstep < alpha[i] )
                alpha_maxstep = alpha[i];
        }

        int imin = 0;
        const int nSTEP = 25;
        double STEPARRAY[nSTEP];
        STEPARRAY[0] = 5.00E+04;
        for ( i = 1; i < nSTEP; i++ )
            STEPARRAY[i] = 4.0*STEPARRAY[i-1];

        /* ---- OPTIMIZATION STARTS HERE -------------------- */

        while ( N < N_MAX ) {

            if ( verbose )
                printf("\nN: %d, METRIC: %6.5e, Ratio, %e\n\n", N, METRIC, METRIC_OLD/METRIC);

            METRIC_OLD = METRIC;


            METRIC_min = 1.00E+40;
            METRIC_prev = 1.00E+40;

            /* ---- LINE SEARCH BEGINS HERE --------- */
            
            ICNTRL[6] = 1; // No adjoint evaluation

            j = 0;
            while ( j < nSTEP ) {
            
                alpha_step = alpha_maxstep/STEPARRAY[j];

                for ( i = 0; i < NVAR; i++ )
                    VAR_RUN[i] = VAR_DIR[i];
                for ( i = 0; i < NOPT; i++ ) {
                    VAR_RUN[ind_OPT[i]] = VAR_DIR[ind_OPT[i]] + alpha_step * dir[i]; // / DIAG[i];
                    VAR_INIT[ind_OPT[i]] = VAR_RUN[ind_OPT[i]];
                }

                /* Ensure that NOy is conserved */
                VAR_RUN[ind_HNO3] -= alpha_step * ( dir[0] + dir[1] ); //( dir[0] / DIAG[0] + dir[1] / DIAG[1] );
                VAR_INIT[ind_HNO3] = VAR_RUN[ind_HNO3];
    
                for ( i = 0; i < NREACT; i++ )
                    RCONST[i] = 0.0E+00;

                for ( i = 0; i < NSPEC; i++ ) {
                    HET[i][0] = 0.0E+00;
                    HET[i][1] = 0.0E+00;
                    HET[i][2] = 0.0E+00;
                }

                Update_RCONST( temperature_K, pressure_Pa, airDens, VAR_RUN[ind_H2O] );
    
                TIME = TSTART;
                IERR = INTEGRATE_ADJ( NADJ, VAR_RUN, Y_adj, TSTART, TEND, ATOL_adj, RTOL_adj, ATOL, RTOL, ICNTRL, RCNTRL, ISTATUS, RSTATUS, STEPMIN );
                
                /* Compute metric */
                METRIC = 0.0;
                for ( i = 0; i < NOPT - 1; i++ ) {
                    if ( i < 1 )
                        METRIC += pow( ( DIAG[0] * ( VAR_RUN[ind_OPT[0]] + VAR_RUN[ind_OPT[1]] - finalPlume[ind_OPT[0]] - finalPlume[ind_OPT[1]] ) ), 2);
                    else
                        METRIC += pow( ( DIAG[i+1] * ( VAR_RUN[ind_OPT[i+1]] - finalPlume[ind_OPT[i+1]] ) ), 2);
                }

                if ( verbose )
                    printf("Step: %d, METRIC: %e\n", j, METRIC);

                if ( METRIC < METRIC_min ) {
                    METRIC_min = METRIC;
                    imin = j;
                }

                if ( METRIC < METRIC_ABS_MIN ) {
                    METRIC_ABS_MIN = METRIC;
                    for ( i = 0; i < NOPT; i++ ) {
                        VAR_OUTPUT[ind_OPT[i]] = VAR_INIT[ind_OPT[i]];
                        VAR_OUTPUT[ind_HNO3]   = VAR_INIT[ind_HNO3];
                    }
                }

                j++;

                if ( ( METRIC > METRIC_prev ) && ( METRIC < 2.0 * METRIC_ABS_MIN ) ){
                    IERR = 0;
                    BREAK = 1;
                    break;
                }

                METRIC_prev = METRIC;

            }
            
            /* ---- LINE SEARCH ENDS HERE ----------- */
            
            if ( ( N > 4 ) && ( ABS( METRIC_OLD / METRIC - 1 ) < 1.00E-04 ) ) {
                IERR = 0;
                BREAK = 1;
                break;
            }

            alpha_step = alpha_maxstep/STEPARRAY[imin];

            for ( i = 0; i < NVAR; i++ )
                VAR_RUN[i] = VAR_DIR[i];
            for ( i = 0; i < NOPT; i++ ) {
                VAR_RUN[ind_OPT[i]]  = VAR_DIR[ind_OPT[i]] + alpha_step * dir[i];
                VAR_INIT[ind_OPT[i]] = VAR_RUN[ind_OPT[i]];
            }
            VAR_RUN[ind_HNO3] -= alpha_step * ( dir[0] + dir[1] );
            VAR_INIT[ind_HNO3] = VAR_RUN[ind_HNO3];

            /* Adjoint computation */
            ICNTRL[6] = 2;

            /* ---- INITIALIZE SENSITIVITIES ---------- */

            for( i = 0; i < NADJ; i++ ) {
                for( j = 0; j < NVAR; j++ )
                    Y_adj[i][j] = (double)0.0;
            }
            for( j = 0; j < NADJ; j++ )
                Y_adj[j][j] = (double)1.0;

            /* ---- INITIALIZE RATES ------------------ */

            for ( i = 0; i < NREACT; i++ )
                RCONST[i] = 0.0E+00;

            for ( i = 0; i < NSPEC; i++ ) {
                HET[i][0] = 0.0E+00;
                HET[i][1] = 0.0E+00;
                HET[i][2] = 0.0E+00;
            }

            Update_RCONST( temperature_K, pressure_Pa, airDens, VAR_RUN[ind_H2O] );

            /* ---- COMPUTE SENSITIVITIES -------------- */

            TIME = TSTART;
            IERR = INTEGRATE_ADJ( NADJ, VAR_RUN, Y_adj, TSTART, TEND, ATOL_adj, RTOL_adj, ATOL, RTOL, ICNTRL, RCNTRL, ISTATUS, RSTATUS, STEPMIN );

            if ( IERR < 0 ) {
                printf(" Adjoint integration failed\n Metric: %e\n", METRIC);
                printf(" [NO](t = t0)    = %f [ppt]\n", VAR_INIT[ind_NO]/airDens*TOPPT);
                printf(" [NO2](t = t0)   = %f [ppt]\n", VAR_INIT[ind_NO2]/airDens*TOPPT);
                printf(" [O3](t = t0)    = %f [ppb]\n", VAR_INIT[ind_O3]/airDens*TOPPB);
                printf(" [HNO3](t = t0)  = %f [ppb]\n", VAR_INIT[ind_HNO3]/airDens*TOPPB);
                return IERR;
            }
            
            // if O3 delta smaller than 0.5 ppt and NOx delta smaller than 0.05 ppt
            if ( ABS((VAR_RUN[ind_O3] - finalPlume[ind_O3])/airDens*TOPPT) < 0.1 && ABS((VAR_RUN[ind_NO] + VAR_RUN[ind_NO2] - finalPlume[ind_NO] - finalPlume[ind_NO2])/airDens*TOPPT) < 0.05 ) {
                IERR = 0;
                BREAK = 1;
                break;
            }

            METRIC = 0.0;
            for ( i = 0; i < NOPT - 1; i++ ) {
                if ( i < 1 ) //NOx
                    METRIC += pow( ( DIAG[0] * ( VAR_RUN[ind_OPT[0]] + VAR_RUN[ind_OPT[1]] - finalPlume[ind_OPT[0]] - finalPlume[ind_OPT[1]] ) ), 2);
                else
                    METRIC += pow( ( DIAG[i+1] * ( VAR_RUN[ind_OPT[i+1]] - finalPlume[ind_OPT[i+1]] ) ), 2);
            }

            if ( verbose ) {
                printf("\n METRIC: %6.5e\n\n",METRIC);
            }

            for ( i = 0; i < NVAR; i++ )
                VAR_DIR[i] = VAR_INIT[i];
        
            DELTAO3f = VAR_RUN[ind_O3] - finalPlume[ind_O3];
//            VAR_DIR[ind_O3] -= DELTAO3f/Y_adj[ind_O3][ind_O3];
            
            for ( i = 0; i < NVAR; i++ )
                VAR_RUN[i] = VAR_DIR[i];
            VAR_RUN[ind_HNO3] -= alpha_step * ( dir[0] + dir[1] );

            /* No adjoint computation */
            ICNTRL[6] = 1;

            /* ---- INITIALIZE RATES ------------------ */

            for ( i = 0; i < NREACT; i++ )
                RCONST[i] = 0.0E+00;

            for ( i = 0; i < NSPEC; i++ ) {
                HET[i][0] = 0.0E+00;
                HET[i][1] = 0.0E+00;
                HET[i][2] = 0.0E+00;
            }

            Update_RCONST( temperature_K, pressure_Pa, airDens, VAR_RUN[ind_H2O] );

            /* ---- COMPUTE SENSITIVITIES -------------- */

            TIME = TSTART;
            IERR = INTEGRATE_ADJ( NADJ, VAR_RUN, Y_adj, TSTART, TEND, ATOL_adj, RTOL_adj, ATOL, RTOL, ICNTRL, RCNTRL, ISTATUS, RSTATUS, STEPMIN );

            if ( IERR < 0 ) {
                printf(" Forward integration failed\n");
                printf(" [NO](t = t0)    = %f [ppt]\n", VAR_INIT[ind_NO]/airDens*TOPPT);
                printf(" [NO2](t = t0)   = %f [ppt]\n", VAR_INIT[ind_NO2]/airDens*TOPPT);
                printf(" [O3](t = t0)    = %f [ppb]\n", VAR_INIT[ind_O3]/airDens*TOPPB);
                printf(" [HNO3](t = t0)  = %f [ppb]\n", VAR_INIT[ind_HNO3]/airDens*TOPPB);
                return IERR;
            }

            /* Assigning new data */
            for ( i = 0; i < NOPT; i++ ) {
                for ( j = 0; j < NOPT; j++ )
                    A_mat[i][j] = (double) Y_adj[ind_OPT[j]][ind_OPT[i]];
                wDiff[i] = DIAG[i] * ( VAR_RUN[ind_OPT[i]] - finalPlume[ind_OPT[i]] );
            }
        
            /* Print reduced sensitivity matrix and weighted difference
             * to console? */
            if ( verbose ) {
                printf(" Printing reduced sensitivity matrix\n");
                printf(" J^T = [");
                for ( i = 0; i < NOPT; i++ ) {
                    if ( i > 0 )
                        printf("       [");
                    for ( j = 0; j < NOPT; j++ )
                        printf("%+e ",A_mat[i][j]);
                    printf("]\n");
                }
                printf(" Weighted difference\n");
                printf(" wDiff = [");
                for ( i = 0; i < NOPT; i++ )
                    printf("%+e ",wDiff[i]);
                printf("]\n");
            }

            /* ---- COMPUTE NEW DIRECTION ----------------- */

            /* Compute the new direction */
            for ( i = 0; i < NOPT; i++ ) {
                dirold[i] = dir[i];
                dir[i] = 0.0E+00;
                for ( j = 0; j < NOPT; j++ )
                    dir[i] -= A_mat[i][j] * wDiff[j];
            }
            
            /* Update the CG parameter */
            beta = 0.0E+00;
            /* Only select ONE of the following update formula!! */
            for ( i = 0; i < NOPT; i++ ) {
                /* Fletcher-Reeves */
//                beta += dir[i] * dir[i] / ( dirold[i] * dirold[i] );
                /* Polak-Ribière */
//                beta += dir[i] * ( dir[i] - dirold[i] ) / ( dirold[i] * dirold[i] );
                /* Hestenes-Stiefel */
                /* ... */
                beta = MAX( 0, beta );
            }

            if ( verbose ) {
                printf(" beta = %e\n", beta);
                printf(" dir = [ %e, %e, %e ]\n",dir[0],dir[1],dir[2]);
            }


            for ( i = 0; i < NOPT; i++ )
                dir[i] += beta * dirold[i];

            if ( verbose )
                printf(" dir = [ %e, %e, %e ]\n",dir[0],dir[1],dir[2]);

            if ( verbose ) {
                printf(" O3: %2.7f, %2.7f, %2.7f\n",VAR_RUN[ind_O3]/airDens*1E9, finalPlume[ind_O3]/airDens*1E9,(VAR_RUN[ind_O3] - finalPlume[ind_O3])/airDens*TOPPT);
                printf(" NOx: %2.7f, %2.7f, %2.7f\n",(VAR_RUN[ind_NO]+VAR_RUN[ind_NO2])/airDens*TOPPT, (finalPlume[ind_NO] + finalPlume[ind_NO2])/airDens*TOPPT, (VAR_RUN[ind_NO]+VAR_RUN[ind_NO2] - finalPlume[ind_NO] - finalPlume[ind_NO2])/airDens*TOPPT);
            }

            alpha_maxstep = -1.0;
            /* Computing step */
            for ( i = 0; i < NOPT; i++ ) {
                alpha[i] = ABS( VAR_DIR[ind_OPT[i]] / dir[i] );
                if ( alpha_maxstep < alpha[i] )
                    alpha_maxstep = alpha[i];
            }

            IERR = 0;
            
            N++;

        }
       
        if ( verbose )
            printf("\n Trying to tweak ozone concentration\n");

        /* Try tweaking the ozone value */
        for ( i = 0; i < NVAR; i++ )
            VAR_RUN[i] = VAR_OUTPUT[i];
        for ( i = 0; i < NOPT; i++ )
            VAR_INIT[ind_OPT[i]] = VAR_RUN[ind_OPT[i]];
        VAR_INIT[ind_HNO3] = VAR_RUN[ind_HNO3];

        /* No adjoint computation */
        ICNTRL[6] = 1;

        /* ---- INITIALIZE RATES ------------------ */

        for ( i = 0; i < NREACT; i++ )
            RCONST[i] = 0.0E+00;

        for ( i = 0; i < NSPEC; i++ ) {
            HET[i][0] = 0.0E+00;
            HET[i][1] = 0.0E+00;
            HET[i][2] = 0.0E+00;
        }

        Update_RCONST( temperature_K, pressure_Pa, airDens, VAR_RUN[ind_H2O] );

        /* ---- COMPUTE SENSITIVITIES -------------- */

        TIME = TSTART;
        IERR = INTEGRATE_ADJ( NADJ, VAR_RUN, Y_adj, TSTART, TEND, ATOL_adj, RTOL_adj, ATOL, RTOL, ICNTRL, RCNTRL, ISTATUS, RSTATUS, STEPMIN );

        if ( IERR < 0 ) {
            printf(" Forward integration failed\n");
            printf(" [NO](t = t0)    = %f [ppt]\n", VAR_INIT[ind_NO]/airDens*TOPPT);
            printf(" [NO2](t = t0)   = %f [ppt]\n", VAR_INIT[ind_NO2]/airDens*TOPPT);
            printf(" [O3](t = t0)    = %f [ppb]\n", VAR_INIT[ind_O3]/airDens*TOPPB);
            printf(" [HNO3](t = t0)  = %f [ppb]\n", VAR_INIT[ind_HNO3]/airDens*TOPPB);
            return IERR;
        }
                
        /* Compute metric */
        METRIC_ABS_MIN = 0.0;
        for ( i = 0; i < NOPT - 1; i++ ) {
            if ( i < 1 )
                METRIC_ABS_MIN += pow( ( DIAG[0] * ( VAR_RUN[ind_OPT[0]] + VAR_RUN[ind_OPT[1]] - finalPlume[ind_OPT[0]] - finalPlume[ind_OPT[1]] ) ), 2);
            else
                METRIC_ABS_MIN += pow( ( DIAG[i+1] * ( VAR_RUN[ind_OPT[i+1]] - finalPlume[ind_OPT[i+1]] ) ), 2);
        }

        if ( verbose )
            printf("METRIC: %e\n", METRIC_ABS_MIN);

            
        /* Try tweaking the ozone value */
        DELTAO3f = VAR_RUN[ind_O3] - finalPlume[ind_O3];
        for ( i = 0; i < NVAR; i++ )
            VAR_RUN[i] = VAR_OUTPUT[i];
        VAR_RUN[ind_O3]   -= DELTAO3f/Y_adj[ind_O3][ind_O3];
        for ( i = 0; i < NOPT; i++ )
            VAR_INIT[ind_OPT[i]] = VAR_RUN[ind_OPT[i]];
        VAR_INIT[ind_HNO3] = VAR_RUN[ind_HNO3];

        /* No adjoint computation */
        ICNTRL[6] = 1;

        /* ---- INITIALIZE RATES ------------------ */

        for ( i = 0; i < NREACT; i++ )
            RCONST[i] = 0.0E+00;

        for ( i = 0; i < NSPEC; i++ ) {
            HET[i][0] = 0.0E+00;
            HET[i][1] = 0.0E+00;
            HET[i][2] = 0.0E+00;
        }

        Update_RCONST( temperature_K, pressure_Pa, airDens, VAR_RUN[ind_H2O] );

        /* ---- COMPUTE SENSITIVITIES -------------- */

        TIME = TSTART;
        IERR = INTEGRATE_ADJ( NADJ, VAR_RUN, Y_adj, TSTART, TEND, ATOL_adj, RTOL_adj, ATOL, RTOL, ICNTRL, RCNTRL, ISTATUS, RSTATUS, STEPMIN );

        if ( IERR < 0 ) {
            printf(" Forward integration failed\n");
            printf(" [NO](t = t0)    = %f [ppt]\n", VAR_INIT[ind_NO]/airDens*TOPPT);
            printf(" [NO2](t = t0)   = %f [ppt]\n", VAR_INIT[ind_NO2]/airDens*TOPPT);
            printf(" [O3](t = t0)    = %f [ppb]\n", VAR_INIT[ind_O3]/airDens*TOPPB);
            printf(" [HNO3](t = t0)  = %f [ppb]\n", VAR_INIT[ind_HNO3]/airDens*TOPPB);
            return IERR;
        }
                
        /* Compute metric */
        METRIC = 0.0;
        for ( i = 0; i < NOPT - 1; i++ ) {
            if ( i < 1 )
                METRIC += pow( ( DIAG[0] * ( VAR_RUN[ind_OPT[0]] + VAR_RUN[ind_OPT[1]] - finalPlume[ind_OPT[0]] - finalPlume[ind_OPT[1]] ) ), 2);
            else
                METRIC += pow( ( DIAG[i+1] * ( VAR_RUN[ind_OPT[i+1]] - finalPlume[ind_OPT[i+1]] ) ), 2);
        }

        if ( verbose )
            printf("METRIC: %e\n", METRIC);

        if ( METRIC < METRIC_ABS_MIN ) {
            METRIC_ABS_MIN = METRIC;
            if ( verbose )
                printf("Tweaking worked!\n");
            for ( i = 0; i < NOPT; i++ )
                VAR_OUTPUT[ind_OPT[i]] = VAR_INIT[ind_OPT[i]];
            VAR_OUTPUT[ind_HNO3]   = VAR_INIT[ind_HNO3];
        }

        /* ---- OPTIMIZATION ENDS HERE ---------------------- */


#pragma omp critical
        {
            #ifdef OMP
                printf("\n ## ON THREAD: %d.", omp_get_thread_num());
            #endif /* OMP */
            printf(" Integration successful (METRIC: %e), break = %d\n", METRIC_ABS_MIN, BREAK);
            printf(" ## O3 Delta : %+f [ppt], %f %% \n",(VAR_RUN[ind_O3] - finalPlume[ind_O3])/airDens*TOPPT, 100.0 * ABS((VAR_RUN[ind_O3] - finalPlume[ind_O3]) / DELTAO3f_0));
            printf(" ## NOx Delta: %+f [ppt], %f %% \n",(VAR_RUN[ind_NO] + VAR_RUN[ind_NO2] - finalPlume[ind_NO] - finalPlume[ind_NO2])/airDens*TOPPT, 100.0 * ABS((VAR_RUN[ind_NO] + VAR_RUN[ind_NO2] - finalPlume[ind_NO] - finalPlume[ind_NO2])/ ( DELTANOf_0 + DELTANO2f_0 )));
        }

    }

    return IERR;

}

/* End of MAIN_ADJ function                                          */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

