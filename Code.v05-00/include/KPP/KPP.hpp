/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/*                                                                  */
/* KPP Header File                                                  */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 10/4/2018                                 */
/* File                 : KPP.hpp                                   */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef KPP_H_INCLUDED
#define KPP_H_INCLUDED

#include "KPP/KPP_Parameters.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

int INTEGRATE( double VAR[] , double FIX[], double TIN   , double TOUT, \
               double ATOL[], double RTOL[], double STEPMIN );
int KPP_Main_ADJ( const double finalPlume[], const double initBackg[],  \
                  const double temperature_K, const double pressure_Pa, \
                  const double airDens, const double timeArray[],       \
                  const unsigned int NT,                                \
                  const double RTOLS, const double ATOLS,               \
                  double VAR_OUTPUT[], double *METRIC,                  \
                  const bool VERBOSE = 0, const bool RETRY = 0 );
int INTEGRATE_ADJ( int NADJ, double Y[], double Lambda[][NVAR],        \
		           double TIN, double TOUT, double ATOL_adj[][NVAR],   \
        	       double RTOL_adj[][NVAR], double ATOL[],             \
                   double RTOL[], int ICNTRL_U[],                      \
		           double RCNTRL_U[], int ISTATUS_U[],                 \
                   double RSTATUS_U[], double STEPMIN );
void Update_RCONST( const double TEMP, const double PRESS,  \
                    const double AIRDENS, const double H2O );
void GC_SETHET( const double TEMP, const double PATM, const double AIRDENS, \
                const double RELHUM, const unsigned int STATE_PSC,          \
                const double SPC[], const double AREA[NAERO],               \
                const double RADI[NAERO], const double IWC,                 \
                const double KHETI_SLA[11], double tropopausePressure);
void Update_JRates ( double JRates[], const double CSZA );
void ComputeFamilies( const double V[], const double F[], const double RCT[], \
                      double familyRates[] );

static int KPP_FAIL    = -1;
// static int KPPADJ_FAIL = -5;
    

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* KPP_H_INCLUDED */
