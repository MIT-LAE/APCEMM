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

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
int KPP_Main( double varArray[], double fixArray[], double currentT, double dt, \
              double airDens, double temperature, double pressure, \
              double sinLAT, double cosLAT, double sinDEC, double cosDEC, \
              double rtols, double atols );
#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* KPP_H_INCLUDED */
