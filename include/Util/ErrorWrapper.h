/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* ErrorWrapper Header File                                         */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 10/17/2018                                */
/* File                 : ErrorWrapper.h                            */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef ERRORWRAPPER_H_INCLUDED
#define ERRORWRAPPER_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

    bool SafeDiv_f( float num, float denom );
    bool SafeDiv_d( double num, double denom );

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* ERRORWRAPPER_H_INCLUDED */
