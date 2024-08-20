/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* MC_Rand Header File                                              */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 1/25/2019                                 */
/* File                 : MC_Rand.hpp                               */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef MC_RAND_H_INCLUDED
#define MC_RAND_H_INCLUDED

#include <cstdlib>

/* Set seed for pseudo-random generator */
void setSeed();

/* Generates a random number of type T between fMin and fMax */
template <typename T>
T fRand(const T fMin, const T fMax);

#endif /* MC_RAND_H_INCLUDED */
