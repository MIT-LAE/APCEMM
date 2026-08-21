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
#include "Core/Input_Mod.hpp"

// Set the seed of the pseudo-random generator and return it.
// Also writes the resolved seed back into input.SIMULATION_SEED_VALUE.
unsigned int setSeed(OptInput& input);

/* Generates a random number of type T between fMin and fMax */
template <typename T>
T fRand(const T fMin, const T fMax);

#endif /* MC_RAND_H_INCLUDED */
