/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* MC_Rand Program File                                             */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 1/25/2019                                 */
/* File                 : MC_Rand.cpp                               */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
#include <iostream>
#include "APCEMM.h"
#include "Util/MC_Rand.hpp"

void setSeed() {

    // Sets seed for pseudo-random generator.
    #ifdef DEBUG
        // With DEBUG compile flag set a constant seed for reproducibility
        std::cout << "Compiled in DEBUG mode: random seed is set to 0 for all simulations" << std::endl;
        srand(0);
    #else
        // Otherwise use the current unix timestamp as our random seed. 
        srand(time(NULL));
    #endif

} /* End of setSeed */

template <typename T>
T fRand(const T fMin, const T fMax) {

    /* Returns a random number between fMin and fMax */

    double f = (double) rand()/RAND_MAX;
    return (T) fMin + f * (fMax - fMin);

} /* End of fRand */

template double fRand(const double fMin, const double fMax);
template float fRand(const float fMin, const float fMax);
template int fRand(const int fMin, const int fMax);
template unsigned int fRand(const unsigned int fMin, const unsigned int fMax);

/* End of MC_Rand.cpp */
