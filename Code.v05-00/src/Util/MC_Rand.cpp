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
#include <ctime>
#include <iostream>
#include "APCEMM.h"
#include "Util/MC_Rand.hpp"
#include "Core/Input_Mod.hpp"

unsigned int setSeed(OptInput& input) {

    // Sets the seed of the pseudo-random generator, stores it in OptInput and returns
    // it, so the caller can record the seed the run actually used.
    #ifdef DEBUG
        // With DEBUG compile flag set a constant seed for reproducibility
        input.SIMULATION_SEED_VALUE = 0;
        std::cout << "Compiled in DEBUG mode: random seed is set to " << input.SIMULATION_SEED_VALUE << " for all simulations" << std::endl;
    #else
        if(!input.SIMULATION_FORCE_SEED){
            // If the seed is not being forced to a value use the current unix timestamp
            // instead. The conversion from time_t to unsigned int wraps (time_t is larger)
            // but that's fine as any unsigned int is a usable seed.
            input.SIMULATION_SEED_VALUE = static_cast<unsigned int>(std::time(nullptr));
        }
        std::cout << "Random seed is set to " << input.SIMULATION_SEED_VALUE << " for all simulations" << std::endl;
    #endif

    srand(input.SIMULATION_SEED_VALUE);
    return input.SIMULATION_SEED_VALUE;
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
