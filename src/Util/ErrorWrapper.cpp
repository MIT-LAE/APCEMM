/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* ErrorWrapper Program File                                        */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 10/17/2018                                */
/* File                 : ErrorWrapper.cpp                          */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "Util/Error.hpp"
#include "Util/ErrorWrapper.h"

using namespace Error;

extern "C" {

    bool SafeDiv( float num, float denom )
    {
     
        return Error::SafeDiv( num, denom );

    } /* End of SafeDiv */

    bool SafeDiv( double num, double denom )
    {
     
        return Error::SafeDiv( num, denom );

    } /* End of SafeDiv */


}

/* End of ErrorWrapper.cpp */
