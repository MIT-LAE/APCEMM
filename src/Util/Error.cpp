/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* Error Program File                                               */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 10/17/2018                                */
/* File                 : Error.cpp                                 */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "Util/Error.hpp"

namespace Error
{

    template<typename T> bool SafeDiv( T num, T denom ) 
    {

        typedef std::numeric_limits<T> limits;

        if ( denom == 0 ) {
            return 0;
        } else {
            if ( ( (num / denom) > limits::max() ) || ( (num/denom) < limits::min() ) )
                return 0;
            else
                return 1;
        }

    } /* End of SafeDiv */

}

/* End of Error.cpp */
