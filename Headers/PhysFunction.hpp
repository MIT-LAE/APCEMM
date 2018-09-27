/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*     Aircraft Plume Chemistry, Emission and Microphysics Model    */
/*                             (APCEMM)                             */
/*                                                                  */
/* PhysFunction Header File                                         */
/*                                                                  */
/* Author               : Thibaud M. Fritz                          */
/* Time                 : 7/26/2018                                 */
/* File                 : PhysFunction.hpp                          */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef PHYSFUNCTION_H_INCLUDED
#define PHYSFUNCTION_H_INCLUDED

#include <cmath>

#include "PhysConstant.hpp"

class physFunc
{

    public:

        physFunc();
        ~physFunc();

        double pSat_H2Ol( double T );
        double pSat_H2Os( double T );
        double pSat_H2SO4( double T );
        double pSat_HNO3( double T , double PPH2O );

};

#endif /* PHYSFUNCTION_H_INCLUDED */
