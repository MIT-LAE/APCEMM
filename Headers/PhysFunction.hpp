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

        /* Saturation pressures [Pa] */
        double pSat_H2Ol( double T );
        double pSat_H2Os( double T );
        double pSat_H2SO4( double T );
        double pSat_HNO3( double T , double PPH2O );

        /* Density of air [kg/m^3] */
        double rhoAir( double T, double P );

        /* Dynamic viscosity [kg/(m s)] */
        double dynVisc( double T );

        /* Kinematic viscosity [m^2/s] */
        double kinVisc( double T, double P );

        /* Thermal speed of an air molecule [m/s] */
        double thermalSpeed( double T );

        /* Mean free path [m] */
        double lambda( double T, double P );

        /* Knudsen number [-] */
        double Kn( double r, double T, double P );

        /* Particle diffusion coefficient [m^2/s] */
        double partDiffCoef( double r, double T, double P );

        /* Cunningham slip-flow correction factor [-] */
        double slip_flowCorrection( double Kn );

};

#endif /* PHYSFUNCTION_H_INCLUDED */
