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

namespace physFunc
{

    /* Defined physical functions */

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

    /* Thermal speed of an air molecule/particle [m/s] */
    double thermalSpeed( double T, double m = physConst::M_Air );

    /* Mean free path in air [m] */
    double lambda( double T, double P );

    /* Mass of spherical particle [kg] */
    double mass_sphere( double r, double rho );
    
    /* Terminal fall speed of a spherical particle [m/s] */
    double vFall( double r, double rho, double T, double P );

    /* Knudsen number [-] */
    double Kn( double r, double T, double P );

    /* Particle diffusion coefficient [m^2/s] */
    double partDiffCoef( double r, double T, double P );

    /* Cunningham slip-flow correction factor [-] */
    double slip_flowCorrection( double Kn );
    
    /* Mean free path of particles in air [m] */
    double lambda_p( double r, double m, double T, double P );
    
    /* Mean distance [m] from the center of a sphere reach by particles leaving */
    double delta_p( double r, double m, double T, double P );

}

#endif /* PHYSFUNCTION_H_INCLUDED */
