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

#include "ForwardsDecl.hpp"
#include "PhysConstant.hpp"

namespace physFunc
{

    /* Defined physical functions */

    /* Saturation pressures [Pa] */
    RealDouble pSat_H2Ol( RealDouble T );
    RealDouble pSat_H2Os( RealDouble T );
    RealDouble pSat_H2SO4( RealDouble T );
    RealDouble pSat_HNO3( RealDouble T , RealDouble PPH2O );

    /* Density of air [kg/m^3] */
    RealDouble rhoAir( RealDouble T, RealDouble P );

    /* Dynamic viscosity [kg/(m s)] */
    RealDouble dynVisc( RealDouble T );

    /* Kinematic viscosity [m^2/s] */
    RealDouble kinVisc( RealDouble T, RealDouble P );

    /* Thermal speed of an air molecule/particle [m/s] */
    RealDouble thermalSpeed( RealDouble T, RealDouble m = physConst::M_Air );

    /* Mean free path in air [m] */
    RealDouble lambda( RealDouble T, RealDouble P );

    /* Mass of spherical particle [kg] */
    RealDouble mass_sphere( RealDouble r, RealDouble rho );
    
    /* Terminal fall speed of a spherical particle [m/s] */
    RealDouble vFall( RealDouble r, RealDouble rho, RealDouble T, RealDouble P );

    /* Knudsen number [-] */
    RealDouble Kn( RealDouble r, RealDouble T, RealDouble P );

    /* Particle diffusion coefficient [m^2/s] */
    RealDouble partDiffCoef( RealDouble r, RealDouble T, RealDouble P );

    /* Cunningham slip-flow correction factor [-] */
    RealDouble slip_flowCorrection( RealDouble Kn );
    
    /* Mean free path of particles in air [m] */
    RealDouble lambda_p( RealDouble r, RealDouble m, RealDouble T, RealDouble P );
    
    /* Mean distance [m] from the center of a sphere reach by particles leaving */
    RealDouble delta_p( RealDouble r, RealDouble m, RealDouble T, RealDouble P );
    
    /* Particle Reynolds number [-] */
    RealDouble Reynolds_p( RealDouble r, RealDouble rho, RealDouble T, RealDouble P );
    
    /* Particle Schmidt number [-] */
    RealDouble Schmidt_p( RealDouble r, RealDouble T, RealDouble P );

    /* Particle Stokes number [-] */
    RealDouble Stokes_p( RealDouble r_1, RealDouble rho_1, RealDouble r_2, RealDouble rho_2, RealDouble T, RealDouble P );
    
    /* Aggregation efficiency for liquid particles [-] */
    RealDouble E_agg( RealDouble r_1, RealDouble rho_1, RealDouble r_2, RealDouble rho_2, RealDouble T, RealDouble P );

}

#endif /* PHYSFUNCTION_H_INCLUDED */
