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

#include "ForwardDecl.hpp"
#include "PhysConstant.hpp"

namespace physFunc
{

    /* Defined physical functions */

    /* Saturation pressures [Pa] */
    RealDouble pSat_H2Ol( const RealDouble T );
    RealDouble pSat_H2Os( const RealDouble T );
    RealDouble pSat_H2SO4( const RealDouble T );
    RealDouble pSat_HNO3( const RealDouble T , const RealDouble PPH2O );

    /* Derivative of the saturation pressures wrt to the temperature [Pa/K] */
    RealDouble dpSat_H2Os( const RealDouble T );

    /* Density of air [kg/m^3] */
    RealDouble rhoAir( const RealDouble T, const RealDouble P );

    /* Dynamic viscosity [kg/(m s)] */
    RealDouble dynVisc( const RealDouble T );

    /* Kinematic viscosity [m^2/s] */
    RealDouble kinVisc( const RealDouble T, const RealDouble P );

    /* Thermal speed of an air molecule/particle [m/s] */
    RealDouble thermalSpeed( const RealDouble T, \
                             const RealDouble m = physConst::M_Air );

    /* Mean free path in air [m] */
    RealDouble lambda( const RealDouble T, const RealDouble P );

    /* Mass of spherical particle [kg] */
    RealDouble mass_sphere( const RealDouble r, const RealDouble rho );
    
    /* Terminal fall speed of a spherical particle [m/s] */
    RealDouble vFall( const RealDouble r, const RealDouble rho, \
                      const RealDouble T, const RealDouble P );

    /* Knudsen number [-] */
    RealDouble Kn( const RealDouble r, const RealDouble T, \
                   const RealDouble P );

    /* Particle diffusion coefficient [m^2/s] */
    RealDouble partDiffCoef( const RealDouble r, const RealDouble T, \
                             const RealDouble P );

    /* Cunningham slip-flow correction factor [-] */
    RealDouble slip_flowCorrection( const RealDouble Kn );
    
    /* Mean free path of particles in air [m] */
    RealDouble lambda_p( const RealDouble r, const RealDouble m, \
                         const RealDouble T, const RealDouble P );
    
    /* Mean distance [m] from the center of a sphere reach by particles leaving */
    RealDouble delta_p( const RealDouble r, const RealDouble m, \
                        const RealDouble T, const RealDouble P );
    
    /* Particle Reynolds number [-] */
    RealDouble Reynolds_p( const RealDouble r, const RealDouble rho, \
                           const RealDouble T, const RealDouble P );
    
    /* Particle Schmidt number [-] */
    RealDouble Schmidt_p( const RealDouble r, const RealDouble T, \
                          const RealDouble P );

    /* Particle Stokes number [-] */
    RealDouble Stokes_p( const RealDouble r_1, const RealDouble rho_1, \
                         const RealDouble r_2, const RealDouble rho_2, \
                         const RealDouble T, const RealDouble P );
    
    /* Aggregation efficiency for liquid particles [-] */
    RealDouble E_agg( const RealDouble r_1, const RealDouble rho_1, \
                      const RealDouble r_2, const RealDouble rho_2, \
                      const RealDouble T, const RealDouble P );
    
    /* H2O gas phase diffusion coefficient in [m^2/s] */
    RealDouble DiffCoef_H2O( const RealDouble T, const RealDouble P );
    
    /* H2SO4 gas phase diffusion coefficient in [m^2/s] */
    RealDouble DiffCoef_H2SO4( const RealDouble T, const RealDouble P );
    
    /* HNO3 gas phase diffusion coefficient in [m^2/s] */
    RealDouble DiffCoef_HNO3( const RealDouble T, const RealDouble P );

    /* Corrected H2O gas phase diffusion coefficient in [m^2/s] */
    RealDouble CorrDiffCoef_H2O( const RealDouble r, const RealDouble T, \
                                 const RealDouble P );
    
    /* Corrected H2SO4 gas phase diffusion coefficient in [m^2/s] */
    RealDouble CorrDiffCoef_H2SO4( const RealDouble r, const RealDouble T, \
                                   const RealDouble P );

    /* Corrected HNO3 gas phase diffusion coefficient in [m^2/s] */
    RealDouble CorrDiffCoef_HNO3( const RealDouble r, const RealDouble T, \
                                  const RealDouble P );
    
    /* Thermal conductivity of dry air in [J/(msK)] */
    RealDouble ThermalCond( const RealDouble r, const RealDouble T, \
                            const RealDouble P );

    /* Latent heat of sublimation of water vapor in [J/kg] */
    RealDouble LHeatSubl_H2O( const RealDouble T );

    /* Kelvin factor [-] */
    RealDouble Kelvin( const RealDouble r );

    /* RH Field */
    Vector_2D RHi_Field(const Vector_2D& H2O, const Vector_2D& T, const Vector_1D& P);

}

#endif /* PHYSFUNCTION_H_INCLUDED */
